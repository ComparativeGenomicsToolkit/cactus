#!/usr/bin/env python3

#Copyright (C) 2009-2021 by Benedict Paten, Joel Armstrong and Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

"""Script strings together all the components to make the basic pipeline for reconstruction.
"""

import os
import copy
import sys
from toil.lib.bioio import system
from toil.lib.bioio import getLogLevelString
from toil.realtimeLogger import RealtimeLogger
from toil.lib.conversions import bytes2human
from sonLib.nxnewick import NXNewick
from cactus.shared.common import makeURL
from cactus.shared.common import cactus_call
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import findRequiredNode, getOptionalAttrib
from cactus.shared.common import cactus_clamp_memory
from cactus.shared.common import cactus_walltime

############################################################
############################################################
############################################################
##The consolidate phase, which runs the setup, caf,
## bar, reference and cactus to hal algorithms in one job
## on a multi-node machine
############################################################
############################################################
############################################################

# bar is 63.5% of cactus_consolidated's time across the 576 VGP alignments, caf 20.3%,
# reference 15.4%.  The parallel part stops improving somewhere around 24 cores -- which is why
# --consCores above that buys memory rather than speed -- so a job given fewer than that, and
# only then, takes proportionally longer.  The fits themselves all ran at 64 cores, on the
# plateau.
CONS_CORE_BASELINE = 24
CONS_PARALLEL_FRACTION = 0.64

def cons_core_scale(cores, baseline=CONS_CORE_BASELINE, parallel=CONS_PARALLEL_FRACTION):
    """ how much longer cactus_consolidated takes when given `cores` rather than a full node """
    if cores and 0 < cores < baseline:
        return (1.0 - parallel) + parallel * (float(baseline) / cores)
    return 1.0

def cactus_cons_with_resources(job, tree, ancestor_event, config_node, seq_id_map, og_map, paf_id,
                               cons_cores = None, cons_memory = None, intermediate_results_url = None, chrom_name = None,
                               cons_retain_pages = None):
    ''' run cactus_consolidated as a child job, requesting resources based on input sizes '''

    cons_node = findRequiredNode(config_node, 'consolidated')
    name = chrom_name if chrom_name else ancestor_event
    outgroups = set(og_map[ancestor_event] if ancestor_event in og_map else [])
    og_size_scale = getOptionalAttrib(cons_node, 'og_size_scale_pct', typeFn=float, default=100.0)

    # compute resource requirements
    total_sequence_size = 0
    for seq_name, seq_id in seq_id_map.items():
        seq_size = seq_id.size
        if seq_name in outgroups:
            seq_size = int(seq_size * (og_size_scale / 100.))
        total_sequence_size += seq_size
        
    disk = 5 * total_sequence_size + 2 * paf_id.size

    # Peak memory saturates: it is set by the largest flowers bar builds, not by input volume.
    # AnuraAnc6 (0.44 GB paf) peaked at 685.9 GiB and salamander Anc3 (19.8 GB paf, 45x larger)
    # at 682.3 GiB.  So the estimate is a power of the input, not linear in it -- a linear term
    # fitted on the VGP range extrapolates to roughly double the truth at salamander scale.
    #
    # Fitted to 575 VGP alignments (64 cores, stock bar parameters), against the old table:
    #
    #   old seq_size table   14.3% under-provisioned, 1.7% still short at 2x, 2.37x allocated
    #   this model           14.3% under-provisioned, 0.7% still short at 2x, 1.85x allocated
    #
    # The middle column is the one that costs a run: --doubleMem retries at twice the request,
    # so an under-estimate only loses the work when 2x misses too.  The old table keyed off
    # total sequence size, which correlates just 0.27 with peak -- two of those alignments
    # peaked within 1% of each other (989 and 980 GiB) and were estimated 2150 and 373 GiB
    # purely on their sequence sizes.  Held-out check: this predicts salamander Anc3 at 1.01x
    # its measured peak, where the table gave 3.15x.
    mem_coef = getOptionalAttrib(cons_node, 'memory_coefficient_gb', typeFn=float, default=220.0)
    mem_exp = getOptionalAttrib(cons_node, 'memory_input_exponent', typeFn=float, default=0.30)
    mem_ramp_gb = getOptionalAttrib(cons_node, 'memory_ramp_gb', typeFn=float, default=0.02)
    input_gb = (paf_id.size + total_sequence_size) / 1e9
    # Below the smallest alignment in the fit (0.027 GB of input) the power law is pure
    # extrapolation and would hand an evolver-sized test tens of GiB, so taper it to zero.
    # The taper ends at 0.02 GB, under every point in the fit, so it changes none of them.
    ramp = min(1.0, input_gb / mem_ramp_gb) if mem_ramp_gb > 0 else 1.0
    mem = int(mem_coef * (input_gb ** mem_exp) * ramp * 2**30) if input_gb > 0 else 0
    # the window and core factors below scale bar's working set.  With retention off the peak is
    # caf's instead -- measured at 350.0 GiB on five salamander runs with different filters,
    # windows and core counts -- so the unscaled figure is what that path uses.
    unscaled_mem = mem

    # Memory is *not* quadratic in the poa window, whatever the poa comment says: halving it from
    # 10000 to 5000 measured 858.4 -> 639.4 GiB on salamander Anc3, a 0.745x ratio, because halving
    # the window also halves the number of windows and peak is bounded by how many run at once.
    #
    # The exponent is deliberately gentler than that ratio.  0.745 is the ratio of two *peaks*, and
    # a peak carries a large component that the window does not touch -- caf's own high water is
    # most of it -- so applying the full ratio to the whole estimate reduces the part that does not
    # scale.  At 0.43 a salamander at 64 cores was estimated 690 GiB against a measured 639 and was
    # OOM-killed in bar; the retry cost a day.  0.20 leaves ~27% instead of ~8%.  Over-estimating
    # costs queue time, under-estimating costs the whole alignment, and these are the jobs where
    # that trade is most lopsided.  Nothing is subtracted for partialOrderAlignmentMaskFilter even
    # though it matters more, because every alignment in the fit ran with it disabled.
    poa_node = findRequiredNode(config_node, 'bar').find('poa')
    poa_window = getOptionalAttrib(poa_node, 'partialOrderAlignmentWindow', typeFn=int, default=10000) if poa_node is not None else 10000

    # Giant genomes get a smaller poa window.  The window is the only bound on the DP once the
    # sequence is long and repeat-rich, and on 22 Gb salamanders halving it took bar's peak from
    # 798.5 to 368.7 GiB and doubled throughput, for 0.06 points of recall on evolver mammals.
    # Any role, not just ingroups.  Bar aligns an outgroup's sequence too, and lungfish Anc2 is the
    # case that shows it: a 40 GB outgroup against two ingroups of 1 and 2 GB, 91% of the input,
    # and the job was OOM-killed.  An ingroup test sees 2 GB there and does nothing.
    big_window = getOptionalAttrib(poa_node, 'partialOrderAlignmentWindowBigGenome', typeFn=int, default=0) if poa_node is not None else 0
    big_threshold = getOptionalAttrib(poa_node, 'partialOrderAlignmentWindowBigGenomeThreshold', typeFn=float, default=0) if poa_node is not None else 0
    if big_window > 0 and big_threshold > 0:
        biggest_genome = max([seq_id.size for seq_id in seq_id_map.values()] or [0])
        if biggest_genome >= big_threshold and big_window < poa_window:
            RealtimeLogger.info('cactus_consolidated({}): largest genome is {}, at or above the {} threshold, so the poa window drops from {} to {}'.format(
                name, bytes2human(biggest_genome), bytes2human(int(big_threshold)), poa_window, big_window))
            poa_window = big_window
    window_exp = getOptionalAttrib(cons_node, 'memory_poa_window_exponent', typeFn=float, default=0.20)
    if poa_window > 0 and poa_window != 10000:
        mem = int(mem * (poa_window / 10000.0) ** window_exp)

    # Scale with the core count in BOTH directions.  Concurrent abPOA instances now equal
    # the core count exactly (bar's nested parallelism is gone), so a job given fewer cores
    # genuinely needs less memory and should say so -- the old form only ever scaled up, so a
    # 24-core job asked for the 32-core figure.  Measured on AnuraAnc6 with the current
    # defaults: 465.5 GiB at 24 cores, 542.2 at 32, 685.9 at 64, which is ~0.65%/core against
    # a 64-core baseline (exact at 32, 9% conservative at 24).  The baseline is 64 because
    # that is where the 575 alignments the model was fitted to were run.
    core_scale_pct = getOptionalAttrib(cons_node, 'memory_core_scale_pct', typeFn=float, default=0.65)
    core_scale_baseline = getOptionalAttrib(cons_node, 'memory_core_scale_baseline', typeFn=int, default=64)
    if cons_cores and cons_cores != core_scale_baseline and core_scale_pct > 0:
        extra_cores = cons_cores - core_scale_baseline
        scale_factor = 1.0 + (extra_cores * core_scale_pct / 100.0)
        scaled_mem = int(mem * scale_factor)
        RealtimeLogger.info('Scaling cactus_consolidated({}) memory by {:.1f}% for {} cores ({:+d} against the {} baseline): {} -> {}'.format(
            chrom_name if chrom_name else ancestor_event, (scale_factor - 1) * 100, cons_cores, extra_cores, core_scale_baseline,
            bytes2human(mem), bytes2human(scaled_mem)))
        mem = scaled_mem

    # abPOA needs a table even for tiny alignments; apply the floor last so neither the window
    # nor the core scaling can push a small job below it
    if getOptionalAttrib(findRequiredNode(config_node, 'bar'), 'partialOrderAlignment', typeFn=bool, default=True):
        mem = max(mem, int(4e9))

    RealtimeLogger.info('Estimating cactus_consolidated({}) memory as {} from {} sequences with total-sequence-size {} and paf-size {} using <consolidated> configuration settings'.format(chrom_name if chrom_name else ancestor_event, bytes2human(mem), len(seq_id_map), bytes2human(total_sequence_size), paf_id.size))

    if cons_memory is not None and cons_memory != mem:
        RealtimeLogger.info('Overriding cactus_conslidated({}) memory estimate of {} with {} value {} from --consMemory'.format(
            chrom_name if chrom_name else ancestor_event, bytes2human(mem), 'greater' if cons_memory > mem else 'lesser', bytes2human(cons_memory)))
        mem = cons_memory

    max_system_memory = ConfigWrapper(config_node).getSystemMemory()

    # Whether cactus_consolidated keeps the pages jemalloc frees (see <consolidated retain_pages>).
    # The memory fit above was made with retention on, so when it is off the estimate is scaled
    # down by memory_retain_ratio.  "auto" keeps the pages unless the retained estimate exceeds
    # what the job can be given: the system memory on a single machine, or --maxMemory.
    retain_pages = cons_retain_pages if cons_retain_pages is not None else getOptionalAttrib(cons_node, 'retain_pages', default='auto')
    retain_pages = str(retain_pages).lower()
    if retain_pages not in ['auto', '0', '1']:
        raise RuntimeError('<consolidated retain_pages> / --consRetainPages must be auto, 0 or 1, not {}'.format(retain_pages))
    retain_ratio = getOptionalAttrib(cons_node, 'memory_retain_ratio', typeFn=float, default=2.0)
    retain_fraction = getOptionalAttrib(cons_node, 'memory_retain_auto_fraction', typeFn=float, default=0.5)
    if retain_pages == 'auto':
        limits = [l for l in [max_system_memory, int(os.environ['CACTUS_MAX_MEMORY']) if 'CACTUS_MAX_MEMORY' in os.environ else None] if l]
        limit = min(limits) if limits else None
        # Retention is only attempted with room to be wrong.  mem here is the retained estimate for
        # this job's core count, and 14% of VGP alignments came in over their estimate -- so asking
        # it to fit in a fraction of what the job can be given means a miss still has somewhere to
        # land.  Without the fraction, any under-estimate at the ceiling OOMs the whole alignment.
        # `budget is not None`, not `budget`: a fraction of 0 means never retain, and a bare
        # truthiness test would make 0.0 falsy and fall through to retaining every time.
        budget = limit * retain_fraction if limit is not None else None
        if budget is not None and mem > budget:
            RealtimeLogger.info('cactus_consolidated({}): the memory estimate of {} with jemalloc page retention exceeds {:g} of the {} the job can be given, so the pages will not be retained'.format(
                name, bytes2human(mem), retain_fraction, bytes2human(limit)))
            retain_pages = '0'
        else:
            retain_pages = '1'
    if retain_pages == '0' and cons_memory is None and retain_ratio > 1:
        RealtimeLogger.info('cactus_consolidated({}): without jemalloc page retention the estimate is {}, from the unscaled {} over {}'.format(
            name, bytes2human(int(unscaled_mem / retain_ratio)), bytes2human(unscaled_mem), retain_ratio))
        mem = int(unscaled_mem / retain_ratio)
    RealtimeLogger.info('cactus_consolidated({}): jemalloc page retention {}'.format(name, 'on' if retain_pages == '1' else 'off'))

    if max_system_memory and mem > max_system_memory:
        RealtimeLogger.info('Clamping cactus_conslidated({}) memory estimate of {} to maximum system memory {}'.format(
            name, bytes2human(mem), bytes2human(max_system_memory)))
        mem = max_system_memory

    # Runtime has two regimes, and which one you are in is decided by <bar bandingLimit> (what
    # --maxLen sets).  They are different in shape, not just in scale:
    #
    #   unbanded, 1 Mb (progressive)   secs = 214 * disk_gb**0.906   r = 0.64
    #   banded,  10 kb (pangenome)     secs = 2522 * disk_gb**0.251  r = 0.26
    #
    # fitted to 576 VGP 577-way alignments and to the 50 chromosome alignments of two HPRC
    # pangenomes.  Banding to 10 kb bounds bar's work per column, so the cost follows the number
    # of reference columns rather than the sequence volume -- which is why the pangenome
    # exponent is nearly flat, and why adding haplotypes barely moves it even though it moves
    # `disk` a great deal.  Applying the progressive fit to a pangenome chromosome overshoots by
    # 13-17x at the median: HPRC chr5 takes 4.1 h and would have been given 312.
    #
    # The coefficients below are the fits scaled so that cactus_walltime()'s factor covers the
    # worst residual, then divided by 4: the 2x cactus_consolidated speedup that has already
    # landed since both sets of logs, and a further 2x that is expected but NOT yet measured
    # here.  If that second 2x underdelivers the cost is small and bounded -- replaying the
    # fits against today's times, 5 of 576 progressive and 3 of 50 pangenome alignments would
    # run over, none by more than 1.5x, so a single --doubleTime retry rescues every one.  The
    # gain is not small: it takes the median progressive request from 15.6 h to 7.8 h and the
    # number of them over 24 h from 181 to 7.  Raise these two attributes if that turns out to
    # be optimistic.  `disk` is the size term because it already combines the sequence and paf
    # sizes in the proportions that drive the work (5:2).
    #
    # Only the two measured banding values are in real use, so this selects between them rather
    # than interpolating a curve through data that does not exist.  Anything between them takes
    # the unbanded model, which is the conservative side.
    banding_limit = getOptionalAttrib(findRequiredNode(config_node, 'bar'), 'bandingLimit', typeFn=int, default=0)
    banding_threshold = getOptionalAttrib(cons_node, 'walltime_banding_threshold', typeFn=int, default=100000)
    if banding_limit and banding_limit < banding_threshold:
        wt_coef = getOptionalAttrib(cons_node, 'walltime_banded_coefficient_secs', typeFn=float, default=550.0)
        wt_exp = getOptionalAttrib(cons_node, 'walltime_banded_exponent', typeFn=float, default=0.35)
    else:
        # the progressive fit's exponent was 0.906 and rounded to 0.95 here, which is near enough
        # to linear that carrying it as a tunable was not worth the knob: over sizes from 0.1 to
        # 50 GB and core counts from 8 to 64, dropping it moves 3 of 36 jobs across a partition
        # boundary, and always upward.  The banded exponent below is a different matter -- 0.35 is
        # genuinely concave, and linearising it either under-provisions the middle of the range by
        # 1.7x or over-provisions the top by 4.4x -- so that one stays.
        wt_coef = getOptionalAttrib(cons_node, 'walltime_coefficient_secs', typeFn=float, default=400.0)
        wt_exp = 1.0
    walltime_secs = wt_coef * ((disk / 1e9) ** wt_exp) if disk > 0 else 0
    # bar is 63.5% of consolidated's time across those 576 alignments, caf 20.3%, reference
    # 15.4%; the parallel part of that stops improving somewhere around 24 cores (which is why
    # --consCores above ~24 buys memory, not speed), so scale up only when a job is given fewer
    # than that.  The fit's own jobs all ran at 64 cores, i.e. already on the plateau.
    # CONS_CORE_BASELINE and CONS_PARALLEL_FRACTION rather than config attributes: they describe
    # how cactus_consolidated parallelises, not anything about a particular alignment, and the
    # three other callers of cons_core_scale already take them from the module.
    walltime_secs *= cons_core_scale(cons_cores)

    cons_job = job.addChildJobFn(cactus_cons, tree, ancestor_event, config_node, seq_id_map, og_map, paf_id,
                                 intermediate_results_url=intermediate_results_url, chrom_name=chrom_name, cores = cons_cores,
                                 memory=cactus_clamp_memory(mem), disk=disk, retain_pages=retain_pages,
                                 poa_window=poa_window, walltime=cactus_walltime(walltime_secs))
    return cons_job.rv()

def cactus_cons(job, tree, ancestor_event, config_node, seq_id_map, og_map, paf_id,
                intermediate_results_url = None, chrom_name = None, retain_pages = None,
                poa_window = None):
    ''' run cactus_consolidated '''

    # cactus_consolidated reads its settings from the config, so the resolved page retention
    # goes into the copy it is given (this job's copy of the node, so nothing else sees it)
    if retain_pages is not None or poa_window is not None:
        config_node = copy.deepcopy(config_node)
        if retain_pages is not None:
            findRequiredNode(config_node, 'consolidated').set('retain_pages', str(retain_pages))
        # the estimator resolved the window (it may have been lowered for a giant genome), and the
        # estimate it produced only holds if cactus_consolidated uses that same value
        poa_node = findRequiredNode(config_node, 'bar').find('poa')
        if poa_window is not None and poa_node is not None:
            poa_node.set('partialOrderAlignmentWindow', str(poa_window))

    # Build up a genome -> fasta map.
    work_dir = job.fileStore.getLocalTempDir()
    seq_path_map = {}
    for event, seq_id in seq_id_map.items():
        seq_path = os.path.join(work_dir, '{}.fa'.format(event))
        job.fileStore.readGlobalFile(seq_id, seq_path)
        seq_path_map[event] = seq_path

    outgroups = og_map[ancestor_event] if ancestor_event in og_map else []

    # Get the alignments file
    paf_path = os.path.join(work_dir, f'{ancestor_event}.paf')
    job.fileStore.readGlobalFile(paf_id, paf_path)

    # Split the alignments file into primary and secondary
    primary_alignment_file = os.path.join(work_dir, f'{ancestor_event}_primary.paf')
    # grep exits 1 when nothing matched, which is fine and has to be tolerated,
    # but it exits 2 when the write itself failed.  `|| true` cannot tell those
    # apart, so a full disk here used to leave a short paf and carry on.
    system(f"grep -v 'tp:A:S' {paf_path} > {primary_alignment_file} || [ $? -eq 1 ]")  # Alignments that are not-secondaries

    # Optionally parse our secondary alignments
    use_secondary_alignments = int(config_node.find("blast").attrib["outputSecondaryAlignments"])  # We should really switch to
    # the empty string being false instead of 0
    assert use_secondary_alignments == 0 or use_secondary_alignments == 1
    if use_secondary_alignments:
        secondary_alignment_file = os.path.join(work_dir, f'{ancestor_event}_secondary.paf')
        system(f"grep 'tp:A:S' {paf_path} > {secondary_alignment_file} || [ $? -eq 1 ]")  # Alignments that are secondaries

    # Optionally copy the alignments to a specified location for debug purposes
    if config_node.find("caf").attrib["writeInputAlignmentsTo"]:
        system(f'cp {paf_path} {config_node.find("caf").attrib["writeInputAlignmentsTo"]}/{ancestor_event}.paf')

    # Temporary place to store the output c2h file
    tmpHal = os.path.join(work_dir, f'{ancestor_event}.c2h')
    tmpFasta = os.path.join(work_dir, f'{ancestor_event}.c2h.fa')
    tmpRef = os.path.join(work_dir, f'{ancestor_event}.ref')

    tmpConfig = os.path.join(work_dir, f'{ancestor_event}.config.xml')
    ConfigWrapper(config_node).writeXML(tmpConfig)

    # We pass the tree and species in with a seqFile
    tmpSeqFilePath = os.path.join(work_dir, f'{ancestor_event}.seqfile')
    with open(tmpSeqFilePath, 'w') as seqFile:
        seqFile.write(f'{NXNewick().writeString(tree)}\n')
        for genome, faPath in list(seq_path_map.items()):
            # note: path must be relative for docker support
            seqFile.write(f'{genome}\t{os.path.basename(faPath)}\n')
            
    args = ["--seqFile", tmpSeqFilePath, "--logLevel", getLogLevelString(),
            "--alignments", primary_alignment_file, "--params", tmpConfig, "--outputFile", tmpHal,
            "--outputHalFastaFile", tmpFasta, "--outputReferenceFile", tmpRef, "--outgroupEvents", " ".join(outgroups),
            "--referenceEvent", ancestor_event, "--threads", str(job.cores)]
    if use_secondary_alignments:  # Optionally add the secondary alignments
        args += ["--secondaryAlignments", secondary_alignment_file]

    # jemalloc reads MALLOC_CONF once, before main, so transparent huge pages cannot be switched on
    # through mallctl the way page retention is.  Prefixing `env` puts it in front of this one
    # process rather than in the worker's own environment, which a follow-on job would inherit --
    # and it travels into the container, which an exported variable does not, since dockerCommand
    # passes no -e.  Measured on salamander Anc3 without page retention, same processor model,
    # 44398 -> 33124 seconds for 0.9% more peak (347.0 -> 350.0 GiB).  An existing MALLOC_CONF wins.
    #
    # Only when the pages are being returned to the OS.  The measurement above is a retention-off
    # pair, and retention-on is where huge pages are most dangerous: nothing is ever handed back,
    # so a partly-used 2 MB page keeps all 2 MB for the life of the process.  A salamander Anc3 run
    # with both on was OOM-killed in bar at an estimate that pre-THP measurements say should have
    # been ample, which is what this guard is for.  Worth revisiting once the pair is measured.
    env_prefix = []
    if getOptionalAttrib(findRequiredNode(config_node, 'consolidated'), 'transparent_huge_pages', typeFn=bool, default=True) \
       and str(retain_pages) == '0' and 'MALLOC_CONF' not in os.environ:
        env_prefix = ['env', 'MALLOC_CONF=thp:always']

    # A way out if the estimate was wrong: bar hands the pages back rather than being OOM-killed.
    # Only with retention on -- with it off there is nothing to hand back.
    guard_pct = getOptionalAttrib(findRequiredNode(config_node, 'consolidated'), 'retain_pages_release_pct', typeFn=float, default=65.0)
    if str(retain_pages) == '1' and guard_pct > 0 and job.memory:
        limit_mb = int(job.memory * guard_pct / 100.0 / 2**20)
        if not env_prefix:
            env_prefix = ['env']
        env_prefix.append('CACTUS_BAR_RETENTION_OFF_MB={}'.format(limit_mb))

    messages = cactus_call(check_output=True, returnStdErr=True,
                           realtimeStderrPrefix=f'cactus_consolidated({chrom_name if chrom_name else ancestor_event})',
                           parameters=env_prefix + ["cactus_consolidated"] + args,
                           work_dir=work_dir,
                           job_memory=job.memory)[1]  # Get just the standard error output

    # if cactus was run with --realTimeLogging, cactus_call will print out conslidated's stderr messages as they happen
    # (and not return anything)
    # otherwise, if run without --realTimeLogging, cactus_call will return (but not print) stderr messages
    if messages:
        job.fileStore.logToMaster(f"cactus_consolidated event:{ancestor_event}\n{messages}")  # Log the messages
    else:
        job.fileStore.logToMaster("Ran cactus consolidated okay")

    # Write the temporary output files to the final output
    # At top level--have the final .c2h file
    halID = job.fileStore.writeGlobalFile(tmpHal)
    fastaID = job.fileStore.writeGlobalFile(tmpFasta)
    referenceID = job.fileStore.writeGlobalFile(tmpRef)

    if intermediate_results_url is not None:
        # The user requested to keep the c2h files in a separate place. Export it there.
        url = intermediate_results_url + ".c2h"
        job.fileStore.exportFile(halID, makeURL(url))

        # The user requested to keep the hal fasta files in a separate place. Export it there.
        url = intermediate_results_url + ".hal.fa"
        job.fileStore.exportFile(fastaID, makeURL(url))

        # The user requested to keep the reference fasta files in a separate place. Export it there.
        url = intermediate_results_url + ".reference.fa"
        job.fileStore.exportFile(referenceID, makeURL(url))

    return (ancestor_event, halID, fastaID, referenceID)


if __name__ == '__main__':
    runCactusWorkflow(sys.argv)
