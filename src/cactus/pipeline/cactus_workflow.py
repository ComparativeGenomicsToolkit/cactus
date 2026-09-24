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

    # One number: the peak with jemalloc page retention OFF.  That peak is caf's, not bar's --
    # 350.0 GiB across five salamander Anc3 runs with different filters, windows and core counts
    # -- so it answers to input size and to nothing bar does.  It is a power of the input, well
    # under linear: the two workloads measured for this quantity under the same configuration
    # are MammalsAnc0 at 122.5 GiB on 4.3 GB of input and salamander Anc3 at 350.0 on 122 GB,
    # which is 28x the input for 2.9x the peak, an exponent of 0.31.
    #
    # Fitted to 575 VGP alignments from the 577way run.  Those ran in Feb-Mar 2026, before page
    # retention existed and with bar's nested parallelism still in, which inflated peaks ~2.3x,
    # so the fitted 220 is scaled down -- to 175, which is where the three retention-off peaks
    # since measured under the shipped configuration (huge pages on) put it: MammalsAnc0 122.5
    # GiB on 4.3 GB of input (2.2x under the estimate), AnuraAnc6 265.1 on 6.2 GB (1.14x) and
    # salamander Anc3 350.0 on 122 GB (2.1x).  AnuraAnc6 is the thin one and is what sets the
    # coefficient: it has flowers of up to 847 ends against MammalsAnc0's 208, so it peaks 2.2x
    # higher on 1.5x the input, and no power of the input can follow that -- only the margin
    # can cover it.  At the previous 147 it was under, 254 against 265, and would have been
    # OOM-killed on its own request.  Conservative otherwise on purpose: an under-estimate costs
    # the alignment, an over-estimate costs queue time.
    mem_coef = getOptionalAttrib(cons_node, 'memory_coefficient_gb', typeFn=float, default=175.0)
    mem_exp = getOptionalAttrib(cons_node, 'memory_input_exponent', typeFn=float, default=0.30)
    input_gb = (paf_id.size + total_sequence_size) / 1e9
    # Below 0.02 GB of input, under every point in the fit, taper to zero rather than hand an
    # evolver-sized test the power law's extrapolation.
    ramp = min(1.0, input_gb / 0.02)
    estimate = int(mem_coef * (input_gb ** mem_exp) * ramp * 2**30) if input_gb > 0 else 0
    # A POA aligner (abPOA or minipoa) needs a table even for tiny alignments
    bar_node = findRequiredNode(config_node, 'bar')
    base_aligner = getOptionalAttrib(bar_node, 'baseAligner', typeFn=str, default=None)
    if base_aligner is None:
        base_aligner = 'abpoa' if getOptionalAttrib(bar_node, 'partialOrderAlignment', typeFn=bool, default=True) else 'pecan'
    if base_aligner != 'pecan':
        estimate = max(estimate, int(4e9))

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

    RealtimeLogger.info('Estimating cactus_consolidated({}) memory without page retention as {} from {} sequences with total-sequence-size {} and paf-size {}'.format(
        name, bytes2human(estimate), len(seq_id_map), bytes2human(total_sequence_size), paf_id.size))

    max_system_memory = ConfigWrapper(config_node).getSystemMemory()

    # Page retention (see <consolidated retain_pages>) makes bar 2-2.5x faster, and its peak
    # follows cumulative churn rather than the working set, so it is not estimated -- it is given
    # room, and the release guard (retain_pages_release_pct) caps it.  The request is
    # memory_retain_multiple times the estimate, which is what retention is allowed to grow into;
    # the guard hands the pages back before it gets there, so the request does not have to predict
    # the peak, only leave the guard somewhere sensible to sit.
    #
    # Retention is only attempted with a whole further estimate of headroom beyond that request.
    # That is what keeps salamander Anc3 out: it would ask 1242 GiB of a 1700 GiB ceiling to buy a
    # caf 10% faster and a guard that fires either way, and a job that would take essentially the
    # machine should not be gambling on a peak that has no ceiling of its own.
    retain_pages = cons_retain_pages if cons_retain_pages is not None else getOptionalAttrib(cons_node, 'retain_pages', default='auto')
    retain_pages = str(retain_pages).lower()
    if retain_pages not in ['auto', '0', '1']:
        raise RuntimeError('<consolidated retain_pages> / --consRetainPages must be auto, 0 or 1, not {}'.format(retain_pages))
    retain_multiple = getOptionalAttrib(cons_node, 'memory_retain_multiple', typeFn=float, default=2.0)
    if retain_pages == 'auto':
        # sys.maxsize is what CACTUS_MAX_MEMORY holds when --maxMemory was not given and the
        # batch system could not be asked (see cactus_clamp_memory's setup): every one except
        # single_machine, and slurm when the node probe fails.  It is truthy, so taking it at
        # face value would make the test below always pass and every job retain.  An unknown
        # ceiling is a reason to decline, not a licence: retention is the optimisation, and a
        # 2x request that no node can satisfy pends forever.
        limits = [l for l in [max_system_memory,
                              int(os.environ['CACTUS_MAX_MEMORY']) if 'CACTUS_MAX_MEMORY' in os.environ else None]
                  if l and l < sys.maxsize]
        limit = min(limits) if limits else None
        if limit is None:
            RealtimeLogger.info('cactus_consolidated({}): how much memory a job can be given is not known here, so the pages will not be retained; pass --maxMemory to enable it'.format(name))
            retain_pages = '0'
        elif (retain_multiple + 1) * estimate > limit:
            RealtimeLogger.info('cactus_consolidated({}): a {:g}x request of {} plus a further {} of headroom does not fit the {} the job can be given, so the pages will not be retained'.format(
                name, retain_multiple, bytes2human(int(retain_multiple * estimate)), bytes2human(estimate), bytes2human(limit)))
            retain_pages = '0'
        else:
            retain_pages = '1'
    mem = int(retain_multiple * estimate) if retain_pages == '1' else estimate
    RealtimeLogger.info('cactus_consolidated({}): jemalloc page retention {}, requesting {}'.format(
        name, 'on' if retain_pages == '1' else 'off', bytes2human(mem)))

    if cons_memory is not None and cons_memory != mem:
        RealtimeLogger.info('Overriding cactus_consolidated({}) memory request of {} with {} value {} from --consMemory'.format(
            name, bytes2human(mem), 'greater' if cons_memory > mem else 'lesser', bytes2human(cons_memory)))
        mem = cons_memory

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
        # only the workflow applies partialOrderAlignmentWindowBigGenome, so the window it resolved
        # has to be written into the config bar reads.  It no longer feeds the memory estimate.
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
    # process rather than in the worker's own environment, and it travels into a container, which
    # an exported variable does not.  An existing MALLOC_CONF wins.
    #
    # Always, in both retention modes.  Retention off: 44398 -> 33124 s on salamander Anc3 for 0.9%
    # more peak.  Retention on: this used to be gated off, on the theory that a partly used 2 MB page
    # held whole for the life of the process was what OOM-killed two salamander runs -- but that was
    # retention's own unbounded growth, which the release guard now handles, and measured directly
    # on MammalsAnc0 the pair is the best cell of its round: 4580 s bar / 269 GiB against 6536 s /
    # 289 GiB without huge pages.  It also fixes what happens after the guard fires: the run then
    # continues in stock decay, and with huge pages that is the ~4100 s mode instead of the 13433 s
    # one that two salamander runs spent 50 hours in.
    env_prefix = []
    if getOptionalAttrib(findRequiredNode(config_node, 'consolidated'), 'transparent_huge_pages', typeFn=bool, default=True) \
       and 'MALLOC_CONF' not in os.environ:
        env_prefix = ['env', 'MALLOC_CONF=thp:always']

    # A way out if the estimate was wrong: the pages are handed back rather than the job being
    # OOM-killed.  Checked from caf's melting rounds and bar's poa loop.  Only with retention on.
    guard_pct = getOptionalAttrib(findRequiredNode(config_node, 'consolidated'), 'retain_pages_release_pct', typeFn=float, default=80.0)
    if str(retain_pages) == '1' and guard_pct > 0 and job.memory:
        limit_mb = int(job.memory * guard_pct / 100.0 / 2**20)
        if not env_prefix:
            env_prefix = ['env']
        env_prefix.append('CACTUS_RETENTION_OFF_MB={}'.format(limit_mb))

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
