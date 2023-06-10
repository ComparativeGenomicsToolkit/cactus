#!/usr/bin/env python3

"""This script will take the output of cactus-align-batch (which was run on the output of cactus-graphmap-split).  This would be
a set of .vg files, one for each chromosome.  It will do the following in order to make a single output graph
- clip unmapped regions out of each graph
- join the ids of the clipped graph
- convert each clipped graph to gfa
- merge the gfas into a whole-genome graph
- convert the gfa into a gbwt/gbwt-graph pair
- export the gbwt to xg
- export the xg/gbwt to vcf

so the final output will be
- set of clipped vg graphs
- graph.gfa.gz
- graph.gbwt
- graph.gg
- graph.xg
- graph.snarls
- graph.vcf.gz
- graph.vcf.gz.tbi
"""
import os, sys
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.common import unzip_gz, write_s3
from cactus.shared.version import cactus_commit
from cactus.preprocessor.fileMasking import get_mask_bed_from_fasta
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from cactus.refmap.cactus_minigraph import check_sample_names

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory, getTempFile, catFiles
    
def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("--vg", required=True, nargs='+',  help = "Input vg files (PackedGraph or HashGraph format)")
    parser.add_argument("--hal", nargs='+', default = [], help = "Input hal files (for merging)")    
    parser.add_argument("--outDir", required=True, type=str, help = "Output directory")
    parser.add_argument("--outName", required=True, type=str, help = "Basename of all output files")
    parser.add_argument("--reference", required=True, nargs='+', type=str, help = "Reference event name(s). The first will be the \"true\" reference and will be left unclipped and uncollapsed. It also should have been used with --reference in all upstream commands. Other names will be promoted to reference paths in vg")

    graphmap_join_options(parser)
        
    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--latest", dest="latest", action="store_true",
                        help="Use the latest version of the docker container "
                        "rather than pulling one matching this version of cactus")
    parser.add_argument("--containerImage", dest="containerImage", default=None,
                        help="Use the the specified pre-built containter image "
                        "rather than pulling one from quay.io")
    parser.add_argument("--binariesMode", choices=["docker", "local", "singularity"],
                        help="The way to run the Cactus binaries", default=None)

    options = parser.parse_args()

    setupBinaries(options)
    set_logging_from_options(options)
    enableDumpStack()

    if options.outDir and not options.outDir.startswith('s3://'):
        if not os.path.isdir(options.outDir):
            os.makedirs(options.outDir)

    #make sure our options make sense and fill in sensible defaults
    options = graphmap_join_validate_options(options)
        
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()
    graphmap_join(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-graphmap-join has finished after {} seconds".format(run_time))

def graphmap_join_options(parser):
    """ we share these options with cactus-pangenome """
    parser.add_argument("--clip", type=int, default=10000, help = "Generate clipped graph by removing anything longer than this amount that is unaligned to the underlying minigraph. Set to 0 to disable (must also set --filter 0 as well). [default=10000]")
    
    parser.add_argument("--filter", type=int, default=2, help = "Generate a frequency filtered graph (from the clipped graph) by removing any sequence present in fewer than this many sequences. Set to 0 to disable. [default=2]")

    parser.add_argument("--gfa", nargs='*', default=None, help = "Produce a GFA for given graph type(s) if specified. Valid types are 'full', 'clip', and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space. [--gfa clip assumed by default]")
    
    parser.add_argument("--gbz", nargs='*', default=None, help = "Generate GBZ/snarls indexes for the given graph type(s) if specified. Valid types are 'full', 'clip' and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space. --giraffe will also produce these (and other) indexes")

    parser.add_argument("--chrom-vg", nargs='*', default=None, help = "Produce a directory of chromosomal graphs is vg format for the graph type(s) specified. Valid typs are 'full', 'clip' and 'filter'. If no type specified 'clip' will be used ('full' used if clipping disabled). Multiple types can be provided separated by a space.  The output will be the <outDir>/<outName>.chroms/ directory")
    
    parser.add_argument("--vcf", nargs='*', default=None, help = "Generate a VCF from the given graph type(s). Valid types are 'full', 'clip' and 'filter'. If no type specified, 'clip' will be used ('full' used if clipping disabled). Multipe types can be provided separated by space")
    parser.add_argument("--vcfReference", nargs='+', default=None, help = "If multiple references were provided with --reference, this option can be used to specify a subset for vcf creation with --vcf. By default, --vcf will create VCFs for the first reference only")
    parser.add_argument("--vcfbub", type=int, default=100000, help = "Use vcfbub to flatten nested sites (sites with reference alleles > this will be replaced by their children)). Setting to 0 will disable, only prudcing full VCF [default=100000].")
    
    parser.add_argument("--giraffe", nargs='*', default=None, help = "Generate Giraffe (.dist, .min) indexes for the given graph type(s). Valid types are 'full', 'clip' and 'filter'. If not type specified, 'filter' will be used (will fall back to 'clip' than full if filtering, clipping disabled, respectively). Multiple types can be provided seperated by a space")
    parser.add_argument("--indexCores", type=int, default=None, help = "cores for general indexing and VCF constructions (defaults to the same as --maxCores)")

def graphmap_join_validate_options(options):
    """ make sure the options make sense and fill in sensible defaults """
    if options.hal and len(options.hal) != len(options.vg):
        raise RuntimeError("If --hal and --vg should specify the same number of files")

    # apply cpu override
    if options.batchSystem.lower() in ['single_machine', 'singleMachine']:
        if not options.indexCores:
            options.indexCores = sys.maxsize
        options.indexCores = min(options.indexCores, max(1, cactus_cpu_count() - 1), int(options.maxCores) if options.maxCores else sys.maxsize)
    else:
        if not options.indexCores:
            raise RuntimeError("--indexCores required run *not* running on single machine batch system")
    
    # sanity check the workflow options and apply defaults
    if options.filter and not options.clip:
        raise RuntimeError('--filter cannot be used without also disabling --clip.')

    # check the reference name suffix
    check_sample_names(options.reference, options.reference[0])

    if not options.gfa:
        options.gfa = ['clip'] if options.clip else ['full']
    options.gfa = list(set(options.gfa))
    for gfa in options.gfa:
        if gfa not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --gfa: {}. Must be one of {clip, filter, full}'.format(gfa))
        if gfa == 'clip' and not options.clip:
            raise RuntimError('--gfa cannot be set to clip since clipping is disabled')
        if gfa == 'filter' and not options.filter:
            raise RuntimeError('--gfa cannot be set to filter since filtering is disabled')

    if options.gbz == []:
        options.gbz = ['clip'] if options.clip else ['full']
    options.gbz = list(set(options.gbz)) if options.gbz else []
    for gbz in options.gbz:
        if gbz not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --gbz: {}. Must be one of {clip, filter, full}'.format(gbz))
        if gbz == 'clip' and not options.clip:
            raise RuntimError('--gbz cannot be set to clip since clipping is disabled')
        if gbz == 'filter' and not options.filter:
            raise RuntimeError('--gbz cannot be set to filter since filtering is disabled')

    if options.chrom_vg == []:
        options.chrom_vg = ['clip'] if options.clip else ['full']
    options.chrom_vg = list(set(options.chrom_vg)) if options.chrom_vg else []
    for chrom_vg in options.chrom_vg:
        if chrom_vg not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --chrom-vg: {}. Must be one of {clip, filter, full}'.format(chrom_vg))
        if chrom_vg == 'clip' and not options.clip:
            raise RuntimError('--chrom-vg cannot be set to clip since clipping is disabled')
        if chrom_vg == 'filter' and not options.filter:
            raise RuntimeError('--chrom-vg cannot be set to filter since filtering is disabled')        
    
    if options.vcf == []:
        options.vcf = ['clip'] if options.clip else ['full']
    options.vcf = list(set(options.vcf)) if options.vcf else []    
    for vcf in options.vcf:
        if vcf not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --vcf: {}. Must be one of {clip, filter, full}'.format(vcf))
        if vcf == 'clip' and not options.clip:
            raise RuntimError('--vcf cannot be set to clip since clipping is disabled')
        if vcf == 'filter' and not options.filter:
            raise RuntimeError('--vcf cannot be set to filter since filtering is disabled')

    # if someone used --vcfReference they probably want --vcf too
    if options.vcfReference and not options.vcf:
        raise RuntimeError('--vcfReference cannot be used without --vcf')
        
    if not options.vcfReference:
        options.vcfReference = [options.reference[0]]
    else:
        for vcfref in options.vcfReference:
            if vcfref not in options.reference:
                raise RuntimeError('--vcfReference {} invalid because {} was not specified as a --reference'.format(vcfref, vcfref))

    if options.giraffe == []:
        options.giraffe = ['filter'] if options.filter else ['clip'] if options.clip else ['full']
    options.giraffe = list(set(options.giraffe)) if options.giraffe else []        
    for giraffe in options.giraffe:
        if giraffe not in ['clip', 'filter', 'full']:
            raise RuntimeError('Unrecognized value for --giraffe: {}. Must be one of {clip, filter, full}'.format(giraffe))
        if giraffe == 'clip' and not options.clip:
            raise RuntimError('--giraffe cannot be set to clip since clipping is disabled')
        if giraffe == 'filter' and not options.filter:
            raise RuntimeError('--giraffe cannot be set to filter since filtering is disabled')

    # Prevent some useless compute due to default param combos
    if options.clip and 'clip' not in options.gfa + options.gbz + options.chrom_vg + options.vcf + options.giraffe \
       and 'filter' not in options.gfa + options.gbz + options.chrom_vg + options.vcf + options.giraffe:
        options.clip = None
    if options.filter and 'filter' not in options.gfa + options.gbz + options.chrom_vg + options.vcf + options.giraffe:
        options.filter = None

    return options
    
def graphmap_join(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            wf_output = toil.restart()
        else:
            options.cactusDir = getTempDirectory()

            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals()
                
            # load up the vgs
            vg_ids = []
            for vg_path in options.vg:
                logger.info("Importing {}".format(vg_path))
                vg_ids.append(toil.importFile(makeURL(vg_path)))
                
            # load up the hals
            hal_ids = []
            for hal_path in options.hal:
                logger.info("Importing {}".format(hal_path))
                hal_ids.append(toil.importFile(makeURL(hal_path)))

            # run the workflow
            wf_output = toil.start(Job.wrapJobFn(graphmap_join_workflow, options, config, vg_ids, hal_ids))
                
        #export the split data
        export_join_data(toil, options, wf_output[0], wf_output[1], wf_output[2], wf_output[3], wf_output[4])

def graphmap_join_workflow(job, options, config, vg_ids, hal_ids):

    root_job = Job()
    job.addChild(root_job)

    # make sure reference doesn't have a haplotype suffix, as it will have been changed upstream
    ref_base, ref_ext = os.path.splitext(options.reference[0])
    assert len(ref_base) > 0
    if ref_ext:
        assert len(ref_ext) > 1 and ref_ext[1:].isnumeric()
        options.reference[0] = ref_base
        if options.vcfReference:
            for i in range(len(options.vcfReference)):
                if options.vcfReference[i] == ref_base + ref_ext:
                    options.vcfReference[i] = ref_base

    # run the "full" phase to do pre-clipping stuff
    full_vg_ids = []
    assert len(options.vg) == len(vg_ids)
    for vg_path, vg_id in zip(options.vg, vg_ids):
        full_job = Job.wrapJobFn(clip_vg, options, config, vg_path, vg_id, 'full',
                                disk=vg_id.size * 2, memory=vg_id.size * 4)
        root_job.addChild(full_job)
        full_vg_ids.append(full_job.rv(0))
    prev_job = root_job
    
    # join the ids
    join_job = prev_job.addFollowOnJobFn(join_vg, options, config, full_vg_ids,
                                         disk=sum([f.size for f in vg_ids]))
    full_vg_ids = [join_job.rv(i) for i in range(len(vg_ids))]
    prev_job = join_job

    # run the "clip" phase to do the clip-vg clipping
    clip_vg_ids = []
    clipped_stats = None 
    if options.clip or options.filter:
        clip_root_job = Job()
        prev_job.addFollowOn(clip_root_job)
        clip_vg_stats = []
        assert len(options.vg) == len(full_vg_ids) == len(vg_ids)
        for vg_path, vg_id, input_vg_id in zip(options.vg, full_vg_ids, vg_ids):
            clip_job = Job.wrapJobFn(clip_vg, options, config, vg_path, vg_id, 'clip',
                                     disk=input_vg_id.size * 2, memory=input_vg_id.size * 4)
            clip_root_job.addChild(clip_job)
            clip_vg_ids.append(clip_job.rv(0))
            clip_vg_stats.append(clip_job.rv(1))
 
        # join the stats
        clipped_stats = clip_root_job.addFollowOnJobFn(cat_stats, clip_vg_stats).rv()
        prev_job = clip_root_job

    # run the "filter" phase to do the vg clip clipping
    filter_vg_ids = []
    if options.filter:
        filter_root_job = Job()
        prev_job.addFollowOn(filter_root_job)
        assert len(options.vg) == len(clip_vg_ids) == len(vg_ids)
        for vg_path, vg_id, input_vg_id in zip(options.vg, clip_vg_ids, vg_ids):
            filter_job = filter_root_job.addChildJobFn(vg_clip_vg, options, config, vg_path, vg_id, 
                                                       disk=input_vg_id.size * 2)
            filter_vg_ids.append(filter_job.rv())
        prev_job = filter_root_job

    # set up our whole-genome output
    out_dicts = []

    # optional hal merge
    if hal_ids:
        hal_merge_job = job.addChildJobFn(merge_hal, options, hal_ids,
                                          cores = 1,
                                          disk=sum(f.size for f in hal_ids) * 2)
        hal_id_dict = hal_merge_job.rv()
        out_dicts.append(hal_id_dict)

    workflow_phases = [('full', full_vg_ids, join_job)]
    if options.clip:
        workflow_phases.append(('clip', clip_vg_ids, clip_root_job))
    if options.filter:
        workflow_phases.append(('filter', filter_vg_ids, filter_root_job))        
    for workflow_phase, phase_vg_ids, phase_root_job in workflow_phases:
        
        # make a gfa for each
        gfa_root_job = Job()
        phase_root_job.addFollowOn(gfa_root_job)
        gfa_ids = []
        current_out_dict = None
        if workflow_phase in options.gfa + options.gbz + options.vcf + options.giraffe:
            assert len(options.vg) == len(phase_vg_ids) == len(vg_ids)
            for vg_path, vg_id, input_vg_id in zip(options.vg, phase_vg_ids, vg_ids):
                gfa_job = gfa_root_job.addChildJobFn(vg_to_gfa, options, config, vg_path, vg_id,
                                                     disk=input_vg_id.size * 5)
                gfa_ids.append(gfa_job.rv())

            # merge up the gfas and make the various vg indexes
            gfa_merge_job = gfa_root_job.addFollowOnJobFn(make_vg_indexes, options, config, gfa_ids,
                                                          tag=workflow_phase + '.',
                                                          cores=options.indexCores,
                                                          disk=sum(f.size for f in vg_ids) * 5)
            out_dicts.append(gfa_merge_job.rv())
            prev_job = gfa_merge_job
            current_out_dict = gfa_merge_job.rv()

        # optional vcf
        if workflow_phase in options.vcf:
            for vcf_ref in options.vcfReference:
                vcftag = vcf_ref + '.' + workflow_phase if vcf_ref != options.reference[0] else workflow_phase
                deconstruct_job = prev_job.addFollowOnJobFn(make_vcf, config, options.outName, vcf_ref, current_out_dict,
                                                            max_ref_allele=options.vcfbub,
                                                            tag=vcftag + '.', ref_tag = workflow_phase + '.',
                                                            cores=options.indexCores,
                                                            disk = sum(f.size for f in vg_ids) * 2)
                out_dicts.append(deconstruct_job.rv())                

        # optional giraffe
        if workflow_phase in options.giraffe:
            giraffe_job = gfa_merge_job.addFollowOnJobFn(make_giraffe_indexes, options, config, current_out_dict,
                                                         tag=workflow_phase + '.',
                                                         cores=options.indexCores,
                                                         disk = sum(f.size for f in vg_ids) * 4)
            out_dicts.append(giraffe_job.rv())

    
    return full_vg_ids, clip_vg_ids, clipped_stats, filter_vg_ids, out_dicts

def clip_vg(job, options, config, vg_path, vg_id, phase):
    """ run clip-vg 
    """
    assert phase in ['full', 'clip']
    
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)

    clipped_path = vg_path + '.clip'
    clipped_bed_path = vg_path + '.clip.bed'

    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")

    # optional gfaffixing -- only happens in full
    if phase == 'full' and getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "gfaffix", typeFn=bool, default=True):
        normalized_path = vg_path + '.gfaffixed'
        gfa_in_path = vg_path + '.gfa'
        gfa_out_path = normalized_path + '.gfa'
        # GFAffix's --dont_collapse option doesn't work on W-lines, so we strip them for now with -W
        cactus_call(parameters=['vg', 'convert', '-W', '-f', vg_path], outfile=gfa_in_path)
        fix_cmd = ['gfaffix', gfa_in_path, '--output_refined', gfa_out_path, '--check_transformation']
        if options.reference:
            fix_cmd += ['--dont_collapse', options.reference[0] + '*']
        cactus_call(parameters=fix_cmd)
        # GFAffix strips the header, until this is fixed we need to add it back (for the RS tag)
        gfa_header = cactus_call(parameters=['head', '-1', gfa_in_path], check_output=True).strip()
        cactus_call(parameters=['sed', '-i', gfa_out_path, '-e', '1s/.*/{}/'.format(gfa_header)])
        # Come back from gfa to vg
        conv_cmd = [['vg', 'convert', '-g', '-p', gfa_out_path]]
        # GFAFfix doesn't unchop, so we do that in vg after
        conv_cmd.append(['vg', 'mod', '-u', '-'])
        cactus_call(parameters=conv_cmd, outfile=normalized_path)        
        vg_path = normalized_path

    # run clip-vg no matter what, but we don't actually remove anything in full
    cmd = []
    clip_vg_cmd = ['clip-vg', vg_path, '-f']

    # don't clip the (first) reference, also enforces it's in forward orientation (for rGFA)
    clip_vg_cmd += ['-e', options.reference[0]]

    # our vg file has minigraph sequences -- we'll filter them out, along with any nodes
    # that don't appear in a non-minigraph path
    clip_vg_cmd += ['-d', graph_event]
    if phase == 'full':
        # but... in the full graph we'll leave the minigraph path fragments that are aligned to anything else in the vg's
        # since we'll need them to run the actual clipping
        clip_vg_cmd += ['-L']

    if phase == 'clip':
        if options.clip:
            clip_vg_cmd += ['-u', str(options.clip)]
            if getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "clipNonMinigraph", typeFn=bool, default=True):
                clip_vg_cmd += ['-a', graph_event]
            clip_vg_cmd += ['-o', clipped_bed_path]

    cmd.append(clip_vg_cmd)
    if options.reference:
        # GFAFfix can leave uncovered nodes with --dont_collapse.  We filter out here so they dont cause trouble later
        # Also: any kind of clipping or path dropping can leave uncovered edges, so we remove them with vg clip        
        clip_cmd = ['vg', 'clip', '-d', '1', '-']
        for ref in options.reference:
            clip_cmd += ['-P', ref]
        cmd.append(clip_cmd)
        if phase == 'clip' and getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "removeStubs", typeFn=bool, default=True):
            # todo: could save a little time by making vg clip smart enough to do two things at once
            stub_cmd = ['vg', 'clip', '-s', '-']
            for ref in options.reference:
                stub_cmd += ['-P', ref]

            # todo: do we want to add the minigraph prefix to keep stubs from minigraph? but I don't think it makes stubs....
            cmd.append(stub_cmd)

    # and we sort by id on the first go-around
    if phase == 'full':
        cmd.append(['vg', 'ids', '-s', '-'])
        
    cactus_call(parameters=cmd, outfile=clipped_path)

    # worth it
    cactus_call(parameters=['vg', 'validate', clipped_path])

    if phase == 'full':
        job.fileStore.deleteGlobalFile(vg_id)
        
    # keep some stats of the chromosomal vg file:
    # Path Lengths
    if phase == 'clip':
        chr_name = os.path.splitext(os.path.basename(vg_path))[0]
        path_stats_path = vg_path + '.path-stats.tsv'
        cactus_call(parameters=[['vg', 'paths', '-E', '-v', clipped_path],
                                ['awk', '{{print "{}\t" $0}}'.format(chr_name)]],
                    outfile=path_stats_path)
        # Nodes, edges and total length
        graph_stats_path = vg_path + '.graph-stats.tsv'
        cactus_call(parameters=[['vg', 'stats', '-l', '-z', clipped_path],
                                ['awk', '{{print "{}\t" $0}}'.format(chr_name)]],
                    outfile=graph_stats_path)
        # Stick the contig identifier onto the clipped regions bed
        clipped_bed_chr_path = clipped_bed_path + '.named'
        cactus_call(parameters=['awk',  '{{print $0 "\t{}"}}'.format(chr_name), clipped_bed_path],
                    outfile=clipped_bed_chr_path)
        sample_stats_path = vg_path + '.sample-stats.tsv'
        sample_stats = {}
        with open(path_stats_path, 'r') as path_stats_file:
            for line in path_stats_file:
                toks = line.split()
                contig_name = toks[1]
                contig_length = int(toks[2])
                sample_name = contig_name if '#' not in contig_name else '#'.join(contig_name.split('#')[:2])
                if sample_name not in sample_stats:
                    sample_stats[sample_name] = contig_length
                else:
                    sample_stats[sample_name] += contig_length
        with open(sample_stats_path, 'w') as sample_stats_file:
            for sample in sorted(sample_stats.keys()):
                sample_stats_file.write("{}\t{}\t{}\n".format(chr_name, sample, sample_stats[sample]))

        out_stats = { 'clip-stats.bed' : job.fileStore.writeGlobalFile(clipped_bed_chr_path),
                      'path-stats.tsv' : job.fileStore.writeGlobalFile(path_stats_path),
                      'graph-stats.tsv' : job.fileStore.writeGlobalFile(graph_stats_path),
                      'sample-stats.tsv' : job.fileStore.writeGlobalFile(sample_stats_path) }
    else:
        out_stats = None
    return job.fileStore.writeGlobalFile(clipped_path), out_stats

def vg_clip_vg(job , options, config, vg_path, vg_id):
    """ run vg clip, chaining multiple invocations if desired.
    """
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    clipped_path = vg_path + '.clip'

    clip_cmd = ['vg', 'clip', vg_path, '-d', str(options.filter)]
    for ref in options.reference:
        clip_cmd += ['-P', ref]
    min_fragment = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "minFilterFragment", typeFn=int, default=None)
    if min_fragment:
        clip_cmd += ['-m', str(min_fragment)]
        # the min-fragment option is kind of ill thought out, and works as post-processing.  therefore it can
        # actually leave uncovered nodes, we clean here but note that it can leave fragments smaller than min-fragment
        clean_cmd = ['vg', 'clip', '-d', '1', '-']
        for ref in options.reference:
            clean_cmd += ['-P', ref]
        clip_cmd = [clip_cmd, clean_cmd]

    if getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "removeStubs", typeFn=bool, default=True):
        # this command can also leave fragments smaller than min-fragment
        stub_cmd = ['vg', 'clip', '-s', '-']
        for ref in options.reference:
            stub_cmd += ['-P', ref]
        if min_fragment:
            clip_cmd.append(stub_cmd)
        else:
            clip_cmd = [clip_cmd, stub_cmd]

    cactus_call(parameters=clip_cmd, outfile=clipped_path)

    # worth it
    cactus_call(parameters=['vg', 'validate', clipped_path])

    return job.fileStore.writeGlobalFile(clipped_path)
    
def join_vg(job, options, config, clipped_vg_ids):
    """ run vg ids -j
    """
    work_dir = job.fileStore.getLocalTempDir()
    vg_paths = []

    assert len(options.vg) == len(clipped_vg_ids)
    for vg_path, vg_id in zip(options.vg, clipped_vg_ids):
        vg_path = os.path.join(work_dir, os.path.basename(vg_path))
        job.fileStore.readGlobalFile(vg_id, vg_path, mutable=True)
        vg_paths.append(vg_path)
        
    cactus_call(parameters=['vg', 'ids', '-j'] + vg_paths)

    return [job.fileStore.writeGlobalFile(f) for f in vg_paths]

def vg_to_gfa(job, options, config, vg_path, vg_id):
    """ run gfa conversion """
    work_dir = job.fileStore.getLocalTempDir()
    vg_path = os.path.join(work_dir, os.path.basename(vg_path))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    out_path = vg_path + '.gfa'

    cmd = ['vg', 'convert', '-f', '-Q', options.reference[0], os.path.basename(vg_path), '-B']
    
    cactus_call(parameters=cmd, outfile=out_path, work_dir=work_dir)

    return job.fileStore.writeGlobalFile(out_path)

def make_vg_indexes(job, options, config, gfa_ids, tag=''):
    """ merge of the gfas, then make gbz / snarls / trans
    """ 
    work_dir = job.fileStore.getLocalTempDir()
    vg_paths = []
    merge_gfa_path = os.path.join(work_dir, '{}merged.gfa'.format(tag))

    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    
    # merge the gfas
    assert len(options.vg) == len(gfa_ids)
    for i, (vg_path, gfa_id) in enumerate(zip(options.vg, gfa_ids)):
        gfa_path = os.path.join(work_dir, os.path.basename(vg_path) +  '.gfa')
        job.fileStore.readGlobalFile(gfa_id, gfa_path, mutable=True)
        cmd = ['grep', '-v', '{}^W     {}'.format('^H\|' if i else '', graph_event), gfa_path]
        # add in the additional references here
        if i == 0 and len(options.reference) > 1:
            # if so, will need a new tool (or perhaps interface on vg paths?)
            cmd = [cmd, ['sed', '-e', '1s/{}/{}/'.format(options.reference[0], ' '.join(options.reference)),
                         '-e', '1s/{}//'.format(graph_event)]]
        cactus_call(parameters=cmd, outfile=merge_gfa_path, outappend=True)
        job.fileStore.deleteGlobalFile(gfa_id)

    # make the gbz
    gbz_path = os.path.join(work_dir, '{}merged.gbz'.format(tag))
    cactus_call(parameters=['vg', 'gbwt', '-G', merge_gfa_path, '--gbz-format', '-g', gbz_path])

    # zip the gfa
    cactus_call(parameters=['bgzip', merge_gfa_path, '--threads', str(job.cores)])
    gfa_path = merge_gfa_path + '.gz'

    # make the snarls
    snarls_path = os.path.join(work_dir, '{}merged.snarls'.format(tag))
    cactus_call(parameters=['vg', 'snarls', gbz_path, '-T', '-t', str(job.cores)], outfile=snarls_path)
                            
    return { '{}gfa.gz'.format(tag) : job.fileStore.writeGlobalFile(gfa_path),
             '{}gbz'.format(tag) : job.fileStore.writeGlobalFile(gbz_path),
             '{}snarls'.format(tag) : job.fileStore.writeGlobalFile(snarls_path) }

def make_vcf(job, config, out_name, vcf_ref, index_dict, tag='', ref_tag='', max_ref_allele=None):
    """ make the vcf
    """ 
    work_dir = job.fileStore.getLocalTempDir()
    gbz_path = os.path.join(work_dir, ref_tag + os.path.basename(out_name) + '.gbz')
    snarls_path = os.path.join(work_dir, ref_tag + os.path.basename(out_name) + '.snarls')
    job.fileStore.readGlobalFile(index_dict['{}gbz'.format(ref_tag)], gbz_path)
    job.fileStore.readGlobalFile(index_dict['{}snarls'.format(ref_tag)], snarls_path)

    # make the vcf
    vcf_path = os.path.join(work_dir, '{}merged.vcf.gz'.format(tag))
    decon_cmd = ['vg', 'deconstruct', gbz_path, '-P', vcf_ref, '-a', '-r', snarls_path, '-t', str(job.cores)]
    if getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "GFANodeIDsInVCF", typeFn=bool, default=True):
        decon_cmd.append('-O')
    cactus_call(parameters=[decon_cmd, ['bgzip', '--threads', str(job.cores)]], outfile=vcf_path)
    try:
        cactus_call(parameters=['tabix', '-p', 'vcf', vcf_path])
        tbi_path = vcf_path + '.tbi'
    except Exception as e:
        # todo: better support for larger chromosomes
        RealtimeLogger.warning('WARNING: tabix failed on VCF with this error {}'.format(str(e)))
        tbi_path = None

    output_dict = { '{}raw.vcf.gz'.format(tag) : job.fileStore.writeGlobalFile(vcf_path) }
    if tbi_path:
        output_dict['{}raw.vcf.gz.tbi'.format(tag)] = job.fileStore.writeGlobalFile(tbi_path)

    # make the filtered vcf
    if max_ref_allele:
        vcfbub_path = os.path.join(work_dir, 'merged.filtered.vcf.gz')
        cactus_call(parameters=[['vcfbub', '--input', vcf_path, '--max-ref-length', str(max_ref_allele), '--max-level', '0'],
                                ['bgzip', '--threads', str(job.cores)]],
                outfile=vcfbub_path)
        try:
            cactus_call(parameters=['tabix', '-p', 'vcf', vcfbub_path])
            tbi_path = vcfbub_path + '.tbi'
        except Exception as e:
            # todo: better support for larger chromosomes
            RealtimeLogger.warning('WARNING: tabix failed on VCF with this error {}'.format(str(e)))
            tbi_path = None            
        output_dict['{}vcf.gz'.format(tag)] = job.fileStore.writeGlobalFile(vcfbub_path)
        if tbi_path:
            output_dict['{}vcf.gz.tbi'.format(tag)] = job.fileStore.writeGlobalFile(vcfbub_path + '.tbi')

    return output_dict
    
def make_giraffe_indexes(job, options, config, index_dict, tag=''):
    """ make giraffe-specific indexes: distance and minimaer """
    work_dir = job.fileStore.getLocalTempDir()
    gbz_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.gbz')
    job.fileStore.readGlobalFile(index_dict['{}gbz'.format(tag)], gbz_path)

    # make the distance index
    dist_path = os.path.join(work_dir, tag + os.path.basename(options.outName) + '.dist')
    cactus_call(parameters=['vg', 'index', '-t', str(job.cores), '-j', dist_path, gbz_path])

    # make the minimizer index
    min_path = os.path.join(work_dir, os.path.basename(options.outName) + '.min')
    min_opts = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap_join"), "minimizerOptions", default='').split()
    cactus_call(parameters=['vg', 'minimizer'] + min_opts + ['-t', str(job.cores), '-d', dist_path, '-o', min_path, gbz_path])
                            
    return { '{}min'.format(tag) : job.fileStore.writeGlobalFile(min_path),
             '{}dist'.format(tag) : job.fileStore.writeGlobalFile(dist_path) }

def merge_hal(job, options, hal_ids):
    """ call halMergeChroms to make one big hal file out of the chromosome hal files """
    work_dir = job.fileStore.getLocalTempDir()
    hal_paths = []
    assert len(options.hal) == len(hal_ids)
    for in_path, hal_id in zip(options.hal, hal_ids):
        hal_path = os.path.join(work_dir, os.path.basename(in_path))
        job.fileStore.readGlobalFile(hal_id, hal_path)
        hal_paths.append(hal_path)

    merged_path = os.path.join(work_dir, '__merged__.hal')
    assert merged_path not in hal_paths

    # note: cactus_call tries to sort out relative paths by itself for docker.  but the comma-separated list
    # will likely throw it off, so we take care to specify it relative manually.
    # also note: most hal commands need --inMemory to run at scale, but the access patterns for chrom
    # merging are linear enough that it shouldn't be needed
    cmd = ['halMergeChroms',
           ','.join([os.path.basename(p) for p in hal_paths]),
           os.path.basename(merged_path),
           '--progress']
    cactus_call(parameters=cmd, work_dir = work_dir)

    return { 'full.hal' : job.fileStore.writeGlobalFile(merged_path) }

def unclip_hal(job, hal_id_dict, seq_id_map):
    """ run halUnclip """
    work_dir = job.fileStore.getLocalTempDir()
    in_hal_path = os.path.join(work_dir, "in.hal")
    job.fileStore.readGlobalFile(hal_id_dict['hal'], in_hal_path)
    out_hal_path = os.path.join(work_dir, "out.hal")
    seqfile_path = os.path.join(work_dir, "seqfile.txt")
    with open(seqfile_path, "w") as seqfile:
        for event, name_id in seq_id_map.items():
            seqfile.write("{}\t{}\n".format(event, job.fileStore.readGlobalFile(name_id[1])))

    cactus_call(parameters=["halUnclip", in_hal_path, seqfile_path, out_hal_path, "--progress", "--validate"])

    return { 'hal' : job.fileStore.writeGlobalFile(out_hal_path) }

def cat_bed_files(job, bed_ids):
    bed_paths = [job.fileStore.readGlobalFile(bed_id) for bed_id in bed_ids]
    cat_bed_path = job.fileStore.getLocalTempFile()
    catFiles(bed_paths, cat_bed_path)
    renamed_bed_path = cat_bed_path + '.renamed'
    # cactus-preprocess will make a bed like "ID=SAMPLE.1|CONTIG"
    # but by the name its been through cactus/hal it'll be SAMPLE.1.CONTIG
    # so we need to make the same transformation in the bed
    # if the id=| prefix is not present, it's up to the user to make sure it's consistent on input
    with open(cat_bed_path, 'r') as ifile, open(renamed_bed_path, 'w') as ofile:
        for line in ifile:
            if line.startswith('ID='):
                p = line.find('|')
                if p > 3:
                    sample = line[3:p]
                    rest = line[p+1:]
                    line = '{}.{}'.format(sample, rest)
            ofile.write(line)    
    return job.fileStore.writeGlobalFile(renamed_bed_path)

def cat_stats(job, stats_dict_list, zip_stats=True):
    work_dir = job.fileStore.getLocalTempDir()
    merged_dict = {}
    for key in stats_dict_list[0].keys():
        merged_dict[key] = os.path.join(work_dir, key)
        catFiles([job.fileStore.readGlobalFile(sd[key]) for sd in stats_dict_list], merged_dict[key])
    if zip_stats:
        cactus_call(parameters=['tar', 'czf', os.path.join(work_dir, 'stats.tgz')] + list(merged_dict.values()))
        return { 'stats.tgz': job.fileStore.writeGlobalFile(os.path.join(work_dir, 'stats.tgz')) }
    else:
        for key in stats_dict_list[0].keys():
            merged_dict[key] = job.fileStore.writeGlobalFile(merged_dict[key])
        return merged_dict

def export_join_data(toil, options, full_ids, clip_ids, clip_stats, filter_ids, idx_maps):
    """ download all the output data
    """

    # make a directory for the chromosomal vgs
    if options.chrom_vg:
        clip_base = os.path.join(options.outDir, '{}.chroms'.format(options.outName))
        if not clip_base.startswith('s3://') and not os.path.isdir(clip_base):
            os.makedirs(clip_base)

        if 'full' in options.chrom_vg:
            # download the "full" vgs
            assert len(options.vg) == len(full_ids)
            for vg_path, full_id,  in zip(options.vg, full_ids):
                name = os.path.splitext(vg_path)[0] + '.full.vg'
                toil.exportFile(full_id, makeURL(os.path.join(clip_base, os.path.basename(name))))

        if 'clip' in options.chrom_vg:
            # download the "clip" vgs
            assert len(options.vg) == len(clip_ids)
            for vg_path, clip_id,  in zip(options.vg, clip_ids):
                name = os.path.splitext(vg_path)[0] + '.vg'
                toil.exportFile(clip_id, makeURL(os.path.join(clip_base, os.path.basename(name))))

        if 'filter' in options.chrom_vg:
            # download the "filter" vgs
            assert len(options.vg) == len(filter_ids)
            for vg_path, filter_id,  in zip(options.vg, filter_ids):
                name = os.path.splitext(vg_path)[0] + '.d{}.vg'.format(options.filter)
                toil.exportFile(filter_id, makeURL(os.path.join(clip_base, os.path.basename(name))))
                
    # download the stats files
    if clip_stats:
        for stats_file in clip_stats.keys():
            toil.exportFile(clip_stats[stats_file], makeURL(os.path.join(options.outDir, '{}.{}'.format(options.outName, stats_file))))
        
    # download everything else
    for idx_map in idx_maps:
        for ext, idx_id in idx_map.items():
            # hacky filtering of intermediate indexes that the user doesn't want
            # ex if someone did --vcf clip --gfa full this would filter out clip.gfa etc.
            is_intermediate = False
            for out_ext, out_sel in [('gfa', options.gfa), ('gbz', options.gbz + options.giraffe), ('snarls', options.gbz)]:
                for phase in ['full', 'clip', 'filter']:
                    if '{}.{}'.format(phase, out_ext) in ext and phase not in out_sel:
                        is_intermediate = True
            if not is_intermediate:
                out_ext = ext
                if 'clip.' in out_ext:
                    out_ext = out_ext.replace('clip.', '')
                if 'filter.' in out_ext:
                    out_ext = out_ext.replace('filter.', 'd{}.'.format(options.filter))
                toil.exportFile(idx_id, makeURL(os.path.join(options.outDir, '{}.{}'.format(options.outName, out_ext))))
        
if __name__ == "__main__":
    main()
