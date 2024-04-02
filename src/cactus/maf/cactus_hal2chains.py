#!/usr/bin/env python3

"""
This is a script to convert a HAL file into a set of pairwise alignment in UCSC chains format.
"""

import os, sys
from argparse import ArgumentParser
import copy
import timeit, time
import math
import xml.etree.ElementTree as ET

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options
from cactus.shared.common import cactus_call
from cactus.shared.common import getOptionalAttrib, findRequiredNode
from cactus.shared.version import cactus_commit
from cactus.progressive.multiCactusTree import MultiCactusTree
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from cactus.shared.common import cactus_clamp_memory
from toil.lib.humanize import bytes2human
from sonLib.bioio import getTempDirectory
from sonLib.nxnewick import NXNewick

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("halFile", help = "HAL file to convert to MAF")
    parser.add_argument("outDir", help = "Output directory")
    parser.add_argument("--targetGenomes", nargs='*',
                        help="name(s) of target genomes (all leaves if empty).",
                        default=[])
    parser.add_argument("--queryGenomes", nargs='*',
                        help="name(s) of query genomes (all leaves if empty).",
                        default=[])
    parser.add_argument("--includeSelfAlignments", action="store_true",
                        help="include genome-vs-self chains (which are trivial one-block alignments)",
                        default=False)
    parser.add_argument("--bigChain", action="store_true",
                        help="output bigChain.bb and bigChain.link.bb files as well",
                        default=False)
    parser.add_argument("--useHalSynteny",
                        help="use halSynteny instead of halLiftover. halLiftover is default because that's what CAT uses",
                        action="store_true")
    parser.add_argument("--maxAnchorDistance", type=int,
                        help="set --maxAnchorDistance in halSynteny")
    parser.add_argument("--minBlockSize", type=int,
                        help="set --minBlockSize in halSynteny")

    parser.add_argument("--inMemory",
                       help="use --inMemory for halLiftover and/or halSynteny",
                       action="store_true")
    
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

    if options.maxAnchorDistance and not options.useHalSynteny:
        raise RuntimeError('--maxAnchorDistance can only be used with --useHalSynteny')
    if options.minBlockSize and not options.useHalSynteny:
        raise RuntimeError('--minBlockSize can only be used with --useHalSynteny')
    if options.includeSelfAlignments and options.useHalSynteny:
        raise RuntimeError('--includeSelfAlignments cannot be used with --useHalSynteny')
    
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    # make sure we have an outdir
    if not options.outDir.startswith('s3://') and not os.path.isdir(options.outDir):
        os.makedirs(options.outDir)

    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            chains_id_dict = toil.restart()
        else:
            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals(options)
            
            logger.info("Importing {}".format(options.halFile))
            hal_id = toil.importFile(options.halFile)            
            chains_id_dict = toil.start(Job.wrapJobFn(hal2chains_workflow, config, options, hal_id))

        #export the chains
        for query_genome in chains_id_dict.keys():
            clean_query_name = query_genome.replace('#', '.').replace(' ', '.')
            for target_genome in chains_id_dict[query_genome].keys():                
                clean_target_name = target_genome.replace('#', '.').replace(' ', '.')
                chain_id = chains_id_dict[query_genome][target_genome]['chains']
                out_base = makeURL(os.path.join(options.outDir, clean_target_name + '_vs_' + clean_query_name))
                toil.exportFile(chain_id, out_base + ".chain.gz")
                if options.bigChain:
                    bigchain_id = chains_id_dict[query_genome][target_genome]['bigChain']
                    toil.exportFile(bigchain_id, out_base + ".bigChain.bb")
                    biglink_id = chains_id_dict[query_genome][target_genome]['bigLink']
                    toil.exportFile(biglink_id, out_base + ".bigChain.link.bb")

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-hal2chains has finished after {} seconds".format(run_time))


def hal2chains_workflow(job, config, options, hal_id):
    root_job = Job()
    job.addChild(root_job)
    check_tools_job = root_job.addChildJobFn(hal2chains_check_tools, options)
    get_genomes_job = check_tools_job.addFollowOnJobFn(hal2chains_get_genomes, config, options, hal_id)
    leaf_genomes = get_genomes_job.rv()
    chrom_info_job = get_genomes_job.addFollowOnJobFn(hal2chains_chrom_info_all, config, options, hal_id, leaf_genomes)
    hal2chains_all_job = chrom_info_job.addFollowOnJobFn(hal2chains_all, config, options, hal_id, chrom_info_job.rv())
    return hal2chains_all_job.rv()

def hal2chains_check_tools(job, options):
    """ make sure we have the required ucsc commands available on the PATH"""
    tools = ['axtChain', 'faToTwoBit', 'pslPosTarget']
    if options.bigChain:
        tools += ['hgLoadChain', 'bedToBigBed']
    for tool in tools:
        try:
            cactus_call(parameters=[tool])
        except Exception as e:
            # hacky way to see if the usage info came out of running the command with no arguments
            if "usage:" in str(e):
                continue
        raise RuntimeError('Required tool, {}, not found in PATH. Please download it with wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/{}'.format(tool, tool))

def hal2chains_get_genomes(job, config, options, hal_id):    
    """ get the genomes in the hal file"""
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))    
    job.fileStore.readGlobalFile(hal_id, hal_path)
    tree_str = cactus_call(parameters=['halStats', hal_path, '--tree'], check_output=True).strip()
    mc_tree = MultiCactusTree(NXNewick().parseString(tree_str, addImpliedRoots=False))
    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    leaf_genomes = []
    anc_genomes = []
    for node in mc_tree.preOrderTraversal():
        if node == graph_event:
            continue
        if mc_tree.isLeaf(node):
            leaf_genomes.append(mc_tree.getName(node))
        else:
            anc_genomes.append(mc_tree.getName(node))

    if options.targetGenomes:
        for input_genome in options.targetGenomes:
            if input_genome not in leaf_genomes and input_genome not in anc_genomes:
                raise RuntimeError('Input --targetGenomes {} not found in HAL file'.format(input_genome))

    if options.queryGenomes:
        for input_genome in options.queryGenomes:
            if input_genome not in leaf_genomes and input_genome not in anc_genomes:
                raise RuntimeError('Input --queryGenomes {} not found in HAL file'.format(input_genome))
            
    return leaf_genomes
    
def hal2chains_chrom_info(job, config, options, hal_id, genome):
    """ get the BED and 2bit file for the genome from the hal """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))    
    job.fileStore.readGlobalFile(hal_id, hal_path)
    RealtimeLogger.info("Computing chromosomes and 2bit ref sequence")

    bed_path = os.path.join(work_dir, genome + '.bed')
    cactus_call(parameters=[['halStats', hal_path, '--chromSizes', genome],
                            ['awk', '{{print $1 "\t0\t" $2}}']],
                outfile=bed_path)
    if options.bigChain:
        sizes_path = os.path.join(work_dir, genome + '.chrom.sizes')
        cactus_call(parameters=['halStats', hal_path, '--chromSizes', genome], outfile=sizes_path)

    tbit_path = os.path.join(work_dir, genome + '.2bit')    
    cactus_call(parameters=[['hal2fasta', hal_path, genome],
                            ['faToTwoBit', 'stdin', tbit_path]])

    out_dict = { 'bed' : job.fileStore.writeGlobalFile(bed_path),
                 '2bit' : job.fileStore.writeGlobalFile(tbit_path) }
    if options.bigChain:
        out_dict['sizes'] = job.fileStore.writeGlobalFile(sizes_path)
    return out_dict

def hal2chains_chrom_info_all(job, config, options, hal_id, genomes):
    """ get the 2bit and chrom sizes of every relevant genome"""

    # TODO: this is gonna fall down on huge hal files because it runs a separate job
    # for each pair, and each job is going to require copying the full hal.  only work-around
    # would be to do something like cacts-hal2maf where jobs are done in multithread batches...

    # default query / target genome sets to all leaves if they weren't input
    if not options.queryGenomes:
        options.queryGenomes = genomes
    if not options.targetGenomes:
        options.targetGenomes = genomes

    chrom_info_dict = {}
    for genome in set(options.queryGenomes + options.targetGenomes):
        chrom_info_dict[genome] = job.addChildJobFn(hal2chains_chrom_info, config, options, hal_id, genome).rv()
    return chrom_info_dict        

def hal2chains_all(job, config, options, hal_id, chrom_info_dict):    
    """ convert all the genomes in parallel """

    # default query / target genome sets to all leaves if they weren't input
    if not options.queryGenomes:
        options.queryGenomes = list(chrom_info_dict.keys())
    if not options.targetGenomes:
        options.targetGenomes = list(chrom_info_dict.keys())

    # TODO: this is gonna fall down on huge hal files because it runs a separate job
    # for each pair, and each job is going to require copying the full hal.  only work-around
    # would be to do something like cacts-hal2maf where jobs are done in multithread batches...
        
    output_dict = {}
    hal_mem = hal_id.size if options.inMemory else int(hal_id.size / 10)
    for query_genome in options.queryGenomes:
        output_dict[query_genome] = {}
        query_info = chrom_info_dict[query_genome]
        for target_genome in options.targetGenomes:
            target_info = chrom_info_dict[target_genome]
            if options.includeSelfAlignments or target_genome != query_genome:
                chains_job = job.addChildJobFn(hal2chains_genome, config, options, hal_id,
                                               query_genome, query_info, target_genome,
                                               disk=int(hal_id.size * 1.2) + query_info['2bit'].size * 10 + target_info['2bit'].size * 10,
                                               memory=cactus_clamp_memory(query_info['2bit'].size * 10 + target_info['2bit'].size * 10 + hal_mem))
                output_dict[query_genome][target_genome] = {}
                output_dict[query_genome][target_genome]['chains'] = chains_job.rv()
                if options.bigChain:
                    bigchains_job = chains_job.addFollowOnJobFn(chain2bigchain, options, query_genome,
                                                                target_genome, chrom_info_dict[target_genome], chains_job.rv(),
                                                                disk = query_info['2bit'].size * 20 + target_info['2bit'].size * 20,
                                                                memory = cactus_clamp_memory(query_info['2bit'].size * 10 + target_info['2bit'].size * 10))
                    output_dict[query_genome][target_genome]['bigChain'] = bigchains_job.rv(0)
                    output_dict[query_genome][target_genome]['bigLink'] = bigchains_job.rv(1)

    return output_dict

    
def hal2chains_genome(job, config, options, hal_id, query_genome, query_info, target_genome):
    """ convert one genome to chains """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))    
    job.fileStore.readGlobalFile(hal_id, hal_path)
    RealtimeLogger.info("Computing chains for query {} to target {}".format(query_genome, target_genome))

    bed_path = os.path.join(work_dir, query_genome + '.bed')
    job.fileStore.readGlobalFile(query_info['bed'], bed_path)
    query_2bit_path = os.path.join(work_dir, query_genome + '.2bit')
    job.fileStore.readGlobalFile(query_info['2bit'], query_2bit_path)

    # we need the target sequence for axtChain
    target_2bit_path = os.path.join(work_dir, target_genome + '.2bit')    
    cactus_call(parameters=[['hal2fasta', hal_path, target_genome],
                            ['faToTwoBit', 'stdin', target_2bit_path]])

    # the output
    query_chain_path = os.path.join(work_dir, target_genome + '_vs_' + query_genome + '.chain.gz')

    # make our hal -> psl -> chain piped command

    cmd = []

    if options.useHalSynteny:
        hal_synteny_cmd = ['halSynteny', hal_path, '/dev/stdout', '--queryGenome', query_genome, '--targetGenome', target_genome]
        if options.inMemory:
            hal_synteny_cmd.append('--inMemory')
        if options.maxAnchorDistance:
            hal_synteny_cmd += ['--maxAnchorDistance', str(options.maxAnchorDistance)]
        if options.minBlockSize:
            hal_synteny_cmd += ['--minBlockSize', str(options.minBlockSize)]
        cmd.append(hal_synteny_cmd)
    else:
        hal_liftover_cmd = ['halLiftover', hal_path, query_genome, bed_path, target_genome, '/dev/stdout', '--outPSL']
        if options.inMemory:
            hal_liftover_cmd.append('--inMemory')
        cmd.append(hal_liftover_cmd)

    cmd.append(['pslPosTarget', '/dev/stdin', '/dev/stdout'])

    cmd.append(['axtChain', '-psl', '-verbose=0', '-linearGap=medium', '/dev/stdin', target_2bit_path, query_2bit_path, '/dev/stdout'])

    cmd.append(['gzip'])

    cactus_call(parameters=cmd, outfile = query_chain_path)

    return job.fileStore.writeGlobalFile(query_chain_path)

# https://genome.ucsc.edu/goldenPath/help/examples/bigChain.as
bigChain_as = \
'''table bigChain
"bigChain pairwise alignment"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name or ID of item, ideally both human readable and unique"
    uint score;         "Score (0-1000)"
    char[1] strand;     "+ or - for strand"
    uint tSize;         "size of target sequence"
    string qName;       "name of query sequence"
    uint qSize;         "size of query sequence"
    uint qStart;        "start of alignment on query sequence"
    uint qEnd;          "end of alignment on query sequence"
    uint chainScore;    "score from chain"
    )
'''

# https://genome.ucsc.edu/goldenPath/help/examples/bigLink.as
bigLink_as = \
'''
table bigLink
"bigLink pairwise alignment"
    (
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name or ID of item, ideally both human readable and unique"
    uint qStart;        "start of alignment on query sequence"
    )
'''

def chain2bigchain(job, options, query_genome, target_genome, target_info, chain_id):
    """ convert the chain to big chain. from https://genome.ucsc.edu/goldenPath/help/bigChain.html """
    work_dir = job.fileStore.getLocalTempDir()
    chain_path = os.path.join(work_dir, target_genome + '_vs_' + query_genome + '.chain.gz')
    job.fileStore.readGlobalFile(chain_id, chain_path)
    target_sizes_path = os.path.join(work_dir, target_genome + '.chrom.sizes')
    job.fileStore.readGlobalFile(target_info['sizes'], target_sizes_path)

    # unzip the chains
    uz_chain_path = chain_path[:-3]
    cactus_call(parameters=['gzip', '-dc', chain_path], outfile=uz_chain_path)

    # make the bigChain.as file
    bigchain_as_path = os.path.join(work_dir, 'bigChain.as')
    with open(bigchain_as_path, 'w') as as_file:
        as_file.write(bigChain_as)

    # Use the hgLoadChain utility to generate the chain.tab and link.tab files needed to create the bigChain file:
    cactus_call(parameters=['hgLoadChain', '-noBin', '-test', target_genome, 'bigChain', uz_chain_path])

    # Create the bigChain file from your input chain file using a combination of sed, awk and the bedToBigBed utility:
    bigchain_path = os.path.join(work_dir, target_genome + '_vs_' + query_genome + '.bigChain')
    cactus_call(parameters=[['sed', 's/\\.000000//', 'chain.tab'],
                            ['awk', 'BEGIN {OFS=\"\\t\"} {print $2, $4, $5, $11, 1000, $8, $3, $6, $7, $9, $10, $1}']],
                outfile=bigchain_path)

    bigchain_bb_path = bigchain_path + '.bb'
    cactus_call(parameters=['bedToBigBed', '-type=bed6+6', '-as={}'.format(bigchain_as_path),
                            '-tab', bigchain_path, target_sizes_path, bigchain_bb_path])

    # make the bigLink.as file
    biglink_as_path = os.path.join(work_dir, 'bigLink.as')
    with open(biglink_as_path, 'w') as as_file:
        as_file.write(bigLink_as)

    # To display your date in the Genome Browser, you must also create a binary indexed link file to accompany your bigChain file:
    biglink_path = bigchain_path + '.link'
    cactus_call(parameters=[['awk', 'BEGIN {OFS=\"\\t\"} {print $1, $2, $3, $5, $4}', 'link.tab'],
                            ['sort', '-k1,1', '-k2,2n']],
                outfile=biglink_path)

    biglink_bb_path = biglink_path + '.bb'
    cactus_call(parameters=['bedToBigBed', '-type=bed4+1', 'as={}'.format(biglink_as_path),
                            '-tab', biglink_path, target_sizes_path, biglink_bb_path])

    return (job.fileStore.writeGlobalFile(bigchain_bb_path), job.fileStore.writeGlobalFile(biglink_bb_path))
    
         


