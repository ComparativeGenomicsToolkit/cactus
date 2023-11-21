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
from toil.job import Job
from toil.common import Toil
from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger
from cactus.shared.common import cactus_cpu_count
from toil.lib.humanize import bytes2human
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("halFile", help = "HAL file to convert to MAF")
    parser.add_argument("outDir", help = "Output directory")

    parser.add_argument("--refGenome", required=True,
                        help="name of reference genome (root if empty)",
                        default=None)
    parser.add_argument("--targetGenomes",
                        help="comma-separated (no spaces) list of target "
                        "genomes (others are excluded) (vist all non-ancestral if empty)",
                        default=None)
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
        for tgt_genome, chain_id in chains_id_dict.items():
            clean_name = tgt_genome.replace('#', '.').replace(' ', '.')
            toil.exportFile(chain_id, makeURL(os.path.join(options.outDir, clean_name + ".chain.gz")))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-hal2chains has finished after {} seconds".format(run_time))


def hal2chains_workflow(job, config, options, hal_id):
    root_job = Job()
    job.addChild(root_job)
    check_tools_job = root_job.addChildJobFn(hal2chains_check_tools)
    ref_info_job = check_tools_job.addFollowOnJobFn(hal2chains_ref_info, config, options, hal_id)
    hal2chains_all_job = ref_info_job.addFollowOnJobFn(hal2chains_all, config, options, hal_id, ref_info_job.rv())
    return hal2chains_all_job.rv()

def hal2chains_check_tools(job):
    """ make sure we have the required ucsc commands available on the PATH"""
    for tool in ['axtChain', 'faToTwoBit', 'pslPosTarget']:
        try:
            cactus_call(parameters=[tool])
        except Exception as e:
            # hacky way to see if the usage info came out of running the command with no arguments
            if "usage:" in str(e):
                continue
        raise RuntimeError('Required tool, {}, not found in PATH. Please download it with wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/{}'.format(tool, tool))
            
def hal2chains_ref_info(job, config, options, hal_id):
    """ get the BED and 2bit file for the reference genome from the hal """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))    
    job.fileStore.readGlobalFile(hal_id, hal_path)
    RealtimeLogger.info("Computing chromosomes and 2bit ref sequence")

    bed_path = os.path.join(work_dir, options.refGenome + '.bed')
    cactus_call(parameters=[['halStats', hal_path, '--chromSizes', options.refGenome],
                            ['awk', '{{print $1 "\t0\t" $2}}']],
                outfile=bed_path)

    tbit_path = os.path.join(work_dir, options.refGenome + '.2bit')    
    cactus_call(parameters=[['hal2fasta', hal_path, options.refGenome],
                            ['faToTwoBit', 'stdin', tbit_path]])

    genomes=cactus_call(parameters=['halStats', '--genomes', hal_path], check_output=True).split()

    graph_event = getOptionalAttrib(findRequiredNode(config.xmlRoot, "graphmap"), "assemblyName", default="_MINIGRAPH_")
    leaf_genomes = []
    for genome in genomes:
        if genome != options.refGenome and genome != graph_event and \
           not cactus_call(parameters=['halStats', hal_path, '--children', genome], check_output=True).strip():
            leaf_genomes.append(genome)
    
    return { 'bed' : job.fileStore.writeGlobalFile(bed_path),
             '2bit' : job.fileStore.writeGlobalFile(tbit_path),
             'genomes' : genomes,
             'leaf_genomes' : leaf_genomes }

def hal2chains_all(job, config, options, hal_id, ref_info):
    """ convert all the genomes in parallel """
    tgt_genomes = []
    hal_genomes = set(ref_info['genomes'])
    if options.targetGenomes:
        for opt_target_genome in options.targetGenomes.split(','):
            if opt_target_genome not in hal_genomes:
                raise RuntimeError('--targetGenome {} not found in HAL'.format(opt_target_genome))
            else:
                tgt_genomes.append(opt_target_genome)
    else:
        tgt_genomes = ref_info['leaf_genomes']
        
    chains_output_dict = {}
    for tgt_genome in tgt_genomes:
        chains_job = job.addChildJobFn(hal2chains_genome, config, options, hal_id, ref_info, tgt_genome,
                                       disk=int(hal_id.size * 1.2),
                                       memory=(ref_info['2bit'].size * 10))
        chains_output_dict[tgt_genome] = chains_job.rv()
    return chains_output_dict
    
def hal2chains_genome(job, config, options, hal_id, ref_info, tgt_genome):
    """ convert one genome to chains """
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))    
    job.fileStore.readGlobalFile(hal_id, hal_path)
    RealtimeLogger.info("Computing chains for {}".format(tgt_genome))

    bed_path = os.path.join(work_dir, options.refGenome + '.bed')
    job.fileStore.readGlobalFile(ref_info['bed'], bed_path)
    ref_2bit_path = os.path.join(work_dir, options.refGenome + '.2bit')
    job.fileStore.readGlobalFile(ref_info['2bit'], ref_2bit_path)

    # we need the target sequence for axtChain
    tgt_2bit_path = os.path.join(work_dir, tgt_genome + '.2bit')    
    cactus_call(parameters=[['hal2fasta', hal_path, tgt_genome],
                            ['faToTwoBit', 'stdin', tgt_2bit_path]])

    # the output
    tgt_chain_path = os.path.join(work_dir, tgt_genome + '.chain.gz')

    # make our hal -> psl -> chain piped command

    cmd = []

    if options.useHalSynteny:
        hal_synteny_cmd = ['halSynteny', hal_path, '/dev/stdout', '--queryGenome', options.refGenome, '--targetGenome', tgt_genome]
        if options.inMemory:
            hal_synteny_cmd.append('--inMemory')
        if options.maxAnchorDistance:
            hal_synteny_cmd += ['--maxAnchorDistance', str(options.maxAnchorDistance)]
        if options.minBlockSize:
            hal_synteny_cmd += ['--minBlockSize', str(options.minBlockSize)]
        cmd.append(hal_synteny_cmd)
    else:
        hal_liftover_cmd = ['halLiftover', hal_path, options.refGenome, bed_path, tgt_genome, '/dev/stdout', '--outPSL']
        if options.inMemory:
            hal_liftover_cmd.append('--inMemory')
        cmd.append(hal_liftover_cmd)

    cmd.append(['pslPosTarget', '/dev/stdin', '/dev/stdout'])

    cmd.append(['axtChain', '-psl', '-verbose=0', '-linearGap=medium', '/dev/stdin', tgt_2bit_path, ref_2bit_path, '/dev/stdout'])

    cmd.append(['gzip'])

    cactus_call(parameters=cmd, outfile = tgt_chain_path)

    return job.fileStore.writeGlobalFile(tgt_chain_path)
            
