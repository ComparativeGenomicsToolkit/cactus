#!/usr/bin/env python3

"""
This is a script to convert a MAF (from cactus-hal2maf) to BigMaf
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
from cactus.maf.cactus_hal2maf import get_sed_rename_scripts
from toil.lib.humanize import bytes2human
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)

    parser.add_argument("mafFile", help = "MAF file to convert to BigMaf (can be gzipped)")
    parser.add_argument("outFile", help = "Output bigMaf file (.bb)")
    parser.add_argument("--refGenome", help = "reference genome to get chrom sizes from", required=True)

    
    parser.add_argument("--halFile", help = "HAL file to use to get chrom sizes from. Will also be used to work around dots in genome names")
    parser.add_argument("--chromSizes", help = "File of chromosome sizes (can be obtained with halStats --chromSizes)")
    
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

    if bool(options.halFile) == bool(options.chromSizes):
        raise RuntimeError('Either --chromSizes or --halFile must be used to get the chromosome sizes, not both')
    if not options.outFile.endswith('.bb'):
        raise RuntimeError('output file path must end with .bb')
    
    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

    logger.info('Cactus Command: {}'.format(' '.join(sys.argv)))
    logger.info('Cactus Commit: {}'.format(cactus_commit))
    start_time = timeit.default_timer()

    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            bigmaf_id_dict = toil.restart()
        else:
            #load cactus config
            configNode = ET.parse(options.configFile).getroot()
            config = ConfigWrapper(configNode)
            config.substituteAllPredefinedConstantsWithLiterals(options)
            
            logger.info("Importing {}".format(options.mafFile))
            maf_id = toil.importFile(options.mafFile)
            chrom_sizes_id = None
            if options.chromSizes:
                logger.info("Importing {}".format(options.chromSizes))
                chrom_sizes_id = toil.importFile(options.chromSizes)
            hal_id = None
            if options.halFile:
                logger.info("Importing {}".format(options.halFile))
                hal_id = toil.importFile(options.halFile)
                
            bigmaf_id_dict = toil.start(Job.wrapJobFn(maf2bigmaf_workflow, config, options, maf_id, chrom_sizes_id, hal_id))

        #export the big maf
        out_bm_path = makeURL(options.outFile)
        logger.info("Exporting {}".format(out_bm_path))
        toil.exportFile(bigmaf_id_dict['bb'], out_bm_path)
        #export the big maf summary
        out_bms_path = makeURL(os.path.splitext(options.outFile)[0] + '.summary.bb')
        logger.info("Exporting {}".format(out_bms_path))
        toil.exportFile(bigmaf_id_dict['summary.bb'], out_bms_path)
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-maf2bigmaf has finished after {} seconds".format(run_time))


def maf2bigmaf_workflow(job, config, options, maf_id, chrom_sizes_id, hal_id):
    root_job = Job()
    job.addChild(root_job)
    check_tools_job = root_job.addChildJobFn(maf2bigmaf_check_tools)
    genomes_list = None
    if not chrom_sizes_id:
        chrom_sizes_job = check_tools_job.addFollowOnJobFn(maf2bigmaf_chrom_sizes, options, hal_id,
                                                           disk=hal_id.size)
        prev_job = chrom_sizes_job
        chrom_sizes_id = chrom_sizes_job.rv(0)
        genomes_list = chrom_sizes_job.rv(1)
    else:
        prev_job = check_tools_job
    if options.mafFile.endswith('.gz'):
        disk_mult = 7
    else:
        disk_mult = 3
    bigmaf_job = prev_job.addFollowOnJobFn(maf2bigmaf, maf_id, chrom_sizes_id, genomes_list, options,
                                           disk=disk_mult * maf_id.size,
                                           memory=min(2**32, max(maf_id.size, 2**36)))
    bigmaf_summary_job = prev_job.addFollowOnJobFn(maf2bigmaf_summary, maf_id, chrom_sizes_id, genomes_list, options,
                                                   disk=2 * maf_id.size,
                                                   memory=min(2**32, max(maf_id.size, 2**36)))

    return { 'bb' : bigmaf_job.rv(),  'summary.bb' : bigmaf_summary_job.rv() }

def maf2bigmaf_check_tools(job):
    """ make sure we have the required ucsc commands available on the PATH"""
    for tool in ['mafToBigMaf', 'bedToBigBed', 'hgLoadMafSummary']:
        try:
            cactus_call(parameters=[tool])
        except Exception as e:
            # hacky way to see if the usage info came out of running the command with no arguments
            if "usage:" in str(e):
                continue
        raise RuntimeError('Required tool, {}, not found in PATH. Please download it with wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/{}'.format(tool, tool))

def maf2bigmaf_chrom_sizes(job, options, hal_id):
    """ get chromsizes from hal, also get a unique token to map underscores and dots to"""
    assert options.refGenome
    
    work_dir = job.fileStore.getLocalTempDir()
    hal_path = os.path.join(work_dir, os.path.basename(options.halFile.replace(' ', '.')))
    RealtimeLogger.info("Reading HAL file from job store to {}".format(hal_path))    
    job.fileStore.readGlobalFile(hal_id, hal_path)
    RealtimeLogger.info("Computing chromosome sizes")

    chrom_sizes = os.path.join(work_dir, options.refGenome + '.chrom_sizes')
    cactus_call(parameters=['halStats', hal_path, '--chromSizes', options.refGenome],
                outfile = chrom_sizes)

    genomes = set(cactus_call(parameters=['halStats', '--genomes', hal_path], check_output=True).split())
    
    return job.fileStore.writeGlobalFile(chrom_sizes), genomes

# https://genome.ucsc.edu/goldenPath/help/examples/bigMaf.as
bigMaf_as = \
'''table bedMaf
"Bed3 with MAF block"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    lstring mafBlock;   "MAF block"
    )
'''

# https://genome.ucsc.edu/goldenPath/help/examples/mafSummary.as
mafSummary_as = \
'''table mafSummary
"Positions and scores for alignment blocks"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string src;        "Sequence name or database of alignment"
    float  score;      "Floating point score."
    char[1] leftStatus;  "Gap/break annotation for preceding block"
    char[1] rightStatus; "Gap/break annotation for following block"
    )
'''

def maf2bigmaf(job, maf_id, chrom_sizes_id, genomes_list, options):
    """ do the conversion from https://genome.ucsc.edu/goldenPath/help/bigMaf.html to
    make the BigMaf bigbed file. This takes a lot of disk since it requires making
    an uncompressed BED file intermediary"""
    work_dir = job.fileStore.getLocalTempDir()
    maf_path = os.path.join(work_dir, os.path.basename(options.mafFile.replace(' ', '.')))    
    RealtimeLogger.info("Reading MAF file from job store to {}".format(maf_path))
    job.fileStore.readGlobalFile(maf_id, maf_path, mutable=True)    
    chrom_sizes_path = os.path.join(work_dir, options.refGenome + '.chrom_sizes')    
    job.fileStore.readGlobalFile(chrom_sizes_id, chrom_sizes_path)
    bigmaf_path = os.path.join(work_dir, os.path.basename(options.outFile))
                               
    # make the bigMaf.as file
    bigmaf_as_path = os.path.join(work_dir, 'bigMaf.as')
    assert(maf_path != bigmaf_as_path)
    with open(bigmaf_as_path, 'w') as as_file:
        as_file.write(bigMaf_as)
        
    # make the sed scripts to handle non-alphanumeric characters
    sed_scripts = get_sed_rename_scripts(work_dir, genomes_list, out_bed=True)

    cat_cmd = ['gzip', '-dc'] if maf_path.endswith('.gz') else ['cat']

    # make the bigmaf bed file
    bigmaf_bed_path = os.path.join(work_dir, 'bigMaf.bed')
    bigmaf_cmd = [cat_cmd + [maf_path],
                  ['mafDuplicateFilter', '-km', '-'],
                  ['mafToBigMaf', options.refGenome, 'stdin', 'stdout'],
                  ['sort', '-k1,1', '-k2,2n']]
    cactus_call(parameters=bigmaf_cmd, outfile=bigmaf_bed_path)

    # convert it to bigbed (can't pipe!)
    cactus_call(parameters=['bedToBigBed', '-type=bed3+1', '-as={}'.format(bigmaf_as_path),
                            '-tab', bigmaf_bed_path, chrom_sizes_path, bigmaf_path])
    
    return job.fileStore.writeGlobalFile(bigmaf_path)

def maf2bigmaf_summary(job, maf_id, chrom_sizes_id, genomes_list, options):
    """ do the conversion from https://genome.ucsc.edu/goldenPath/help/bigMaf.html to 
    make the BigMafSummary BigBed """
    work_dir = job.fileStore.getLocalTempDir()
    maf_path = os.path.join(work_dir, os.path.basename(options.mafFile.replace(' ', '.')))    
    RealtimeLogger.info("Reading MAF file from job store to {}".format(maf_path))
    job.fileStore.readGlobalFile(maf_id, maf_path)    
    chrom_sizes_path = os.path.join(work_dir, options.refGenome + '.chrom_sizes')    
    job.fileStore.readGlobalFile(chrom_sizes_id, chrom_sizes_path)
    bigmaf_path = os.path.join(work_dir, os.path.basename(options.outFile))
    summary_path = os.path.splitext(bigmaf_path)[0] + '.summary.bb'
                                       
    # make the mafSummary.as file
    mafsummary_as_path = os.path.join(work_dir, 'mafSummary.as')
    assert(maf_path != mafsummary_as_path)
    with open(mafsummary_as_path, 'w') as as_file:
        as_file.write(mafSummary_as)

    # make the sed scripts to handle non-alphanumeric characters
    sed_scripts = get_sed_rename_scripts(work_dir, genomes_list, out_bed=True)

    cat_cmd = ['gzip', '-dc'] if maf_path.endswith('.gz') else ['cat']

    # make the bigmaf summary file
    ref_name = options.refGenome
    summary_bed_path = os.path.join(work_dir, 'summary.bed')
    summary_cmd = [cat_cmd + [maf_path]]
    if sed_scripts:
        summary_cmd += [['sed', '-f', sed_scripts[0]]]
        if options.refGenome in sed_scripts[2]:
            ref_name = sed_scripts[2][options.refGenome]
    summary_cmd += [['mafDuplicateFilter', '-km', '-'],
                    ['hgLoadMafSummary', '-minSeqSize=1', '-test', ref_name, 'bigMafSummary', 'stdin']]
    cactus_call(parameters=summary_cmd)
    sort_cmd = [['cut', '-f2-', 'bigMafSummary.tab']]                
    if sed_scripts:
        sort_cmd += [['sed', '-f', sed_scripts[1]]]
    sort_cmd += [['sort', '-k1,1', '-k2,2n']]
    cactus_call(parameters=sort_cmd, outfile=summary_bed_path)    

    # convert it to bigbed (can't pipe!)
    cactus_call(parameters=['bedToBigBed', '-type=bed3+4', '-as={}'.format(mafsummary_as_path),
                            '-tab', summary_bed_path, chrom_sizes_path, summary_path])

    return job.fileStore.writeGlobalFile(summary_path)
