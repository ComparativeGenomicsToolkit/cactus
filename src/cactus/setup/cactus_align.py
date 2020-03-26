#!/usr/bin/env python3

#Released under the MIT license, see LICENSE.txt

"""Run the multiple alignment on pairwise alignment input (ie cactus_setup_phase and beyond)
   
"""
import os
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
import timeit

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import setupBinaries, importSingularityImage
from cactus.progressive.cactus_progressive import exportHal
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.progressive.schedule import Schedule
from cactus.progressive.projectWrapper import ProjectWrapper
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import CactusTrimmingBlastPhase
from cactus.pipeline.cactus_workflow import CactusSetupCheckpoint
from cactus.pipeline.cactus_workflow import prependUniqueIDs
from cactus.blast.blast import calculateCoverage
from cactus.shared.common import makeURL
from toil.realtimeLogger import RealtimeLogger

from toil.job import Job
from toil.common import Toil
from toil.lib.bioio import logger
from toil.lib.bioio import setLoggingFromOptions

from sonLib.nxnewick import NXNewick
from sonLib.bioio import getTempDirectory

def main():
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    addCactusWorkflowOptions(parser)

    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("blastOutput", type=str, help = "Blast output (from cactus-blast)")
    parser.add_argument("outputHal", type=str, help = "Output HAL file")

    #Progressive Cactus Options
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--root", dest="root", help="Name of ancestral node (which"
                        " must appear in NEWICK tree in <seqfile>) to use as a "
                        "root for the alignment.  Any genomes not below this node "
                        "in the tree may be used as outgroups but will never appear"
                        " in the output.  If no root is specifed then the root"
                        " of the tree is used. ", default=None, required=True)
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
    setLoggingFromOptions(options)

    options.database = 'kyoto_tycoon'

    options.buildHal = True
    options.buildFasta = True
    
    # Mess with some toil options to create useful defaults.

    # Caching generally slows down the cactus workflow, plus some
    # methods like readGlobalFileStream don't support forced
    # reads directly from the job store rather than from cache.
    options.disableCaching = True
    # Job chaining breaks service termination timing, causing unused
    # databases to accumulate and waste memory for no reason.
    options.disableChaining = True
    if options.retryCount is None:
        # If the user didn't specify a retryCount value, make it 5
        # instead of Toil's default (1).
        options.retryCount = 5

    start_time = timeit.default_timer()
    runCactusAfterBlastOnly(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-blast has finished after {} seconds".format(run_time))

def runCactusAfterBlastOnly(options):
    with Toil(options) as toil:
        importSingularityImage(options)
        #Run the workflow
        if options.restart:
            alignmentID = toil.restart()
        else:
            options.cactusDir = getTempDirectory()
            
            #Create the progressive cactus project (as we do in runCactusProgressive)
            projWrapper = ProjectWrapper(options, options.configFile, ignoreSeqPaths=options.root)
            projWrapper.writeXml()

            pjPath = os.path.join(options.cactusDir, ProjectWrapper.alignmentDirName,
                                  '%s_project.xml' % ProjectWrapper.alignmentDirName)
            assert os.path.exists(pjPath)

            project = MultiCactusProject()

            if not os.path.isdir(options.cactusDir):
                os.makedirs(options.cactusDir)

            project.readXML(pjPath)

            # open up the experiment (as we do in ProgressiveUp.run)
            # note that we copy the path into the options here
            experimentFile = project.expMap[options.root]
            expXml = ET.parse(experimentFile).getroot()
            experiment = ExperimentWrapper(expXml)
            configPath = experiment.getConfigPath()
            configXml = ET.parse(configPath).getroot()
            
            seqIDMap = dict()
            tree = MultiCactusTree(experiment.getTree()).extractSubTree(options.root)
            leaves = [tree.getName(leaf) for leaf in tree.getLeaves()]
            outgroups = experiment.getOutgroupGenomes()
            genome_set = set(leaves + outgroups)

            # import the outgroups
            outgroupIDs = []
            cactus_blast_input = True
            for i, outgroup in enumerate(outgroups):
                try:
                    outgroupID = toil.importFile(makeURL(options.blastOutput) + '.og_fragment_{}'.format(i))
                    outgroupIDs.append(outgroupID)
                    experiment.setSequenceID(outgroup, outgroupID)
                except:
                    # we assume that input is not coming from cactus blast, so we'll treat output
                    # sequences normally and not go looking for fragments
                    outgroupIDs = []
                    cactus_blast_input = False
                    break

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            for genome, seq in list(project.inputSequenceMap.items()):
                if genome in leaves or (not cactus_blast_input and genome in outgroups):
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    
                    experiment.setSequenceID(genome, toil.importFile(seq))
                    
            if not cactus_blast_input:
                outgroupIDs = [experiment.getSequenceID(outgroup) for outgroup in outgroups]
    
            # write back the experiment, as CactusWorkflowArguments wants a path
            experiment.writeXML(experimentFile)

            #import cactus config
            if options.configFile:
                cactusConfigID = toil.importFile(makeURL(options.configFile))
            else:
                cactusConfigID = toil.importFile(makeURL(project.getConfigPath()))
            project.setConfigID(cactusConfigID)

            project.syncToFileStore(toil)
            configNode = ET.parse(project.getConfigPath()).getroot()
            configWrapper = ConfigWrapper(configNode)
            configWrapper.substituteAllPredefinedConstantsWithLiterals()            
            
            workFlowArgs = CactusWorkflowArguments(options, experimentFile=experimentFile, configNode=configNode, seqIDMap = project.inputSequenceIDMap)

            #import the files that cactus-blast made
            workFlowArgs.alignmentsID = toil.importFile(makeURL(options.blastOutput))
            try:
                workFlowArgs.secondaryAlignmentsID = toil.importFile(makeURL(options.blastOutput) + '.secondary')
            except:
                workFlowArgs.secondaryAlignmentsID = None
            workFlowArgs.outgroupFragmentIDs = outgroupIDs
            workFlowArgs.ingroupCoverageIDs = []
            if cactus_blast_input and len(outgroups) > 0:
                for i in range(len(leaves)):
                    workFlowArgs.ingroupCoverageIDs.append(toil.importFile(makeURL(options.blastOutput) + '.ig_coverage_{}'.format(i)))

            halID = toil.start(Job.wrapJobFn(run_cactus_align, configWrapper, workFlowArgs, project, cactus_blast_input))

        # export the hal
        toil.exportFile(halID, makeURL(options.outputHal))
        
def run_cactus_align(job, configWrapper, cactusWorkflowArguments, project, cactus_blast_input):

    if cactus_blast_input:
        # cactus expects unique ids, but they weren't saved from cactus-blast
        first_job = job.addChildJobFn(run_prepend_unique_ids, cactusWorkflowArguments, project,
                                    #todo disk=
        )
    else:
        # we only have cigar input, so compute the ingroup coverage bed files now
        first_job = job.addChildJobFn(run_ingroup_coverage, cactusWorkflowArguments, project)
    cactusWorkflowArguments = first_job.rv()        
    
    # run cactus setup all the way through cactus2hal generation
    setup_job = first_job.addFollowOnJobFn(run_setup_phase, cactusWorkflowArguments)

    # set up the project
    prepare_hal_export_job = setup_job.addFollowOnJobFn(run_prepare_hal_export, project, setup_job.rv())

    # create the hal
    hal_export_job = prepare_hal_export_job.addFollowOnJobFn(exportHal, prepare_hal_export_job.rv(0), event=prepare_hal_export_job.rv(1),
                                                             memory=configWrapper.getDefaultMemory(),
                                                             disk=configWrapper.getExportHalDisk(),
                                                             preemptable=False)
    return hal_export_job.rv()

def run_prepend_unique_ids(job, cactusWorkflowArguments, project):
    """ prepend the unique ids on the input fasta.  this is required for cactus to work (would be great to relax it though"""

    # note, there is an order dependence to everything where we have to match what was done in cactus_workflow
    # (so the code is pasted exactly as it is there)
    # this is horrible and needs to be fixed via drastic interface refactor
    exp = cactusWorkflowArguments.experimentWrapper
    ingroupsAndOriginalIDs = [(g, exp.getSequenceID(g)) for g in exp.getGenomesWithSequence() if g not in exp.getOutgroupGenomes()]
    sequences = [job.fileStore.readGlobalFile(id) for id in map(itemgetter(1), ingroupsAndOriginalIDs)]
    cactusWorkflowArguments.totalSequenceSize = sum(os.stat(x).st_size for x in sequences)
    renamedInputSeqDir = job.fileStore.getLocalTempDir()
    uniqueFas = prependUniqueIDs(sequences, renamedInputSeqDir)
    uniqueFaIDs = [job.fileStore.writeGlobalFile(seq, cleanup=True) for seq in uniqueFas]
    # Set the uniquified IDs for the ingroups and outgroups
    ingroupsAndNewIDs = list(zip(list(map(itemgetter(0), ingroupsAndOriginalIDs)), uniqueFaIDs[:len(ingroupsAndOriginalIDs)]))
    for event, sequenceID in ingroupsAndNewIDs:
        cactusWorkflowArguments.experimentWrapper.setSequenceID(event, sequenceID)
    return cactusWorkflowArguments

def run_ingroup_coverage(job, cactusWorkflowArguments, project):
    """ for every ingroup genome, make a bed file by computing its coverge vs the outgroups """
    work_dir=job.fileStore.getLocalTempDir()
    exp = cactusWorkflowArguments.experimentWrapper
    ingroupsAndOriginalIDs = [(g, exp.getSequenceID(g)) for g in exp.getGenomesWithSequence() if g not in exp.getOutgroupGenomes()]
    outgroups = [job.fileStore.readGlobalFile(id) for id in cactusWorkflowArguments.outgroupFragmentIDs]
    sequences = [job.fileStore.readGlobalFile(id) for id in map(itemgetter(1), ingroupsAndOriginalIDs)]
    cactusWorkflowArguments.totalSequenceSize = sum(os.stat(x).st_size for x in sequences)
    ingroups = map(itemgetter(0), ingroupsAndOriginalIDs)
    cigar = job.fileStore.readGlobalFile(cactusWorkflowArguments.alignmentsID)
    # should we parallelize with child jobs?
    for ingroup, sequence in zip(ingroups, sequences):
        coverage_path = os.path.join(work_dir, '{}.coverage'.format(sequence))
        calculateCoverage(sequence, cigar, coverage_path, fromGenome=outgroups, work_dir=work_dir)
        cactusWorkflowArguments.ingroupCoverageIDs.append(job.fileStore.writeGlobalFile(coverage_path))
    return cactusWorkflowArguments
    
def run_setup_phase(job, cactusWorkflowArguments):
    # needs to be its own job to resovolve the workflowargument promise
    return job.addChild(CactusSetupCheckpoint(cactusWorkflowArguments=cactusWorkflowArguments, phaseName="setup")).rv()

def run_prepare_hal_export(job, project, experiment):
    """ hack up the given project into something that gets exportHal() to do what we want """
    event = experiment.getRootGenome()
    exp_path = os.path.join(job.fileStore.getLocalTempDir(), event + '_experiment.xml')
    experiment.writeXML(exp_path)
    project.expMap = {event : experiment}
    project.expIDMap = {event : job.fileStore.writeGlobalFile(exp_path)}
    return project, event
    

if __name__ == '__main__':
    main()
