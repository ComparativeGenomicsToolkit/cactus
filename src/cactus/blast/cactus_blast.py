#!/usr/bin/env python3

#Released under the MIT license, see LICENSE.txt

"""Run the pairwise blast stage only
   
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
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.progressive.schedule import Schedule
from cactus.progressive.projectWrapper import ProjectWrapper
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.pipeline.cactus_workflow import CactusWorkflowArguments
from cactus.pipeline.cactus_workflow import addCactusWorkflowOptions
from cactus.pipeline.cactus_workflow import CactusTrimmingBlastPhase
from cactus.shared.common import makeURL

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
    parser.add_argument("outputFile", type=str, help = "Output pairwise alignment file")

    #Progressive Cactus Options
    parser.add_argument("--database", dest="database",
                      help="Database type: tokyo_cabinet or kyoto_tycoon"
                      " [default: %(default)s]",
                      default="kyoto_tycoon")
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
    runCactusBlastOnly(options)
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    logger.info("cactus-blast has finished after {} seconds".format(run_time))

def runCactusBlastOnly(options):
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
            logger.info("Experiment {}".format(ET.tostring(expXml)))
            experiment = ExperimentWrapper(expXml)
            configPath = experiment.getConfigPath()
            configXml = ET.parse(configPath).getroot()
            
            seqIDMap = dict()
            tree = MultiCactusTree(experiment.getTree()).extractSubTree(options.root)
            leaves = tree.getChildNames(tree.getRootName())
            outgroups = experiment.getOutgroupGenomes()
            genome_set = set(leaves + outgroups)
            logger.info("Genomes in blastonly, {}: {}".format(options.root, list(genome_set)))

            #import the sequences (that we need to align for the given event, ie leaves and outgroups)
            for genome, seq in list(project.inputSequenceMap.items()):
                if genome in genome_set:
                    if os.path.isdir(seq):
                        tmpSeq = getTempFile()
                        catFiles([os.path.join(seq, subSeq) for subSeq in os.listdir(seq)], tmpSeq)
                        seq = tmpSeq
                    seq = makeURL(seq)
                    project.inputSequenceIDMap[genome] = toil.importFile(seq)
                else:
                    # out-of-scope sequences will only cause trouble later on
                    del project.inputSequenceMap[genome]

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

            outWorkFlowArgs = toil.start(CactusTrimmingBlastPhase(standAlone=True, cactusWorkflowArguments=workFlowArgs, phaseName="trimBlast"))

        # export the alignments
        toil.exportFile(outWorkFlowArgs.alignmentsID, makeURL(options.outputFile))
        outWorkFlowArgs.experimentWrapper.writeXML(options.outputFile + '.exp.xml')
        # optional secondary alignments
        if outWorkFlowArgs.secondaryAlignmentsID:
            toil.exportFile(outWorkFlowArgs.secondaryAlignmentsID, makeURL(options.outputFile) + '.secondary')
        # outgroup fragments and coverage are necessary for cactus-align, as the sequence names got changed in the above alignemnts
        for i, outgroupFragmentID in enumerate(outWorkFlowArgs.outgroupFragmentIDs):
            toil.exportFile(outgroupFragmentID, makeURL(options.outputFile) + '.og_fragment_{}'.format(i))
        # cactus-align can recompute coverage on the fly, but we save them because we have them 
        for i, ingroupCoverageID in enumerate(outWorkFlowArgs.ingroupCoverageIDs):
            toil.exportFile(ingroupCoverageID, makeURL(options.outputFile) + '.ig_coverage_{}'.format(i))
        
        
if __name__ == '__main__':
    main()
