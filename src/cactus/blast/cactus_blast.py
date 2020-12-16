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
from cactus.shared.common import makeURL, catFiles
from cactus.shared.common import enableDumpStack
from cactus.shared.common import cactus_override_toil_options

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
    parser.add_argument("--pathOverrides", nargs="*", help="paths (multiple allowed) to override from seqFile")
    parser.add_argument("--pathOverrideNames", nargs="*", help="names (must be same number as --pathOverrides) of path overrides")

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
    options.database = "kyoto_tycoon"

    setupBinaries(options)
    setLoggingFromOptions(options)
    enableDumpStack()

    if (options.pathOverrides or options.pathOverrideNames):
        if not options.pathOverrides or not options.pathOverrideNames or \
           len(options.pathOverrideNames) != len(options.pathOverrides):
            raise RuntimeError('same number of values must be passed to --pathOverrides and --pathOverrideNames')

    # Mess with some toil options to create useful defaults.
    cactus_override_toil_options(options)

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

            # apply path overrides.  this was necessary for wdl which doesn't take kindly to
            # text files of local paths (ie seqfile).  one way to fix would be to add support
            # for s3 paths and force wdl to use it.  a better way would be a more fundamental
            # interface shift away from files of paths throughout all of cactus
            if options.pathOverrides:
                seqFile = SeqFile(options.seqFile)
                configNode = ET.parse(options.configFile).getroot()
                config = ConfigWrapper(configNode)
                tree = MultiCactusTree(seqFile.tree)
                tree.nameUnlabeledInternalNodes(prefix = config.getDefaultInternalNodePrefix())
                for name, override in zip(options.pathOverrideNames, options.pathOverrides):
                    seqFile.pathMap[name] = override
                override_seq = os.path.join(options.cactusDir, 'seqFile.override')
                with open(override_seq, 'w') as out_sf:
                    out_sf.write(str(seqFile))
                options.seqFile = override_seq

            #to be consistent with all-in-one cactus, we make sure the project
            #isn't limiting itself to the subtree (todo: parameterize so root can
            #be passed through from prepare to blast/align)
            proj_options = copy.deepcopy(options)
            proj_options.root = None
            #Create the progressive cactus project (as we do in runCactusProgressive)
            projWrapper = ProjectWrapper(proj_options, proj_options.configFile, ignoreSeqPaths=options.root)
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

            print (str(project.inputSequenceMap))


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
