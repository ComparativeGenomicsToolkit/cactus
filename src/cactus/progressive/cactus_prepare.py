#!/usr/bin/env python3

#Released under the MIT license, see LICENSE.txt

"""Set up a given seqfile to include proprocessed/ancestral sequences, as well as to

"""
import os
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy
from sonLib.bioio import getTempDirectory

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import cactusRootPath
from cactus.shared.common import enableDumpStack
from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.schedule import Schedule
from cactus.progressive.projectWrapper import ProjectWrapper

def main():
    parser = ArgumentParser()
    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("outSeqDir", help='Directory where the processed leaf sequence and ancestral sequences will be placed')
    parser.add_argument("outSeqFile", help = "Path for annotated Seq file output")
    parser.add_argument("outputHal", type=str, help = "Output HAL file")
    parser.add_argument("--configFile", default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))
    parser.add_argument("--preprocessBatchSize", type=int, default=3, help="size (number of genomes) of suggested preprocessing jobs")
    parser.add_argument("--jobStore", type=str, default="$JOBSTORE", help="jobstore to use in suggested commands")
    parser.add_argument("--halOptions", type=str, default="--hdf5InMemory", help="options for every hal command")
    parser.add_argument("--cactusOptions", type=str, default="--realTimeLogging --logInfo", help="options for every cactus command")
    parser.add_argument("--preprocessOnly", action="store_true", help="only decompose into preprocessor and cactus jobs")

    options = parser.parse_args()
    options.database = 'kyoto_tycoon'
    #todo support root option
    options.root = None

    # need to go through this garbage (copied from the main() in progressive_cactus) to
    # come up with the project
    options.cactusDir = getTempDirectory()
    #Create the progressive cactus project
    projWrapper = ProjectWrapper(options, options.configFile)
    projWrapper.writeXml()

    pjPath = os.path.join(options.cactusDir, ProjectWrapper.alignmentDirName,
                          '%s_project.xml' % ProjectWrapper.alignmentDirName)
    assert os.path.exists(pjPath)

    project = MultiCactusProject()

    if not os.path.isdir(options.cactusDir):
        os.makedirs(options.cactusDir)

    project.readXML(pjPath)

    enableDumpStack()
    cactusPrepare(options, project)

def cactusPrepare(options, project):
    """ annotate a SeqFile with ancestral names as well as paths for output sequences."""

    # read the input
    seqFile = SeqFile(options.seqFile)
    configNode = ET.parse(options.configFile).getroot()
    config = ConfigWrapper(configNode)

    # prepare output sequence directory
    # todo: support remote (ie s3) output directory
    try:
        os.makedirs(options.outSeqDir)
    except:
        pass
    if not os.path.isdir(options.outSeqDir):
        raise RuntimeError('Unable to create output sequence directory \'{}\''.format(options.outSeqDir))
    if not os.access(options.outSeqDir, os.W_OK):
        logger.warning('Output sequence directory is not writeable: \'{}\''.format(options.outSeqDir))

    # hack the configfile to skip preprocessing and write it to the output dir
    if options.preprocessOnly:
        config.removePreprocessors()
        options.configFile = os.path.join(options.outSeqDir, 'config.xml')
        config.writeXML(options.configFile)
        
    # pass through the config file to the options
    # todo (don't like second hard-code check of .xml path)
    if options.configFile != os.path.join(cactusRootPath(), "cactus_progressive_config.xml"):
        options.cactusOptions += ' --configFile {}'.format(options.configFile)

    # get the ancestor names
    tree = MultiCactusTree(seqFile.tree)
    tree.nameUnlabeledInternalNodes(prefix = config.getDefaultInternalNodePrefix())

    # make the output
    outSeqFile = SeqFile()
    outSeqFile.tree= tree
    outSeqFile.pathMap = seqFile.pathMap
    outSeqFile.outgroups = seqFile.outgroups

    # update paths for preprocessed leaves or inferred ancestors
    for node in outSeqFile.tree.breadthFirstTraversal():
        name = outSeqFile.tree.getName(node)
        leaf = outSeqFile.tree.isLeaf(node)
        if leaf or (not leaf and name not in seqFile.pathMap and not options.preprocessOnly):
            out_basename = seqFile.pathMap[name] if name in seqFile.pathMap else '{}.fa'.format(name)
            outSeqFile.pathMap[name] = os.path.join(options.outSeqDir, os.path.basename(out_basename))

    # write the output
    with open(options.outSeqFile, 'w') as out_sf:
        out_sf.write(str(outSeqFile))

    # write the instructions
    print(get_plan(options, project, outSeqFile))

def get_plan(options, project, outSeqFile):

    plan = ''

    # preprocessing
    plan += '\n## Preprocessor\n'
    leaves = [outSeqFile.tree.getName(leaf) for leaf in outSeqFile.tree.getLeaves()]
    for i in range(0, len(leaves), options.preprocessBatchSize):
        pre_batch = leaves[i:i+options.preprocessBatchSize]
        plan += 'cactus-preprocess {} {} {} --inputNames {} {}\n'.format(
            options.jobStore, options.seqFile, options.outSeqFile, ' '.join(pre_batch),
            options.cactusOptions)

    if options.preprocessOnly:
        # if we're only making preprocess jobs, we can end early
        plan += '\n## Cactus\n'
        plan += 'cactus {} {} {} {}\n'.format(options.jobStore, options.outSeqFile,
                                              options.outputHal, options.cactusOptions)
        return plan

    # shedule up the alignments
    schedule = Schedule()
    schedule.loadProject(project)
    schedule.compute()

    # set of all jobs, as genome names from the (fully resolved, output) seqfile
    events = set(outSeqFile.pathMap.keys()) - set(leaves)
    resolved = set(leaves)

    # convert follow-ons to dependencies
    follow_on_deps = {}
    for event in events:
        fo = schedule.followOn(event)
        if fo:
            follow_on_deps[fo] = event

    def get_deps(event):
        deps = set(schedule.deps(event))
        if event in follow_on_deps:
            deps = deps.union(set(follow_on_deps[event]))
        # I don't know why the schedule doesn't always give the children
        # todo: understand!
        try:
            has_name = outSeqFile.tree.getNodeId(event) is not None
        except:
            has_name = False
        if has_name:
            for node in outSeqFile.tree.getChildren(outSeqFile.tree.getNodeId(event)):
                if not outSeqFile.tree.isLeaf(node):
                    deps.add(outSeqFile.tree.getName(node))
        return deps

    events_and_virtuals = set()
    for event in events:
        events_and_virtuals.add(event)
        if get_deps(event):
            events_and_virtuals = events_and_virtuals.union(get_deps(event))

    # group jobs into rounds.  where all jobs of round i can be run in parallel
    groups = []
    while len(events_and_virtuals) > 0:
        group = []
        to_remove = []
        added = 0
        for event in events_and_virtuals:
            if all([dep in resolved for dep in get_deps(event)]):
                if not schedule.isVirtual(event):
                    group.append(event)
                to_remove.append(event)
                added += 1
        if added == 0:
            sys.stderr.write("schedule deadlock:\n")
            for event in events_and_virtuals:
                sys.stderr.write("{} has deps {}\b".format(event, get_deps(event)))
            sys.exit(1)
        for tr in to_remove:
            resolved.add(tr)
            events_and_virtuals.remove(tr)
        groups.append(group)

    def halPath(event):
        if event == project.mcTree.getRootName():
            return options.outputHal
        else:
            return os.path.join(options.outSeqDir, event + '.hal')
    def cigarPath(event):
        return os.path.join(options.outSeqDir, event + '.cigar')

    # alignment groups
    plan += '\n## Alignment\n'
    for i, group in enumerate(groups):
        plan += '\n### Round {}'.format(i)
        for event in sorted(group):
            # todo: support cactus interface (it's easy enough here, but cactus_progressive.py needs changes to handle)
            plan += '\n'
            plan += 'cactus-blast {} {} {} --root {} {}\n'.format(
                options.jobStore, options.outSeqFile, cigarPath(event), event, options.cactusOptions)
            plan += 'cactus-align {} {} {} {} --root {} {}\n'.format(
                options.jobStore, options.outSeqFile, cigarPath(event), halPath(event), event, options.cactusOptions)
            # todo: just output the fasta in cactus-align.
            plan += 'hal2fasta {} {} {} > {}\n'.format(halPath(event), event, options.halOptions, outSeqFile.pathMap[event])

    # stitch together the final tree
    plan += '\n## HAL merging\n'
    root = project.mcTree.getRootName()
    for group in reversed(groups):
        for event in group:
            if event != root:
                plan += 'halAppendSubtree {} {} {} {} --merge {}\n'.format(
                    halPath(root), halPath(event), event, event, options.halOptions)

    return plan

if __name__ == '__main__':
    main()
