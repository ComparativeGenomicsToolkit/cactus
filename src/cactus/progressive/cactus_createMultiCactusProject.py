#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

"""Create the multi_cactus xml and directory structure from a workflow template
"""
import os
import xml.etree.ElementTree as ET
import copy

from operator import itemgetter

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.outgroup import GreedyOutgroup, DynamicOutgroup
from sonLib.nxnewick import NXNewick

def createMCProject(tree, experiment, config, options):
    """
    Creates a properly initialized MultiCactusProject.

    TODO: This should really all be in the constructor for MultiCactusProject.
    """
    mcTree = MultiCactusTree(tree)
    mcTree.nameUnlabeledInternalNodes(config.getDefaultInternalNodePrefix())
    mcTree.computeSubtreeRoots()
    mcProj = MultiCactusProject()
    for genome in experiment.getGenomesWithSequence():
        mcProj.inputSequenceMap[genome] = experiment.getSequenceID(genome)
    mcProj.mcTree = mcTree
    if config.getDoSelfAlignment():
        mcTree.addSelfEdges()
    for name in mcProj.mcTree.getSubtreeRootNames():
        expPath = "%s/%s/%s_experiment.xml" % (options.path, name, name)
        mcProj.expMap[name] = os.path.abspath(expPath)
    alignmentRootId = mcProj.mcTree.getRootId()
    if options.root is not None:
        try:
            alignmentRootId = mcProj.mcTree.getNodeId(options.root)
        except:
            raise RuntimeError("Specified root name %s not found in tree" % options.root)

    fillInOutgroups(mcProj, options.outgroupNames, config, alignmentRootId)

    # if necessary, we reroot the tree at the specified alignment root id.  all leaf genomes
    # that are no longer in the tree, but still used as outgroups, are moved into special fields
    # so that we can remember to, say, get their paths for preprocessing.
    specifyAlignmentRoot(mcProj, experiment, alignmentRootId)
    return mcProj

def fillInOutgroups(mcProj, outgroupNames, config, alignmentRootId):
    """
    Determines the outgroups for a MultiCactusProject using the strategy from the config.
    """
    mcProj.outgroup = None
    if config.getOutgroupStrategy() == 'greedy':
        # use the provided outgroup candidates, or use all outgroups
        # as candidates if none are given
        mcProj.outgroup = GreedyOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree, alignmentRootId)
        mcProj.outgroup.greedy(threshold=config.getOutgroupThreshold(),
                               candidateSet=outgroupNames,
                               candidateChildFrac=config.getOutgroupAncestorQualityFraction(),
                               maxNumOutgroups=config.getMaxNumOutgroups())
    elif config.getOutgroupStrategy() == 'greedyLeaves':
        # use all leaves as outgroups, unless outgroup candidates are given
        mcProj.outgroup = GreedyOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree, alignmentRootId)
        mcProj.outgroup.greedy(threshold=config.getOutgroupThreshold(),
                               candidateSet=outgroupNames,
                               candidateChildFrac=2.0,
                               maxNumOutgroups=config.getMaxNumOutgroups())
    elif config.getOutgroupStrategy() == 'greedyPreference':
        # prefer the provided outgroup candidates, if any, but use
        # other nodes as "filler" if we can't find enough.
        mcProj.outgroup = GreedyOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree, alignmentRootId)
        mcProj.outgroup.greedy(threshold=config.getOutgroupThreshold(),
                               candidateSet=outgroupNames,
                               candidateChildFrac=config.getOutgroupAncestorQualityFraction(),
                               maxNumOutgroups=config.getMaxNumOutgroups())
        mcProj.outgroup.greedy(threshold=config.getOutgroupThreshold(),
                               candidateSet=None,
                               candidateChildFrac=config.getOutgroupAncestorQualityFraction(),
                               maxNumOutgroups=config.getMaxNumOutgroups())
    elif config.getOutgroupStrategy() == 'dynamic':
        # dynamic programming algorithm that exactly optimizes probability
        # that base in target node aligns to at least one base in the
        # outgroup set.  Caveats are that it only returns leaves, and
        # the model used for optimization is super naive. Still, it does
        # some things better than greedy approaches such as properly account
        # for phylogenetic redundancy, as well as try to factor assembly
        # size/quality automatically.
        mcProj.outgroup = DynamicOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree, mcProj.inputSequenceMap, alignmentRootId,
                                   candidateSet=outgroupNames)
        mcProj.outgroup.compute(maxNumOutgroups=config.getMaxNumOutgroups())
    elif config.getOutgroupStrategy() != 'none':
        raise RuntimeError("Could not understand outgroup strategy %s" % config.getOutgroupStrategy())

# it is possible that we start with a much bigger tree than we actually want to align
# (controlled by the --root option in cactus_createMultiCactusProject.py).  We use
# the entire tree when selecting outgroups, but right afterward have no use for
# genomes that are neither outgroups, nor in the alignment.  We especially don't
# want to waste time preprocessing them.  This function, reroots the tree at the
# alignment root then tacks on all the outgroups from ooutside the new tree
# (which must be leaves) as children of the root
def specifyAlignmentRoot(mcProj, expTemplate, alignmentRootId):
    # ugly hack to keep track of external outgroups for root experiment (yuck)
    mcProj.externalOutgroupNames = set()

    # keep around the entire un-rerooted tree so that we can calculate
    # the induced species trees for each node correctly -- gross!
    mcProj.entireTree = copy.deepcopy(mcProj.mcTree)

    # nothing to do
    if alignmentRootId == mcProj.mcTree.getRootId():
        return

    # dig out every outgroup
    outGroupNames = set()
    eventsInSubtree = set(mcProj.mcTree.getName(i) for i in mcProj.mcTree.postOrderTraversal(alignmentRootId))
    if mcProj.outgroup is not None:
        for event, ogNameDistList in list(mcProj.outgroup.ogMap.items()):
            if event in eventsInSubtree:
                for og, dist in ogNameDistList:
                    outGroupNames.add(og)

    # find outgroups we want to extract
    keptNodes = set(mcProj.mcTree.postOrderTraversal(alignmentRootId))
    deadNodes = []
    for node in mcProj.mcTree.postOrderTraversal():
        if node not in keptNodes:
            deadNodes.append(node)
            name = mcProj.mcTree.getName(node)
            if name in outGroupNames:
                mcProj.externalOutgroupNames.add(name)

    # reroot the tree!
    mcProj.mcTree.reroot(alignmentRootId)

    # add the outgroups to the tree (and sequence map)
    # computing distance to new root for each one
    for ogName in mcProj.externalOutgroupNames:
        ogId = mcProj.mcTree.getNodeId(ogName)
        dist = 0.
        x = ogId
        while mcProj.mcTree.hasParent(x):
            d = mcProj.mcTree.getWeight(mcProj.mcTree.getParent(x), x)
            if d is None:
                dist = None
                break
            else:
                dist += d
            x = mcProj.mcTree.getParent(x)
        mcProj.mcTree.addOutgroup(ogName, dist)

    # remove any experiment directories that have become invalid
    for event in list(mcProj.expMap.keys()):
        if mcProj.mcTree.getNodeId(event) in deadNodes:
            del mcProj.expMap[event]

    # flush out all unused nodes, set the new root, and update the
    # experiment template tree to match the new project tree
    for node in deadNodes:
        assert mcProj.mcTree.hasParent(node)
        mcProj.mcTree.removeEdge(mcProj.mcTree.getParent(node), node)
    expTemplate.setTree(mcProj.mcTree)

    # reset input sequences to only contain genomes in tree
    genomesInTree = [mcProj.mcTree.getName(node) for node in mcProj.mcTree.postOrderTraversal() if mcProj.mcTree.hasName(node)]
    genomesNotInTree = []
    for genome in mcProj.inputSequenceMap:
        if genome not in genomesInTree:
            genomesNotInTree.append(genome)
    for genome in genomesNotInTree:
        del mcProj.inputSequenceMap[genome]

# Make the subdirs for each subproblem:  name/ and name/name_DB
# and write the experiment files
# and copy over a config with updated reference field
def createFileStructure(mcProj, expTemplate, configTemplate, options):
    if not os.path.exists(options.path):
        os.makedirs(options.path)
    mcProj.writeXML(os.path.join(options.path, "%s_project.xml" % options.name))

    for name, expPath in list(mcProj.expMap.items()):
        path = os.path.join(options.path, name)
        children = mcProj.entireTree.getChildNames(name)

        # Get outgroups
        outgroups = []
        if configTemplate.getOutgroupStrategy() != 'none' \
        and name in mcProj.outgroup.ogMap:
            # Outgroup name is the first element of the ogMap tuples
            outgroups.extend(list(map(itemgetter(0), mcProj.outgroup.ogMap[name])))

        subtree = mcProj.entireTree.extractSpanningTree(children + [name] + outgroups)
        exp = ExperimentWrapper.createExperimentWrapper(NXNewick().writeString(subtree),
                                                        children + [name] + outgroups,
                                                        databaseConf=expTemplate.confElem)

        exp.setRootGenome(name)
        exp.setOutgroupGenomes(outgroups)

        if not os.path.exists(path):
            os.makedirs(path)
        config = ConfigWrapper(copy.deepcopy(configTemplate.xmlRoot))
        if expTemplate.getSequenceID(name):
            exp.setRootReconstructed(False)
            exp.setSequenceID(name, expTemplate.getSequenceID(name))
        else:
            exp.setRootReconstructed(True)
        exp.writeXML(expPath)

def checkInputSequencePaths(exp):
    for seq in exp.getSequences():
        if not os.path.exists(seq):
            raise RuntimeError("sequence path %s does not exist\n" % seq)
        elif os.path.isdir(seq):
            contents = os.listdir(seq)
            size = 0
            for i in contents:
                if i[0] != '.':
                    size += 1
            if size == 0:
                raise RuntimeError("Sequence path %s is an empty directory\n" % seq)

class CreateMultiCactusProjectOptions:
    def __init__(self, expFile, projectFile,
                 outgroupNames, root, overwrite):
        self.expFile = expFile
        self.path = projectFile
        self.name = os.path.basename(self.path)

        self.outgroupNames = outgroupNames
        self.root = root
        self.overwrite = overwrite

def runCreateMultiCactusProject(expFile, projectFile,
            outgroupNames=None, root=None, overwrite=False):
    options = CreateMultiCactusProjectOptions(expFile, projectFile,
            outgroupNames=outgroupNames, root=root, overwrite=overwrite)

    expTemplate = ExperimentWrapper(ET.parse(options.expFile).getroot())
    configPath = expTemplate.getConfigPath()
    confTemplate = ConfigWrapper(ET.parse(configPath).getroot())
    tree = expTemplate.getTree()
    if options.outgroupNames is not None:
        options.outgroupNames = set(options.outgroupNames)
        projNames = set([tree.getName(x) for x in tree.getLeaves()])
        for outgroupName in options.outgroupNames:
            if outgroupName not in projNames:
                raise RuntimeError("Specified outgroup %s not found in tree" % outgroupName)

    mcProj = createMCProject(tree, expTemplate, confTemplate, options)

    expTemplate.setTree(mcProj.mcTree)

    # Now do the file tree creation
    createFileStructure(mcProj, expTemplate, confTemplate, options)
