#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Create the multi_cactus xml and directory structure from a workflow template
""" 
import os
import sys
import xml.etree.ElementTree as ET
import copy

from optparse import OptionParser
from operator import itemgetter

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.outgroup import GreedyOutgroup, DynamicOutgroup
from sonLib.nxnewick import NXNewick
from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor

def createMCProject(tree, sequences, experiment, config, options):
    """
    Creates a properly initialized MultiCactusProject.

    TODO: This should really all be in the constructor for MultiCactusProject.
    """
    mcTree = MultiCactusTree(tree)
    mcTree.nameUnlabeledInternalNodes(config.getDefaultInternalNodePrefix())
    mcTree.computeSubtreeRoots()
    mcProj = MultiCactusProject()
    mcProj.mcTree = mcTree
    mcProj.inputSequences = sequences
    mcProj.outputSequenceDir = experiment.getOutputSequenceDir()
    if config.getDoSelfAlignment():
        mcTree.addSelfEdges()
    for name in mcProj.mcTree.getSubtreeRootNames():
        expPath = "%s/%s/%s_experiment.xml" % (options.path, name, name)
        mcProj.expMap[name] = os.path.abspath(expPath)
    alignmentRootId = mcProj.mcTree.getRootId()
    if options.root is not None:
        try:
            alignmentRootId = mcProj.mcTree.getNodeId(options.root)
        except Exception as e:
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
        mcProj.outgroup.importTree(mcProj.mcTree, mcProj.getInputSequenceMap(), alignmentRootId,
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
    if mcProj.outgroup is not None:
        for event, ogNameDistList in mcProj.outgroup.ogMap.items():
            for og, dist in ogNameDistList:
                outGroupNames.add(og)

    # find outgroups we want to extract
    keptNodes = set([x for x in mcProj.mcTree.postOrderTraversal(alignmentRootId)])
    deadNodes = []
    extractOutgroupMap = dict()
    for node in mcProj.mcTree.postOrderTraversal():
        if node not in keptNodes:
            deadNodes.append(node)
            name = mcProj.mcTree.getName(node)
            if name in outGroupNames:
                assert name in expTemplate.getGenomesWithSequence()
                extractOutgroupMap[name] = expTemplate.getSequencePath(name)
                mcProj.externalOutgroupNames.add(name)

    # reroot the tree!
    oldParent = mcProj.mcTree.getParent(alignmentRootId)
    mcProj.mcTree.reroot(alignmentRootId)

    # add the outgroups to the tree (and sequence map)
    # computing distance to new root for each one
    for ogName, ogPath in extractOutgroupMap.items():
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
    for event in mcProj.expMap.keys():
        if mcProj.mcTree.getNodeId(event) in deadNodes:
            del mcProj.expMap[event]
            
    # flush out all unused nodes, and set the new root
    for node in deadNodes:
        assert mcProj.mcTree.hasParent(node)
        mcProj.mcTree.removeEdge(mcProj.mcTree.getParent(node), node)

    # reset input sequences to only contain genomes in tree
    genomesInTree = [mcProj.mcTree.getName(node) for node in mcProj.mcTree.postOrderTraversal() if mcProj.mcTree.hasName(node)]
    seqsInTree = set(expTemplate.getSequencePath(genome) for genome in genomesInTree)
    mcProj.inputSequences = [seq for seq in mcProj.inputSequences if seq in seqsInTree]

# Make the subdirs for each subproblem:  name/ and name/name_DB
# and write the experiment files
# and copy over a config with updated reference field
def createFileStructure(mcProj, expTemplate, configTemplate, options):
    if not os.path.exists(options.path):
        os.makedirs(options.path)
    mcProj.writeXML(os.path.join(options.path, "%s_project.xml" % options.name))
    portOffset = 0

    for name, expPath in mcProj.expMap.items():
        path = os.path.join(options.path, name)
        children = mcProj.entireTree.getChildNames(name)

        # Get outgroups
        outgroups = []
        if configTemplate.getOutgroupStrategy() != 'none' \
        and name in mcProj.outgroup.ogMap:
            # Outgroup name is the first element of the ogMap tuples
            outgroups.extend(map(itemgetter(0), mcProj.outgroup.ogMap[name]))

        subtree = mcProj.entireTree.extractSpanningTree(children + outgroups)
        exp = ExperimentWrapper.createExperimentWrapper(NXNewick().writeString(subtree),
                                                        expTemplate.getOutputSequenceDir(),
                                                        databaseConf=expTemplate.confElem)

        exp.setRootGenome(name)

        # Get subtree connecting children + outgroups
        assert len(children) > 0

        for genome in children + outgroups:
            seqPath = expTemplate.getSequencePath(genome)
            if seqPath is None:
                # Ancestral sequence yet to be reconstructed. Point it
                # at where it will be once this problem runs.
                subPath = os.path.join(options.path, genome)
                seqPath = os.path.join(subPath, genome + '.fa')
            exp.setSequencePath(genome, seqPath)

        dbBase = path
        if expTemplate.getDbDir() is not None:
            dbBase = os.path.abspath(expTemplate.getDbDir())
        exp.setDbDir(os.path.join(dbBase, name, "%s_DB" % name))
        if expTemplate.getDbType() == "kyoto_tycoon" and \
            os.path.splitext(name)[1] != ".kch":
            exp.setDbName("%s.kch" % name)
        else:
            exp.setDbName(name)
        if expTemplate.getDbType() == "kyoto_tycoon":
            exp.setDbPort(expTemplate.getDbPort() + portOffset)
            portOffset += 1
            host = expTemplate.getDbHost()
            if host is not None:
                exp.setDbHost(host)
        exp.setReferencePath(os.path.join(path, name + '.fa'))
        if configTemplate.getBuildHal() == True:
            exp.setHALPath(os.path.join(path, "%s_hal.c2h" % name))
        if configTemplate.getBuildFasta() == True:
            exp.setHALFastaPath(os.path.join(path, "%s_hal.fa" % name))
        exp.setTree(subtree)
        exp.setOutgroupGenomes(outgroups)
        exp.setConfigPath(os.path.join(path, "%s_config.xml" % name))
        if not os.path.exists(exp.getDbDir()):
            os.makedirs(exp.getDbDir())
        if not os.path.exists(path):
            os.makedirs(path)
        config = ConfigWrapper(copy.deepcopy(configTemplate.xmlRoot))
        if expTemplate.getSequencePath(name):
            exp.setRootReconstructed(False)
        else:
            exp.setRootReconstructed(True)
        exp.writeXML(expPath)
        config.writeXML(exp.getConfigPath())

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

def main():
    usage = "usage: %prog [options] <experiment> <output project path>"
    description = "Setup a multi-cactus project using an experiment xml as template"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--fixNames",
                      help="This function is now a no-op")
    parser.add_option("--outgroupNames", dest="outgroupNames",  default = None, 
                      help="comma-separated names of high quality assemblies to use as outgroups [default=everything]")
    parser.add_option("--root", dest="root", type=str,
                      help="name of alignment root (must be labeled ancestral node in tree in input experiment).  Useful "
                      "for allowing the tree to contain nodes that won't be in the alignment but can still be used for "
                      "outgroups.",
                      default=None)
    parser.add_option("--overwrite", action="store_true", help="Overwrite existing experiment files", default=False)

    options, args = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    options.expFile = args[0]    
    options.path = os.path.abspath(args[1])
    options.name = os.path.basename(options.path)

    if (os.path.isdir(options.path) and not options.overwrite) or os.path.isfile(options.path):
        raise RuntimeError("Output project path %s exists\n" % options.path)

    expTemplate = ExperimentWrapper(ET.parse(options.expFile).getroot())
    configPath = expTemplate.getConfigPath()
    confTemplate = ConfigWrapper(ET.parse(configPath).getroot())
    checkInputSequencePaths(expTemplate)
    tree = expTemplate.getTree()

    # Check that the tree is sensible (root has at least 1 child)
    if len(tree.getChildren(tree.getRootId())) == 0:
        raise RuntimeError("Input species tree has only one node.")

    # Validate the outgroups
    if options.outgroupNames is not None:
        projNames = set([tree.getName(x) for x in tree.getLeaves()])
        options.outgroupNames = set(options.outgroupNames.split(","))
        for outgroupName in options.outgroupNames:
            if outgroupName not in projNames:
                raise RuntimeError("Specified outgroup %s not found in tree" % outgroupName)

    genomes = expTemplate.getGenomesWithSequence()
    sequences = []
    for genome in genomes:
        sequences.append(expTemplate.getSequencePath(genome))

    mcProj = createMCProject(tree, sequences, expTemplate, confTemplate, options)

    # Replace the sequences with output sequences
    newSequences = CactusPreprocessor.getOutputSequenceFiles(sequences, expTemplate.getOutputSequenceDir())
    for genome, newSequence in zip(genomes, newSequences):
        expTemplate.setSequencePath(genome, newSequence)

    # Now do the file tree creation
    createFileStructure(mcProj, expTemplate, confTemplate, options)

if __name__ == '__main__':    
    main()
