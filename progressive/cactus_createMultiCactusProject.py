#!/usr/bin/env python

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt
#!/usr/bin/env python

"""Create the multi_cactus xml and directory structure from a workflow template
""" 
import os
import sys
from optparse import OptionParser
import xml.etree.ElementTree as ET
import copy

from cactus.progressive.multiCactusProject import MultiCactusProject
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.experimentWrapper import ExperimentWrapper
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.outgroup import GreedyOutgroup, DynamicOutgroup
from sonLib.nxnewick import NXNewick
from cactus.preprocessor.cactus_preprocessor import CactusPreprocessor

def createMCProject(tree, experiment, config, options):
    mcTree = MultiCactusTree(tree, config.getSubtreeSize())
    mcTree.nameUnlabeledInternalNodes(config.getDefaultInternalNodePrefix())
    mcTree.computeSubtreeRoots()
    mcProj = MultiCactusProject()
    mcProj.mcTree = mcTree
    mcProj.inputSequences = experiment.getSequences()[:] 
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
    mcProj.outgroup = None
    if config.getOutgroupStrategy() == 'greedy':
        # use the provided outgroup candidates, or use all outgroups
        # as candidates if none are given
        mcProj.outgroup = GreedyOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree, alignmentRootId)
        mcProj.outgroup.greedy(threshold=config.getOutgroupThreshold(),
                               candidateSet=options.outgroupNames,
                               candidateChildFrac=config.getOutgroupAncestorQualityFraction(),
                               maxNumOutgroups=config.getMaxNumOutgroups())
    elif config.getOutgroupStrategy() == 'greedyLeaves':
        # use all leaves as outgroups, unless outgroup candidates are given
        mcProj.outgroup = GreedyOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree, alignmentRootId)
        ogSet = options.outgroupNames
        if ogSet is None:
            ogSet = set([mcProj.mcTree.getName(x) for x in mcProj.mcTree.getLeaves()])
        mcProj.outgroup.greedy(threshold=config.getOutgroupThreshold(),
                               candidateSet=ogSet,
                               candidateChildFrac=2.0,
                               maxNumOutgroups=config.getMaxNumOutgroups())
    elif config.getOutgroupStrategy() == 'greedyPreference':
        # prefer the provided outgroup candidates, if any, but use
        # other nodes as "filler" if we can't find enough.
        mcProj.outgroup = GreedyOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree, alignmentRootId)
        mcProj.outgroup.greedy(threshold=config.getOutgroupThreshold(),
                               candidateSet=options.outgroupNames,
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
                                   candidateSet=options.outgroupNames)
        mcProj.outgroup.compute(maxNumOutgroups=config.getMaxNumOutgroups())
    elif config.getOutgroupStrategy() != 'none':
        raise RuntimeError("Could not understand outgroup strategy %s" % config.getOutgroupStrategy())

    # if necessary, we reroot the tree at the specified alignment root id.  all leaf genomes
    # that are no longer in the tree, but still used as outgroups, are moved into special fields
    # so that we can remember to, say, get their paths for preprocessing. 
    specifyAlignmentRoot(mcProj, alignmentRootId)
    return mcProj

# it is possible that we start with a much bigger tree than we actually want to align
# (controlled by the --root option in cactus_createMultiCactusProject.py).  We use
# the entire tree when selecting outgroups, but right afterward have no use for
# genomes that are neither outgroups, nor in the alignment.  We especially don't
# want to waste time preprocessing them.  This function, reroots the tree at the
# alignment root then tacks on all the outgroups from ooutside the new tree
# (which must be leaves) as children of the root
def specifyAlignmentRoot(mcProj, alignmentRootId):
    # ugly hack to keep track of external outgroups for root experiment (yuck)
    mcProj.externalOutgroupNames = set()
    # nothing to do
    if alignmentRootId == mcProj.mcTree.getRootId():
        return
    
    # dig out every outgroup
    outGroupNames = set()
    if mcProj.outgroup is not None:
        for event, ogNameDistList in mcProj.outgroup.ogMap.items():
            for og, dist in ogNameDistList:
                outGroupNames.add(og)

    # i don't like this but we have to assume here that the sequences are
    # written in postorder (as in experiments)
    allLeafMap = dict()
    idx = 0
    for node in mcProj.mcTree.postOrderTraversal():
        if mcProj.mcTree.isLeaf(node) is True:
            allLeafMap[mcProj.mcTree.getName(node)] = mcProj.inputSequences[idx]
            idx += 1

    # find outgroups we want to extract
    keptNodes = set([x for x in mcProj.mcTree.postOrderTraversal(alignmentRootId)])
    deadNodes = []
    extractOutgroupMap = dict()
    for node in mcProj.mcTree.postOrderTraversal():
        if node not in keptNodes:
            deadNodes.append(node)
            name = mcProj.mcTree.getName(node)
            if name in outGroupNames:
                assert name in allLeafMap
                extractOutgroupMap[name] = allLeafMap[name]
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
        allLeafMap[ogName] = ogPath

    # remove any experiment directories that have become invalid
    for event in mcProj.expMap.keys():
        if mcProj.mcTree.getNodeId(event) in deadNodes:
            del mcProj.expMap[event]
            
    # flush out all unused nodes, and set the new root
    for node in deadNodes:
        assert mcProj.mcTree.hasParent(node)
        mcProj.mcTree.removeEdge(mcProj.mcTree.getParent(node), node)
    from sonLib.nxnewick import NXNewick

    # reset input sequences to only contain genomes in tree (making sure to
    # keep in postorder)
    mcProj.inputSequences = []
    for node in mcProj.mcTree.postOrderTraversal():
        if mcProj.mcTree.isLeaf(node):
            mcProj.inputSequences.append(allLeafMap[mcProj.mcTree.getName(node)])


# go through the tree (located in the template experimet)
# and make sure event names are unique up unitil first dot
def cleanEventTree(experiment):
    tree = MultiCactusTree(experiment.getTree())
    tree.nameUnlabeledInternalNodes()
    for node in tree.breadthFirstTraversal():
        if tree.hasName(node):
            name = tree.getName(node)
            if '.' in name:
                newName = name.replace('.', '_')
                sys.stderr.write('WARNING renaming event %s to %s\n' %(name, newName))
                tree.setName(node, newName)
                name = newName
            parent = tree.getParent(node)
            if parent is not None:
                weight = tree.getWeight(parent, node)
                if weight is None:
                    raise RuntimeError('Missing branch length in species_tree tree')
    redoPrefix = True
    newSuffix = 0
    while redoPrefix is True:
        redoPrefix = False
        for node1 in tree.breadthFirstTraversal():
            name1 = tree.getName(node1)
            for node2 in tree.breadthFirstTraversal():
                name2 = tree.getName(node2)
                if node1 != node2 and name1 == name2:
                    newName = "%s%i" % (name2, newSuffix)
                    newSuffix += 1
                    tree.setName(node2, newName)
                    sys.stderr.write('WARNING renaming event %s to %s\n' % (
                        name2, newName))
                    redoPrefix = True
                    
    experiment.xmlRoot.attrib["species_tree"] = NXNewick().writeString(tree)
    experiment.seqMap = experiment.buildSequenceMap()

# Make the subdirs for each subproblem:  name/ and name/name_DB
# and write the experiment files
# and copy over a config with updated reference field
def createFileStructure(mcProj, expTemplate, configTemplate, options):
    if not os.path.exists(options.path):
        os.makedirs(options.path)
    mcProj.writeXML(os.path.join(options.path, "%s_project.xml" % options.name))
    seqMap = expTemplate.seqMap
    portOffset = 0
    for name, expPath in mcProj.expMap.items():
        path = os.path.join(options.path, name)
        seqMap[name] = os.path.join(path, name + '.fa')
    for name, expPath in mcProj.expMap.items():
        path = os.path.join(options.path, name)
        subtree = mcProj.mcTree.extractSubTree(name)
        # hack the outgroups back out of the root (yuck)
        if name == mcProj.mcTree.getRootName():
            for rootChild in subtree.getChildren(subtree.getRootId()):
                if subtree.getName(rootChild) in mcProj.externalOutgroupNames:
                    subtree.removeEdge(subtree.getRootId(), rootChild)
        exp = copy.deepcopy(expTemplate)
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
        exp.updateTree(subtree, seqMap)
        exp.setConfigPath(os.path.join(path, "%s_config.xml" % name))
        if configTemplate.getOutgroupStrategy() != 'none' \
        and name in mcProj.outgroup.ogMap:
            for og, ogDist in mcProj.outgroup.ogMap[name]:
                if og in expTemplate.seqMap:
                    ogPath = expTemplate.seqMap[og]
                else:
                    ogPath = os.path.join(options.path, og)
                    ogPath = os.path.join(ogPath, refFileName(og))
                exp.addOutgroupSequence(og, ogDist, ogPath)
        elif name == mcProj.mcTree.getRootName() \
        and options.rootOutgroupPaths is not None:
            exp.xmlRoot.attrib["outgroup_events"] = " ".join([ogName for (ogName, _) in mcProj.outgroup.ogMap[name]])
        if not os.path.exists(exp.getDbDir()):
            os.makedirs(exp.getDbDir())
        if not os.path.exists(path):
            os.makedirs(path)
        exp.writeXML(expPath)
        config = ConfigWrapper(copy.deepcopy(configTemplate.xmlRoot))
        config.setReferenceName(name)
        config.verifyMinBlockDegree(exp)
        config.writeXML(exp.getConfigPath())
        
def checkInputSequencePaths(exp):
    for event, seq in exp.seqMap.items():
        if not os.path.exists(seq):
            sys.stderr.write("WARNING: sequence path %s does not exist\n" % seq)
        elif os.path.isdir(seq):
            contents = os.listdir(seq)
            size = 0
            for i in contents:
                if i[0] != '.':
                    size += 1
            if size == 0:
                sys.stderr.write("WARNING: sequence path %s is an empty directory\n" % seq)
            
def main():
    usage = "usage: %prog [options] <experiment> <output project path>"
    description = "Setup a multi-cactus project using an experiment xml as template"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--fixNames", dest="fixNames",  default = "True", 
                      help="try to make sequence and event names MAF-compliant [default=true]")
    parser.add_option("--outgroupNames", dest="outgroupNames",  default = None, 
                      help="comma-separated names of high quality assemblies to use as outgroups [default=everything]")
    parser.add_option("--rootOutgroupPaths", dest="rootOutgroupPaths", type=str,
                      help="comma-separated root outgroup paths (other root outgroup options " +
                      "must be given as well)", default=None)
    parser.add_option("--rootOutgroupDists", dest="rootOutgroupDists", type=str,
                      help="comma-separated root outgroup distances (other root outgroup " +
                      "options must be given as well", default=None)
    parser.add_option("--root", dest="root", type=str,
                      help="name of alignment root (must be labeled ancestral node in tree in input experiment).  Useful "
                      "for allowing the tree to contain nodes that won't be in the alignment but can still be used for "
                      "outgroups.  NOTE: ADDED TO REPLACE --rootOutgroupPaths and --rootOutgroupDists.",
                      default=None)
    parser.add_option("--overwrite", action="store_true", help="Overwrite existing experiment files", default=False)

    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    options.expFile = args[0]    
    options.path = os.path.abspath(args[1])
    options.name = os.path.basename(options.path)
    options.fixNames = not options.fixNames.lower() == "false"

    if (options.rootOutgroupDists is not None) \
    ^ (options.rootOutgroupPaths is not None):
        parser.error("--rootOutgroupDists and --rootOutgroupPaths must be " +
                         "provided together")

    if options.rootOutgroupDists is not None and len(options.rootOutgroupDists.split(",")) != len(options.rootOutgroupPaths.split(",")):
        parser.error("root outgroup options must have the same number of entries")

    if (os.path.isdir(options.path) and not options.overwrite) or os.path.isfile(options.path):
        raise RuntimeError("Output project path %s exists\n" % options.path)
    
    expTemplate = ExperimentWrapper(ET.parse(options.expFile).getroot())
    configPath = expTemplate.getConfigPath()
    confTemplate = ConfigWrapper(ET.parse(configPath).getroot())
    if options.fixNames:
        cleanEventTree(expTemplate)
    checkInputSequencePaths(expTemplate)
    tree = expTemplate.getTree()
    if options.outgroupNames is not None:
        projNames = set([tree.getName(x) for x in tree.getLeaves()])
        options.outgroupNames = set(options.outgroupNames.split(","))
        for outgroupName in options.outgroupNames:
            if outgroupName not in projNames:
                raise RuntimeError("Specified outgroup %s not found in tree" % outgroupName)
    mcProj = createMCProject(tree, expTemplate, confTemplate, options)
    #Replace the sequences with output sequences
    expTemplate.updateTree(mcProj.mcTree, expTemplate.buildSequenceMap())
    expTemplate.setSequences(CactusPreprocessor.getOutputSequenceFiles(mcProj.inputSequences, expTemplate.getOutputSequenceDir()))
    if options.rootOutgroupPaths is not None:
        # hacky -- add the root outgroup to the tree after everything
        # else.  This ends up being ugly, but avoids running into
        # problems with assumptions about tree structure
        #
        # [this hack is hopefully made redundant by the --root option]
        #
        options.rootOutgroupPaths = [os.path.abspath(path) for path in options.rootOutgroupPaths.split(",")]
        options.rootOutgroupDists = [float(i) for i in options.rootOutgroupDists.split(",")]
        mcProj.inputSequences.extend(options.rootOutgroupPaths)
        # replace the root outgroup sequence by post-processed output
        for i, (outgroupPath, outgroupDist) in enumerate(zip(options.rootOutgroupPaths, options.rootOutgroupDists)):
            outgroupPath = CactusPreprocessor.getOutputSequenceFiles(expTemplate.getSequences() + options.rootOutgroupPaths[:i+1], expTemplate.getOutputSequenceDir())[-1]
            rootOutgroupName = "rootOutgroup%d" % i
            expTemplate.seqMap[rootOutgroupName] = outgroupPath
            # Add to tree
            mcProj.mcTree.addOutgroup(rootOutgroupName, outgroupDist)
            mcProj.mcTree.computeSubtreeRoots()
            if mcProj.mcTree.getRootName() not in mcProj.outgroup.ogMap:
                # initialize ogMap entry
                mcProj.outgroup.ogMap[mcProj.mcTree.getRootName()] = []
            mcProj.outgroup.ogMap[mcProj.mcTree.getRootName()].append((rootOutgroupName, outgroupDist))
    #Now do the file tree creation
    createFileStructure(mcProj, expTemplate, confTemplate, options)
   # mcProj.check()
    return 0
    

if __name__ == '__main__':    
    main()
