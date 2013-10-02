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
from cactus.progressive.outgroup import GreedyOutgroup
from sonLib.nxnewick import NXNewick

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
    if config.getOutgroupStrategy() == 'greedy':
        mcProj.outgroup = GreedyOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree)
        mcProj.outgroup.greedy(threshold=config.getOutgroupThreshold(),
                               candidateSet=options.outgroupNames,
                               candidateChildFrac=config.getOutgroupAncestorQualityFraction())
    elif config.getOutgroupStrategy() == 'greedyLeaves':
        mcProj.outgroup = GreedyOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree)
        ogSet = options.outgroupNames
        if ogSet is None:
            ogSet = set([mcProj.mcTree.getName(x) for x in mcProj.mcTree.getLeaves()])
        mcProj.outgroup.greedy(threshold=config.getOutgroupThreshold(),
                               candidateSet=ogSet,
                               candidateChildFrac=2.0)
    return mcProj

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
                if node1 != node2 and name1.find(name2) == 0 and\
                 name1 != name2 + tree.self_suffix:
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
            og, ogDist = mcProj.outgroup.ogMap[name]
            if og in expTemplate.seqMap:
                ogPath = expTemplate.seqMap[og]
            else:
                ogPath = os.path.join(options.path, og)
                ogPath = os.path.join(ogPath, refFileName(og))
            exp.addTheOutgroupSequence(og, ogDist, ogPath)
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
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    options.expFile = args[0]    
    options.path = os.path.abspath(args[1])
    options.name = os.path.basename(options.path)
    options.fixNames = not options.fixNames.lower() == "false"

    if os.path.isdir(options.path) or os.path.isfile(options.path):
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
    expTemplate.setSequences([ os.path.join(expTemplate.getOutputSequenceDir(), os.path.split(i)[-1]) for i in expTemplate.getSequences() ])
    #Now do the file tree creation
    createFileStructure(mcProj, expTemplate, confTemplate, options)
   # mcProj.check()
    return 0
    

if __name__ == '__main__':    
    main()
