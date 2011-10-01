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
from cactus.progressive.experimentWrapper import ExperimentWrapper
from cactus.progressive.outgroup import GreedyOutgroup
from sonLib.nxnewick import NXNewick

prepCmdLine = 'cactus_addFastaHeaderDots.py TARGET_FILE OUT_FILE --event EVENT_STRING'

def createMCProject(tree, options):
    mcTree = MultiCactusTree(tree, options.subtreeSize)
    mcTree.nameUnlabeledInternalNodes(options.prefix)
    mcTree.computeSubtreeRoots()
    mcProj = MultiCactusProject()
    mcProj.mcTree = mcTree
    if options.selfAlign:
        mcTree.addSelfEdges()
    for name in mcProj.mcTree.getSubtreeRootNames():
        expPath = "%s/%s/%s_experiment.xml" % (options.path, name, name)
        mcProj.expMap[name] = os.path.abspath(expPath)
    if options.outgroup or options.outgroupLeaves:
        mcProj.outgroup = GreedyOutgroup()
        mcProj.outgroup.importTree(mcProj.mcTree)
        mcProj.outgroup.greedy(options.outgroupLeaves)
    return mcProj

# progressive alignment relies on input sequences having dots and being 
# unique.  This can be enforced using the preprocessor
def addPreprocessor(configElem):
    found = False
    prepNodes = configElem.findall("preprocessor")
    for prep in prepNodes:
        if "preprocessorString" in prep.attrib:
            pString = prep.attrib["preprocessorString"]
            if pString.find(prepCmdLine) >= 0:
                found = True
                break
    if found == False:
        prep = ET.SubElement(configElem, "preprocessor")
        prep.attrib["chunkSize"] = "10000"
        prep.attrib["preprocessorString"] = prepCmdLine

# go through the tree (located in the template experimet)
# and make sure event names are unique up unitil first dot
def cleanEventTree(experiment):
    tree = experiment.getTree()
    eventIds = set()
    for node in tree.breadthFirstTraversal():
        if tree.hasName(node):
            name = tree.getName(node)
            if '.' in name:
                newName = name.replace('.', '_')
                sys.stderr.write('WARNING renaming event %s to %s\n' %(name, newName))
                tree.setName(node, newName)
                name = newName
            if name in eventIds:
                 raise RuntimeException('Duplicate event in tree: %s' % name)
            eventIds.add(name)
            parent = tree.getParent(node)
            if parent is not None:
                weight = tree.getWeight(parent, node)
            if weight is None:
                raise RuntimeException('Missing branch length in species_tree tree')
    experiment.xmlRoot.attrib["species_tree"] = NXNewick().writeString(tree)
    experiment.seqMap = experiment.buildSequenceMap()

# Make the subdirs for each subproblem:  name/ and name/name_DB
# and write the experiment files
# and copy over a config with updated reference field
def createFileStructure(mcProj, expTemplate, options):
    os.makedirs(options.path)
    mcProj.writeXML(os.path.join(options.path, "%s_project.xml" % options.name))
    baseConfig = expTemplate.getConfigPath()
    baseConfigXML = ET.parse(baseConfig).getroot()
    seqMap = expTemplate.seqMap
    portOffset = 0
    for name, expPath in mcProj.expMap.items():
        path = os.path.join(options.path, name)
        seqMap[name] = os.path.join(path, name + '.fa')
    for name, expPath in mcProj.expMap.items():
        path = os.path.join(options.path, name)
        subtree = mcProj.mcTree.extractSubTree(name)
        exp = copy.deepcopy(expTemplate)
        exp.setDbDir(os.path.join(path, "%s_DB" % name))
        if expTemplate.getDbType() == "kyoto_tycoon" and \
            os.path.splitext(name)[1] != ".kch":
            exp.setDbName("%s.kch" % name)
        else:
            exp.setDbName(name)
        if expTemplate.getDbType() == "kyoto_tycoon":
            exp.setDbPort(expTemplate.getDbPort() + portOffset)
            portOffset += 1
        exp.setReferencePath(os.path.join(path, name + '.fa'))
        exp.setMAFPath(os.path.join(path, "%s.maf" % name))
        exp.updateTree(subtree, seqMap)
        exp.setConfigPath(os.path.join(path, "%s_config.xml" % name))
        if options.outgroup and name in mcProj.outgroup.ogMap:
            og, ogDist = mcProj.outgroup.ogMap[name]
            if og in expTemplate.seqMap:
                ogPath = expTemplate.seqMap[og]
            else:
                ogPath = os.path.join(options.path, og)
                ogPath = os.path.join(ogPath, refFileName(og))
            exp.addOutgroup(og, ogPath, ogDist)
        os.makedirs(exp.getDbDir())
        exp.writeXML(expPath)
        configElem = copy.deepcopy(baseConfigXML)
        refElem = configElem.find("reference")
        refElem.attrib["reference"] = name
        if options.fixNames:
            addPreprocessor(configElem)
        ET.ElementTree(configElem).write(exp.getConfigPath()) 
        
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
    parser.add_option("--subtreeSize", dest="subtreeSize", type="int", 
                      help="Max number of sequences to align at a time [default=2]", 
                      default=2)
    parser.add_option("--ancestorPrefix", dest="prefix", type="string",
                      help="Name to assign unlabeled tree nodes [default=\"Anc\"]",
                      default="Anc")
    parser.add_option("--outgroup", dest="outgroup", action="store_true", 
                      default = False, help="Assign outgroups")
    parser.add_option("--outgroupLeaves", dest="outgroupLeaves", action="store_true",
                      default = False, help="Assign (only leaves as) outgroups")
    parser.add_option("--self", dest="selfAlign", action="store_true", 
                      default = False, help="Align each sequence to itself")
    parser.add_option("--fixNames", dest="fixNames",  default = "True", 
                      help="try to make sequence and event names MAF-compliant [default=true]")
    
    options, args = parser.parse_args()
    
    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    options.expFile = args[0]    
    options.path = os.path.abspath(args[1])
    options.name = os.path.basename(options.path)
    options.fixNames = not options.fixNames.lower() == "false"
    assert options.subtreeSize > 1

    if os.path.isdir(options.path) or os.path.isfile(options.path):
        raise RuntimeError("Output project path %s exists\n" % options.path)
    
    expTemplate = ExperimentWrapper(ET.parse(options.expFile).getroot())
    if options.fixNames:
        cleanEventTree(expTemplate)
    checkInputSequencePaths(expTemplate) 
    mcProj = createMCProject(expTemplate.getTree(), options)
    createFileStructure(mcProj, expTemplate, options)
   # mcProj.check()
    return 0
    

if __name__ == '__main__':    
    main()
