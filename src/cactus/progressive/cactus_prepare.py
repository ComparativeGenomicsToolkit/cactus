#!/usr/bin/env python3

#Released under the MIT license, see LICENSE.txt

"""Set up a given seqfile to include proprocessed/ancestral sequences, as well as to 
   
"""
import os
from argparse import ArgumentParser
import xml.etree.ElementTree as ET
import copy

from operator import itemgetter

from cactus.progressive.seqFile import SeqFile
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper

def main():
    parser = ArgumentParser()
    parser.add_argument("seqFile", help = "Seq file")
    parser.add_argument("outputSequenceDir", help='Directory where the processed leaf sequence and ancestral sequences will be placed')
    parser.add_argument("outSeqFile", help = "Path for annotated Seq file output")
    parser.add_argument("--configFile", default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))

    options = parser.parse_args()

    cactusPrepare(options.seqFile, options.outputSequenceDir, options.outSeqFile, options.configFile)

def cactusPrepare(seqFilePath, outSeqDir, outSeqFilePath, configFilePath):
    """ annotate a SeqFile with ancestral names as well as paths for output sequences."""

    # read the input
    seqFile = SeqFile(seqFilePath)
    configNode = ET.parse(configFilePath).getroot()
    config = ConfigWrapper(configNode)

    # prepare output sequence directory
    # todo: support remote (ie s3) output directory
    try:
        os.makedirs(outSeqDir)
    except:
        pass
    if not os.path.isdir(outSeqDir):
        raise RuntimeError('Unable to create output sequence directory \'{}\''.format(outSeqDir))
    if not os.access(outSeqDir, os.W_OK):
        logger.warning('Output sequence directory is not writeable: \'{}\''.format(outSeqDir))

    # get the ancestor names
    tree = MultiCactusTree(seqFile.tree)
    tree.nameUnlabeledInternalNodes(prefix = config.getDefaultInternalNodePrefix())

    # make the output
    outSeqFile = SeqFile()
    outSeqFile.tree= tree
    outSeqFile.pathMap = seqFile.pathMap
    outSeqFile.outgroups = seqFile.outgroups

    # update paths for preprocessed leaves or inferred ancestors
    preprocess = len(configNode.findall('preprocessor')) > 0
    for node in outSeqFile.tree.breadthFirstTraversal():
        name = outSeqFile.tree.getName(node)
        leaf = outSeqFile.tree.isLeaf(node)
        if (leaf and preprocess) or (not leaf and name not in seqFile.pathMap):
            out_basename = seqFile.pathMap[name] if name in seqFile.pathMap else '{}.fa'.format(name)
            outSeqFile.pathMap[name] = os.path.join(outSeqDir, os.path.basename(out_basename))

    # write the output
    with open(outSeqFilePath, 'w') as out_sf:
        out_sf.write(str(outSeqFile))

    # todo: print out some kind of alignment plan based on a tree decomposition that the user
    # can run on their own

if __name__ == '__main__':
    main()
