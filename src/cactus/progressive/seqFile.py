#!/usr/bin/env python3

# Progressive Cactus Package
# Copyright (C) 2009-2012 by Glenn Hickey (hickey@soe.ucsc.edu)
# and Benedict Paten (benedictpaten@gmail.com)

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


import os
import sys
import xml.etree.ElementTree as ET
import math
import time
import random
import copy
from optparse import OptionParser
from optparse import OptionGroup
import string
import re

from sonLib.bioio import absSymPath
from sonLib.nxtree import NXTree
from sonLib.nxnewick import NXNewick
from sonLib.bioio import popenCatch

# parse the input seqfile for progressive cactus.  this file is in the
# format of:
#   newick tree
#   name sequencePath
#   name sequencePath
#   ...
# example:
#
# (human, (dog, cat));
# human /hive/seq/human
# dog /users/name/fasta/dog.fa
# cat /tmp/cat/
#
# now added * to signify high quality addembly. ex, use human
# and cat as outgroups but never dog:
#
# (human, (dog, cat));
# *human /hive/seq/human
# dog /users/name/fasta/dog.fa
# *cat /tmp/cat/

class SeqFile:
    branchLen = 1
    def __init__(self, path=None):
        if path is not None:
            self.parseFile(path)

    def parseFile(self, path):
        if not os.path.isfile(path):
            raise RuntimeError("File not found: %s" % path)
        self.tree = None
        self.pathMap = dict()
        self.outgroups = []
        seqFile = open(path, "r")
        for l in seqFile:
            line = l.strip()
            if line:
                if line[0] == "#":
                    continue
                tokens = line.split()
                if self.tree is None and (len(tokens) == 1 or line[0] == '('):
                    newickParser = NXNewick()
                    try:
                        self.tree = newickParser.parseString(line)
                    except:
                        raise RuntimeError("Failed to parse newick tree: %s" %
                                           line)
                elif len(tokens) > 0 and tokens[0] == '*':
                    sys.stderr.write("Skipping line %s\n" % l)
                elif line[0] != '(' and len(tokens) >= 2:
                    name = tokens[0]
                    if name[0] == '*':
                        name = name[1:]
                        self.outgroups.append(name)
                    path = " ".join(tokens[1:])
                    if name in self.pathMap:
                        raise RuntimeError("Duplicate name found: %s" % name)
                    self.pathMap[name] = path
                elif len(tokens) > 0:
                    sys.stderr.write("Skipping line %s\n" % l)

        if self.tree is None:
            self.starTree()
        self.cleanTree()
        self.validate()

    def starTree(self):
        self.tree = NXTree()
        label = 0
        self.tree.nxDg.add_node(label)
        self.tree.rootId = label
        for name in list(self.pathMap.keys()):
            label += 1
            self.tree.nxDg.add_edge(0, label)
            self.tree.setName(label, name)
            self.tree.setWeight(0, label, SeqFile.branchLen)

    def validate(self):
        if len([i for i in self.tree.postOrderTraversal()]) <= 2:
            raise RuntimeError("At least two valid leaf genomes required in"
                               " input tree")
        for node in self.tree.postOrderTraversal():
            if self.tree.isLeaf(node):
                name = self.tree.getName(node)
                if name not in self.pathMap:
                    raise RuntimeError("No sequence specified for %s" % name)
                else:
                    path = self.pathMap[name]
                    #if not os.path.exists(path):
                    #    raise RuntimeError("Sequence path not found: %s" % path)
                    #self.sanityCheckSequence(path)

    def sanityCheckSequence(self, path):
        """Warns the user about common problems with the input sequences."""
        # Relies on cactus_analyseAssembly output staying in the
        # format it's currently in.
        cmdline = "cactus_analyseAssembly"
        if os.path.isdir(path):
            cmdline = "cat %s/* | %s -" % (path, cmdline)
        else:
            cmdline += " %s" % path
        output = popenCatch(cmdline)
        try:
            repeatMaskedFrac = float(re.search(r'Proportion-repeat-masked: ([0-9.]*)', output).group(1))
            nFrac = float(re.search(r'ProportionNs: ([0-9.]*)', output).group(1))
        except ValueError:
            # This can happen if the genome has 0 length, making the fractions NaN.
            # We warn the user but return afterwards, as the rest of the checks are
            # dependent on the fraction values.
            sys.stderr.write("WARNING: sequence path %s has 0 length. Consider "
                             "removing it from your input file.\n\n" % path)
            return
        # These thresholds are pretty arbitrary, but should be good for
        # badly- to well-assembled vertebrate genomes.
        if repeatMaskedFrac > 0.70:
            sys.stderr.write("WARNING: sequence path %s has an extremely high "
                             "proportion of masked bases: %f. progressiveCactus"
                             " expects a soft-masked genome, i.e. all lowercase"
                             " characters are considered masked. The process "
                             "will proceed normally, but make sure you haven't "
                             "accidentally provided an all-lowercase genome, "
                             "in which case nothing will be aligned to "
                             "it!\n\n" % (path, repeatMaskedFrac))
        if nFrac > 0.30:
            sys.stderr.write("WARNING: sequence path %s has an extremely high "
                             "proportion of 'N' bases: %f. The process will "
                             "proceed normally, but make sure your genome "
                             "isn't hard-masked! Alignments to hard-masked "
                             "genomes are much worse than to soft-masked "
                             "genomes. If the genome just has a lot of "
                             "poorly assembled regions, feel free to "
                             "ignore this message.\n\n" % (path, nFrac))

    # remove leaves that do not have sequence data associated with them
    def cleanTree(self):
        numLeaves = 0
        removeList = []
        for node in self.tree.postOrderTraversal():
            if self.tree.isLeaf(node):
                name = self.tree.getName(node)
                if name not in self.pathMap:
                    removeList.append(node)
                numLeaves += 1
        if numLeaves < 2:
            raise RuntimeError("At least two valid leaf genomes required in"
                               " input tree")
        if len(removeList) == numLeaves:
            raise RuntimeError("No sequence path specified for any leaves in the tree")
        for leaf in removeList:
             sys.stderr.write("No sequence path found for %s: skipping\n" % (
                 self.tree.getName(leaf)))
             self.tree.removeLeaf(leaf)

        for node in self.tree.postOrderTraversal():
            if self.tree.hasParent(node):
                parent = self.tree.getParent(node)
                if self.tree.getWeight(parent, node) is None:
                    sys.stderr.write(
                        "No branch length for %s: setting to %d\n" % (
                            self.tree.getName(node), SeqFile.branchLen))
                    self.tree.setWeight(parent, node, SeqFile.branchLen)


    # create the cactus_workflow_experiment xml element which serves as
    # the root node of the experiment template file needed by
    # cactus_createMultiCactusProject.  Note the element is incomplete
    # until the cactus_disk child element has been added
    def toXMLElement(self):
        assert self.tree is not None
        elem = ET.Element("cactus_workflow_experiment")
        for node in self.tree.postOrderTraversal():
            name = self.tree.getName(node)
            if name in self.pathMap:
                path = self.pathMap[name]
                genomeNode = ET.SubElement(elem, "genome")
                genomeNode.attrib['name'] = name
                genomeNode.attrib['sequence'] = path
        elem.attrib["species_tree"] = NXNewick().writeString(self.tree)
        elem.attrib["config"] = "defaultProgressive"
        return elem
