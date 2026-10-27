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
from cactus.shared.common import normalize_pangenome_name, pangenome_name_to_pansn_prefix, split_pangenome_name
from cactus.shared.common import pangenome_file_name, split_dot_haplotype

from toil.statsAndLogging import logger
from toil.statsAndLogging import set_logging_from_options
from toil.realtimeLogger import RealtimeLogger

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
    def __init__(self, path=None, defaultBranchLen=1.0, root_name = None, pangenome = False, dot_haplotypes = None):
        self.branchLen = defaultBranchLen        
        if path is not None:
            self.parseFile(path, root_name, pangenome, dot_haplotypes)

    def parseFile(self, path, root_name = None, pangenome = False, dot_haplotypes = None):
        """ with pangenome, the genome names are PanSN (see normalize_pangenome_name), and
        dot_haplotypes (<graphmap dotHaplotypes>, which has to be given) rewrites the names Cactus
        used to give a haplotype with (HG002.1) to it (HG002#1) """
        assert not pangenome or dot_haplotypes is not None, 'a pangenome seqfile is read with its config\'s dotHaplotypes'
        if not os.path.isfile(path):
            raise RuntimeError("File not found: %s" % path)
        self.tree = None
        self.pathMap = dict()
        self.outgroups = []
        self.seqOrder = []
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
                    self.seqOrder.append(name)
                elif len(tokens) > 0:
                    sys.stderr.write("Skipping line %s\n" % l)

        if self.tree is None:
            self.starTree()
        self.cleanTree(root_name)
        self.validate()
        if pangenome:
            self.normalizePangenomeNames(dot_haplotypes)

    def normalizePangenomeNames(self, dot_haplotypes):
        """ rewrite every genome name to its normalize_pangenome_name() spelling, warning about any
        read the old way (HG002.1 as HG002#1), and make sure no two names are the same PanSN genome
        (HG002.1 and HG002#1, or CHM13 and CHM13#0) or would give the same file names
        (GCA_054883195.1#0 and GCA_054883195#1).  ancestor names in a given tree are left alone """
        renamed = {}
        by_prefix = {}
        by_file_name = {}
        for name in self.pathMap:
            new_name = normalize_pangenome_name(name, dot_haplotypes)
            if new_name != name:
                renamed[name] = new_name
            prefix = pangenome_name_to_pansn_prefix(new_name)
            if prefix in by_prefix:
                raise RuntimeError('Genome names {} and {} are both haplotype {} of sample {}: '
                                   'each sample#haplotype can only be given once'.format(
                                       by_prefix[prefix], name, *split_pangenome_name(new_name)[::-1]))
            by_prefix[prefix] = name
            file_name = pangenome_file_name(new_name)
            if file_name in by_file_name:
                raise RuntimeError('Genome names {} and {} cannot both be used: files are named after '
                                   'their genomes with any "#" made a ".", which would name both {}.  '
                                   'Please rename one of them'.format(by_file_name[file_name], name, file_name))
            by_file_name[file_name] = name
        dot_renamed = {name : new_name for name, new_name in renamed.items() if '#' not in name}
        if dot_renamed:
            logger.warning('Haplotypes are now given in the PanSN style, SAMPLE#HAPLOTYPE, and names '
                           'are no longer read as SAMPLE.HAPLOTYPE except for this backward compatible '
                           'rewrite, which is also how they will be spelled in the output (HAL, MAF '
                           'etc.): {}.  Please use the new names in the seqfile.  A sample whose name '
                           'really ends in .N, like a versioned accession, must give its haplotype '
                           'explicitly, ex GCA_054883195.1#0, or set <graphmap dotHaplotypes="0"> in '
                           'the config to take every name without a "#" as it is'.format(
                               ', '.join('{} -> {}'.format(k, v) for k, v in dot_renamed.items())))
        elif not dot_haplotypes:
            # with dotHaplotypes="0" (some day the default), names in the old spelling are samples
            # in their own right, which is quietly wrong for the haplotypes of a diploid
            by_stem = {}
            for name in self.pathMap:
                if split_dot_haplotype(name):
                    by_stem.setdefault(split_dot_haplotype(name)[0], []).append(name)
            lookalikes = [names for names in by_stem.values() if len(names) > 1]
            if lookalikes:
                logger.warning('Genome names {} are each a haploid sample of their own, as <graphmap '
                               'dotHaplotypes="0"> takes a name without a "#" as it is.  If they are '
                               'the haplotypes of one sample, name them SAMPLE#1, SAMPLE#2 etc.'.format(
                                   ', '.join(' '.join(names) for names in lookalikes)))
        if not renamed:
            return
        self.pathMap = {renamed.get(name, name) : path for name, path in self.pathMap.items()}
        self.seqOrder = [renamed.get(name, name) for name in self.seqOrder]
        self.outgroups = [renamed.get(name, name) for name in self.outgroups]
        for node in self.tree.postOrderTraversal():
            if self.tree.isLeaf(node) and self.tree.getName(node) in renamed:
                self.tree.setName(node, renamed[self.tree.getName(node)])

    def starTree(self):
        self.tree = NXTree()
        label = 0
        self.tree.nxDg.add_node(label)
        self.tree.rootId = label
        for name in list(self.pathMap.keys()):
            label += 1
            self.tree.nxDg.add_edge(0, label)
            self.tree.setName(label, name)
            self.tree.setWeight(0, label, self.branchLen)

    def validate(self):
        if len([i for i in self.tree.postOrderTraversal()]) <= 2:
            raise RuntimeError("At least two valid leaf genomes required in"
                               " input tree")
        nameSet = set()
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
                nameSet.add(name)
            elif self.tree.hasName(node):
                name = self.tree.getName(node)
                if name.isnumeric():
                    sys.stderr.write("WARNING: ignoring numeric ancestor name {} in input tree\n".format(name))
                    self.tree.setName(node, None)
                elif name in nameSet:
                    sys.stderr.write("WARNING: ignoring duplicate ancestor name {} in input tree\n".format(name))
                    self.tree.setName(node, None)
                else:
                    nameSet.add(name)

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
    def cleanTree(self, root_name = None):
        numLeaves = 0
        removeList = []
        rootID = self.tree.getRootId()
        if root_name:
            for node in self.tree.breadthFirstTraversal():
                if self.tree.hasName(node) and self.tree.getName(node) == root_name:
                    rootID = node
                    break

        # we will never clean inside our root
        mustKeepSet = set()
        for node in self.tree.breadthFirstTraversal(root = rootID):
            mustKeepSet.add(node)
            
        for node in self.tree.postOrderTraversal():
            if self.tree.isLeaf(node):
                name = self.tree.getName(node)
                if name not in self.pathMap and node not in mustKeepSet:
                    removeList.append(node)
                numLeaves += 1
        if numLeaves < 2:
            raise RuntimeError("At least two valid leaf genomes required in"
                               " input tree")
        if len(removeList) == numLeaves:
            raise RuntimeError("No sequence path specified for any leaves in the tree")
        for leaf in removeList:
             logger.warning("No sequence path found for species %s: removing from tree\n" % (
                 self.tree.getName(leaf)))             
             self.tree.removeLeaf(leaf)

        for node in self.tree.postOrderTraversal():
            if self.tree.hasParent(node):
                parent = self.tree.getParent(node)
                if self.tree.getWeight(parent, node) is None:
                    logger.warning(
                        "No branch length for %s: setting to %d\n" % (
                            self.tree.getName(node), self.branchLen))
                    self.tree.setWeight(parent, node, self.branchLen)

    # create the cactus_workflow_experiment xml element which serves as
    # the root node of the experiment template file needed by
    # cactus_createMultiCactusProject.  Note the element is incomplete
    # until the cactus_disk child element has been added
    def toXMLElement(self, ignoreSeqPaths=[]):
        assert self.tree is not None
        elem = ET.Element("cactus_workflow_experiment")
        for node in self.tree.postOrderTraversal():
            name = self.tree.getName(node)
            if name in self.pathMap:
                path = self.pathMap[name]
                genomeNode = ET.SubElement(elem, "genome")
                genomeNode.attrib['name'] = name
                if name not in ignoreSeqPaths:
                    genomeNode.attrib['sequence'] = path
        elem.attrib["species_tree"] = NXNewick().writeString(self.tree)
        elem.attrib["config"] = "defaultProgressive"
        return elem

    # convert to string
    def __str__(self):
        og_set = set(self.outgroups)
            
        s = NXNewick().writeString(self.tree)
        s += '\n'
        for name, path in self.pathMap.items():
            if name in og_set:
                s += '*'
            s += '{}\t{}\n'.format(name, path)

        return s
            
