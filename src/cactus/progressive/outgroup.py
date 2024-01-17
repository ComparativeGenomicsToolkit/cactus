#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Compute outgroups for subtrees by greedily finding
the nearest valid candidate.  has option to only assign
leaves as outgroups

"""

import os
import sys
import math
import copy
import itertools
import networkx as NX
import xml.etree.ElementTree as ET
from collections import defaultdict, namedtuple
from optparse import OptionParser
from sonLib.nxnewick import NXNewick
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import cactus_call
from toil.statsAndLogging import logger

class GreedyOutgroup(object):
    def __init__(self):
        self.dag = None
        self.dm = None
        self.dmDirected = None
        self.root = None
        self.ogMap = None
        self.mcTree = None
        self.chrom_map = None

    # add edges from sonlib tree to self.dag
    # compute self.dm: an undirected distance matrix
    def importTree(self, mcTree, rootId = None):
        # input tree
        self.mcTree = mcTree
        # initialize digraph with input tree
        self.dag = mcTree.nxDg.copy()
        # root of input tree
        self.root = mcTree.rootId
        # remove nodes that aren't events (probably unused in practice)
        self.stripNonEvents(self.root, mcTree.subtreeRoots)
        # distance matrix of tree (directed)
        self.dmDirected = dict(NX.algorithms.shortest_paths.weighted.\
        all_pairs_dijkstra_path_length(self.dag))
        # nodes that are outside scope and should just be ignored
        self.invalidSet = self.getInvalid(rootId)
        # undirected graph
        graph = NX.Graph(self.dag)
        # undirected distance matrix
        self.dm = dict(NX.algorithms.shortest_paths.weighted.\
        all_pairs_dijkstra_path_length(graph))
        # the results: map tree events to otugroups
        self.ogMap = defaultdict(list)

    # return set of ancestral nodes that aren't below alignment root
    # they can't be outgroups as they are effective "out of project"
    def getInvalid(self, rootId):
        invalid = []
        if rootId is None or rootId == self.mcTree.getRootId():
            return invalid
        good = set([x for x in self.mcTree.postOrderTraversal(rootId)])
        for node in self.mcTree.postOrderTraversal():
            if not self.mcTree.isLeaf(node) and node not in good:
                invalid.append(node)
        return set(invalid)

    # get rid of any node that's not an event
    def stripNonEvents(self, id, subtreeRoots):
        children = []
        for outEdge in self.dag.out_edges(id):
            children.append(outEdge[1])
        parentCount = len(self.dag.in_edges(id))
        assert parentCount <= 1
        parent = None
        if parentCount == 1:
            parent = self.dag.in_edges[0][0]
        if id not in subtreeRoots and len(children) >= 1:
            assert parentCount == 1
            self.dag.remove_node(id)
            for child in children:
                self.dag.add_edge(parent, child)
                self.stripNonEvents(child, subtreeRoots)

    # are source and sink son same path to to root? if so,
    # they shouldn't be outgroups of each other.
    def onSamePath(self, source, sink):
        if source in self.dmDirected:
            if sink in self.dmDirected[source]:
                return True
        if sink in self.dmDirected:
            if source in self.dmDirected[sink]:
                return True
        return False

    # fill up a dictionary of node id -> height in tree where
    # leaves have height = 0
    def heightTable(self):
        htable = dict()
        def rec(node):
            children = [x[1] for x in self.dag.out_edges(node)]

            # Update the table for those children not already in the
            # table.
            for child in children:
                if child not in htable:
                    rec(child)

            if len(children) == 0:
                htable[node] = 0
            else:
                htable[node] = max([htable[i] for i in children]) + 1
        rec(self.root)
        return htable
    
    # check the candidate using the set and and fraction
    def inCandidateSet(self, node, candidateChildFrac):
        if self.candidateMap is None or len(self.candidateMap) == 0:
            return True
        if self.mcTree.getName(node) in self.candidateMap:
            return self.candidateMap[self.mcTree.getName(node)]
        children = self.mcTree.breadthFirstTraversal(node)
        leaves = []
        for child in children:
            if self.mcTree.isLeaf(child):
                leaves.append(child)
        candidateLeaves = 0
        for leaf in leaves:
            if self.mcTree.getName(leaf) in self.candidateMap:
                candidateLeaves += 1
        if len(leaves) == 0:
            self.candidateMap[self.mcTree.getName(node)] = False
            return False
        frac = float(candidateLeaves) / float(len(leaves))
        if frac >= candidateChildFrac:
            self.candidateMap[self.mcTree.getName(node)] = True
            return True
        self.candidateMap[self.mcTree.getName(node)] = False
        return False

    def loadChromInfo(self, chrom_info_name):
        """
        load the chrom info.  it's like a seqfile, but maps genome names
        to lists of chromosomes.  not all genomes need to be in it
        """
        self.chrom_map = {}
        
        if not os.path.isfile(chrom_info_name):
            raise RuntimeError('Unable to open chromosome info file {}'.format(chrom_info_name))

        with open(chrom_info_name, 'r') as chrom_info_file:
            for line in chrom_info_file:
                toks = line.rstrip().split()
                if len(toks) == 2:
                    genome_name = toks[0]
                    chroms = toks[1]
                    self.chrom_map[genome_name] = chroms.split(',')
                elif len(toks):
                    raise RuntimeError('Unable to parse line in {}, expecting 2 columns: {}'.format(chrom_info_name, line))
                    
    def check_chrom_satisfied(self, node_id, id_to_chroms):
        """
        see if the current outgroup selection for the given node has at least
        one copy of every chromosome in/under that node.  this check is only
        relevant if input chroms have been specified...

        these chromosomes are strings like "X", "Y" etc.  that can
        be used to denote the sex of certain leaf assemblies.  This
        will then be used to guide outgroup selection, where outgroups
        are chosen so that at least two copies of each chromosome
        are present in the ingroup+outgroup set.  

        """
        node_chrom_set = id_to_chroms[node_id]
        if not node_chrom_set:
            return True
        outgroup_names = self.ogMap[self.mcTree.getName(node_id)]
        og_ids = [self.mcTree.getNodeId(og_name[0]) for og_name in outgroup_names]
        og_chroms = set()
        for og_id in og_ids:
            for og_chrom in id_to_chroms[og_id]:
                og_chroms.add(og_chrom)
        return node_chrom_set.issubset(og_chroms)

    def refine_og_chroms(self, node_chroms, max_outgroups):
        """
        go over the outgroup assignment and try to refine the selections to
        best represent the chromosome specificaitons
        """
        for source_name, og_names_dists in self.ogMap.items():
            source_id = self.mcTree.getNodeId(source_name)
            source_chroms = node_chroms[source_id]            
            # only need to do anything if we have chromosomes
            if source_chroms:
                chrom_set = set()
                essential_ogs = []
                for og_name, og_dist in og_names_dists:
                    og_id = self.mcTree.getNodeId(og_name)
                    og_chroms = node_chroms[og_id]
                    is_essential = False
                    for og_chrom in og_chroms:
                        if og_chrom in source_chroms and og_chrom not in chrom_set:
                            # this outgroup brings with it a needed chromosome we haven't seen yet,
                            # so we mark it as "essential"
                            is_essential = True
                            chrom_set.add(og_chrom)
                    if is_essential:
                        essential_ogs.append((og_name, og_dist))

                # todo: we can/should use some kind of threshold to make sure we
                # don't take rediculously distance outgroups just because they
                # have the right chromosomes
                if len(essential_ogs) > max_outgroups:
                    logger.warning("Warning: Limiting outgroups for {} to {}, which means leaving out {}".format(
                        source_name, max_outgroups, ",".join([dog[0] for dog in essential_ogs[max_outgroups:]])))
                    # todo: rank by how many chromosomes they have!
                    essential_ogs = essential_ogs[:max_outgroups]
                
                if len(chrom_set) < len(source_chroms):
                    logger.warning("Warning: Unable to fulfull chromosome requirements when selecting outgroups for {}".format(source_name) +
                                   ". In particular, these chromosomes could not be found: {}".format(",".join(source_chroms - chrom_set)))
                # if we have room for outgroups that aren't listed "essential" for their chromosome content,
                # then pad by distance
                if len(essential_ogs) < max_outgroups:
                    padded_ogs = []
                    padding = max_outgroups - len(essential_ogs)
                    for og_name, og_dist in og_names_dists:
                        is_essential = (og_name, og_dist) in essential_ogs
                        if is_essential or padding > 0:
                            padded_ogs.append((og_name, og_dist))
                            if not is_essential:
                                padding -= 1
                    essential_ogs = padded_ogs

                # update the outgroupss
                self.ogMap[source_name] = essential_ogs                    

    # greedily assign closest possible valid outgroups
    # If some outgroups are already assigned, keep the existing
    # assignments but attempt to add more, if possible.
    # all outgroups are stored in self.ogMap
    # edges between leaves ARE NOT kept in the dag
    # the threshold parameter specifies how much parallelism can
    # be sacrificed by the selection of an outgroup
    # threshold = None : just greedy with no constraints
    # threshold = 0 : depth of schedule guaranteed to be unaffected by outgroup
    # threshold = k : depth increases by at most k per outgroup
    # candidateSet : names of valid outgroup genomes. (all if None)
    # candidateChildFrac : min fraction of children of ancestor in
    # candidateSet in order for the ancestor to be an outgroup candidate
    # if > 1, then only members of the candidate set and none of their
    # ancestors are chosen
    # maxNumOutgroups : max number of outgroups to put in each entry of self.ogMap
    # maxNumChromOutgroups : max number of outgroups when trying to make chromosome assignment
    # should be bigger than maxNumOutgroups so that the genomes with rare chromosomes can
    # get the extra help they need.
    def greedy(self, threshold = None, candidateSet = None,
               candidateChildFrac = 2., maxNumOutgroups = 1,
               maxNumChromOutgroups = 1):
        # sort the (undirected) distance map
        orderedPairs = []
        for source, sinks in list(self.dm.items()):
            for sink, dist in list(sinks.items()):
                if source != self.root and sink != self.root:
                    orderedPairs.append((dist, (source, sink)))
        orderedPairs.sort(key = lambda x: x[0])
        finished = set()
        self.candidateMap = dict()
        if candidateSet is not None:
            assert isinstance(candidateSet, set)
            for candidate in candidateSet:
                self.candidateMap[candidate] = True

        htable = self.heightTable()

        # convert the input (leaf) chroms to id-space
        node_chroms = defaultdict(set)
        if self.chrom_map:
            for node_name, chr_list in self.chrom_map.items():
                node_chroms[self.mcTree.getNodeId(node_name)] = set(chr_list)
                
        # fill in the chromosomes of the ancestors, where they are (maybe naively)
        # just the union of all chroms below
        # these are the chroms that outgroups will need to have if possible.
        for node, chr_list in list(node_chroms.items()):
            parent = node
            while (self.mcTree.hasParent(parent)):
                parent = self.mcTree.getParent(parent)
                prev_len = len(node_chroms[parent])
                node_chroms[parent] = node_chroms[parent].union(chr_list)
                if len(node_chroms[parent]) == prev_len:
                    # no sense continuing traveral upwards if we're not adding anything
                    break

        # fast check to see if chromosomes of a given source are satisfied
        source_satisfied = {}
        for node in self.mcTree.postOrderTraversal():
            if node in node_chroms and len(node_chroms[node]) > 0 and not self.mcTree.isLeaf(node):
                # ancestral nodes with descendant chromosome specificaitons need to be satisfied
                source_satisfied[node] = False
            else:
                source_satisfied[node] = True

        # sort the orderedPairs by source
        ordered_pairs_by_source = defaultdict(list)
        for candidate in orderedPairs:
            # source is the ancestor we're trying to find outgroups for
            source = candidate[1][0]
            ordered_pairs_by_source[source].append(candidate)

        # visit the tree bottom up
        for node in self.mcTree.postOrderTraversal():
            # visit the candidates in order of increasing distance
            orderedPairs = ordered_pairs_by_source[node]
            for candidate in orderedPairs:
                # source is the ancestor we're trying to find outgroups for
                source = candidate[1][0]
                # sink it the candidate outgroup
                sink = candidate[1][1]
                sourceName = self.mcTree.getName(source)
                sinkName = self.mcTree.getName(sink)
                dist = candidate[0]

                # skip leaves (as sources)
                if len(self.dag.out_edges(source)) == 0:
                    finished.add(source)

                # skip nodes that were already finished in a previous run
                if sourceName in self.ogMap and source_satisfied[source] and len(self.ogMap[sourceName]) >= maxNumOutgroups:
                    finished.add(source)

                # skip invalid outgroups
                if sink in self.invalidSet:
                    continue

                # skip nodes that aren't in the candidate set (if specified)
                # or don't have enough candidate children
                if not self.inCandidateSet(sink, candidateChildFrac):
                    continue

                # candidate (ancestral) sink is too high up the tree compared to the source, so we skip
                if (threshold is not None and
                    not self.mcTree.isLeaf(sink) and
                    # "+ 1" because we want a threshold of 0 to indicate we
                    # probably won't have to wait for the outgroup (it will
                    # already have been finished, assuming all subproblems
                    # proceed at the same rate)
                    htable[sink] - htable[source] + 1 > threshold):
                    continue

                # Don't use any outgroups that are a child of another node
                # already in the outgroup set
                if any([self.onSamePath(x, sink) for x in self.dag.successors(source)]):
                    continue

                if source not in finished and \
                not self.onSamePath(source, sink):
                    self.dag.add_edge(source, sink, weight=dist, info='outgroup')
                    if NX.is_directed_acyclic_graph(self.dag):
                        htable[source] = max(htable[source], htable[sink] + 1)
                        existingOutgroups = [i[0] for i in self.ogMap[sourceName]]
                        if sinkName in existingOutgroups:
                            # This outgroup was already assigned to this source in a previous run
                            # Sanity check that the distance is equal
                            existingOutgroupDist = dict(self.ogMap[sourceName])
                            assert existingOutgroupDist[sinkName] == dist
                            continue
                        self.ogMap[sourceName].append((sinkName, dist))
                        source_satisfied[source] = self.check_chrom_satisfied(source, node_chroms)
                        if len(self.ogMap[sourceName]) >= maxNumOutgroups and source_satisfied[source]:
                            finished.add(source)
                    else:
                        self.dag.remove_edge(source, sink)

        # Since we could be adding to the ogMap instead of creating
        # it, sort the outgroups by distance again. Sorting the
        # outgroups is critical for the multiple-outgroups code to
        # work well.
        for node, outgroups in list(self.ogMap.items()):
            self.ogMap[node] = sorted(outgroups, key=lambda x: x[1])

        # the chromosome specification logic trumps the maximum number of outgroups
        # so we do a second pass to reconcile them (greedily) as best as possible
        if node_chroms:
            self.refine_og_chroms(node_chroms, max(maxNumOutgroups, maxNumChromOutgroups))


def main():
    usage = "usage: %prog <seqfile> <output graphviz .dot file>"
    description = "TEST: draw the outgroup DAG"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--justLeaves", dest="justLeaves", action="store_true",
                      default = False, help="Assign only leaves as outgroups")
    parser.add_option("--threshold", dest="threshold", type='int',
                      default = None, help="greedy threshold")
    parser.add_option("--numOutgroups", dest="maxNumOutgroups",
                      help="Maximum number of outgroups to provide", type=int)
    parser.add_option("--chromInfo", dest="chromInfo",
                      help="File mapping genomes to sex chromosome lists")
    parser.add_option("--configFile", dest="configFile",
                      help="Specify cactus configuration file",
                      default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))    
    options, args = parser.parse_args()

    options.binariesMode = 'local'
    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    # load up the seqfile and figure out the outgroups and schedule
    # can't use progressive_decomposition.py for circular deps...
    config_node = ET.parse(options.configFile).getroot()
    config_wrapper = ConfigWrapper(config_node)
    config_wrapper.substituteAllPredefinedConstantsWithLiterals(options)
    seq_file = SeqFile(args[0])
    mc_tree = MultiCactusTree(seq_file.tree)
    mc_tree.nameUnlabeledInternalNodes(config_wrapper.getDefaultInternalNodePrefix())
    mc_tree.computeSubtreeRoots()

    outgroup = GreedyOutgroup()
    outgroup.importTree(mc_tree)
    if options.chromInfo:
        outgroup.loadChromInfo(options.chromInfo)
    if options.justLeaves:
        candidates = set([mc_tree.getName(x)
                        for x in mc_tree.getLeaves()])
    else:
        candidates = None
    outgroup.greedy(threshold=options.threshold, candidateSet=candidates,
                    candidateChildFrac=1.1,
                    maxNumOutgroups=options.maxNumOutgroups)

    try:
        NX.drawing.nx_agraph.write_dot(outgroup.dag, args[1])
    except Exception as e:
        print(("NetworkX failed: %s" % str(e)))
        print("Writing ogMap in non-graphviz format")
        print("Tree %s" % NXNewick().writeString(mc_tree))
        with open(args[1], "w") as f:
            for node, ogs in list(outgroup.ogMap.items()):
                f.write("%s -> %s\n" % (node, str(ogs)))

    return 0

if __name__ == '__main__':
    main()
