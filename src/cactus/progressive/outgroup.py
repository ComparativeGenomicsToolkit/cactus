#!/usr/bin/env python3

#Copyright (C) 2011 by Glenn Hickey
#
#Released under the MIT license, see LICENSE.txt

""" Compute outgroups for subtrees by greedily finding
the nearest valid candidate.  has option to only assign
leaves as outgroups

"""

import os
import math
import copy
import itertools
import networkx as NX
from collections import defaultdict, namedtuple
from optparse import OptionParser
from cactus.shared.common import cactusRootPath
from cactus.shared.configWrapper import ConfigWrapper
from cactus.progressive.multiCactusTree import MultiCactusTree
from cactus.progressive.seqFile import SeqFile
from cactus.shared.common import cactus_call

class GreedyOutgroup(object):
    def __init__(self):
        self.dag = None
        self.dm = None
        self.dmDirected = None
        self.root = None
        self.ogMap = None
        self.mcTree = None

    # add edges from sonlib tree to self.dag
    # compute self.dm: an undirected distance matrix
    def importTree(self, mcTree, rootId = None):
        self.mcTree = mcTree
        self.dag = mcTree.nxDg.copy()
        self.root = mcTree.rootId
        self.stripNonEvents(self.root, mcTree.subtreeRoots)
        self.dmDirected = dict(NX.algorithms.shortest_paths.weighted.\
        all_pairs_dijkstra_path_length(self.dag))
        self.invalidSet = self.getInvalid(rootId)
        graph = NX.Graph(self.dag)
        self.dm = dict(NX.algorithms.shortest_paths.weighted.\
        all_pairs_dijkstra_path_length(graph))
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
    def greedy(self, threshold = None, candidateSet = None,
               candidateChildFrac = 2., maxNumOutgroups = 1):
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

        for candidate in orderedPairs:
            source = candidate[1][0]
            sink = candidate[1][1]
            sourceName = self.mcTree.getName(source)
            sinkName = self.mcTree.getName(sink)
            dist = candidate[0]

            # skip leaves (as sources)
            if len(self.dag.out_edges(source)) == 0:
                finished.add(source)

            # skip nodes that were already finished in a previous run
            if sourceName in self.ogMap and len(self.ogMap[sourceName]) >= maxNumOutgroups:
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
                    if len(self.ogMap[sourceName]) >= maxNumOutgroups:
                        finished.add(source)
                else:
                    self.dag.remove_edge(source, sink)

        # Since we could be adding to the ogMap instead of creating
        # it, sort the outgroups by distance again. Sorting the
        # outgroups is critical for the multiple-outgroups code to
        # work well.
        for node, outgroups in list(self.ogMap.items()):
            self.ogMap[node] = sorted(outgroups, key=lambda x: x[1])


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
    parser.add_argument("--configFile", dest="configFile",
                        help="Specify cactus configuration file",
                        default=os.path.join(cactusRootPath(), "cactus_progressive_config.xml"))    
    options, args = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    # load up the seqfile and figure out the outgroups and schedule
    # can't use progressive_decomposition.py for circular deps...
    config_node = ET.parse(options.configFile).getroot()
    config_wrapper = ConfigWrapper(config_node)
    config_wrapper.substituteAllPredefinedConstantsWithLiterals()
    seq_file = SeqFile(args[0])
    mc_tree = MultiCactusTree(seq_file.tree)
    mc_tree.nameUnlabeledInternalNodes(config_wrapper.getDefaultInternalNodePrefix())
    mc_tree.computeSubtreeRoots()

    outgroup = GreedyOutgroup()
    outgroup.importTree(mc_tree)
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
        with open(args[1], "w") as f:
            for node, ogs in list(outgroup.ogMap.items()):
                f.write("%s -> %s\n" % (node, str(ogs)))

    return 0

if __name__ == '__main__':
    main()
