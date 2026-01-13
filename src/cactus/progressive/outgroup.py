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
        self.lcaMap = None
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
        # LCA map: (node1, node2) -> lowest common ancestor
        self.lcaMap = dict(NX.tree_all_pairs_lowest_common_ancestor(self.dag, root=self.root))
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

    def computeNovelty(self, source, candidate, existing_outgroups):
        """Compute the novelty fraction of a candidate outgroup relative to existing outgroups.

        Novelty measures how much of the candidate's path to source is NOT shared with
        existing outgroups' paths. Returns a value between 0 and 1.

        For source X, existing outgroup A, and candidate B:
        - novel_length = (dist(X,B) - dist(X,A) + dist(A,B)) / 2
        - novelty_fraction = novel_length / dist(X,B)

        If there are multiple existing outgroups, returns the minimum novelty
        (candidate must be novel relative to ALL existing outgroups).
        """
        if not existing_outgroups:
            return 1.0  # First outgroup is always 100% novel

        dist_source_candidate = self.dm[source].get(candidate, float('inf'))
        if dist_source_candidate == 0 or dist_source_candidate == float('inf'):
            return 0.0

        min_novelty = 1.0
        for existing_og, _ in existing_outgroups:
            existing_id = self.mcTree.getNodeId(existing_og)
            dist_source_existing = self.dm[source].get(existing_id, float('inf'))
            dist_existing_candidate = self.dm[existing_id].get(candidate, float('inf'))

            if dist_existing_candidate == float('inf') or dist_source_existing == float('inf'):
                continue

            # novel_length = (dist(X,B) - dist(X,A) + dist(A,B)) / 2
            novel_length = (dist_source_candidate - dist_source_existing + dist_existing_candidate) / 2.0
            novelty_fraction = novel_length / dist_source_candidate

            # Clamp to [0, 1] in case of numerical issues
            novelty_fraction = max(0.0, min(1.0, novelty_fraction))
            min_novelty = min(min_novelty, novelty_fraction)

        return min_novelty

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

        nodes_in_tree = set()
        if self.mcTree:
            for node_id in self.mcTree.breadthFirstTraversal():
                nodes_in_tree.add(self.mcTree.getName(node_id))            

        with open(chrom_info_name, 'r') as chrom_info_file:
            for line in chrom_info_file:
                toks = line.rstrip().split()
                if len(toks) <= 2:
                    genome_name = toks[0]
                    if genome_name in self.chrom_map:
                        RuntimeError('Duplicate genome, {}, found in chromInfo file {}'.format(genome, chrom_info_name))
                    if nodes_in_tree and genome_name not in nodes_in_tree:
                        RuntimeError('Genome name, {}, from chromInfo file {}, not found in tree'.format(genome, chrom_info_name))
                    if len(toks) > 1:
                        chroms = toks[1].split(',')
                    else:
                        chroms = []
                    self.chrom_map[genome_name] = chroms
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

    def refine_og_chroms(self, node_chroms, max_outgroups, extra_chrom_outgroups):
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

                # cut off the outgroups list according to our input cap (if it's not set to auto with -1)
                if extra_chrom_outgroups >= 0 and len(essential_ogs) > max_outgroups + extra_chrom_outgroups:
                    cutoff = max_outgroups + extra_chrom_outgroups
                    logger.warning("Warning: Limiting outgroups for {} to {}, which means leaving out {}".format(
                        source_name, cutoff, ",".join([dog[0] for dog in essential_ogs[cutoff:]])))
                    # todo: rank by how many chromosomes they have!
                    essential_ogs = essential_ogs[:cutoff]
                
                if len(chrom_set) < len(source_chroms):
                    logger.warning("Warning: Unable to fulfull chromosome requirements when selecting outgroups for {}".format(source_name) +
                                   ". In particular, these chromosomes could not be found: {}".format(",".join(source_chroms - chrom_set)))
                # this is the amount we add no matter what, since we want to get to max_outgroups 
                required_padding = max(0, max_outgroups - len(essential_ogs))
                # this is the amount we add if the nearest-by-distance outgroups aren't included.
                # for example if we're allowed 1 extra outgroup and our max_outgroups is 3 and
                # we've gotten 2 essenital outgroups to satisfy the chromosomes, then we can
                # have optional_padding of 1.  And we'd add to it if either of the top-2-by distance
                # outgroups aren't in the set
                extra_cutoff = extra_chrom_outgroups if extra_chrom_outgroups >= 0 else len(essential_ogs)                
                optional_padding = required_padding + extra_cutoff
                
                if optional_padding > 0:
                    padded_ogs = []
                    for og_name, og_dist in og_names_dists:
                        is_essential = (og_name, og_dist) in essential_ogs
                        # always add essential
                        to_add = is_essential
                        # pad out as required
                        if not is_essential and required_padding > 0:
                            to_add = True
                            required_padding -= 1
                        # pad out if we have room and top greedy choices aren't essential
                        elif not is_essential and optional_padding > 0:
                            to_add = True
                        # tick down optional padding no matter what, as we only want it to apply to top choices
                        optional_padding -= 1
                        if to_add:
                            padded_ogs.append((og_name, og_dist))
                            required_padding 
                    essential_ogs = padded_ogs

                # update the outgroups
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
    # extraChromOutgroups : number of extra outgroups that can be added to attemp
    # to satisfy chromosomes.
    # topological : if True, sort outgroup candidates by (LCA_height, distance)
    # instead of just distance. This prefers topologically closer outgroups
    # (those sharing a more recent common ancestor) over those that happen
    # to have short branch lengths.
    # minNovelty : minimum novelty fraction (0-1) for additional outgroups.
    # Novelty measures how much of an outgroup's path is NOT shared with
    # existing outgroups. If no candidate meets minNovelty, falls back to
    # best available candidate (soft constraint).
    def greedy(self, threshold = None, candidateSet = None,
               candidateChildFrac = 2., maxNumOutgroups = 1,
               extraChromOutgroups = -1, topological = False,
               minNovelty = 0.0):
        # compute height table early (needed for topological sorting and threshold)
        htable = self.heightTable()

        # sort the (undirected) distance map
        orderedPairs = []
        for source, sinks in list(self.dm.items()):
            for sink, dist in list(sinks.items()):
                if source != self.root and sink != self.root:
                    orderedPairs.append((dist, (source, sink)))

        if topological:
            # Sort by (LCA_height, distance) tuple to prefer topologically
            # closer outgroups over those with short branch lengths.
            # Lower LCA height means a more recent common ancestor.
            def topo_key(x):
                dist, (source, sink) = x
                # lcaMap may have either (source, sink) or (sink, source) as key
                lca = self.lcaMap.get((source, sink), self.lcaMap.get((sink, source)))
                return (htable[lca], dist)
            orderedPairs.sort(key=topo_key)
        else:
            orderedPairs.sort(key = lambda x: x[0])

        finished = set()
        self.candidateMap = dict()
        if candidateSet is not None:
            assert isinstance(candidateSet, set)
            for candidate in candidateSet:
                self.candidateMap[candidate] = True

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

        # track fallback candidates when novelty requirement not met
        # key: source node id, value: list of (sink, sinkName, dist) tuples
        novelty_fallbacks = defaultdict(list)

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

                        # Check novelty requirement for additional outgroups
                        if minNovelty > 0 and len(self.ogMap[sourceName]) > 0:
                            novelty = self.computeNovelty(source, sink, self.ogMap[sourceName])
                            if novelty < minNovelty:
                                # Store as fallback and undo the edge addition
                                novelty_fallbacks[source].append((sink, sinkName, dist))
                                self.dag.remove_edge(source, sink)
                                htable[source] = max([htable[i] for i in self.dag.successors(source)] + [0]) + 1 if self.dag.out_edges(source) else htable[source]
                                continue

                        self.ogMap[sourceName].append((sinkName, dist))
                        source_satisfied[source] = self.check_chrom_satisfied(source, node_chroms)
                        if len(self.ogMap[sourceName]) >= maxNumOutgroups and source_satisfied[source]:
                            finished.add(source)
                    else:
                        self.dag.remove_edge(source, sink)

        # Fallback pass: for sources that didn't get enough outgroups due to
        # novelty requirement, use the best available fallbacks
        if minNovelty > 0:
            for source, fallbacks in novelty_fallbacks.items():
                sourceName = self.mcTree.getName(source)
                if len(self.ogMap[sourceName]) < maxNumOutgroups:
                    for sink, sinkName, dist in fallbacks:
                        if len(self.ogMap[sourceName]) >= maxNumOutgroups:
                            break
                        # Re-check if this fallback is still valid
                        if any([self.onSamePath(x, sink) for x in self.dag.successors(source)]):
                            continue
                        if sinkName in [og[0] for og in self.ogMap[sourceName]]:
                            continue
                        # Try to add the edge
                        self.dag.add_edge(source, sink, weight=dist, info='outgroup')
                        if NX.is_directed_acyclic_graph(self.dag):
                            self.ogMap[sourceName].append((sinkName, dist))
                            htable[source] = max(htable[source], htable[sink] + 1)
                        else:
                            self.dag.remove_edge(source, sink)

        # Since we could be adding to the ogMap instead of creating
        # it, sort the outgroups again. Sorting the outgroups is critical
        # for the multiple-outgroups code to work well.
        if topological:
            # Sort by (LCA_height, distance) to maintain topological preference
            def topo_sort_key(source_name, og_name, og_dist):
                source_id = self.mcTree.getNodeId(source_name)
                sink_id = self.mcTree.getNodeId(og_name)
                lca = self.lcaMap.get((source_id, sink_id), self.lcaMap.get((sink_id, source_id)))
                return (htable[lca], og_dist)
            for node, outgroups in list(self.ogMap.items()):
                self.ogMap[node] = sorted(outgroups, key=lambda x: topo_sort_key(node, x[0], x[1]))
        else:
            for node, outgroups in list(self.ogMap.items()):
                self.ogMap[node] = sorted(outgroups, key=lambda x: x[1])

        # the chromosome specification logic trumps the maximum number of outgroups
        # so we do a second pass to reconcile them (greedily) as best as possible
        if node_chroms:
            self.refine_og_chroms(node_chroms, maxNumOutgroups, extraChromOutgroups)


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
    parser.add_option("--extraChromOutgroups", dest="extraChromOutgroups",
                      help="Maximum number of outgroups to provide if necessitated by --chromInfo", type=int, default=-1)    
    parser.add_option("--chromInfo", dest="chromInfo",
                      help="File mapping genomes to sex chromosome lists")
    parser.add_option("--topological", dest="topological", action="store_true",
                      default=False, help="Sort outgroup candidates by (LCA_height, distance) "
                      "instead of just distance. This prefers topologically closer outgroups "
                      "(sharing a more recent common ancestor) over those with short branch lengths.")
    parser.add_option("--minOutgroupNovelty", dest="minNovelty", type='float',
                      default=0.0, help="Minimum novelty fraction (0-1) for additional outgroups. "
                      "Novelty measures how much of an outgroup's path is NOT shared with existing "
                      "outgroups. Higher values ensure more diverse/informative outgroups. "
                      "Falls back to best available if no candidate meets threshold. (default: 0)")
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
                    maxNumOutgroups=options.maxNumOutgroups,
                    extraChromOutgroups=options.extraChromOutgroups,
                    topological=options.topological,
                    minNovelty=options.minNovelty)

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
