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

from cactus.progressive.multiCactusProject import MultiCactusProject

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
            print((self.mcTree.getName(node), self.mcTree.isLeaf(node), node in good))
            if not self.mcTree.isLeaf(node) and node not in good:
                invalid.append(node)
        print([self.mcTree.getName(i) for i in invalid])
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


# First stab at better outgroup selection.  Uses estimated fraction of
# orthologous bases between two genomes, along with dynamic programming
# to select outgroups to best create ancestors.
#
# Only works with leaves for now (ie will never choose ancestor as outgroup)
#
class DynamicOutgroup(GreedyOutgroup):
    def __init__(self):
        self.SeqInfo = namedtuple("SeqInfo", "count totalLen umLen n50 umN50")
        self.sequenceInfo = None
        self.numOG = 1
        assert self.numOG is not None
        self.defaultBranchLength = 0.05
        # following weights control relative contributions of the three types
        # of "orthology loss" when computing branch scores
        self.lossFac = 1.
        self.fragFac = 1.
        self.mutFac = 1.
        # distance to sequence endpoint where we consider bases unalignable due
        # to fragmentation.
        self.edgeLen = 100

    # create map of leaf id -> sequence stats by scanninf the FASTA
    # files.  will be used to determine assembly quality for each input
    # genome (in a very crude manner, at least to start).
    #
    # for internal nodes, we store the stats of the max leaf underneath
    def importTree(self, mcTree, seqMap, rootId = None, candidateSet = None,
                   candidateBoost = 1.5):
        super(DynamicOutgroup, self).importTree(mcTree, rootId)
        self.candidateSet = candidateSet
        if candidateSet is not None and len(candidateSet) == 0:
            self.candidateSet = None
        self.candidateBoost = candidateBoost
        assert seqMap is not None
        # map name to (numSequences, totalLength)
        self.sequenceInfo = dict()
        for event, inPath in list(seqMap.items()):
            node = self.mcTree.getNodeId(event)
            if os.path.isdir(inPath):
                fastaPaths = [os.path.join(inPath, f) for
                              f in os.listdir(inPath)]
            else:
                fastaPaths = [inPath]
            totalFaInfo = self.__getSeqInfo(fastaPaths, event)
            self.sequenceInfo[node] = totalFaInfo

            # propagate leaf stats up to the root
            # can speed this up by O(N) but not sure if necessary..
            # we are conservative here in that we assume that the
            # ancestor has the longest, least fragmented genome possible
            # when judging from its descendants.
            x = node
            while self.mcTree.hasParent(x):
                x = self.mcTree.getParent(x)
                if x not in self.sequenceInfo:
                    self.sequenceInfo[x] = totalFaInfo
                else:
                    self.sequenceInfo[x] = self.SeqInfo(
                        min(totalFaInfo.count, self.sequenceInfo[x].count),
                        max(totalFaInfo.totalLen, self.sequenceInfo[x].totalLen),
                        max(totalFaInfo.umLen, self.sequenceInfo[x].umLen),
                        max(totalFaInfo.n50, self.sequenceInfo[x].n50),
                        max(totalFaInfo.umN50, self.sequenceInfo[x].umN50))

            #for node, info in self.sequenceInfo.items():
            #    print self.mcTree.getName(node), info

    # run the dynamic programming algorithm on each internal node
    def compute(self, maxNumOutgroups,
                mutationWeight = 1,
                sequenceLossWeight = .5):
        self.mutFac = mutationWeight
        self.lossFac = sequenceLossWeight
        self.numOG = maxNumOutgroups
        self.ogMap = dict()

        for node in self.mcTree.breadthFirstTraversal(self.root):
            if self.mcTree.isLeaf(node) or not self.mcTree.hasParent(node):
                continue
            self.__dpInit(node)
            self.__dpRun(node)
            nodeName = self.mcTree.getName(node)
            bestK = 0
            # we look for highest k with non-zero solution.
            # (can swap >= 0.0 with bestScore below to get the global best
            # not sure we'd want fewer outgroups..)
            for i in range(self.numOG + 1):
                if self.dpTable[node][i].score > 0.0:
                    bestK = i

            # we rank solution based on individual conservation score
            # of each outgroup vis-a-vis the target
            #rankFn = lambda x : 1. - self.__computeBranchConservation(
            #    x, self.dpTree.getRootId())
            # scratch that, we just use distance:
            rankFn = self.__getOgDist
            rankedSolution = sorted(self.dpTable[node][bestK].solution,
                                    key = rankFn)
            # convert to EventName,Dist format.  Note that distance
            # here is not necessarily what we're ranking on, and we include
            # it for consistency only.
            self.ogMap[nodeName] = [(self.dpTree.getName(x),
                                     self.__getOgDist(x))
                                     for x in rankedSolution]
            for og, dist in self.ogMap[nodeName]:
                self.dag.add_edge(node, self.mcTree.getNodeId(og),
                                  weight=dist, info="outgroup")
                #print self.dpTree.getName(node), "-->", og

    # initialize dynamic programming table
    def __dpInit(self, ancestralNodeId):
        self.dpTree = copy.deepcopy(self.mcTree)
        self.rootSeqInfo = None
        self.branchProbs = dict()
        self.DPEntry = namedtuple("DPEntry", "score solution")
        self.leafDistances = dict()
        # map .node id to [(score, solution)]
        # where list is for 0, 1, 2, ... k (ie best score for solution
        # of size k)
        self.dpTable = dict()

        # make a new tree rooted at the target ancestor with everything
        # below it, ie invalid outgroups, cut out
        for child in self.dpTree.getChildren(ancestralNodeId):
            self.dpTree.removeEdge(ancestralNodeId, child)
        self.dpTree.reroot(ancestralNodeId)
        self.rootSeqInfo = self.sequenceInfo[self.dpTree.getRootId()]

        # compute all the branch conservation probabilities
        for node in self.dpTree.preOrderTraversal():
            if self.dpTree.hasParent(node):
                self.branchProbs[node] = self.__computeBranchConservation(node)

        # set table to 0
        for node in self.dpTree.preOrderTraversal():
            self.dpTable[node] = []
            for i in range(self.numOG + 1):
                self.dpTable[node].append(self.DPEntry(0.0, []))

    # compute score for given node from its children using the dynamic
    # programming table
    def __dpNode(self, node):
        children = self.dpTree.getChildren(node)
        numChildren = len(children)
        # special case for leaf
        if numChildren == 0:
            self.dpTable[node][1] = self.DPEntry(1.0, [node])
        else:
            # iterate all possible combinations of child solution sizes
            # (very inefficeint since we only want unique solutions with
            # sum <= numOG, but assume numbers are small enough so doesn't
            # matter for now)
            cset = [x for x in range(0, self.numOG + 1)]
            if math.pow(len(cset), numChildren) > 1e6:
                raise RuntimeError("Dynamic Outgroup selection error for"
                                   " %s: degree limit exceeded.  Need to fix "
                                   "this!!" % self.dpTree.getName(node))
            for scoreAlloc in itertools.product(*[cset] * numChildren):
                csetK = sum(scoreAlloc)
                if  csetK > self.numOG:
                    continue
                # we compute the probability that a base is lost along
                # all the branches (so will be a product of of complement
                # of conservations along each branch)
                lossProb = 1.0
                solution = []
                for childNo, childId in enumerate(children):
                    childK = scoreAlloc[childNo]
                    childCons = self.dpTable[childId][childK].score
                    lossProb *= (1. - self.branchProbs[childId] * childCons)
                    solution += self.dpTable[childId][childK].solution
                # overall conservation is 1 - loss
                consProb = 1. - lossProb
                assert consProb >= 0. and consProb <= 1.
                assert len(solution) <= csetK
                if consProb > self.dpTable[node][csetK].score and \
                  len(solution) == csetK:
                    self.dpTable[node][csetK] = self.DPEntry(consProb, solution)

    # get the dynamic programming solution (for a single ancestor set in
    # __dpInit...)
    def __dpRun(self, node):
        for child in self.dpTree.getChildren(node):
            self.__dpRun(child)
        self.__dpNode(node)

    # compute the probability that a base is not "lost" on a branch
    # from given node to its parent
    # ancestor parameter allows us to consider a path from the node
    # to any ancestor in the tree as a single "branch" for purposes
    # of computation.  if none, then the immediate parent is used.
    def __computeBranchConservation(self, node, ancestor=None):
        if ancestor is None:
            ancestor = self.dpTree.getParent(node)
        nodeInfo = self.sequenceInfo[node]
        ancInfo = self.sequenceInfo[ancestor]

        # Loss probablity models alignment lost due to assembly quality.  We
        # use proportion of N50 (minus Ns) as crude proxy
        if nodeInfo.umN50 >= ancInfo.umN50:
            pLoss = 0.
        else:
            pLoss = 1. - (float(nodeInfo.umN50) / float(ancInfo.umN50))
        pLoss *= self.lossFac

        # Mutation probability is proportional to branch length.  We use
        # Jukes-Cantor model
        branchLength = 0.
        x = node
        while x != ancestor:
            weight = self.dpTree.getWeight(self.dpTree.getParent(x), x, None)
            if weight is None or weight < 0 or weight >= 1:
                # some kind of warning should happen here
                weight = self.defaultBranchLength
            branchLength += weight
            x = self.dpTree.getParent(x)
            assert x is not None
        jcMutProb = .75 - .75 * math.exp(-branchLength)
        jcMutProb *= self.mutFac

        conservationProb = (1. - pLoss) *  (1. - jcMutProb)
        assert conservationProb >= 0. and conservationProb <= 1.
        return conservationProb

    def __getOgDist(self, node):
        dist = 0.
        x = node
        while x != self.dpTree.getRootId():
            dist += self.dpTree.getWeight(self.dpTree.getParent(x), x,
                                          self.defaultBranchLength)
            x = self.dpTree.getParent(x)
            assert x != None
        return dist

    # use cactus_analyseAssembly to get some very basic stats about the
    # length and fragmentation of an assembly.  there is certainly
    # room for investigation of more sophisticated stats...
    #
    # Warning: we are tightly coupled to the output format of
    # cactus_analyseAssembly here... any change may cause an
    # assertion error
    def __getSeqInfo(self, faPaths, event):
        cmdLine = ["cactus_analyseAssembly"]
        for faPath in faPaths:
            if not os.path.isfile(faPath):
                raise RuntimeError("Unable to open sequence file %s" % faPath)
            cmdLine += [faPath]
        isCandidate = False
        if self.candidateSet is not None and event in self.candidateSet:
            isCandidate = True
        analyseOutput = cactus_call(parameters=cmdLine, work_dir=os.path.dirname(faPath), check_output=True)
        tsIdx = analyseOutput.index("Total-sequences:")
        assert tsIdx >= 0 and tsIdx < len(analyseOutput) - 1
        numSequences = int(analyseOutput[tsIdx + len("Total-sequences:") + 1])
        tlIdx = analyseOutput.index("Total-length:")
        assert tlIdx >= 0 and tlIdx < len(analyseOutput) - 1
        totalLength = int(analyseOutput[tlIdx + len("Total-length:") + 1])
        nsIdx = analyseOutput.index("ProportionNs:")
        assert nsIdx >= 0 and nsIdx < len(analyseOutput) - 1
        nsPct = float(analyseOutput[nsIdx + len("ProportionNs:") + 1])
        rmIdx = analyseOutput.index("Proportion-repeat-masked:")
        assert rmIdx >= 0 and rmIdx < len(analyseOutput) - 1
        rmPct = float(analyseOutput[rmIdx + len("Proportion-repeat-masked:") + 1])
        assert rmPct <= 1. and rmPct >= 0.
        n50Idx = analyseOutput.index("N50:")
        assert n50Idx >= 0 and n50Idx < len(analyseOutput) - 1
        n50 = int(analyseOutput[n50Idx + 1])

        if isCandidate is True:
            totalLength *= self.candidateBoost
            n50 *= self.candidateBoost

        umLength = max(0, totalLength * (1. - nsPct))
        umN50 = max(0, n50 * (1. - nsPct))

        return self.SeqInfo(numSequences, totalLength, umLength, n50, umN50)


def main():
    usage = "usage: %prog <project> <output graphviz .dot file>"
    description = "TEST: draw the outgroup DAG"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("--justLeaves", dest="justLeaves", action="store_true",
                      default = False, help="Assign only leaves as outgroups")
    parser.add_option("--threshold", dest="threshold", type='int',
                      default = None, help="greedy threshold")
    parser.add_option("--numOutgroups", dest="maxNumOutgroups",
                      help="Maximum number of outgroups to provide", type=int)
    parser.add_option("--dynamic", help="Use new dynamic programming"
                      " algorithm", action="store_true", default=False)
    options, args = parser.parse_args()

    if len(args) != 2:
        parser.print_help()
        raise RuntimeError("Wrong number of arguments")

    proj = MultiCactusProject()
    proj.readXML(args[0])
    if not options.dynamic:
        outgroup = GreedyOutgroup()
        outgroup.importTree(proj.mcTree)
        if options.justLeaves:
            candidates = set([proj.mcTree.getName(x)
                            for x in proj.mcTree.getLeaves()])
        else:
            candidates = None
        outgroup.greedy(threshold=options.threshold, candidateSet=candidates,
                        candidateChildFrac=1.1,
                        maxNumOutgroups=options.maxNumOutgroups)
    else:
        outgroup = DynamicOutgroup()
        outgroup.importTree(proj.mcTree, proj.getInputSequenceMap())
        outgroup.compute(options.maxNumOutgroups)

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
