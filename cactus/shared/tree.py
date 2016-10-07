#!/usr/bin/env python

#Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

import sys
import os
import re
import math
import random

from cactus.shared.misc import close
#import bioio

#########################################################
#########################################################
#########################################################
#basic tree datastructures
#########################################################
#########################################################
#########################################################

MIN_TREE_DISTANCE = 0.00001
 
class BinaryTree:
    def __init__(self, distance, internal, left, right, iD):
        self.distance = distance
        self.internal = internal
        self.left = left
        self.right = right
        self.iD = iD
        
class TraversalID:
    """
    tree traversal numbers, used as nodeIDs for identifying
    orders in the tree
    """
    def __init__(self, midStart, mid, midEnd):
        self.midStart = midStart
        self.mid = mid
        self.midEnd = midEnd

#########################################################
#########################################################
#########################################################
#tree functions
#########################################################
#########################################################
#########################################################

def binaryTree_depthFirstNumbers(binaryTree, labelTree=True, dontStopAtID=True):
    """
    get mid-order depth first tree numbers
    """
    traversalIDs = {}
    def traverse(binaryTree, mid=0, leafNo=0):
        if binaryTree.internal and (dontStopAtID or binaryTree.iD is None):
            midStart = mid
            j, leafNo = traverse(binaryTree.left, mid, leafNo)
            mid = j
            j, leafNo = traverse(binaryTree.right, j+1, leafNo)
            traversalIDs[binaryTree] = TraversalID(midStart, mid, j)
            return j, leafNo
        traversalID = TraversalID(mid, mid, mid+1)
        traversalID.leafNo = leafNo
        #thus nodes must be unique
        traversalIDs[binaryTree] = traversalID
        return mid+1, leafNo+1
    traverse(binaryTree)
    if labelTree:
        for binaryTree in traversalIDs.keys():
            binaryTree.traversalID = traversalIDs[binaryTree]
    return traversalIDs

def binaryTree_nodeNames(binaryTree):
    """
    creates names for the leave and internal nodes
    of the newick tree from the leaf labels
    """
    def fn(binaryTree, labels):
        if binaryTree.internal:
            fn(binaryTree.left, labels)
            fn(binaryTree.right, labels)
            labels[binaryTree.traversalID.mid] = labels[binaryTree.left.traversalID.mid] + "_" + labels[binaryTree.right.traversalID.mid]
            return labels[binaryTree.traversalID.mid]
        else:
            labels[binaryTree.traversalID.mid] = str(binaryTree.iD)
            return labels[binaryTree.traversalID.mid]
    labels = [None]*binaryTree.traversalID.midEnd
    fn(binaryTree, labels)
    return labels

def getBinaryTreeNodes(binaryTree, l):
    if binaryTree is not None:
        getBinaryTreeNodes(binaryTree.left, l)
        l.append(binaryTree)
        getBinaryTreeNodes(binaryTree.right, l)
        
def binaryTree_leafNo(binaryTree):
    return (binaryTree.traversalID.midEnd - binaryTree.traversalID.midStart + 1)/2

def getDistanceMatrix(tree):
    m = {}
    def fn(tree):
        if tree.internal:
            leftNodes = fn(tree.left)
            rightNodes = fn(tree.right)
            for i, d1 in leftNodes:
                for j, d2 in rightNodes:
                    m[(i, j)] = d1 + d2
                    m[(j, i)] = d1 + d2
            j = tree.traversalID.mid
            for i, d1 in leftNodes + rightNodes:
                m[(i, j)] = d1
                m[(j, i)] = d1
            return [ (i[0], i[1]+tree.distance) for i in leftNodes + rightNodes + [(j, 0.0)] ]
        return [ (tree.traversalID.mid, tree.distance) ]
    fn(tree)       
    return m


def makeRandomBinaryTree(leafNodeNumber=None):
    """Creates a random binary tree.
    """
    while True:
        nodeNo = [-1]
        def fn():
            nodeNo[0] += 1
            if random.random() > 0.6:
                i = str(nodeNo[0])
                return BinaryTree(0.00001 + random.random()*0.8, True, fn(), fn(), i)
            else:
                return BinaryTree(0.00001 + random.random()*0.8, False, None, None, str(nodeNo[0]))
        tree = fn()
        def fn2(tree):
            if tree.internal:
                return fn2(tree.left) + fn2(tree.right)
            return 1
        if leafNodeNumber is None or fn2(tree) == leafNodeNumber:
            return tree

def getRandomBinaryTreeLeafNode(binaryTree):
    """Get random binary tree node.
    """
    if binaryTree.internal == True:
        if random.random() > 0.5:
            return getRandomBinaryTreeLeafNode(binaryTree.left)
        else:
            return getRandomBinaryTreeLeafNode(binaryTree.right)
    else:
        return binaryTree

#########################################################
#########################################################
#########################################################
#substition functions and felsensteins algorithm
#########################################################
#########################################################
#########################################################

def transformByDistance(wV, subModel, alphabetSize=4):
    """
    transform wV by given substitution matrix
    """
    nc = [0.0]*alphabetSize
    for i in xrange(0, alphabetSize):
        j = wV[i]
        k = subModel[i]
        for l in xrange(0, alphabetSize):
            nc[l] += j * k[l]
    return nc

def multiplyWV(wVX, wVY, alphabetSize=4):
    return [ wVX[i] * wVY[i] for i in xrange(0, alphabetSize) ]

def sumWV(wVX, wVY, alphabetSize=4):
    return [ wVX[i] + wVY[i] for i in xrange(0, alphabetSize) ]
    
def normaliseWV(wV, normFac=1.0):
    """
    make char probs divisible by one
    """
    f = sum(wV) / normFac
    return [ i/f for i in wV ]

def sumWVA(wVA, alphabetSize=4):
    totals = [0.0]*alphabetSize
    for wV in wVA:
        for i in xrange(0, alphabetSize):
            totals[i] += wV[i]
    return totals

def felsensteins(binaryTree, subMatrices, ancestorProbs, leaves, alphabetSize):
    """
    calculates the un-normalised probabilties of each non-gap residue position
    """
    l = {}
    def upPass(binaryTree):
        if binaryTree.internal: #is internal binaryTree
            i = branchUp(binaryTree.left)
            j = branchUp(binaryTree.right)
            k = multiplyWV(i, j, alphabetSize)
            l[binaryTree.traversalID.mid] = (k, i, j)
            return k
        l[binaryTree.traversalID.mid] = leaves[binaryTree.traversalID.leafNo]
        return leaves[binaryTree.traversalID.leafNo]
    def downPass(binaryTree, ancestorProbs):
        if binaryTree.internal: #is internal binaryTree
            i = l[binaryTree.traversalID.mid]
            l[binaryTree.traversalID.mid] = multiplyWV(ancestorProbs, i[0], alphabetSize)
            branchDown(binaryTree.left, multiplyWV(ancestorProbs, i[2], alphabetSize))
            branchDown(binaryTree.right, multiplyWV(ancestorProbs, i[1], alphabetSize))
    def branchUp(binaryTree):
        return transformByDistance(upPass(binaryTree), subMatrices[binaryTree.traversalID.mid], alphabetSize)
    def branchDown(binaryTree, ancestorProbs):
        downPass(binaryTree, transformByDistance(ancestorProbs, subMatrices[binaryTree.traversalID.mid], alphabetSize))
    upPass(binaryTree)
    downPass(binaryTree, ancestorProbs)
    return l

def calculateCharacterFrequencies(seq, map, alphabetSize):
    counts = [0.0]*alphabetSize
    for i in seq:
        counts[map(i)] += 1
    return counts
    
#########################################################
#########################################################
#########################################################
#distance matrix to tree building functions
#########################################################
#########################################################
#########################################################

class DistancePair:
    def __init__(self, distance, leaf1, leafNo1, leaf2, leafNo2):
        self.distance = distance
        self.leaf1 = leaf1
        self.leaf2 = leaf2
        self.leafNo1 = leafNo1
        self.leafNo2 = leafNo2
    
    def __cmp__(self, distancePair):
        if self.distance < distancePair.distance:
            return -1
        if self.distance > distancePair.distance:
            return 1
        return 0 #don't care
        #doesn't wort for floats return self.distance.__cmp__(distancePair.distance)
        
def correctTreeDistances(tree):
    if tree is not None:
        if tree.distance < MIN_TREE_DISTANCE:
            tree.distance = MIN_TREE_DISTANCE
        correctTreeDistances(tree.left)
        correctTreeDistances(tree.right)
        
def calculateDNADistanceMatrix(seqNo, fastaIter, transitionTransversionRatio=2.0):
    transitions = [0.1]*seqNo*seqNo
    transversions = [0.1]*seqNo*seqNo
    counts = [1.0]*seqNo*seqNo
    for column in fastaIter:
        for i in xrange(0, seqNo):
            if column[i] in [ 'A', 'C', 'T', 'G' ]:
                for j in xrange(i+1, seqNo):
                    if column[j] in [ 'A', 'C', 'T', 'G' ]:
                        counts[i*seqNo + j] += 1
                        if column[i] != column[j]:
                            if column[i] in [ 'A', 'G' ]:
                                if column[j] in [ 'C', 'T' ]:
                                    transversions[i*seqNo + j] += 1
                                else:
                                    transitions[i*seqNo + j] += 1
                            else:
                                if column[j] in [ 'A', 'G' ]:
                                    transversions[i*seqNo + j] += 1
                                else:
                                    transitions[i*seqNo + j] += 1
    distanceMatrix = [ [None]*seqNo for i in xrange(0, seqNo) ]
    for i in xrange(0, seqNo*seqNo):
        for j in xrange(i+1, seqNo):
            k = i * seqNo + j
            distanceMatrix[i][j] = -0.75*math.log(1 - (4/3)*((transitions[k]+transversions[k])/counts[k])) #jukes cantor correction
            distanceMatrix[j][i] = distanceMatrix[i][j]
            #print "boo", i, j, distanceMatrix[i][j], (transitions[k]+transversions[k])/counts[k]
        #distanceMatrix[i] = -0.5*math.log(1 - 2*P - Q)-0.25*math.log(1 - 2*Q)
    return distanceMatrix

def makeDistancePairs(distanceMatrix, iDs, seqNo):
    binaryTrees = [ BinaryTree(0.0, False, None, None, iDs[i]) for i in xrange(0, seqNo) ]
    distancePairs = []
    for i in xrange(0, seqNo):
        for j in xrange(i+1, seqNo): 
            distancePairs.append(DistancePair(distanceMatrix[i][j], binaryTrees[i], 1, binaryTrees[j], 1))
            distancePairs.append(DistancePair(distanceMatrix[i][j], binaryTrees[j], 1, binaryTrees[i], 1))
    return distancePairs

def upgma(distanceMatrix, iDs, leafNo):
    binaryTree = upgmaI(makeDistancePairs(distanceMatrix, iDs, leafNo), leafNo)
    def fn(tree):
        if tree.internal:
            tree.distance -= tree.left.distance
            fn(tree.left)
            fn(tree.right)
    fn(binaryTree)
    binaryTree.distance = MIN_TREE_DISTANCE
    correctTreeDistances(binaryTree)
    return binaryTree

def upgmaI(distancePairs, leafNo):
    #get min pair
    distancePairs.sort()
    distancePair = distancePairs[0]
    #calculate shared distance
    distancePair.leaf1.distance = distancePair.distance/2
    distancePair.leaf2.distance = distancePair.distance/2
    newLeaf = BinaryTree(0.0, True, distancePair.leaf1, distancePair.leaf2, None)
    if leafNo-1 == 1:
        return newLeaf
    #replace references
    holder1 = {}
    holder2 = {}
    newDistances = []
    for i in distancePairs:
        if i.leaf1 == distancePair.leaf1 and i.leaf2 != distancePair.leaf2:
            holder1[i.leaf2] = i
        if i.leaf1 == distancePair.leaf2 and i.leaf2 != distancePair.leaf1:
            holder2[i.leaf2] = i
    assert len(holder1.keys()) == leafNo-2
    assert len(holder2.keys()) == leafNo-2
    assert set(holder1.keys()) == set(holder2.keys())
    for i in holder1.keys():
        j = holder1[i]
        k = holder2[i]
        newDistance = (j.distance*j.leafNo1 + k.distance*k.leafNo1)/(j.leafNo1 + k.leafNo1)
        newDistances.append(DistancePair(newDistance, j.leaf2, j.leafNo2, newLeaf, j.leafNo1 + k.leafNo1))
        newDistances.append(DistancePair(newDistance, newLeaf, j.leafNo1 + k.leafNo1, j.leaf2, j.leafNo2))
    distancePairs = [ i for i in distancePairs if (i.leaf1 != distancePair.leaf1 and i.leaf1 != distancePair.leaf2 and i.leaf2 != distancePair.leaf1 and i.leaf2 != distancePair.leaf2) ] + newDistances
    return upgmaI(distancePairs, leafNo-1)

def nj(distanceMatrix, iDs, leafNo):
    binaryTree = njI(makeDistancePairs(distanceMatrix, iDs, leafNo), leafNo)
    correctTreeDistances(binaryTree)
    return binaryTree

def getMinPair(distancePairs, rValues, leafNo):
    j = None
    k = sys.maxint
    for i in distancePairs:
        adjustD = i.distance - (rValues[i.leaf1] + rValues[i.leaf2])/(leafNo-2)
        #print "the adjusted value ", adjustD, i.distance, rValues[i.leaf1]/(leafNo-2), rValues[i.leaf2]/(leafNo-2)
        if adjustD < k:
            k = adjustD
            j = i
    #print "value is ", k, j.distance
    return j

def calculateRValues(distancePairs, leafNo):
    j = {}
    for i in distancePairs:
        if j.has_key(i.leaf1):
            j[i.leaf1] += i.distance
        else:
            j[i.leaf1] = i.distance
    assert len(j.keys()) == leafNo
    return j

def njI(distancePairs, leafNo):
    assert leafNo >= 2 
    if leafNo == 2:
        assert len(distancePairs) == 2
        distancePair = distancePairs[0]
        distancePair.leaf1.distance = distancePair.distance*0.5
        distancePair.leaf2.distance = distancePair.distance*0.5
        return BinaryTree(MIN_TREE_DISTANCE, True, distancePair.leaf1, distancePair.leaf2, None)
    #calculate r values
    rValues = calculateRValues(distancePairs, leafNo)
    #get min pair
    distancePair = getMinPair(distancePairs, rValues, leafNo)
    #distance, internal, left, right
    distancePair.leaf1.distance = 0.5*(distancePair.distance + (rValues[distancePair.leaf1] - rValues[distancePair.leaf2])/(leafNo-2))
    distancePair.leaf2.distance = distancePair.distance - distancePair.leaf1.distance
    newLeaf = BinaryTree(0.0, True, distancePair.leaf1, distancePair.leaf2, None)
    #replace references
    holder1 = {}
    holder2 = {}
    newDistances = []
    for i in distancePairs:
        if i.leaf1 == distancePair.leaf1 and i.leaf2 != distancePair.leaf2:
            holder1[i.leaf2] = i
        if i.leaf1 == distancePair.leaf2 and i.leaf2 != distancePair.leaf1:
            holder2[i.leaf2] = i
    assert len(holder1.keys()) == leafNo-2
    assert len(holder2.keys()) == leafNo-2
    assert set(holder1.keys()) == set(holder2.keys())
    for i in holder1.keys():
        j = holder1[i]
        k = holder2[i]
        assert j.leaf2 == k.leaf2
        #print "the leaf is ", j.leaf2
        newDistance = 0.5*(j.distance + k.distance - distancePair.distance)
        #print "now a new distance", newDistance
        newDistances.append(DistancePair(newDistance, j.leaf2, 0, newLeaf, 0)) #leaf numbers are un important, and hence omitted
        newDistances.append(DistancePair(newDistance, newLeaf, 0, j.leaf2, 0)) 
    distancePairs = [ i for i in distancePairs if (i.leaf1 != distancePair.leaf1 and i.leaf1 != distancePair.leaf2 and i.leaf2 != distancePair.leaf1 and i.leaf2 != distancePair.leaf2) ] + newDistances
    return njI(distancePairs, leafNo-1)

#########################################################
#########################################################
#########################################################
#substitution matrix functions
#########################################################
#########################################################
#########################################################

def checkMatrix(m, fV, AS=4, reversible=True):
    #print m
    for i in xrange(0, AS):
        j = sum(m[i])
        #print "AAAAA", j
        assert j <= 1.0001
        assert j >= 0.9999
        if reversible:
            for k in xrange(0, AS):
                #print "comp2", (fV[i] * m[i][k]), (fV[k] * m[k][i] )
                assert close(fV[i] * m[i][k], fV[k] * m[k][i], 0.00001)
    
    wV = fV
    wV2 = fV
    wV3 = transformByDistance(wV, m, AS)
    wV4 = transformByDistance(wV2, m, AS)
    i = sum(multiplyWV(wV2, wV3, AS))
    j = sum(multiplyWV(wV, wV4, AS))
    #print i, j
    assert close(i, j, 0.00001)
    
def reverseSubMatrix(m, AS=4):
    k = [ [None]*AS for i in xrange(0, AS) ]
    for i in xrange(0, AS):
        for j in xrange(0, AS):
            k[j][i] = m[i][j]
    return k
    
def subMatrix_jukesCantor(d):
    i = 0.25 + 0.75*math.exp(-(4.0/3.0)*d)
    j = 0.25 - 0.25*math.exp(-(4.0/3.0)*d)
    return [ [i, j, j, j], [j, i, j, j], [j, j, i, j], [j, j, j, i] ]

"""
def distanceTamureiNei(aF, cF, gF, tF, a2G, t2C, tV):
    rF = aF + gF
    yF = cF + tF
    
    xx = 1.0 - a2G * rF/(2.0 * aF * gF) - tV/(2.0 * rF)
    yy = 1.0 - t2C * yF/(2.0 * tF * cF) - tV/(2.0 * yF)
    i = -(2.0 * aF * gF / rF) * math.log(1.0 - a2G * rF/(2.0 * aF * gF) - tV/(2.0 * rF))
    j = -(2.0 * tF * cF / yF) * math.log(1.0 - t2C * yF/(2.0 * tF * cF) - tV/(2.0 * yF))
    k = -2.0 * (rF * yF - (aF * gF * yF / rF) - (tF * cF * rF / yF)) * math.log(1.0 - tV/(2.0 * rF * yF))
    print i, j, k, xx, yy
    d = -(2.0 * aF * gF / rF) * math.log(1.0 - a2G * rF/(2.0 * aF * gF) - tV/(2.0 * rF)) \
        -(2.0 * tF * cF / yF) * math.log(1.0 - t2C * yF/(2.0 * tF * cF) - tV/(2.0 * yF)) \
        -2.0 * (rF * yF - (aF * gF * yF / rF) - (tF * cF * rF / yF)) * math.log(1.0 - tV/(2.0 * rF * yF))
    return d
"""

def subMatrix_TamuraNei(d, fA, fC, fG, fT, alphaPur, alphaPyr, beta):
    i =  fA + fC + fG + fT
    assert i < 1.00001
    assert i > 0.99999
    assert d >= 0.0
    #assert alphaPur >= 0.0
    #assert alphaPyr >= 0.0
    #assert beta >= 0
    
    AS = 4
    freq = ( fA, fC, fG, fT )
    alpha = ( alphaPur, alphaPyr, alphaPur, alphaPyr )
    matrix = [ [ 0.0 ]*AS for i in xrange(0, AS) ]
    #see page 203 of Felsenstein's Inferring Phylogenies for explanation of calculations
    def watKro(j, k):
        if (j % 2) == (k % 2):
            return 1.0
        return 0.0
    def kroenickerDelta(i, j):
        if i == j:
            return 1.0
        return 0.0
    for i in xrange(0, AS): #long winded, totally unoptimised method for calculating matrix
        for j in xrange(0, AS):
            l = 0.0
            for k in xrange(0, AS):
                l += watKro(j, k) * freq[k]
            matrix[i][j] =\
            math.exp(-(alpha[i] + beta) * d) * kroenickerDelta(i, j) + \
            math.exp(-beta*d) * (1.0 - math.exp(-alpha[i]*d)) * (freq[j] * watKro(i, j) / l) + \
            (1.0 - math.exp(-beta * d)) * freq[j]
    checkMatrix(matrix, (fA, fC, fG, fT))
    return matrix

def subMatrix_HKY(d, fA, fC, fG, fT, transitionTransversionR):
    i =  fA + fC + fG + fT
    assert i < 1.00001
    assert i > 0.99999
    
    fPur = fA + fG
    fPyr = fC + fT
    p = fPur/fPyr #makes like HKY
    
    beta = 1.0 / (2.0 * fPur * fPyr * (1.0 + transitionTransversionR))
    alphaPyr = ((fPur * fPyr * transitionTransversionR) - (fA * fG) - (fC * fT)) \
                / (2.0 * (1.0 + transitionTransversionR) * (fPyr * fA * fG * p + fPur * fC * fT))
    alphaPur = p * alphaPyr
    return subMatrix_TamuraNei(d, fA, fC, fG, fT, alphaPur, alphaPyr, beta)

def subMatrix_HalpernBruno(d, freqColumn, subMatrix, AS=4):
    #return subMatrix_HKY(d, freqColumn[0], freqColumn[1], freqColumn[2], freqColumn[3], 2.0)
    #return subMatrix
    matrix = [ [ 0.0 ]*AS for i in xrange(0, AS) ]
    for i in xrange(0, AS):
        for j in xrange(0, AS):
            a = freqColumn[i] * subMatrix[i][j]
            b = freqColumn[j] * subMatrix[j][i]
            if not close(a, b, 0.0001):
                matrix[i][j] = subMatrix[i][j] * (math.log(b/a) / (1 - (a/b)))
            else:
                matrix[i][j] = subMatrix[i][j]
    #for i in xrange(0, AS):
    #    #print matrix[i][i], sum(matrix[i])
    #    matrix[i][i] -= sum(matrix[i]) - 1.0
    #    assert matrix[i][i] >= 0
    #checkMatrix(matrix, freqColumn)
    return matrix

#########################################################
#########################################################
#########################################################
#misc tree functions
#########################################################
#########################################################
#########################################################

def annotateTree(bT, fn):
    """
    annotate a tree in an external array using the given function
    """
    l = [None]*bT.traversalID.midEnd
    def fn2(bT):
        l[bT.traversalID.mid] = fn(bT)
        if bT.internal:
            fn2(bT.left)
            fn2(bT.right)
    fn2(bT)
    return l

def mapTraversalIDsBetweenTrees(oldTree, newTree):
    map = {}
    leafMap = {}
    internalMap = {}
    def fn(i):
        j = i.traversalID.mid
        if j == oldTree.traversalID.mid or (oldTree.internal and j == oldTree.left.traversalID.mid or j == oldTree.right.traversalID.mid):
            return oldTree.traversalID.mid
        return j
    def fn2(oldTree):
        if oldTree.internal:
            fn2(oldTree.left)
            fn2(oldTree.right)
        else:
            leafMap[oldTree.iD] = fn(oldTree)
    fn2(oldTree)
    def fn3(oldTree):
        if oldTree.internal:
            fn3(oldTree.left)
            fn3(oldTree.right)
            internalMap[(fn(oldTree.left), fn(oldTree.right))] = fn(oldTree)
            internalMap[(fn(oldTree.right), fn(oldTree.left))] = fn(oldTree)
            internalMap[(fn(oldTree), fn(oldTree.right))] = fn(oldTree.left)
            internalMap[(fn(oldTree.right), fn(oldTree))] = fn(oldTree.left)
            internalMap[(fn(oldTree), fn(oldTree.left))] = fn(oldTree.right)
            internalMap[(fn(oldTree.left), fn(oldTree))] = fn(oldTree.right)
    fn3(oldTree)
    print leafMap
    print internalMap
    def fn4(newTree):
        if newTree.internal:
            fn4(newTree.left)
            fn4(newTree.right)
            map[newTree.traversalID.mid] = internalMap[(map[newTree.left.traversalID.mid], map[newTree.right.traversalID.mid])]
        else:
            map[newTree.traversalID.mid] = leafMap[newTree.iD]
    fn4(newTree)
    return map    

def remodelTreeRemovingRoot(root, node):
    """
    Node is mid order number
    """
    import bioio
    assert root.traversalID.mid != node
    hash = {}
    def fn(bT):
        if bT.traversalID.mid == node:
            assert bT.internal == False
            return [ bT ]
        elif bT.internal:
            i = fn(bT.left)
            if i is None:
                i = fn(bT.right)
            if i is not None:
                hash[i[-1]]= bT
                i.append(bT)
            return  i
        return None
    l = fn(root)
    def fn2(i, j):
        if i.left == j:
            return i.right
        assert i.right == j
        return i.left
    def fn3(bT):
        if hash[bT] == root:
            s = '(' + bioio.printBinaryTree(fn2(hash[bT], bT), bT, True)[:-1] + ')'
        else:
            s = '(' + bioio.printBinaryTree(fn2(hash[bT], bT), bT, True)[:-1] + ',' + fn3(hash[bT]) + ')'
        return s + ":" + str(bT.distance)
    s = fn3(l[0]) + ';'
    t = bioio.newickTreeParser(s)
    return t

def moveRoot(root, branch):
    """
    Removes the old root and places the new root at the mid point along the given branch
    """
    import bioio
    if root.traversalID.mid == branch:
        return bioio.newickTreeParser(bioio.printBinaryTree(root, True))
    def fn2(tree, seq):
        if seq is not None:
            return '(' + bioio.printBinaryTree(tree, True)[:-1] + ',' + seq + ')'
        return bioio.printBinaryTree(tree, True)[:-1]
    def fn(tree, seq):
        if tree.traversalID.mid == branch:
            i = tree.distance
            tree.distance /= 2
            seq = '(' + bioio.printBinaryTree(tree, True)[:-1] + ',(' + seq + ('):%s' % tree.distance) + ');'
            tree.distance = i
            return seq
        if tree.internal:
            if branch < tree.traversalID.mid:
                seq = fn2(tree.right, seq)
                return fn(tree.left, seq)
            else:
                assert branch > tree.traversalID.mid
                seq = fn2(tree.left, seq)
                return fn(tree.right, seq)
        else:
            return bioio.printBinaryTree(tree, True)[:-1]
    s = fn(root, None)
    return bioio.newickTreeParser(s)

def checkGeneTreeMatchesSpeciesTree(speciesTree, geneTree, processID):
    """
    Function to check ids in gene tree all match nodes in species tree
    """
    def fn(tree, l):
        if tree.internal:
            fn(tree.left, l)
            fn(tree.right, l)
        else:
            l.append(processID(tree.iD))
    l = []
    fn(speciesTree, l)
    l2 = []
    fn(geneTree, l2)
    for i in l2:
        #print "node", i, l
        assert i in l

def calculateDupsAndLossesByReconcilingTrees(speciesTree, geneTree, processID):
    """
    Reconciles the given gene tree with the species tree and
    report the number of needed duplications and losses
    """  
    checkGeneTreeMatchesSpeciesTree(speciesTree, geneTree, processID)   
    def fn(tree, m):  
        if tree.internal:
            nodes = fn(tree.left, m)
            nodes = nodes.union(fn(tree.right, m))
            m[tree.traversalID.mid] = nodes
        else:
            m[tree.traversalID.mid] = set((processID(tree.iD),))
        return m[tree.traversalID.mid]
    a = {}
    fn(speciesTree, a)
    b = {}
    fn(geneTree, b)
    def fn2(nodes, speciesTree):
        assert nodes.issubset(a[speciesTree.traversalID.mid])
        if speciesTree.internal:
            if nodes.issubset(a[speciesTree.left.traversalID.mid]):
                return fn2(nodes, speciesTree.left)
            if nodes.issubset(a[speciesTree.right.traversalID.mid]):
                return fn2(nodes, speciesTree.right)
        return speciesTree.traversalID.mid
    for iD in b.keys():
        nodes = b[iD]
        b[iD] = fn2(nodes, speciesTree)
    dups = []
    def fn3(geneTree):
        if geneTree.internal:
            i = b[geneTree.traversalID.mid]
            if b[geneTree.left.traversalID.mid] == i or b[geneTree.right.traversalID.mid] == i:
                dups.append(geneTree.traversalID.mid)
            fn3(geneTree.left)
            fn3(geneTree.right)
    fn3(geneTree)
    lossMap = {}
    def fn4(speciesTree):
        nodes = [(speciesTree.traversalID.mid, -1)]
        lossMap[(speciesTree.traversalID.mid, speciesTree.traversalID.mid)] = 0
        if speciesTree.internal:
            for node, losses in fn4(speciesTree.left) + fn4(speciesTree.right):
                lossMap[(speciesTree.traversalID.mid, node)] = losses+1
                nodes.append((node, losses+1))
        return nodes
    for node, losses in fn4(speciesTree):
        lossMap[(sys.maxint, node)] = losses+1
    losses = [0]
    def fn5(geneTree, ancestor):
        if geneTree.internal:
            i = b[geneTree.traversalID.mid]
            if geneTree.traversalID.mid in dups:
                losses[0] += lossMap[(ancestor, b[geneTree.left.traversalID.mid])]
                losses[0] += lossMap[(ancestor, b[geneTree.right.traversalID.mid])]
            else:
                losses[0] += lossMap[(i, b[geneTree.left.traversalID.mid])]
                losses[0] += lossMap[(i, b[geneTree.right.traversalID.mid])]
            fn5(geneTree.left, i)
            fn5(geneTree.right, i)
    ancestorHolder = [None]
    def fn6(speciesTree, ancestor, node):
        if speciesTree.traversalID.mid == node:
            ancestorHolder[0] = ancestor
        if speciesTree.internal:
            fn6(speciesTree.left, speciesTree.traversalID.mid, node)
            fn6(speciesTree.right, speciesTree.traversalID.mid, node)
    ancestor = fn6(speciesTree, sys.maxint, b[geneTree.traversalID.mid])
    assert ancestorHolder[0] is not None
    fn5(geneTree, ancestorHolder[0])
    return len(dups), losses[0]

def calculateProbableRootOfGeneTree(speciesTree, geneTree, processID=lambda x : x):
    """
    Goes through each root possible branch making it the root. 
    Returns tree that requires the minimum number of duplications.
    """
    #get all rooted trees
    #run dup calc on each tree
    #return tree with fewest number of dups
    if geneTree.traversalID.midEnd <= 3:
        return (0, 0, geneTree)
    checkGeneTreeMatchesSpeciesTree(speciesTree, geneTree, processID)
    l = []
    def fn(tree):
        if tree.traversalID.mid != geneTree.left.traversalID.mid and tree.traversalID.mid != geneTree.right.traversalID.mid:
            newGeneTree = moveRoot(geneTree, tree.traversalID.mid)
            binaryTree_depthFirstNumbers(newGeneTree)
            dupCount, lossCount = calculateDupsAndLossesByReconcilingTrees(speciesTree, newGeneTree, processID)
            l.append((dupCount, lossCount, newGeneTree))
        if tree.internal:
            fn(tree.left)
            fn(tree.right)
    fn(geneTree)
    l.sort()
    return l[0][2], l[0][0], l[0][1]
              
#add traversalID.mid to each node name
#print tree
#parse tree
#remove names, and add them to traversalID

def main():
    pass 

def _test():
    import doctest      
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()
