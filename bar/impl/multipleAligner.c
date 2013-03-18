/*
 * Copyright (C) 2013 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "multipleAligner.h"
#include "sonLib.h"
#include "stPosetAlignment.h"
#include "pairwiseAligner.h"
#include <stdlib.h>
#include <math.h>
#include "stGraph.h"

///////////////////////////////
///////////////////////////////
///////////////////////////////
//Stuff to represent the multiple alignment columns and pairwise weights between the columns
///////////////////////////////
///////////////////////////////
///////////////////////////////

/*
 * Structure represents an alignment of sequence positions.
 */
typedef struct _column Column;
struct _column {
    int32_t seqName;
    int32_t position;
    Column *nColumn;
};

static void column_destruct(Column *c) {
    while (c->nColumn != NULL) {
        Column *c2 = c->nColumn;
        free(c);
        c = c2;
    }
    free(c);
}

static int column_cmp(Column *c, Column *c2) {
    //compare by sequence, then position
    return c->seqName > c2->seqName ? 1 : (c->seqName < c2->seqName ? -1 : (c->position > c2->position ? 1
            : (c->position < c2->position ? -1 : 0)));
}

static uint32_t column_hashFn(Column *c) {
    return c->seqName + c->position;
}

static int column_equalsFn(Column *c, Column *c2) {
    return column_cmp(c, c2) == 0;
}

stSet *makeColumns(stList *sequences) {
    /*
     * Makes a set of columns, each containing one sequence position. Represents
     * initially unaligned state of sequence positions.
     */
    stSet *columns = stSet_construct3((uint32_t(*)(const void *)) column_hashFn, (int(*)(const void *, const void *)) column_equalsFn,
            (void(*)(void *)) column_destruct);
    for (int32_t seq = 0; seq < stList_length(sequences); seq++) {
        int32_t seqLength = strlen(stList_get(sequences, seq));
        for (int32_t pos = 0; pos < seqLength; pos++) {
            Column *c = st_malloc(sizeof(Column));
            c->seqName = seq;
            c->position = pos;
            c->nColumn = NULL;
            stSet_insert(columns, c);
        }
    }
    return columns;
}

/*
 * Structure represents pairwise alignment matches between positions in two columns.
 */
typedef struct _alignmentWeight AlignmentWeight;
struct _alignmentWeight {
    Column *column;
    double avgWeight;
    double numberOfWeights;
    AlignmentWeight *rWeight;
};

static int alignmentWeight_cmpByPosition(AlignmentWeight *aW, AlignmentWeight *aW2) {
    /*
     * Cmp by position only
     */
    return column_cmp(aW->column, aW2->column);
}

static int alignmentWeight_cmpByWeight(AlignmentWeight *aW, AlignmentWeight *aW2) {
    /*
     * Cmp primarily by weight, then position.
     */
    return aW->avgWeight > aW2->avgWeight ? 1 : (aW->avgWeight < aW2->avgWeight ? -1 : alignmentWeight_cmpByPosition(aW, aW2));
}

static void insertWeight(AlignmentWeight *aW, stHash *alignmentWeightAdjLists) {
    /*
     * Inserts an alignment weight to an adjacency list.
     */
    stSortedSet *aWs = stHash_search(alignmentWeightAdjLists, aW->column);
    if (aWs == NULL) {
        aWs = stSortedSet_construct3((int(*)(const void *, const void *)) alignmentWeight_cmpByPosition, free);
        stHash_insert(alignmentWeightAdjLists, aW->rWeight->column, aWs);
    }
    stSortedSet_insert(aWs, aW);
}

stHash *makeAlignmentWeightAdjacencyLists(stSet *columns, stList *multipleAlignedPairs) {
    /*
     * Make set of adjacency lists for the trivial columns and the weights linking them.
     */
    stHash *alignmentWeightAdjLists = stHash_construct2(NULL, (void(*)(void *)) stSortedSet_destruct);
    for (int32_t i = 0; i < stList_length(multipleAlignedPairs); i++) {
        /*Tuple of score, seq1, pos1, seq2, pos2 */
        stIntTuple *aP = stList_get(multipleAlignedPairs, i);
        AlignmentWeight *aW = st_malloc(sizeof(AlignmentWeight));
        aW->rWeight = st_malloc(sizeof(AlignmentWeight));
        Column c;
        c.seqName = stIntTuple_getPosition(aP, 1);
        c.position = stIntTuple_getPosition(aP, 2);
        aW->column = stSet_search(columns, &c);
        c.seqName = stIntTuple_getPosition(aP, 3);
        c.position = stIntTuple_getPosition(aP, 4);
        aW->rWeight->column = stSet_search(columns, &c);
        aW->numberOfWeights = 1;
        aW->rWeight->numberOfWeights = 1;
        aW->avgWeight = stIntTuple_getPosition(aP, 0);
        aW->rWeight->avgWeight = aW->avgWeight;
        insertWeight(aW, alignmentWeightAdjLists);
        insertWeight(aW->rWeight, alignmentWeightAdjLists);
    }
    return alignmentWeightAdjLists;
}

stSortedSet *makeOrderedSetOfAlignmentWeights(stHash *alignmentWeightAdjLists) {
    /*
     * Convert the adjacency lists into an ordered set, ordered by weight.
     */
    stHashIterator *it = stHash_getIterator(alignmentWeightAdjLists);
    stSortedSet *alignmentWeightsOrderedByWeight = stSortedSet_construct3((int(*)(const void *, const void *)) alignmentWeight_cmpByWeight, NULL);
    stSortedSet *aWs;
    while ((aWs = stHash_getNext(it)) != NULL) {
        stSortedSetIterator *aWIt = stSortedSet_getIterator(aWs);
        AlignmentWeight *aW;
        while ((aW = stSortedSet_getNext(aWIt)) != NULL) {
            if (aW->column < aW->rWeight->column) {
                stSortedSet_insert(alignmentWeightsOrderedByWeight, aW);
            }
        }
    }
    return alignmentWeightsOrderedByWeight;
}

static void removeAlignmentFromSortedAlignmentWeights(AlignmentWeight *aW, stSortedSet *alignmentWeightsOrderedByWeight) {
    /*
     * Removes the weight from the ordered set.
     */
    stSortedSet_remove(alignmentWeightsOrderedByWeight, aW->column < aW->rWeight->column ? aW : aW->rWeight);
}

static void alignmentWeight_destruct(AlignmentWeight *aW) {
    /*
     * Cleans up a weight and its reverse.
     */
    free(aW->rWeight);
    free(aW);
}

static void mergeColumnsP(AlignmentWeight *aW, stSet *columns, stHash *alignmentWeightAdjLists,
        stSortedSet *alignmentWeightsOrderedByWeight) {
    Column *c1 = aW->column, *c2 = aW->rWeight->column;
    stSortedSet *aWs1 = stHash_search(alignmentWeightAdjLists, c1);
    stSortedSet *aWs2 = stHash_remove(alignmentWeightAdjLists, c2);
    assert(stSortedSet_size(aWs1) >= stSortedSet_size(aWs2));
    //Merge the columns
    Column *c = c1;
    while (c->nColumn != NULL) {
        c = c->nColumn;
    }
    c->nColumn = c2;
    stSet_remove(columns, c2);
    //Cleanup the merging weight
    stSortedSet_remove(aWs1, aW);
    stSortedSet_remove(aWs2, aW->rWeight);
    alignmentWeight_destruct(aW);
    //Now merge the remaining weights
    while (stSortedSet_size(aWs2) > 0) {
        AlignmentWeight *aW2 = stSortedSet_getFirst(aWs2);
        stSortedSet_remove(aWs2, aW2);
        assert(aW2->column != c1);
        aW = stSortedSet_search(aWs1, aW2);
        if (aW != NULL) { //Merge the weight
            removeAlignmentFromSortedAlignmentWeights(aW2, alignmentWeightsOrderedByWeight); //Remove old weights from alignmentWeights ordered
            removeAlignmentFromSortedAlignmentWeights(aW, alignmentWeightsOrderedByWeight); //Remove old weights from alignmentWeights ordered
            aW->avgWeight = ((aW->avgWeight * aW->numberOfWeights) + (aW2->avgWeight * aW2->numberOfWeights)) / (aW->numberOfWeights
                    + aW2->numberOfWeights);
            aW->numberOfWeights += aW2->numberOfWeights;
            aW->rWeight->avgWeight = aW->avgWeight;
            aW->rWeight->numberOfWeights = aW->numberOfWeights;
            alignmentWeight_destruct(aW2);
            stSortedSet_insert(alignmentWeightsOrderedByWeight, aW->column < aW->rWeight->column ? aW : aW->rWeight);
        } else { //No merging,just go for it.
            aW2->rWeight->column = c1;
            stSortedSet_insert(aWs1, aW2);
        }
    }
    stSortedSet_destruct(aWs2);
}

void mergeColumns(AlignmentWeight *aW, stSet *columns, stHash *alignmentWeightAdjLists, stSortedSet *alignmentWeightsOrderedByWeight) {
    /*
     * Merges two columns together, modifying all the associate datastructures, including edge weights, to reflect the merge.
     */
    stSortedSet_size(stHash_search(alignmentWeightAdjLists, aW->column)) <= stSortedSet_size(
            stHash_search(alignmentWeightAdjLists, aW->rWeight->column)) ? mergeColumnsP(aW, columns, alignmentWeightAdjLists,
            alignmentWeightsOrderedByWeight) : mergeColumns(aW->rWeight, columns, alignmentWeightAdjLists, alignmentWeightsOrderedByWeight);
}

stSet *getMultipleSequenceAlignment(stList *sequences, stList *multipleAlignedPairs, double gapGamma) {
    stSet *columns = makeColumns(sequences);
    stHash *alignmentWeightAdjLists = makeAlignmentWeightAdjacencyLists(columns, multipleAlignedPairs);
    stSortedSet *alignmentWeightsOrderedByWeight = makeOrderedSetOfAlignmentWeights(alignmentWeightAdjLists);
    stPosetAlignment *posetAlignment = stPosetAlignment_construct(stList_length(sequences));
    while (stSortedSet_size(alignmentWeightsOrderedByWeight) > 0) {
        AlignmentWeight *aW = stSortedSet_getFirst(alignmentWeightsOrderedByWeight);
        stSortedSet_remove(alignmentWeightsOrderedByWeight, aW);
        Column *c = aW->column, *c2 = aW->rWeight->column;
        if (aW->avgWeight >= gapGamma && stPosetAlignment_add(posetAlignment, c->seqName, c->position, c2->seqName, c2->position)) {
            mergeColumns(aW, columns, alignmentWeightAdjLists, alignmentWeightsOrderedByWeight);
        } else {
            stSortedSet_remove(stHash_search(alignmentWeightAdjLists, aW->column), aW);
            stSortedSet_remove(stHash_search(alignmentWeightAdjLists, aW->rWeight->column), aW->rWeight);
            alignmentWeight_destruct(aW);
        }
    }
    stSortedSet_destruct(alignmentWeightsOrderedByWeight);
    stHash_destruct(alignmentWeightAdjLists);
    stPosetAlignment_destruct(posetAlignment);
    return columns;
}

static Column *getColumn(stHash *columns, int32_t seq, int32_t pos) {
    Column c;
    c.seqName = seq;
    c.position = pos;
    return stHash_search(columns, &c);
}

stList *filterMultipleAlignedPairs(stSet *columns, stList *multipleAlignedPairs) {
    /*
     * Processes the list of multipleAlignedPairs and places those that align pairs within the same column in a list which is
     * returned. Pairs that do not make the list are cleaned up, as is the input list.
     */
    //Build hash of positions to columns
    stSetIterator *it = stSet_getIterator(columns);
    Column *c;
    stHash *positionsToColumns = stHash_construct3((uint32_t (*)(const void *))column_hashFn, (int (*)(const void *, const void *))column_equalsFn, NULL, NULL);
    while ((c = stSet_getNext(it)) != NULL) {
        Column *c2 = c;
        do {
            stHash_insert(positionsToColumns, c2, c);
            c2 = c2->nColumn;
        } while (c2 != NULL);
    }
    //Now walk through pairs
    stList *filteredMultipleAlignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    while (stList_length(multipleAlignedPairs) > 0) {
        stIntTuple *mAP = stList_pop(multipleAlignedPairs);
        if (getColumn(positionsToColumns, stIntTuple_getPosition(mAP, 1), stIntTuple_getPosition(mAP, 2)) == getColumn(positionsToColumns,
                stIntTuple_getPosition(mAP, 3), stIntTuple_getPosition(mAP, 4))) {
            stList_append(filteredMultipleAlignedPairs, mAP);
        } else {
            stIntTuple_destruct(mAP);
        }
    }
    //Cleanup
    stHash_destruct(positionsToColumns);
    stList_destruct(multipleAlignedPairs);
    return filteredMultipleAlignedPairs;
}

///////////////////////////////
///////////////////////////////
///////////////////////////////
//Functions used for selecting which pairwise alignments to compute.
///////////////////////////////
///////////////////////////////
///////////////////////////////

stList *getSpanningTree(stList *sequences) {
    /*
     * Returns a set of pairs such that all the sequences are in one connected component, picking pairs
     * such that the pairs of sequences being aligned have similar lengths.
     */
    stList *l = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(sequences); i++) {
        stIntTuple_construct(2, strlen(stList_get(sequences, i)), i);
    }
    stList_sort(l, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    stList *chosenPairsOfSequencesToAlign = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(sequences); i++) {
        stIntTuple *j = stList_get(sequences, i), *k = stList_get(sequences, i + 1);
        stList_append(chosenPairsOfSequencesToAlign, stIntTuple_construct(2, stIntTuple_getPosition(j, 1), stIntTuple_getPosition(k, 1)));
    }
    stList_destruct(l);
    return chosenPairsOfSequencesToAlign;
}

/*
 * Following functions used to for distance matrix calculated from MSA, represented as a set of columns.
 */

int32_t *getValue(int32_t seq1, int32_t seq2, int32_t *distanceCounts, int32_t sequenceNo) {
    return &distanceCounts[seq1 * sequenceNo + seq2];
}

int32_t *getSubs(int32_t seq1, int32_t seq2, int32_t *distanceCounts, int32_t sequenceNo) {
    return seq1 > seq2 ? getValue(seq2, seq1, distanceCounts, sequenceNo) : getValue(seq1, seq2, distanceCounts, sequenceNo);
}

int32_t *getNonSubs(int32_t seq1, int32_t seq2, int32_t *distanceCounts, int32_t sequenceNo) {
    return seq1 > seq2 ? getValue(seq2, seq1, distanceCounts, sequenceNo) : getValue(seq1, seq2, distanceCounts, sequenceNo);
}

double subsPerSite(int32_t seq1, int32_t seq2, int32_t *distanceCounts, int32_t sequenceNo) {
    /*
     * Returns subs/(subs + identity sites).
     */
    int32_t subs = *getSubs(seq1, seq2, distanceCounts, sequenceNo);
    int32_t identities = *getNonSubs(seq1, seq2, distanceCounts, sequenceNo);
    return subs/((double)subs + identities);
}

int32_t *getDistanceMatrix(stSet *columns, stList *sequences, int64_t maxPairsToConsider) {
    /*
     * Builds a distance matrix from the deepest 'columnsToConsider' columns.
     * Return matrix has for each pair of sequences number of substitutions observed in the alignment
     * and number of identity sights, i.e. sights that have remained the same.
     * Stops calculating pairs when number of pairs of sequence positions compared exceeds maxPairsToConsider
     */
    int32_t *distanceCounts = st_calloc((int64_t)stList_length(sequences) * stList_length(sequences), sizeof(int32_t));
    stSetIterator *setIt = stSet_getIterator(columns);
    Column *c;
    int32_t seqNo = stList_length(sequences);
    int64_t pairsConsidered=0;
    while ((c = stSet_getNext(setIt)) != NULL && pairsConsidered < maxPairsToConsider) {
        Column *c1 = c;
        int32_t seq1 = c1->seqName;
        char base1 = ((char *) stList_get(sequences, c1->seqName))[c1->position];
        do {
            Column *c2 = c1->nColumn;
            while (c2 != NULL) {
                int32_t seq2 = c2->seqName;
                char base2 = ((char *) stList_get(sequences, seq2))[c2->position];
                (*(int32_t *)(base1 == base2 ? getSubs : getNonSubs)(seq1, seq2, distanceCounts, seqNo))++;
                c2 = c2->nColumn;
                pairsConsidered++;
            }
            c1 = c1->nColumn;
        } while (c1 != NULL);
    }
    stSet_destructIterator(setIt);
    return distanceCounts;
}

stGraph *makeAdjacencyList(int32_t *distanceCounts, int32_t sequenceNo, stList *chosenPairsOfSequencesToAlign) {
    /*
     * Creates a graph in which the sequences are the vertices and the edges are the chosen pairwise alignments, with
     * weights that are equal to the number of subs per site.
     */
    stGraph *g = stGraph_construct(sequenceNo);
    for(int32_t i=0; i<stList_length(chosenPairsOfSequencesToAlign); i++) {
        stIntTuple *pairToAlign = stList_get(chosenPairsOfSequencesToAlign, i);
        int32_t seq1 = stIntTuple_getPosition(pairToAlign, 0);
        int32_t seq2 = stIntTuple_getPosition(pairToAlign, 1);
        stGraph_addEdge(g, seq1, seq2, subsPerSite(seq1, seq2, distanceCounts, sequenceNo));
    }
    return g;
}

int32_t getNextBestPair(int32_t seq1, int32_t *distanceCounts, int32_t sequenceNo, stList *chosenPairsOfSequencesToAlign) {
    /*
     * Selects the best next pairwise alignment for seq1 to compute, which we define as the alignment where the difference
     * between the current path of alignments and the predicted pairwise alignment distance is greatest.
     */
    //Do dijkstra's first
    stGraph *graph = makeAdjacencyList(distanceCounts, sequenceNo, chosenPairsOfSequencesToAlign);
    double *distances = stGraph_shortestPaths(graph, seq1);
    double maxGain = -INT64_MAX;
    int32_t maxGainSeq = INT32_MAX;
    for(int32_t seq2=0; seq2<sequenceNo; seq2++) {
        if(seq1 != seq2) {
            double gain = distances[seq2] - subsPerSite(seq1, seq2, distanceCounts, sequenceNo);
            if(gain > maxGain) {
                maxGain = gain;
                maxGainSeq = seq2;
            }
        }
    }
    stGraph_destruct(graph);
    free(distances);
    return maxGainSeq;
}

/*static double getAlignmentScore(stList *alignedPairs, int32_t seqLength1, int32_t seqLength2) {
    int64_t alignmentScore = 0;
    for (int32_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        assert(stIntTuple_length(alignedPair) == 3);
        alignmentScore += stIntTuple_getPosition(alignedPair, 0);
    }
    int64_t j = seqLength1 < seqLength2 ? seqLength1 : seqLength2;
    j = j == 0 ? 1 : j;
    double d = (double) alignmentScore / (j * PAIR_ALIGNMENT_PROB_1);
    d = d > 1.0 ? 1.0 : d;
    d = d < 0.0 ? 0.0 : d;
    return d;
}*/

void addMultipleAlignedPairs(int32_t sequence1, int32_t sequence2, stList *sequences, stList *multipleAlignedPairs,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Computes a pairwise alignment and returns the pairwise match probabilities as tuples of (score, seq1, pos1, seq2, pos2).
     */
    char *string1 = stList_get(sequences, sequence1);
    char *string2 = stList_get(sequences, sequence2);
    stList *alignedPairs = getAlignedPairs(string1, string2, pairwiseAlignmentBandingParameters);
    //double alignmentScore = getAlignmentScore(alignedPairs, strlen(string1), strlen(string2));
    while (stList_length(alignedPairs) > 0) {
        stIntTuple *alignedPair = (stIntTuple *) stList_pop(alignedPairs);
        stList_append(multipleAlignedPairs, stIntTuple_construct(5,
        /* score */stIntTuple_getPosition(alignedPair, 0),
        /* seq 1 */sequence1, stIntTuple_getPosition(alignedPair, 1),
        /* seq 2 */sequence2, stIntTuple_getPosition(alignedPair, 2)));
    }
}

stList *makeAlignment(stList *sequences, int32_t spanningTrees, int64_t maxPairsToConsider, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Computes an MSA, making up to "spanningTrees"*no of sequences pairwise alignments.
     */
    if (spanningTrees == 0) {
        return stList_construct();
    }
    stList *chosenPairsOfSequencesToAlign = getSpanningTree(sequences);
    stList *multipleAlignedPairs = stList_construct();
    for (int32_t i = 0; i < stList_length(chosenPairsOfSequencesToAlign); i++) {
        stIntTuple *pairToAlign = stList_get(chosenPairsOfSequencesToAlign, i);
        addMultipleAlignedPairs(stIntTuple_getPosition(pairToAlign, 0), stIntTuple_getPosition(pairToAlign, 1), sequences,
                multipleAlignedPairs, pairwiseAlignmentBandingParameters);
    }
    int32_t iteration = 1;
    while (1) {
        stSet *columns = getMultipleSequenceAlignment(sequences, multipleAlignedPairs, gapGamma);
        if (iteration++ >= spanningTrees) {
            stList_destruct(chosenPairsOfSequencesToAlign);
            multipleAlignedPairs = filterMultipleAlignedPairs(columns, multipleAlignedPairs);
            stSet_destruct(columns);
            return multipleAlignedPairs;
        }
        int32_t *distanceCounts = getDistanceMatrix(columns, sequences, maxPairsToConsider);
        for (int32_t seq = 0; seq < stList_length(sequences); seq++) {
            int32_t otherSeq = getNextBestPair(seq, distanceCounts, stList_length(sequences), chosenPairsOfSequencesToAlign);
            if (otherSeq != INT32_MAX) {
                addMultipleAlignedPairs(seq, otherSeq, sequences, multipleAlignedPairs, pairwiseAlignmentBandingParameters);
                stList_append(chosenPairsOfSequencesToAlign, stIntTuple_construct(2, seq, otherSeq));
            }
        }
        free(distanceCounts);
    }
}

