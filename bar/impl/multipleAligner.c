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
#include <inttypes.h>

///////////////////////////////
///////////////////////////////
///////////////////////////////
//Stuff to represent a sequence, allowing for uncertainty at each end.
///////////////////////////////
///////////////////////////////
///////////////////////////////

SeqFrag *seqFrag_construct(const char *seq, bool missingLeftEnd, bool missingRightEnd) {
    SeqFrag *seqFrag = st_malloc(sizeof(SeqFrag));
    seqFrag->seq = stString_copy(seq);
    seqFrag->length = strlen(seq);
    seqFrag->missingLeftEnd = missingLeftEnd;
    seqFrag->missingRightEnd = missingRightEnd;
    return seqFrag;
}

void seqFrag_destruct(SeqFrag *seqFrag) {
    free(seqFrag->seq);
    free(seqFrag);
}

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

void column_destruct(Column *c) {
    while (c->nColumn != NULL) {
        Column *c2 = c->nColumn;
        free(c);
        c = c2;
    }
    free(c);
}

int column_cmp(Column *c, Column *c2) {
    //compare by sequence, then position
    return c->seqName > c2->seqName ? 1 : (c->seqName < c2->seqName ? -1 : (c->position > c2->position ? 1
            : (c->position < c2->position ? -1 : 0)));
}

uint64_t column_hashFn(Column *c) {
    return c->seqName + c->position;
}

int column_equalsFn(Column *c, Column *c2) {
    return column_cmp(c, c2) == 0;
}

stSet *makeColumns(stList *seqFrags) {
    /*
     * Makes a set of columns, each containing one sequence position. Represents
     * initially unaligned state of sequence positions.
     */
    stSet *columns = stSet_construct3((uint64_t(*)(const void *)) column_hashFn, (int(*)(const void *, const void *)) column_equalsFn,
            (void(*)(void *)) column_destruct);
    for (int64_t seq = 0; seq < stList_length(seqFrags); seq++) {
        int64_t seqLength = ((SeqFrag *)(stList_get(seqFrags, seq)))->length;
        for (int64_t pos = 0; pos < seqLength; pos++) {
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

int alignmentWeight_cmpByPosition(AlignmentWeight *aW, AlignmentWeight *aW2) {
    /*
     * Cmp by position only
     */
    return column_cmp(aW->column, aW2->column);
}

int alignmentWeight_cmpByWeight(AlignmentWeight *aW, AlignmentWeight *aW2) {
    /*
     * Cmp primarily by weight, then position.
     */
    if (aW->avgWeight == aW2->avgWeight) {
        return aW < aW2 ? -1 : (aW > aW2 ? 1 : 0);
    }
    return aW->avgWeight > aW2->avgWeight ? 1 : -1;
}

void insertWeight(AlignmentWeight *aW, stHash *alignmentWeightAdjLists) {
    /*
     * Inserts an alignment weight to an adjacency list.
     */
    stSortedSet *aWs = stHash_search(alignmentWeightAdjLists, aW->rWeight->column);
    if (aWs == NULL) {
        aWs = stSortedSet_construct3((int(*)(const void *, const void *)) alignmentWeight_cmpByPosition, free);
        stHash_insert(alignmentWeightAdjLists, aW->rWeight->column, aWs);
    }
    stSortedSet_insert(aWs, aW);
}

static Column *getColumn(stSet *columns, int64_t seqName, int64_t position) {
    Column c;
    c.seqName = seqName;
    c.position = position;
    return stSet_search(columns, &c);
}

static AlignmentWeight *makeAlignmentWeight(stSet *columns, int64_t score, int64_t seqName, int64_t position) {
    AlignmentWeight *aW = st_malloc(sizeof(AlignmentWeight));
    aW->column = getColumn(columns, seqName, position);
    assert(aW->column != NULL);
    aW->numberOfWeights = 1;
    aW->avgWeight = ((double) score) / PAIR_ALIGNMENT_PROB_1 + st_random() * 0.00001; //This randomness avoids nasty types of unbalanced trees and doesn't really affect accuracy
    return aW;
}

stHash *makeAlignmentWeightAdjacencyLists(stSet *columns, stList *multipleAlignedPairs) {
    /*
     * Make set of adjacency lists for the trivial columns and the weights linking them.
     */
    stHash *alignmentWeightAdjLists = stHash_construct2(NULL, (void(*)(void *)) stSortedSet_destruct);
    for (int64_t i = 0; i < stList_length(multipleAlignedPairs); i++) {
        /*Tuple of score, seq1, pos1, seq2, pos2 */
        stIntTuple *aP = stList_get(multipleAlignedPairs, i);
        assert(stIntTuple_length(aP) == 5);
        AlignmentWeight *aW = makeAlignmentWeight(columns, stIntTuple_get(aP, 0), stIntTuple_get(aP, 1),
                stIntTuple_get(aP, 2));
        aW->rWeight = makeAlignmentWeight(columns, stIntTuple_get(aP, 0), stIntTuple_get(aP, 3),
                stIntTuple_get(aP, 4));
        aW->rWeight->rWeight = aW;
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
    stSortedSet *alignmentWeightsOrderedByWeight = stSortedSet_construct3((int(*)(const void *, const void *)) alignmentWeight_cmpByWeight,
            NULL);
    Column *c;
    while ((c = stHash_getNext(it)) != NULL) {
        stSortedSet *aWs = stHash_search(alignmentWeightAdjLists, c);
        stSortedSetIterator *aWIt = stSortedSet_getIterator(aWs);
        AlignmentWeight *aW;
        while ((aW = stSortedSet_getNext(aWIt)) != NULL) {
            assert(aW->column != aW->rWeight->column);
            if (aW->column < aW->rWeight->column) {
                stSortedSet_insert(alignmentWeightsOrderedByWeight, aW);
            }
        }
        stSortedSet_destructIterator(aWIt);
    }
    stHash_destructIterator(it);
    return alignmentWeightsOrderedByWeight;
}

static void removeAlignmentFromSortedAlignmentWeights(AlignmentWeight *aW, stSortedSet *alignmentWeightsOrderedByWeight) {
    /*
     * Removes the weight from the ordered set.
     */
    stSortedSet_remove(alignmentWeightsOrderedByWeight, aW->column < aW->rWeight->column ? aW : aW->rWeight);
}

static void insertAlignmentIntoSortedAlignmentWeights(AlignmentWeight *aW, stSortedSet *alignmentWeightsOrderedByWeight) {
    /*
     * Adds weight to the ordered set.
     */
    stSortedSet_insert(alignmentWeightsOrderedByWeight, aW->column < aW->rWeight->column ? aW : aW->rWeight);
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
    assert(c1 != c2);
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
    stSortedSet_remove(aWs1, aW->rWeight);
    stSortedSet_remove(aWs2, aW);
    alignmentWeight_destruct(aW);
    //Now merge the remaining weights
    while (stSortedSet_size(aWs2) > 0) {
        AlignmentWeight *aW2 = stSortedSet_getFirst(aWs2);
        stSortedSet_remove(aWs2, aW2);
        stSortedSet *aWs3 = stHash_search(alignmentWeightAdjLists, aW2->column);
        stSortedSet_remove(aWs3, aW2->rWeight); //Remember to remove it from the other list!
        aW = stSortedSet_search(aWs1, aW2);
        removeAlignmentFromSortedAlignmentWeights(aW2, alignmentWeightsOrderedByWeight); //Remove old weights from alignmentWeights ordered, must do this for this weight, as we'll be altering it's column, which is used for search
        if (aW != NULL) { //Merge the weight
            removeAlignmentFromSortedAlignmentWeights(aW, alignmentWeightsOrderedByWeight); //Remove old weights from alignmentWeights ordered
            aW->avgWeight = ((aW->avgWeight * aW->numberOfWeights) + (aW2->avgWeight * aW2->numberOfWeights)) / (aW->numberOfWeights
                    + aW2->numberOfWeights);
            aW->numberOfWeights += aW2->numberOfWeights;
            aW->rWeight->avgWeight = aW->avgWeight;
            aW->rWeight->numberOfWeights = aW->numberOfWeights;
            alignmentWeight_destruct(aW2); //The old weight can now be disposed of.
            insertAlignmentIntoSortedAlignmentWeights(aW, alignmentWeightsOrderedByWeight);
        } else { //Transferring an edge across to the new vertex.
            aW2->rWeight->column = c1;
            stSortedSet_insert(aWs1, aW2);
            stSortedSet_insert(aWs3, aW2->rWeight);
            insertAlignmentIntoSortedAlignmentWeights(aW2, alignmentWeightsOrderedByWeight);
        }
    }
    assert(stSortedSet_size(aWs2) == 0);
    stSortedSet_destruct(aWs2);
}

void mergeColumns(AlignmentWeight *aW, stSet *columns, stHash *alignmentWeightAdjLists, stSortedSet *alignmentWeightsOrderedByWeight) {
    /*
     * Merges two columns together, modifying all the associate datastructures, including edge weights, to reflect the merge.
     */
    if (stSortedSet_size(stHash_search(alignmentWeightAdjLists, aW->column)) < stSortedSet_size(
            stHash_search(alignmentWeightAdjLists, aW->rWeight->column))) {
        aW = aW->rWeight;
    }
    mergeColumnsP(aW, columns, alignmentWeightAdjLists, alignmentWeightsOrderedByWeight);
}

stSet *getMultipleSequenceAlignment(stList *seqFrags, stList *multipleAlignedPairs, double gapGamma, bool checkConsistency) {
    stSet *columns = makeColumns(seqFrags);
    stHash *alignmentWeightAdjLists = makeAlignmentWeightAdjacencyLists(columns, multipleAlignedPairs);
    stSortedSet *alignmentWeightsOrderedByWeight = makeOrderedSetOfAlignmentWeights(alignmentWeightAdjLists);
    stPosetAlignment *posetAlignment = checkConsistency ? stPosetAlignment_construct(stList_length(seqFrags)) : NULL;
    while (stSortedSet_size(alignmentWeightsOrderedByWeight) > 0) {
        AlignmentWeight *aW = stSortedSet_getLast(alignmentWeightsOrderedByWeight);
        if (aW->avgWeight < gapGamma) {
            break;
        }
        removeAlignmentFromSortedAlignmentWeights(aW, alignmentWeightsOrderedByWeight);
        Column *c = aW->column, *c2 = aW->rWeight->column;
        if (!checkConsistency || (c->seqName != c2->seqName && stPosetAlignment_add(posetAlignment, c->seqName, c->position, c2->seqName, c2->position))) {
            mergeColumns(aW, columns, alignmentWeightAdjLists, alignmentWeightsOrderedByWeight);
        } else {
            stSortedSet_remove(stHash_search(alignmentWeightAdjLists, aW->rWeight->column), aW);
            stSortedSet_remove(stHash_search(alignmentWeightAdjLists, aW->column), aW->rWeight);
            alignmentWeight_destruct(aW);
        }
    }
    stSortedSet_destruct(alignmentWeightsOrderedByWeight);
    stHash_destruct(alignmentWeightAdjLists);
    if(checkConsistency) {
        stPosetAlignment_destruct(posetAlignment);
    }
    return columns;
}

/*
 * Alternative way to build a multiple alignment, by progressive alignment. Scales linearly with number
 * of sequences.
 */

typedef struct _columnPair ColumnPair;
struct _columnPair {
    int64_t xIndex, yIndex;
    int64_t score;
    ColumnPair *pPair;
    int64_t refCount;
};

ColumnPair *columnPair_construct(int64_t xIndex, int64_t yIndex, int64_t score, ColumnPair *pPair) {
    ColumnPair *cP = st_malloc(sizeof(ColumnPair));
    cP->xIndex = xIndex;
    cP->yIndex = yIndex;
    cP->score = score;
    cP->pPair = pPair;
    cP->refCount = 0;
    return cP;
}

void columnPair_destruct(ColumnPair *cP) {
    if(cP->pPair != NULL) {
        cP->pPair->refCount--;
        assert(cP->pPair->refCount >= 0);
        if(cP->pPair->refCount == 0) {
            columnPair_destruct(cP->pPair);
        }
    }
    free(cP);
}

int columnPair_cmpByYIndex(const void *c1, const void *c2) {
    int64_t i = ((ColumnPair *)c1)->yIndex;
    int64_t j = ((ColumnPair *)c2)->yIndex;
    return i > j ? 1 : (i < j ? -1 : 0);
}

stList *pairwiseAlignColumns(stList *seqXColumns, stList *seqYColumns, stHash *alignmentWeightAdjLists) {
    //Use indices of columns in list, have index --> column (obviously), but need to build column --> index.
    stHash *columnToIndexHash = stHash_construct2(NULL, (void (*)(void *))stIntTuple_destruct);
    for(int64_t i=0; i<stList_length(seqYColumns); i++) {
        stHash_insert(columnToIndexHash, stList_get(seqYColumns, i), stIntTuple_construct1(i));
    }

    //Best scoring pairs
    stSortedSet *bestScoringAlignments = stSortedSet_construct3(columnPair_cmpByYIndex, (void (*)(void *))columnPair_destruct);
    //Add in buffering first and last pairs
    stSortedSet_insert(bestScoringAlignments, columnPair_construct(-1, -1, 0, NULL));
    stSortedSet_insert(bestScoringAlignments, columnPair_construct(INT64_MAX, INT64_MAX, INT64_MAX, NULL));

    //For each column in X.
    for(int64_t i=0; i<stList_length(seqXColumns); i++) {
        Column *cX = stList_get(seqXColumns, i);
        //For each weight involving column X.
        stSortedSet *aWsX = stHash_search(alignmentWeightAdjLists, cX);
        stSortedSetIterator *aWXIt = stSortedSet_getIterator(aWsX);
        AlignmentWeight *aWX;
        while((aWX = stSortedSet_getNext(aWXIt)) != NULL) {
            //Locate index of other column
            ColumnPair cP;
            cP.yIndex = stIntTuple_get(stHash_search(columnToIndexHash, aWX->column), 0);
            //Search for highest scoring point up to but less than that index.
            ColumnPair *cPP = stSortedSet_searchLessThan(bestScoringAlignments, &cP);
            assert(cPP != NULL);
            //New score
            cP.score = cPP->score + aWX->avgWeight;
            //Find point that is equal or to the right of j
            ColumnPair *cPN = stSortedSet_searchGreaterThanOrEqual(bestScoringAlignments, &cP);
            assert(cPN != NULL);
            if(cP.score >= cPN->score || cPN->yIndex > cP.yIndex) {
                //Remove points that overlap or are to the right that score more poorly and clean them up.
                while(cP.score >= cPN->score) {
                    ColumnPair *cPNN = stSortedSet_searchGreaterThan(bestScoringAlignments, cPN);
                    assert(cPNN != NULL);
                    stSortedSet_remove(bestScoringAlignments, cPN);
                    columnPair_destruct(cPN);
                    cPN = cPNN;
                }
                cPP->refCount++; //This is for memory management
                //Insert new point.
                stSortedSet_insert(bestScoringAlignments, columnPair_construct(i, cP.yIndex, cP.score, cPP));
            }
        }
    }

    //Now traceback from highest scoring point to generate the alignment
    ColumnPair *maxPair = stSortedSet_getLast(bestScoringAlignments);
    assert(maxPair != NULL);
    stList *alignment = stList_construct();
    //For each alignment pair
    while(1) {
        //Add any unaligned Y columns
        while(--maxPair->yIndex > maxPair->pPair->yIndex) {
            stList_append(alignment, stList_get(seqYColumns, maxPair->yIndex));
        }
        //Add any unaligned X columns
        while(--maxPair->xIndex > maxPair->pPair->xIndex) {
            stList_append(alignment, stList_get(seqXColumns, maxPair->xIndex));
        }
        //Now move to previous pair
        maxPair = maxPair->pPair;
        //If this is the final pair we're done
        if(maxPair->pPair == NULL) {
            break;
        }
        //Merge two columns.
        Column *mergedColumn;
        //Add to array.
        stList_append(alignment, mergedColumn);
    }
    //Make the list of columns left-to-right
    stList_reverse(alignment);

    //Cleanup
    stSortedSet_destruct(bestScoringAlignments);
    stHash_destruct(columnToIndexHash);

    return alignment;
}

stSet *getMultipleSequenceAlignmentProgressive(stList *seqFrags, stList *multipleAlignedPairs, double gapGamma) {
    //Convert each sequence into sequence of columns.

    //Get set of weights
    //While there are multiple sequences/alignments left
    //Pick pair of alignments that is on average closest
    //Do alignment
    //Merge alignments
    //Return final set of columns.
    return NULL;
}

/*
 * Methods to extract consistent pairs.
 */

static Column *getColumn2(stHash *columns, int64_t seq, int64_t pos) {
    Column c;
    c.seqName = seq;
    c.position = pos;
    return stHash_search(columns, &c);
}

static stList *filterMultipleAlignedPairs(stSet *columns, stList *multipleAlignedPairs) {
    /*
     * Processes the list of multipleAlignedPairs and places those that align pairs within the same column in a list which is
     * returned. Pairs that do not make the list are cleaned up, as is the input list.
     */
    //Build hash of positions to columns
    stSetIterator *it = stSet_getIterator(columns);
    Column *c;
    stHash *positionsToColumns = stHash_construct3((uint64_t(*)(const void *)) column_hashFn,
            (int(*)(const void *, const void *)) column_equalsFn, NULL, NULL);
    while ((c = stSet_getNext(it)) != NULL) {
        Column *c2 = c;
        do {
            stHash_insert(positionsToColumns, c2, c);
            c2 = c2->nColumn;
        } while (c2 != NULL);
    }
    stSet_destructIterator(it);
    //Now walk through pairs
    stList *filteredMultipleAlignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    while (stList_length(multipleAlignedPairs) > 0) {
        stIntTuple *mAP = stList_pop(multipleAlignedPairs);
        if (getColumn2(positionsToColumns, stIntTuple_get(mAP, 1), stIntTuple_get(mAP, 2)) == getColumn2(
                positionsToColumns, stIntTuple_get(mAP, 3), stIntTuple_get(mAP, 4))) {
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

/*static double getAlignmentScore(stList *alignedPairs, int64_t seqLength1, int64_t seqLength2) {
 int64_t alignmentScore = 0;
 for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
 stIntTuple *alignedPair = stList_get(alignedPairs, i);
 alignmentScore += stIntTuple_getPosition(alignedPair, 0);
 }
 int64_t j = seqLength1 < seqLength2 ? seqLength1 : seqLength2;
 j = j == 0 ? 1 : j;
 double d = (double) alignmentScore / (j * PAIR_ALIGNMENT_PROB_1);
 d = d > 1.0 ? 1.0 : d;
 d = d < 0.0 ? 0.0 : d;
 return d;
 }*/

static void convertToMultipleAlignedPairs(stList *alignedPairs, stList *multipleAlignedPairs, int64_t sequence1, int64_t sequence2) {
    /*
     * Converts pairwise alignment matches to pairs with sequence indices, destroying the old pairs in the process.
     */
    while (stList_length(alignedPairs) > 0) {
        stIntTuple *alignedPair = (stIntTuple *) stList_pop(alignedPairs);
        stList_append(multipleAlignedPairs, stIntTuple_construct5(
        /* score */stIntTuple_get(alignedPair, 0), //(int64_t)(stIntTuple_getPosition(alignedPair, 0) * alignmentScore),
                /* seq 1 */sequence1, stIntTuple_get(alignedPair, 1),
                /* seq 2 */sequence2, stIntTuple_get(alignedPair, 2)));
        stIntTuple_destruct(alignedPair);
    }
    stList_destruct(alignedPairs);
}

static void addMultipleAlignedPairs(int64_t sequence1, int64_t sequence2, stList *seqFrags, stList *multipleAlignedPairs,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Computes a pairwise alignment and returns the pairwise match probabilities as tuples of (score, seq1, pos1, seq2, pos2).
     */
    SeqFrag *seqFrag1 = stList_get(seqFrags, sequence1);
    SeqFrag *seqFrag2 = stList_get(seqFrags, sequence2);
    stList *alignedPairs = getAlignedPairs(seqFrag1->seq, seqFrag2->seq, pairwiseAlignmentBandingParameters, seqFrag1->missingLeftEnd || seqFrag2->missingLeftEnd, seqFrag1->missingRightEnd || seqFrag2->missingRightEnd);
    convertToMultipleAlignedPairs(alignedPairs, multipleAlignedPairs, sequence1, sequence2);
}

static void addMultipleAlignedPairsFilteringInconsistentPairs(int64_t sequence1, int64_t sequence2, stList *seqFrags,
        stList *multipleAlignedPairs, stList *discardedMultipleAlignedPairs,
        float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Get pairwise alignment between sequences, filtering pairs in to those in a consistent pairwise alignment and those that are inconsistent with the first alignment.
     */
    SeqFrag *seqFrag1 = stList_get(seqFrags, sequence1);
    SeqFrag *seqFrag2 = stList_get(seqFrags, sequence2);
    stList *alignedPairs = getAlignedPairs(seqFrag1->seq, seqFrag2->seq, pairwiseAlignmentBandingParameters, seqFrag1->missingLeftEnd || seqFrag2->missingLeftEnd, seqFrag1->missingRightEnd || seqFrag2->missingRightEnd);
    stList *discardedAlignedPairs = filterPairwiseAlignmentToMakePairsOrdered(alignedPairs, gapGamma);
    convertToMultipleAlignedPairs(alignedPairs, multipleAlignedPairs, sequence1, sequence2);
    convertToMultipleAlignedPairs(discardedAlignedPairs, discardedMultipleAlignedPairs, sequence1, sequence2);
}

stList *makeAllPairwiseAlignments(stList *seqFrags, PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Generate the set of pairwise alignments between the sequences.
     */
    stList *multipleAlignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    int64_t seqNo = stList_length(seqFrags);
    for (int64_t seq1 = 0; seq1 < seqNo; seq1++) {
        for (int64_t seq2 = seq1 + 1; seq2 < seqNo; seq2++) {
            addMultipleAlignedPairs(seq1, seq2, seqFrags, multipleAlignedPairs, pairwiseAlignmentBandingParameters);
        }
    }
    return multipleAlignedPairs;
}

stList *makeAlignmentUsingAllPairs(stList *seqFrags, float gapGamma, PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Generate a multiple alignment considering all pairs of sequences.
     */
    stList *multipleAlignedPairs = makeAllPairwiseAlignments(seqFrags, pairwiseAlignmentBandingParameters);
    stSet *columns = getMultipleSequenceAlignment(seqFrags, multipleAlignedPairs, gapGamma, 1);
    multipleAlignedPairs = filterMultipleAlignedPairs(columns, multipleAlignedPairs);
    stSet_destruct(columns);
    return multipleAlignedPairs;
}

///////////////////////////////
///////////////////////////////
///////////////////////////////
//Functions used for selecting which pairwise alignments to compute, so that fewer than n choose 2 pairs need
//be calculated to create a good alignment.
///////////////////////////////
///////////////////////////////
///////////////////////////////

static stIntTuple *makePairToAlign(int64_t seq1, int64_t seq2) {
    assert(seq1 != seq2);
    return seq1 < seq2 ? stIntTuple_construct2( seq1, seq2) : stIntTuple_construct2( seq2, seq1);
}

stList *getReferencePairwiseAlignments(stList *seqFrags) {
    /*
     * Picks a reference sequence ensures everyone is aligned.
     */
    stList *l = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(seqFrags); i++) {
        stList_append(l, stIntTuple_construct2(((SeqFrag *)stList_get(seqFrags, i))->length, i));
    }
    stList_sort(l, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    stList *chosenPairsOfSequencesToAlign = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    if (stList_length(l) > 0) {
        stIntTuple *k = stList_get(l, stList_length(l) / 2);
        for (int64_t i = 0; i < stList_length(l); i++) {
            stIntTuple *j = stList_get(l, i);
            if (k != j) {
                stList_append(chosenPairsOfSequencesToAlign,
                        makePairToAlign(stIntTuple_get(j, 1), stIntTuple_get(k, 1)));
            }
        }
    }
    stList_destruct(l);
    return chosenPairsOfSequencesToAlign;
}

stSortedSet *getReferencePairwiseAlignments2(stList *seqFrags) {
    /*
     * Same as above, but returns sorted set.
     */
    stList *l = getReferencePairwiseAlignments(seqFrags);
    stSortedSet *s = stList_getSortedSet(l, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    stList_setDestructor(l, NULL);
    stList_destruct(l);
    stSortedSet_setDestructor(s, (void (*)(void *))stIntTuple_destruct);
    return s;
}

/*
 * Following functions used to for distance matrix calculated from MSA, represented as a set of columns.
 */

int64_t *getValue(int64_t seq1, int64_t seq2, int64_t *distanceCounts, int64_t seqNo) {
    return &distanceCounts[seq1 * seqNo + seq2];
}

int64_t *getSubs(int64_t seq1, int64_t seq2, int64_t *distanceCounts, int64_t seqNo) {
    return seq1 > seq2 ? getValue(seq1, seq2, distanceCounts, seqNo) : getValue(seq2, seq1, distanceCounts, seqNo);
}

int64_t *getNonSubs(int64_t seq1, int64_t seq2, int64_t *distanceCounts, int64_t seqNo) {
    return seq1 < seq2 ? getValue(seq1, seq2, distanceCounts, seqNo) : getValue(seq2, seq1, distanceCounts, seqNo);
}

double subsPerSite(int64_t seq1, int64_t seq2, int64_t *distanceCounts, int64_t seqNo) {
    /*
     * Returns subs/(subs + identity sites).
     */
    int64_t subs = *getSubs(seq1, seq2, distanceCounts, seqNo);
    int64_t identities = *getNonSubs(seq1, seq2, distanceCounts, seqNo);
    return (subs + identities == 0) ? 0.0 : subs / ((double) subs + identities);
}

int64_t *getDistanceMatrix(stSet *columns, stList *seqFrags, int64_t maxPairsToConsider) {
    /*
     * Builds a distance matrix from the deepest 'columnsToConsider' columns.
     * Return matrix has for each pair of seqs number of substitutions observed in the alignment
     * and number of identity sights, i.e. sights that have remained the same.
     * Stops calculating pairs when number of pairs of sequence positions compared exceeds maxPairsToConsider
     */
    int64_t *distanceCounts = st_calloc((int64_t) stList_length(seqFrags) * stList_length(seqFrags), sizeof(int64_t));
    stSetIterator *setIt = stSet_getIterator(columns);
    Column *c;
    int64_t seqNo = stList_length(seqFrags);
    int64_t pairsConsidered = 0;
    while ((c = stSet_getNext(setIt)) != NULL && pairsConsidered < maxPairsToConsider) {
        Column *c1 = c;
        do {
            int64_t seq1 = c1->seqName;
            char base1 = ((SeqFrag *) stList_get(seqFrags, seq1))->seq[c1->position];
            Column *c2 = c1->nColumn;
            while (c2 != NULL) {
                int64_t seq2 = c2->seqName;
                char base2 = ((SeqFrag *) stList_get(seqFrags, seq2))->seq[c2->position];
                (*(int64_t *) (base1 == base2 ? getNonSubs : getSubs)(seq1, seq2, distanceCounts, seqNo))++;
                c2 = c2->nColumn;
                pairsConsidered++;
            }
            c1 = c1->nColumn;
        } while (c1 != NULL);
    }
    stSet_destructIterator(setIt);
    return distanceCounts;
}

stGraph *makeAdjacencyList(int64_t *distanceCounts, int64_t seqNo, stSortedSet *chosenPairsOfSequencesToAlign) {
    /*
     * Creates a graph in which the seqs are the vertices and the edges are the chosen pairwise alignments, with
     * weights that are equal to the number of subs per site.
     */
    stGraph *g = stGraph_construct(seqNo);
    stSortedSetIterator *pairIt = stSortedSet_getIterator(chosenPairsOfSequencesToAlign);
    stIntTuple *pairToAlign;
    while ((pairToAlign = stSortedSet_getNext(pairIt)) != NULL) {
        int64_t seq1 = stIntTuple_get(pairToAlign, 0);
        int64_t seq2 = stIntTuple_get(pairToAlign, 1);
        stGraph_addEdge(g, seq1, seq2, subsPerSite(seq1, seq2, distanceCounts, seqNo));
    }
    stSortedSet_destructIterator(pairIt);
    return g;
}

int64_t getNextBestPair(int64_t seq1, int64_t *distanceCounts, int64_t seqNo, stSortedSet *chosenPairsOfSequencesToAlign) {
    /*
     * Selects the best next pairwise alignment for seq1 to compute, which we define as the alignment where the difference
     * between the current path of alignments and the predicted pairwise alignment distance is greatest.
     */
    //Do dijkstra's first
    stGraph *graph = makeAdjacencyList(distanceCounts, seqNo, chosenPairsOfSequencesToAlign);
    double *distances = stGraph_shortestPaths(graph, seq1);
    double maxGain = INT64_MIN;
    int64_t maxGainSeq = INT64_MAX;
    for (int64_t seq2 = 0; seq2 < seqNo; seq2++) {
        if (seq1 != seq2) {
            double gain = distances[seq2] - subsPerSite(seq1, seq2, distanceCounts, seqNo);
            if (gain > maxGain) {
                maxGain = gain;
                maxGainSeq = seq2;
            }
        }
    }
    stGraph_destruct(graph);
    free(distances);
    return maxGainSeq;
}

stList *makeAlignment(stList *seqFrags, int64_t spanningTrees, int64_t maxPairsToConsider, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Computes an MSA, making up to "spanningTrees"*no of seqs pairwise alignments.
     */
    if (spanningTrees == 0) {
        return stList_construct();
    }
    int64_t seqNo = stList_length(seqFrags);
    if (spanningTrees * (seqNo - 1) >= (seqNo * (seqNo - 1)) / 2) { //Do all pairs if we can
        return makeAlignmentUsingAllPairs(seqFrags, gapGamma, pairwiseAlignmentBandingParameters);
    }
    stList *multipleAlignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct); //pairwise alignment pairs, with sequence indices
    stList *discardedMultipleAlignedPairs = stList_construct3(0, (void (*)(void *))stIntTuple_destruct); //pairwise alignment pairs that were inconsistent with some of the pairs in the multipleAlignedPairs set
    stSortedSet *chosenPairsOfSequencesToAlign = getReferencePairwiseAlignments2(seqFrags);
    stSortedSetIterator *pairIt = stSortedSet_getIterator(chosenPairsOfSequencesToAlign);
    stIntTuple *pairToAlign;
    while ((pairToAlign = stSortedSet_getNext(pairIt)) != NULL) {
        //We get pairwise alignments, for this first alignment we filter the pairs greedily to make them consistent
        addMultipleAlignedPairsFilteringInconsistentPairs(stIntTuple_get(pairToAlign, 0), stIntTuple_get(pairToAlign, 1),
                seqFrags, multipleAlignedPairs, discardedMultipleAlignedPairs, gapGamma, pairwiseAlignmentBandingParameters);
    }
    stSortedSet_destructIterator(pairIt);
    if (spanningTrees == 1) { //As there is no further alignments to compute we exit, as all the pairs are already consistent.
        stSortedSet_destruct(chosenPairsOfSequencesToAlign);
        stList_destruct(discardedMultipleAlignedPairs);
        //We are done
        return multipleAlignedPairs;
    }
    int64_t iteration = 1;
    bool checkConsistency = 0; //The first alignment of multiple aligned pairs is already consistent
    while (1) {
        stSet *columns = getMultipleSequenceAlignment(seqFrags, multipleAlignedPairs, gapGamma, checkConsistency);
        if (iteration++ >= spanningTrees) {
            stSortedSet_destruct(chosenPairsOfSequencesToAlign);
            multipleAlignedPairs = filterMultipleAlignedPairs(columns, multipleAlignedPairs);
            stSet_destruct(columns);
            stList_destruct(discardedMultipleAlignedPairs);
            return multipleAlignedPairs;
        }
        int64_t *distanceCounts = getDistanceMatrix(columns, seqFrags, maxPairsToConsider);
        stSet_destruct(columns);
        for (int64_t seq = 0; seq < stList_length(seqFrags); seq++) {
            int64_t otherSeq = getNextBestPair(seq, distanceCounts, seqNo, chosenPairsOfSequencesToAlign);
            if (otherSeq != INT64_MAX) {
                assert(seq != otherSeq);
                stIntTuple *pairToAlign = makePairToAlign(seq, otherSeq);
                if (stSortedSet_search(chosenPairsOfSequencesToAlign, pairToAlign) == NULL) {
                    addMultipleAlignedPairs(seq, otherSeq, seqFrags, multipleAlignedPairs, pairwiseAlignmentBandingParameters);
                    stSortedSet_insert(chosenPairsOfSequencesToAlign, pairToAlign);
                } else {
                    stIntTuple_destruct(pairToAlign);
                }
            }
        }
        if (!checkConsistency) {
            while (stList_length(discardedMultipleAlignedPairs) > 0) {
                stList_append(multipleAlignedPairs, stList_pop(discardedMultipleAlignedPairs)); //Add back the filtered pairs, which may be next set of alignments gathered.
            }
            checkConsistency = 1;
        }
        free(distanceCounts);
    }
}
