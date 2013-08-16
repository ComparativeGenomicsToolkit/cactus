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

SeqFrag *seqFrag_construct(const char *seq, int64_t leftEndId, int64_t rightEndId) {
    SeqFrag *seqFrag = st_malloc(sizeof(SeqFrag));
    seqFrag->seq = stString_copy(seq);
    seqFrag->length = strlen(seq);
    seqFrag->leftEndId = leftEndId;
    seqFrag->rightEndId = rightEndId;
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
        int64_t seqLength = ((SeqFrag *) (stList_get(seqFrags, seq)))->length;
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
        AlignmentWeight *aW = makeAlignmentWeight(columns, stIntTuple_get(aP, 0), stIntTuple_get(aP, 1), stIntTuple_get(aP, 2));
        aW->rWeight = makeAlignmentWeight(columns, stIntTuple_get(aP, 0), stIntTuple_get(aP, 3), stIntTuple_get(aP, 4));
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

static Column *mergeColumnsP(AlignmentWeight *aW, stSet *columns, stHash *alignmentWeightAdjLists,
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
    return c1;
}

Column *mergeColumns(AlignmentWeight *aW, stSet *columns, stHash *alignmentWeightAdjLists, stSortedSet *alignmentWeightsOrderedByWeight) {
    /*
     * Merges two columns together, modifying all the associate datastructures, including edge weights, to reflect the merge.
     */
    if (stSortedSet_size(stHash_search(alignmentWeightAdjLists, aW->column)) < stSortedSet_size(
            stHash_search(alignmentWeightAdjLists, aW->rWeight->column))) {
        aW = aW->rWeight;
    }
    return mergeColumnsP(aW, columns, alignmentWeightAdjLists, alignmentWeightsOrderedByWeight);
}

stSet *getMultipleSequenceAlignment(stList *seqFrags, stList *multipleAlignedPairs, double gapGamma) {
    stSet *columns = makeColumns(seqFrags);
    stHash *alignmentWeightAdjLists = makeAlignmentWeightAdjacencyLists(columns, multipleAlignedPairs);
    stSortedSet *alignmentWeightsOrderedByWeight = makeOrderedSetOfAlignmentWeights(alignmentWeightAdjLists);
    stPosetAlignment *posetAlignment = stPosetAlignment_construct(stList_length(seqFrags));
    while (stSortedSet_size(alignmentWeightsOrderedByWeight) > 0) {
        AlignmentWeight *aW = stSortedSet_getLast(alignmentWeightsOrderedByWeight);
        if (aW->avgWeight < gapGamma) {
            break;
        }
        removeAlignmentFromSortedAlignmentWeights(aW, alignmentWeightsOrderedByWeight);
        Column *c = aW->column, *c2 = aW->rWeight->column;
        if (c->seqName != c2->seqName && stPosetAlignment_add(posetAlignment, c->seqName, c->position, c2->seqName,
                c2->position)) {
            mergeColumns(aW, columns, alignmentWeightAdjLists, alignmentWeightsOrderedByWeight);
        } else {
            stSortedSet_remove(stHash_search(alignmentWeightAdjLists, aW->rWeight->column), aW);
            stSortedSet_remove(stHash_search(alignmentWeightAdjLists, aW->column), aW->rWeight);
            alignmentWeight_destruct(aW);
        }
    }
    stSortedSet_destruct(alignmentWeightsOrderedByWeight);
    stHash_destruct(alignmentWeightAdjLists);
    stPosetAlignment_destruct(posetAlignment);
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
    AlignmentWeight *aW;
    int64_t refCount;
};

static ColumnPair *columnPair_construct(int64_t xIndex, int64_t yIndex, int64_t score, ColumnPair *pPair, AlignmentWeight *aW) {
    ColumnPair *cP = st_malloc(sizeof(ColumnPair));
    cP->xIndex = xIndex;
    cP->yIndex = yIndex;
    cP->score = score;
    cP->pPair = pPair;
    cP->aW = aW;
    cP->refCount = 1;
    if(pPair != NULL) {
        pPair->refCount++;
    }
    return cP;
}

static void columnPair_destruct(ColumnPair *cP) {
    cP->refCount--;
    assert(cP->refCount >= 0);
    if(cP->refCount == 0) {
        if (cP->pPair != NULL) {
            columnPair_destruct(cP->pPair);
        }
        free(cP);
    }
}

static int columnPair_cmpByYIndex(const void *c1, const void *c2) {
    int64_t i = ((ColumnPair *) c1)->yIndex;
    int64_t j = ((ColumnPair *) c2)->yIndex;
    return i > j ? 1 : (i < j ? -1 : 0);
}

int64_t getTotalWeights(stList *seqColumns, stHash *alignmentWeightAdjLists) {
    int64_t totalWeights = 0;
    for(int64_t i=0; i<stList_length(seqColumns); i++) {
        stSortedSet *aWs = stHash_search(alignmentWeightAdjLists, stList_get(seqColumns, i));
        if(aWs != NULL) {
            totalWeights += stSortedSet_size(aWs);
        }
    }
    return totalWeights;
}

stList *pairwiseAlignColumns(stList *seqXColumns, stList *seqYColumns, stHash *alignmentWeightAdjLists, stSet *columns,
        stSortedSet *alignmentWeightsOrderedByWeight, double gapGamma) {
    //Switch seqX and seqY if seqX has more alignment weights associated with it. This is critical to ensure linear scaling,
    //else worse case performance is quadratic
    int64_t totalXWeights = getTotalWeights(seqXColumns, alignmentWeightAdjLists);
    int64_t totalYWeights = getTotalWeights(seqYColumns, alignmentWeightAdjLists);
    if(totalXWeights > totalYWeights) {
        stList *l = seqYColumns;
        seqYColumns = seqXColumns;
        seqXColumns = l;
    }

    //Use indices of columns in list, have index --> column (obviously), but need to build column --> index.
    stHash *columnToIndexHash = stHash_construct2(NULL, (void(*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(seqYColumns); i++) {
        stHash_insert(columnToIndexHash, stList_get(seqYColumns, i), stIntTuple_construct1(i));
    }

    //Best scoring pairs
    stSortedSet *bestScoringAlignments = stSortedSet_construct3(columnPair_cmpByYIndex, (void(*)(void *)) columnPair_destruct);
    //Add in buffering first and last pairs
    ColumnPair *minPair = columnPair_construct(-1, -1, 0, NULL, NULL);
    stSortedSet_insert(bestScoringAlignments, minPair);
    stSortedSet_insert(bestScoringAlignments, columnPair_construct(stList_length(seqXColumns), stList_length(seqYColumns), INT64_MAX, minPair, NULL));

    //For each column in X.
    for (int64_t i = 0; i < stList_length(seqXColumns); i++) {
        Column *cX = stList_get(seqXColumns, i);
        //For each weight involving column X.
        stSortedSet *aWsX = stHash_search(alignmentWeightAdjLists, cX);
        if(aWsX != NULL) {
            //We first get all the valid new column pairs.
            stList *l = stList_construct();
            stSortedSetIterator *aWXIt = stSortedSet_getIterator(aWsX);
            AlignmentWeight *aWX;
            while ((aWX = stSortedSet_getPrevious(aWXIt)) != NULL) {
                //Add pair if exceeds the gap gamma.
                if(aWX->avgWeight >= gapGamma) {
                    //Locate index of other column
                    ColumnPair cP;
                    //The column weight may point to a column not in the Y column sequence, if so ignore.
                    if(stHash_search(columnToIndexHash, aWX->column) != NULL) {
                        cP.yIndex = stIntTuple_get(stHash_search(columnToIndexHash, aWX->column), 0);
                        //Search for highest scoring point up to but less than that index.
                        ColumnPair *cPP = stSortedSet_searchLessThan(bestScoringAlignments, &cP);
                        assert(cPP != NULL);
                        assert(i - cPP->xIndex > 0);
                        assert(cP.yIndex - cPP->yIndex > 0);
                        stList_append(l, columnPair_construct(i, cP.yIndex, /*new score */ cPP->score + aWX->avgWeight * aWX->numberOfWeights, cPP, aWX)); //Make first to increase ref-count of previous position.
                    }
                }
            }
            stSortedSet_destructIterator(aWXIt);
            //We now work through the new column pairs, from right-to-left along Y.
            stList_sort(l, columnPair_cmpByYIndex);
            while(stList_length(l) > 0) {
                ColumnPair *cP = stList_pop(l);
                //Find point that is equal or to the right of cP->yIndex
                ColumnPair *cPN = stSortedSet_searchGreaterThanOrEqual(bestScoringAlignments, cP);
                assert(cPN != NULL);
                if (cP->score >= cPN->score || cPN->yIndex > cP->yIndex) {
                    //Remove points that overlap or are to the right that score more poorly and clean them up.
                    while (cP->score >= cPN->score) {
                        ColumnPair *cPNN = stSortedSet_searchGreaterThan(bestScoringAlignments, cPN);
                        assert(cPNN != NULL);
                        stSortedSet_remove(bestScoringAlignments, cPN);
                        columnPair_destruct(cPN);
                        cPN = cPNN;
                    }
                    //Insert new point.
                    assert(stSortedSet_search(bestScoringAlignments, cP) == NULL);
                    stSortedSet_insert(bestScoringAlignments, cP);
                }
                else { //The new cP is redundant.
                    columnPair_destruct(cP);
                }
            }
            stList_destruct(l);
        }
    }

    //Link the right-most Y pair to the next-right most.
    ColumnPair *maxPair = stSortedSet_getLast(bestScoringAlignments);
    assert(maxPair != NULL);
    stSortedSet_remove(bestScoringAlignments, maxPair);
    assert(stSortedSet_getLast(bestScoringAlignments) != NULL);
    maxPair->pPair = stSortedSet_getLast(bestScoringAlignments);
    maxPair->pPair->refCount++;

    //Now traceback from highest scoring point to generate the alignment
    ColumnPair *cP = maxPair;
    stList *alignment = stList_construct();
    //For each alignment pair
    int64_t merges = 0;
    while (1) {
        assert(cP->pPair != NULL);
        //Add any unaligned Y columns
        assert(cP->yIndex > cP->pPair->yIndex);
        while (--cP->yIndex > cP->pPair->yIndex) {
            stList_append(alignment, stList_get(seqYColumns, cP->yIndex));
        }
        //Add any unaligned X columns
        assert(cP->xIndex > cP->pPair->xIndex);
        while (--cP->xIndex > cP->pPair->xIndex) {
            stList_append(alignment, stList_get(seqXColumns, cP->xIndex));
        }
        //Now move to previous pair
        cP = cP->pPair;
        //If this is the final pair we're done
        if (cP == minPair) {
            break;
        }
        //Merge two columns.
        Column *mergedColumn = mergeColumns(cP->aW, columns, alignmentWeightAdjLists, alignmentWeightsOrderedByWeight);
        merges++;
        //Add to array.
        stList_append(alignment, mergedColumn);
    }
    assert(stList_length(alignment) + merges == stList_length(seqXColumns) + stList_length(seqYColumns));
    //Make the list of columns left-to-right
    stList_reverse(alignment);

    //Cleanup
    assert(maxPair->refCount == 1);
    columnPair_destruct(maxPair);
    stSortedSet_destruct(bestScoringAlignments);
    assert(minPair->refCount == 1);
    columnPair_destruct(minPair);
    stHash_destruct(columnToIndexHash);
    stList_destruct(seqXColumns);
    stList_destruct(seqYColumns);

    return alignment;
}

stList *makeColumnSequences(stList *seqFrags, stSet *columns) {
    /*
     * Converts each seqFrag into a sequence of columns.
     */
    stList *columnSequences = stList_construct();
    for (int64_t seq = 0; seq < stList_length(seqFrags); seq++) {
        int64_t seqLength = ((SeqFrag *) (stList_get(seqFrags, seq)))->length;
        assert(seqLength >= 0);
        stList *columnSequence = stList_construct();
        for (int64_t pos = 0; pos < seqLength; pos++) {
            assert(getColumn(columns, seq, pos) != NULL);
            stList_append(columnSequence, getColumn(columns, seq, pos));
        }
        stList_append(columnSequences, columnSequence);
    }
    return columnSequences;
}

stSet *getMultipleSequenceAlignmentProgressive(stList *seqFrags, stList *multipleAlignedPairs, double gapGamma, stList *seqPairSimilarityScores) {
    //Get the data-structures needed for the pairwise alignments
    stSet *columns = makeColumns(seqFrags);
    stHash *alignmentWeightAdjLists = makeAlignmentWeightAdjacencyLists(columns, multipleAlignedPairs);
    stSortedSet *alignmentWeightsOrderedByWeight = makeOrderedSetOfAlignmentWeights(alignmentWeightAdjLists);

    //sort list of pairwise distances
    seqPairSimilarityScores = stList_copy(seqPairSimilarityScores, NULL);
    stList_sort(seqPairSimilarityScores, (int(*)(const void *, const void *)) stIntTuple_cmpFn);

    //get list of column-sequences
    stList *columnSequences = makeColumnSequences(seqFrags, columns);

    //do n-1 merges
    while (stList_length(seqPairSimilarityScores) > 0) {
        //get next closest pair
        stIntTuple *seqPair = stList_pop(seqPairSimilarityScores);

        //if not same column then align and update hash from seq indices to column-sequences
        int64_t seqX = stIntTuple_get(seqPair, 1);
        int64_t seqY = stIntTuple_get(seqPair, 2);
        stList *seqXColumns = stList_get(columnSequences, seqX);
        stList *seqYColumns = stList_get(columnSequences, seqY);
        if (seqXColumns != seqYColumns) {
            stList *seqColumns = pairwiseAlignColumns(seqXColumns, seqYColumns, alignmentWeightAdjLists, columns, alignmentWeightsOrderedByWeight, gapGamma);
            for(int64_t i=0; i<stList_length(columnSequences); i++) { //Replace instances of seqXColumns and seqYColumns with seqColumns
                stList *j = stList_get(columnSequences, i);
                if(j == seqXColumns || j == seqYColumns) {
                    stList_set(columnSequences, i, seqColumns);
                }
            }
        }
    }

    //Clean up
    stSortedSet_destruct(alignmentWeightsOrderedByWeight);
    stHash_destruct(alignmentWeightAdjLists);
    if(stList_length(columnSequences) > 0) {
        stList_destruct(stList_peek(columnSequences)); //This is because we have repeated copies of the same column-sequence in the list
    }
    stList_destruct(columnSequences);
    stList_destruct(seqPairSimilarityScores);
    //Return final set of columns.
    return columns;
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

stList *filterMultipleAlignedPairs(stSet *columns, stList *multipleAlignedPairs) {
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
        if (getColumn2(positionsToColumns, stIntTuple_get(mAP, 1), stIntTuple_get(mAP, 2)) == getColumn2(positionsToColumns,
                stIntTuple_get(mAP, 3), stIntTuple_get(mAP, 4))) {
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

static int64_t getAlignmentScore(stList *alignedPairs, int64_t seqLength1, int64_t seqLength2) {
    /*
     * Gets the normalised average posterior probability that a position in the shorter of the two sequences is aligned to a match.
     */
    int64_t alignmentScore = 0;
    for (int64_t i = 0; i < stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        alignmentScore += stIntTuple_get(alignedPair, 0);
    }
    int64_t j = seqLength1 < seqLength2 ? seqLength1 : seqLength2;
    j = j == 0 ? 1 : j;
    double d = (double) alignmentScore / (j * PAIR_ALIGNMENT_PROB_1);
    d = d > 1.0 ? 1.0 : d;
    d = d < 0.0 ? 0.0 : d;
    return d * PAIR_ALIGNMENT_PROB_1;
}

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

static int64_t addMultipleAlignedPairs(int64_t sequence1, int64_t sequence2, stList *seqFrags, stList *multipleAlignedPairs,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Computes a pairwise alignment and returns the pairwise match probabilities as tuples of (score, seq1, pos1, seq2, pos2).
     */
    SeqFrag *seqFrag1 = stList_get(seqFrags, sequence1);
    SeqFrag *seqFrag2 = stList_get(seqFrags, sequence2);
    stList *alignedPairs = getAlignedPairs(seqFrag1->seq, seqFrag2->seq, pairwiseAlignmentBandingParameters,
            seqFrag1->leftEndId != seqFrag2->leftEndId, seqFrag1->rightEndId != seqFrag2->rightEndId);
    int64_t distance = getAlignmentScore(alignedPairs, seqFrag1->length, seqFrag2->length);
    convertToMultipleAlignedPairs(alignedPairs, multipleAlignedPairs, sequence1, sequence2);
    return distance;
}

stList *makeAllPairwiseAlignments(stList *seqFrags, PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, stList **seqPairSimilarityScores) {
    /*
     * Generate the set of pairwise alignments between the sequences.
     */
    *seqPairSimilarityScores = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    stList *multipleAlignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    int64_t seqNo = stList_length(seqFrags);
    for (int64_t seq1 = 0; seq1 < seqNo; seq1++) {
        for (int64_t seq2 = seq1 + 1; seq2 < seqNo; seq2++) {
            stList_append(*seqPairSimilarityScores, stIntTuple_construct3(addMultipleAlignedPairs(seq1, seq2, seqFrags, multipleAlignedPairs, pairwiseAlignmentBandingParameters), seq1, seq2));
        }
    }
    return multipleAlignedPairs;
}

MultipleAlignment *makeAlignmentUsingAllPairs(stList *seqFrags, float gapGamma, PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Generate a multiple alignment considering all pairs of sequences.
     */
    MultipleAlignment *mA = st_calloc(1, sizeof(MultipleAlignment));
    mA->alignedPairs = makeAllPairwiseAlignments(seqFrags, pairwiseAlignmentBandingParameters, &mA->chosenPairwiseAlignments);
    mA->columns = getMultipleSequenceAlignment(seqFrags, mA->alignedPairs, gapGamma);
    mA->alignedPairs = filterMultipleAlignedPairs(mA->columns, mA->alignedPairs);
    return mA;
}

void multipleAlignment_destruct(MultipleAlignment *mA) {
    stList_destruct(mA->alignedPairs);
    stSet_destruct(mA->columns);
    stList_destruct(mA->chosenPairwiseAlignments);
    free(mA);
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
    return seq1 < seq2 ? stIntTuple_construct2(seq1, seq2) : stIntTuple_construct2(seq2, seq1);
}

stList *getReferencePairwiseAlignments(stList *seqFrags) {
    /*
     * Picks a reference sequence ensures everyone is aligned.
     */
    stList *l = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(seqFrags); i++) {
        stList_append(l, stIntTuple_construct2(((SeqFrag *) stList_get(seqFrags, i))->length, i));
    }
    stList_sort(l, (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    stList *chosenPairsOfSequencesToAlign = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    if (stList_length(l) > 0) {
        stIntTuple *k = stList_get(l, stList_length(l) / 2);
        for (int64_t i = 0; i < stList_length(l); i++) {
            stIntTuple *j = stList_get(l, i);
            if (k != j) {
                stList_append(chosenPairsOfSequencesToAlign, makePairToAlign(stIntTuple_get(j, 1), stIntTuple_get(k, 1)));
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
    stSortedSet_setDestructor(s, (void(*)(void *)) stIntTuple_destruct);
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

int64_t getNextBestPair(int64_t seq1, int64_t *distanceCounts, int64_t seqNo,
        stSortedSet *chosenPairsOfSequencesToAlign,
        bool chooseSequenceWithCommonRightEnd, stList *seqFrags) {
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
            if(!chooseSequenceWithCommonRightEnd ||
                    ((SeqFrag *)stList_get(seqFrags, seq1))->rightEndId == ((SeqFrag *)stList_get(seqFrags, seq2))->rightEndId) {
                double gain = distances[seq2] - subsPerSite(seq1, seq2, distanceCounts, seqNo);
                if (gain > maxGain || (gain == maxGain && st_random() > 0.5) /*lame attempt to reduce always picking the same seq*/) {
                    stIntTuple *pairToAlign = makePairToAlign(seq1, seq2);
                    if (stSortedSet_search(chosenPairsOfSequencesToAlign, pairToAlign) == NULL) { //So that any pair is unique
                        maxGain = gain;
                        maxGainSeq = seq2;
                    }
                    stIntTuple_destruct(pairToAlign);
                }
            }
        }
    }
    stGraph_destruct(graph);
    free(distances);
    return maxGainSeq;
}

MultipleAlignment *makeAlignment(stList *seqFrags, int64_t spanningTrees, int64_t maxPairsToConsider,
        int64_t maximumNumberOfSequencesBeforeSwitchingToFast, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    /*
     * Computes an MSA, making up to "spanningTrees"*no of seqs pairwise alignments.
     */
    int64_t seqNo = stList_length(seqFrags);
    if (spanningTrees * (seqNo - 1) >= (seqNo * (seqNo - 1)) / 2) { //Do all pairs if we can
        return makeAlignmentUsingAllPairs(seqFrags, gapGamma, pairwiseAlignmentBandingParameters);
    }
    MultipleAlignment *mA = st_calloc(1, sizeof(MultipleAlignment));
    mA->alignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct); //pairwise alignment pairs, with sequence indices
    stSortedSet *chosenPairwiseAlignmentsSet = getReferencePairwiseAlignments2(seqFrags);
    stSortedSetIterator *pairIt = stSortedSet_getIterator(chosenPairwiseAlignmentsSet);
    stIntTuple *pairToAlign;
    mA->chosenPairwiseAlignments = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    while ((pairToAlign = stSortedSet_getNext(pairIt)) != NULL) {
        int64_t seqX = stIntTuple_get(pairToAlign, 0);
        int64_t seqY = stIntTuple_get(pairToAlign, 1);
        //We get pairwise alignments, for this first alignment we filter the pairs greedily to make them consistent
        stList_append(mA->chosenPairwiseAlignments, stIntTuple_construct3(addMultipleAlignedPairs(seqX, seqY, seqFrags,
                mA->alignedPairs, pairwiseAlignmentBandingParameters),
                seqX, seqY));
    }
    stSortedSet_destructIterator(pairIt);

    int64_t iteration = 0;
    //The first alignment of multiple aligned pairs is already consistent
    while (1) {
        mA->columns = (stList_length(seqFrags) == 2 || stList_length(seqFrags) > maximumNumberOfSequencesBeforeSwitchingToFast)
                ? getMultipleSequenceAlignmentProgressive(seqFrags, mA->alignedPairs, gapGamma, mA->chosenPairwiseAlignments)
                : getMultipleSequenceAlignment(seqFrags, mA->alignedPairs, gapGamma);
        if (++iteration >= spanningTrees) {
            stSortedSet_destruct(chosenPairwiseAlignmentsSet);
            mA->alignedPairs = filterMultipleAlignedPairs(mA->columns, mA->alignedPairs);
            return mA;
        }
        int64_t *distanceCounts = getDistanceMatrix(mA->columns, seqFrags, maxPairsToConsider);
        stSet_destruct(mA->columns);
        for (int64_t seq = 0; seq < stList_length(seqFrags); seq++) {
            int64_t otherSeq = getNextBestPair(seq, distanceCounts, seqNo, chosenPairwiseAlignmentsSet, iteration == 1, seqFrags);
            if (otherSeq != INT64_MAX) {
                assert(seq != otherSeq);
                stIntTuple *pairToAlign = makePairToAlign(seq, otherSeq);
                assert(stSortedSet_search(chosenPairwiseAlignmentsSet, pairToAlign) == NULL);
                stList_append(mA->chosenPairwiseAlignments, stIntTuple_construct3(addMultipleAlignedPairs(seq, otherSeq, seqFrags, mA->alignedPairs, pairwiseAlignmentBandingParameters), seq, otherSeq));
                stSortedSet_insert(chosenPairwiseAlignmentsSet, pairToAlign);
            }
        }
        free(distanceCounts);
    }
}
