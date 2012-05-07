/*
 * stPinchIterator.c
 *
 *  Created on: 21 Mar 2012
 *      Author: benedictpaten
 */

#include <stdlib.h>
#include "sonLib.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"

struct _stPinchIterator {
    void *alignmentArg;
    stPinch *(*getNextAlignment)(void *);
    void *(*startAlignmentStack)(void *);
    void (*destructAlignmentArg)(void *);
    void (*cleanupAlignment)(void *);
};

stPinch *stPinchIterator_getNext(stPinchIterator *pinchIterator) {
    return pinchIterator->getNextAlignment(pinchIterator->alignmentArg);
}

void stPinchIterator_reset(stPinchIterator *pinchIterator) {
    pinchIterator->alignmentArg = pinchIterator->startAlignmentStack(pinchIterator->alignmentArg);
}

void stPinchIterator_destruct(stPinchIterator *pinchIterator) {
    pinchIterator->destructAlignmentArg(pinchIterator->alignmentArg);
    free(pinchIterator);
}

void stPinchIterator_destructAlignment(stPinchIterator *pinchIterator, stPinch *stPinch) {
    if (pinchIterator->cleanupAlignment != NULL) {
        pinchIterator->cleanupAlignment(stPinch);
    }
}

static FILE *startAlignmentStackForFile(FILE *fileHandle) {
    fseek(fileHandle, 0, SEEK_SET);
    return fileHandle;
}

typedef struct _pairwiseAlignmentToPinch {
    FILE *alignmentFile;
    struct PairwiseAlignment *pairwiseAlignment;
    int64_t alignmentIndex, xCoordinate, yCoordinate, xName, yName;
} PairwiseAlignmentToPinch;

static void fillOutPinch(stPinch *pinch, struct AlignmentOperation *op, )

static stPinch *getNextPinchFromAlignment(PairwiseAlignmentToPinch *pA) {
    static stPinch pinch;
    if (pA->pairwiseAlignment == NULL) {
        pA->pairwiseAlignment = cigarRead(pA->alignmentFile);
        if (pA->pairwiseAlignment == NULL) {
            return NULL;
        }
        pA->alignmentIndex = 0;
        pA->xCoordinate = pA->pairwiseAlignment->start1;
        pA->yCoordinate = pA->pairwiseAlignment->start1;
        pA->xName = cactusMisc_stringToName(pA->contig1);
        pA->yName = cactusMisc_stringToName(pA->contig2);
    }
    while (pA->alignmentIndex < pA->pairwiseAlignment->operationList->length) {
        struct AlignmentOperation *op = pA->pairwiseAlignment->operationList->list[pA->alignmentIndex++];
        if (op->opType == PAIRWISE_MATCH) {
            if (op->length >= 1) { //deal with the possibility of a zero length match (strange, but not illegal)
                if (pA->strand1) {
                    stPinch_fillOut(&pinch, pA->xName, pA->yName, pA->xCoordinate, pA->yCoordinate, op->length, 1);
                    piece_recycle(&pinch, contig1, j, j + op->length - 1);
                } else {
                    piece_recycle(&piece1, contig1, -(j - 1), -(j - op->length));
                }
                if (pA->strand2) {
                    piece_recycle(&piece2, contig2, k, k + op->length - 1);
                } else {
                    piece_recycle(&piece2, contig2, -(k - 1), -(k - op->length));
                }
                addFunction(graph, &piece1, &piece2, vertexToAdjacencyComponents, extraParameter);
            }
        }
        if (op->opType != PAIRWISE_INDEL_Y) {
            pA->xCoordinate += pA->strand1 ? op->length : -op->length;
        }
        if (op->opType != PAIRWISE_INDEL_X) {
            pA->yCoordinate += pA->strand2 ? op->length : -op->length;
        }
    }
    pA->pairwiseAlignment = NULL;
    assert(pA->xCoordinate == pA->end1);
    assert(pA->yCoordinate == pA->end2);
    return getNextPinchFromAlignment(pA);
}

stPinchIterator *stPinchIterator_constructFromFile(const char *alignmentFile) {
    stPinchIterator *pinchIterator = st_malloc(sizeof(stPinchIterator));
    pinchIterator->alignmentArg = fopen(alignmentFile, "r");
    pinchIterator->getNextAlignment = (stPinch *(*)(void *)) cigarRead;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) fclose;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) startAlignmentStackForFile;
    pinchIterator->cleanupAlignment = (void(*)(void *)) destructstPinch;
    return pinchIterator;
}

static stListIterator *startAlignmentStackForList(stListIterator *listIt) {
    while (stList_getPrevious(listIt) != NULL)
        ;
    return listIt;
}

stPinchIterator *stPinchIterator_constructFromList(stList *alignmentsList) {
    stPinchIterator *pinchIterator = st_malloc(sizeof(stPinchIterator));
    pinchIterator->alignmentArg = stList_getIterator(alignmentsList);
    pinchIterator->getNextAlignment = (stPinch *(*)(void *)) stList_getNext;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) stList_destructIterator;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) startAlignmentStackForList;
    pinchIterator->cleanupAlignment = NULL;
    return stPinchIterator;
}

stSortedSetIterator *startAlignmentStackForAlignedPairs(stSortedSetIterator *it) {
    while (stSortedSet_getPrevious(it) != NULL)
        ;
    return it;
}

stPinchIterator *stPinchIterator_constructFromAlignedPairs(stSortedSet *alignedPairs,
        stPinch *(*getNextAlignedPairAlignment)(stSortedSetIterator *)) {
    stPinchIterator *pinchIterator = st_malloc(sizeof(stPinchIterator));
    pinchIterator->alignmentArg = stSortedSet_getIterator(alignedPairs);
    pinchIterator->getNextAlignment = (stPinch *(*)(void *)) getNextAlignedPairAlignment;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) stSortedSet_destructIterator;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) startAlignmentStackForAlignedPairs;
    pinchIterator->cleanupAlignment = NULL;
    return pinchIterator;
}
