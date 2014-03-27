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
#include "pairwiseAlignment.h"
#include "cactus.h"

struct _stPinchIterator {
    int64_t alignmentTrim;
    void *alignmentArg;
    stPinch *(*getNextAlignment)(void *);
    void *(*startAlignmentStack)(void *);
    void (*destructAlignmentArg)(void *);
};

stPinch *stPinchIterator_getNext(stPinchIterator *pinchIterator) {
    stPinch *pinch = NULL;
    while (1) {
        pinch = pinchIterator->getNextAlignment(pinchIterator->alignmentArg);
        if (pinch == NULL || pinchIterator->alignmentTrim <= 0) {
            break;
        }
        pinch->start1 += pinchIterator->alignmentTrim;
        pinch->start2 += pinchIterator->alignmentTrim;
        pinch->length -= 2 * pinchIterator->alignmentTrim;
        if (pinch->length > 0) {
            break;
        }
    }
    return pinch;
}

void stPinchIterator_reset(stPinchIterator *pinchIterator) {
    pinchIterator->alignmentArg = pinchIterator->startAlignmentStack(pinchIterator->alignmentArg);
}

void stPinchIterator_destruct(stPinchIterator *pinchIterator) {
    pinchIterator->destructAlignmentArg(pinchIterator->alignmentArg);
    free(pinchIterator);
}

typedef struct _pairwiseAlignmentToPinch {
    void *alignmentArg;
    struct PairwiseAlignment *(*getPairwiseAlignment)(void *);
    struct PairwiseAlignment *pairwiseAlignment;
    int64_t alignmentIndex, xCoordinate, yCoordinate, xName, yName;
    bool freeAlignments;
} PairwiseAlignmentToPinch;

static PairwiseAlignmentToPinch *pairwiseAlignmentToPinch_construct(void *alignmentArg,
        struct PairwiseAlignment *(*getPairwiseAlignment)(void *), bool freeAlignments) {
    PairwiseAlignmentToPinch *pairwiseAlignmentToPinch = st_calloc(1, sizeof(PairwiseAlignmentToPinch));
    pairwiseAlignmentToPinch->alignmentArg = alignmentArg;
    pairwiseAlignmentToPinch->getPairwiseAlignment = getPairwiseAlignment;
    pairwiseAlignmentToPinch->freeAlignments = freeAlignments;
    return pairwiseAlignmentToPinch;
}

static stPinch *pairwiseAlignmentToPinch_getNext(PairwiseAlignmentToPinch *pA) {
    static stPinch pinch;
    while (1) {
        if (pA->pairwiseAlignment == NULL) {
            pA->pairwiseAlignment = pA->getPairwiseAlignment(pA->alignmentArg);
            if (pA->pairwiseAlignment == NULL) {
                return NULL;
            }
            pA->alignmentIndex = 0;
            pA->xCoordinate = pA->pairwiseAlignment->start1;
            pA->yCoordinate = pA->pairwiseAlignment->start2;
            pA->xName = cactusMisc_stringToName(pA->pairwiseAlignment->contig1);
            pA->yName = cactusMisc_stringToName(pA->pairwiseAlignment->contig2);
        }
        while (pA->alignmentIndex < pA->pairwiseAlignment->operationList->length) {
            struct AlignmentOperation *op = pA->pairwiseAlignment->operationList->list[pA->alignmentIndex++];
            if (op->opType == PAIRWISE_MATCH && op->length >= 1) { //deal with the possibility of a zero length match (strange, but not illegal)
                if (pA->pairwiseAlignment->strand1) {
                    if (pA->pairwiseAlignment->strand2) {
                        stPinch_fillOut(&pinch, pA->xName, pA->yName, pA->xCoordinate, pA->yCoordinate, op->length, 1);
                        pA->yCoordinate += op->length;
                    } else {
                        pA->yCoordinate -= op->length;
                        stPinch_fillOut(&pinch, pA->xName, pA->yName, pA->xCoordinate, pA->yCoordinate, op->length, 0);
                    }
                    pA->xCoordinate += op->length;
                } else {
                    pA->xCoordinate -= op->length;
                    if (pA->pairwiseAlignment->strand2) {
                        stPinch_fillOut(&pinch, pA->xName, pA->yName, pA->xCoordinate, pA->yCoordinate, op->length, 0);
                        pA->yCoordinate += op->length;
                    } else {
                        pA->yCoordinate -= op->length;
                        stPinch_fillOut(&pinch, pA->xName, pA->yName, pA->xCoordinate, pA->yCoordinate, op->length, 1);
                    }
                }
                return &pinch;
            }
            if (op->opType != PAIRWISE_INDEL_Y) {
                pA->xCoordinate += pA->pairwiseAlignment->strand1 ? op->length : -op->length;
            }
            if (op->opType != PAIRWISE_INDEL_X) {
                pA->yCoordinate += pA->pairwiseAlignment->strand2 ? op->length : -op->length;
            }
        }
        assert(pA->xCoordinate == pA->pairwiseAlignment->end1);
        assert(pA->yCoordinate == pA->pairwiseAlignment->end2);
        if (pA->freeAlignments) {
            destructPairwiseAlignment(pA->pairwiseAlignment);
        }
        pA->pairwiseAlignment = NULL;
    }
    return NULL;
}

static PairwiseAlignmentToPinch *pairwiseAlignmentToPinch_resetForFile(PairwiseAlignmentToPinch *pA) {
    fseek(pA->alignmentArg, 0, SEEK_SET);
    pA->pairwiseAlignment = NULL;
    return pA;
}

static void pairwiseAlignmentToPinch_destructForFile(PairwiseAlignmentToPinch *pA) {
    fclose(pA->alignmentArg);
    free(pA);
}

stPinchIterator *stPinchIterator_constructFromFile(const char *alignmentFile) {
    stPinchIterator *pinchIterator = st_calloc(1, sizeof(stPinchIterator));
    pinchIterator->alignmentArg = pairwiseAlignmentToPinch_construct(fopen(alignmentFile, "r"),
            (struct PairwiseAlignment *(*)(void *)) cigarRead, 1);
    pinchIterator->getNextAlignment = (stPinch *(*)(void *)) pairwiseAlignmentToPinch_getNext;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) pairwiseAlignmentToPinch_destructForFile;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) pairwiseAlignmentToPinch_resetForFile;
    return pinchIterator;
}

static PairwiseAlignmentToPinch *pairwiseAlignmentToPinch_resetForList(PairwiseAlignmentToPinch *pA) {
    while (stList_getPrevious(pA->alignmentArg) != NULL)
        ;
    pA->pairwiseAlignment = NULL;
    return pA;
}

void pairwiseAlignmentToPinch_destructForList(PairwiseAlignmentToPinch *pA) {
    stList_destructIterator(pA->alignmentArg);
    free(pA);
}

stPinchIterator *stPinchIterator_constructFromList(stList *alignmentsList) {
    stPinchIterator *pinchIterator = st_calloc(1, sizeof(stPinchIterator));
    pinchIterator->alignmentArg = pairwiseAlignmentToPinch_construct(stList_getIterator(alignmentsList),
            (struct PairwiseAlignment *(*)(void *)) stList_getNext, 0);
    pinchIterator->getNextAlignment = (stPinch *(*)(void *)) pairwiseAlignmentToPinch_getNext;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) pairwiseAlignmentToPinch_destructForList;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) pairwiseAlignmentToPinch_resetForList;
    return pinchIterator;
}

stSortedSetIterator *startAlignmentStackForAlignedPairs(stSortedSetIterator *it) {
    while (stSortedSet_getPrevious(it) != NULL) {
        ;
    }
    return it;
}

stPinchIterator *stPinchIterator_constructFromAlignedPairs(stSortedSet *alignedPairs,
        stPinch *(*getNextAlignedPairAlignment)(stSortedSetIterator *)) {
    stPinchIterator *pinchIterator = st_calloc(1, sizeof(stPinchIterator));
    pinchIterator->alignmentArg = stSortedSet_getIterator(alignedPairs);
    pinchIterator->getNextAlignment = (stPinch *(*)(void *)) getNextAlignedPairAlignment;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) stSortedSet_destructIterator;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) startAlignmentStackForAlignedPairs;
    return pinchIterator;
}

void stPinchIterator_setTrim(stPinchIterator *pinchIterator, int64_t alignmentTrim) {
    pinchIterator->alignmentTrim = alignmentTrim;
}
