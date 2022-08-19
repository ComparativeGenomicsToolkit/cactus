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
#include "paf.h"
#include "cactus.h"

stPinch *stPinchIterator_getNext(stPinchIterator *pinchIterator, stPinch *pinchToFillOut) {
    stPinch *pinch;
    while (1) {
        pinch = pinchIterator->getNextAlignment(pinchIterator->alignmentArg, pinchToFillOut);
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
    Paf *(*getPairwiseAlignment)(void *);
    Paf *paf;
    int64_t xCoordinate, yCoordinate, xName, yName;
    Cigar *op;
    bool freeAlignments;
} PairwiseAlignmentToPinch;

static PairwiseAlignmentToPinch *pairwiseAlignmentToPinch_construct(void *alignmentArg,
        Paf *(*getPairwiseAlignment)(void *), bool freeAlignments) {
    PairwiseAlignmentToPinch *pairwiseAlignmentToPinch = st_calloc(1, sizeof(PairwiseAlignmentToPinch));
    pairwiseAlignmentToPinch->alignmentArg = alignmentArg;
    pairwiseAlignmentToPinch->getPairwiseAlignment = getPairwiseAlignment;
    pairwiseAlignmentToPinch->freeAlignments = freeAlignments;
    return pairwiseAlignmentToPinch;
}

static stPinch *pairwiseAlignmentToPinch_getNext(PairwiseAlignmentToPinch *pA, stPinch *pinchToFillOut) {
    while (1) {
        if (pA->paf == NULL) {
            pA->paf = pA->getPairwiseAlignment(pA->alignmentArg);
            if (pA->paf == NULL) {
                return NULL;
            }
            pA->op = pA->paf->cigar;
            pA->xCoordinate = pA->paf->same_strand ? pA->paf->query_start : pA->paf->query_end;
            pA->yCoordinate = pA->paf->target_start;
            pA->xName = cactusMisc_stringToName(pA->paf->query_name);
            pA->yName = cactusMisc_stringToName(pA->paf->target_name);
        }
        while (pA->op != NULL) {
            assert(pA->op->length >= 1);
            if (pA->op->op == match) { //deal with the possibility of a zero length match (strange, but not illegal)
                if (pA->paf->same_strand) {
                    stPinch_fillOut(pinchToFillOut, pA->xName, pA->yName, pA->xCoordinate, pA->yCoordinate, pA->op->length,
                                    1);
                    pA->xCoordinate += pA->op->length;
                } else {
                    pA->xCoordinate -= pA->op->length;
                    stPinch_fillOut(pinchToFillOut, pA->xName, pA->yName, pA->xCoordinate, pA->yCoordinate, pA->op->length,
                                    0);
                }
                pA->yCoordinate += pA->op->length;
                pA->op = pA->op->next;
                return pinchToFillOut;
            }
            if (pA->op->op != query_delete) {
                pA->xCoordinate += pA->paf->same_strand ? pA->op->length : -pA->op->length;
            }
            if (pA->op->op != query_insert) {
                pA->yCoordinate += pA->op->length;
            }
            pA->op = pA->op->next;
        }
        if (pA->paf->same_strand) {
            assert(pA->xCoordinate == pA->paf->query_end);
        }
        else {
            assert(pA->xCoordinate == pA->paf->query_start);
        }
        assert(pA->yCoordinate == pA->paf->target_end);
        if (pA->freeAlignments) {
            paf_destruct(pA->paf);
        }
        pA->paf = NULL;
    }
    return NULL;
}

static PairwiseAlignmentToPinch *pairwiseAlignmentToPinch_resetForFile(PairwiseAlignmentToPinch *pA) {
    fseek(pA->alignmentArg, 0, SEEK_SET);
    pA->paf = NULL;
    return pA;
}

static void pairwiseAlignmentToPinch_destructForFile(PairwiseAlignmentToPinch *pA) {
    fclose(pA->alignmentArg);
    free(pA);
}

stPinchIterator *stPinchIterator_constructFromFile(const char *alignmentFile) {
    stPinchIterator *pinchIterator = st_calloc(1, sizeof(stPinchIterator));
    pinchIterator->alignmentArg = pairwiseAlignmentToPinch_construct(fopen(alignmentFile, "r"),
            (Paf *(*)(void *)) paf_read, 1);
    pinchIterator->getNextAlignment = (stPinch *(*)(void *, stPinch *)) pairwiseAlignmentToPinch_getNext;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) pairwiseAlignmentToPinch_destructForFile;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) pairwiseAlignmentToPinch_resetForFile;
    return pinchIterator;
}

stSortedSetIterator *startAlignmentStackForAlignedPairs(stSortedSetIterator *it) {
    while (stSortedSet_getPrevious(it) != NULL) {
        ;
    }
    return it;
}

stPinchIterator *stPinchIterator_constructFromAlignedPairs(stSortedSet *alignedPairs,
        stPinch *(*getNextAlignedPairAlignment)(stSortedSetIterator *, stPinch *)) {
    stPinchIterator *pinchIterator = st_calloc(1, sizeof(stPinchIterator));
    pinchIterator->alignmentArg = stSortedSet_getIterator(alignedPairs);
    pinchIterator->getNextAlignment = (stPinch *(*)(void *, stPinch *)) getNextAlignedPairAlignment;
    pinchIterator->destructAlignmentArg = (void(*)(void *)) stSortedSet_destructIterator;
    pinchIterator->startAlignmentStack = (void *(*)(void *)) startAlignmentStackForAlignedPairs;
    return pinchIterator;
}

void stPinchIterator_setTrim(stPinchIterator *pinchIterator, int64_t alignmentTrim) {
    pinchIterator->alignmentTrim = alignmentTrim;
}
