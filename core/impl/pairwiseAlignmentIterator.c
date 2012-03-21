/*
 * pairwiseAlignmentIterator.c
 *
 *  Created on: 21 Mar 2012
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "pairwiseAlignment.h"
#include "pairwiseAlignmentIterator.h"
#include "cactus.h"

struct _pairwiseAlignmentIterator {
    void *alignmentArg;
    struct PairwiseAlignment *(*getNextAlignment)(void *);
    void *(*startAlignmentStack)(void *);
    void (*destructAlignmentArg)(void *);
    void (*cleanupAlignment)(void *);
};

struct PairwiseAlignment *pairwiseAlignmentIterator_getNext(
        PairwiseAlignmentIterator *pairwiseAlignmentIterator) {
    return pairwiseAlignmentIterator->getNextAlignment(
            pairwiseAlignmentIterator->alignmentArg);
}

void pairwiseAlignmentIterator_reset(
        PairwiseAlignmentIterator *pairwiseAlignmentIterator) {
    pairwiseAlignmentIterator->alignmentArg
            = pairwiseAlignmentIterator->startAlignmentStack(
                    pairwiseAlignmentIterator->alignmentArg);
}

void pairwiseAlignmentIterator_destruct(
        PairwiseAlignmentIterator *pairwiseAlignmentIterator) {
    pairwiseAlignmentIterator->destructAlignmentArg(
            pairwiseAlignmentIterator->alignmentArg);
    free(pairwiseAlignmentIterator);
}

static FILE *startAlignmentStackForFile(FILE *fileHandle) {
    fseek(fileHandle, 0, SEEK_SET);
    return fileHandle;
}

void pairwiseAlignmentIterator_destructAlignment(PairwiseAlignmentIterator *pairwiseAlignmentIterator, struct PairwiseAlignment *pairwiseAlignment) {
    if(pairwiseAlignmentIterator->cleanupAlignment != NULL) {
        pairwiseAlignmentIterator->cleanupAlignment(pairwiseAlignment);
    }
}

PairwiseAlignmentIterator *pairwiseAlignmentIterator_constructFromFile(
        const char *alignmentFile) {
    PairwiseAlignmentIterator *pairwiseAlignmentIterator = st_malloc(
            sizeof(PairwiseAlignmentIterator));
    pairwiseAlignmentIterator->alignmentArg = fopen(alignmentFile, "r");
    pairwiseAlignmentIterator->getNextAlignment
            = (struct PairwiseAlignment *(*)(void *)) cigarRead;
    pairwiseAlignmentIterator->destructAlignmentArg = (void(*)(void *)) fclose;
    pairwiseAlignmentIterator->startAlignmentStack
            = (void *(*)(void *)) startAlignmentStackForFile;
    pairwiseAlignmentIterator->cleanupAlignment = (void (*)(void *))destructPairwiseAlignment;
    return pairwiseAlignmentIterator;
}

static stListIterator *startAlignmentStackForList(stListIterator *listIt) {
    while (stList_getPrevious(listIt) != NULL)
        ;
    return listIt;
}

PairwiseAlignmentIterator *pairwiseAlignmentIterator_constructFromList(
        stList *alignmentsList) {
    PairwiseAlignmentIterator *pairwiseAlignmentIterator = st_malloc(
            sizeof(PairwiseAlignmentIterator));
    pairwiseAlignmentIterator->alignmentArg
            = stList_getIterator(alignmentsList);
    pairwiseAlignmentIterator->getNextAlignment
            = (struct PairwiseAlignment *(*)(void *)) stList_getNext;
    pairwiseAlignmentIterator->destructAlignmentArg
            = (void(*)(void *)) stList_destructIterator;
    pairwiseAlignmentIterator->startAlignmentStack
            = (void *(*)(void *)) startAlignmentStackForList;
    return pairwiseAlignmentIterator;
}

static void startAlignmentStackForAlignedPairs(stSortedSetIterator *it) {
    while(stSortedSet_getPrevious(it) != NULL);
}

PairwiseAlignmentIterator *pairwiseAlignmentIterator_constructFromAlignedPairs(
        stSortedSet *alignedPairs, struct PairwiseAlignment *(*getNextAlignedPairAlignment)(stSortedSetIterator *)) {
    PairwiseAlignmentIterator *pairwiseAlignmentIterator = st_malloc(
            sizeof(PairwiseAlignmentIterator));
    pairwiseAlignmentIterator->alignmentArg
            = stSortedSet_getIterator(alignedPairs);
    pairwiseAlignmentIterator->getNextAlignment
            = (struct PairwiseAlignment *(*)(void *)) getNextAlignedPairAlignment;
    pairwiseAlignmentIterator->destructAlignmentArg
            = (void(*)(void *)) stSortedSet_destructIterator;
    pairwiseAlignmentIterator->startAlignmentStack
            = (void *(*)(void *)) startAlignmentStackForAlignedPairs;
    return pairwiseAlignmentIterator;
}
