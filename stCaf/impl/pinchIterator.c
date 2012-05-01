/*
 * stPinchIterator.c
 *
 *  Created on: 21 Mar 2012
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "stPinchGraphs.h"
#include "stPinchIterator.h"
#include "cactus.h"

struct _stPinchIterator {
    void *alignmentArg;
    stPinch *(*getNextAlignment)(void *);
    void *(*startAlignmentStack)(void *);
    void (*destructAlignmentArg)(void *);
    void (*cleanupAlignment)(void *);
};

stPinch *stPinchIterator_getNext(
        stPinchIterator *stPinchIterator) {
    return stPinchIterator->getNextAlignment(
            stPinchIterator->alignmentArg);
}

void stPinchIterator_reset(
        stPinchIterator *stPinchIterator) {
    stPinchIterator->alignmentArg
            = stPinchIterator->startAlignmentStack(
                    stPinchIterator->alignmentArg);
}

void stPinchIterator_destruct(
        stPinchIterator *stPinchIterator) {
    stPinchIterator->destructAlignmentArg(
            stPinchIterator->alignmentArg);
    free(stPinchIterator);
}

void stPinchIterator_destructAlignment(stPinchIterator *stPinchIterator, stPinch *stPinch) {
    if(stPinchIterator->cleanupAlignment != NULL) {
        stPinchIterator->cleanupAlignment(stPinch);
    }
}

/*static FILE *startAlignmentStackForFile(FILE *fileHandle) {
    fseek(fileHandle, 0, SEEK_SET);
    return fileHandle;
}

stPinchIterator *stPinchIterator_constructFromFile(
        const char *alignmentFile) {
    stPinchIterator *stPinchIterator = st_malloc(
            sizeof(stPinchIterator));
    stPinchIterator->alignmentArg = fopen(alignmentFile, "r");
    stPinchIterator->getNextAlignment
            = (stPinch *(*)(void *)) cigarRead;
    stPinchIterator->destructAlignmentArg = (void(*)(void *)) fclose;
    stPinchIterator->startAlignmentStack
            = (void *(*)(void *)) startAlignmentStackForFile;
    stPinchIterator->cleanupAlignment = (void (*)(void *))destructstPinch;
    return stPinchIterator;
}

static stListIterator *startAlignmentStackForList(stListIterator *listIt) {
    while (stList_getPrevious(listIt) != NULL)
        ;
    return listIt;
}

stPinchIterator *stPinchIterator_constructFromList(
        stList *alignmentsList) {
    stPinchIterator *stPinchIterator = st_malloc(
            sizeof(stPinchIterator));
    stPinchIterator->alignmentArg
            = stList_getIterator(alignmentsList);
    stPinchIterator->getNextAlignment
            = (stPinch *(*)(void *)) stList_getNext;
    stPinchIterator->destructAlignmentArg
            = (void(*)(void *)) stList_destructIterator;
    stPinchIterator->startAlignmentStack
            = (void *(*)(void *)) startAlignmentStackForList;
    stPinchIterator->cleanupAlignment = NULL;
    return stPinchIterator;
}*/

static stSortedSetIterator *startAlignmentStackForAlignedPairs(stSortedSetIterator *it) {
    while(stSortedSet_getPrevious(it) != NULL);
    return it;
}

stPinchIterator *stPinchIterator_constructFromAlignedPairs(
        stSortedSet *alignedPairs, stPinch *(*getNextAlignedPairAlignment)(stSortedSetIterator *)) {
    stPinchIterator *stPinchIterator = st_malloc(
            sizeof(stPinchIterator));
    stPinchIterator->alignmentArg
            = stSortedSet_getIterator(alignedPairs);
    stPinchIterator->getNextAlignment
            = (stPinch *(*)(void *)) getNextAlignedPairAlignment;
    stPinchIterator->destructAlignmentArg
            = (void(*)(void *)) stSortedSet_destructIterator;
    stPinchIterator->startAlignmentStack
            = (void *(*)(void *)) startAlignmentStackForAlignedPairs;
    stPinchIterator->cleanupAlignment = NULL;
    return stPinchIterator;
}
