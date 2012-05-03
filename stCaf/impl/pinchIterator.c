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

stPinch *stPinchIterator_getNext(
        stPinchIterator *pinchIterator) {
    return pinchIterator->getNextAlignment(
            pinchIterator->alignmentArg);
}

void stPinchIterator_reset(
        stPinchIterator *pinchIterator) {
    pinchIterator->alignmentArg
            = pinchIterator->startAlignmentStack(
                    pinchIterator->alignmentArg);
}

void stPinchIterator_destruct(
        stPinchIterator *pinchIterator) {
    pinchIterator->destructAlignmentArg(
            pinchIterator->alignmentArg);
    free(pinchIterator);
}

void stPinchIterator_destructAlignment(stPinchIterator *pinchIterator, stPinch *stPinch) {
    if(pinchIterator->cleanupAlignment != NULL) {
        pinchIterator->cleanupAlignment(stPinch);
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

stSortedSetIterator *startAlignmentStackForAlignedPairs(stSortedSetIterator *it) {
    while(stSortedSet_getPrevious(it) != NULL);
    return it;
}

stPinchIterator *stPinchIterator_constructFromAlignedPairs(
        stSortedSet *alignedPairs, stPinch *(*getNextAlignedPairAlignment)(stSortedSetIterator *)) {
    stPinchIterator *pinchIterator = st_malloc(
            sizeof(stPinchIterator));
    pinchIterator->alignmentArg
            = stSortedSet_getIterator(alignedPairs);
    pinchIterator->getNextAlignment
            = (stPinch *(*)(void *)) getNextAlignedPairAlignment;
    pinchIterator->destructAlignmentArg
            = (void(*)(void *)) stSortedSet_destructIterator;
    pinchIterator->startAlignmentStack
            = (void *(*)(void *)) startAlignmentStackForAlignedPairs;
    pinchIterator->cleanupAlignment = NULL;
    return pinchIterator;
}
