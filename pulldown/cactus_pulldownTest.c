/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "cactus.h"
#include "cactusFlower.h"
#include "cactusSegment.h"
#include "cactusSequence.h"
#include "cactusMetaSequence.h"
#include "cactusTestCommon.h"

#include "cactus_pulldown.h"

CactusDisk *cactusDisk;

void destroyExistingFlower() {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(cactusDisk);
        cactusDisk = NULL;
    }
    else {
        fprintf(stderr, "destroyExistingFlower(): No cactus disk.\n");
        exit(1);
    }
}

// Create the simplest flower possible.
Flower *buildSimpleFlower() {
    if (cactusDisk == NULL) {
        cactusDisk = testCommon_getTemporaryCactusDisk();
        Flower *flower = flower_construct(cactusDisk);
        EventTree *eventTree = eventTree_construct2(flower);
        reference_construct(flower);

        Block *block = block_construct(10, flower);
        //MetaSequence *rootMetaSequence = 
            //metaSequence_construct( 0, 10, "ACTGACTGAC", ">one",
                                    //event_getName(eventTree_getRootEvent(eventTree)), 
                                    //cactusDisk);
        //Sequence *rootSequence = sequence_construct(rootMetaSequence, flower);
        Segment *rootSegment 
            = segment_construct(block, eventTree_getRootEvent(eventTree));

        // I guess ends are getting constructed automatically?
        End *stub = end_construct(1, flower);
        Cap *stubCap = cap_construct(stub, eventTree_getRootEvent(eventTree));
        End *stub2 = end_construct(1, flower);
        Cap *stubCap2 = cap_construct(stub2, eventTree_getRootEvent(eventTree));

        Cap *_5Cap = segment_get5Cap(rootSegment);
        Cap *_3Cap = segment_get3Cap(rootSegment);

        assert(cap_getAdjacency(_5Cap) == NULL && cap_getAdjacency(_3Cap) == NULL);

        cap_makeAdjacent(stubCap, _5Cap);
        cap_makeAdjacent(stubCap2, _3Cap);

        // Root thread now goes <stub>-<segment>-<stub2>.

        return flower;
    }
    else {
        fprintf(stderr, "buildSimpleFlower(): cactus disk is not NULL (flower exists)\n");
        exit(1);
    }
}

void pulldown_simpleFlowerTest() {
    Flower *flower = buildSimpleFlower();
    assert(flower != NULL);

    Flower_CapIterator *capIt = flower_getCapIterator(flower);
    Cap *cap;
    for (cap = flower_getNextCap(capIt);
            cap != NULL; cap = flower_getNextCap(capIt)) {
        int32_t capBadness = pulldown_getCapAdjacencyBadness(cap);
        assert(capBadness == 0);
    }
    flower_destructCapIterator(capIt);

    Flower_SegmentIterator *segIt = flower_getSequenceIterator(flower);
    Segment *segment;
    for (segment = flower_getNextSegment(segIt);
            segment != NULL; flower_getNextSegment(segIt)) {
        int32_t segmentBadness = pulldown_getSegmentBadness(segment);
        assert(segmentBadness == 0);
    }
    flower_destructSegmentIterator(segIt);
    destroyExistingFlower();
}

int main (int argc, char *argv[]) {
    fprintf(stderr, "Entering pulldown test.\n");
    pulldown_simpleFlowerTest();
    fprintf(stderr, "Exiting pulldown test.\n");
    return 0;
}
