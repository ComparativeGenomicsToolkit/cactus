/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
static Flower *flower;
static EventTree *eventTree;
static Sequence *sequence;
static Event *rootEvent;
static Event *leafEvent;

static Block *block;
static Segment *rootSegment;
static Segment *leaf1Segment;
static Segment *leaf2Segment;

static void cactusBlocksTestSharedTeardown(const char *testName) {
    if (cactusDisk != NULL) {
        cactusDisk_destruct(cactusDisk);
        cactusDisk = NULL;
    }
}

static void cactusBlocksTestSharedSetup(const char *testName) {
    cactusBlocksTestSharedTeardown(testName);
    cactusDisk = cactusDisk_construct();
    flower = flower_construct(cactusDisk);

    eventTree = eventTree_construct2(cactusDisk);

    rootEvent = eventTree_getRootEvent(eventTree);
    leafEvent = event_construct3("LEAF1", 0.2, rootEvent, eventTree);

    sequence = sequence_construct(1, 10, "ACTGACTGAC", ">one",
            leafEvent, cactusDisk);
    flower_addSequence(flower, sequence);

    block = block_construct(3, flower);
    leaf2Segment = segment_construct2(block_getReverse(block), 4, 0, sequence);
    leaf1Segment = segment_construct2(block, 2, 1, sequence);
    rootSegment = segment_construct(block_getReverse(block), rootEvent);
}
