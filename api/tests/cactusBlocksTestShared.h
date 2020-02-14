/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
static Flower *flower;
static EventTree *eventTree;
static MetaSequence *metaSequence;
static Sequence *sequence;
static Event *rootEvent;
static Event *leafEvent;

static Block *block;
static Segment *rootSegment;
static Segment *leaf1Segment;
static Segment *leaf2Segment;

static void cactusBlocksTestSharedTeardown(const char *testName) {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testName, cactusDisk);
        cactusDisk = NULL;
    }
}

static void cactusBlocksTestSharedSetup(const char *testName) {
    cactusBlocksTestSharedTeardown(testName);
    cactusDisk = testCommon_getTemporaryCactusDisk(testName);
    flower = flower_construct(cactusDisk);

    eventTree = eventTree_construct2(cactusDisk);

    rootEvent = eventTree_getRootEvent(eventTree);
    leafEvent = event_construct3("LEAF1", 0.2, rootEvent, eventTree);

    metaSequence = metaSequence_construct(1, 10, "ACTGACTGAC", ">one",
            event_getName(leafEvent), cactusDisk);
    sequence = sequence_construct(metaSequence, flower);

    block = block_construct(3, flower);
    rootSegment = segment_construct(block_getReverse(block), rootEvent);
    leaf1Segment = segment_construct2(block, 2, 1, sequence);
    leaf2Segment = segment_construct2(block_getReverse(block), 4, 0, sequence);
    segment_makeParentAndChild(rootSegment, leaf1Segment);
    segment_makeParentAndChild(rootSegment, leaf2Segment);
    block_setRootInstance(block, rootSegment);
}
