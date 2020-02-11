/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;
static Flower *flower;
static EventTree *eventTree;
static MetaSequence *metaSequence;
static Sequence *sequence;
static End *end;

static Event *rootEvent;
static Event *leafEvent;

static Cap *rootCap;
static Cap *leaf1Cap;
static Cap *leaf2Cap;
static Cap *leaf3Cap;

static void cactusEndsTestSharedTeardown(const char *testName) {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testName, cactusDisk);
        cactusDisk = NULL;
    }
}

static void cactusEndsTestSharedSetup(const char *testName) {
    cactusEndsTestSharedTeardown(testName);
    cactusDisk = testCommon_getTemporaryCactusDisk(testName);
    flower = flower_construct(cactusDisk);

    eventTree = eventTree_construct2(cactusDisk);

    rootEvent = eventTree_getRootEvent(eventTree);
    leafEvent = event_construct3("LEAF2", 0.2, rootEvent, eventTree);

    metaSequence = metaSequence_construct(0, 10, "ACTGACTGAC", ">one",
            event_getName(leafEvent), cactusDisk);
    sequence = sequence_construct(metaSequence, flower);

    end = end_construct(1, flower);
    rootCap = cap_construct(end_getReverse(end), rootEvent);
    leaf1Cap = cap_construct2(end_getReverse(end), 4, 1, sequence);
    leaf2Cap = cap_construct2(end, 6, 0, sequence);
    leaf3Cap = cap_construct2(end_getReverse(end), 7, 0, sequence);
    cap_makeParentAndChild(rootCap, leaf1Cap);
    cap_makeParentAndChild(rootCap, leaf2Cap);
    cap_makeParentAndChild(rootCap, leaf3Cap);
    end_setRootInstance(end, rootCap);
}
