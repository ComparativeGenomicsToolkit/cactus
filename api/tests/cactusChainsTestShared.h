/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
static Flower *flower;
static Flower *nestedFlower1;
static Flower *nestedFlower2;
static End *end1;
static Block *block;
static End *end2;
static Group *group1;
static Group *group2;
static Chain *chain;
static Link *link1;
static Link *link2;

static Block *block2;
static Block *block3;
static Group *group3;
static Group *group4;
static Chain *chain2;
static Link *link4;
static Segment *segment1;
static Segment *segment2;

//Circular chain
static Block *block4;
static Group *group5;
static Chain *chain3;
static Link *link5;
static Segment *segment3;


static void cactusChainsSharedTestTeardown(const char *testName) {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testName, cactusDisk);
        cactusDisk = NULL;
    }
}

static void cactusChainsSharedTestSetup(const char *testName) {
    cactusChainsSharedTestTeardown(testName);
    cactusDisk = testCommon_getTemporaryCactusDisk(testName);
    eventTree_construct2(cactusDisk);
    flower = flower_construct(cactusDisk);
    nestedFlower1 = flower_construct(cactusDisk);
    nestedFlower2 = flower_construct(cactusDisk);
    end1 = end_construct2(0, 0, flower);
    end2 = end_construct(0, flower);
    block = block_construct(2, flower);
    end_copyConstruct(end1, nestedFlower1);
    end_copyConstruct(block_get5End(block), nestedFlower1);
    end_copyConstruct(block_get3End(block), nestedFlower2);
    end_copyConstruct(end2, nestedFlower2);
    group1 = group_construct(flower, nestedFlower1);
    group2 = group_construct(flower, nestedFlower2);
    chain = chain_construct(flower);
    link1 = link_construct(end1, block_get5End(block), group1, chain);
    link2 = link_construct(block_get3End(block), end2, group2, chain);
    //Stuff for the making trivial chains
    block2 = block_construct(1, flower);
    block3 = block_construct(1, flower);
    group3 = group_construct2(flower);
    group4 = group_construct2(flower);
    chain2 = chain_construct(flower);
    end_setGroup(block_get5End(block2), group3);
    end_setGroup(block_get3End(block3), group3);
    end_setGroup(block_get3End(block2), group4);
    end_setGroup(block_get5End(block3), group4);
    link4 = link_construct(block_get3End(block2), block_get5End(block3), group4, chain2);
    Event *event = eventTree_getRootEvent(flower_getEventTree(flower));
    MetaSequence *metaSequence = metaSequence_construct(1, 2, "AA", NULL, event_getName(event), cactusDisk);
    Sequence *sequence = sequence_construct(metaSequence, flower);
    segment1 = segment_construct2(block2, 1, 1, sequence);
    segment2 = segment_construct2(block3, 2, 1, sequence);
    assert(cap_getEnd(segment_get3Cap(segment1)) == block_get3End(block2));
    assert(cap_getEnd(segment_get5Cap(segment2)) == block_get5End(block3));
    cap_makeAdjacent(segment_get3Cap(segment1), segment_get5Cap(segment2));
    group_makeNestedFlower(group4);
    flower_setBuiltBlocks(flower, 1);
    //Circular chain
    block4 = block_construct(2, flower);
    group5 = group_construct2(flower);
    chain3 = chain_construct(flower);
    end_setGroup(block_get3End(block4), group5);
    end_setGroup(block_get5End(block4), group5);
    link5 = link_construct(block_get3End(block4), block_get5End(block4), group5, chain3);
    MetaSequence *metaSequence2 = metaSequence_construct(1, 2, "CT", NULL, event_getName(event), cactusDisk);
    Sequence *sequence2 = sequence_construct(metaSequence2, flower);
    segment3 = segment_construct2(block4, 1, 1, sequence2);
}
