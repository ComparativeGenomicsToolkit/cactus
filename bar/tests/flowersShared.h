/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "cactus.h"
#include <stdlib.h>
#include <string.h>
#include "pairwiseAligner.h"

//Statemachine
StateMachine *stateMachine;

//Basic flower.
static CactusDisk *cactusDisk;
static Flower *flower;

//Event tree...
static EventTree *eventTree;
static Event *rootEvent;
static Event *leafEvent;

//Sequences
static MetaSequence *metaSequence1;
static Sequence *sequence1;

static MetaSequence *metaSequence2;
static Sequence *sequence2;

static MetaSequence *metaSequence3;
static Sequence *sequence3;

static MetaSequence *metaSequence4;
static Sequence *sequence4;

//Ends
static End *end1;
static End *end2;
static End *end3;

//Caps
static Cap *cap1;
static Cap *cap2;
static Cap *cap3;
static Cap *cap4;
static Cap *cap5;
static Cap *cap6;
static Cap *cap7;
static Cap *cap8;
static Cap *cap9;
static Cap *cap10;
static Cap *cap11;
static Cap *cap12;

static PairwiseAlignmentParameters *pairwiseParameters;

static void teardown(CuTest* testCase) {
    if (cactusDisk != NULL) {
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
        cactusDisk = NULL;
        pairwiseAlignmentBandingParameters_destruct(pairwiseParameters);
        stateMachine_destruct(stateMachine);
    }
}

static void setup(CuTest* testCase) {
    teardown(testCase);
    cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
    flower = flower_construct(cactusDisk);
    stateMachine = stateMachine5_construct(fiveState);

    //Event tree
    eventTree = eventTree_construct2(cactusDisk);
    rootEvent = eventTree_getRootEvent(eventTree);
    leafEvent = event_construct3("LEAF1", 0.2, rootEvent, eventTree);

    //Sequences
    metaSequence1 = metaSequence_construct(1, 10, "ACTGACTGAC", ">one",
            event_getName(leafEvent), cactusDisk);
    sequence1 = sequence_construct(metaSequence1, flower);

    metaSequence2 = metaSequence_construct(1, 8, "AACCGGAA", ">two",
            event_getName(leafEvent), cactusDisk);
    sequence2 = sequence_construct(metaSequence2, flower);

    metaSequence3 = metaSequence_construct(1, 4, "CGGG", ">three",
                event_getName(leafEvent), cactusDisk);
    sequence3 = sequence_construct(metaSequence3, flower);

    metaSequence4 = metaSequence_construct(1, 1, "C", ">four",
                    event_getName(leafEvent), cactusDisk);
    sequence4 = sequence_construct(metaSequence4, flower);

    //Ends
    end1 = end_construct2(0, 1, flower);
    end2 = end_construct2(1, 1, flower);
    end3 = end_construct2(1, 1, flower);

    //Caps
    //ACTG, Seq 1, 1:4 (End 1 to 2)
    cap1 = cap_construct2(end1, 0, 1, sequence1);
    cap2 = cap_construct2(end2, 5, 1, sequence1);
    cap_makeAdjacent(cap1, cap2);

    //CTGAC, Seq 1, 6:10 (End 1 to 2)
    cap3 = cap_construct2(end1, 5, 1, sequence1);
    cap4 = cap_construct2(end2, 11, 1, sequence1);
    cap_makeAdjacent(cap3, cap4);

    //T, Seq 2, -8:-8 (End 1 to 3)
    cap5 = cap_construct2(end1, 9, 0, sequence2);
    cap6 = cap_construct2(end3, 7, 0, sequence2);
    cap_makeAdjacent(cap5, cap6);

    //CCGGTT, Seq 2, -6:-1 (End 3 to 1)
    cap7 = cap_construct2(end_getReverse(end3), 7, 0, sequence2);
    cap8 = cap_construct2(end_getReverse(end1), 0, 0, sequence2);
    cap_makeAdjacent(cap7, cap8);

    //CGGG, Seq 3, 1:4 (End 1 to 1) (self loop)
    cap9 = cap_construct2(end1, 0, 1, sequence3);
    cap10 = cap_construct2(end_getReverse(end1), 5, 1, sequence3);
    cap_makeAdjacent(cap9, cap10);

    //Seq 4, (between 1 and 2) (End 1 to 1) (zero length)
    cap11 = cap_construct2(end_getReverse(end2), 1, 1, sequence4);
    cap12 = cap_construct2(end3, 2, 1, sequence4);
    cap_makeAdjacent(cap11, cap12);

    pairwiseParameters = pairwiseAlignmentBandingParameters_construct();
}

