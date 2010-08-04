#include "CuTest.h"
#include "sonLib.h"
#include "cactus.h"
#include <stdlib.h>
#include <string.h>

//Basic net.
static CactusDisk *cactusDisk;
static Net *net;

//Event tree...
static EventTree *eventTree;
static MetaEvent *rootMetaEvent;
static MetaEvent *leafMetaEvent;
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

static void teardown() {
    if (cactusDisk != NULL) {
        cactusDisk_destruct(cactusDisk);
        testCommon_deleteTemporaryNetDisk();
        cactusDisk = NULL;
    }
}

static void setup() {
    teardown();
    cactusDisk = cactusDisk_construct(testCommon_getTemporaryNetDisk());
    net = net_construct(cactusDisk);

    //Event tree
    rootMetaEvent = metaEvent_construct("ROOT", cactusDisk);
    leafMetaEvent = metaEvent_construct("LEAF1", cactusDisk);
    eventTree = eventTree_construct(rootMetaEvent, net);
    rootEvent = eventTree_getRootEvent(eventTree);
    leafEvent = event_construct(leafMetaEvent, 0.2, rootEvent, eventTree);

    //Sequences
    metaSequence1 = metaSequence_construct(1, 10, "ACTGACTGAC", ">one",
            metaEvent_getName(leafMetaEvent), cactusDisk);
    sequence1 = sequence_construct(metaSequence1, net);

    metaSequence2 = metaSequence_construct(1, 8, "AACCGGAA", ">two",
            metaEvent_getName(leafMetaEvent), cactusDisk);
    sequence2 = sequence_construct(metaSequence2, net);

    metaSequence3 = metaSequence_construct(1, 4, "CGGG", ">three",
                metaEvent_getName(leafMetaEvent), cactusDisk);
    sequence3 = sequence_construct(metaSequence3, net);

    metaSequence4 = metaSequence_construct(1, 1, "C", ">four",
                    metaEvent_getName(leafMetaEvent), cactusDisk);
    sequence4 = sequence_construct(metaSequence4, net);

    //Ends
    end1 = end_construct2(0, 1, net);
    end2 = end_construct2(1, 1, net);
    end3 = end_construct2(1, 1, net);

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
}

