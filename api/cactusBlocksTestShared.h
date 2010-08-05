#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
static Flower *flower;
static EventTree *eventTree;
static MetaSequence *metaSequence;
static Sequence *sequence;
static MetaEvent *rootMetaEvent;
static MetaEvent *leafMetaEvent;
static Event *rootEvent;
static Event *leafEvent;

static Block *block;
static Segment *rootSegment;
static Segment *leaf1Segment;
static Segment *leaf2Segment;

static void cactusBlocksTestSharedTeardown() {
    if (cactusDisk != NULL) {
        cactusDisk_destruct(cactusDisk);
        testCommon_deleteTemporaryNetDisk();
        cactusDisk = NULL;
    }
}

static void cactusBlocksTestSharedSetup() {
    cactusBlocksTestSharedTeardown();
    cactusDisk = cactusDisk_construct(testCommon_getTemporaryNetDisk());
    flower = flower_construct(cactusDisk);

    rootMetaEvent = metaEvent_construct("ROOT", cactusDisk);
    leafMetaEvent = metaEvent_construct("LEAF1", cactusDisk);

    eventTree = eventTree_construct(rootMetaEvent, flower);

    rootEvent = eventTree_getRootEvent(eventTree);
    leafEvent = event_construct(leafMetaEvent, 0.2, rootEvent, eventTree);

    metaSequence = metaSequence_construct(1, 10, "ACTGACTGAC", ">one",
            metaEvent_getName(leafMetaEvent), cactusDisk);
    sequence = sequence_construct(metaSequence, flower);

    block = block_construct(3, flower);
    rootSegment = segment_construct(block_getReverse(block), rootEvent);
    leaf1Segment = segment_construct2(block, 2, 1, sequence);
    leaf2Segment = segment_construct2(block_getReverse(block), 4, 0, sequence);
    segment_makeParentAndChild(rootSegment, leaf1Segment);
    segment_makeParentAndChild(rootSegment, leaf2Segment);
    block_setRootInstance(block, rootSegment);
}
