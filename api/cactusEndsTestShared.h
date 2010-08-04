#include "cactusGlobalsPrivate.h"

static NetDisk *netDisk = NULL;
static Net *net;
static EventTree *eventTree;
static MetaSequence *metaSequence;
static Sequence *sequence;
static End *end;

static MetaEvent *rootMetaEvent;
static MetaEvent *leafMetaEvent;

static Event *rootEvent;
static Event *leafEvent;

static Cap *rootCap;
static Cap *leaf1Cap;
static Cap *leaf2Cap;
static Cap *leaf3Cap;

static void cactusEndsTestSharedTeardown() {
    if (netDisk != NULL) {
        netDisk_destruct(netDisk);
        testCommon_deleteTemporaryNetDisk();
        netDisk = NULL;
    }
}

static void cactusEndsTestSharedSetup() {
    cactusEndsTestSharedTeardown();
    netDisk = netDisk_construct(testCommon_getTemporaryNetDisk());
    net = net_construct(netDisk);

    rootMetaEvent = metaEvent_construct("ROOT", netDisk);
    leafMetaEvent = metaEvent_construct("LEAF1", netDisk);

    eventTree = eventTree_construct(rootMetaEvent, net);

    rootEvent = eventTree_getRootEvent(eventTree);
    leafEvent = event_construct(leafMetaEvent, 0.2, rootEvent, eventTree);

    metaSequence = metaSequence_construct(0, 10, "ACTGACTGAC", ">one",
            metaEvent_getName(leafMetaEvent), netDisk);
    sequence = sequence_construct(metaSequence, net);

    end = end_construct(1, net);
    rootCap = cap_construct(end_getReverse(end), rootEvent);
    leaf1Cap = cap_construct2(end_getReverse(end), 4, 1, sequence);
    leaf2Cap = cap_construct2(end, 6, 0, sequence);
    leaf3Cap = cap_construct2(end_getReverse(end), 7, 0, sequence);
    cap_makeParentAndChild(rootCap, leaf1Cap);
    cap_makeParentAndChild(rootCap, leaf2Cap);
    cap_makeParentAndChild(rootCap, leaf3Cap);
    end_setRootInstance(end, rootCap);
}
