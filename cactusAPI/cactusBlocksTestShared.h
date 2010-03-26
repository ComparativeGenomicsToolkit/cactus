#include "cactusGlobalsPrivate.h"

static NetDisk *netDisk;
static Net *net;
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
	if(netDisk != NULL) {
		netDisk_destruct(netDisk);
		testCommon_deleteTemporaryNetDisk();
		netDisk = NULL;
	}
}

static void cactusBlocksTestSharedSetup() {
	cactusBlocksTestSharedTeardown();
	netDisk = netDisk_construct(testCommon_getTemporaryNetDisk());
	net = net_construct(netDisk);

	rootMetaEvent = metaEvent_construct("ROOT", netDisk);
	leafMetaEvent = metaEvent_construct("LEAF1", netDisk);

	eventTree = eventTree_construct(rootMetaEvent, net);

	rootEvent = eventTree_getRootEvent(eventTree);
	leafEvent = event_construct(leafMetaEvent, 0.2, rootEvent, eventTree);

	metaSequence = metaSequence_construct(1, 10, "ACTGACTGAC",
			">one", metaEvent_getName(leafMetaEvent), netDisk);
	sequence = sequence_construct(metaSequence, net);

	block = block_construct(3, net);
	rootSegment = segment_construct(block_getReverse(block), rootEvent);
	leaf1Segment = segment_construct2(block, 2, 1, sequence);
	leaf2Segment = segment_construct2(block_getReverse(block), 4, 0, sequence);
	segment_makeParentAndChild(rootSegment, leaf1Segment);
	segment_makeParentAndChild(rootSegment, leaf2Segment);
}
