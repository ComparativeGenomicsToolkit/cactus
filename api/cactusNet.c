#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int net_constructSequencesP(const void *o1, const void *o2) {
	return netMisc_nameCompare(sequence_getName((Sequence *)o1), sequence_getName((Sequence *)o2));
}

int net_constructCapsP(const void *o1, const void *o2) {
	return netMisc_nameCompare(cap_getName((Cap *)o1), cap_getName((Cap *)o2));
}

int net_constructEndsP(const void *o1, const void *o2) {
	return netMisc_nameCompare(end_getName((End *)o1), end_getName((End *)o2));
}

int net_constructSegmentsP(const void *o1, const void *o2) {
	return netMisc_nameCompare(segment_getName((Segment *)o1), segment_getName((Segment *)o2));
}

int net_constructBlocksP(const void *o1, const void *o2) {
	return netMisc_nameCompare(block_getName((Block *)o1), block_getName((Block *)o2));
}

int net_constructGroupsP(const void *o1, const void *o2) {
	return netMisc_nameCompare(group_getName((Group *)o1),
			group_getName((Group *)o2));
}

int net_constructChainsP(const void *o1, const void *o2) {
	return netMisc_nameCompare(chain_getName((Chain *)o1), chain_getName((Chain *)o2));
}

int net_constructFacesP(const void *o1, const void *o2) {
	assert(o1 != NULL);
	assert(o2 != NULL);
	return o1 == o2 ? 0 : o1 > o2 ? 1 : -1;
}

Net *net_construct(NetDisk *netDisk) {
	return net_construct2(netDisk_getUniqueID(netDisk), netDisk);
}

Net *net_construct2(Name name, NetDisk *netDisk) {
	Net *net;
	net = st_malloc(sizeof(Net));

	net->name = name;

	net->sequences = stSortedSet_construct3(net_constructSequencesP, NULL);
	net->caps = stSortedSet_construct3(net_constructCapsP, NULL);
	net->ends = stSortedSet_construct3(net_constructEndsP, NULL);
	net->segments = stSortedSet_construct3(net_constructSegmentsP, NULL);
	net->blocks = stSortedSet_construct3(net_constructBlocksP, NULL);
	net->groups = stSortedSet_construct3(net_constructGroupsP, NULL);
	net->chains = stSortedSet_construct3(net_constructChainsP, NULL);
	net->faces = stSortedSet_construct3(net_constructFacesP, NULL);
	net->reference = NULL;

	net->parentNetName = NULL_NAME;
	net->netDisk = netDisk;
	net->faceIndex = 0;
	net->chainIndex = 0;

	net->builtBlocks = 0;
	net->builtFaces = 0;
	net->builtTrees = 0;

	netDisk_addNet(net->netDisk, net);

	//Do this bit last.. so the netdisk relationship is established
	net->eventTree = NULL;
	eventTree_construct2(net);

	return net;
}


void net_destruct(Net *net, int32_t recursive) {
	Net_GroupIterator *iterator;
	Sequence *sequence;
	End *end;
	Block *block;
	Group *group;
	Chain *chain;
	Net *nestedNet;

	if(recursive) {
		iterator = net_getGroupIterator(net);
		while((group = net_getNextGroup(iterator)) != NULL) {
			nestedNet = group_getNestedNet(group);
			if(nestedNet != NULL) {
				net_destruct(nestedNet, recursive);
			}
		}
		net_destructGroupIterator(iterator);
	}

	netDisk_unloadNet(net->netDisk, net);

	net_destructFaces(net);
	stSortedSet_destruct(net->faces);

	assert(net_getEventTree(net) != NULL);
	eventTree_destruct(net_getEventTree(net));

	while((sequence = net_getFirstSequence(net)) != NULL) {
		sequence_destruct(sequence);
	}
	stSortedSet_destruct(net->sequences);

	if(net_getReference(net) != NULL) {
		reference_destruct(net_getReference(net));
	}

	while((chain = net_getFirstChain(net)) != NULL) {
		chain_destruct(chain);
	}
	stSortedSet_destruct(net->chains);

	while((end = net_getFirstEnd(net)) != NULL) {
		end_destruct(end);
	}
	stSortedSet_destruct(net->caps);
	stSortedSet_destruct(net->ends);

	while((block = net_getFirstBlock(net)) != NULL) {
		block_destruct(block);
	}
	stSortedSet_destruct(net->segments);
	stSortedSet_destruct(net->blocks);

	while((group = net_getFirstGroup(net)) != NULL) {
		group_destruct(group);
	}
	stSortedSet_destruct(net->groups);



	free(net);
}

Name net_getName(Net *net) {
	return net->name;
}

NetDisk *net_getNetDisk(Net *net) {
	return net->netDisk;
}

EventTree *net_getEventTree(Net *net) {
	return net->eventTree;
}

Sequence *net_getFirstSequence(Net *net) {
	return stSortedSet_getFirst(net->sequences);
}

Sequence *net_getSequence(Net *net, Name name) {
	Sequence *sequence;
	sequence = sequence_getStaticNameWrapper(name);
	return stSortedSet_search(net->sequences, sequence);
}

int32_t net_getSequenceNumber(Net *net) {
	return stSortedSet_size(net->sequences);
}

Net_SequenceIterator *net_getSequenceIterator(Net *net) {
	return stSortedSet_getIterator(net->sequences);
}

Sequence *net_getNextSequence(Net_SequenceIterator *sequenceIterator) {
	return stSortedSet_getNext(sequenceIterator);
}

Sequence *net_getPreviousSequence(Net_SequenceIterator *sequenceIterator) {
	return stSortedSet_getPrevious(sequenceIterator);
}

Net_SequenceIterator *net_copySequenceIterator(Net_SequenceIterator *sequenceIterator) {
	return stSortedSet_copyIterator(sequenceIterator);
}

void net_destructSequenceIterator(Net_SequenceIterator *sequenceIterator) {
	stSortedSet_destructIterator(sequenceIterator);
}

Cap *net_getFirstCap(Net *net) {
	return stSortedSet_getFirst(net->caps);
}

Cap *net_getCap(Net *net, Name name) {
	Cap *cap;
	cap = cap_getStaticNameWrapper(name);
	return stSortedSet_search(net->caps, cap);
}

int32_t net_getCapNumber(Net *net) {
	return stSortedSet_size(net->caps);
}

Net_CapIterator *net_getCapIterator(Net *net) {
	return stSortedSet_getIterator(net->caps);
}

Cap *net_getNextCap(Net_CapIterator *capIterator) {
	return stSortedSet_getNext(capIterator);
}

Cap *net_getPreviousCap(Net_CapIterator *capIterator) {
	return stSortedSet_getPrevious(capIterator);
}

Net_CapIterator *net_copyCapIterator(Net_CapIterator *capIterator) {
	return stSortedSet_copyIterator(capIterator);
}

void net_destructCapIterator(Net_CapIterator *capIterator) {
	stSortedSet_destructIterator(capIterator);
}

End *net_getFirstEnd(Net *net) {
	return stSortedSet_getFirst(net->ends);
}

End *net_getEnd(Net *net, Name name) {
	End *end;
	end = end_getStaticNameWrapper(name);
	return stSortedSet_search(net->ends, end);
}

int32_t net_getEndNumber(Net *net) {
	return stSortedSet_size(net->ends);
}

int32_t net_getBlockEndNumber(Net *net) {
	return net_getBlockNumber(net) * 2;
}

int32_t net_getStubEndNumber(Net *net) {
	return net_getEndNumber(net) - net_getBlockEndNumber(net);
}

int32_t net_getFreeStubEndNumber(Net *net) {
	End *end;
	Net_EndIterator *iterator = net_getEndIterator(net);
	int32_t i = 0;
	while((end = net_getNextEnd(iterator)) != NULL) {
		if(end_isStubEnd(end) && end_isFree(end)) {
			i++;
		}
	}
	net_destructEndIterator(iterator);
	return i;
}

int32_t net_getAttachedStubEndNumber(Net *net) {
	return net_getStubEndNumber(net) - net_getFreeStubEndNumber(net);
}

Net_EndIterator *net_getEndIterator(Net *net) {
	return stSortedSet_getIterator(net->ends);
}

End *net_getNextEnd(Net_EndIterator *endIterator) {
	return stSortedSet_getNext(endIterator);
}

End *net_getPreviousEnd(Net_EndIterator *endIterator) {
	return stSortedSet_getPrevious(endIterator);
}

Net_EndIterator *net_copyEndIterator(Net_EndIterator *endIterator) {
	return stSortedSet_copyIterator(endIterator);
}

void net_destructEndIterator(Net_EndIterator *endIterator) {
	stSortedSet_destructIterator(endIterator);
}

Segment *net_getFirstSegment(Net *net) {
	return stSortedSet_getFirst(net->segments);
}

Segment *net_getSegment(Net *net, Name name) {
	Segment *segment;
	segment = segment_getStaticNameWrapper(name);
	return stSortedSet_search(net->segments, segment);
}

int32_t net_getSegmentNumber(Net *net) {
	return stSortedSet_size(net->segments);
}

Net_SegmentIterator *net_getSegmentIterator(Net *net) {
	return stSortedSet_getIterator(net->segments);
}

Segment *net_getNextSegment(Net_SegmentIterator *segmentIterator) {
	return stSortedSet_getNext(segmentIterator);
}

Segment *net_getPreviousSegment(Net_SegmentIterator *segmentIterator) {
	return stSortedSet_getPrevious(segmentIterator);
}

Net_SegmentIterator *net_copySegmentIterator(Net_SegmentIterator *segmentIterator) {
	return stSortedSet_copyIterator(segmentIterator);
}

void net_destructSegmentIterator(Net_SegmentIterator *segmentIterator) {
	stSortedSet_destructIterator(segmentIterator);
}

Block *net_getFirstBlock(Net *net) {
	return stSortedSet_getFirst(net->blocks);
}

Block *net_getBlock(Net *net, Name name) {
	Block *block;
	block = block_getStaticNameWrapper(name);
	return stSortedSet_search(net->blocks, block);
}

int32_t net_getBlockNumber(Net *net) {
	return stSortedSet_size(net->blocks);
}

Net_BlockIterator *net_getBlockIterator(Net *net) {
	return stSortedSet_getIterator(net->blocks);
}

Block *net_getNextBlock(Net_BlockIterator *blockIterator) {
	return stSortedSet_getNext(blockIterator);
}

Block *net_getPreviousBlock(Net_BlockIterator *blockIterator) {
	return stSortedSet_getPrevious(blockIterator);
}

Net_BlockIterator *net_copyBlockIterator(Net_BlockIterator *blockIterator) {
	return stSortedSet_copyIterator(blockIterator);
}

void net_destructBlockIterator(Net_BlockIterator *blockIterator) {
	stSortedSet_destructIterator(blockIterator);
}

Group *net_getFirstGroup(Net *net) {
	return stSortedSet_getFirst(net->groups);
}

Group *net_getGroup(Net *net, Name netName) {
	Group *group= group_getStaticNameWrapper(netName);
	return stSortedSet_search(net->groups, group);
}

int32_t net_getGroupNumber(Net *net) {
	return stSortedSet_size(net->groups);
}

Net_GroupIterator *net_getGroupIterator(Net *net) {
	return stSortedSet_getIterator(net->groups);
}

Group *net_getNextGroup(Net_GroupIterator *groupIterator) {
	return stSortedSet_getNext(groupIterator);
}

Group *net_getPreviousGroup(Net_GroupIterator *groupIterator) {
	return stSortedSet_getPrevious(groupIterator);
}

Net_GroupIterator *net_copyGroupIterator(Net_GroupIterator *groupIterator) {
	return stSortedSet_copyIterator(groupIterator);
}

void net_destructGroupIterator(Net_GroupIterator *groupIterator) {
	stSortedSet_destructIterator(groupIterator);
}

Group *net_getParentGroup(Net *net) {
	if(net->parentNetName == NULL_NAME) {
		return NULL;
	}
	Net *net2 = netDisk_getNet(net_getNetDisk(net), net->parentNetName);
	assert(net2 != NULL);
	return net_getGroup(net2, net_getName(net));
}

Chain *net_getFirstChain(Net *net) {
	return stSortedSet_getFirst(net->chains);
}

Chain *net_getChain(Net *net, Name name) {
	Chain *chain = chain_getStaticNameWrapper(name);
	return stSortedSet_search(net->chains, chain);
}

int32_t net_getChainNumber(Net *net) {
	return stSortedSet_size(net->chains);
}

Net_ChainIterator *net_getChainIterator(Net *net) {
	return stSortedSet_getIterator(net->chains);
}

Chain *net_getNextChain(Net_ChainIterator *chainIterator) {
	return stSortedSet_getNext(chainIterator);
}

Chain *net_getPreviousChain(Net_ChainIterator *chainIterator) {
	return stSortedSet_getPrevious(chainIterator);
}

Net_ChainIterator *net_copyChainIterator(Net_ChainIterator *chainIterator) {
	return stSortedSet_copyIterator(chainIterator);
}

void net_destructChainIterator(Net_ChainIterator *chainIterator) {
	stSortedSet_destructIterator(chainIterator);
}

Face *net_getFirstFace(Net *net) {
	return stSortedSet_getFirst(net->faces);
}

int32_t net_getFaceNumber(Net *net) {
	return stSortedSet_size(net->faces);
}

Net_FaceIterator *net_getFaceIterator(Net *net) {
	return stSortedSet_getIterator(net->faces);
}

Face *net_getNextFace(Net_FaceIterator *faceIterator) {
	return stSortedSet_getNext(faceIterator);
}

Face *net_getPreviousFace(Net_FaceIterator *faceIterator) {
	return stSortedSet_getPrevious(faceIterator);
}

Net_FaceIterator *net_copyFaceIterator(Net_FaceIterator *faceIterator) {
	return stSortedSet_copyIterator(faceIterator);
}

void net_destructFaceIterator(Net_FaceIterator *faceIterator) {
	stSortedSet_destructIterator(faceIterator);
}

int64_t net_getTotalBaseLength(Net *net) {
	/*
	 * The implementation of this fubction is very like that in group_getTotalBaseLength, with a few differences. Consider merging them.
	 */
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	int64_t totalLength = 0;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(!end_isBlockEnd(end)) {
			End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
			Cap *cap;
			while((cap = end_getNext(instanceIterator)) != NULL) {
				cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
				if(!cap_getSide(cap) && cap_getSequence(cap) != NULL) {
					Cap *cap2 = cap_getAdjacency(cap);
					assert(cap2 != NULL);
					while(end_isBlockEnd(cap_getEnd(cap2))) {
						Segment *segment = cap_getSegment(cap2);
						assert(segment != NULL);
						assert(segment_get5Cap(segment) == cap2);
						cap2 = cap_getAdjacency(segment_get3Cap(segment));
						assert(cap2 != NULL);
						assert(cap_getStrand(cap2));
						assert(cap_getSide(cap2));
					}
					assert(cap_getStrand(cap2));
					assert(cap_getSide(cap2));
					int32_t length = cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1;
					assert(length >= 0);
					totalLength += length;
				}
			}
			end_destructInstanceIterator(instanceIterator);
		}
	}
	net_destructEndIterator(endIterator);
	return totalLength;
}

Reference *net_getReference(Net *net) {
	return net->reference;
}

/*Net *net_mergeNets(Net *net1, Net *net2) {
	if(net_getParentGroup(net1) == NULL) { //We are merging two top level reconstructions!
		assert(net_getParentGroup(net2) == NULL);
		net_mergeNetsP(net1, net2);
	}
	else { //We are merging two sister nets, merge there parent nets, which in turn will merge the child nets.
		group_mergeGroups(net_getParentGroup(net1), net_getParentGroup(net2));
	}
	return net2;
}*/

void net_check(Net *net) {
	eventTree_check(net_getEventTree(net));

	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		group_check(group);
	}
	net_destructGroupIterator(groupIterator);

	Net_ChainIterator *chainIterator = net_getChainIterator(net);
	Chain *chain;
	while((chain = net_getNextChain(chainIterator)) != NULL) {
		chain_check(chain);
	}
	net_destructCapIterator(chainIterator);

	if(net_getReference(net) != NULL) {
		reference_check(net_getReference(net));
	}

	//We check built trees in here.
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		end_check(end);
		end_check(end_getReverse(end)); //We will test everything backwards also.
	}
	net_destructEndIterator(endIterator);

	if(net_builtFaces(net)) {
		Net_FaceIterator *faceIterator = net_getFaceIterator(net);
		Face *face;
		while((face = net_getNextFace(faceIterator)) != NULL) {
			face_check(face);
		}
		net_destructFaceIterator(faceIterator);
		face_checkFaces(net);
	}
	else {
		assert(net_getFaceNumber(net) == 0);
	}

	if(net_builtBlocks(net)) { //Note that a net for which the blocks are not yet built must be a leaf.
		Net_BlockIterator *blockIterator = net_getBlockIterator(net);
		Block *block;
		while((block = net_getNextBlock(blockIterator)) != NULL) {
			block_check(block);
			block_check(block_getReverse(block)); //We will test everything backwards also.
		}
		net_destructBlockIterator(blockIterator);
	}
	else {
		assert(net_isLeaf(net)); //Defensive
		assert(net_isTerminal(net)); //Checks that a net without built blocks is a leaf and does not
		//contain any blocks.
	}

	Net_SequenceIterator *sequenceIterator = net_getSequenceIterator(net);
	Sequence *sequence;
	while((sequence = net_getNextSequence(sequenceIterator)) != NULL) {
		sequence_check(sequence);
	}
	net_destructSequenceIterator(sequenceIterator);
}

bool net_builtBlocks(Net *net) {
	return net->builtBlocks;
}

void net_setBuiltBlocks(Net *net, bool b) {
	net->builtBlocks = b;
}

bool net_builtTrees(Net *net) {
	return net->builtTrees;
}

void net_setBuiltTrees(Net *net, bool b) {
	net->builtTrees = b;
}

bool net_builtFaces(Net *net) {
	return net->builtFaces;
}

void net_setBuildFaces(Net *net, bool b) {
	net->builtFaces = b;
	if(net_builtFaces(net)) {
		net_reconstructFaces(net);
	}
}

bool net_isLeaf(Net *net) {
	Group *group;
	Net_GroupIterator *iterator = net_getGroupIterator(net);
	while((group = net_getNextGroup(iterator)) != NULL) {
		if(!group_isLeaf(group)) {
			return 0;
		}
	}
	net_destructGroupIterator(iterator);
	return 1;
}

bool net_isTerminal(Net *net) {
	return net_isLeaf(net) && net_getStubEndNumber(net) == net_getEndNumber(net);
}

/*
 * Private functions
 */

void net_setEventTree(Net *net, EventTree *eventTree) {
	if(net_getEventTree(net) != NULL) {
		eventTree_destruct(net_getEventTree(net));
	}
	net->eventTree = eventTree;
}

void net_removeEventTree(Net *net, EventTree *eventTree) {
	assert(net_getEventTree(net) == eventTree);
	net->eventTree = NULL;
}

void net_addSequence(Net *net, Sequence *sequence) {
	assert(stSortedSet_search(net->sequences, sequence) == NULL);
	stSortedSet_insert(net->sequences, sequence);
}

void net_removeSequence(Net *net, Sequence *sequence) {
	assert(stSortedSet_search(net->sequences, sequence) != NULL);
	stSortedSet_remove(net->sequences, sequence);
}

void net_addCap(Net *net, Cap *cap) {
	cap = cap_getPositiveOrientation(cap);
	assert(stSortedSet_search(net->caps, cap) == NULL);
	stSortedSet_insert(net->caps, cap);
}

void net_removeCap(Net *net, Cap *cap) {
	cap = cap_getPositiveOrientation(cap);
	assert(stSortedSet_search(net->caps, cap) != NULL);
	stSortedSet_remove(net->caps, cap);
}

void net_addEnd(Net *net, End *end) {
	end = end_getPositiveOrientation(end);
	assert(stSortedSet_search(net->ends, end) == NULL);
	stSortedSet_insert(net->ends, end);
}

void net_removeEnd(Net *net, End *end) {
	end = end_getPositiveOrientation(end);
	assert(stSortedSet_search(net->ends, end) != NULL);
	stSortedSet_remove(net->ends, end);
}

void net_addSegment(Net *net, Segment *segment) {
	segment = segment_getPositiveOrientation(segment);
	assert(stSortedSet_search(net->segments, segment) == NULL);
	stSortedSet_insert(net->segments, segment);
}

void net_removeSegment(Net *net, Segment *segment) {
	segment = segment_getPositiveOrientation(segment);
	assert(stSortedSet_search(net->segments, segment) != NULL);
	stSortedSet_remove(net->segments, segment);
}

void net_addBlock(Net *net, Block *block) {
	block = block_getPositiveOrientation(block);
	assert(stSortedSet_search(net->blocks, block) == NULL);
	stSortedSet_insert(net->blocks, block);
}

void net_removeBlock(Net *net, Block *block) {
	block = block_getPositiveOrientation(block);
	assert(stSortedSet_search(net->blocks, block) != NULL);
	stSortedSet_remove(net->blocks, block);
}

void net_addChain(Net *net, Chain *chain) {
	assert(stSortedSet_search(net->chains, chain) == NULL);
	stSortedSet_insert(net->chains, chain);
}

void net_removeChain(Net *net, Chain *chain) {
	assert(stSortedSet_search(net->chains, chain) != NULL);
	stSortedSet_remove(net->chains, chain);
}

void net_addGroup(Net *net, Group *group) {
	assert(stSortedSet_search(net->groups, group) == NULL);
	stSortedSet_insert(net->groups, group);
}

void net_removeGroup(Net *net, Group *group) {
	assert(stSortedSet_search(net->groups, group) != NULL);
	stSortedSet_remove(net->groups, group);
}

void net_setParentGroup(Net *net, Group *group) {
	//assert(net->parentNetName == NULL_NAME); we can change this if merging the parent nets, so this no longer applies.
	net->parentNetName = net_getName(group_getNet(group));
}

void net_addFace(Net *net, Face *face) {
	assert(stSortedSet_search(net->faces, face) == NULL);
	stSortedSet_insert(net->faces, face);
}

void net_removeFace(Net *net, Face *face) {
	assert(stSortedSet_search(net->faces, face) != NULL);
	stSortedSet_remove(net->faces, face);
}

void net_setReference(Net *net, Reference *reference) {
	if(net_getReference(net) != NULL) {
		reference_destruct(net_getReference(net));
	}
	net->reference = reference;
}

void net_removeReference(Net *net, Reference *reference) {
	assert(net_getReference(net) == reference);
	net->reference = NULL;
}


/*void net_mergeNetsP(Net *net1, Net *net2) {
	//Check the build settings match
	assert(net_builtBlocks(net1) == net_builtBlocks(net2));
	assert(net_builtTrees(net1) == net_builtTrees(net2));
	assert(net_builtFaces(net1) == net_builtFaces(net2));

	//Transfers the events not in event tree 1 into event tree 2.
	EventTree *eventTree1 = net_getEventTree(net1);
	EventTree *eventTree2 = net_getEventTree(net2);
	EventTree_Iterator *eventIterator = eventTree_getIterator(eventTree1);
	Event *event;
	while((event = eventTree_getNext(eventIterator)) != NULL) {
		if(eventTree_getEvent(eventTree2, event_getName(event)) == NULL) {
			eventTree_addSiblingUnaryEvent(eventTree2, event);
		}
	}
	eventTree_destructIterator(eventIterator);

	//Transfers the sequences not in net2 from net1.
	Sequence *sequence;
	Net_SequenceIterator *sequenceIterator = net_getSequenceIterator(net1);
	while((sequence = net_getNextSequence(sequenceIterator)) != NULL) {
		if(net_getSequence(net2, sequence_getName(sequence)) == NULL) {
			sequence_construct(sequence_getMetaSequence(sequence), net2);
		}
	}
	net_destructSequenceIterator(sequenceIterator);

	//This is the difficult bit.. we're going to try and replace all the references
	//to the objects left in net1 to those in net2.

	//Ensure caps and segments have event in the second event tree..
	Net_CapIterator *capIterator = net_getCapIterator(net1);
	Cap *cap;
	while((cap = net_getNextCap(capIterator)) != NULL) {
		net_addCap(net2, cap);
		cap_setEvent(cap, eventTree_getEvent(eventTree2, event_getName(cap_getEvent(cap))));
		if(cap_getSequence(cap) != NULL) {
			cap_setSequence(cap, net_getSequence(net2, sequence_getName(cap_getSequence(cap))));
		}
	}
	net_destructCapIterator(capIterator);

	//Segments currently use caps to get events, but we must include them from the
	//net
	Net_SegmentIterator *segmentIterator = net_getSegmentIterator(net1);
	Segment *segment;
	while((segment = net_getNextSegment(segmentIterator)) != NULL) {
		net_addSegment(net2, segment);
	}
	net_destructSegmentIterator(segmentIterator);

	while(net_getEndNumber(net1) > 0) {
		end_setNet(net_getFirstEnd(net1), net2);
	}

	while(net_getBlockNumber(net1) > 0) {
		block_setNet(net_getFirstBlock(net1), net2);
	}

	while(net_getFaceNumber(net1) > 0) {
		face_setNet(net_getFirstFace(net1), net2);
	}

	while(net_getGroupNumber(net1) > 0) {
		group_setNet(net_getFirstGroup(net1), net2);
	}

	while(net_getChainNumber(net1) > 0) {
		chain_setNet(net_getFirstChain(net1), net2);
	}

	while(net_getReferenceNumber(net1) > 0) {
		reference_setNet(net_getFirstReference(net1), net2);
	}

	//Now destroy the first net.
	Name netName1 = net_getName(net1);
	net_destruct(net1, 0);
	//ensure net1 is not in the netdisk..
	netDisk_deleteNetFromDisk(net_getNetDisk(net2), netName1);
	assert(netDisk_getNet(net_getNetDisk(net2), netName1) == NULL);
}*/

/*
 * Serialisation functions.
 */


void net_writeBinaryRepresentation(Net *net, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Net_SequenceIterator *sequenceIterator;
	Net_EndIterator *endIterator;
	Net_BlockIterator *blockIterator;
	Net_GroupIterator *groupIterator;
	Net_ChainIterator *chainIterator;
	Sequence *sequence;
	End *end;
	Block *block;
	Group *group;
	Chain *chain;

	binaryRepresentation_writeElementType(CODE_NET, writeFn);
	binaryRepresentation_writeName(net_getName(net), writeFn);
	binaryRepresentation_writeBool(net_builtBlocks(net), writeFn);
	binaryRepresentation_writeBool(net_builtTrees(net), writeFn);
	binaryRepresentation_writeName(net->parentNetName, writeFn);

	if(net_getEventTree(net) != NULL) {
		eventTree_writeBinaryRepresentation(net_getEventTree(net), writeFn);
	}

	sequenceIterator = net_getSequenceIterator(net);
	while((sequence = net_getNextSequence(sequenceIterator)) != NULL) {
		sequence_writeBinaryRepresentation(sequence, writeFn);
	}
	net_destructSequenceIterator(sequenceIterator);

	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		end_writeBinaryRepresentation(end, writeFn);
	}
	net_destructEndIterator(endIterator);

	blockIterator = net_getBlockIterator(net);
	while((block = net_getNextBlock(blockIterator)) != NULL) {
		block_writeBinaryRepresentation(block, writeFn);
	}
	net_destructBlockIterator(blockIterator);

	if(net_getReference(net) != NULL) {
		reference_writeBinaryRepresentation(net_getReference(net), writeFn);
	}

	groupIterator = net_getGroupIterator(net);
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		group_writeBinaryRepresentation(group, writeFn);
	}
	net_destructGroupIterator(groupIterator);

	chainIterator = net_getChainIterator(net);
	while((chain = net_getNextChain(chainIterator)) != NULL) {
		chain_writeBinaryRepresentation(chain, writeFn);
	}
	net_destructChainIterator(chainIterator);

	binaryRepresentation_writeBool(net_builtFaces(net), writeFn);
	binaryRepresentation_writeElementType(CODE_NET, writeFn); //this avoids interpretting things wrong.
}

Net *net_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk) {
	Net *net = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_NET) {
		binaryRepresentation_popNextElementType(binaryString);
		net = net_construct2(binaryRepresentation_getName(binaryString), netDisk);
		net_setBuiltBlocks(net, binaryRepresentation_getBool(binaryString));
		net_setBuiltTrees(net, binaryRepresentation_getBool(binaryString));
		net->parentNetName = binaryRepresentation_getName(binaryString);
		eventTree_loadFromBinaryRepresentation(binaryString, net);
		while(sequence_loadFromBinaryRepresentation(binaryString, net) != NULL);
		while(end_loadFromBinaryRepresentation(binaryString, net) != NULL);
		while(block_loadFromBinaryRepresentation(binaryString, net) != NULL);
		reference_loadFromBinaryRepresentation(binaryString, net);
		while(group_loadFromBinaryRepresentation(binaryString, net) != NULL);
		while(chain_loadFromBinaryRepresentation(binaryString, net) != NULL);
		net_setBuildFaces(net, binaryRepresentation_getBool(binaryString));
		assert(binaryRepresentation_popNextElementType(binaryString) == CODE_NET);
	}
	return net;
}
