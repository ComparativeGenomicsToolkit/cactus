#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic net functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t net_constructSequencesP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(sequence_getName((Sequence *)o1), sequence_getName((Sequence *)o2));
}

int32_t net_constructCapsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(cap_getName((Cap *)o1), cap_getName((Cap *)o2));
}

int32_t net_constructEndsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(end_getName((End *)o1), end_getName((End *)o2));
}

int32_t net_constructSegmentsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(segment_getName((Segment *)o1), segment_getName((Segment *)o2));
}

int32_t net_constructBlocksP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(block_getName((Block *)o1), block_getName((Block *)o2));
}

int32_t net_constructGroupsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(group_getName((Group *)o1),
			group_getName((Group *)o2));
}

int32_t net_constructChainsP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(chain_getName((Chain *)o1), chain_getName((Chain *)o2));
}

int32_t net_constructFacesP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(face_getName((Face *)o1), face_getName((Face *)o2));
}

int32_t net_constructReferencesP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(reference_getName((Reference *)o1), reference_getName((Reference *)o2));
}

Net *net_construct(NetDisk *netDisk) {
	return net_construct2(netDisk_getUniqueID(netDisk), netDisk);
}

Net *net_construct2(Name name, NetDisk *netDisk) {
	Net *net;
	net = malloc(sizeof(Net));

	net->name = name;

	net->sequences = sortedSet_construct(net_constructSequencesP);
	net->caps = sortedSet_construct(net_constructCapsP);
	net->ends = sortedSet_construct(net_constructEndsP);
	net->segments = sortedSet_construct(net_constructSegmentsP);
	net->blocks = sortedSet_construct(net_constructBlocksP);
	net->groups = sortedSet_construct(net_constructGroupsP);
	net->chains = sortedSet_construct(net_constructChainsP);
	net->faces = sortedSet_construct(net_constructFacesP);
	net->references = sortedSet_construct(net_constructReferencesP);
	net->eventTree = NULL;

	net->parentNetName = NULL_NAME;
	net->netDisk = netDisk;
	net->faceIndex = 0;
	net->chainIndex = 0;

	netDisk_addNet(net->netDisk, net);

	return net;
}


void net_destruct(Net *net, int32_t recursive) {
	Net_GroupIterator *iterator;
	Sequence *sequence;
	End *end;
	Block *block;
	Group *group;
	Chain *chain;
	Face *face;
	Reference *reference;
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

	while((face = net_getFirstFace(net)) != NULL) {
		face_destruct(face);
	}
	sortedSet_destruct(net->faces, NULL);

	if(net_getEventTree(net) != NULL) {
		eventTree_destruct(net_getEventTree(net));
	}

	while((sequence = net_getFirstSequence(net)) != NULL) {
		sequence_destruct(sequence);
	}
	sortedSet_destruct(net->sequences, NULL);

	while((end = net_getFirstEnd(net)) != NULL) {
		end_destruct(end);
	}
	sortedSet_destruct(net->caps, NULL);
	sortedSet_destruct(net->ends, NULL);

	while((block = net_getFirstBlock(net)) != NULL) {
		block_destruct(block);
	}
	sortedSet_destruct(net->segments, NULL);
	sortedSet_destruct(net->blocks, NULL);

	while((group = net_getFirstGroup(net)) != NULL) {
		group_destruct(group);
	}
	sortedSet_destruct(net->groups, NULL);

	while((chain = net_getFirstChain(net)) != NULL) {
		chain_destruct(chain);
	}
	sortedSet_destruct(net->chains, NULL);

	while((reference = net_getFirstReference(net)) != NULL) {
		reference_destruct(reference);
	}
	sortedSet_destruct(net->references, NULL);

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
	return sortedSet_getFirst(net->sequences);
}

Sequence *net_getSequence(Net *net, Name name) {
	Sequence *sequence;
	sequence = sequence_getStaticNameWrapper(name);
	return sortedSet_find(net->sequences, sequence);
}

int32_t net_getSequenceNumber(Net *net) {
	return sortedSet_getLength(net->sequences);
}

Net_SequenceIterator *net_getSequenceIterator(Net *net) {
	return iterator_construct(net->sequences);
}

Sequence *net_getNextSequence(Net_SequenceIterator *sequenceIterator) {
	return iterator_getNext(sequenceIterator);
}

Sequence *net_getPreviousSequence(Net_SequenceIterator *sequenceIterator) {
	return iterator_getPrevious(sequenceIterator);
}

Net_SequenceIterator *net_copySequenceIterator(Net_SequenceIterator *sequenceIterator) {
	return iterator_copy(sequenceIterator);
}

void net_destructSequenceIterator(Net_SequenceIterator *sequenceIterator) {
	iterator_destruct(sequenceIterator);
}

Cap *net_getFirstCap(Net *net) {
	return sortedSet_getFirst(net->caps);
}

Cap *net_getCap(Net *net, Name name) {
	Cap *cap;
	cap = cap_getStaticNameWrapper(name);
	return sortedSet_find(net->caps, cap);
}

int32_t net_getCapNumber(Net *net) {
	return sortedSet_getLength(net->caps);
}

Net_CapIterator *net_getCapIterator(Net *net) {
	return iterator_construct(net->caps);
}

Cap *net_getNextCap(Net_CapIterator *capIterator) {
	return iterator_getNext(capIterator);
}

Cap *net_getPreviousCap(Net_CapIterator *capIterator) {
	return iterator_getPrevious(capIterator);
}

Net_CapIterator *net_copyCapIterator(Net_CapIterator *capIterator) {
	return iterator_copy(capIterator);
}

void net_destructCapIterator(Net_CapIterator *capIterator) {
	iterator_destruct(capIterator);
}

End *net_getFirstEnd(Net *net) {
	return sortedSet_getFirst(net->ends);
}

End *net_getEnd(Net *net, Name name) {
	End *end;
	end = end_getStaticNameWrapper(name);
	return sortedSet_find(net->ends, end);
}

int32_t net_getEndNumber(Net *net) {
	return sortedSet_getLength(net->ends);
}

Net_EndIterator *net_getEndIterator(Net *net) {
	return iterator_construct(net->ends);
}

End *net_getNextEnd(Net_EndIterator *endIterator) {
	return iterator_getNext(endIterator);
}

End *net_getPreviousEnd(Net_EndIterator *endIterator) {
	return iterator_getPrevious(endIterator);
}

Net_EndIterator *net_copyEndIterator(Net_EndIterator *endIterator) {
	return iterator_copy(endIterator);
}

void net_destructEndIterator(Net_EndIterator *endIterator) {
	iterator_destruct(endIterator);
}

Block *net_getFirstBlock(Net *net) {
	return sortedSet_getFirst(net->blocks);
}

Block *net_getBlock(Net *net, Name name) {
	Block *block;
	block = block_getStaticNameWrapper(name);
	return sortedSet_find(net->blocks, block);
}

Segment *net_getSegment(Net *net, Name completeName) {
	Segment *segment;
	segment = segment_getStaticNameWrapper(completeName);
	return sortedSet_find(net->segments, segment);
}

int32_t net_getBlockNumber(Net *net) {
	return sortedSet_getLength(net->blocks);
}

Net_BlockIterator *net_getBlockIterator(Net *net) {
	return iterator_construct(net->blocks);
}

Block *net_getNextBlock(Net_BlockIterator *blockIterator) {
	return iterator_getNext(blockIterator);
}

Block *net_getPreviousBlock(Net_BlockIterator *blockIterator) {
	return iterator_getPrevious(blockIterator);
}

Net_BlockIterator *net_copyBlockIterator(Net_BlockIterator *blockIterator) {
	return iterator_copy(blockIterator);
}

void net_destructBlockIterator(Net_BlockIterator *blockIterator) {
	iterator_destruct(blockIterator);
}

Group *net_getFirstGroup(Net *net) {
	return sortedSet_getFirst(net->groups);
}

Group *net_getGroup(Net *net, Name netName) {
	Group *group= group_getStaticNameWrapper(netName);
	return sortedSet_find(net->groups, group);
}

int32_t net_getGroupNumber(Net *net) {
	return sortedSet_getLength(net->groups);
}

Net_GroupIterator *net_getGroupIterator(Net *net) {
	return iterator_construct(net->groups);
}

Group *net_getNextGroup(Net_GroupIterator *groupIterator) {
	return iterator_getNext(groupIterator);
}

Group *net_getPreviousGroup(Net_GroupIterator *groupIterator) {
	return iterator_getPrevious(groupIterator);
}

Net_GroupIterator *net_copyGroupIterator(Net_GroupIterator *groupIterator) {
	return iterator_copy(groupIterator);
}

void net_destructGroupIterator(Net_GroupIterator *groupIterator) {
	iterator_destruct(groupIterator);
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
	return sortedSet_getFirst(net->chains);
}

Chain *net_getChain(Net *net, Name name) {
	Chain *chain = chain_getStaticNameWrapper(name);
	return sortedSet_find(net->chains, chain);
}

int32_t net_getChainNumber(Net *net) {
	return sortedSet_getLength(net->chains);
}

Net_ChainIterator *net_getChainIterator(Net *net) {
	return iterator_construct(net->chains);
}

Chain *net_getNextChain(Net_ChainIterator *chainIterator) {
	return iterator_getNext(chainIterator);
}

Chain *net_getPreviousChain(Net_ChainIterator *chainIterator) {
	return iterator_getPrevious(chainIterator);
}

Net_ChainIterator *net_copyChainIterator(Net_ChainIterator *chainIterator) {
	return iterator_copy(chainIterator);
}

void net_destructChainIterator(Net_ChainIterator *chainIterator) {
	iterator_destruct(chainIterator);
}

Face *net_getFirstFace(Net *net) {
	return sortedSet_getFirst(net->faces);
}

Face *net_getFace(Net *net, Name name) {
	Face *face;
	face = face_getStaticNameWrapper(name);
	return sortedSet_find(net->faces, face);
}

int32_t net_getFaceNumber(Net *net) {
	return sortedSet_getLength(net->faces);
}

Net_FaceIterator *net_getFaceIterator(Net *net) {
	return iterator_construct(net->faces);
}

Face *net_getNextFace(Net_FaceIterator *faceIterator) {
	return iterator_getNext(faceIterator);
}

Face *net_getPreviousFace(Net_FaceIterator *faceIterator) {
	return iterator_getPrevious(faceIterator);
}

Net_FaceIterator *net_copyFaceIterator(Net_FaceIterator *faceIterator) {
	return iterator_copy(faceIterator);
}

void net_destructFaceIterator(Net_FaceIterator *faceIterator) {
	iterator_destruct(faceIterator);
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
				if(!cap_getSide(cap)) {
					Cap *cap2 = cap_getAdjacency(cap);
					while(end_isBlockEnd(cap_getEnd(cap2))) {
						Segment *segment = cap_getSegment(cap2);
						assert(segment != NULL);
						assert(segment_get5Cap(segment) == cap2);
						cap2 = cap_getAdjacency(segment_get3Cap(segment));
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

Reference *net_getFirstReference(Net *net) {
	return sortedSet_getFirst(net->references);
}

Reference *net_getReference(Net *net, Name name) {
	Reference *reference;
	reference = reference_getStaticNameWrapper(name);
	return sortedSet_find(net->references, reference);
}

int32_t net_getReferenceNumber(Net *net) {
	return sortedSet_getLength(net->references);
}

Net_ReferenceIterator *net_getReferenceIterator(Net *net) {
	return iterator_construct(net->references);
}

Reference *net_getNextReference(Net_ReferenceIterator *referenceIterator) {
	return iterator_getNext(referenceIterator);
}

Reference *net_getPreviousReference(Net_ReferenceIterator *referenceIterator) {
	return iterator_getPrevious(referenceIterator);
}

Net_ReferenceIterator *net_copyReferenceIterator(Net_ReferenceIterator *referenceIterator) {
	return iterator_copy(referenceIterator);
}

void net_destructReferenceIterator(Net_ReferenceIterator *referenceIterator) {
	iterator_destruct(referenceIterator);
}

void net_mergeNets(Net *net1, Net *net2) {
	if(net_getParentGroup(net1) == NULL) { //We are merging two top level reconstructions!
		assert(net_getParentGroup(net2) == NULL);
		net_mergeNetsP(net1, net2);
	}
	else { //We are merging two sister nets, merge there parent nets, which in turn will merge the child nets.
		group_mergeGroups(net_getParentGroup(net1), net_getParentGroup(net2));
	}
}


/*
 * Private functions
 */

void net_addEventTree(Net *net, EventTree *eventTree) {
	net->eventTree = eventTree;
}

void net_addSequence(Net *net, Sequence *sequence) {
	assert(sortedSet_find(net->sequences, sequence) == NULL);
	sortedSet_insert(net->sequences, sequence);
}

void net_removeSequence(Net *net, Sequence *sequence) {
	assert(sortedSet_find(net->sequences, sequence) != NULL);
	sortedSet_delete(net->sequences, sequence);
}

void net_addCap(Net *net, Cap *cap) {
	cap = cap_getPositiveOrientation(cap);
	assert(sortedSet_find(net->caps, cap) == NULL);
	sortedSet_insert(net->caps, cap);
}

void net_removeCap(Net *net, Cap *cap) {
	cap = cap_getPositiveOrientation(cap);
	assert(sortedSet_find(net->caps, cap) != NULL);
	sortedSet_delete(net->caps, cap);
}

void net_addEnd(Net *net, End *end) {
	end = end_getPositiveOrientation(end);
	assert(sortedSet_find(net->ends, end) == NULL);
	sortedSet_insert(net->ends, end);
}

void net_removeEnd(Net *net, End *end) {
	end = end_getPositiveOrientation(end);
	assert(sortedSet_find(net->ends, end) != NULL);
	sortedSet_delete(net->ends, end);
}

void net_addSegment(Net *net, Segment *segment) {
	segment = segment_getPositiveOrientation(segment);
	assert(sortedSet_find(net->segments, segment) == NULL);
	sortedSet_insert(net->segments, segment);
}

void net_removeSegment(Net *net, Segment *segment) {
	segment = segment_getPositiveOrientation(segment);
	assert(sortedSet_find(net->segments, segment) != NULL);
	sortedSet_delete(net->segments, segment);
}

void net_addBlock(Net *net, Block *block) {
	block = block_getPositiveOrientation(block);
	assert(sortedSet_find(net->blocks, block) == NULL);
	sortedSet_insert(net->blocks, block);
}

void net_removeBlock(Net *net, Block *block) {
	block = block_getPositiveOrientation(block);
	assert(sortedSet_find(net->blocks, block) != NULL);
	sortedSet_delete(net->blocks, block);
}

void net_addChain(Net *net, Chain *chain) {
	assert(sortedSet_find(net->chains, chain) == NULL);
	sortedSet_insert(net->chains, chain);
}

void net_removeChain(Net *net, Chain *chain) {
	assert(sortedSet_find(net->chains, chain) != NULL);
	sortedSet_delete(net->chains, chain);
}

void net_addGroup(Net *net, Group *group) {
	assert(sortedSet_find(net->groups, group) == NULL);
	sortedSet_insert(net->groups, group);
}

void net_removeGroup(Net *net, Group *group) {
	assert(sortedSet_find(net->groups, group) != NULL);
	sortedSet_delete(net->groups, group);
}

void net_setParentGroup(Net *net, Group *group) {
	assert(net->parentNetName == NULL_NAME);
	net->parentNetName = net_getName(group_getNet(group));
}

void net_addFace(Net *net, Face *face) {
	assert(sortedSet_find(net->faces, face) == NULL);
	sortedSet_insert(net->faces, face);
}

void net_removeFace(Net *net, Face *face) {
	assert(sortedSet_find(net->faces, face) != NULL);
	sortedSet_delete(net->faces, face);
}

void net_addReference(Net *net, Reference *reference) {
	assert(sortedSet_find(net->references, reference) == NULL);
	sortedSet_insert(net->references, reference);
}

void net_removeReference(Net *net, Reference *reference) {
	assert(sortedSet_find(net->references, reference) != NULL);
	sortedSet_delete(net->references, reference);
}

void net_mergeNetsP(Net *net1, Net *net2) {
	//Make binary strings for the two nets
	//Destruct the two nets and any children that are loaded.
	//Load the net using both strings again.


	//Write the e


	Sequence *sequence;

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

	while((sequence = net_getFirstSequence(net1)) != NULL) {
		if(net_getSequence(net2, sequence_getName(sequence)) == NULL) {
			sequence_setNet(sequence, net2);
			//sequence_setEvent(sequence, eventTree_getEvent(eventTree2, event_getName(sequence_getEvent(sequence)))); //ensures it has the right event.
		}
		else {
			sequence_destruct(sequence);
		}
	}

	//Ensure caps, segments and sequences have event in the second event tree..
	Net_CapIterator *capIterator = net_getCapIterator(net1);
	Cap *cap;
	while((cap = net_getNextCap(capIterator)) != NULL) {
		//cap_setEvent(cap, eventTree_getEvent(eventTree2, event_getName(cap_getEvent(cap))));
		//cap_setSequence(cap, eventTree_getEvent(eventTree2, event_getName(cap_getEvent(cap))));
	}
	net_destructCapIterator(capIterator);
	//Segments currently use caps to get events.

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
	net_destruct(net1, 0);
}

/*
 * Serialisation functions.
 */


void net_writeBinaryRepresentation(Net *net, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Net_SequenceIterator *sequenceIterator;
	Net_EndIterator *endIterator;
	Net_BlockIterator *blockIterator;
	Net_GroupIterator *groupIterator;
	Net_ChainIterator *chainIterator;
	Net_FaceIterator *faceIterator;
	Net_ReferenceIterator *referenceIterator;
	Sequence *sequence;
	End *end;
	Block *block;
	Group *group;
	Chain *chain;
	Face *face;
	Reference *reference;

	binaryRepresentation_writeElementType(CODE_NET, writeFn);
	binaryRepresentation_writeName(net_getName(net), writeFn);
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

	faceIterator = net_getFaceIterator(net);
	while((face = net_getNextFace(faceIterator)) != NULL) {
		face_writeBinaryRepresentation(face, writeFn);
	}
	net_destructFaceIterator(faceIterator);

	referenceIterator = net_getReferenceIterator(net);
	while((reference = net_getNextReference(referenceIterator)) != NULL) {
		reference_writeBinaryRepresentation(reference, writeFn);
	}
	net_destructReferenceIterator(referenceIterator);

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

	binaryRepresentation_writeElementType(CODE_NET, writeFn); //this avoids interpretting things wrong.
}

Net *net_loadFromBinaryRepresentation(void **binaryString, NetDisk *netDisk) {
	Net *net = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_NET) {
		binaryRepresentation_popNextElementType(binaryString);
		net = net_construct2(binaryRepresentation_getName(binaryString), netDisk);
		net->parentNetName = binaryRepresentation_getName(binaryString);
		eventTree_loadFromBinaryRepresentation(binaryString, net);
		while(sequence_loadFromBinaryRepresentation(binaryString, net) != NULL);
		while(end_loadFromBinaryRepresentation(binaryString, net) != NULL);
		while(block_loadFromBinaryRepresentation(binaryString, net) != NULL);
		while(face_loadFromBinaryRepresentation(binaryString, net) != NULL);
		while(reference_loadFromBinaryRepresentation(binaryString, net) != NULL);
		while(group_loadFromBinaryRepresentation(binaryString, net) != NULL);
		while(chain_loadFromBinaryRepresentation(binaryString, net) != NULL);
		assert(binaryRepresentation_popNextElementType(binaryString) == CODE_NET);
	}
	return net;
}
