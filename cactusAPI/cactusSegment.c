#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic segment functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Segment *segment_construct(Block *block, Event *event) {
	return segment_construct3(
			netDisk_getUniqueID(net_getNetDisk(block_getNet(block))),
			block,
			cap_construct(block_getLeftEnd(block), event),
			cap_construct(block_getRightEnd(block), event));
}

Segment *segment_construct2(Block *block,
		int32_t startCoordinate, bool strand, Sequence *sequence) {
	assert(startCoordinate >= sequence_getStart(sequence));
	assert(startCoordinate + block_getLength(block) <= sequence_getStart(sequence) + sequence_getLength(sequence));

	int32_t i, j;
	i = startCoordinate;
	j = startCoordinate + block_getLength(block) - 1;
	if(!strand) {
		i = j;
		j = startCoordinate;
	}
	return segment_construct3(
		netDisk_getUniqueID(net_getNetDisk(block_getNet(block))),
		block,
		cap_construct2(block_getLeftEnd(block),
				i, strand, 1, sequence),
		cap_construct2(block_getRightEnd(block),
					j, strand, 0, sequence));
}

Segment *segment_construct3(Name name, Block *block,
		Cap *_5Cap, Cap *_3Cap) {
	Segment *segment;
	segment = malloc(sizeof(Segment));
	segment->rInstance = malloc(sizeof(Segment));
	segment->rInstance->rInstance = segment;
	segment->name = name;
	segment->rInstance->name = name;
	segment->block = block;
	segment->rInstance->block = block_getReverse(block);
	segment->_5Cap = _5Cap;
	segment->rInstance->_5Cap = cap_getReverse(_3Cap);
	cap_setSegment(_5Cap, segment);
	cap_setSegment(_3Cap, segment);
	block_addInstance(block, segment);
	net_addSegment(block_getNet(block), segment);
	return segment;
}

void segment_destruct(Segment *segment) {
	block_removeInstance(segment_getBlock(segment), segment);
	net_removeSegment(block_getNet(segment_getBlock(segment)), segment);
	free(segment->rInstance);
	free(segment);
}

Block *segment_getBlock(Segment *segment) {
	return segment->block;
}

Name segment_getName(Segment *segment) {
	return segment->name;
}

bool segment_getOrientation(Segment *segment) {
	return block_getOrientation(segment_getBlock(segment));
}

Segment *segment_getPositiveOrientation(Segment *segment) {
	return segment_getOrientation(segment) ? segment : segment_getReverse(segment);
}

Segment *segment_getReverse(Segment *segment) {
	return segment->rInstance;
}

Event *segment_getEvent(Segment *segment) {
	return cap_getEvent(segment_get5Cap(segment));
}

int32_t segment_getStart(Segment *segment) {
	return cap_getCoordinate(segment_get5Cap(segment));
}

bool segment_getStrand(Segment *segment) {
	return cap_getStrand(segment_get5Cap(segment));
}

int32_t segment_getLength(Segment *segment) {
	return block_getLength(segment_getBlock(segment));
}

Sequence *segment_getSequence(Segment *segment) {
	return cap_getSequence(segment_get5Cap(segment));
}

char *segment_getString(Segment *segment) {
	Sequence *sequence = segment_getSequence(segment);
	return sequence == NULL ? NULL : sequence_getString(sequence,
			segment_getStart(segment_getStrand(segment) ? segment : segment_getReverse(segment)),
			segment_getLength(segment),
			segment_getStrand(segment));
}

Cap *segment_get5Cap(Segment *segment) {
	return segment->_5Cap;
}

Cap *segment_get3Cap(Segment *segment) {
	return cap_getReverse(segment->rInstance->_5Cap);
}

Segment *segment_getParent(Segment *segment) {
	Cap *cap;
	Segment *segment2;
	cap = segment_get5Cap(segment);
	while((cap = cap_getParent(cap)) != NULL) {
		if((segment2 = cap_getSegment(cap)) != NULL) {
			return segment2;
		}
	}
	return NULL;
}

int32_t segment_getChildNumber(Segment *segment) {
	return cap_getChildNumber(segment_get5Cap(segment));
}

Segment *segment_getChild(Segment *segment, int32_t index) {
	Cap *cap;
	Segment *segment2;
	cap = cap_getChild(segment_get5Cap(segment), index);
	while(cap_isAugmented(cap)) {
		assert(cap_getChildNumber(cap) == 1);
		cap = cap_getChild(cap, 0);
	}
	segment2 = cap_getSegment(cap);
	assert(segment_getOrientation(segment) == segment_getOrientation(segment2));
	return segment2;
}

void segment_makeParentAndChild(Segment *segmentParent, Segment *segmentChild) {
	segmentParent = segment_getPositiveOrientation(segmentParent);
	segmentChild = segment_getPositiveOrientation(segmentChild);
	cap_makeParentAndChild(segment_get5Cap(segmentParent), segment_get5Cap(segmentChild));
	cap_makeParentAndChild(segment_get3Cap(segmentParent), segment_get3Cap(segmentChild));
}

/*
 * Serialisation functions.
 */

void segment_writeBinaryRepresentation(Segment *segment, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	assert(segment_getOrientation(segment));
	binaryRepresentation_writeElementType(CODE_ATOM_INSTANCE, writeFn);
	binaryRepresentation_writeName(segment_getName(segment), writeFn);
	binaryRepresentation_writeName(cap_getName(segment_get5Cap(segment)), writeFn);
	binaryRepresentation_writeName(cap_getName(segment_get3Cap(segment)), writeFn);
}

Segment *segment_loadFromBinaryRepresentation(void **binaryString, Block *block) {
	Name name, _5InstanceName, _3InstanceName;
	Segment *segment;

	segment = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ATOM_INSTANCE) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		_5InstanceName = binaryRepresentation_getName(binaryString);
		_3InstanceName = binaryRepresentation_getName(binaryString);
		segment = segment_construct3(name, block,
				end_getInstance(block_getLeftEnd(block), _5InstanceName),
				end_getInstance(block_getRightEnd(block), _3InstanceName));
	}
	return segment;
}

Segment *segment_getStaticNameWrapper(Name name) {
	static Segment segment;
	segment.name = name;
	return &segment;
}

