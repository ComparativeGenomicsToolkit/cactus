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
			cap_construct5(event, block_get5End(block), 1),
			cap_construct5(event, block_get3End(block), 0));
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
		cap_construct2(block_get5End(block),
				i, strand, 1, sequence),
		cap_construct2(block_get3End(block),
					j, strand, 0, sequence));
}

Segment *segment_construct3(Name name, Block *block,
		Cap *_5Cap, Cap *_3Cap) {
	Segment *segment;
	segment = st_malloc(sizeof(Segment));
	segment->rInstance = st_malloc(sizeof(Segment));
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
	while(cap_getSegment(cap) == NULL) {
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

void segment_check(Segment *segment) {
	//Check segment is properly linked to block.
	Block *block = segment_getBlock(segment);
	assert(block_getInstance(block, segment_getName(segment)) == segment);
	//Orientations consistent.
	assert(segment_getOrientation(segment) == block_getOrientation(block));
	//Check lengths are consistent
	assert(segment_getLength(segment) == block_getLength(block));

	//Checks the two ends have caps.
	Cap *_5Cap = segment_get5Cap(segment);
	Cap *_3Cap = segment_get3Cap(segment);
	assert(_5Cap != NULL); //check segment has other ends.
	assert(_3Cap != NULL);
	assert(cap_getOtherSegmentCap(_5Cap) == _3Cap); //check we can get the other end
	assert(cap_getOtherSegmentCap(_3Cap) == _5Cap);

	//Checks the coordinates of the caps are consistent with the segment.
	assert(cap_getOrientation(_5Cap) == segment_getOrientation(segment)); //check orientations consistent
	assert(cap_getOrientation(_3Cap) == segment_getOrientation(segment));
	assert(cap_getSide(_5Cap)); //check sides correctly configured
	assert(!cap_getSide(_3Cap));
	assert(segment_getStrand(segment) == cap_getStrand(_5Cap)); //Check strand is consistent.
	assert(segment_getStrand(segment) == cap_getStrand(_3Cap));
	assert(segment_getSequence(segment) == cap_getSequence(_5Cap)); //Check sequences are common (may be null).
	assert(segment_getSequence(segment) == cap_getSequence(_3Cap));
	assert(segment_getStart(segment) == cap_getCoordinate(_5Cap)); //Check 5End coordinate is same as start, may both be int32_max.
	assert(segment_getLength(segment) == block_getLength(block)); //Check coordinate length is consistent with block
	if(segment_getStart(segment) != INT32_MAX) { //check _3End coordinate is consistent
		if(segment_getStrand(segment)) {
			assert(segment_getStart(segment) + segment_getLength(segment) - 1 == cap_getCoordinate(_3Cap));
		}
		else {
			assert(segment_getStart(segment) - segment_getLength(segment) + 1 == cap_getCoordinate(_3Cap));
		}
	}
	else {
		assert(cap_getCoordinate(_3Cap) == INT32_MAX);
	}

	//Checks the the segment has a parent, unless the root.
	if(block_getRootInstance(block) == NULL) {
		assert(segment_getParent(segment) == NULL);
	}
	else {
		if(block_getRootInstance(block) == segment) {
			assert(segment_getParent(segment) == NULL);
		}
		else { //Check the parent-child links are correct.
			Segment *ancestorSegment = segment_getParent(segment); //Check the parent / child is consistent.
			assert(ancestorSegment != NULL);
			assert(event_isAncestor(segment_getEvent(segment), segment_getEvent(ancestorSegment)));
			assert(segment_getOrientation(segment) == segment_getOrientation(ancestorSegment));

			int32_t i;
			for(i=0; i<segment_getChildNumber(segment); i++) {
				Segment *childSegment = segment_getChild(segment, i);
				assert(childSegment != NULL);
				assert(segment_getParent(childSegment) == segment);
			}
		}
	}

	//Check the reverse..
	Segment *rSegment = segment_getReverse(segment);
	assert(rSegment != NULL);
	assert(segment_getReverse(rSegment) == segment);
	assert(block == block_getReverse(segment_getBlock(rSegment)));
	assert(segment_getOrientation(segment) == !segment_getOrientation(rSegment));
	assert(segment_getName(segment) == segment_getName(rSegment));
	assert(segment_getEvent(segment) == segment_getEvent(rSegment));
	assert(segment_getSequence(segment) == segment_getSequence(rSegment));
	assert(segment_getStrand(segment) == !segment_getStrand(rSegment));
	assert(segment_getStart(segment) == cap_getCoordinate(segment_get3Cap(rSegment)));
	assert(segment_getStart(rSegment) == cap_getCoordinate(segment_get3Cap(segment)));
	assert(segment_getLength(segment) == segment_getLength(rSegment));
	assert(segment_get5Cap(segment) == cap_getReverse(segment_get3Cap(rSegment)));
	assert(segment_get3Cap(segment) == cap_getReverse(segment_get5Cap(rSegment)));
	if(segment_getParent(segment) == NULL) {
		assert(segment_getParent(rSegment) == NULL);
	}
	else {
		assert(segment_getParent(segment) == segment_getReverse(segment_getParent(rSegment)));
	}
	assert(segment_getChildNumber(segment) == segment_getChildNumber(rSegment));
	int32_t i;
	for(i=0; i<segment_getChildNumber(segment); i++) {
		assert(segment_getChild(segment, i) == segment_getReverse(segment_getChild(rSegment, i)));
	}
}

/*
 * Serialisation functions.
 */

void segment_writeBinaryRepresentation(Segment *segment, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	assert(segment_getOrientation(segment));
	binaryRepresentation_writeElementType(CODE_SEGMENT, writeFn);
	binaryRepresentation_writeName(segment_getName(segment), writeFn);
	binaryRepresentation_writeName(cap_getName(segment_get5Cap(segment)), writeFn);
	binaryRepresentation_writeName(cap_getName(segment_get3Cap(segment)), writeFn);
}

Segment *segment_loadFromBinaryRepresentation(void **binaryString, Block *block) {
	Name name, _5InstanceName, _3InstanceName;
	Segment *segment;

	segment = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_SEGMENT) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		_5InstanceName = binaryRepresentation_getName(binaryString);
		_3InstanceName = binaryRepresentation_getName(binaryString);
		segment = segment_construct3(name, block,
				end_getInstance(block_get5End(block), _5InstanceName),
				end_getInstance(block_get3End(block), _3InstanceName));
	}
	return segment;
}

Segment *segment_getStaticNameWrapper(Name name) {
	static Segment segment;
	segment.name = name;
	return &segment;
}

