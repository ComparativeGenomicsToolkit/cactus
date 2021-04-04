/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic block functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

BlockEndContents *block_getContents(Block *block) {
    assert(end_isBlock(block));
    return end_getBlockEndContents(block);
}

int blockConstruct_constructP(const void *o1, const void *o2) {
	return cactusMisc_nameCompare(segment_getName((Segment *)o1), segment_getName((Segment *)o2));
}

Block *block_construct(int64_t length, Flower *flower) {
    assert(flower != NULL);

    Name name = cactusDisk_getUniqueIDInterval(flower_getCactusDisk(flower), 3);

	Block *block = st_calloc(1, 6*sizeof(Block) + sizeof(BlockEndContents));
    // Bits: (0) orientation / (1) part_of_block / (2) is_block / (3) left / (4) is_attached / (5) side
    (block+0)->bits = 0x2B; // binary: 101011
    (block+1)->bits = 0xA; // binary: 001010
    (block+2)->bits = 0x7; // binary: 000111
    (block+3)->bits = 0x6; // binary: 000110
    (block+4)->bits = 0x3; // binary: 000011
    (block+5)->bits = 0x22; // binary: 100010
    block = block+2;

	block_getContents(block)->name = name;
	block_getContents(block)->length = length;
	block_getContents(block)->flower = flower;

	flower_addEnd(flower, block_get5End(block));
    flower_addEnd(flower, block_get3End(block));

	// Some non-exhaustive checks
    assert(block_getName(block) == name+1);
    assert(block_getLength(block) == length);
    assert(block_getFlower(block) == flower);
    assert(block_getReverse(block_getReverse(block)) == block);
    assert(block_get5End(block) == block-2);
    assert(block_get3End(block) == block+2);
    assert(end_getName(block_get5End(block)) == name);
    assert(end_getName(block_get3End(block)) == name+2);
    assert(end_getName(end_getReverse(block_get5End(block))) == name);
    assert(end_getName(end_getReverse(block_get3End(block))) == name+2);
    assert(!end_isAttached(block_get5End(block)));
    assert(!end_isAttached(block_get3End(block)));
    assert(end_getSide(block_get5End(block)));
    assert(!end_getSide(block_get3End(block)));
    assert(!end_getSide(end_getReverse(block_get5End(block))));
    assert(end_getSide(end_getReverse(block_get3End(block))));
    assert(flower_getBlock(flower, name+1) == block);
    assert(flower_getEnd(flower, name) == block_get5End(block));
    assert(flower_getEnd(flower, name+2) == block_get3End(block));

	return block;
}

bool block_getOrientation(Block *block) {
    assert(end_isBlock(block));
	return end_getOrientation(block);
}

Block *block_getPositiveOrientation(Block *block) {
	return block_getOrientation(block) ? block : block_getReverse(block);
}

Block *block_getReverse(Block *block) {
	return block_getOrientation(block) ? block+1 : block-1;
}

Name block_getName(Block *block) {
	return block_getContents(block)->name+1;
}

int64_t block_getLength(Block *block) {
	return block_getContents(block)->length;
}

Flower *block_getFlower(Block *block) {
	return block_getContents(block)->flower;
}

End *block_get5End(Block *block) {
    assert(end_isBlock(block));
	return block_getOrientation(block) ? block-2 : block+2;
}

End *block_get3End(Block *block) {
    assert(end_isBlock(block));
    return block_getOrientation(block) ? block+2 : block-2;
}

int64_t block_getInstanceNumber(Block *block) {
    Segment *segment = block_getContents(block)->firstSegment;
    int64_t totalSegments=0;
    while(segment != NULL) {
        totalSegments++; segment = segment_getContents(segment)->nSegment;
    }
    return totalSegments;
}

Segment *block_getInstanceP(Block *block, Segment *connectedSegment) {
	return block_getOrientation(block) || connectedSegment == NULL ? connectedSegment : segment_getReverse(connectedSegment);
}

Segment *block_getInstance(Block *block, Name name) {
    Segment *segment = block_getContents(block)->firstSegment;
    while(segment != NULL) {
        if(segment_getName(segment) == name) {
            return block_getInstanceP(block, segment);
        }
        segment = segment_getContents(segment)->nSegment; ;
    }
    return NULL;
}

Segment *block_getFirst(Block *block) {
	return block_getInstanceP(block, block_getContents(block)->firstSegment);
}

Block_InstanceIterator *block_getInstanceIterator(Block *block) {
	Block_InstanceIterator *iterator;
	iterator = st_malloc(sizeof(struct _block_instanceIterator));
	iterator->block = block;
	iterator->segment = block_getContents(block)->firstSegment;
	return iterator;
}

Segment *block_getNext(Block_InstanceIterator *iterator) {
    Segment *segment = iterator->segment;
    if(segment == NULL) {
        return NULL;
    }
    iterator->segment = segment_getContents(segment)->nSegment;
    return block_getInstanceP(iterator->block, segment);
}

void block_destructInstanceIterator(Block_InstanceIterator *iterator) {
	free(iterator);
}

Chain *block_getChain(Block *block) {
	Link *link;
	Chain *chain1, *chain2;
	Group *group = end_getGroup(block_get5End(block));
	chain1 = (group != NULL && (link = group_getLink(group)) != NULL) ? link_getChain(link) : NULL;
	group = end_getGroup(block_get3End(block));
	chain2 = (group != NULL && (link = group_getLink(group)) != NULL) ? link_getChain(link) : NULL;
	if(chain1 != NULL && chain2 != NULL) {
		assert(chain1 == chain2); //block should not be in more than one chain!
	}
	return chain1 != NULL ? chain1 : chain2;
}

Segment *block_getSegmentForEvent(Block *block, Name eventName) {
    /*
     * Get the segment for a given event.
     */
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment;
    while ((segment = block_getNext(it)) != NULL) {
        if (event_getName(segment_getEvent(segment)) == eventName) {
            block_destructInstanceIterator(it);
            return segment;
        }
    }
    block_destructInstanceIterator(it);
    return NULL;
}

void block_check(Block *block) {
	//Check is connected to flower properly
	assert(flower_getBlock(block_getFlower(block), block_getName(block)) == block_getPositiveOrientation(block));
	//Check we have actually set built blocks for the flower..
	assert(flower_builtBlocks(block_getFlower(block)));

	//Checks the two ends are block ends.
	End *_5End = block_get5End(block);
	End *_3End = block_get3End(block);
	assert(end_isBlockEnd(_5End));
	assert(end_isBlockEnd(_3End));
	assert(end_getOrientation(_5End) == block_getOrientation(block));
	assert(end_getOrientation(_3End) == block_getOrientation(block));
	assert(end_getBlock(_5End) == block);
	assert(end_getBlock(_3End) == block);
	assert(end_getSide(_5End)); //Check the sides of the ends are consistent.
	assert(!end_getSide(_3End));

	assert(block_getLength(block) > 0); //check block has non-zero length

	//Check reverse
	Block *rBlock = block_getReverse(block);
	assert(rBlock != NULL);
	assert(block_getReverse(block) == rBlock);
	assert(block_getOrientation(block) == !block_getOrientation(rBlock));
	assert(block_getLength(block) == block_getLength(rBlock));
	assert(block_get5End(block) == end_getReverse(block_get3End(rBlock)));
	assert(block_get3End(block) == end_getReverse(block_get5End(rBlock)));
	assert(block_getInstanceNumber(block) == block_getInstanceNumber(rBlock));
	if(block_getInstanceNumber(block) > 0) {
		assert(block_getFirst(block) == segment_getReverse(block_getFirst(rBlock)));
		if(block_getRootInstance(block) == NULL) {
			assert(block_getRootInstance(rBlock) == NULL);
		}
		else {
			assert(block_getRootInstance(block) == segment_getReverse(block_getRootInstance(rBlock)));
		}
	}

	//For each segment calls segment_check.
	Block_InstanceIterator *iterator = block_getInstanceIterator(block);
	Segment *segment;
	while((segment = block_getNext(iterator)) != NULL) {
		segment_check(segment);
	}
	block_destructInstanceIterator(iterator);
}

bool block_isTrivialChain(Block *block) {
	End *_5End = block_get5End(block);
	End *_3End = block_get3End(block);
	assert(end_getGroup(_5End) != NULL);
	assert(end_getGroup(_3End) != NULL);
	return group_getLink(end_getGroup(_5End)) == NULL &&
			group_getLink(end_getGroup(_3End)) == NULL;
}

/*
 * Private functions.
 */

void block_addInstance(Block *block, Segment *segment) {
    assert(end_isBlock(block));
    segment = segment_getPositiveOrientation(segment);
    assert(segment_getContents(segment)->nSegment == NULL);
    segment_getContents(segment)->nSegment = block_getContents(block)->firstSegment;
    block_getContents(block)->firstSegment = segment;
}

void block_removeInstance(Block *block, Segment *segment) {
    assert(end_isBlock(block));
    Segment **segmentP = &(block_getContents(block)->firstSegment);
    while(*segmentP != NULL) {
        if(segment_getName(segment) == segment_getName(*segmentP)) {
            (*segmentP) = segment_getContents(*segmentP)->nSegment; // Splice it out
            return;
        }
        segmentP = &(segment_getContents(*segmentP)->nSegment);
    }
}
