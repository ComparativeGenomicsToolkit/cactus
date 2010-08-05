#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic block functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int blockConstruct_constructP(const void *o1, const void *o2) {
	return cactusMisc_nameCompare(segment_getName((Segment *)o1), segment_getName((Segment *)o2));
}

Block *block_construct(int32_t length, Flower *flower) {
	return block_construct2(cactusDisk_getUniqueID(flower_getCactusDisk(flower)), length,
			end_construct3(cactusDisk_getUniqueID(flower_getCactusDisk(flower)), 0, 0, 1, flower),
			end_construct3(cactusDisk_getUniqueID(flower_getCactusDisk(flower)), 0, 0, 0, flower), flower);
}

Block *block_construct2(Name name, int32_t length,
		End *leftEnd, End *rightEnd,
		Flower *flower) {
	Block *block;
	block = st_malloc(sizeof(Block));
	block->rBlock = st_malloc(sizeof(Block));
	block->rBlock->rBlock = block;
	block->blockContents = st_malloc(sizeof(BlockContents));
	block->rBlock->blockContents = block->blockContents;

	block->orientation = 1;
	block->rBlock->orientation = 0;

	block->blockContents->name = name;
	block->blockContents->segments = stSortedSet_construct3(blockConstruct_constructP, NULL);
	block->blockContents->length = length;
	block->blockContents->flower = flower;

	block->leftEnd = leftEnd;
	end_setBlock(leftEnd, block);
	block->rBlock->leftEnd = end_getReverse(rightEnd);
	end_setBlock(rightEnd, block);

	flower_addBlock(flower, block);

	return block;
}

void block_destruct(Block *block) {
	Segment *segment;
	//remove from flower.
	flower_removeBlock(block_getFlower(block), block);

	//remove instances
	while((segment = block_getFirst(block)) != NULL) {
		segment_destruct(segment);
	}
	//now the actual instances.
	stSortedSet_destruct(block->blockContents->segments);

	free(block->rBlock);
	free(block->blockContents);
	free(block);
}

bool block_getOrientation(Block *block) {
	return block->orientation;
}

Block *block_getPositiveOrientation(Block *block) {
	return block_getOrientation(block) ? block : block_getReverse(block);
}

Block *block_getReverse(Block *block) {
	return block->rBlock;
}

Name block_getName(Block *block) {
	return block->blockContents->name;
}

int32_t block_getLength(Block *block) {
	return block->blockContents->length;
}

Flower *block_getFlower(Block *block) {
	return block->blockContents->flower;
}

End *block_get5End(Block *block) {
	return block->leftEnd;
}

End *block_get3End(Block *block) {
	return end_getReverse(block->rBlock->leftEnd);
}

int32_t block_getInstanceNumber(Block *block) {
	return stSortedSet_size(block->blockContents->segments);
}

Segment *block_getInstanceP(Block *block, Segment *connectedSegment) {
	return block_getOrientation(block) || connectedSegment == NULL ? connectedSegment : segment_getReverse(connectedSegment);
}

Segment *block_getInstance(Block *block, Name name) {
	Segment *segment = segment_getStaticNameWrapper(name);
	return block_getInstanceP(block, stSortedSet_search(block->blockContents->segments, segment));
}

Segment *block_getFirst(Block *block) {
	return block_getInstanceP(block, stSortedSet_getFirst(block->blockContents->segments));
}

Segment *block_getRootInstance(Block *block) {
	Cap *cap = end_getRootInstance(block_get5End(block));
	return cap != NULL ? cap_getSegment(cap) : NULL;
}

void block_setRootInstance(Block *block, Segment *segment) {
	block = block_getPositiveOrientation(block);
	segment = segment_getPositiveOrientation(segment);
	assert(block_getInstance(block, segment_getName(segment)) == segment);
	end_setRootInstance(block_get5End(block), segment_get5Cap(segment));
	end_setRootInstance(block_get3End(block), segment_get3Cap(segment));
}

Block_InstanceIterator *block_getInstanceIterator(Block *block) {
	Block_InstanceIterator *iterator;
	iterator = st_malloc(sizeof(struct _block_instanceIterator));
	iterator->block = block;
	iterator->iterator = stSortedSet_getIterator(block->blockContents->segments);
	return iterator;
}

Segment *block_getNext(Block_InstanceIterator *iterator) {
	return block_getInstanceP(iterator->block, stSortedSet_getNext(iterator->iterator));
}

Segment *block_getPrevious(Block_InstanceIterator *iterator) {
	return block_getInstanceP(iterator->block, stSortedSet_getPrevious(iterator->iterator));
}

Block_InstanceIterator *block_copyInstanceIterator(Block_InstanceIterator *iterator) {
	Block_InstanceIterator *iterator2;
	iterator2 = st_malloc(sizeof(struct _block_instanceIterator));
	iterator2->block = iterator->block;
	iterator2->iterator = stSortedSet_copyIterator(iterator->iterator);
	return iterator2;
}

void block_destructInstanceIterator(Block_InstanceIterator *iterator) {
	stSortedSet_destructIterator(iterator->iterator);
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

Segment *block_splitP(Segment *segment,
		Block *leftBlock, Block *rightBlock) {
	Segment *leftSegment = segment_getSequence(segment) != NULL
				? segment_construct2(leftBlock, segment_getStart(segment), segment_getStrand(segment), segment_getSequence(segment))
				: segment_construct(leftBlock, segment_getEvent(segment));
				Segment *rightSegment = segment_getSequence(segment) != NULL
				? segment_construct2(rightBlock, segment_getStart(segment) + block_getLength(leftBlock), segment_getStrand(segment), segment_getSequence(segment))
						: segment_construct(rightBlock, segment_getEvent(segment));
	//link together.
	cap_makeAdjacent(segment_get3Cap(leftSegment), segment_get5Cap(rightSegment));
	//update adjacencies.
	Cap *_5Cap = segment_get5Cap(segment);
	Cap *new5Cap = segment_get5Cap(leftSegment);
	Cap *_3Cap = segment_get3Cap(segment);
	Cap *new3Cap = segment_get3Cap(rightSegment);
	if(cap_getAdjacency(_5Cap) != NULL) {
		cap_makeAdjacent(cap_getAdjacency(_5Cap), new5Cap);
	}
	if(cap_getAdjacency(_3Cap) != NULL) {
		cap_makeAdjacent(cap_getAdjacency(_3Cap), new3Cap);
	}
	return leftSegment;
}

static void block_splitP2(Segment *segment,
		Segment *parentLeftSegment,
		Segment *parentRightSegment,
		Block *leftBlock, Block *rightBlock) {
	Segment *leftSegment = block_splitP(segment, leftBlock, rightBlock);
	Segment *rightSegment = cap_getSegment(cap_getAdjacency(segment_get3Cap(leftSegment)));
	if(parentLeftSegment != NULL) {
		assert(parentRightSegment != NULL);
		segment_makeParentAndChild(parentLeftSegment, leftSegment);
		segment_makeParentAndChild(parentRightSegment, rightSegment);
	}
	else {
		assert(parentRightSegment == NULL);
		block_setRootInstance(leftBlock, leftSegment);
		block_setRootInstance(rightBlock, rightSegment);
	}
	int32_t i;
	for(i=0; i<segment_getChildNumber(segment); i++) {
		block_splitP2(segment_getChild(segment, i), leftSegment, rightSegment, leftBlock, rightBlock);
	}
}

void block_split(Block *block, int32_t splitPoint, Block **leftBlock, Block **rightBlock) {
	assert(splitPoint > 0);
	assert(splitPoint < block_getLength(block));
	*leftBlock = block_construct(splitPoint, block_getFlower(block));
	*rightBlock = block_construct(block_getLength(block) - splitPoint, block_getFlower(block));

	Segment *segment = block_getRootInstance(block);
	if(segment != NULL) {
		block_splitP2(segment, NULL, NULL, *leftBlock, *rightBlock);
	}
	else {
		Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
		while((segment = block_getNext(instanceIterator)) != NULL) {
			block_splitP(segment, *leftBlock, *rightBlock);
		}
		block_destructInstanceIterator(instanceIterator);
	}
	block_destruct(block);
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

char *block_makeNewickStringP(Segment *segment, int32_t includeInternalNames, int32_t includeUnaryEvents) {
	if(!includeUnaryEvents && segment_getChildNumber(segment) == 1) {
		return block_makeNewickStringP(segment_getChild(segment, 0), includeInternalNames, includeUnaryEvents);
	}
	if(segment_getChildNumber(segment) > 0) {
		char *left = stString_print("(");
		int32_t i;
		int32_t comma = 0;
		for(i=0; i<segment_getChildNumber(segment); i++) {
			Segment *childSegment = segment_getChild(segment, i);
			char *cA = block_makeNewickStringP(childSegment, includeInternalNames, includeUnaryEvents);
			char *cA2 = stString_print(comma ? "%s,%s" : "%s%s", left, cA);
			free(cA);
			left = cA2;
			comma = 1;
		}
		char *final = includeInternalNames ?
				stString_print("%s)%s", left, cactusMisc_nameToStringStatic(segment_getName(segment))) :
				stString_print("%s)", left);
		free(left);
		return final;
	}
	return cactusMisc_nameToString(segment_getName(segment));
}

char *block_makeNewickString(Block *block, int32_t includeInternalNames, int32_t includeUnaryEvents) {
	Segment *segment = block_getRootInstance(block);
	if(segment != NULL) {
		assert(segment != NULL);
		char *cA = block_makeNewickStringP(segment, includeInternalNames, includeUnaryEvents);
		char *cA2 = stString_print("%s;", cA);
		free(cA);
		return cA2;
	}
	return NULL;
}

/*
 * Private functions.
 */

void block_addInstance(Block *block, Segment *segment) {
	stSortedSet_insert(block->blockContents->segments, segment_getPositiveOrientation(segment));
}

void block_removeInstance(Block *block, Segment *segment) {
	stSortedSet_remove(block->blockContents->segments, segment_getPositiveOrientation(segment));
}

void block_setFlower(Block *block, Flower *flower) {
	flower_removeBlock(block_getFlower(block), block);
	block->blockContents->flower = flower;
	flower_addBlock(flower, block);
}

/*
 * Serialisation functions.
 */

void block_writeBinaryRepresentation(Block *block, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Block_InstanceIterator *iterator;
	Segment *segment;

	assert(block_getOrientation(block));
	binaryRepresentation_writeElementType(CODE_BLOCK, writeFn);
	binaryRepresentation_writeName(block_getName(block), writeFn);
	binaryRepresentation_writeInteger(block_getLength(block), writeFn);
	binaryRepresentation_writeName(end_getName(block_get5End(block)), writeFn);
	binaryRepresentation_writeName(end_getName(block_get3End(block)), writeFn);
	iterator = block_getInstanceIterator(block);
	while((segment = block_getNext(iterator)) != NULL) {
		segment_writeBinaryRepresentation(segment, writeFn);
	}
	block_destructInstanceIterator(iterator);
}

Block *block_loadFromBinaryRepresentation(void **binaryString, Flower *flower) {
	Block *block;
	Name name, leftEndName, rightEndName;
	int32_t length;

	block = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_BLOCK) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		length = binaryRepresentation_getInteger(binaryString);
		leftEndName = binaryRepresentation_getName(binaryString);
		rightEndName = binaryRepresentation_getName(binaryString);
		block = block_construct2(name, length, flower_getEnd(flower, leftEndName), flower_getEnd(flower, rightEndName), flower);
		while(segment_loadFromBinaryRepresentation(binaryString, block) != NULL);
	}
	return block;
}

Block *block_getStaticNameWrapper(Name name) {
	static Block block;
	static BlockContents blockContents;
	block.blockContents = &blockContents;
	blockContents.name = name;
	return &block;
}
