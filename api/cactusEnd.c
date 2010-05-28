#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int32_t end_constructP(const void *o1, const void *o2) {
	return netMisc_nameCompare(cap_getName((Cap *)o1), cap_getName((Cap *)o2));
}

End *end_construct(bool isAttached, Net *net) {
	return end_construct2(netDisk_getUniqueID(net_getNetDisk(net)), 1, isAttached, 1, net);
}

End *end_construct2(Name name, int32_t isStub, int32_t isAttached, int32_t side, Net *net) {
	End *end;
	end = st_malloc(sizeof(End));
	end->rEnd = st_malloc(sizeof(End));
	end->rEnd->rEnd = end;
	end->endContents = st_malloc(sizeof(EndContents));
	end->rEnd->endContents = end->endContents;

	end->orientation = 1;
	end->rEnd->orientation = 0;

	end->side = side;
	end->rEnd->side = !side;

	if(!isStub) {
		assert(!isAttached);
	}

	end->endContents->isStub = isStub;
	end->endContents->isAttached = isAttached;

	end->endContents->rootInstance = NULL;
	end->endContents->name = name;
	end->endContents->caps = stSortedSet_construct3(end_constructP, NULL);
	end->endContents->attachedBlock = NULL;
	end->endContents->group = NULL;
	end->endContents->net = net;
	net_addEnd(net, end);
	return end;
}

End *end_copyConstruct(End *end, Net *newNet) {
	End *end2;
	End_InstanceIterator *iterator;
	Cap *cap;
	Cap *cap2;

	assert(net_getEnd(newNet, end_getName(end)) == NULL);

	end2 = end_construct2(end_getName(end), 1, end_isBlockEnd(end) ? 1 : end_isAttached(end), end_getSide(end), newNet);
	//Copy the instances.
	iterator = end_getInstanceIterator(end);
	while((cap = end_getNext(iterator)) != NULL) {
		cap_copyConstruct(end2, cap);
	}
	end_destructInstanceIterator(iterator);

	//Copy any parent child links.
	iterator = end_getInstanceIterator(end);
	while((cap = end_getNext(iterator)) != NULL) {
		if((cap2 = cap_getParent(cap)) != NULL) {
			cap_makeParentAndChild(end_getInstance(end2, cap_getName(cap2)),
										   end_getInstance(end2, cap_getName(cap)));
		}
	}
	end_destructInstanceIterator(iterator);

	//Copy root.
	if(end_getRootInstance(end) != NULL) {
		end_setRootInstance(end2, end_getInstance(end, cap_getName(end_getRootInstance(end))));
	}
	return end2;
}

void end_destruct(End *end) {
	Cap *cap;
	//remove from net.
	net_removeEnd(end_getNet(end), end);

	//remove from group.
	end_setGroup(end, NULL);

	//remove instances
	while((cap = end_getFirst(end)) != NULL) {
		cap_destruct(cap);
	}
	//now the actual instances.
	stSortedSet_destruct(end->endContents->caps);

	free(end->endContents);
	free(end->rEnd);
	free(end);
}

void end_setBlock(End *end, Block *block) {
	assert(end_getOrientation(end));
	assert(block_getOrientation(block));
	end->endContents->attachedBlock = block;
}

Name end_getName(End *end) {
	return end->endContents->name;
}

bool end_getOrientation(End *end) {
	return end->orientation;
}

End *end_getPositiveOrientation(End *end) {
	return end_getOrientation(end) ? end : end_getReverse(end);
}

End *end_getReverse(End *end) {
	return end->rEnd;
}

bool end_getSide(End *end) {
	return end->side;
}

Net *end_getNet(End *end) {
	return end->endContents->net;
}

Block *end_getBlock(End *end) {
	Block *a = end->endContents->attachedBlock;
	return a == NULL || end_getOrientation(end) ? a : block_getReverse(a);
}

End *end_getOtherBlockEnd(End *end) {
	if(!end_isBlockEnd(end)) {
		return NULL; //the end must be block end to return the other end of a block!
	}
	Block *block = end_getBlock(end);
	assert(block != NULL);
	End *otherEnd = end_getSide(end) ? block_get3End(block) :  block_get5End(block);
	assert(end_getOrientation(end) == end_getOrientation(otherEnd));
	assert(end != otherEnd);
	return otherEnd;
}

Group *end_getGroup(End *end) {
	return end->endContents->group;
}

int32_t end_getInstanceNumber(End *end) {
	return stSortedSet_size(end->endContents->caps);
}

Cap *end_getInstanceP(End *end, Cap *connectedCap) {
	return connectedCap == NULL ? connectedCap : (end_getOrientation(end) ? connectedCap : cap_getReverse(connectedCap));
}

Cap *end_getInstance(End *end, Name name) {
	Cap *cap = cap_getStaticNameWrapper(name);
	return end_getInstanceP(end, stSortedSet_search(end->endContents->caps, cap));
}

Cap *end_getFirst(End *end) {
	return end_getInstanceP(end, stSortedSet_getFirst(end->endContents->caps));
}

Cap *end_getRootInstance(End *end) {
	return end_getInstanceP(end, end->endContents->rootInstance);
}

void end_setRootInstance(End *end, Cap *cap) {
	end->endContents->rootInstance = cap_getOrientation(cap) ? cap : cap_getReverse(cap);
}

End_InstanceIterator *end_getInstanceIterator(End *end) {
	End_InstanceIterator *iterator;
	iterator = st_malloc(sizeof(struct _end_instanceIterator));
	iterator->end = end;
	iterator->iterator = stSortedSet_getIterator(end->endContents->caps);
	return iterator;
}

Cap *end_getNext(End_InstanceIterator *iterator) {
	return end_getInstanceP(iterator->end, stSortedSet_getNext(iterator->iterator));
}

Cap *end_getPrevious(End_InstanceIterator *iterator) {
	return end_getInstanceP(iterator->end, stSortedSet_getPrevious(iterator->iterator));
}

End_InstanceIterator *end_copyInstanceIterator(End_InstanceIterator *iterator) {
	End_InstanceIterator *iterator2;
	iterator2 = st_malloc(sizeof(struct _end_instanceIterator));
	iterator2->end = iterator->end;
	iterator2->iterator = stSortedSet_copyIterator(iterator->iterator);
	return iterator2;
}

void end_destructInstanceIterator(End_InstanceIterator *iterator) {
	stSortedSet_destructIterator(iterator->iterator);
	free(iterator);
}

bool end_isBlockEnd(End *end) {
	return !end_isStubEnd(end);
}

bool end_isStubEnd(End *end) {
	return end->endContents->isStub;
}

bool end_isAttached(End *end) {
	return end->endContents->isAttached;
}

bool end_isFree(End *end) {
	return !end_isAttached(end);
}

void end_setGroup(End *end, Group *group) {
	if(end_getGroup(end) != NULL) {
		group_removeEnd(end_getGroup(end), end);
	}
	end->endContents->group = group;
	if(group != NULL) {
		group_addEnd(group, end);
	}
}

void end_check(End *end) {
	//Check is connected to net properly
	assert(net_getEnd(end_getNet(end), end_getName(end)) == end_getPositiveOrientation(end));

	//check end is part of group..
	Group *group = end_getGroup(end);
	assert(group != NULL);
	assert(group_getEnd(group, end_getName(end)) == end_getPositiveOrientation(end));

	if(end_isBlockEnd(end)) {
		assert(!end_isStubEnd(end));
		assert(end_isFree(end));
		//Check block..
		Block *block = end_getBlock(end);
		assert(block != NULL);
		assert(block_getOrientation(block) == end_getOrientation(end));
		//check not attached
		assert(end_isFree(end));
		assert(!end_isAttached(end));
		//Check sides correspond..
		if(end_getSide(end)) {
			assert(block_get5End(block) == end);
		}
		else {
			assert(block_get3End(block) == end);
		}
	}
	else {
		assert(end_isStubEnd(end)); //Is stub end:
		//there must be no attached block.
		assert(end_getBlock(end) == NULL);
		Group *parentGroup = net_getParentGroup(end_getNet(end));
		if(parentGroup != NULL) {
			// if attached the is inherited from a parent net to the containing net.
			End *parentEnd = group_getEnd(parentGroup, end_getName(end));
			if(end_isAttached(end)) {
				assert(parentEnd != NULL);
			}
			if(parentEnd != NULL) {
				assert(end_getSide(end_getPositiveOrientation(parentEnd)) == end_getSide(end_getPositiveOrientation(end)));
			}
		}
	}

	//Check reverse, not comprehensively, perhaps.
	End *rEnd = end_getReverse(end);
	assert(rEnd != NULL);
	assert(end_getReverse(rEnd) == end);
	assert(end_getOrientation(end) == !end_getOrientation(rEnd));
	assert(end_getSide(end) == !end_getSide(rEnd));
	assert(end_getName(end) == end_getName(rEnd));
	assert(end_getInstanceNumber(end) == end_getInstanceNumber(rEnd));
	assert(end_isAttached(end) == end_isAttached(rEnd));
	assert(end_isStubEnd(end) == end_isStubEnd(rEnd));
	if(end_getRootInstance(end) == NULL) {
		assert(end_getRootInstance(rEnd) == NULL);
	}
	else {
		assert(end_getRootInstance(end) == cap_getReverse(end_getRootInstance(rEnd)));
	}
	if(end_getInstanceNumber(end) > 0) {
		assert(end_getFirst(end) == cap_getReverse(end_getFirst(rEnd)));
	}

	//Check has tree if built_trees set
	if(net_builtTrees(end_getNet(end)) && end_getInstanceNumber(end) > 0) {
		assert(end_getRootInstance(end) != NULL);
	}

	//For each segment calls segment_check.
	End_InstanceIterator *iterator = end_getInstanceIterator(end);
	Cap *cap;
	while((cap = end_getNext(iterator)) != NULL) {
		cap_check(cap);
	}
	end_destructInstanceIterator(iterator);
}


/*
 * Private functions
 */

void end_addInstance(End *end, Cap *cap) {
	stSortedSet_insert(end->endContents->caps, cap_getPositiveOrientation(cap));
}

void end_removeInstance(End *end, Cap *cap) {
	stSortedSet_remove(end->endContents->caps, cap);
}

void end_setNet(End *end, Net *net) {
	net_removeEnd(end_getNet(end), end);
	end->endContents->net = net;
	net_addEnd(net, end);
}

/*
 * Serialisation functions.
 */


void end_writeBinaryRepresentationP(Cap *cap, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int32_t i;
	cap_writeBinaryRepresentation(cap, writeFn);
	for(i=0; i<cap_getChildNumber(cap); i++) {
		end_writeBinaryRepresentationP(cap_getChild(cap, i), writeFn);
	}
}

void end_writeBinaryRepresentation(End *end, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	End_InstanceIterator *iterator;
	Cap *cap;

	assert(end_getOrientation(end));
	cap = end_getRootInstance(end);
	binaryRepresentation_writeElementType(cap == NULL ?
			CODE_END_WITHOUT_PHYLOGENY : CODE_END_WITH_PHYLOGENY, writeFn);
	binaryRepresentation_writeName(end_getName(end), writeFn);
	binaryRepresentation_writeBool(end_isStubEnd(end), writeFn);
	binaryRepresentation_writeBool(end_isAttached(end), writeFn);
	binaryRepresentation_writeBool(end_getSide(end), writeFn);

	if(cap == NULL) {
		iterator = end_getInstanceIterator(end);
		while((cap = end_getNext(iterator)) != NULL) {
			assert(cap_getParent(cap) == NULL);
			cap_writeBinaryRepresentation(cap, writeFn);
		}
		end_destructInstanceIterator(iterator);
	}
	else {
		end_writeBinaryRepresentationP(cap, writeFn);
	}
}

End *end_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	End *end;
	Name name;
	int32_t isStub;
	int32_t isAttached;
	int32_t side;

	end = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_END_WITHOUT_PHYLOGENY) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		isStub = binaryRepresentation_getBool(binaryString);
		isAttached = binaryRepresentation_getBool(binaryString);
		side = binaryRepresentation_getBool(binaryString);
		end = end_construct2(name, isStub, isAttached, side, net);
		while(cap_loadFromBinaryRepresentation(binaryString, end) != NULL);
	}
	else {
		if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_END_WITH_PHYLOGENY) {
			binaryRepresentation_popNextElementType(binaryString);
			name = binaryRepresentation_getName(binaryString);
			isStub = binaryRepresentation_getBool(binaryString);
			isAttached = binaryRepresentation_getBool(binaryString);
			side = binaryRepresentation_getBool(binaryString);
			end = end_construct2(name, isStub, isAttached, side, net);
			end_setRootInstance(end, cap_loadFromBinaryRepresentation(binaryString, end));
			while(cap_loadFromBinaryRepresentation(binaryString, end) != NULL);
		}
	}

	return end;
}

End *end_getStaticNameWrapper(Name name) {
	static End end;
	static EndContents endContents;
	end.endContents = &endContents;
	endContents.name = name;
	end.orientation = 1;
	return &end;
}

uint32_t end_hashKey(const void *o) {
	return end_getName((End *)o);
}

int end_hashEqualsKey(const void *o, const void *o2) {
	End *end1 = (End *)o;
	End *end2 = (End *)o2;
	return end_getName(end1) == end_getName(end2) && end_getOrientation(end1) == end_getOrientation(end2);
}
