#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic end functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t end_constructP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(cap_getName((Cap *)o1), cap_getName((Cap *)o2));
}

End *end_construct(bool isAttached, Net *net) {
	return end_construct2(netDisk_getUniqueID(net_getNetDisk(net)), 1, isAttached, net);
}

End *end_construct2(Name name, int32_t isStub, int32_t isAttached, Net *net) {
	End *end;
	end = malloc(sizeof(End));
	end->rEnd = malloc(sizeof(End));
	end->rEnd->rEnd = end;
	end->endContents = malloc(sizeof(EndContents));
	end->rEnd->endContents = end->endContents;

	end->orientation = 1;
	end->rEnd->orientation = 0;

	if(!isStub) {
		assert(!isAttached);
	}

	end->endContents->isStub = isStub;
	end->endContents->isAttached = isAttached;

	end->endContents->rootInstance = NULL;
	end->endContents->name = name;
	end->endContents->caps = sortedSet_construct(end_constructP);
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

	end2 = end_construct2(end_getName(end), 1, end_isBlockEnd(end) ? 1 : end_isAttached(end), newNet);
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

	//remove instances
	while((cap = end_getFirst(end)) != NULL) {
		cap_destruct(cap);
	}
	//now the actual instances.
	sortedSet_destruct(end->endContents->caps, NULL);

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

Net *end_getNet(End *end) {
	return end->endContents->net;
}

Block *end_getBlock(End *end) {
	Block *a = end->endContents->attachedBlock;
	return a == NULL || end_getOrientation(end) ? a : block_getReverse(a);
}

Group *end_getGroup(End *end) {
	return end->endContents->group;
}

int32_t end_getInstanceNumber(End *end) {
	return sortedSet_getLength(end->endContents->caps);
}

Cap *end_getInstanceP(End *end, Cap *connectedCap) {
	return connectedCap == NULL ? connectedCap : (end_getOrientation(end) ? connectedCap : cap_getReverse(connectedCap));
}

Cap *end_getInstance(End *end, Name name) {
	Cap *cap = cap_getStaticNameWrapper(name);
	return end_getInstanceP(end, sortedSet_find(end->endContents->caps, cap));
}

Cap *end_getFirst(End *end) {
	return end_getInstanceP(end, sortedSet_getFirst(end->endContents->caps));
}

Cap *end_getRootInstance(End *end) {
	return end_getInstanceP(end, end->endContents->rootInstance);
}

void end_setRootInstance(End *end, Cap *cap) {
	end->endContents->rootInstance = cap_getOrientation(cap) ? cap : cap_getReverse(cap);
}

End_InstanceIterator *end_getInstanceIterator(End *end) {
	End_InstanceIterator *iterator;
	iterator = malloc(sizeof(struct _end_instanceIterator));
	iterator->end = end;
	iterator->iterator = iterator_construct(end->endContents->caps);
	return iterator;
}

Cap *end_getNext(End_InstanceIterator *iterator) {
	return end_getInstanceP(iterator->end, iterator_getNext(iterator->iterator));
}

Cap *end_getPrevious(End_InstanceIterator *iterator) {
	return end_getInstanceP(iterator->end, iterator_getPrevious(iterator->iterator));
}

End_InstanceIterator *end_copyInstanceIterator(End_InstanceIterator *iterator) {
	End_InstanceIterator *iterator2;
	iterator2 = malloc(sizeof(struct _end_instanceIterator));
	iterator2->end = iterator->end;
	iterator2->iterator = iterator_copy(iterator->iterator);
	return iterator2;
}

void end_destructInstanceIterator(End_InstanceIterator *iterator) {
	iterator_destruct(iterator->iterator);
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


/*
 * Private functions
 */

void end_addInstance(End *end, Cap *cap) {
	sortedSet_insert(end->endContents->caps, cap_getPositiveOrientation(cap));
}

void end_removeInstance(End *end, Cap *cap) {
	sortedSet_delete(end->endContents->caps, cap);
}

void end_setGroup(End *end, Group *group) {
	//argument may be NULL
	end->endContents->group = group;
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

	end = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_END_WITHOUT_PHYLOGENY) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		isStub = binaryRepresentation_getBool(binaryString);
		isAttached = binaryRepresentation_getBool(binaryString);
		end = end_construct2(name, isStub, isAttached, net);
		while(cap_loadFromBinaryRepresentation(binaryString, end) != NULL);
	}
	else {
		if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_END_WITH_PHYLOGENY) {
			binaryRepresentation_popNextElementType(binaryString);
			name = binaryRepresentation_getName(binaryString);
			isStub = binaryRepresentation_getBool(binaryString);
			isAttached = binaryRepresentation_getBool(binaryString);
			end = end_construct2(name, isStub, isAttached, net);
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
	return &end;
}
