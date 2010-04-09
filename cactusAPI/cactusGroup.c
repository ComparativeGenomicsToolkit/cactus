#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic group functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int32_t group_constructP(const void *o1, const void *o2, void *a) {
	assert(a == NULL);
	return netMisc_nameCompare(end_getName((End *)o1), end_getName((End *)o2));
}

Group *group_construct(Net *net, Net *nestedNet) {
	Group *group;

	group = group_construct3(net, net_getName(nestedNet), 0);
	group_updateContainedEnds(group);
	net_setParentGroup(nestedNet, group);
	return group;
}

Group *group_construct2(Net *net) {
	Group *group;

	group = group_construct3(net, netDisk_getUniqueID(net_getNetDisk(net)), 1);
	return group;
}

bool group_isTerminal(Group *group) {
	return group->terminalGroup;
}

static int32_t returnsTrue(Event *event) {
	assert(event != NULL);
	return 1;
}

void group_makeNonTerminal(Group *group) {
	assert(group_isTerminal(group));
	group->terminalGroup = 0;
	Net *nestedNet = net_construct2(group_getName(group), net_getNetDisk(group_getNet(group)));
	net_setParentGroup(nestedNet, group);
	if(net_getEventTree(group_getNet(group)) != NULL) {
		eventTree_copyConstruct(net_getEventTree(group_getNet(group)), nestedNet, returnsTrue);
	}
	Group *nestedGroup = group_construct2(nestedNet);
	//Add the ends to the nested net.
	Group_EndIterator *endIterator = group_getEndIterator(group);
	End *end;
	while((end = group_getNextEnd(endIterator)) != NULL) {
		end_setGroup(end_copyConstruct(end, nestedNet), nestedGroup);
	}
	group_destructEndIterator(endIterator);
	netMisc_addAdjacenciesToLeafCaps(nestedNet);
	assert(group_getTotalBaseLength(group) == net_getTotalBaseLength(nestedNet));
}

void group_updateContainedEnds(Group *group) {
	assert(!group_isTerminal(group));
	Net *net;
	Net_EndIterator *iterator;
	End *end;
	End *end2;
	//wipe the slate clean.
	while(group_getEndNumber(group) != 0) {
		end_setGroup(group_getFirstEnd(group), NULL);
	}
	sortedSet_destruct(group->ends, NULL);
	group->ends = sortedSet_construct(group_constructP);
	//now calculate the ends
	net = group_getNet(group);
	iterator = net_getEndIterator(group_getNestedNet(group));
	while((end = net_getNextEnd(iterator)) != NULL) {
		if((end2 = net_getEnd(net, end_getName(end))) != NULL) {
			end_setGroup(end2, group);
		}
	}
	net_destructEndIterator(iterator);
}

void group_addEnd(Group *group, End *end) {
	end = end_getPositiveOrientation(end);
	sortedSet_insert(group->ends, end);
}

void group_destruct(Group *group) {
	//Detach from the parent net.
	net_removeGroup(group_getNet(group), group);
	while(group_getEndNumber(group) != 0) {
		end_setGroup(group_getFirstEnd(group), NULL);
	}
	sortedSet_destruct(group->ends, NULL);
	//Free the memory
	free(group);
}

Net *group_getNet(Group *group) {
	return group->net;
}

Name group_getName(Group *group) {
	return group->name;
}

Net *group_getNestedNet(Group *group) {
	return group_isTerminal(group) ? NULL : netDisk_getNet(net_getNetDisk(group_getNet(group)), group->name);
}

Link *group_getLink(Group *group) {
	return group->link;
}

End *group_getFirstEnd(Group *group) {
	return sortedSet_getFirst(group->ends);
}

End *group_getEnd(Group *group, Name name) {
	static End end;
	static EndContents endContents;
	end.endContents = &endContents;
	endContents.name = name;
	return sortedSet_find(group->ends, &end);
}

int32_t group_getEndNumber(Group *group) {
	return sortedSet_getLength(group->ends);
}

Group_EndIterator *group_getEndIterator(Group *group) {
	return iterator_construct(group->ends);
}

End *group_getNextEnd(Group_EndIterator *endIterator) {
	return iterator_getNext(endIterator);
}

End *group_getPreviousEnd(Group_EndIterator *endIterator) {
	return iterator_getPrevious(endIterator);
}

Group_EndIterator *group_copyEndIterator(Group_EndIterator *endIterator) {
	return iterator_copy(endIterator);
}

void group_destructEndIterator(Group_EndIterator *endIterator) {
	iterator_destruct(endIterator);
}

int64_t group_getTotalBaseLength(Group *group) {
	Group_EndIterator *endIterator = group_getEndIterator(group);
	End *end;
	int64_t totalLength = 0;
	while((end = group_getNextEnd(endIterator)) != NULL) {
		End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
		Cap *cap;
		while((cap = end_getNext(instanceIterator)) != NULL) {
			cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
			if(!cap_getSide(cap)) {
				Cap *cap2 = cap_getAdjacency(cap);
				assert(cap_getStrand(cap2));
				assert(cap_getSide(cap2));
				assert(end_getGroup(cap_getEnd(cap2)) == group);
				int32_t length = cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1;
				assert(length >= 0);
				totalLength += length;
			}
		}
		end_destructInstanceIterator(instanceIterator);
	}
	group_destructEndIterator(endIterator);
	return totalLength;
}

Group *group_mergeGroups(Group *group1, Group *group2) {
	//Check they are in the same net..
	assert(group_getNet(group1) == group_getNet(group2));
	assert(group1 != group2);

	if(!group_isTerminal(group1) || !group_isTerminal(group2)) { //We must first merge the nested nets
		if(group_isTerminal(group1)) { //Need to make a nested net to merge with the other
			group_makeNonTerminal(group1);
		}
		if(group_isTerminal(group2)) { //Need to make a nested net to merge with the other
			group_makeNonTerminal(group2);
		}
		assert(group_getNestedNet(group1) != NULL);
		assert(group_getNestedNet(group2) != NULL);
		net_mergeNetsP(group_getNestedNet(group1), group_getNestedNet(group2));
	}
	End *end;
	while((end = group_getFirstEnd(group1)) != NULL) {
		end_setGroup(end, group2);
	}
	group_destruct(group1);

#ifdef BEN_DEBUG
	if(!group_isTerminal(group2)) {
		assert(net_getParentGroup(group_getNestedNet(group2)) == group2);
	}
#endif

	return group2;
}

/*
 * Private functions.
 */

Group *group_construct3(Net *net, Name name, bool terminalGroup) {
	Group *group;
	group = malloc(sizeof(Group));

	group->net = net;
	group->link = NULL;
	group->name = name;
	group->ends = sortedSet_construct(group_constructP);
	group->terminalGroup = terminalGroup;
	net_addGroup(net, group);

	return group;
}

void group_setLink(Group *group, Link *link) {
	//argument may be NULL
	group->link = link;
	assert(group_getEnd(group, end_getName(link_getLeft(link))) == link_getLeft(link));
	assert(group_getEnd(group, end_getName(link_getRight(link))) == link_getRight(link));
}

void group_removeEnd(Group *group, End *end) {
	assert(group_getEnd(group, end_getName(end)) == end);
	sortedSet_delete(group->ends, end);
}

void group_setNet(Group *group, Net *net) {
	net_removeGroup(group_getNet(group), group);
	group->net = net;
	net_addGroup(net, group);
	Net *nestedNet = group_getNestedNet(group);
	if(nestedNet != NULL) { //we re-do this link, because the parent net has changed.
		net_setParentGroup(nestedNet, group);
	}
}

/*
 * Serialisation functions
 */

void group_writeBinaryRepresentation(Group *group, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	End *end;
	Group_EndIterator *iterator;

	binaryRepresentation_writeElementType(CODE_ADJACENCY_COMPONENT, writeFn);
	binaryRepresentation_writeBool(group_isTerminal(group), writeFn);
	binaryRepresentation_writeName(group_getName(group), writeFn);
	iterator = group_getEndIterator(group);
	while((end = group_getNextEnd(iterator)) != NULL) {
		binaryRepresentation_writeElementType(CODE_ADJACENCY_COMPONENT_END, writeFn);
		binaryRepresentation_writeName(end_getName(end), writeFn);
	}
	group_destructEndIterator(iterator);
}

Group *group_loadFromBinaryRepresentation(void **binaryString, Net *net) {
	Group *group;

	group = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ADJACENCY_COMPONENT) {
		binaryRepresentation_popNextElementType(binaryString);
		bool terminalGroup = binaryRepresentation_getBool(binaryString);
		Name name = binaryRepresentation_getName(binaryString);
		group = group_construct3(net, name, terminalGroup);
		while(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ADJACENCY_COMPONENT_END) {
			binaryRepresentation_popNextElementType(binaryString);
			end_setGroup(net_getEnd(net, binaryRepresentation_getName(binaryString)), group);
		}
	}
	return group;
}

Group *group_getStaticNameWrapper(Name netName) {
	static Group group;
	group.name = netName;
	return &group;
}

