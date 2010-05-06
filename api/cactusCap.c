#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic cap functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Cap *cap_construct(End *end, Event *event) {
	return cap_construct3(netDisk_getUniqueID(net_getNetDisk(end_getNet(end))), event, end, 1);
}

Cap *cap_construct3(Name instance, Event *event, End *end, int32_t side) {
	assert(end != NULL);
	assert(event != NULL);
	assert(instance != NULL_NAME);
	Cap *cap;

	cap = malloc(sizeof(Cap));
	cap->capContents = malloc(sizeof(CapContents));
	cap->rCap = malloc(sizeof(Cap));
	cap->rCap->rCap = cap;
	cap->rCap->capContents = cap->capContents;

	cap->end = end;
	cap->rCap->end = end_getReverse(end);

	cap->capContents->instance = instance;
	cap->capContents->coordinate = INT32_MAX;
	cap->capContents->sequence = NULL;
	cap->capContents->adjacency = NULL;
	cap->capContents->adjacency2 = NULL;
	cap->capContents->face = NULL;
	cap->capContents->segment = NULL;
	cap->capContents->parent = NULL;
	cap->capContents->children = constructEmptyList(0, NULL);
	cap->capContents->event = event;
	cap->capContents->strand = end_getOrientation(end);
	cap->capContents->side = end_getOrientation(end) ? side : !side;

	end_addInstance(end, cap);
	net_addCap(end_getNet(end), cap);
	return cap;
}

Cap *cap_construct2(End *end,
		int32_t coordinate, bool strand, bool side, Sequence *sequence) {
	return cap_construct4(netDisk_getUniqueID(net_getNetDisk(end_getNet(end))),
			end, coordinate, strand, side, sequence);
}

Cap *cap_construct4(Name instance, End *end,
		int32_t coordinate, int32_t strand, int32_t side, Sequence *sequence) {
	Cap *cap;
	cap = cap_construct3(instance, sequence_getEvent(sequence), end, side);
	cap->capContents->coordinate = coordinate;
	cap->capContents->strand = cap_getOrientation(cap) ? strand : !strand;
	cap->capContents->sequence = sequence;
	return cap;
}

Cap *cap_construct5(Event *event, End *end, int32_t side) {
	return cap_construct3(netDisk_getUniqueID(net_getNetDisk(end_getNet(end))), event, end, side);
}

Cap *cap_copyConstruct(End *end, Cap *cap) {
	assert(end_getName(cap_getEnd(cap)) == end_getName(end));
	Event *event;
	Name sequenceName;
	Sequence *sequence;

	Net *net = end_getNet(end);
	if(cap_getCoordinate(cap) != INT32_MAX) {
		sequenceName = sequence_getName(cap_getSequence(cap));
		sequence = net_getSequence(net, sequenceName);
		if(sequence == NULL) { //add sequence to the net.
			sequence = sequence_construct(netDisk_getMetaSequence(net_getNetDisk(net), sequenceName), net);
			assert(sequence != NULL);
		}
		return cap_construct4(cap_getName(cap), end,
				cap_getCoordinate(cap), cap_getStrand(cap),
				cap_getSide(cap), sequence);
	}
	else {
		event = eventTree_getEvent(net_getEventTree(net), event_getName(cap_getEvent(cap)));
		assert(event != NULL);
		return cap_construct3(cap_getName(cap), event, end, cap_getSide(cap));
	}
}

void cap_destruct(Cap *cap) {
	//Remove from end.
	end_removeInstance(cap_getEnd(cap), cap);
	net_removeCap(end_getNet(cap_getEnd(cap)), cap);

	destructList(cap->capContents->children);
	free(cap->rCap);
	free(cap->capContents);
	free(cap);
}

Name cap_getName(Cap *cap) {
	return cap->capContents->instance;
}

bool cap_getOrientation(Cap *cap) {
	return end_getOrientation(cap_getEnd(cap));
}

Cap *cap_getPositiveOrientation(Cap *cap) {
	return cap_getOrientation(cap) ? cap : cap_getReverse(cap);
}


Cap *cap_getReverse(Cap *cap) {
	return cap->rCap;
}

Event *cap_getEvent(Cap *cap) {
	return cap->capContents->event;
}

End *cap_getEnd(Cap *cap) {
	return cap->end;
}

Segment *cap_getSegment(Cap *cap) {
	return cap_getOrientation(cap) ?
			cap->capContents->segment :
	(cap->capContents->segment != NULL ? segment_getReverse(cap->capContents->segment) : NULL);
}

Cap *cap_getOtherSegmentCap(Cap *cap) {
	if(!end_isBlockEnd(cap_getEnd(cap))) {
		assert(cap_getSegment(cap) == NULL);
		return NULL;
	}
	Segment *segment = cap_getSegment(cap);
	assert(segment != NULL);
	Cap *otherCap = cap_getSide(cap) ? segment_get3Cap(segment) : segment_get5Cap(segment);
	assert(cap != otherCap);
	return otherCap;
}

int32_t cap_getCoordinate(Cap *cap) {
	return cap->capContents->coordinate;
}

bool cap_getStrand(Cap *cap) {
	return cap_getOrientation(cap) ? cap->capContents->strand : !cap->capContents->strand;
}

bool cap_getSide(Cap *cap) {
	return cap_getOrientation(cap) ? cap->capContents->side : !cap->capContents->side;
}

Sequence *cap_getSequence(Cap *cap) {
	return cap->capContents->sequence;
}

#ifdef BEN_DEBUG
static void cap_checkProposedAdjacency(Cap *cap, Cap *cap2) {
	if(cap_getCoordinate(cap) != INT32_MAX && cap_getCoordinate(cap2) != INT32_MAX) {
		assert(cap_getSequence(cap) == cap_getSequence(cap2));
		assert(cap_getStrand(cap) == cap_getStrand(cap2));
		assert(cap_getCoordinate(cap) != cap_getCoordinate(cap2));
		if(cap_getCoordinate(cap) < cap_getCoordinate(cap2)) {
			assert(!cap_getSide(cap));
			assert(cap_getSide(cap2));
		}
		else {
			assert(cap_getSide(cap));
			assert(!cap_getSide(cap2));
		}
	}
}
#endif

void cap_makeAdjacent(Cap *cap, Cap *cap2) {
	//We put them both on the same strand, as the strand is not important in the pairing
	//logDebug(">>>> Making adjacency %p -- %p\n", cap, cap2);
	cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
	cap2 = cap_getStrand(cap2) ? cap2 : cap_getReverse(cap2);
#ifdef BEN_DEBUG
	cap_checkProposedAdjacency(cap, cap2); //checks, if they have coordinates, that they lie along a single range.
	assert(cap_getEvent(cap) == cap_getEvent(cap2));
	if (cap_getParent(cap) && cap_getParent(cap2)) {
		assert(event_isDescendant(cap_getEvent(cap_getParent(cap)), cap_getEvent(cap2)));
		assert(event_isDescendant(cap_getEvent(cap_getParent(cap2)), cap_getEvent(cap)));
	}
#endif
	cap_breakAdjacency(cap);
	cap_breakAdjacency(cap2);
	//we ensure we have them right with respect there orientation.
	cap->capContents->adjacency = cap_getOrientation(cap) ? cap2 : cap_getReverse(cap2);
	cap2->capContents->adjacency = cap_getOrientation(cap2) ? cap : cap_getReverse(cap);
}

Cap *cap_getP(Cap *cap, Cap *connectedCap) {
	return connectedCap == NULL ? NULL :
			cap_getOrientation(cap) ? connectedCap : cap_getReverse(connectedCap);
}

Cap *cap_getAdjacency(Cap *cap) {
	return cap_getP(cap, cap->capContents->adjacency);
}

Face *cap_getFace(Cap *cap) {
	return cap->capContents->face;
}

Cap *cap_getParent(Cap *cap) {
	Cap *e = cap->capContents->parent;
	return e == NULL ? NULL :
		cap_getOrientation(cap) ? e : cap_getReverse(e);
}

int32_t cap_getChildNumber(Cap *cap) {
	return cap->capContents->children->length;
}

Cap *cap_getChild(Cap *cap, int32_t index) {
#ifdef BEN_DEBUG
	assert(cap_getChildNumber(cap) > index);
	assert(index >= 0);
#endif
	return cap_getP(cap, cap->capContents->children->list[index]);
}

void cap_makeParentAndChild(Cap *capParent, Cap *capChild) {
	capParent = cap_getPositiveOrientation(capParent);
	capChild = cap_getPositiveOrientation(capChild);
	assert(capChild->capContents->parent == NULL);


	if(!listContains(capParent->capContents->children, capChild)) { //defensive, means second calls will have no effect.
		//logDebug("YYYYYYYYYYY New parent %p->%p\n", capParent, capChild);
#ifdef BEN_DEBUG
		assert(event_isDescendant(cap_getEvent(capParent), cap_getEvent(capChild)));
#endif
		listAppend(capParent->capContents->children, capChild);
	}
	capChild->capContents->parent = capParent;
}

void cap_changeParentAndChild(Cap* newCapParent, Cap* capChild) {
	Cap * oldCapParent = capChild->capContents->parent;
	newCapParent = cap_getPositiveOrientation(newCapParent);
	capChild = cap_getPositiveOrientation(capChild);

	assert(oldCapParent);

	logInfo("Parent without child %p\n", oldCapParent);

	if(!listContains(newCapParent->capContents->children, capChild)) { //defensive, means second calls will have no effect.
		listAppend(newCapParent->capContents->children, capChild);
	}
	listRemove(oldCapParent->capContents->children, capChild);
	capChild->capContents->parent = newCapParent;
}

int32_t cap_isInternal(Cap *cap) {
	return cap_getChildNumber(cap) > 0;
}

void cap_check(Cap *cap) {
	End *end = cap_getEnd(cap);
	assert(end_getInstance(end, cap_getName(cap)) == cap);
	assert(cap_getOrientation(cap) == end_getOrientation(end));

	//If we've built the trees
	if(net_builtTrees(end_getNet(cap_getEnd(cap)))) {
		// checks the cap has a parent which has an ancestral event to the caps event, unless it is the root.
		if(end_getRootInstance(end) == cap) {
			assert(cap_getParent(cap) == NULL);
		}
		else {
			Cap *ancestorCap = cap_getParent(cap);
			assert(ancestorCap != NULL);
			event_isAncestor(cap_getEvent(cap), cap_getEvent(ancestorCap));
			assert(cap_getOrientation(cap) == cap_getOrientation(ancestorCap));
		}
		//Check caps ancestor/descendant links are proper.
		int32_t i;
		for(i=0; i<cap_getChildNumber(cap); i++) {
			Cap *childCap = cap_getChild(cap, i);
			assert(childCap != NULL);
			assert(cap_getParent(childCap) == cap);
		}
	}
	else {
		assert(cap_getParent(cap) == NULL); //No root, so no tree.
	}

	//If stub end checks, there is no attached segment.
	if(end_isStubEnd(end)) {
		assert(cap_getSegment(cap) == NULL);
	}
	else {
		Segment *segment = cap_getSegment(cap);
		if(segment != NULL) {
			assert(cap_getOrientation(cap) == segment_getOrientation(segment));
		}
	}

	//Checks adjacencies are properly linked and have consistent coordinates and the same group
	Cap *cap2 = cap_getAdjacency(cap);
	if(cap2 != NULL) {
		assert(end_getGroup(cap_getEnd(cap2)) == end_getGroup(end)); //check they have the same group.
		assert(cap_getAdjacency(cap2) == cap); //reciprocal connection
		assert(cap_getEvent(cap) == cap_getEvent(cap2)); //common event
		assert(cap_getSequence(cap) == cap_getSequence(cap2)); //common sequence (which may be null)
		assert(cap_getStrand(cap) == cap_getStrand(cap2)); //common strand

		if(cap_getCoordinate(cap) != INT32_MAX) { //if they have a coordinate
			assert(cap_getSide(cap) != cap_getSide(cap2)); //they have to represent an interval
			if(cap_getStrand(cap)) {
				if(!cap_getSide(cap)) {
					assert(cap_getCoordinate(cap) < cap_getCoordinate(cap2));
				}
				else {
					assert(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
				}
			}
			else {
				if(cap_getSide(cap)) {
					assert(cap_getCoordinate(cap) < cap_getCoordinate(cap2));
				}
				else {
					assert(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
				}
			}
		}
	}

	//Checks the reverse
	Cap *rCap = cap_getReverse(cap);
	assert(rCap != NULL);
	assert(cap_getReverse(rCap) == cap);
	assert(cap_getOrientation(cap) == !cap_getOrientation(rCap));
	assert(cap_getEnd(cap) == end_getReverse(cap_getEnd(rCap)));
	assert(cap_getName(cap) == cap_getName(rCap));
	assert(cap_getEvent(cap) == cap_getEvent(rCap));
	assert(cap_getEnd(cap) == end_getReverse(cap_getEnd(rCap)));
	if(cap_getSegment(cap) == NULL) {
		assert(cap_getSegment(rCap) == NULL);
	}
	else {
		assert(cap_getSegment(cap) == segment_getReverse(cap_getSegment(rCap)));
	}
	assert(cap_getSide(cap) == !cap_getSide(rCap));
	assert(cap_getCoordinate(cap) == cap_getCoordinate(rCap));
	assert(cap_getSequence(cap) == cap_getSequence(rCap));
	assert(cap_getStrand(cap) == !cap_getStrand(rCap));
	if(cap_getAdjacency(cap) == NULL) {
		assert(cap_getAdjacency(rCap) == NULL);
	}
	else {
		assert(cap_getReverse(cap_getAdjacency(rCap)) == cap_getAdjacency(cap));
	}
	assert(cap_getFace(cap) == cap_getFace(rCap));
	if(cap_getParent(cap) == NULL) {
		assert(cap_getParent(rCap) == NULL);
	}
	else {
		assert(cap_getParent(cap) == cap_getReverse(cap_getParent(rCap)));
	}
	assert(cap_isInternal(cap) == cap_isInternal(rCap));
	assert(cap_getChildNumber(cap) == cap_getChildNumber(rCap));
	int32_t i;
	for(i=0; i<cap_getChildNumber(cap); i++) {
		assert(cap_getChild(cap, i) == cap_getReverse(cap_getChild(rCap, i)));
	}

	//it is consistent with any copy of end in the nested net, in terms of events, connections and coordinates.
	Net *nestedNet = group_getNestedNet(end_getGroup(end));
	if(nestedNet != NULL) {
		End *childEnd = net_getEnd(nestedNet, end_getName(end));
		assert(childEnd != NULL); //End must be present in child
		//Cap *childCap = end_getInstance(childEnd, cap_getName(cap));
	}
}

/*
 * Private functions.
 */

void cap_setSegment(Cap *cap, Segment *segment) {
	cap->capContents->segment = cap_getOrientation(cap) ? segment : segment_getReverse(segment);
}

void cap_setFace(Cap *cap, Face *face) {
	cap->capContents->face = face;
}

void cap_breakAdjacency(Cap *cap) {
	Cap *cap2;
	cap2 = cap_getAdjacency(cap);
	if(cap2 != NULL) {
		cap2->capContents->adjacency = NULL;
		cap->capContents->adjacency = NULL;
	}
}

/*
 * Serialisation functions.
 */


void cap_writeBinaryRepresentationP(Cap *cap2, int32_t elementType, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(elementType, writeFn);
	binaryRepresentation_writeName(cap_getName(cap2), writeFn);
}

void cap_writeBinaryRepresentation(Cap *cap, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	Cap *cap2;
	if(cap_getCoordinate(cap) == INT32_MAX) {
		binaryRepresentation_writeElementType(CODE_CAP, writeFn);
		binaryRepresentation_writeName(cap_getName(cap), writeFn);
		binaryRepresentation_writeName(event_getName(cap_getEvent(cap)), writeFn);
		binaryRepresentation_writeBool(cap_getSide(cap), writeFn);
	}
	else {
		binaryRepresentation_writeElementType(CODE_CAP_WITH_COORDINATES, writeFn);
		binaryRepresentation_writeName(cap_getName(cap), writeFn);
		binaryRepresentation_writeInteger(cap_getCoordinate(cap), writeFn);
		binaryRepresentation_writeBool(cap_getStrand(cap), writeFn);
		binaryRepresentation_writeBool(cap_getSide(cap), writeFn);
		binaryRepresentation_writeName(sequence_getName(cap_getSequence(cap)), writeFn);
	}
	if((cap2 = cap_getAdjacency(cap)) != NULL) {
		cap_writeBinaryRepresentationP(cap2, CODE_ADJACENCY, writeFn);
	}
	if((cap2 = cap_getParent(cap)) != NULL) {
		cap_writeBinaryRepresentationP(cap2, CODE_PARENT, writeFn);
	}
}

int32_t cap_loadFromBinaryRepresentationP(Cap *cap, void **binaryString, void (*linkFn)(Cap *, Cap *)) {
	Cap *cap2;
	binaryRepresentation_popNextElementType(binaryString);
	cap2 = net_getCap(end_getNet(cap_getEnd(cap)), binaryRepresentation_getName(binaryString));
	if(cap2 != NULL) { //if null we'll make the adjacency when the other end is parsed.
		linkFn(cap2, cap);
		return 0;
	}
	return 1;
}

Cap *cap_loadFromBinaryRepresentation(void **binaryString, End *end) {
	Cap *cap;
	Name name;
	Event *event;
	int32_t coordinate;
	int32_t strand;
	int32_t side;
	Sequence *sequence;

	cap = NULL;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_CAP) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		event = eventTree_getEvent(net_getEventTree(end_getNet(end)), binaryRepresentation_getName(binaryString));
		side = binaryRepresentation_getBool(binaryString);
		cap = cap_construct3(name, event, end, side);
	}
	else if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_CAP_WITH_COORDINATES) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		coordinate = binaryRepresentation_getInteger(binaryString);
		strand = binaryRepresentation_getBool(binaryString);
		side = binaryRepresentation_getBool(binaryString);
		sequence = net_getSequence(end_getNet(end), binaryRepresentation_getName(binaryString));
		cap = cap_construct4(name, end, coordinate, strand, side, sequence);
	}
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_ADJACENCY) {
		cap_loadFromBinaryRepresentationP(cap, binaryString, cap_makeAdjacent);
	}
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_PARENT) {
		assert(cap_loadFromBinaryRepresentationP(cap, binaryString, cap_makeParentAndChild) == 0);
	}

	return cap;
}

Cap *cap_getStaticNameWrapper(Name name) {
	static Cap cap;
	static CapContents capContents;
	cap.capContents = &capContents;
	cap.capContents->instance = name;
	return &cap;
}

void cap_setEvent(Cap *cap, Event *event) {
	cap->capContents->event = event;
}

void cap_setSequence(Cap *cap, Sequence *sequence) {
	cap->capContents->sequence = sequence;
}
