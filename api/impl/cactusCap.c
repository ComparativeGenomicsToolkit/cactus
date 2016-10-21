/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic cap functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

Cap *cap_construct(End *end, Event *event) {
    return cap_construct3(cactusDisk_getUniqueID(flower_getCactusDisk(end_getFlower(end))), event, end);
}

Cap *cap_construct3(Name instance, Event *event, End *end) {
    assert(end != NULL);
    assert(event != NULL);
    assert(instance != NULL_NAME);
    Cap *cap;

    cap = st_malloc(sizeof(Cap));
    cap->capContents = st_malloc(sizeof(CapContents));
    cap->rCap = st_malloc(sizeof(Cap));
    cap->rCap->rCap = cap;
    cap->rCap->capContents = cap->capContents;

    cap->end = end;
    cap->rCap->end = end_getReverse(end);

    cap->capContents->instance = instance;
    cap->capContents->coordinate = INT64_MAX;
    cap->capContents->sequence = NULL;
    cap->capContents->adjacency = NULL;
    cap->capContents->adjacency2 = NULL;
    cap->capContents->face = NULL;
    cap->capContents->segment = NULL;
    cap->capContents->parent = NULL;
    cap->capContents->children = constructEmptyList(0, NULL);
    cap->capContents->event = event;
    cap->capContents->strand = end_getOrientation(end);

    end_addInstance(end, cap);
    flower_addCap(end_getFlower(end), cap);
    return cap;
}

Cap *cap_construct2(End *end, int64_t coordinate, bool strand, Sequence *sequence) {
    return cap_construct4(cactusDisk_getUniqueID(flower_getCactusDisk(end_getFlower(end))), end, coordinate, strand,
            sequence);
}

Cap *cap_construct4(Name instance, End *end, int64_t coordinate, bool strand, Sequence *sequence) {
    Cap *cap;
    cap = cap_construct3(instance, sequence_getEvent(sequence), end);
    cap->capContents->coordinate = coordinate;
    cap->capContents->strand = cap_getOrientation(cap) ? strand : !strand;
    cap->capContents->sequence = sequence;
    return cap;
}

Cap *cap_construct5(Event *event, End *end) {
    return cap_construct3(cactusDisk_getUniqueID(flower_getCactusDisk(end_getFlower(end))), event, end);
}

void cap_setCoordinates(Cap *cap, int64_t coordinate, bool strand, Sequence *sequence) {
    cap->capContents->coordinate = coordinate;
    cap->capContents->strand = cap_getOrientation(cap) ? strand : !strand;
    cap->capContents->sequence = sequence;
}

Cap *cap_copyConstruct(End *end, Cap *cap) {
    assert(end_getName(cap_getEnd(cap)) == end_getName(end));
    assert(end_getSide(end) == cap_getSide(cap));
    Event *event;
    Name sequenceName;
    Sequence *sequence;

    Flower *flower = end_getFlower(end);
    if (cap_getSequence(cap) != NULL) { //cap_getCoordinate(cap) != INT64_MAX) {
        sequenceName = sequence_getName(cap_getSequence(cap));
        sequence = flower_getSequence(flower, sequenceName);
        if (sequence == NULL) { //add sequence to the flower.
            sequence = sequence_construct(cactusDisk_getMetaSequence(flower_getCactusDisk(flower), sequenceName),
                    flower);
            assert(sequence != NULL);
        }
        return cap_construct4(cap_getName(cap), end, cap_getCoordinate(cap), cap_getStrand(cap), sequence);
    } else {
        event = eventTree_getEvent(flower_getEventTree(flower), event_getName(cap_getEvent(cap)));
        assert(event != NULL);
        Cap *cap2 = cap_construct3(cap_getName(cap), event, end);
        cap_setCoordinates(cap2, cap_getCoordinate(cap), cap_getStrand(cap), NULL);
        return cap2;
    }
}

void cap_destruct(Cap *cap) {
    //Remove from end.
    end_removeInstance(cap_getEnd(cap), cap);
    flower_removeCap(end_getFlower(cap_getEnd(cap)), cap);

    // Remove parent->child link from parent (if any).
    Cap *capParent = cap->capContents->parent;
    if (capParent != NULL) {
        listRemove(capParent->capContents->children, cap);
    }

    // Remove child->parent link from children (if any).
    for (int64_t i = 0; i < cap->capContents->children->length; i++) {
        Cap *capChild = cap->capContents->children->list[i];
        capChild->capContents->parent = NULL;
    }

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
    return cap_getOrientation(cap) ? cap->capContents->segment
            : (cap->capContents->segment != NULL ? segment_getReverse(cap->capContents->segment) : NULL);
}

Cap *cap_getOtherSegmentCap(Cap *cap) {
    if (!end_isBlockEnd(cap_getEnd(cap))) {
        assert(cap_getSegment(cap) == NULL);
        return NULL;
    }
    Segment *segment = cap_getSegment(cap);
    assert(segment != NULL);
    Cap *otherCap = cap_getSide(cap) ? segment_get3Cap(segment) : segment_get5Cap(segment);
    assert(cap != otherCap);
    return otherCap;
}

int64_t cap_getCoordinate(Cap *cap) {
    return cap->capContents->coordinate;
}

bool cap_getStrand(Cap *cap) {
    return cap_getOrientation(cap) ? cap->capContents->strand : !cap->capContents->strand;
}

bool cap_getSide(Cap *cap) {
    return end_getSide(cap_getEnd(cap));
}

Sequence *cap_getSequence(Cap *cap) {
    return cap->capContents->sequence;
}

void cap_makeAdjacent(Cap *cap, Cap *cap2) {
    //We put them both on the same strand, as the strand is important in the pairing
    cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
    cap2 = cap_getStrand(cap2) ? cap2 : cap_getReverse(cap2);
    assert(cap != cap2);
    assert(cap_getEvent(cap) == cap_getEvent(cap2));
    cap_breakAdjacency(cap);
    cap_breakAdjacency(cap2);
    //we ensure we have them right with respect there orientation.
    cap->capContents->adjacency = cap_getOrientation(cap) ? cap2 : cap_getReverse(cap2);
    cap2->capContents->adjacency = cap_getOrientation(cap2) ? cap : cap_getReverse(cap);
}

Cap *cap_getP(Cap *cap, Cap *connectedCap) {
    return connectedCap == NULL ? NULL : cap_getOrientation(cap) ? connectedCap : cap_getReverse(connectedCap);
}

Cap *cap_getAdjacency(Cap *cap) {
    return cap_getP(cap, cap->capContents->adjacency);
}

Cap *cap_getTopCap(Cap *cap) {
    if (cap_getAdjacency(cap) == NULL || end_getRootInstance(cap_getEnd(cap)) == cap) {
        return NULL;
    }
    Cap *cap2 = cap_getParent(cap);
    assert(cap2 != NULL);
    while (1) {
        if (cap_getAdjacency(cap2) != NULL) {
            return cap2;
        }
        if (cap_getParent(cap2) == NULL) {
            assert(end_getRootInstance(cap_getEnd(cap2)) == cap2);
            return cap2;
        }
        cap2 = cap_getParent(cap2);
    }
    return NULL;
}

Face *cap_getTopFace(Cap *cap) {
    return cap->capContents->face;
}

FaceEnd *cap_getTopFaceEnd(Cap *cap) {
    Face *face = cap_getTopFace(cap);
    if (face != NULL) {
        FaceEnd *faceEnd = face_getFaceEndForTopNode(face, cap);
        assert(faceEnd != NULL);
        return faceEnd;
    }
    return NULL;
}

FaceEnd *cap_getBottomFaceEnd(Cap *cap) {
    Cap *topNode = cap_getTopCap(cap);
    if (topNode != NULL) {
        return cap_getTopFaceEnd(topNode);
    }
    return NULL;
}

Cap *cap_getParent(Cap *cap) {
    Cap *e = cap->capContents->parent;
    return e == NULL ? NULL : cap_getOrientation(cap) ? e : cap_getReverse(e);
}

int64_t cap_getChildNumber(Cap *cap) {
    return cap->capContents->children->length;
}

Cap *cap_getChild(Cap *cap, int64_t index) {
    assert(cap_getChildNumber(cap) > index);
    assert(index >= 0);
    return cap_getP(cap, cap->capContents->children->list[index]);
}

void cap_makeParentAndChild(Cap *capParent, Cap *capChild) {
    capParent = cap_getPositiveOrientation(capParent);
    capChild = cap_getPositiveOrientation(capChild);
    assert(capChild->capContents->parent == NULL);

    if (!listContains(capParent->capContents->children, capChild)) { //defensive, means second calls will have no effect.
        assert(event_isDescendant(cap_getEvent(capParent), cap_getEvent(capChild)));
        listAppend(capParent->capContents->children, capChild);
    }
    capChild->capContents->parent = capParent;
}

void cap_changeParentAndChild(Cap* newCapParent, Cap* capChild) {
    Cap * oldCapParent = cap_getParent(capChild);
    newCapParent = cap_getPositiveOrientation(newCapParent);
    capChild = cap_getPositiveOrientation(capChild);
    assert(oldCapParent);
    if (!listContains(newCapParent->capContents->children, capChild)) { //defensive, means second calls will have no effect.
        listAppend(newCapParent->capContents->children, capChild);
    }
    listRemove(oldCapParent->capContents->children, capChild);
    capChild->capContents->parent = newCapParent;
}

bool cap_isInternal(Cap *cap) {
    return cap_getChildNumber(cap) > 0;
}

void cap_check(Cap *cap) {
    End *end = cap_getEnd(cap);
    cactusCheck(end_getInstance(end, cap_getName(cap)) == cap);
    cactusCheck(cap_getOrientation(cap) == end_getOrientation(end));
    cactusCheck(end_getSide(end) == cap_getSide(cap)); //This is critical, it ensures
    //that we have a consistently oriented set of caps in an end.

    //If we've built the trees
    if (flower_builtTrees(end_getFlower(cap_getEnd(cap)))) {
        // checks the cap has a parent which has an ancestral event to the caps event, unless it is the root.
        cactusCheck(end_getRootInstance(end) != NULL);
        if (end_getRootInstance(end) == cap) {
            cactusCheck(cap_getParent(cap) == NULL);
        } else {
            Cap *ancestorCap = cap_getParent(cap);
            cactusCheck(ancestorCap != NULL);
            event_isAncestor(cap_getEvent(cap), cap_getEvent(ancestorCap));
            cactusCheck(cap_getOrientation(cap) == cap_getOrientation(ancestorCap));
        }
        //Check caps ancestor/descendant links are proper.
        int64_t i;
        for (i = 0; i < cap_getChildNumber(cap); i++) {
            Cap *childCap = cap_getChild(cap, i);
            cactusCheck(childCap != NULL);
            cactusCheck(cap_getParent(childCap) == cap);
        }
    } else {
        cactusCheck(cap_getParent(cap) == NULL); //No root, so no tree.
    }

    //If stub end checks, there is no attached segment.
    if (end_isStubEnd(end)) {
        cactusCheck(cap_getSegment(cap) == NULL);
    } else {
        Segment *segment = cap_getSegment(cap);
        if (segment != NULL) {
            cactusCheck(cap_getOrientation(cap) == segment_getOrientation(segment));
        }
    }

    //Checks adjacencies are properly linked and have consistent coordinates and the same group
    Cap *cap2 = cap_getAdjacency(cap);
    if (cap2 != NULL) {
        cactusCheck(end_getGroup(cap_getEnd(cap2)) == end_getGroup(end)); //check they have the same group.
        cactusCheck(cap_getAdjacency(cap2) == cap); //reciprocal connection
        cactusCheck(cap_getEvent(cap) == cap_getEvent(cap2)); //common event
        cactusCheck(cap_getStrand(cap) == cap_getStrand(cap2)); //common strand
        cactusCheck(cap_getSequence(cap) == cap_getSequence(cap2)); //common sequence (which may be null)

        if (cap_getCoordinate(cap) != INT64_MAX) { //if they have a coordinate
            cactusCheck(cap_getSide(cap) != cap_getSide(cap2)); //they have to represent an interval
            if (cap_getStrand(cap)) {
                if (!cap_getSide(cap)) {
                    cactusCheck(cap_getCoordinate(cap) < cap_getCoordinate(cap2));
                } else {
                    cactusCheck(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
                }
            } else {
                if (cap_getSide(cap)) {
                    cactusCheck(cap_getCoordinate(cap) < cap_getCoordinate(cap2));
                } else {
                    cactusCheck(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
                }
            }
        } else {
            cactusCheck(cap_getCoordinate(cap2) == INT64_MAX);
        }
    }

    //Checks the reverse
    Cap *rCap = cap_getReverse(cap);
    cactusCheck(rCap != NULL);
    cactusCheck(cap_getReverse(rCap) == cap);
    cactusCheck(cap_getOrientation(cap) == !cap_getOrientation(rCap));
    cactusCheck(cap_getEnd(cap) == end_getReverse(cap_getEnd(rCap)));
    cactusCheck(cap_getName(cap) == cap_getName(rCap));
    cactusCheck(cap_getEvent(cap) == cap_getEvent(rCap));
    cactusCheck(cap_getEnd(cap) == end_getReverse(cap_getEnd(rCap)));
    if (cap_getSegment(cap) == NULL) {
        cactusCheck(cap_getSegment(rCap) == NULL);
    } else {
        cactusCheck(cap_getSegment(cap) == segment_getReverse(cap_getSegment(rCap)));
    }
    cactusCheck(cap_getSide(cap) == !cap_getSide(rCap));
    cactusCheck(cap_getCoordinate(cap) == cap_getCoordinate(rCap));
    cactusCheck(cap_getSequence(cap) == cap_getSequence(rCap));
    cactusCheck(cap_getStrand(cap) == !cap_getStrand(rCap));
    if (cap_getAdjacency(cap) == NULL) {
        cactusCheck(cap_getAdjacency(rCap) == NULL);
    } else {
        cactusCheck(cap_getReverse(cap_getAdjacency(rCap)) == cap_getAdjacency(cap));
    }
    cactusCheck(cap_getTopFace(cap) == cap_getTopFace(rCap));
    if (cap_getParent(cap) == NULL) {
        cactusCheck(cap_getParent(rCap) == NULL);
    } else {
        cactusCheck(cap_getParent(cap) == cap_getReverse(cap_getParent(rCap)));
    }
    cactusCheck(cap_isInternal(cap) == cap_isInternal(rCap));
    cactusCheck(cap_getChildNumber(cap) == cap_getChildNumber(rCap));
    int64_t i;
    for (i = 0; i < cap_getChildNumber(cap); i++) {
        cactusCheck(cap_getChild(cap, i) == cap_getReverse(cap_getChild(rCap, i)));
    }

    //it is consistent with any copy of end in the nested flower, in terms of events, connections and coordinates.
    Flower *nestedFlower = group_getNestedFlower(end_getGroup(end));
    if (nestedFlower != NULL) {
        End *childEnd = flower_getEnd(nestedFlower, end_getName(end));
        cactusCheck(childEnd != NULL); //End must be present in child
        //Cap *childCap = end_getInstance(childEnd, cap_getName(cap));
    }
}

/*
 * Private functions.
 */

void cap_setSegment(Cap *cap, Segment *segment) {
    cap->capContents->segment = cap_getOrientation(cap) ? segment : segment_getReverse(segment);
}

void cap_setTopFace(Cap *cap, Face *face) {
    cap->capContents->face = face;
}

void cap_breakAdjacency(Cap *cap) {
    Cap *cap2;
    cap2 = cap_getAdjacency(cap);
    if (cap2 != NULL) {
        cap2->capContents->adjacency = NULL;
        cap->capContents->adjacency = NULL;
    }
}

/*
 * Serialisation functions.
 */

void cap_writeBinaryRepresentationP(Cap *cap2, int64_t elementType,
        void(*writeFn)(const void * ptr, size_t size, size_t count)) {
    binaryRepresentation_writeElementType(elementType, writeFn);
    binaryRepresentation_writeName(cap_getName(cap2), writeFn);
}

void cap_writeBinaryRepresentation(Cap *cap, void(*writeFn)(const void * ptr, size_t size, size_t count)) {
    Cap *cap2;
    if (cap_getCoordinate(cap) == INT64_MAX) {
        binaryRepresentation_writeElementType(CODE_CAP, writeFn);
        binaryRepresentation_writeName(cap_getName(cap), writeFn);
        binaryRepresentation_writeBool(cap_getStrand(cap), writeFn);
        binaryRepresentation_writeName(event_getName(cap_getEvent(cap)), writeFn);
    } else if (cap_getSequence(cap) != NULL) {
        binaryRepresentation_writeElementType(CODE_CAP_WITH_COORDINATES, writeFn);
        binaryRepresentation_writeName(cap_getName(cap), writeFn);
        binaryRepresentation_writeInteger(cap_getCoordinate(cap), writeFn);
        binaryRepresentation_writeBool(cap_getStrand(cap), writeFn);
        binaryRepresentation_writeName(sequence_getName(cap_getSequence(cap)), writeFn);
    } else {
        binaryRepresentation_writeElementType(CODE_CAP_WITH_COORDINATES_BUT_NO_SEQUENCE, writeFn);
        binaryRepresentation_writeName(cap_getName(cap), writeFn);
        binaryRepresentation_writeInteger(cap_getCoordinate(cap), writeFn);
        binaryRepresentation_writeBool(cap_getStrand(cap), writeFn);
        binaryRepresentation_writeName(event_getName(cap_getEvent(cap)), writeFn);
    }
    if ((cap2 = cap_getAdjacency(cap)) != NULL) {
        cap_writeBinaryRepresentationP(cap2, CODE_ADJACENCY, writeFn);
    }
    if ((cap2 = cap_getParent(cap)) != NULL) {
        cap_writeBinaryRepresentationP(cap2, CODE_PARENT, writeFn);
    }
}

bool cap_loadFromBinaryRepresentationP(Cap *cap, void **binaryString, void(*linkFn)(Cap *, Cap *)) {
    Cap *cap2;
    binaryRepresentation_popNextElementType(binaryString);
    cap2 = flower_getCap(end_getFlower(cap_getEnd(cap)), binaryRepresentation_getName(binaryString));
    if (cap2 != NULL) { //if null we'll make the adjacency when the other end is parsed.
        linkFn(cap2, cap);
        return 0;
    }
    return 1;
}

void cap_loadFromBinaryRepresentationP2(void **binaryString, Cap *cap) {
    if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_ADJACENCY) {
        cap_loadFromBinaryRepresentationP(cap, binaryString, cap_makeAdjacent);
    }
    if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_PARENT) {
        bool i = cap_loadFromBinaryRepresentationP(cap, binaryString, cap_makeParentAndChild);
        (void) i;
        assert(i == 0);
    }
}

Cap *cap_loadFromBinaryRepresentation(void **binaryString, End *end) {
    Cap *cap;
    Name name;
    Event *event;
    int64_t coordinate;
    int64_t strand;
    Sequence *sequence;

    cap = NULL;
    if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_CAP) {
        binaryRepresentation_popNextElementType(binaryString);
        name = binaryRepresentation_getName(binaryString);
        strand = binaryRepresentation_getBool(binaryString);
        event = eventTree_getEvent(flower_getEventTree(end_getFlower(end)), binaryRepresentation_getName(binaryString));
        cap = cap_construct3(name, event, end);
        cap_setCoordinates(cap, INT64_MAX, strand, NULL); //Hacks
        cap_loadFromBinaryRepresentationP2(binaryString, cap);
    } else if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_CAP_WITH_COORDINATES) {
        binaryRepresentation_popNextElementType(binaryString);
        name = binaryRepresentation_getName(binaryString);
        coordinate = binaryRepresentation_getInteger(binaryString);
        strand = binaryRepresentation_getBool(binaryString);
        sequence = flower_getSequence(end_getFlower(end), binaryRepresentation_getName(binaryString));
        cap = cap_construct4(name, end, coordinate, strand, sequence);
        cap_loadFromBinaryRepresentationP2(binaryString, cap);
    } else if (binaryRepresentation_peekNextElementType(*binaryString) == CODE_CAP_WITH_COORDINATES_BUT_NO_SEQUENCE) {
        binaryRepresentation_popNextElementType(binaryString);
        name = binaryRepresentation_getName(binaryString);
        coordinate = binaryRepresentation_getInteger(binaryString);
        strand = binaryRepresentation_getBool(binaryString);
        event = eventTree_getEvent(flower_getEventTree(end_getFlower(end)), binaryRepresentation_getName(binaryString));
        cap = cap_construct3(name, event, end);
        cap_setCoordinates(cap, coordinate, strand, NULL);
        cap_loadFromBinaryRepresentationP2(binaryString, cap);
    }

    return cap;
}

void cap_setEvent(Cap *cap, Event *event) {
    cap->capContents->event = event;
}

void cap_setSequence(Cap *cap, Sequence *sequence) {
    cap->capContents->sequence = sequence;
}
