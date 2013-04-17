/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"

#include "cactus_buildFaces.h"

/*
 * Simplify a given face
 */
void buildFaces_simplify(Face * face, Flower * flower) {
    // TODO
    // Here lies NP-completeness
    assert(face_isSimple(face));
}

/*
 * Simplify all the faces in the flower
 */
void buildFaces_simplifyFaces(Flower * flower) {
    Flower_FaceIterator *iter = flower_getFaceIterator(flower);
    Face *face;

    st_logInfo("Simplifying faces\n");

    while ((face = flower_getNextFace(iter)))
	buildFaces_simplify(face, flower);

    flower_destructFaceIterator(iter);
}

/*
 * Create an interpolation between parent and child caps at event
 */
static Cap *buildFaces_interpolateCaps(Cap * parentCap, Cap * childCap,
	Event * event) {
    Cap * newCap = cap_construct(cap_getEnd(childCap), event);
    cap_changeParentAndChild(newCap, childCap);
    cap_makeParentAndChild(parentCap, newCap);
    return newCap;
}

/*
 * Creates an inteprolation half way on the branch between two events
 */
static Event *buildFaces_interpolateEvents(Event* parentEvent,
	Event* childEvent) {
    EventTree * eventTree = event_getEventTree(parentEvent);
    float branchLength = event_getBranchLength(childEvent) / 2;

    return event_construct4(NULL, branchLength, parentEvent, childEvent,
	    eventTree);
}

/*
 * Creates the interpolation of a top node
 */
static Cap *buildFaces_interpolateTopNode(Face * face, int64_t topIndex) {
    Cap *topNode = face_getTopNode(face, topIndex);
    FaceEnd * faceEnd = cap_getTopFaceEnd(topNode); 
    FaceEnd_BottomNodeIterator * bottomNodeIterator = faceEnd_getBottomNodeIterator(faceEnd);
    Cap * bottomNode;
    Cap *adjacency = cap_getPositiveOrientation(cap_getAdjacency(topNode));
    Event *topEvent, *bottomEvent, *interpolatedEvent;
    Cap * derivedEdgeBottomNode = NULL;

    // Search for bottom node which generated the derived edge
    while ((bottomNode = faceEnd_getNextBottomNode(bottomNodeIterator))) {
	if (cap_getTopCap(cap_getPositiveOrientation(bottomNode)) != adjacency) {
	    derivedEdgeBottomNode = bottomNode;
	    break;
	}
    }

    // If none found
    if (derivedEdgeBottomNode == NULL)
	return NULL;

    // Go to the appropriate descent edge
    while (cap_getParent(derivedEdgeBottomNode) != topNode)
	derivedEdgeBottomNode = cap_getParent(derivedEdgeBottomNode);

    // Event interpolation
    topEvent = cap_getEvent(topNode);
    bottomEvent = cap_getEvent(derivedEdgeBottomNode);
    while (event_getParent(bottomEvent) != topEvent)
	bottomEvent = event_getParent(bottomEvent);
    interpolatedEvent = buildFaces_interpolateEvents(topEvent, bottomEvent);

    // Cap interpolation
    return buildFaces_interpolateCaps(topNode, derivedEdgeBottomNode,
	    interpolatedEvent);
}

/*
 * Produces array of interpolations for the top nodes of a face
 */
static Cap **buildFaces_interpolateTopNodes(Face * face) {
    Cap **interpolations = st_calloc(face_getCardinal(face), sizeof(Cap *));
    uint64_t topIndex;

    for (topIndex = 0; topIndex < face_getCardinal(face); topIndex++)
	interpolations[topIndex]
		= buildFaces_interpolateTopNode(face, topIndex);

    return interpolations;
}

/*
 * Connects an interpolated node within a face
 */
static void buildFaces_connectInterpolatedNode(Cap ** interpolations,
	int64_t nodeIndex, Face * face) {
    Cap *node = interpolations[nodeIndex];
    Cap *topNode = face_getTopNode(face, nodeIndex);
    Cap *adjacentNode = NULL;
    int64_t adjacentIndex;

    if (cap_getAdjacency(topNode))
	adjacentNode = cap_getPositiveOrientation(cap_getAdjacency(topNode));

    // If not interpolated or previously connected or top node disconnected
    if (node == NULL || cap_getAdjacency(node) || adjacentNode == NULL)
	return;

    // Look for index of adjacent node
    for (adjacentIndex = 0; adjacentIndex < face_getCardinal(face); adjacentIndex++)
	if (face_getTopNode(face, adjacentIndex) == adjacentNode)
	    break;

#ifndef NDEBUG
    // What if adjacent node is not a top node??
    if (adjacentIndex == face_getCardinal(face))
	abort();
#endif

    // Tie the knot
    cap_makeAdjacent(node, interpolations[adjacentIndex]);
}

/*
 * Connects the interpolated nodes of a face
 */
static void buildFaces_connectInterpolatedNodes(Cap ** interpolations,
	Face * face) {
    int64_t topIndex;

    for (topIndex = 0; topIndex < face_getCardinal(face); topIndex++)
	buildFaces_connectInterpolatedNode(interpolations, topIndex, face);
}

/*
 * Create face from interpolated nodes
 */
static void buildFaces_createInterpolatedFace(Face * face,
	Cap ** interpolations, Flower * flower) {
    int64_t nodeIndex;

    // Create "lower" canonical face
    for (nodeIndex = 0; nodeIndex < face_getCardinal(face); nodeIndex++)
	if (interpolations[nodeIndex] != NULL && cap_getTopFace(interpolations[nodeIndex]) == NULL)
	    buildFaces_reconstructFromCap(interpolations[nodeIndex], flower);

    // Create "upper" duplication faces
    for (nodeIndex = 0; nodeIndex < face_getCardinal(face); nodeIndex++)
	if (cap_getTopFace(face_getTopNode(face, nodeIndex)) == face)
	    buildFaces_reconstructFromCap(face_getTopNode(face, nodeIndex), flower);

}

/*
 * Isolates into regular and trivial faces
 */
void buildFaces_isolate(Face * face, Flower * flower) {
    Cap **interpolations;

    // If uncessary
    if (face_isRegular(face))
	return;

    // Interpolate top nodes
    interpolations = buildFaces_interpolateTopNodes(face);

    // Connect intepolated nodes
    buildFaces_connectInterpolatedNodes(interpolations, face);

    // Create rearrangment face
    buildFaces_createInterpolatedFace(face, interpolations, flower);

    // Cleaning up
    face_destruct(face);
    free(interpolations);

#ifndef NDEBUG
    assert(face_isRegular(face));
#endif
}

/*
 * Isolates all the faces in the flower
 */
void buildFaces_isolateFaces(Flower * flower) {
    Flower_FaceIterator *iter = flower_getFaceIterator(flower);
    Face *face;

    st_logInfo("Isolating faces\n");

    while ((face = flower_getNextFace(iter)))
	buildFaces_isolate(face, flower);

    flower_destructFaceIterator(iter);
}

/*
 * Builds a stub node and creates at the same time an ad hoc end
 */
static Cap * buildFaces_constructStub(Cap * adjacentCap) {
    //Construct the new stub and the new cap..
    Flower * flower = end_getFlower(cap_getEnd(adjacentCap));
    End *newFreeStubEnd = createNewFreeStubEnd(flower);
    Cap *cap = cap_construct(newFreeStubEnd, cap_getEvent(adjacentCap));

    //Now set the group of the new stub end (they should be in the same group)
    End *adjacentCapEnd = cap_getEnd(adjacentCap);
    Group *group = end_getGroup(adjacentCapEnd);
    end_setGroup(newFreeStubEnd, group);

    //Now make adjacent
    cap_makeAdjacent(cap, adjacentCap);

    return cap;

    /*
     The 0 argument to the end constructor is a bool saying the stub end is 'free',
     ie. not necessarily inherited from the parent (though it can be), and not part
     of the reference structure. I.e. some stub ends are 'attached' - opposite of free,
     representing the case where we know what happened to the other end of the stub
     (i.e. if it is a block end at a higher level or the defined end of a sequence defined
     at the top level).
     */
}

/*
 * Adds a top node with a simgle bottom node
 */
/*
 * Engineers a node so that a regular face has an even number of
 * top/bottom node pairs
 */
static void buildFaces_engineerCaps(Face * face, Flower * flower) {
    Cap *nonAdjacent = NULL;
    Cap * nonDerivedBottomNode = NULL;
    Cap *X, *Xprime;
    Cap *parent;
    Face_FaceEndIterator * faceEndIterator;
    FaceEnd * faceEnd = NULL;
    FaceEnd_BottomNodeIterator * bottomNodeIterator;
    Cap *bottomNode = NULL, *adjacency = NULL;

    // Look for terminal nodes in face according to their nature
    faceEndIterator = face_getFaceEndIterator(face);
    while ((faceEnd = face_getNextFaceEnd(faceEndIterator))) {
	if (!cap_getAdjacency(faceEnd_getTopNode(faceEnd))) {
	    // Found top node that has no adjacency
	    nonAdjacent = faceEnd_getTopNode(faceEnd);
	    if (nonDerivedBottomNode)
		break;
	} else {
	    // Look for top node with node derived lifted edge
	    adjacency = cap_getPositiveOrientation(cap_getAdjacency(faceEnd_getTopNode(faceEnd)));
	    bottomNodeIterator = faceEnd_getBottomNodeIterator(faceEnd);
	    while((bottomNode = faceEnd_getNextBottomNode(bottomNodeIterator)))
		if (cap_getTopCap(cap_getPositiveOrientation(cap_getAdjacency(bottomNode))) != adjacency)
		    break;

	    faceEnd_destructBottomNodeIterator(bottomNodeIterator);

	    if (!bottomNode) {
		bottomNodeIterator = faceEnd_getBottomNodeIterator(faceEnd);
		nonDerivedBottomNode = faceEnd_getNextBottomNode(bottomNodeIterator);
		faceEnd_destructBottomNodeIterator(bottomNodeIterator);
		if (nonAdjacent)
		    break;
	    }
	}
    }
    face_destructFaceEndIterator(faceEndIterator);

    // Engineer appropriate nodes
    X = buildFaces_constructStub(nonAdjacent);
    Xprime = cap_construct(cap_getEnd(X), cap_getEvent(bottomNode));
    cap_makeAdjacent(Xprime, nonDerivedBottomNode);
    cap_makeParentAndChild(X, Xprime);

    if (cap_getParent(nonAdjacent)) {
	parent = cap_construct(cap_getEnd(X), cap_getEvent(cap_getParent(
		nonAdjacent)));
	cap_makeParentAndChild(parent, X);
    }

    // Correct face:
    buildFaces_reconstructFromCap(nonAdjacent, flower);
    face_destruct(face);
}

/*
 * Connects the ends of an even length regular path to from
 * a canonical cycle
 */
static void buildFaces_close(Face * face) {
    FaceEnd * faceEnd;
    Face_FaceEndIterator * faceEndIterator;
    FaceEnd_BottomNodeIterator * bottomNodeIterator;
    Cap * adjacency, *bottomNode;
    Cap *nonAdjacent1 = NULL;
    Cap *nonAdjacent2 = NULL;
    Cap *nonDerived1 = NULL;
    Cap *nonDerived2 = NULL;

    // Look for terminal nodes in face according to their nature
    faceEndIterator = face_getFaceEndIterator(face);
    while ((faceEnd = face_getNextFaceEnd(faceEndIterator))) {
	if (!cap_getAdjacency(faceEnd_getTopNode(faceEnd))) {
	    // Found top node that has no adjacency
	    if (nonAdjacent1 == NULL)
		nonAdjacent1 = faceEnd_getTopNode(faceEnd);
	    else {
		assert(!nonAdjacent2);
		nonAdjacent2 = faceEnd_getTopNode(faceEnd);
	    }
	} else {
	    // Look for top node with node derived lifted edge
	    adjacency = cap_getPositiveOrientation(cap_getAdjacency(faceEnd_getTopNode(faceEnd)));
	    bottomNodeIterator = faceEnd_getBottomNodeIterator(faceEnd);
	    while((bottomNode = faceEnd_getNextBottomNode(bottomNodeIterator)))
		if (cap_getTopCap(cap_getPositiveOrientation(cap_getAdjacency(bottomNode))) != adjacency)
		    break;

	    faceEnd_destructBottomNodeIterator(bottomNodeIterator);

	    if (!bottomNode) {
		bottomNodeIterator = faceEnd_getBottomNodeIterator(faceEnd);
		if (nonDerived1)
		    nonDerived1 = faceEnd_getNextBottomNode(bottomNodeIterator);
		else {
		    assert(!nonDerived2);
		    nonDerived2 = faceEnd_getNextBottomNode(bottomNodeIterator);

		}
		faceEnd_destructBottomNodeIterator(bottomNodeIterator);
	    }
	}
    }
    face_destructFaceEndIterator(faceEndIterator);

    // What if some nodes with no adjacency and others with no derived edge?
    assert(!nonAdjacent1 || !nonDerived1);

    // If the terminal nodes lack an adjacency edge
    if (nonAdjacent2)
	cap_makeAdjacent(nonAdjacent1, nonAdjacent2);
    else if (nonDerived2)
	cap_makeAdjacent(nonDerived1, nonDerived2);
}

/*
 * Canonizes face into a regular cycle
 */
void buildFaces_canonize(Face * face, Flower * flower) {
    if (face_getCardinal(face) % 2 == 0)
	buildFaces_engineerCaps(face, flower);

#ifndef NDEBUG
    assert(face_isCanonical(face));
#endif

    buildFaces_close(face);
}

/*
 * Canonizes all the faces in the flower
 */
void buildFaces_canonizeFaces(Flower * flower) {
    Flower_FaceIterator *iter = flower_getFaceIterator(flower);
    Face *face;

    st_logInfo("Canonizing faces\n");

    while ((face = flower_getNextFace(iter)))
	if (!face_isCanonical(face))
	    buildFaces_canonize(face, flower);

    flower_destructFaceIterator(iter);
}

/*
 * The works: create, regularize and canonize faces in flower
 */
void buildFaces_buildAndProcessFaces(Flower * flower) {
    buildFaces_simplifyFaces(flower);
    buildFaces_isolateFaces(flower);
    buildFaces_canonizeFaces(flower);
}

//Todo - fine me a home.
End *createNewFreeStubEnd(Flower *flower) {
    End *newFreeStubEnd = end_construct(0, flower);
    assert(flower_getGroupNumber(flower) == 1);
    end_setGroup(newFreeStubEnd, flower_getFirstGroup(flower));
    return newFreeStubEnd;
}

