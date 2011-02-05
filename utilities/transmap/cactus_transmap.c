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

#ifndef true
#define true 1
#endif
#ifndef false
#define false 0
#endif
#ifndef maybe
#define maybe -1
#endif

// Sub for recursion
static int32_t transmap_getTotalDistanceBetweenCaps(Cap * A, Cap * B);

/*
 * Computes the current distance between two adjacent
 * cap in the nested sub-graph
 */
static int32_t transmap_getTotalDistanceAtAdjacency(Cap * A) {
	Flower * flower = group_getNestedFlower(end_getGroup(cap_getEnd(A)));
	Cap * B = cap_getAdjacency(A);
	Cap * childA, * childB;

	if (!flower || !B)
		return 0;

	childA = flower_getCap(flower, cap_getName(A));
	childB = flower_getCap(flower, cap_getName(B));

#ifdef BEN_DEBUG
	assert(childA);
	assert(childB);
#endif

	return transmap_getTotalDistanceBetweenCaps(childA, childB);

}

/* 
 * Computes the current distance between two caps
 * Returns -1 if no simple path was found
 */
static int32_t transmap_getTotalDistanceBetweenCaps(Cap * A, Cap * B) {
	int32_t remaining_length;
	int32_t length;
	
	// End of recursion scenarios
	if (A == B)
		return 0;
	if (A == NULL)
		return -1;

	// Step forward
	if (cap_getSide(A)) {
		length = segment_getLength(cap_getSegment(A));
		A = cap_getOtherSegmentCap(A);
	} else {  
		length = transmap_getTotalDistanceAtAdjacency(A);
		A = cap_getAdjacency(A);
	}

	// Recursion to next step
	remaining_length = transmap_getTotalDistanceBetweenCaps(A, B);

	if (remaining_length == -1)
		return -1;
	else 
		return  segment_getLength(cap_getSegment(A)) + remaining_length;
}

// Stub for recursion
static int32_t transmap_sampleOrderAndOrientationAtEvent(Event * E, Cap * A, Cap * B, int32_t * allowed_distance);

/* 
 * Projects the stochastic search to a nested flower
 */
static int32_t transmap_sampleOrderAndOrientationAtNestedFlower(Event * E, Cap * A, int * allowed_distance) {
	Flower * flower = group_getNestedFlower(end_getGroup(cap_getEnd(A)));
	Cap * B, * childA, * childB;
	Event * childE;

	if (!flower)
		return true;
	if (!A) 
		return false;
	B = cap_getAdjacency(A);
	if (!B)
		return false;

	childA = flower_getCap(flower, cap_getName(A));
	childB = flower_getCap(flower, cap_getName(B));
	childE = eventTree_getEvent(flower_getEventTree(flower), event_getName(E));

#ifdef BEN_DEBUG
	assert(childA);
	assert(childB);
	assert(childE);
#endif

	switch (transmap_sampleOrderAndOrientationAtEvent(childE, childA, childB, allowed_distance)) {
	case true:
		return true;
	case false:
		return false;
	default:
		return (transmap_sampleOrderAndOrientationAtEvent(childE, cap_getReverse(childB), cap_getReverse(childA), allowed_distance) == maybe);
	}
}

/* 
 * Stochastic decision whether to jump a face or not
 */
static int32_t transmap_goByAncestralPath(Face * face, Event * E) {
	Event * topEvent, * bottomEvent, * tmpEvent;
	float totalTimeSpan = 0;
	float partialTimeSpan = 0;
	int32_t index;

	// Get the framing events of the face
	for (index = 0; index < face_getCardinal(face); index++) {
		tmpEvent = cap_getEvent(face_getTopNode(face, index));
		if (index == 0 || event_isAncestor(topEvent, tmpEvent))
			topEvent = tmpEvent;

#ifdef BEN_DEBUG
		// Face should be canonical by now...
		assert(face_getBottomNodeNumber(face, index) == 1);
#endif

		tmpEvent = cap_getEvent(face_getBottomNode(face, index, 0));
		if (index == 0 || event_isAncestor(tmpEvent, bottomEvent))
			bottomEvent = tmpEvent;
	}

	// If event out of bounds	
	if (event_isAncestor(E, topEvent) || E == topEvent)
		return true;
	if (event_isAncestor(bottomEvent, E) || E == bottomEvent)
		return false;

	// Measure time spans
	for (tmpEvent = bottomEvent; event_isAncestor(E, tmpEvent); tmpEvent = event_getParent(tmpEvent))
		partialTimeSpan += event_getBranchLength(tmpEvent);
	for (tmpEvent = bottomEvent; event_isAncestor(topEvent, tmpEvent); tmpEvent = event_getParent(tmpEvent))
		totalTimeSpan += event_getBranchLength(tmpEvent);

	// Random decision
	if (totalTimeSpan == 0)
		return false;
	else
		return (rand() / RAND_MAX) > (partialTimeSpan / totalTimeSpan);
}

/*
 * Returns the number of highest attached caps below cap A
 */
static int32_t transmap_findBottomCapNumberBelow(Cap * A) {
	int32_t total = 0;
	int32_t childIndex;

	// End of recursion
	if (cap_getAdjacency(A))
		return 1;

	// Recursion
	for (childIndex = 0; childIndex < cap_getChildNumber(A); childIndex++)
		total += transmap_findBottomCapNumberBelow(cap_getChild(A, childIndex));

	return total;
}

/*
 * Returns the number of highest attached cap below or at event E
 */
static int32_t transmap_findBottomCapNumber(Event * E, Cap * A) {
	Cap * tmp = A;

	// Try upwards
	for (tmp = A; tmp && !event_isAncestor(cap_getEvent(tmp), E); tmp = cap_getParent(tmp))
		if (cap_getAdjacency(tmp))
			return 1;

	// Still no luck, OK, going down... (risk of multiple hits)
	return transmap_findBottomCapNumberBelow(A);
}

/*
 * Reurns pointer to lowest attached cap above or at event E
 */
static Cap* transmap_findTopCap(Event * E, Cap * A) {
	Cap * tmp;

	for (tmp = A; tmp; tmp = cap_getParent(tmp))
		if (cap_getAdjacency(tmp) && !event_isAncestor(E, cap_getEvent(tmp)))
			return tmp;

	return NULL;
}

/*
 * Returns a pointer to one of the highest attached caps below cap A
 * Normally you only run this function when you know there is only one
 * such cap.
 */
static Cap* transmap_findBottomCapBelow(Cap * A) {
	int32_t childIndex;
	Cap * result;

	// End of recursion
	if (cap_getAdjacency(A))
		return A;

	// Recursion
	for (childIndex = 0; childIndex < cap_getChildNumber(A); childIndex++)
		if ((result = transmap_findBottomCapBelow(cap_getChild(A, childIndex))));
			return result;

	return NULL;
}

/*
 * Reurns pointer to highest attached cap below E
 */
static Cap* transmap_findBottomCap(Event * E, Cap * A) {
	Cap * tmp;

	// Try upwards
	for (tmp = A; tmp && !event_isAncestor(cap_getEvent(tmp), E); tmp = cap_getParent(tmp))
		if (cap_getAdjacency(tmp))
			return tmp;

	// Did not work? That's all right, look around
	if (cap_getAdjacency(A))
		return A;


	return NULL;
}

/*
 * Tries to connect two caps at the time of event
 * Returns true if success
 */
static int32_t transmap_sampleOrderAndOrientationAtEvent(Event * E, Cap * A, Cap * B, int32_t * allowed_distance) {
	Face * face;
	int32_t index;
	Cap * topCap;
	int32_t bottomCapNumber;

	// Termination of recursion
	if (cap_getEnd(A) == cap_getEnd(B))
		return true;
	if (A == NULL)
		return false;

	// Synchronicity
#ifdef BEN_DEBUG
	assert(event_isAncestor(cap_getEvent(A), E) 
	       || event_isAncestor(E, cap_getEvent(A))
	       || cap_getEvent(A) == E);
#endif

	if (event_isAncestor(cap_getEvent(A), E)) {
#ifdef BEN_DEBUG
		assert(cap_getSide(A));
		assert(cap_getTopFace(A));
#endif
		face = cap_getTopFace(A);
		for (index = 0; index < face_getCardinal(face); index++)
			if (face_getTopNode(face, index) == cap_getPositiveOrientation(A))
				break;
#ifdef BEN_DEBUG
		assert(index < face_getCardinal(face));
		assert(face_getBottomNodeNumber(face, index) == 1);
#endif
		if (cap_getOrientation(A))
			A = face_getBottomNode(face, index, 0);
		else 
			A = cap_getReverse(face_getBottomNode(face, index, 0));
#ifdef BEN_DEBUG
		assert(event_isAncestor(E, cap_getEvent(A)));
#endif
	} 

	// Step Forward
	if (cap_getSide(A)) {
		// If 5' cap
		*allowed_distance -= segment_getLength(cap_getSegment(A));
		if (*allowed_distance < 0)
			return false;
		A = cap_getOtherSegmentCap(A);
	} else {
		// If 3' cap
		if (cap_getEvent(A) != E) {
			bottomCapNumber = transmap_findBottomCapNumber(E,A);

			if (bottomCapNumber == 0)
				return false;
			else if (bottomCapNumber > 1) {
#ifdef BEN_DEBUG
				assert(transmap_findTopCap(E,A) == NULL);
#endif
				return false;
			} else {
				topCap = transmap_findTopCap(E,A);
				if (topCap && cap_getEvent(topCap) == E)
					A = topCap;
				if (topCap && cap_getTopFace(topCap)
				    && transmap_goByAncestralPath(cap_getTopFace(topCap), E)) 
					A = topCap;				
				else
					A = transmap_findBottomCap(E,A);
			}
		}

		// Recursion down into nested flower
		if (!transmap_sampleOrderAndOrientationAtNestedFlower(E, A, allowed_distance))
			return false;

		A = cap_getAdjacency(A);
		if (!cap_getOtherSegmentCap(A))
			return maybe;
	}

	// Recursion
	return transmap_sampleOrderAndOrientationAtEvent(E, A, B, allowed_distance);
}

/*
 * Wrapper to the transmap_sampleOrderAndOrientationAtEvent function
 */
static bool transmap_samplePathAtEvent(Event * E, Cap * A, Cap * B, int32_t distance) {
	switch (transmap_sampleOrderAndOrientationAtEvent(E, A, B, &distance)) {
	case true:
		return true;
	case false:
		return false;
	default:
		return (transmap_sampleOrderAndOrientationAtEvent(E, cap_getReverse(B), cap_getReverse(A), &distance) == maybe);
	}
}

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal thread containing
 * (A, B) at event E.
 *
 * All the *below* functions can be created by logic on this function.
 */
bool transmap_connectivityOrderAndOrientationWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier) {
	int index;
	int result = 0;
	int32_t distance = distance_multiplier * transmap_getTotalDistanceBetweenCaps(A,B);

	srand(time(NULL));

	if (distance == -1)
		return false;

	for (index = 0; index < sample_size; index++) 
		if (transmap_samplePathAtEvent(E, A, B, distance) 
		    && ++result > result_cutoff)
			return true;

	return false;
}

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal thread containing
 * (A, B) and/or (A, -B) (and their mirrors) at event E.
 */
bool transmap_connectivityAndOrientationWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier) {
	int index;
	int result = 0;
	int32_t distance = distance_multiplier * transmap_getTotalDistanceBetweenCaps(A,B);

	srand(time(NULL));

	if (distance == -1)
		return false;

	for (index = 0; index < sample_size; index++) 
		if (transmap_samplePathAtEvent(E, A, B, distance)
		    || transmap_samplePathAtEvent(E, B, A, distance))
			if (++result > result_cutoff)
				return true;

	return false;
}

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal thread containing
 * (A, B) and/or (B, A) (and their mirrors) at event E.
 */
bool transmap_connectivityAndOrderWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier) {
	int index;
	int result = 0;
	int32_t distance = distance_multiplier * transmap_getTotalDistanceBetweenCaps(A,B);

	srand(time(NULL));

	if (distance == -1)
		return false;

	for (index = 0; index < sample_size; index++) 
		if (transmap_samplePathAtEvent(E, A, B, distance)
		    || transmap_samplePathAtEvent(E, A, cap_getReverse(B), distance))
			if (++result > result_cutoff)
				return true;

	return false;

}

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal terminal thread containing
 * (A, B), (B, A), (A, -B), and/or (-B, A) (and their mirrors) at event E.
 */
bool transmap_connectivityWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier) {
	int index;
	int result = 0;
	int32_t distance = distance_multiplier * transmap_getTotalDistanceBetweenCaps(A,B);

	srand(time(NULL));

	if (distance == -1)
		return false;

	for (index = 0; index < sample_size; index++) 
		if (transmap_samplePathAtEvent(E, A, B, distance)
		    || transmap_samplePathAtEvent(E, A, cap_getReverse(B), distance)
		    || transmap_samplePathAtEvent(E, B, A, distance)
		    || transmap_samplePathAtEvent(E, cap_getReverse(B), A, distance))
			if (++result > result_cutoff)
				return true;

	return false;
}
