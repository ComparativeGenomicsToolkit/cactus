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

/* 
 * Computes the current distance between two caps
 * Returns -1 if no simple path was found
 */
static int32_t transmap_getTotalDistanceBetweenCaps(Cap * A, Cap * B) {
	int32_t length = 0;
	int32_t remaining_length;
	
	if (A == B)
		return 0;

	if (A == NULL)
		return -1;

	if (cap_getSide(A)) {
		A = cap_getOtherSegmentCap(A);
		length += segment_getLength(cap_getSegment(A));
		if (A == B)
			return length;
		if (A == NULL)
			return -1;
	}

	remaining_length = transmap_getTotalDistanceBetweenCaps(cap_getAdjacency(A), B);

	if (remaining_length == -1)
		return -1;
	else 
		return length + remaining_length;
}

// Stub for recursion
static bool transmap_sampleOrderAndOrientationAtEvent(Event * E, Cap * A, Cap * B, int32_t * allowed_distance);

/* 
 * Projects the stochastic search to a nested net
 */
static bool transmap_sampleOrderAndOrientationAtNestedNet(Event * E, Cap * A, int * allowed_distance) {
	Net * net = group_getNestedNet(end_getGroup(cap_getEnd(A))); 
	Cap * B, * childA, * childB;
	Event * childE;

	if (!net)
		return true;

	B = cap_getAdjacency(A);

#ifdef BEN_DEBUG
	assert(B);
#endif

	childA = net_getCap(net, cap_getName(A));
	childB = net_getCap(net, cap_getName(B));
	childE = eventTree_getEvent(net_getEventTree(net), event_getName(E));

#ifdef BEN_DEBUG
	assert(childA);
	assert(childB);
	assert(childE);
#endif

	return transmap_sampleOrderAndOrientationAtEvent(childE, childA, childB, allowed_distance);
}

/* 
 * Stochastic decision whether to jump a face or not
 */
static bool transmap_goByAncestralPath(Cap * A, Event * E) {
	Cap * parent = cap_getParent(A);
	Event * topEvent = cap_getEvent(parent);
	Event * bottomEvent;
	Event * tmpEvent;
	Face * face = cap_getFace(parent);
	float totalTimeSpan = 0;
	float partialTimeSpan = 0;
	float timeRatio;
	int32_t index;

	// Get the framing events of the face
	for (index = 0; index < face_getCardinal(face); index++) {
		tmpEvent = cap_getEvent(face_getTopNode(face, index));
		if (event_isAncestor(topEvent, tmpEvent))
			topEvent = tmpEvent;

#ifdef BEN_DEBUG
		assert(face_getBottomNodeNumber(face, index) == 1);
#endif
		tmpEvent = cap_getEvent(face_getBottomNode(face, index, 0));
		if (event_isAncestor(tmpEvent, bottomEvent))
			bottomEvent = tmpEvent;
	}

	// If event out of bounds	
	if (event_isAncestor(E, topEvent) || E == topEvent)
		return true;
	if (event_isAncestor(bottomEvent, E) || E == bottomEvent)
		return false;

	// Measure time span
	for (tmpEvent = bottomEvent; event_isAncestor(E, tmpEvent); tmpEvent = event_getParent(tmpEvent))
		partialTimeSpan += event_getBranchLength(tmpEvent);
	for (tmpEvent = bottomEvent; event_isAncestor(topEvent, tmpEvent); tmpEvent = event_getParent(tmpEvent))
		partialTimeSpan += event_getBranchLength(tmpEvent);

	if (partialTimeSpan == 0)
		return false;

	timeRatio = partialTimeSpan / totalTimeSpan;	

	// Random decision
	srand(time(NULL));
	return (rand() / RAND_MAX) > timeRatio;
}

/*
 * Tries to connect two caps at the time of event
 * Returns true if success
 */
static bool transmap_sampleOrderAndOrientationAtEvent(Event * E, Cap * A, Cap * B, int32_t * allowed_distance) {
	Face * face;
	int32_t index;

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
		assert(cap_getFace(A));
#endif
		face = cap_getFace(A);
		for (index = 0; index < face_getCardinal(face); index++)
			if (face_getTopNode(face, index) == A)
				break;
#ifdef BEN_DEBUG
		assert(index < face_getCardinal(face));
		assert(face_getBottomNodeNumber(face, index) == 1);
#endif
		A = face_getBottomNode(face, index, 0);
	} 

	if (event_isAncestor(E, cap_getEvent(A))) {
		while (cap_getParent(A)
		       && (event_isAncestor(E, cap_getEvent(cap_getParent(A))) || cap_getEvent(cap_getParent(A)) == E))
			A = cap_getParent(A);
	}

	// Recursion
	if (cap_getSide(A)) {
		// If 5' cap
		A = cap_getOtherSegmentCap(A);
		if (A == NULL)
			return false;
		*allowed_distance -= segment_getLength(cap_getSegment(A));
		if (*allowed_distance < 0)
			return false;
	} else {
		// If 3' cap
		if (cap_getEvent(A) != E
		    && cap_getFace(cap_getParent(A)) 
		    && transmap_goByAncestralPath(A, E))
			A = cap_getParent(A);				

		if (!cap_getAdjacency(A))
			return false;
		if (!transmap_sampleOrderAndOrientationAtNestedNet(E, A, allowed_distance))
			return false;
	}

	return transmap_sampleOrderAndOrientationAtEvent(E, A, B, allowed_distance);
}

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal thread containing
 * (A, B) at event E.
 *
 * All the below functions can be created by logic on this function.
 */
bool transmap_connectivityOrderAndOrientationWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier) {
	int index;
	int result = 0;
	int32_t distance = transmap_getTotalDistanceBetweenCaps(A,B);
	int32_t allowed_distance;

	if (distance == -1)
		return false;

	for (index = 0; index < sample_size; index++) {
		allowed_distance = distance * distance_multiplier;
		if (transmap_sampleOrderAndOrientationAtEvent(E, A, B, &distance))
			if (++result > result_cutoff)
				return true;
	}

	return false;
}

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal thread containing
 * (A, B) and/or (A, -B) (and their mirrors) at event E.
 */
bool transmap_connectivityAndOrientationWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier) {
	int index;
	int result = 0;
	int32_t distance = transmap_getTotalDistanceBetweenCaps(A,B);
	int32_t allowed_distance, allowed_distance2;

	if (distance == -1)
		return false;

	for (index = 0; index < sample_size; index++) {
		allowed_distance = distance * distance_multiplier;
		allowed_distance2 = distance * distance_multiplier;
		if (transmap_sampleOrderAndOrientationAtEvent(E, A, B, &distance)
		    || transmap_sampleOrderAndOrientationAtEvent(E, B, A, &distance))
			if (++result > result_cutoff)
				return true;
	}

	return false;
}

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal thread containing
 * (A, B) and/or (B, A) (and their mirrors) at event E.
 */
bool transmap_connectivityAndOrderWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier) {
	int index;
	int result = 0;
	int32_t distance = transmap_getTotalDistanceBetweenCaps(A,B);
	int32_t allowed_distance, allowed_distance2;

	if (distance == -1)
		return false;

	for (index = 0; index < sample_size; index++) {
		allowed_distance = distance * distance_multiplier;
		allowed_distance2 = distance * distance_multiplier;
		if (transmap_sampleOrderAndOrientationAtEvent(E, A, B, &distance)
		    || transmap_sampleOrderAndOrientationAtEvent(E, A, cap_getReverse(B), &distance))
			if (++result > result_cutoff)
				return true;
	}

	return false;

}

/*
 * Returns non-zero if and only if there exists a terminal thread or a virtual terminal terminal thread containing
 * (A, B), (B, A), (A, -B), and/or (-B, A) (and their mirrors) at event E.
 */
bool transmap_connectivityWasPresentAtEvent(Event *E, Cap *A, Cap *B, int sample_size, int result_cutoff, int distance_multiplier) {
	int index;
	int result = 0;
	int32_t distance = transmap_getTotalDistanceBetweenCaps(A,B);
	int32_t allowed_distance, allowed_distance2, allowed_distance3, allowed_distance4;

	if (distance == -1)
		return false;

	for (index = 0; index < sample_size; index++) {
		allowed_distance = distance * distance_multiplier;
		allowed_distance2 = distance * distance_multiplier;
		allowed_distance3 = distance * distance_multiplier;
		allowed_distance4 = distance * distance_multiplier;
		if (transmap_sampleOrderAndOrientationAtEvent(E, A, B, &distance)
		    || transmap_sampleOrderAndOrientationAtEvent(E, A, cap_getReverse(B), &distance)
		    || transmap_sampleOrderAndOrientationAtEvent(E, B, A, &distance)
		    || transmap_sampleOrderAndOrientationAtEvent(E, cap_getReverse(B), A, &distance))
			if (++result > result_cutoff)
				return true;
	}

	return false;
}
