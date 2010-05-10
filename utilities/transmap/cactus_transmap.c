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

// Sub for recursion
static int32_t transmap_getTotalDistanceBetweenCaps(Cap * A, Cap * B);

/*
 * Computes the current distance between two adjacent
 * cap in the nested sub-graph
 */
static int32_t transmap_getTotalDistanceAtAdjacency(Cap * A) {
	Net * net = group_getNestedNet(end_getGroup(cap_getEnd(A))); 
	Cap * B = cap_getAdjacency(A);
	Cap * childA, * childB;

	if (!net || !B)
		return 0;

	childA = net_getCap(net, cap_getName(A));
	childB = net_getCap(net, cap_getName(B));

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
	if (!A) 
		return false;
	B = cap_getAdjacency(A);
	if (!B)
		return false;

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
static bool transmap_goByAncestralPath(Face * face, Event * E) {
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
static bool transmap_sampleOrderAndOrientationAtEvent(Event * E, Cap * A, Cap * B, int32_t * allowed_distance) {
	Face * face;
	int32_t index;
	Cap * topCap, *bottomCap;
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
			if (face_getTopNode(face, index) == A)
				break;
#ifdef BEN_DEBUG
		assert(index < face_getCardinal(face));
		assert(face_getBottomNodeNumber(face, index) == 1);
#endif
		A = face_getBottomNode(face, index, 0);
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
				topCap = transmap_findTopCap(E,A);
				assert(topCap == NULL);
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

		// Recursion down into nested net
		if (!transmap_sampleOrderAndOrientationAtNestedNet(E, A, allowed_distance))
			return false;

		A = cap_getAdjacency(A);
	}

	// Recursion
	return transmap_sampleOrderAndOrientationAtEvent(E, A, B, allowed_distance);
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
	int32_t distance = transmap_getTotalDistanceBetweenCaps(A,B);
	int32_t allowed_distance;

	srand(time(NULL));

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

	srand(time(NULL));

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

	srand(time(NULL));

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

	srand(time(NULL));

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
