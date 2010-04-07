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
 * Check that a given direct adjacency (i.e. no nested net)
 * was preserved from the present back to the specified event (over_time == true), 
 * or at the specified time point (over_time == false)
 */
static int32_t transmap_directAdjacencyWasBroken(Cap * cap1, Cap * cap2, Event * event, int32_t over_time) {
	Event * event1 = cap_getEvent(cap1);
	Event * event2 = cap_getEvent(cap2);

	// TODO: This function assumes that the event DAG is a tree, other wise the whole metric 
	// makes no sense... 

	// If the chosen caps do not coincide in time
	if (event1 != event2) {
		// if cap2 is the daddy, switch caps
		if (event_isDescendant(event1, event2)) {
			Cap * tmp = cap2;
			cap2 = cap1;
			cap1 = tmp;
		}

		if (over_time 
		    || event2 == event
		    || event_isDescendant(event2, event)) {
			if (!transmap_isDescendant(cap_getAdjacency(cap1), cap2) 
			    || !cap_getAdjacency(cap2)
			    || cap_getParent(cap2)
			    || cap_getSegment(cap_getAdjacency(cap2)))
				return true;
			else 
				return transmap_directAdjacencyWasConserved(cap1, cap_getParent(cap2), event);	
		}
	}
	// Termination of recursion: if the process has reached an event anterior to the target event
	else if (event1 == event || event_isDescendant(event1, event))
		return cap_getAdjacency(cap1) == cap2; 
	// Termination of recursion: if the process hit a break
	else if (over_time && cap_getAdjacency(cap1) != cap2)
		return true;
	// Recursion: go up one generation
	else 
		return transmap_syntenyWasConserved(cap_getParent(cap1) , cap_getParent(cap2), event);
}

// Stub for recursive call in transmap_adjacencyWasBroken 
static int32_t transmap_syntenyWasBroken(Cap * cap1, Cap * cap2, Event * event, int32_t over_time);

/* 
 * Returns a boolean depending on whether there is an actual 
 * conserved nested path from cap to its adjacent cap
 */
static int32_t transmap_adjacencyWasBroken(Cap * cap, Event * event, int32_t over_time) {
	Cap * target = cap_getAdjacency(cap);
	Net * nestedNet = NULL;

#ifdef BEN_DEBUG
	assert(end_getGroup(cap_getEnd(cap)));
#endif

	nestedNet = group_getNestedNet(end_getGroup(cap_getEnd(cap)));

	// End of recursion: we hit rock bottom
	if (nestedNet == NULL)
		return transmap_directAdjacencyWasBroken(cap, target, event, over_time);

	// Project into nested net
	cap = net_getCap(nestedNet, cap_getName(cap));
	target = net_getCap(nestedNet, cap_getName(target));

#ifdef BEN_DEBUG
	assert(cap);
	assert(target);
#endif

	// Recursion: step into nested net
	return transmap_syntenyWasBroken(cap, target, event, over_time);
}


/*
 * Returns a boolean depending on whether there is a continuous
 * conserved path from cap1 to cap2
 */
static int32_t transmap_syntenyWasBroken(Cap * cap1, Cap * cap2, Event * event, int32_t over_time) {
	// Termination of recursion
	if (cap1 == cap2)
		return false;
	// If the caps have no official neighbours
	else if (!cap_getAdjacency(cap1) || !cap_getAdjacency(cap2))
		return true;
	// If the nested path is broken between cap1 and its immediate neighbour
	else if (transmap_adjacencyWasBroken(cap1, event, over_time)) 
		return true;
	// Recursion: test the synteny one block closer to cap2
	else 
		return transmap_syntenyWasBroken(cap_getReverse(cap_getAdjacency(cap1)), cap2, event, over_time);
	
}

/* 
 * Checks the breaks in a sequence of blocks from the present back to a 
 * specified event. 
 * Returns an array of booleans which is 1 element shorter than the array
 * of blocks (fence post principle)
 */
int32_t * transmap_syntenyWasConservedSinceEvent(Block ** blocks, int32_t blockCount, Event * event) {
	// TODO
}

/* 
 * Checks the breaks in a sequence of blocks at a given time point
 * Returns an array of booleans which is 1 element shorter than the array
 * of blocks (fence post principle)
 */
int32_t * transmap_syntenyWasConservedAtEvent(Block ** blocks, int32_t blockCount, Event * event) {
	// TODO
}
