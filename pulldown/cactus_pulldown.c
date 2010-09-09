#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>

#include "cactus.h"
#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"

#include "cactus_pulldown.h"

/*
 * Take a flower, and determine if any pulldowns are necessary. 
 * Returns the Name of either a Segment, or a Cap w/ adjacency, 
 *  that needs a pulldown if such an object exists.
 * Returns NULL_NAME if no pulldowns are necessary.
 * */
Name pulldown_getFlowerPulldownName(Flower *flower) {
    return NULL_NAME;
}

// Calculate the weight/badness for a given segment in some Flower.
int32_t pulldown_getSegmentWeight(Segment *segment) {
    // Fail if given a NULL segment.
    if (segment == NULL) {
        fprintf(stderr, "pulldown_getSegmentWeight() passed NULL argument");
        exit(1);
    }
    int32_t wt = 0;

    // Orient the segment correctly and get the raw sequence (string).
    Sequence *sequence = segment_getSequence(segment_getPositiveOrientation(segment));
    char *rawSeq = sequence_getString(
            sequence,
            sequence_getStart(sequence),
            sequence_getLength(sequence), 0);

    Segment *child;
    Sequence *childSequence;
    char *childRawSeq;
    for (int32_t i = 0; i < segment_getChildNumber(segment); ++i) {
        child = segment_getChild(segment, i);

        // Get a possitively oriented sequence for child i.
        childSequence = segment_getSequence(segment_getPositiveOrientation(child));

        // Increase the weight and continue looping if the sequence lengths differ.
        if (sequence_getLength(sequence) != sequence_getLength(childSequence)) {
            ++wt; 
            continue;
        }

        // Get the raw sequence data for the oriented child sequence.
        childRawSeq = sequence_getString(
                childSequence,
                sequence_getStart(childSequence),
                sequence_getLength(childSequence), 0);

        // Increase the weight if the sequences are not the same.
        int unequalSequences = strncmp(rawSeq, childRawSeq, sequence_getLength(sequence));
        if (unequalSequences) {
            ++wt;
        }
        
        free(childRawSeq);
    }
    free(rawSeq);

    return wt;
}

int32_t pulldown_getSegmentBadness(Segment *segment) {
    int32_t wt = pulldown_getSegmentBadness(segment);
    assert(wt >= 0);

    // When the weight is not zero, the badness is one less.
    if (wt) {
        return wt - 1;
    }

    // The weight is zero otherwise.
    return 0;
}

void pulldown_getMostAncientAttachedDescendants_(Cap *cap, stList *acc) {
    Cap *child;
    for (int32_t i = 0; i < cap_getChildNumber(cap); i++) {
        child = cap_getChild(cap, i);
        if (cap_getAdjacency(child) != NULL) {
            stList_append(acc, child);
        }
        else { // The child is not attached. So dive deeper.
            pulldown_getMostAncientAttachedDescendants_(child, acc);
        }
    }
}

stList *pulldown_getMostAncientAttachedDescendants(Cap * cap) {
    stList *attachedDescendants = stList_construct();

    if (cap_getChildNumber(cap) > 0) 
        pulldown_getMostAncientAttachedDescendants_(cap, attachedDescendants);

    return attachedDescendants;
}

// Returns true when a descendant cap's adjacency lifts to an ancestor's
// adjacency.
int32_t pulldown_adjacencyLiftsConsistently(Cap *ancestor, Cap *descendant) {
    Cap *realAdjacency = cap_getAdjacency(ancestor);
    Cap *liftedAdjacency = cap_getAdjacency(descendant);

    assert(realAdjacency != NULL);
    assert(liftedAdjacency != NULL);

    // Climb up ancestry until reaching the real adjacency, or a source node.
    for (   ;
            liftedAdjacency != NULL;
            liftedAdjacency = cap_getParent(liftedAdjacency)) {

        // Stop climbing when we find the real adjacency.
        if (realAdjacency == liftedAdjacency)
            break;

        // Stop climing if we are otherwise done lifting the edge.
        if (cap_getAdjacency(liftedAdjacency) != NULL)
            break;
    }

    // Did we find the real adjacency?
    return realAdjacency == liftedAdjacency ? 1 : 0;
}

// Calculate the weight/badness for the adjacency of a Cap in some flower.
// 0 is returned if the given cap has no adjacency; we don't need to
// touch it.
int32_t pulldown_getCapAdjacencyWeight(Cap *cap) {
    int32_t wt = 0;
    Cap *adjacent = cap_getAdjacency(cap);

    if (adjacent == NULL) {
        return 0;
    }

    // Loop through each of the cap's descendants.
    stList *attachedDescendants = pulldown_getMostAncientAttachedDescendants(cap);
    //fprintf(stderr, "got attached %d descendants.\n", stList_length(attachedDescendants));
    Cap *descendant;
    stListIterator *it = stList_getIterator(attachedDescendants);
    for (   descendant = stList_getNext(it);
            descendant != NULL; descendant = stList_getNext(it)) {
        // Increase weight if the adjacency does not pull up correctly.
        if (!pulldown_adjacencyLiftsConsistently(cap, descendant)) 
            ++wt;
    }
    stList_destructIterator(it);
    stList_destruct(attachedDescendants);

    // Loop through each of the adjacent cap's descendants.
    stList *adjacentAttachedDescendants
        = pulldown_getMostAncientAttachedDescendants(adjacent);
    //fprintf(stderr, "got attached descendants of adjacent.\n");
    it = stList_getIterator(adjacentAttachedDescendants);
    for (   descendant = (Cap *)stList_getNext(it);
            descendant != NULL; descendant = (Cap *)stList_getNext(it)) {
        // Increase weight if the adjacency does not pull up correctly.
        if (!pulldown_adjacencyLiftsConsistently(adjacent, descendant)) 
            ++wt;
    }
    stList_destructIterator(it);
    stList_destruct(adjacentAttachedDescendants);
    it = NULL;

    //Cap *child;
    //for (int32_t i = 0; i < cap_getChildNumber(cap); ++i) {
    //    child = cap_getChild(cap,i);
    //}

    // Return the total number of inconsistent lifted adjacency edges.
    return wt;
}

int32_t pulldown_getCapAdjacencyBadness(Cap *cap) {
    int32_t wt = pulldown_getCapAdjacencyWeight(cap);
    assert(wt >= 0);

    // When the weight is not zero, the badness is one less.
    if (wt) {
        return wt - 1;
    }

    // The weight is zero otherwise.
    return 0;
}

/*
 * Take a flower, and perform pulldowns as necessary to insure junctions
 * are consistently attached.
 * */
void pulldown_doFlowerPulldowns(Flower * flower) {
}
