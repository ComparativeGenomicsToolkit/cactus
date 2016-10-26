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
 * Globals definitions
 */
typedef struct adjacency_vote_st AdjacencyVote;
typedef struct adjacency_vote_table_st AdjacencyVoteTable;

/*
 * Global variables
 */
static int VOTE_CUTOFF = 16;

///////////////////////////////////////////////
// Adjacency vote :
//      message passing structure for ambiguous
//      partner configurations
///////////////////////////////////////////////
struct adjacency_vote_st {
    int length;
    Cap **candidates;
};

/*
 * Basic empty constructor
 */
static AdjacencyVote *adjacencyVote_construct(int64_t length) {
    AdjacencyVote *vote = st_calloc(1, sizeof(AdjacencyVote));

    vote->length = length;
    vote->candidates = st_calloc(length, sizeof(Cap *));

    return vote;
}

/*
 * Destructor
 */
static void adjacencyVote_destruct(AdjacencyVote * vote) {
    free(vote);
}

/*
 * Copy function
 */
static AdjacencyVote *adjacencyVote_copy(AdjacencyVote * vote) {
    AdjacencyVote *copy = adjacencyVote_construct(vote->length);
    int index;

    for (index = 0; index < vote->length; index++)
        copy->candidates[index] = vote->candidates[index];

    return copy;
}

static AdjacencyVote *adjacencyVoteTable_getVote(Cap * cap,
        AdjacencyVoteTable * table);
/*
 * Searches a Cap's ancestors in search of given node
 */
static int64_t adjacencyVote_isDescendantOf(Cap * descendant, Cap * ancestor,
        AdjacencyVoteTable * table) {
    Cap *current = descendant;

    if (current == ancestor)
        return true;

    if (cap_getParent(current) == ancestor)
        return true;

    while ((current = cap_getParent(current)) && adjacencyVoteTable_getVote(
            current, table)) {
        st_logInfo("Testing %p\n", current);
        if (cap_getParent(current) == ancestor)
            return true;
    }

    return false;
}

/*
 * Search for particular value
 */
static Cap *adjacencyVote_containsDescent(AdjacencyVote * vote, Cap * cap,
        AdjacencyVoteTable * table) {
    int index;
    End *end = cap_getEnd(cap);

    for (index = 0; index < vote->length; index++) {
        if (cap_getEnd(vote->candidates[index]) != end)
            continue;
        else if (adjacencyVote_isDescendantOf(vote->candidates[index],
                cap_getPositiveOrientation(cap), table))
            return vote->candidates[index];
    }

    return NULL;
}

/*
 * Return number of candidates
 */
static int adjacencyVote_getLength(AdjacencyVote * vote) {
    return vote->length;
}

/*
 * Return first value
 */
static Cap *adjacencyVote_getWinner(AdjacencyVote * vote) {
    return vote->candidates[0];
}

/*
 * Checks if vote points to NULL
 */
static bool adjacencyVote_isBlank(AdjacencyVote * vote) {
    return vote->length == 1 && vote->candidates[0] == NULL;
}

/*
 * Determine what adjacencies are feasible for a given node depending on the adjacencies of children
 */
static AdjacencyVote *adjacencyVote_processVotes(AdjacencyVote * vote_1,
        AdjacencyVote * vote_2) {
    AdjacencyVote *merged_vote;
    int index_merged = 0;
    int index_1 = 0;
    int index_2 = 0;
    bool blank_1 = adjacencyVote_isBlank(vote_1);
    bool blank_2 = adjacencyVote_isBlank(vote_2);

    if (blank_1 && blank_2)
        return adjacencyVote_copy(vote_1);
    else if (blank_1)
        return adjacencyVote_copy(vote_2);
    else if (blank_2)
        return adjacencyVote_copy(vote_1);

    merged_vote = adjacencyVote_construct(vote_1->length + vote_2->length);

    // Computing intersection first
    while (index_1 < vote_1->length && index_2 < vote_2->length) {
        if (vote_1->candidates[index_1] == vote_2->candidates[index_2]) {
            merged_vote->candidates[index_merged++]
                    = vote_1->candidates[index_1++];
            index_2++;
        } else if (cap_getName(vote_1->candidates[index_1]) < cap_getName(
                vote_2->candidates[index_2]))
            index_1++;
        else if (cap_getName(vote_1->candidates[index_1]) > cap_getName(
                vote_2->candidates[index_2]))
            index_2++;
    }

    if (index_merged > 0) {
        merged_vote->length = index_merged;
        return merged_vote;
    }

    index_merged = 0;
    index_1 = 0;
    index_2 = 0;

    // If intersection failed: union
    while (index_1 < vote_1->length || index_2 < vote_2->length) {
        if (index_1 < vote_1->length && index_2 < vote_2->length
                && vote_1->candidates[index_1] == vote_2->candidates[index_2]) {
            merged_vote->candidates[index_merged++]
                    = vote_1->candidates[index_1++];
            index_2++;
        } else if (index_2 == vote_2->length)
            merged_vote->candidates[index_merged++]
                    = vote_1->candidates[index_1++];
        else if (index_1 == vote_1->length)
            merged_vote->candidates[index_merged++]
                    = vote_2->candidates[index_2++];
        else if (cap_getName(vote_1->candidates[index_1]) < cap_getName(
                vote_2->candidates[index_2]))
            merged_vote->candidates[index_merged++]
                    = vote_1->candidates[index_1++];
        else if (cap_getName(vote_1->candidates[index_1]) > cap_getName(
                vote_2->candidates[index_2]))
            merged_vote->candidates[index_merged++]
                    = vote_2->candidates[index_2++];
    }

    merged_vote->length = index_merged;
    return merged_vote;
}

/*
 * Obtains the end instance which is contemporary to given event, or closest after event
 */
static Cap *adjacencyVote_getWitnessAncestor(Cap * cap, Event * event) {
    Cap *current = cap;
    Cap *parent;

    // TODO:
    // If not event DAG, then have to choose parent
    // Tomorrow's an other day...
    while ((parent = cap_getParent(current)) && event_isDescendant(event,
            cap_getEvent(parent)))
        current = parent;

    if (parent && cap_getEvent(parent) == event)
        return parent;

    return current;
}

/*
 * Replaces the candidates of a vote with their parents
 */
static void adjacencyVote_goToParents(AdjacencyVote * vote, Cap * generation) {
    int index;
    Event *event = cap_getEvent(generation);

    for (index = 0; index < vote->length; index++)
        vote->candidates[index]
                = cap_getPositiveOrientation(adjacencyVote_getWitnessAncestor(
                        vote->candidates[index], event));

}

///////////////////////////////////////////////
// Adjacency vote table :
//      simple lookup table of AdjacencyVotes
///////////////////////////////////////////////

struct adjacency_vote_table_st {
    struct hashtable *table;
    struct List * computationFront;
};

static uint64_t hash_from_key_fn(const void *key) {
    return (uint64_t) cap_getName((Cap *) key);
}

static int keys_equal_fn(const void *key1, const void *key2) {
    return key1 == key2;
}

static void valueFree(void *value) {
    adjacencyVote_destruct((AdjacencyVote *) value);
}

// Basic constructor
static AdjacencyVoteTable *adjacencyVoteTable_construct() {
    AdjacencyVoteTable *table = st_calloc(1, sizeof(AdjacencyVoteTable));
    table->table = create_hashtable(16, hash_from_key_fn, keys_equal_fn, NULL,
            valueFree);
    table->computationFront = constructZeroLengthList(1000, NULL);

    return table;
}

// Basic destructor
static void adjacencyVoteTable_destruct(AdjacencyVoteTable * table) {
    hashtable_destroy(table->table, true, false);
    destructList(table->computationFront);
    free(table);
}

// Get votes for Cap
static AdjacencyVote *adjacencyVoteTable_getVote(Cap * cap,
        AdjacencyVoteTable * table) {
    if (cap)
        return (AdjacencyVote *) hashtable_search(table->table,
                (void *) cap_getPositiveOrientation(cap));
    else
        return NULL;
}

/*
 * Insert vote for Cap
 */
static void adjacencyVoteTable_recordVote(AdjacencyVoteTable * table,
        Cap * cap, AdjacencyVote * vote) {
    AdjacencyVote * old_vote;

    if (!cap_getOrientation(cap))
        cap = cap_getReverse(cap);
    old_vote = hashtable_remove(table->table, (void *) cap, false);

    assert(vote->length == 0 || vote->candidates);
    if (listContains(table->computationFront, cap))
        listRemove(table->computationFront, cap);

    adjacencyVote_destruct(old_vote);
    hashtable_insert(table->table, (void *) cap, (void *) vote);
}

/*
 * Check if node actually votes for something
 */
static bool adjacencyVoteTable_doesNotVote(Cap * cap,
        AdjacencyVoteTable * table) {
    AdjacencyVote * vote = adjacencyVoteTable_getVote(cap, table);

    return (vote == NULL || vote->length == 0);
}

///////////////////////////////////////////////
// Filling in method
///////////////////////////////////////////////

/*
 * If possible register parent node into compuation front
 */
static void fillingIn_registerParent(Cap * cap, AdjacencyVoteTable * table) {
    Cap *parent = cap_getParent(cap);
    int64_t childIndex, childNumber;

    st_logInfo("Registering parent node %p\n", cap);

    if (parent == NULL)
        return;
    if (cap_getAdjacency(parent)) {
        fillingIn_registerParent(parent, table);
        return;
    }

    assert(cap_getAdjacency(parent) == NULL);

    childNumber = cap_getChildNumber(parent);

    for (childIndex = 0; childIndex < childNumber; childIndex++)
        if (adjacencyVoteTable_doesNotVote(cap_getChild(parent, childIndex),
                table))
            return;

    st_logInfo("New node in coputation front: %p\n", parent);
    listAppend(table->computationFront, cap_getPositiveOrientation(parent));
}

/*
 * Initiate propagation by counting bottom nodes as decided
 */
static void fillingIn_registerLeafCap(Cap * cap, AdjacencyVoteTable * table) {
    AdjacencyVote *vote;

    st_logInfo("Visiting leaf %p\n", cap);

    // Mark as decided
    vote = adjacencyVote_construct(1);
    vote->candidates[0] = cap_getPositiveOrientation(cap_getAdjacency(cap));
    adjacencyVoteTable_recordVote(table, cap, vote);

    // Register parent
    fillingIn_registerParent(cap, table);
}

static void fillingIn_propagateAdjacencyDownwards(Cap * cap, Cap * partner,
        AdjacencyVoteTable * table);

/*
 * Propagate adjacency to all children
 */
static void fillingIn_propagateToChildren(Cap * cap, Cap * partner,
        AdjacencyVoteTable * table) {
    int64_t childIndex;
    int64_t childNumber = cap_getChildNumber(cap);
    Cap ** array = calloc(cap_getChildNumber(cap), sizeof(Cap*));

    // This playing around is necessary to avoid list corruption between steps
    for (childIndex = 0; childIndex < childNumber; childIndex++)
        array[childIndex] = cap_getChild(cap, childIndex);

    // Propagate to children if necessary
    for (childIndex = 0; childIndex < childNumber; childIndex++) {
        st_logInfo("Child %" PRIi64 ": %p\n", childIndex, array[childIndex]);
        fillingIn_propagateAdjacencyDownwards(array[childIndex], partner, table);
    }

    free(array);
}

/*
 * Register decision to remain without adjacency
 */
static void fillingIn_giveUp(Cap * cap, AdjacencyVoteTable * table) {
    AdjacencyVote *blank_vote;
    Cap * partner;

    if (cap == NULL || !adjacencyVoteTable_getVote(cap, table))
        return;

    // If forcing surrender through ancestral links
    if ((partner = cap_getAdjacency(cap))) {
        partner = cap_getPositiveOrientation(partner);
        cap_breakAdjacency(cap);
        blank_vote = adjacencyVote_construct(0);
        adjacencyVoteTable_recordVote(table, partner, blank_vote);
        fillingIn_giveUp(partner, table);
    }

    // Create blank vote
    blank_vote = adjacencyVote_construct(0);
    adjacencyVoteTable_recordVote(table, cap, blank_vote);
}

static void
        fillingIn_processChildrenVote(Cap * cap, AdjacencyVoteTable * table);

/*
 * When joining a node to a partner, ensure that partner's parents are aware of this
 */
static void fillingIn_propagateToParents(Cap * child,
        AdjacencyVoteTable * table) {
    Cap *parent = cap_getParent(child);

    assert(child != NULL);

    // Already at root
    if (parent == NULL)
        return;

    // Parent already happily paired up
    if (cap_getAdjacency(parent))
        return;

    // TODO: how about converting an old abstentionist?
    if (adjacencyVoteTable_doesNotVote(parent, table))
        return;

    // Update parent's vote
    fillingIn_processChildrenVote(parent, table);
}

/*
 * Create an interpolation between parent and child caps at event
 */
static Cap *fillingIn_interpolateCaps(Cap * parentCap, Cap * childCap,
        Event * event) {
    Cap * newCap = cap_construct(cap_getEnd(childCap), event);

    if (!cap_getOrientation(newCap))
        newCap = cap_getReverse(newCap);

    st_logInfo("Interpolating %p\n", newCap);
    assert(event_isDescendant(cap_getEvent(parentCap), cap_getEvent(newCap)));
    assert(event_isDescendant(cap_getEvent(newCap), cap_getEvent(childCap)));
    assert(event_isDescendant(cap_getEvent(parentCap), cap_getEvent(childCap)));
    cap_changeParentAndChild(newCap, childCap);
    cap_makeParentAndChild(parentCap, newCap);

    return newCap;
}

/*
 * Make two end instances adjacent, and record the decision
 */
static void fillingIn_uniteCaps(Cap * cap, Cap * partner, AdjacencyVote * vote,
        AdjacencyVoteTable * table) {
    AdjacencyVote * partnerVote;

    if (cap_getEvent(cap) != cap_getEvent(partner))
        partner = fillingIn_interpolateCaps(cap_getParent(partner), partner,
                cap_getEvent(cap));

    cap_makeAdjacent(cap, partner);

    partnerVote = adjacencyVote_construct(1);
    partnerVote->candidates[0] = cap_getPositiveOrientation(cap);
    adjacencyVoteTable_recordVote(table, partner, partnerVote);

    // End instance's vote
    vote = adjacencyVote_construct(1);
    vote->candidates[0] = cap_getPositiveOrientation(partner);
    adjacencyVoteTable_recordVote(table, cap, vote);

    // Tell all the family
    st_logInfo("UNITE CAPS: PARNTERS CHILDREN\n");
    fillingIn_propagateToChildren(partner, cap, table);
    st_logInfo("UNITE CAPS: PARTNERS PARENTS\n");
    fillingIn_propagateToParents(partner, table);
    st_logInfo("UNITE CAPS: CHILDREN of %p\n", cap);
    fillingIn_propagateToChildren(cap, partner, table);
}

/*
 * Tests if a partner agrees to be paired up
 */
static bool fillingIn_partnerConsents(Cap * partner, Cap * cap,
        AdjacencyVoteTable * table) {
    AdjacencyVote *partner_vote = adjacencyVoteTable_getVote(partner, table);
    Cap * partnerPartner = cap_getAdjacency(partner);

    if (!partner_vote)
        return false;

    if (partner == cap)
        return false;

    if (partnerPartner)
        partnerPartner = cap_getPositiveOrientation(partnerPartner);

    if (partner == cap || adjacencyVote_isDescendantOf(partner, cap, table))
        return false;

    if (partnerPartner && partnerPartner != cap
            && !adjacencyVote_isDescendantOf(partnerPartner, cap, table))
        return false;

    return adjacencyVote_containsDescent(partner_vote, cap, table);
}

/*
 * Creates an inteprolation half way on the branch between two events
 */
static Event *fillingIn_interpolateEvents(Event* parentEvent, Event* childEvent) {
    float branchLength = 0.0;
    EventTree * eventTree = event_getEventTree(parentEvent);
    //Flower * flower = eventTree_getFlower(eventTree);
    Event * ptr = childEvent;
    Event * result;

    assert(event_isDescendant(parentEvent, childEvent));

    // Compute total branch length
    for (ptr = childEvent; ptr != parentEvent; ptr = event_getParent(ptr))
        branchLength += event_getBranchLength(ptr);

    if (branchLength <= 1) {
        return event_construct4("interpolation", 0, event_getParent(childEvent),
                childEvent, eventTree);
    }

    // Compute desired branch length _from the bottom_
    branchLength /= 2;

    for (ptr = childEvent; ptr != parentEvent; ptr = event_getParent(ptr)) {
        if (event_getBranchLength(ptr) == branchLength) {
            assert(event_isDescendant(parentEvent, event_getParent(ptr)));
            return event_getParent(ptr);
        }
        if (event_getBranchLength(ptr) > branchLength)
            break;

        branchLength -= event_getBranchLength(ptr);
        assert(ptr);
    }

    result = event_construct4("interpolation", event_getBranchLength(ptr)
            - branchLength, event_getParent(ptr), ptr, eventTree);
    assert(event_getBranchLength(result) >= 0);
    assert(event_getBranchLength(ptr) >= 0);
    assert(event_getBranchLength(event_getParent(ptr)) >= 0);
    assert(event_isDescendant(parentEvent, ptr));
    assert(event_isDescendant(parentEvent, result));
    return result;
}

/*
 * Builds a free stub end with the given group. The new end contains a cap with the given
 * event, and if the given event is not the root event, a root cap to attach the new
 * cap to the event.
 */
static Cap * fillingIn_constructStub(Event *event, Group *group) {
    EventTree *eventTree = event_getEventTree(event);
    End *newFreeStubEnd = end_construct(0, group_getFlower(group));
    end_setGroup(newFreeStubEnd, group);
    Cap *cap = cap_construct(newFreeStubEnd, event);
    st_logInfo("Constructing %p\n", cap);
    Event *rootEvent = eventTree_getRootEvent(eventTree);
    if (event == rootEvent) {
        end_setRootInstance(newFreeStubEnd, cap);
    } else {
        end_setRootInstance(newFreeStubEnd, cap_construct(newFreeStubEnd,
                rootEvent));
        st_logInfo("Constructing %p\n", newFreeStubEnd);
        cap_makeParentAndChild(end_getRootInstance(newFreeStubEnd), cap);
    }
    return cap;
}

/*
 * Pair up a node with a null stub if all else fails
 */
static void fillingIn_pairUpToNullStub(Cap * cap, AdjacencyVoteTable * table) {
    Cap * stub = fillingIn_constructStub(cap_getEvent(cap), end_getGroup(
            cap_getEnd(cap)));
    AdjacencyVote * vote, *stubVote;

    cap_makeAdjacent(cap, stub);

    if (table) {
        // Create new votes bulletins to stuff the electoral box
        vote = adjacencyVote_construct(1);
        vote->candidates[0] = cap_getPositiveOrientation(stub);
        adjacencyVoteTable_recordVote(table, cap, vote);

        stubVote = adjacencyVote_construct(1);
        stubVote->candidates[0] = cap_getPositiveOrientation(cap);
        adjacencyVoteTable_recordVote(table, stub, stubVote);

        // Send info down
        fillingIn_propagateToChildren(cap, stub, table);
    }
}

/*
 * Register decision to pick a candidate at random
 */
static void fillingIn_chooseAtRandom(Cap * cap, AdjacencyVote * vote,
        AdjacencyVoteTable * table) {
    Cap *partner;
    int64_t index;

    for (index = 0; index < vote->length; index++) {
        partner = vote->candidates[index];
        partner = adjacencyVote_getWitnessAncestor(partner, cap_getEvent(cap));
        if (fillingIn_partnerConsents(partner, cap, table)) {
            st_logInfo("Found consenting partner %p to %p\n", partner, cap);
            fillingIn_uniteCaps(cap, partner, vote, table);
            return;
        }
    }
}

/*
 * Determine adjacencies from the module pointed by the given Cap
 */
static void fillingIn_processChildrenVote(Cap * cap, AdjacencyVoteTable * table) {
    int64_t childrenNumber = cap_getChildNumber(cap);
    Cap *first_child, *second_child;
    AdjacencyVote *first_child_vote, *second_child_vote, *merged_vote = NULL;
    AdjacencyVote *partner_vote = NULL;
    Cap *partner = NULL;

    //Already paired up
    if (cap_getAdjacency(cap))
        return;

    // Consult children
    switch (childrenNumber) {
        case 0:
            // NOTE: This should not happen in theory since only
            // children can register their parents
            assert(false);
            break;
        case 1:
            first_child = cap_getChild(cap, 0);
            merged_vote = adjacencyVote_copy(adjacencyVoteTable_getVote(
                    first_child, table));
            adjacencyVote_goToParents(merged_vote, cap);
            break;
        case 2:
            first_child = cap_getChild(cap, 0);
            second_child = cap_getChild(cap, 1);

            first_child_vote = adjacencyVote_copy(adjacencyVoteTable_getVote(
                    first_child, table));
            second_child_vote = adjacencyVote_copy(adjacencyVoteTable_getVote(
                    second_child, table));
            adjacencyVote_goToParents(first_child_vote, cap);
            adjacencyVote_goToParents(second_child_vote, cap);
            merged_vote = adjacencyVote_processVotes(first_child_vote,
                    second_child_vote);
            adjacencyVote_destruct(first_child_vote);
            adjacencyVote_destruct(second_child_vote);
            break;
    }

    // Make decision
    if (merged_vote->length == 1) {
        partner = adjacencyVote_getWinner(merged_vote);

        if (partner == NULL || !(partner_vote = adjacencyVoteTable_getVote(
                partner, table))) {
            // If all children unattached or
            // Partner undecided... quizas, quizas, quizas
            st_logInfo("1\n");
            fillingIn_propagateToParents(cap, table);
            adjacencyVoteTable_recordVote(table, cap, merged_vote);
        } else if (fillingIn_partnerConsents(partner, cap, table)) {
            // Partner agrees!
            st_logInfo("2\n");
            fillingIn_uniteCaps(cap, partner, merged_vote, table);
        } else {
            // Partner has someone else in their life...
            st_logInfo("3\n");
            fillingIn_pairUpToNullStub(cap, table);
            adjacencyVote_destruct(merged_vote);
        }
    }

    // If still undecided record alternatives
    else if (adjacencyVote_getLength(merged_vote) < VOTE_CUTOFF) {
        st_logInfo("4\n");
        fillingIn_propagateToParents(cap, table);
        adjacencyVoteTable_recordVote(table, cap, merged_vote);
    }

    // If necessary to make a decision because junction node
    else if (cap_getChildNumber(cap) > 1) {
        st_logInfo("5\n");
        fillingIn_chooseAtRandom(cap, merged_vote, table);
        // If failure on junction: force to null stub
        if (!cap_getAdjacency(cap)) {
            partner = adjacencyVote_getWinner(merged_vote);
            fillingIn_pairUpToNullStub(cap, table);
        }
        adjacencyVote_destruct(merged_vote);
    }

    // If too lazy to make a decision
    else {
        st_logInfo("6\n");
        fillingIn_giveUp(cap, table);
        adjacencyVote_destruct(merged_vote);
    }
}

/*
 * Recursively apply second Fitch rule downwards when an ambiguity is resolved
 */
static void fillingIn_propagateAdjacencyDownwards(Cap * cap, Cap * partner,
        AdjacencyVoteTable * table) {
    Cap *decision;
    AdjacencyVote *vote = adjacencyVoteTable_getVote(cap, table);

    assert(vote->length >= 1 || cap_getAdjacency(cap));
    st_logInfo("Going downwards to %p\n", cap);

    // If Cap does not exist (for NULL partners) or child already decided
    if (!cap || vote->length == 0 || cap_getAdjacency(cap)) {
        st_logInfo("Vote? %" PRIi64 "\n", vote->length);
        if (vote->length)
            st_logInfo("destination %p\n", vote->candidates[0]);
        return;
    }

    assert(partner != NULL);

    decision = adjacencyVote_containsDescent(vote, partner, table);
    if (decision)
        decision
                = adjacencyVote_getWitnessAncestor(decision, cap_getEvent(cap));

    if (!decision || !fillingIn_partnerConsents(decision, cap, table)) {
        // If partner does not exist or refuses
        st_logInfo("A\n");
        fillingIn_chooseAtRandom(cap, vote, table);
    } else {
        // Determined partner accepts
        st_logInfo("B\n");
        fillingIn_uniteCaps(cap, decision, vote, table);
    }

    // If failure on junction: force to null stub
    if (!cap_getAdjacency(cap)) {
        if (cap_getChildNumber(cap) > 1) {
            if (decision && adjacencyVoteTable_getVote(decision, table)) {
                st_logInfo("C\n");
                fillingIn_pairUpToNullStub(cap, table);
            } else {
                st_logInfo("D\n");
                adjacencyVote_goToParents(vote, cap);
                decision = adjacencyVote_getWinner(vote);
                fillingIn_pairUpToNullStub(cap, table);
            }
        }
    }
}

/*
 * Proagate along compution front
 */
static void fillingIn_stepForward(Cap * cap, AdjacencyVoteTable * table) {
    st_logInfo("Stepping into %p\n", cap);
    // Assess node decision
    fillingIn_processChildrenVote(cap, table);

    // Add parents to computationFront
    fillingIn_registerParent(cap, table);
}

/*
 * Force the top junction node in tree to come to a decision
 */
static void fillingIn_forceIndecisiveTree(Cap * cap, AdjacencyVoteTable * table) {
    AdjacencyVote * vote;
    //Cap * partner;

    st_logInfo("Starting at tree top %p\n", cap);

    while (cap_getChildNumber(cap) == 1)
        cap = cap_getChild(cap, 0);

    // If non branched tree
    if (cap_getChildNumber(cap) == 0)
        return;

    // If attached, job done
    if (cap_getAdjacency(cap))
        return;

    // If not done, sort it out!!
    vote = adjacencyVoteTable_getVote(cap, table);
    st_logInfo("AA\n");
    fillingIn_chooseAtRandom(cap, vote, table);
    // If failure on junction: force to null stub
    if (!cap_getAdjacency(cap)) {
        st_logInfo("BB\n");
        adjacencyVote_goToParents(vote, cap);
        adjacencyVote_getWinner(vote);
        fillingIn_pairUpToNullStub(cap, table);
    }
}

/*
 * Go up to the first registered and attached ancestor (if exists)
 */
static Cap *fillingIn_getAttachedAncestor(Cap * cap) {
    Cap *current = cap_getParent(cap);
    Cap *parent;

    while (current && !cap_getAdjacency(cap) && (parent
            = cap_getParent(current)))
        current = parent;

    return current;
}

/*
 * Tests if two event are comparable
 */
static bool fillingIn_eventsAreComparable(Event * A, Event * B) {
    return A == B || event_isAncestor(A, B) || event_isAncestor(B, A);
}

/*
 * Returns the more ancetral of two comparable events
 */
static Event * fillingIn_oldestEvent(Event * A, Event * B) {
#ifndef NDEBUG
    assert(fillingIn_eventsAreComparable(A,B));
#endif
    return (A == B || event_isDescendant(A, B)) ? A : B;
}

/*
 * Stub for recursion
 */
static void fillingIn_testForPullDown(Cap * cap);

/*
 * Amend graph to remove inconsistent adjacencies 
 */
static void fillingIn_pullDown(Cap * A) {
    Cap * B = cap_getAdjacency(A);
    Cap * childA = NULL, *childB =NULL;
    Cap * interpolA = NULL, *interpolB = NULL;
    Event * interpolEvent = NULL;
    int64_t indexA, indexB;
    Event * eventChildA = NULL, *eventChildB = NULL, *childEvent = NULL;
    bool result = false;

    assert(B);
    assert(cap_getChildNumber(A));
    assert(cap_getEvent(A) == cap_getEvent(B));

    // If B has no children (presumably a stub)
    if (!cap_getChildNumber(B)) {
        childA = cap_getChild(A, 0);
        eventChildA = cap_getEvent(childA);
        childEvent = eventChildA;
        interpolEvent
                = fillingIn_interpolateEvents(cap_getEvent(A), childEvent);
        interpolA = fillingIn_interpolateCaps(A, childA, interpolEvent);
        interpolB = cap_construct(cap_getEnd(B), interpolEvent);
        cap_makeParentAndChild(B, interpolB);
        cap_makeAdjacent(interpolA, interpolB);
        puts("C");
        return;
    }

    // Find two children which are comparable:
    for (indexA = 0; indexA < cap_getChildNumber(A); indexA++) {
        childA = cap_getChild(A, indexA);
        eventChildA = cap_getEvent(childA);
        for (indexB = 0; indexB < cap_getChildNumber(B); indexB++) {
            childB = cap_getChild(B, indexB);
            eventChildB = cap_getEvent(childB);
            if ((result = fillingIn_eventsAreComparable(eventChildA,
                    eventChildB)))
                break;
        }

        if (result)
            break;
    }

    if (!result) {
        // Ouch, weird speciation event, got to clean up!!
        fillingIn_pairUpToNullStub(A, NULL);
        fillingIn_pairUpToNullStub(B, NULL);
        fillingIn_testForPullDown(fillingIn_getAttachedAncestor(A));
        fillingIn_testForPullDown(fillingIn_getAttachedAncestor(B));
        puts("A");
        return;
    }

    childEvent = fillingIn_oldestEvent(eventChildA, eventChildB);
    interpolEvent = fillingIn_interpolateEvents(cap_getEvent(A), childEvent);

    interpolA = fillingIn_interpolateCaps(A, childA, interpolEvent);
    interpolB = fillingIn_interpolateCaps(B, childB, interpolEvent);

    cap_makeAdjacent(interpolA, interpolB);
    puts("B");
}

/*
 * Tests if a node and its adjacency should be pulled down:
 */
static void fillingIn_testForPullDown(Cap * cap) {
    Cap *child;
    Cap *partner, *adjacencyAncestor;
    int64_t childIndex;
    int64_t inconsistencyCount = 0;

    assert(cap_getAdjacency(cap));

    // Go through children
    for (childIndex = 0; childIndex < cap_getChildNumber(cap); childIndex++) {
        child = cap_getChild(cap, childIndex);

        // Go down to attached child (if any)
        while (child && !cap_getAdjacency(child)) {
            assert(cap_getChildNumber(child) < 2);
            if (cap_getChildNumber(child))
                child = cap_getChild(child, 0);
            else
                child = NULL;
        }

        // Fringe branch => no worries!
        if (child == NULL)
            break;

        partner = cap_getAdjacency(child);
        partner = cap_getPositiveOrientation(partner);

        // ... lift
        adjacencyAncestor = fillingIn_getAttachedAncestor(partner);

        // Self loop
        if (adjacencyAncestor == cap) {
            fillingIn_pullDown(cap);
            break;
        }

        // Derived branch:
        if (adjacencyAncestor != cap_getAdjacency(cap) && ++inconsistencyCount
                > 1) {
            fillingIn_pullDown(cap);
            break;
        }
    }
}

/*
 * Remove all lifted self loops and inconsistent adjacencies from flower
 */
static void fillingIn_pullDowns(Flower * flower) {
    Flower_CapIterator *iter = flower_getCapIterator(flower);
    Cap *cap;

    // Iterate through potential top nodes
    while ((cap = flower_getNextCap(iter)))
        if (cap_getChildNumber(cap) > 1)
            fillingIn_testForPullDown(cap);

    flower_destructCapIterator(iter);
}

/*
 * Fill adjacencies in flower using descent trees and event tree
 * Correct non-AVG anomalies
 */
void fillingIn_fillAdjacencies(Flower * flower) {
    Cap *cap;
    AdjacencyVoteTable *table = adjacencyVoteTable_construct();
    Flower_CapIterator *iter = flower_getCapIterator(flower);

    //////////////////////////////////////////////////////////////
    // Define a computation front
    //////////////////////////////////////////////////////////////
    st_logInfo("Registering leaf caps\n");
    while ((cap = flower_getNextCap(iter))) {
        st_logInfo("Testing possible leaf %p\n", cap);
        if (cap_getChildNumber(cap) == 0)
            fillingIn_registerLeafCap(cap, table);
    }
    flower_destructCapIterator(iter);

    //////////////////////////////////////////////////////////////
    // Compute greedily
    //////////////////////////////////////////////////////////////
    st_logInfo("Propagation\n");
    while (table->computationFront->length > 0)
        fillingIn_stepForward(listRemoveFirst(table->computationFront), table);

    //////////////////////////////////////////////////////////////
    // Force decision for higher nodes
    //////////////////////////////////////////////////////////////
    st_logInfo("Force top of tree decisions\n");
    iter = flower_getCapIterator(flower);
    while ((cap = flower_getNextCap(iter)))
        if (!cap_getParent(cap))
            fillingIn_forceIndecisiveTree(cap, table);
    flower_destructCapIterator(iter);

    //////////////////////////////////////////////////////////////
    // Remove self loops
    //////////////////////////////////////////////////////////////
    st_logInfo("Pull downs\n");
    fillingIn_pullDowns(flower);

    //////////////////////////////////////////////////////////////
    // Clean up
    //////////////////////////////////////////////////////////////
    st_logInfo("Clean up\n");
    adjacencyVoteTable_destruct(table);
}

/*
 * Prints out information
 */
static void usage() {
    fprintf(
            stderr,
            "cactus_fillAdjacencies [flower-names, ordered by order they should be processed], version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-e --tempDirRoot : The temp file root directory\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

/*
 * Command line wrapper
 */
int main(int argc, char ** argv) {
    /*
     * The script builds trees.
     *
     * (1) It reads the flower.
     *
     * (2) It fill adjacencies
     *
     * (3) It copies the relevant augmented events into the event trees of its descendants.
     *
     */
    CactusDisk *cactusDisk;
    Flower *flower;
    int64_t startTime;
    int64_t j;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' }, { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    //assert(logLevelString == NULL || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
    assert(cactusDiskDatabaseString != NULL);
    //assert(tempFileRootDirectory != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    //////////////////////////////////////////////
    //For each flower do tree building..
    //////////////////////////////////////////////

    for (j = optind; j < argc; j++) {
        const char *flowerName = argv[j];
        st_logInfo("Processing the flower named: %s\n", flowerName);

        ///////////////////////////////////////////////////////////////////////////
        // Parse the basic reconstruction problem
        ///////////////////////////////////////////////////////////////////////////

        flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(
                flowerName));
        assert(flower_builtTrees(flower)); //we must have already built the trees for the problem at this stage
        st_logInfo("Parsed the flower to be refined\n");

        ///////////////////////////////////////////////////////////////////////////
        // Do nothing if we have already built the faces.
        ///////////////////////////////////////////////////////////////////////////

        if (flower_builtFaces(flower)) {
            st_logInfo("We have already built faces for flower %s\n",
                    flowerName);
            continue;
        }

        ///////////////////////////////////////////////////////////////////////////
        // Do nothing if not a leaf.
        ///////////////////////////////////////////////////////////////////////////

        if (!flower_isLeaf(flower)) {
            st_logInfo(
                    "We currently only build flowers for terminal problems: %s\n",
                    flowerName);
            continue;
        }
        assert(flower_isTerminal(flower));
        assert(flower_getBlockNumber(flower) == 0); //this should be true of the terminal problems.

        ///////////////////////////////////////////////////////////////////////////
        // Fill adjencencies
        ///////////////////////////////////////////////////////////////////////////

        startTime = time(NULL);
        fillingIn_fillAdjacencies(flower);
        buildFaces_buildAndProcessFaces(flower);
        st_logInfo("Processed the flowers in: %" PRIi64 " seconds\n", time(NULL)
                - startTime);

        ///////////////////////////////////////////////////////////////////////////
        //Set the faces in the flower to 'built' status, which triggers the building
        //of faces for the flower.
        ///////////////////////////////////////////////////////////////////////////

        assert(!flower_builtFaces(flower));
        flower_setBuildFaces(flower, 1);
    }

    ///////////////////////////////////////////////////////////////////////////
    // (9) Write the flower to disk.
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);
    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk in: %" PRIi64 " seconds\n", time(NULL)
            - startTime);

    ///////////////////////////////////////////////////////////////////////////
    //(15) Clean up.
    ///////////////////////////////////////////////////////////////////////////

    //Destruct stuff
    startTime = time(NULL);
    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    st_logInfo("Cleaned stuff up and am finished in: %" PRIi64 " seconds\n", time(NULL)
            - startTime);
    return 0;

}
