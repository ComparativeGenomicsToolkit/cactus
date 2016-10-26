/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

typedef struct _liftedEdge LiftedEdge;

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Ancestors function
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Utility function for the lifted edge hashtable
 */
static uint64_t buildFaces_hashfunction(const void *ptr) {
    Cap *key = (Cap *) ptr;
    return (uint64_t) cap_getName(key);
}

/*
 * Utility function for the lifted edge hashtable
 */
static int buildFaces_key_eq_fn(const void *ptrA, const void *ptrB) {
    return ptrA == ptrB;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Lifted edges functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct _liftedEdge {
    Cap *destination;
    Cap *bottomNode;
};

/*
 * Lifted edge destructor
 */
static void buildFaces_destructLiftedEdge(LiftedEdge * liftedEdge) {
    free(liftedEdge);
}
/*
 * Utility function for the lifted edge hashtable
 */
static void buildFaces_destructValue(void *ptr) {
    stList_destruct((stList *) ptr);
}

/*
 * Utility function for List struct
 */
static void buildFaces_stList_destructElem(void *ptr) {
    buildFaces_destructLiftedEdge((LiftedEdge *) ptr);
}

/*
 * Fill in a hashtable which to every node associates
 * alist of lifted edges
 */
static stHash *buildFaces_computeLiftedEdges(Flower * flower) {
    stHash *liftedEdgesTable = stHash_construct3(buildFaces_hashfunction,
            buildFaces_key_eq_fn, NULL, buildFaces_destructValue);
    Flower_CapIterator *iter = flower_getCapIterator(flower);
    Cap *cap, *attachedAncestor;
    Cap *adjacency, *adjacencyAncestor;
    stList *liftedEdges;
    LiftedEdge *liftedEdge;

    // Iterate through potential bottom nodes
    while ((cap = flower_getNextCap(iter))) {
        // ... check if connected
        if ((adjacency = cap_getAdjacency(cap))) {
            // ... lift
            attachedAncestor = cap_getTopCap(cap);
            adjacencyAncestor = cap_getTopCap(cap_getPositiveOrientation(
                    adjacency));

#ifndef NDEBUG
            assert((attachedAncestor && adjacencyAncestor) || (!attachedAncestor && !adjacencyAncestor));
#endif

            // If root node
            if (attachedAncestor == NULL)
                continue;

            // ... create lifted edge
            liftedEdge = st_malloc(sizeof(LiftedEdge));
            liftedEdge->destination = adjacencyAncestor;
            liftedEdge->bottomNode = cap;

#ifndef NDEBUG
            // Self loop
            if (adjacencyAncestor == attachedAncestor)
                abort();
#endif

            // ... add it to the hashtable
            if ((liftedEdges
                    = stHash_search(liftedEdgesTable, attachedAncestor))) {
                stList_append(liftedEdges, liftedEdge);
            } else {
                liftedEdges = stList_construct3(2,
                        buildFaces_stList_destructElem);
                stList_append(liftedEdges, liftedEdge);
                stHash_insert(liftedEdgesTable, attachedAncestor, liftedEdges);
            }
        }
    }

    flower_destructCapIterator(iter);
    return liftedEdgesTable;
}

/*
 * Recursive function which fills a givenlist with the
 * connected nodes within a module
 */
static void buildFaces_fillTopNodeList(Cap * cap, stList *list,
        stHash *liftedEdgesTable) {
    stList *liftedEdges;
    int64_t index;

    // Limit of recursion
    if (stList_contains(list, cap))
        return;

    // Actual filling
    st_logInfo("Adding cap %p to face\n", cap);
    stList_append(list, cap);

    // Recursion through lifted edges
    if ((liftedEdges = stHash_search(liftedEdgesTable, cap)))
        for (index = 0; index < stList_length(liftedEdges); index++)
            buildFaces_fillTopNodeList(
                    ((LiftedEdge *) stList_get(liftedEdges, index))->destination,
                   list, liftedEdgesTable);

    // Recursion through adjacency
    if (cap_getAdjacency(cap))
        buildFaces_fillTopNodeList(cap_getAdjacency(cap),list,
                liftedEdgesTable);
}

/*
 * Constructs a face from a given Cap
 */
static void buildFaces_constructFromCap(Cap * startingCap,
        stHash *liftedEdgesTable, Flower * flower) {
    Face *face = face_construct(flower);
    stList *topNodes = stList_construct3(16, NULL);
    stList *liftedEdges;
    Cap *cap, *bottomNode, *ancestor;
    int64_t index, index2;

    printf("Constructing new face");

    // Establishlist of top nodes
    buildFaces_fillTopNodeList(startingCap, topNodes, liftedEdgesTable);

#ifndef NDEBUG
    // What, no top nodes!?
    if (stList_length(topNodes) == 0)
        abort();
#endif

    // Initialize data structure
    face_allocateSpace(face, stList_length(topNodes));

    // For every top node
    for (index = 0; index < stList_length(topNodes); index++) {
        cap = stList_get(topNodes, index);
        face_setTopNode(face, index, cap);
        liftedEdges = stHash_search(liftedEdgesTable, cap);

        if (!liftedEdges) {
            face_setBottomNodeNumber(face, index, 0);
            continue;
        }

        face_setBottomNodeNumber(face, index, stList_length(liftedEdges));
        // For every bottom node of that top node
        for (index2 = 0; index2 < stList_length(liftedEdges); index2++) {
            bottomNode
                    = ((LiftedEdge *) stList_get(liftedEdges, index2))->bottomNode;
            face_addBottomNode(face, index, bottomNode);
            ancestor = cap_getTopCap(cap_getPositiveOrientation(
                    cap_getAdjacency(bottomNode)));
            if (cap_getAdjacency(cap) != ancestor)
                face_setDerivedDestination(face, index, index2, ancestor);
            else
                face_setDerivedDestination(face, index, index2, NULL);

#ifndef NDEBUG
            // If bottom nodes part of top nodes
            assert(!stList_contains(topNodes, cap_getPositiveOrientation(
                    ((LiftedEdge*) stList_get(liftedEdges, index2))->bottomNode)));
#endif
        }
    }

    // Clean up
    stList_destruct(topNodes);
}

/*
 * Fill uplist with bottom nodes for top node
 */
static void buildFaces_computeLiftedEdgesAtTopNode(Cap * cap, stList * liftedEdges) {
    int64_t childIndex;
    LiftedEdge * liftedEdge;

    // Termination
    if (cap_getAdjacency(cap)) {
	liftedEdge = st_malloc(sizeof(LiftedEdge));
	liftedEdge->bottomNode = cap_getPositiveOrientation(cap);
	liftedEdge->destination = cap_getPositiveOrientation(cap_getTopCap(cap_getAdjacency(cap)));
	stList_append(liftedEdges, liftedEdge);
	return;
    }
    
    // Recursion through children
    for (childIndex = 0; childIndex < cap_getChildNumber(cap); childIndex++) 
	buildFaces_computeLiftedEdgesAtTopNode(cap_getChild(cap, childIndex), liftedEdges);
}

/*
 * Recursive function which fills a givenlist with the
 * connected nodes within a module and fills their lifted
 * edges in the same pass
 */
static void buildFaces_fillTopNodeList2(Cap * cap, stList *list,
        stHash *liftedEdgesTable) {
    stList *liftedEdges = stList_construct3(2,
                        buildFaces_stList_destructElem);
    int64_t index;

    // Orientation check
    cap = cap_getPositiveOrientation(cap);

    // Limit of recursion
    if (stList_contains(list, cap))
        return;

    // Actual filling
    st_logInfo("Adding cap %p to face\n", cap);
    stList_append(list, cap);

    // Compute lifted edges
    for (index = 0; index < cap_getChildNumber(cap); index++) 
	buildFaces_computeLiftedEdgesAtTopNode(cap_getChild(cap, index), liftedEdges);

    // If emptylist...
    if (stList_length(liftedEdges) == 0) 
	stList_destruct(liftedEdges);
    // Recursion through lifted edges
    else {
	stHash_insert(liftedEdgesTable, cap, liftedEdges);
        for (index = 0; index < stList_length(liftedEdges); index++)
            buildFaces_fillTopNodeList2(
                    ((LiftedEdge *) stList_get(liftedEdges, index))->destination,
                   list, liftedEdgesTable);
    }

    // Recursion through adjacency
    if (cap_getAdjacency(cap))
        buildFaces_fillTopNodeList2(cap_getAdjacency(cap),list,
                liftedEdgesTable);
}


/*
 * Constructs a face locally from a given Cap but without precomputed liftedEdges
 */
void buildFaces_reconstructFromCap(Cap * startingCap, Flower * flower) {
    Face *face = face_construct(flower);
    stList * liftedEdges;
    stList *topNodes = stList_construct3(16, NULL);
    stHash *liftedEdgesTable = stHash_construct3(buildFaces_hashfunction,
            buildFaces_key_eq_fn, NULL, buildFaces_destructValue);
    Cap *cap, *bottomNode, *ancestor;
    int64_t index, index2;

    printf("Constructing new face");

    // Establishlist of top nodes and fill liftedEdges table
    buildFaces_fillTopNodeList2(startingCap, topNodes, liftedEdgesTable);

#ifndef NDEBUG
    // What, no top nodes!?
    assert(stList_length(topNodes));
#endif

    // Initialize data structure
    face_allocateSpace(face, stList_length(topNodes));

    // For every top node
    for (index = 0; index < stList_length(topNodes); index++) {
        cap = stList_get(topNodes, index);
        face_setTopNode(face, index, cap);
        liftedEdges = stHash_search(liftedEdgesTable, cap);

        if (!liftedEdges) {
            face_setBottomNodeNumber(face, index, 0);
            continue;
        }

        face_setBottomNodeNumber(face, index, stList_length(liftedEdges));
        // For every bottom node of that top node
        for (index2 = 0; index2 < stList_length(liftedEdges); index2++) {
            bottomNode
                    = ((LiftedEdge *) stList_get(liftedEdges, index2))->bottomNode;
            face_addBottomNode(face, index, bottomNode);

            assert(cap_getAdjacency(bottomNode));
            ancestor = cap_getTopCap(cap_getPositiveOrientation(
                    cap_getAdjacency(bottomNode)));
            if (cap_getAdjacency(cap) != ancestor)
                face_setDerivedDestination(face, index, index2, ancestor);
            else
                face_setDerivedDestination(face, index, index2, NULL);

#ifndef NDEBUG
            // If bottom nodes part of top nodes
            if (stList_contains(topNodes, cap_getPositiveOrientation(
                    ((LiftedEdge*) stList_get(liftedEdges, index2))->bottomNode)))
                abort();
#endif
        }
    }

    // Clean up
    stList_destruct(topNodes);
    stHash_destruct(liftedEdgesTable);
}

void flower_reconstructFaces(Flower * flower) {
    flower_destructFaces(flower);
    stHash *liftedEdgesTable = buildFaces_computeLiftedEdges(flower);
    Flower_CapIterator *iter = flower_getCapIterator(flower);
    stList *liftedEdges;
    Cap *current;

    while ((current = flower_getNextCap(iter))) {
        if ((liftedEdges = stHash_search(liftedEdgesTable, current))
                && (stList_length(liftedEdges) >= 1)) {
            buildFaces_constructFromCap(current, liftedEdgesTable, flower);
        }
    }
    stHash_destruct(liftedEdgesTable);
    flower_destructCapIterator(iter);
}

void flower_destructFaces(Flower *flower) {
    Face *face;
    while ((face = flower_getFirstFace(flower)) != NULL) {
        face_destruct(face);
    }
}
