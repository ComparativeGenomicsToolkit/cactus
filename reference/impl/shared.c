
#include "sonLib.h"

/*
 * Miscellaneous basic functions used to deal with nodes and edges in the reference code.
 */

////////////////////////////////////
////////////////////////////////////
//Edge functions
////////////////////////////////////
////////////////////////////////////

stIntTuple *constructEdge(int32_t node1, int32_t node2) {
    assert(node1 != node2);
    if (node1 > node2) {
        return constructEdge(node2, node1);
    }
    return stIntTuple_construct(2, node1, node2);
}

stIntTuple *constructWeightedEdge(int32_t node1, int32_t node2, int32_t weight) {
    assert(node1 != node2);
    if (node1 > node2) {
        return constructWeightedEdge(node2, node1, weight);
    }
    return stIntTuple_construct(3, node1, node2, weight);
}

int compareEdgesByWeight(const void *edge, const void *edge2) {
    int32_t i = stIntTuple_getPosition((stIntTuple *) edge, 2);
    int32_t j = stIntTuple_getPosition((stIntTuple *) edge2, 2);
    return i > j ? 1 : (i < j ? -1 : 0);
}

int32_t getOtherPosition(stIntTuple *edge, int32_t node) {
    /*
     * Get the other position to the given node in the edge.
     */
    if (stIntTuple_getPosition(edge, 0) == node) {
        return stIntTuple_getPosition(edge, 1);
    } else {
        assert(stIntTuple_getPosition(edge, 1) == node);
        return stIntTuple_getPosition(edge, 0);
    }
}

////////////////////////////////////
////////////////////////////////////
//Node functions
////////////////////////////////////
////////////////////////////////////

static void getNodeSetOfEdgesP(stSortedSet *nodes, int32_t node) {
    stIntTuple *i = stIntTuple_construct(1, node);
    if (stSortedSet_search(nodes, i) == NULL) {
        stSortedSet_insert(nodes, i);
    } else {
        stIntTuple_destruct(i);
    }
}

stSortedSet *getNodeSetOfEdges(stList *edges) {
    /*
     * Returns a sorted set of the nodes in the given list of edges.
     */
    stSortedSet *nodes = stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        getNodeSetOfEdgesP(nodes, stIntTuple_getPosition(edge, 0));
        getNodeSetOfEdgesP(nodes, stIntTuple_getPosition(edge, 1));
    }
    return nodes;
}

bool nodeInSet(stSortedSet *nodes, int32_t node) {
    /*
     * Returns non zero iff the node is in the set.
     */
    stIntTuple *i = stIntTuple_construct(1, node);
    bool b = stSortedSet_search(nodes, i) != NULL;
    stIntTuple_destruct(i);
    return b;
}

void addNodeToSet(stSortedSet *nodes, int32_t node) {
    /*
     * Adds a node to the set.
     */
    assert(!nodeInSet(nodes, node));
    stIntTuple *i = stIntTuple_construct(1, node);
    stSortedSet_insert(nodes, i);
}

void *getItemForNode(int32_t node, stHash *nodesToItems) {
    /*
     * Get edges for given node, or NULL if none.
     */
    stIntTuple *node1 = stIntTuple_construct(1, node);
    void *o = stHash_search(nodesToItems, node1);
    stIntTuple_destruct(node1);
    return o;
}

////////////////////////////////////
////////////////////////////////////
//Basic node and edge functions
////////////////////////////////////
////////////////////////////////////

static stIntTuple *getWeightedEdgeFromSetP(int32_t node1, int32_t node2,
        stSortedSet *allAdjacencyEdges) {
    stIntTuple *edge = stIntTuple_construct(3, node1, node2, 0);
    stIntTuple *edge2 = stSortedSet_searchGreaterThanOrEqual(allAdjacencyEdges,
            edge);
    stIntTuple_destruct(edge);
    return edge2 != NULL && stIntTuple_getPosition(edge2, 0) == node1
            && stIntTuple_getPosition(edge2, 1) == node2 ? edge2 : NULL;
}

stIntTuple *getWeightedEdgeFromSet(int32_t node1, int32_t node2,
        stSortedSet *allAdjacencyEdges) {
    /*
     * Gets the edge in the set connecting node1 and node2. Returns NULL if doesn't exist.
     */
    assert(node1 != node2);
    stIntTuple *edge1 =
            getWeightedEdgeFromSetP(node1, node2, allAdjacencyEdges);
    stIntTuple *edge2 =
            getWeightedEdgeFromSetP(node2, node1, allAdjacencyEdges);
    assert(edge1 == NULL || edge2 == NULL);
    return edge1 != NULL ? edge1 : edge2;
}


static void getNodesToEdgesHashP(stHash *nodesToEdges, stIntTuple *edge,
        int32_t position) {
    stIntTuple *node = stIntTuple_construct(1,
            stIntTuple_getPosition(edge, position));
    stList *edges;
    if ((edges = stHash_search(nodesToEdges, node)) == NULL) {
        edges = stList_construct();
        stHash_insert(nodesToEdges, node, edges);
    } else {
        stIntTuple_destruct(node);
    }
    stList_append(edges, edge);
}

stHash *getNodesToEdgesHash(stList *edges) {
    /*
     * Build a hash of nodes to edges.
     */
    stHash *nodesToEdges = stHash_construct3(
            (uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn,
            (void(*)(void *)) stIntTuple_destruct,
            (void(*)(void *)) stList_destruct);

    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        getNodesToEdgesHashP(nodesToEdges, edge, 0);
        getNodesToEdgesHashP(nodesToEdges, edge, 1);
    }

    return nodesToEdges;
}

stIntTuple *getEdgeForNodes(int32_t node1, int32_t node2,
        stHash *nodesToAdjacencyEdges) {
    /*
     * Searches for any edge present in nodesToAdjacencyEdges that links node and node2.
     */
    stList *edges = getItemForNode(node1, nodesToAdjacencyEdges);
    if (edges == NULL) {
        return NULL;
    }
    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        if (getOtherPosition(edge, node1) == node2) {
            return edge;
        }
    }
    return NULL;
}

bool edgeInSet(stSortedSet *edges, int32_t node1, int32_t node2) {
    /*
     * Returns non-zero iff the edge linking node1 and node2 is in the set of edges. Works with addEdgeToSet.
     * Assumes that that node1 < node2
     */
    if (node1 > node2) {
        return edgeInSet(edges, node2, node1);
    }
    stIntTuple *edge = stIntTuple_construct(2, node1, node2);
    bool b = stSortedSet_search(edges, edge) != NULL;
    stIntTuple_destruct(edge);
    return b;
}

void addEdgeToSet(stSortedSet *edges, int32_t node1, int32_t node2) {
    /*
     * Adds an edge to the set.
     */
    if (node1 > node2) {
        addEdgeToSet(edges, node2, node1);
        return;
    }
    assert(!edgeInSet(edges, node1, node2));
    assert(!edgeInSet(edges, node2, node1));
    stIntTuple *edge = stIntTuple_construct(2, node1, node2);
    stSortedSet_insert(edges, edge);
    assert(edgeInSet(edges, node1, node2));
}
