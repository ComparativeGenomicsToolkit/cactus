
#include "sonLib.h"

/*
 * Miscellaneous basic functions used to deal with nodes and edges in the reference code.
 */

////////////////////////////////////
////////////////////////////////////
//Functions for building sets containing nodes and edges
////////////////////////////////////
////////////////////////////////////

stSortedSet *getEmptyNodeOrEdgeSetWithCleanup() {
    return stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
}

stSortedSet *getEmptyNodeOrEdgeSetWithoutCleanup() {
    return stSortedSet_construct3(
            (int(*)(const void *, const void *)) stIntTuple_cmpFn, NULL);
}

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
    assert(stIntTuple_getPosition((stIntTuple *)edge, 0) < stIntTuple_getPosition((stIntTuple *)edge, 1));
    assert(stIntTuple_getPosition((stIntTuple *)edge2, 0) < stIntTuple_getPosition((stIntTuple *)edge2, 1));
    int32_t i = stIntTuple_getPosition((stIntTuple *) edge, 2);
    int32_t j = stIntTuple_getPosition((stIntTuple *) edge2, 2);
    return i > j ? 1 : (i < j ? -1 : 0);
}

int32_t getOtherPosition(stIntTuple *edge, int32_t node) {
    /*
     * Get the other position to the given node in the edge.
     */
    assert(stIntTuple_getPosition(edge, 0) < stIntTuple_getPosition(edge, 1));
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
        assert(stIntTuple_getPosition(edge, 0) < stIntTuple_getPosition(edge, 1));
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

void addWeightedEdgeToList(int32_t node1, int32_t node2, int32_t weight,
        stList *edges) {
    stList_append(edges, constructWeightedEdge(node1, node2, weight));
}

void addEdgeToList(int32_t node1, int32_t node2, stList *edges) {
    stList_append(edges, constructEdge(node1, node2));
}

stIntTuple *getWeightedEdgeFromSet(int32_t node1, int32_t node2,
        stSortedSet *allAdjacencyEdges) {
    if(node1 > node2) {
        return getWeightedEdgeFromSet(node2, node1, allAdjacencyEdges);
    }
    stIntTuple *edge = constructEdge(node1, node2);
    stIntTuple *edge2 = stSortedSet_searchGreaterThanOrEqual(allAdjacencyEdges,
            edge);
    stIntTuple_destruct(edge);
    assert(edge2 == NULL || stIntTuple_getPosition(edge2, 0) < stIntTuple_getPosition(edge2, 1));
    return edge2 != NULL && stIntTuple_getPosition(edge2, 0) == node1
            && stIntTuple_getPosition(edge2, 1) == node2 ? edge2 : NULL;
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
        assert(stIntTuple_getPosition(edge, 0) < stIntTuple_getPosition(edge, 1));
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
        assert(stIntTuple_getPosition(edge, 0) < stIntTuple_getPosition(edge, 1));
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
    stIntTuple *edge = constructEdge(node1, node2);
    bool b = stSortedSet_search(edges, edge) != NULL;
    stIntTuple_destruct(edge);
    return b;
}

void addEdgeToSet(stSortedSet *edges, int32_t node1, int32_t node2) {
    /*
     * Adds an edge to the set.
     */
    assert(!edgeInSet(edges, node1, node2));
    stIntTuple *edge = constructEdge(node1, node2);
    stSortedSet_insert(edges, edge);
    assert(edgeInSet(edges, node1, node2));
}

static bool hasGreaterThan0Weight(stIntTuple *edge) {
    return stIntTuple_getPosition(edge, 2) > 0;
}

stList *getEdgesWithGreaterThanZeroWeight(stList *adjacencyEdges) {
    /*
     * Gets all edges with weight greater than 0.
     */
    return stList_filter(adjacencyEdges, (bool(*)(void *)) hasGreaterThan0Weight);
}

void logEdges(stList *edges, const char *edgesName) {
    for (int32_t i = 0; i < stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        st_logDebug(
                "Edge, type: %s, node1: %i, node2: %i, weight: %i\n",
                edgesName,
                stIntTuple_getPosition(edge, 0),
                stIntTuple_getPosition(edge, 1),
                (stIntTuple_length(edge) == 3 ? stIntTuple_getPosition(edge, 2)
                        : INT32_MAX));
    }
}
