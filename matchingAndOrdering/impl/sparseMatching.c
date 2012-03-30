#include "sonLib.h"
#include "stCheckEdges.h"
#include "shared.h"
#include "stMatchingAlgorithms.h"

static stHash *rebaseNodes(stSortedSet *nodes) {
    /*
     * Renumber the nodes from 0 contiguously.
     */
    stHash *nodesToRebasedNodes = stHash_construct3(
            (uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn,
            NULL,
            (void(*)(void *)) stIntTuple_destruct);
    int32_t i=0;
    stSortedSetIterator *it = stSortedSet_getIterator(nodes);
    stIntTuple *node;
    while((node = stSortedSet_getNext(it)) != NULL) {
        stHash_insert(nodesToRebasedNodes, node, stIntTuple_construct(1, i++));
    }
    stSortedSet_destructIterator(it);
    return nodesToRebasedNodes;
}

static stList *translateEdges(stList *edges, stHash *nodesToRebasedNodes) {
    /*
     * Translate the edges.
     */
    stList *rebasedEdges = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for(int32_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        stIntTuple *node1 = getItemForNode(stIntTuple_getPosition(edge, 0), nodesToRebasedNodes);
        stIntTuple *node2 = getItemForNode(stIntTuple_getPosition(edge, 1), nodesToRebasedNodes);
        assert(node1 != NULL);
        assert(node2 != NULL);
        stList_append(rebasedEdges, constructWeightedEdge(stIntTuple_getPosition(node1, 0), stIntTuple_getPosition(node2, 0), stIntTuple_getPosition(edge, 2)));
    }
    return rebasedEdges;
}

static stList *translateEdges2(stList *rebasedEdges, stHash *rebasedNodesToNodes, stList *originalEdges) {
    /*
     * Translate the edges.
     */
    stSortedSet *originalEdgesSet = stList_getSortedSet(originalEdges, (int (*)(const void *, const void *))stIntTuple_cmpFn);
    stList *edges = stList_construct();
    for(int32_t i=0; i<stList_length(rebasedEdges); i++) {
        stIntTuple *rebasedEdge = stList_get(rebasedEdges, i);
        stIntTuple *node1 = getItemForNode(stIntTuple_getPosition(rebasedEdge, 0), rebasedNodesToNodes);
        stIntTuple *node2 = getItemForNode(stIntTuple_getPosition(rebasedEdge, 1), rebasedNodesToNodes);
        assert(node1 != NULL);
        assert(node2 != NULL);
        stIntTuple *edge = getWeightedEdgeFromSet(stIntTuple_getPosition(node1, 0), stIntTuple_getPosition(node2, 0), originalEdgesSet);
        assert(edge != NULL);
        stList_append(edges, edge);
    }
    stSortedSet_destruct(originalEdgesSet);
    return edges;
}

stList *getSparseMatching(stSortedSet *nodes,
        stList *adjacencyEdges,
        stList *(*matchingAlgorithm)(stList *edges, int32_t nodeNumber)) {
    checkEdges(adjacencyEdges, nodes, 0, 0);

    if (stSortedSet_size(nodes) == 0) { //Some of the following functions assume there are at least 2 nodes.
        return stList_construct();
    }

    /*
     * First calculate the optimal matching.
     */
    stHash *nodesToRebasedNodes = rebaseNodes(nodes);
    stHash *rebasedNodesToNodes = stHash_invert(nodesToRebasedNodes, (uint32_t(*)(const void *)) stIntTuple_hashKey,
            (int(*)(const void *, const void *)) stIntTuple_equalsFn, NULL, NULL);
    stList *rebasedEdges = translateEdges(adjacencyEdges, nodesToRebasedNodes);
    stList *chosenRebasedEdges = matchingAlgorithm(rebasedEdges, stSortedSet_size(nodes));
    stList *chosenEdges = translateEdges2(chosenRebasedEdges, rebasedNodesToNodes, adjacencyEdges);

    /*
     * Clean up
     */
    stList_destruct(chosenRebasedEdges);
    stList_destruct(rebasedEdges);
    stHash_destruct(nodesToRebasedNodes);
    stHash_destruct(rebasedNodesToNodes);

    st_logDebug(
            "Chosen a sparse matching with %i edges, %i cardinality and %i weight\n",
            stList_length(chosenEdges), matchingCardinality(chosenEdges),
            matchingWeight(chosenEdges));

    return chosenEdges;
}

