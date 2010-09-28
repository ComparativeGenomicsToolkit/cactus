#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "pinchGraph.h"
#include "cactusGraph.h"
#include "cactus.h"
#include "cactusFlowerFunctions.h"
#include "avl.h"
#include "sonLib.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions to construct pinch graphs from flowers.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchEdge *hookUpEdge(struct Piece *piece,
        struct PinchGraph *pinchGraph, struct PinchVertex *vertex2,
        struct PinchVertex *vertex3) {
    struct PinchEdge *edge;

    edge = constructPinchEdge(piece);

    //Connect up each end of the black edge.
    edge->from = vertex2;
    edge->rEdge->to = vertex2;
    insertBlackEdge(vertex2, edge);

    edge->to = vertex3;
    edge->rEdge->from = vertex3;
    insertBlackEdge(vertex3, edge->rEdge);

    //Now add pieces connected to edges to the graph.
    avl_insert(pinchGraph->edges, edge);
    avl_insert(pinchGraph->edges, edge->rEdge);

    return edge;
}

struct PinchGraph *constructPinchGraph(Flower *flower) {
    struct PinchGraph *graph;
    Flower_EndIterator *endIterator;
    End_InstanceIterator *instanceIterator;
    End *end;
    Cap *cap;
    Cap *cap2;
    Sequence *sequence;
    struct PinchVertex *sourceVertex;
    struct PinchVertex *pinchVertex;
    struct PinchVertex *pinchVertex2;
    struct PinchEdge *leftCapEdge;
    struct PinchEdge *edge = NULL;
    struct PinchEdge *rightCapEdge;
    struct hashtable *hash;
    struct hashtable *hash2;
    int32_t start;
    int32_t stop;
    int32_t length;

    assert(!flower_builtBlocks(flower));
    assert(flower_isTerminal(flower));

    //make basic object.
    graph = pinchGraph_construct();
    sourceVertex = graph->vertices->list[0];

    //make hashes for ends to vertices
    hash = create_hashtable(flower_getEndNumber(flower) * 2,
            hashtable_stringHashKey, hashtable_stringEqualKey, free, NULL);
    hash2 = create_hashtable(flower_getEndNumber(flower) * 2,
            hashtable_stringHashKey, hashtable_stringEqualKey, free, NULL);

    //for each cap, build a pair of vertices
    endIterator = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        assert(!end_isBlockEnd(end));
        pinchVertex = constructPinchVertex(graph, -1, 0, 1);
        pinchVertex2 = constructPinchVertex(graph, -1, 1, 0);
        //connect to source.
        if (end_isAttached(end)) {
            connectVertices(sourceVertex, pinchVertex);
        }
        hashtable_insert(hash, cactusMisc_nameToString(end_getName(end)),
                pinchVertex);
        hashtable_insert(hash2, cactusMisc_nameToString(end_getName(end)),
                pinchVertex2);
    }
    flower_destructEndIterator(endIterator);

    endIterator = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        instanceIterator = end_getInstanceIterator(end);
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            cap2 = cap_getAdjacency(cap);
            sequence = cap_getSequence(cap);

            assert(cap2 != NULL);
            assert(cap_getStrand(cap2));
            assert(sequence == cap_getSequence(cap2));

            //if(length >= 0)  {
            if (!cap_getSide(cap)) {
                assert(cap_getSide(cap2));

                start = cap_getCoordinate(cap);
                stop = cap_getCoordinate(cap2);
                length = stop - start - 1;
                assert(length >= 0);

                //Make black edges for caps/stubs on left end
                leftCapEdge = hookUpEdge(constructPiece(cap_getName(cap),
                        start, start), graph, hashtable_search(hash,
                        (void *) cactusMisc_nameToStringStatic(end_getName(
                                cap_getEnd(cap)))), hashtable_search(hash2,
                        (void *) cactusMisc_nameToStringStatic(end_getName(
                                cap_getEnd(cap)))));

                //Construct the middle sequence, if not zero length.
                if (length > 0) {
                    edge = hookUpEdge(constructPiece(
                            sequence_getName(sequence), start + 1, stop - 1),
                            graph, constructPinchVertex(graph, -1, 0, 0),
                            constructPinchVertex(graph, -1, 0, 0));
                }

                //Construct the right cap/stub
                rightCapEdge = hookUpEdge(constructPiece(cap_getName(cap2),
                        stop, stop), graph, hashtable_search(hash2,
                        (void *) cactusMisc_nameToStringStatic(end_getName(
                                cap_getEnd(cap2)))), hashtable_search(hash,
                        (void *) cactusMisc_nameToStringStatic(end_getName(
                                cap_getEnd(cap2)))));

                //Connect the edges
                if (length > 0) {
                    connectVertices(leftCapEdge->to, edge->from);
                    connectVertices(edge->to, rightCapEdge->from);
                } else {
                    connectVertices(leftCapEdge->to, rightCapEdge->from);
                }
            }
        }
        end_destructInstanceIterator(instanceIterator);
    }
    flower_destructEndIterator(endIterator);
    //Cleanup the hashes
    hashtable_destroy(hash, FALSE, TRUE);
    hashtable_destroy(hash2, FALSE, TRUE);
    return graph;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions to construct flowers.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////


struct CactusEdge *getNonDeadEndOfStubCactusEdge(struct CactusEdge *edge,
        struct PinchGraph *pinchGraph) {
    struct PinchEdge *pinchEdge;
    pinchEdge = cactusEdgeToFirstPinchEdge(edge, pinchGraph);
    assert(isAStubCactusEdge(edge, pinchGraph));
    assert(vertex_isDeadEnd(pinchEdge->from) || vertex_isDeadEnd(pinchEdge->to));
    return vertex_isDeadEnd(pinchEdge->from) ? edge->rEdge : edge;
}

Name cactusEdgeToEndName(struct CactusEdge *edge,
        struct hashtable *endNamesHash, struct PinchGraph *pinchGraph) {
    struct PinchEdge *pinchEdge = cactusEdgeToFirstPinchEdge(edge, pinchGraph);
    assert(pinchEdge != NULL);
    char *cA = (char *) hashtable_search(endNamesHash, pinchEdge->from);
    assert(cA != NULL);
    return cactusMisc_stringToName(cA);
}

Sequence *copySequence(Flower *flower, Name name) {
    Sequence *sequence = flower_getSequence(flower, name);
    if (sequence == NULL) {
        sequence = sequence_construct(cactusDisk_getMetaSequence(
                flower_getCactusDisk(flower), name), flower);
    }
    return sequence;
}

static Block *constructBlockFromCactusEdge(struct CactusEdge *edge,
        Flower *flower) {
    /*
     * Constructs a block and two connected ends.
     */
    int32_t i;
    Block *block;
    Sequence *sequence;
    struct Piece *piece;
    piece = edge->pieces->list[0];
    block = block_construct(piece->end - piece->start + 1, flower);
    for (i = 0; i < edge->pieces->length; i++) {
        piece = edge->pieces->list[i];
        sequence = copySequence(flower, piece->contig);
        segment_construct2(block,
                piece->start > 0 ? piece->start : -piece->end,
                piece->start > 0, sequence);
    }
    return block;
}

static struct List *addEnvelopedStubEnds(Flower *flower, int32_t addToFlower) {
    /*
     * For each flower contained within a link in a chain, adds the encompassing ends
     * of the chain to the nested flower.
     */
    int32_t j;
    End *end, *end2;
    Flower *flower2;
    struct List *list;
    Flower_EndIterator *endIterator;
    Group_EndIterator *adjacencyIterator;
    Group *group;

    adjacencyIterator = flower_getGroupIterator(flower);
    while ((group = flower_getNextGroup(adjacencyIterator)) != NULL) {
        flower2 = group_getNestedFlower(group);
        if (flower2 != NULL) {
            list = addEnvelopedStubEnds(flower2, 1);
            for (j = 0; j < list->length; j++) {
                end = list->list[j];
                if (addToFlower && flower_getEnd(flower, end_getName(end))
                        == NULL) {
                    end_setGroup(end_copyConstruct(end, flower), group);
                } else {
                    end2 = flower_getEnd(flower, end_getName(end));
                    assert(end2 != NULL);
                    if (end_getGroup(end2) == NULL) {
                        end_setGroup(end2, group);
                    } else {
                        assert(end_getGroup(end2) == group);
                    }
                }
            }
            destructList(list);
        }
    }
    flower_destructGroupIterator(adjacencyIterator);

    list = constructEmptyList(0, NULL);
    endIterator = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        if (end_isStubEnd(end)) {
            listAppend(list, end);
        }
    }
    flower_destructEndIterator(endIterator);
    return list;
}

static int addAdjacenciesToLeafCapsP(Cap **cap1, Cap **cap2) {
    assert(cap_getStrand(*cap1) && cap_getStrand(*cap2));
    Sequence *sequence1 = cap_getSequence(*cap1);
    Sequence *sequence2 = cap_getSequence(*cap2);
    int32_t i = cactusMisc_nameCompare(sequence_getName(sequence1),
            sequence_getName(sequence2));
    if (i == 0) {
        int32_t j = cap_getCoordinate(*cap1);
        int32_t k = cap_getCoordinate(*cap2);
        i = j - k;
        if (i == 0) {
            assert(cap_getSegment(*cap1) == cap_getSegment(*cap2));
            j = cap_getSide(*cap1);
            k = cap_getSide(*cap2);
            assert((j && !k) || (!j && k));
            i = j ? -1 : 1;
        }
    }
    return i;
}

void addAdjacenciesToLeafCaps(Flower *flower) {
    End *end;
    Cap *cap;
    Cap *cap2;
    Flower_EndIterator *endIterator;
    End_InstanceIterator *instanceIterator;
    struct List *list;
    int32_t i;

    /*
     * Build a list of caps, then sort them.
     */
    list = constructEmptyList(0, NULL);
    endIterator = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        instanceIterator = end_getInstanceIterator(end);
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            if (!cap_isInternal(cap)) {
                if (!cap_getStrand(cap)) {
                    cap = cap_getReverse(cap);
                }
                listAppend(list, cap);
            }
        }
        end_destructInstanceIterator(instanceIterator);
    }
    flower_destructEndIterator(endIterator);
    assert((list->length % 2) == 0);

    /*
     * Sort the caps.
     */
    qsort(list->list, list->length, sizeof(void *), (int(*)(const void *v,
            const void *)) addAdjacenciesToLeafCapsP);

    /*
     * Now make the adjacencies.
     */
    for (i = 1; i < list->length; i += 2) {
        cap = list->list[i - 1];
        cap2 = list->list[i];
        cap_makeAdjacent(cap, cap2);
    }

    /*
     * Clean up.
     */
    destructList(list);
}

void addAdjacenciesToEnds(Flower *flower) {
    /*
     * Add the adjacencies between leaf caps of a flower. Does not touch internal instances.
     * All ends must be in the flower before this function is called.
     */
    addAdjacenciesToLeafCaps(flower);
    Group_EndIterator *adjacencyIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(adjacencyIterator)) != NULL) {
        if (group_getNestedFlower(group) != NULL) {
            addAdjacenciesToEnds(group_getNestedFlower(group));
        }
    }
    flower_destructGroupIterator(adjacencyIterator);
}

static int32_t returnsTrue(Event *event) {
    assert(event != NULL);
    return 1;
}

bool groupIsZeroLength(struct List *endNames, Flower *flower) {
    int32_t i;

    for (i = 0; i < endNames->length; i++) {
        End *end = flower_getEnd(flower, cactusMisc_stringToName(
                endNames->list[i]));
        End_InstanceIterator *iterator = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(iterator)) != NULL) {
            Cap *cap2 = cap_getAdjacency(cap);
            assert(cap2 != NULL);
            uint32_t j = abs(cap_getCoordinate(cap) - cap_getCoordinate(cap2));
            assert(j != 0);
            if (j > 1) {
                return 0;
            }
        }
        end_destructInstanceIterator(iterator);
    }
    return 1;
}

void addGroupsP(End *end, Group *group) {
    End_InstanceIterator *capIt = end_getInstanceIterator(end);
    Cap *cap;
    assert(end_getGroup(end) == NULL);
    end_setGroup(end, group);
    assert(end_getGroup(end) == group);
    while ((cap = end_getNext(capIt)) != NULL) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        End *adjacentEnd = cap_getEnd(adjacentCap);
        if (end_getGroup(adjacentEnd) == NULL) {
            addGroupsP(adjacentEnd, group);
        } else {
            assert(end_getGroup(adjacentEnd) == group);
        }
    }
    end_destructInstanceIterator(capIt);
}

void addGroups(Flower *flower) {
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        Group *group = end_getGroup(end);
        if (group == NULL) {
            group = group_construct2(flower);
            addGroupsP(end, group);
        }
    }
    flower_destructEndIterator(endIt);

    //Now call recursively and make any little chains
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        //Call recursively, if necessary
        if (!group_isLeaf(group)) {
            addGroups(group_getNestedFlower(group));
        }
        //Make a little chain for the ends, if the group is empty.
        group_constructChainForLink(group);
    }
    flower_destructGroupIterator(groupIt);
}

static int32_t *vertexDiscoveryTimes;

#ifdef BEN_DEBUG
void checkBiConnectedComponent(struct List *biConnnectedComponent) {
    int32_t i, j;
    struct CactusEdge *cactusEdge;
    struct CactusEdge *cactusEdge2;
    assert(biConnnectedComponent->length> 0);
    j = INT32_MIN;
    for (i = 0; i < biConnnectedComponent->length; i++) {
        cactusEdge = biConnnectedComponent->list[i];
        assert(j < vertexDiscoveryTimes[cactusEdge->from->vertexID]);
        j = vertexDiscoveryTimes[cactusEdge->from->vertexID];
    }
    cactusEdge = biConnnectedComponent->list[0];
    cactusEdge2
            = biConnnectedComponent->list[biConnnectedComponent->length - 1];
    assert(vertexDiscoveryTimes[cactusEdge->from->vertexID] == vertexDiscoveryTimes[cactusEdge2->to->vertexID]);
}
#endif

int fillOutFlowerFromInputsP2(struct List **biConnectedComponent1,
        struct List **biConnectedComponent2) {
    struct CactusEdge *cactusEdge1;
    struct CactusEdge *cactusEdge2;
    int32_t i, j;
    cactusEdge1 = (*biConnectedComponent1)->list[0];
    cactusEdge2 = (*biConnectedComponent2)->list[0];
    i = vertexDiscoveryTimes[cactusEdge1->from->vertexID];
    j = vertexDiscoveryTimes[cactusEdge2->from->vertexID];
    return i - j;
}

void mergeCactusVertices(struct CactusEdge *cactusEdge,
        int32_t *mergedVertexIDs, int32_t j, struct List *biConnectedComponent) {
    int32_t k;
    struct CactusVertex *cactusVertex;
    //merge vertices
    if (vertexDiscoveryTimes[cactusEdge->from->vertexID]
            < vertexDiscoveryTimes[cactusEdge->to->vertexID]) {
        mergedVertexIDs[cactusEdge->to->vertexID]
                = mergedVertexIDs[cactusEdge->from->vertexID];
    } else if (vertexDiscoveryTimes[cactusEdge->from->vertexID]
            > vertexDiscoveryTimes[cactusEdge->to->vertexID]) {
        assert(j == biConnectedComponent->length-1);
        for (k = 0; k <= j; k++) {
            cactusVertex
                    = ((struct CactusEdge *) biConnectedComponent->list[k])->from;
            if (mergedVertexIDs[cactusVertex->vertexID]
                    == mergedVertexIDs[cactusEdge->from->vertexID]) {
                mergedVertexIDs[cactusVertex->vertexID]
                        = mergedVertexIDs[cactusEdge->to->vertexID];
            }
        }
    }
}

static void setBlocksBuilt(Flower *flower) {
    /*
     * Sets the 'built-blocks flag' for all the flowers in the subtree, including the given flower.
     */
    assert(!flower_builtBlocks(flower));
    flower_setBuiltBlocks(flower, 1);
    assert(flower_builtBlocks(flower));
    Flower_GroupIterator *iterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(iterator)) != NULL) {
        if (!group_isLeaf(group)) {
            setBlocksBuilt(group_getNestedFlower(group));
        }
    }
    flower_destructGroupIterator(iterator);
}

static bool getOrientationP(struct CactusEdge *cactusEdge,
        struct PinchGraph *pinchGraph, Flower *flower,
        struct hashtable *endNamesHash) {
    cactusEdge = getNonDeadEndOfStubCactusEdge(cactusEdge, pinchGraph);
    End *end = flower_getEnd(flower, cactusEdgeToEndName(cactusEdge,
            endNamesHash, pinchGraph));
    assert(end != NULL);
    return end_getSide(end);
}

void reverseComponent(struct List *biConnectedComponent) {
    listReverse(biConnectedComponent);
    for (int32_t i = 0; i < biConnectedComponent->length; i++) {
        struct CactusEdge *cactusEdge = biConnectedComponent->list[i];
        biConnectedComponent->list[i] = cactusEdge->rEdge;
    }
}

void getOrientation(struct List *biConnectedComponent,
        struct PinchGraph *pinchGraph, Flower *flower,
        struct hashtable *endNamesHash) {
    //Make the blocks and ends
    assert(biConnectedComponent->length > 0);
    struct CactusEdge *cactusEdge = biConnectedComponent->list[0];
    if (isAStubCactusEdge(cactusEdge, pinchGraph)) {
        if (getOrientationP(cactusEdge, pinchGraph, flower, endNamesHash)) { //Reverse the list
            reverseComponent(biConnectedComponent);
        }
    } else {
        cactusEdge = biConnectedComponent->list[biConnectedComponent->length
                - 1];
        if (isAStubCactusEdge(cactusEdge, pinchGraph)) {
            if (!getOrientationP(cactusEdge, pinchGraph, flower, endNamesHash)) {
                reverseComponent(biConnectedComponent);
            }
        }
    }
}

void fillOutFlowerFromInputs(Flower *parentFlower,
        struct CactusGraph *cactusGraph, struct PinchGraph *pinchGraph,
        stSortedSet *chosenBlocks) {
    Flower *flower;
    Flower *nestedFlower;
    End *end;
    End *end2;
    Block *block;
    Flower_EndIterator *endIterator;
    Cap *cap;
    Chain *chain;
    Group *group;
    //struct CactusVertex *cactusVertex;
    struct CactusEdge *cactusEdge;
    struct CactusEdge *cactusEdge2;
    struct List *biConnectedComponent;
    struct List *list; // *list2;
    struct List *biConnectedComponents;
    void **flowers;
    void **parentFlowers;
    int32_t *mergedVertexIDs;
    int32_t i, j; //, k;
    struct hashtable *endNamesHash;
    struct PinchEdge *pinchEdge;
    struct Piece *piece;

    ////////////////////////////////////////////////
    //Get sorted bi-connected components (sorted as in ordered from root of vertex)
    ////////////////////////////////////////////////

    biConnectedComponents = computeSortedBiConnectedComponents(cactusGraph);

    st_logDebug("Built the bi-connected components: %i\n",
            biConnectedComponents->length);

    ////////////////////////////////////////////////
    //Get DFS numbering on cactus vertices
    ////////////////////////////////////////////////

    vertexDiscoveryTimes = getDFSDiscoveryTimes(cactusGraph);
#ifdef BEN_DEBUG
    for (i = 0; i < biConnectedComponents->length; i++) { //checks the discovery times are as expected.
        checkBiConnectedComponent(biConnectedComponents->list[i]);
    }
#endif
    st_logDebug("Got the vertex discovery times\n");

    ////////////////////////////////////////////////
    //Sort the biconnected components by their start time
    ////////////////////////////////////////////////

    qsort(biConnectedComponents->list, biConnectedComponents->length,
            sizeof(void *),
            (int(*)(const void *v, const void *)) fillOutFlowerFromInputsP2);
    st_logDebug("Sorted the biconnected components by vertex discovery time\n");

    ////////////////////////////////////////////////
    //Build end names hash
    ////////////////////////////////////////////////

    endNamesHash = create_hashtable(stSortedSet_size(chosenBlocks),
            hashtable_key, hashtable_equalKey, NULL, free);
    endIterator = flower_getEndIterator(parentFlower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        cap = end_getFirst(end);
        pinchEdge = getContainingBlackEdge(pinchGraph, cap_getName(cap),
                cap_getCoordinate(cap));
        assert(pinchEdge != NULL);
        if (vertex_isEnd(pinchEdge->from)) {
            assert(vertex_isDeadEnd(pinchEdge->to));
            hashtable_insert(endNamesHash, pinchEdge->from,
                    cactusMisc_nameToString(end_getName(end)));
        } else {
            assert(vertex_isEnd(pinchEdge->to));
            assert(vertex_isDeadEnd(pinchEdge->from));
            hashtable_insert(endNamesHash, pinchEdge->to,
                    cactusMisc_nameToString(end_getName(end)));
        }
    }
    flower_destructEndIterator(endIterator);
    st_logDebug("Built the end names hash\n");

    ////////////////////////////////////////////////
    //Prune the cactus graph to include only those edges relevant to the desired flower.
    ////////////////////////////////////////////////

    mergedVertexIDs
            = st_malloc(sizeof(int32_t) * cactusGraph->vertices->length);
    for (i = 0; i < cactusGraph->vertices->length; i++) {
        mergedVertexIDs[i]
                = ((struct CactusVertex *) cactusGraph->vertices->list[i])->vertexID;
    }

    for (i = 0; i < biConnectedComponents->length; i++) {
        biConnectedComponent = biConnectedComponents->list[i];
        list = constructEmptyList(0, NULL);
        for (j = 0; j < biConnectedComponent->length; j++) {
            cactusEdge = biConnectedComponent->list[j];
            if (isAStubCactusEdge(cactusEdge, pinchGraph)) {
#ifdef BEN_DEBUG
                cactusEdge2 = getNonDeadEndOfStubCactusEdge(cactusEdge,
                        pinchGraph);
                assert(cactusEdge2 != NULL);
                Name name = cactusEdgeToEndName(cactusEdge2, endNamesHash,
                        pinchGraph);
                end = flower_getEnd(parentFlower, name);
                assert(end != NULL);
                assert(end_isStubEnd(end));
                if (end_isFree(end)) { //we don't want free stubs in chains
                    assert(isAFreeStubCactusEdge(cactusEdge, pinchGraph, parentFlower));
                    assert(biConnectedComponent->length == 1);
                } else {
                    assert(end_isAttached(end));
                    assert(!isAFreeStubCactusEdge(cactusEdge, pinchGraph, parentFlower));
                }
#endif
                listAppend(list, cactusEdge);
            } else if (stSortedSet_search(chosenBlocks, cactusEdge) == NULL) { //is a non stub not in the chosen list.
                assert(cactusEdge->pieces->length == 1);
                //merge vertices
                mergeCactusVertices(cactusEdge, mergedVertexIDs, j,
                        biConnectedComponent);
            } else { //is a non stub in the chosen list, so we'll add it.
                listAppend(list, cactusEdge);
            }
        }
        destructList(biConnectedComponent);
        biConnectedComponents->list[i] = list;
    }

    st_logDebug("Built the chosen blocks hash\n");

    ////////////////////////////////////////////////
    //Blocks and ends for each flower.
    ////////////////////////////////////////////////

    flowers = st_malloc(sizeof(void *) * cactusGraph->vertices->length);
    for (i = 1; i < cactusGraph->vertices->length; i++) {
        flowers[i] = NULL;
    }
    flowers[0] = parentFlower;
    parentFlowers = st_malloc(sizeof(void *) * biConnectedComponents->length);
    for (i = 0; i < biConnectedComponents->length; i++) {
        biConnectedComponent = biConnectedComponents->list[i];
        if (biConnectedComponent->length > 0) {
            cactusEdge = biConnectedComponent->list[0];
            //Get the flower.
            flower = flowers[mergedVertexIDs[cactusEdge->from->vertexID]];
            if (flower == NULL) {
                flower = flower_construct(flower_getCactusDisk(parentFlower));
                eventTree_copyConstruct(flower_getEventTree(parentFlower),
                        flower, returnsTrue);
                flowers[mergedVertexIDs[cactusEdge->from->vertexID]] = flower;
            }
            parentFlowers[i] = flower;
            //Make the blocks and ends
            getOrientation(biConnectedComponent, pinchGraph, parentFlower,
                    endNamesHash);
            for (j = 0; j < biConnectedComponent->length; j++) {
                cactusEdge = biConnectedComponent->list[j];
                piece = cactusEdge->pieces->list[0];
                if (!isAStubCactusEdge(cactusEdge, pinchGraph)) {
                    block = constructBlockFromCactusEdge(cactusEdge, flower);
                    pinchEdge = cactusEdgeToFirstPinchEdge(cactusEdge,
                            pinchGraph);
                    hashtable_insert(endNamesHash, pinchEdge->from,
                            cactusMisc_nameToString(end_getName(block_get5End(
                                    block))));
                    hashtable_insert(endNamesHash, pinchEdge->to,
                            cactusMisc_nameToString(end_getName(block_get3End(
                                    block))));
                    assert(cactusEdgeToEndName(cactusEdge, endNamesHash, pinchGraph) == end_getName(block_get5End(block)));
                    assert(cactusEdgeToEndName(cactusEdge->rEdge, endNamesHash, pinchGraph) == end_getName(block_get3End(block)));
                } else {
                    assert(j == 0 || j == biConnectedComponent->length-1);
                    cactusEdge2 = getNonDeadEndOfStubCactusEdge(cactusEdge,
                            pinchGraph);
                    end = flower_getEnd(parentFlower, cactusEdgeToEndName(
                            cactusEdge2, endNamesHash, pinchGraph));
                    assert(end != NULL);
                    if (flower != parentFlower) {
                        end_copyConstruct(end, flower);
                    }
                }
            }
        }
    }
    st_logDebug("Constructed blocks and flowers for the cycle.\n");

    ////////////////////////////////////////////////
    //Link flowers to parent flowers and construct chains.
    ////////////////////////////////////////////////

    for (i = 0; i < biConnectedComponents->length; i++) {
        biConnectedComponent = biConnectedComponents->list[i];
        flower = parentFlowers[i];
        if (biConnectedComponent->length > 1) {
            assert(flower != NULL);
            chain = chain_construct(flower);
            for (j = 1; j < biConnectedComponent->length; j++) {
                cactusEdge = biConnectedComponent->list[j - 1];
                cactusEdge2 = biConnectedComponent->list[j];
                nestedFlower
                        = flowers[mergedVertexIDs[cactusEdge->to->vertexID]];
                assert(cactusEdge->to->vertexID != 0);
                if (nestedFlower == NULL) { //construct a terminal group.
                    group = group_construct2(flower);
                    end
                            = flower_getEnd(
                                    flower,
                                    cactusEdgeToEndName(
                                            isAStubCactusEdge(cactusEdge,
                                                    pinchGraph) ? getNonDeadEndOfStubCactusEdge(
                                                    cactusEdge, pinchGraph)
                                                    : cactusEdge->rEdge,
                                            endNamesHash, pinchGraph));
                    assert(end != NULL);
                    end_setGroup(end, group);

                    end2
                            = flower_getEnd(
                                    flower,
                                    cactusEdgeToEndName(
                                            isAStubCactusEdge(cactusEdge2,
                                                    pinchGraph) ? getNonDeadEndOfStubCactusEdge(
                                                    cactusEdge2, pinchGraph)
                                                    : cactusEdge2,
                                            endNamesHash, pinchGraph));
                    assert(end2 != NULL);
                    end_setGroup(end2, group);
                } else { //construct a link between two existing chains.
                    assert(flower_getEndNumber(nestedFlower)> 0);
                    end
                            = flower_getEnd(
                                    flower,
                                    cactusEdgeToEndName(
                                            isAStubCactusEdge(cactusEdge,
                                                    pinchGraph) ? getNonDeadEndOfStubCactusEdge(
                                                    cactusEdge, pinchGraph)
                                                    : cactusEdge->rEdge,
                                            endNamesHash, pinchGraph));
                    assert(end != NULL);
                    end_copyConstruct(end, nestedFlower);
                    end2
                            = flower_getEnd(
                                    flower,
                                    cactusEdgeToEndName(
                                            isAStubCactusEdge(cactusEdge2,
                                                    pinchGraph) ? getNonDeadEndOfStubCactusEdge(
                                                    cactusEdge2, pinchGraph)
                                                    : cactusEdge2,
                                            endNamesHash, pinchGraph));
                    assert(end2 != NULL);
                    end_copyConstruct(end2, nestedFlower);
                    group = group_construct(flower, nestedFlower);
                    flowers[mergedVertexIDs[cactusEdge->to->vertexID]] = NULL;
                }
                //Make link chain
                link_construct(end, end2, group, chain);
            }
        }
    }
    flower = flowers[0];
#ifdef BEN_DEBUG
    for (i = 1; i < cactusGraph->vertices->length; i++) {
        assert(flowers[i] == NULL);
    }
#endif
    st_logDebug("Constructed the chains and linked together the flowers\n");

    ////////////////////////////////////////////////
    //Add nested ends to flowers.
    ////////////////////////////////////////////////

    destructList(addEnvelopedStubEnds(parentFlower, 0));
    st_logDebug("Added the nested ends to the parent flowers\n");

    ////////////////////////////////////////////////
    //Add adjacencies between ends.
    ////////////////////////////////////////////////

    addAdjacenciesToEnds(flower);
    st_logDebug("Added the adjacencies between the ends\n");

    ////////////////////////////////////////////////
    //Add groups.
    ////////////////////////////////////////////////

    addGroups(flower);
    st_logDebug("Added the tangle groups\n");

    ////////////////////////////////////////////////
    //Set blocks for each flower to 'built'
    ////////////////////////////////////////////////

    assert(flower == parentFlower);
    setBlocksBuilt(flower);

    ////////////////////////////////////////////////
    //Clean up
    ////////////////////////////////////////////////

    free(flowers);
    free(mergedVertexIDs);
    free(vertexDiscoveryTimes);
    free(parentFlowers);
    destructList(biConnectedComponents);
    hashtable_destroy(endNamesHash, TRUE, FALSE);
}

void copyEndTreePhylogenies(Flower *parentFlower, Flower *flower) {
    /*
     * For each end in the parent flower that is also in the flower, we copy
     * the the phylogefloweric information across.
     */
    End *end1;
    End *end2;
    Cap *cap1;
    Cap *cap2;
    Cap *cap3;
    Cap *cap4;
    Flower_EndIterator *endIterator;
    End_InstanceIterator *instanceIterator;

    endIterator = flower_getEndIterator(flower);
    while ((end1 = flower_getNextEnd(endIterator)) != NULL) {
        end2 = flower_getEnd(parentFlower, end_getName(end1));
        assert(end2 != NULL);
        instanceIterator = end_getInstanceIterator(end1);
        while ((cap1 = end_getNext(instanceIterator)) != NULL) {
            assert(cap_getParent(cap1) == NULL);
            cap2 = end_getInstance(end2, cap_getName(cap1));
            assert(cap2 != NULL);
            if ((cap3 = cap_getParent(cap2)) != NULL) {
                cap4 = end_getInstance(end1, cap_getName(cap3));
                assert(cap4 != NULL);
                cap_makeParentAndChild(cap4, cap1);
            } else {
                assert(end_getRootInstance(end2) == cap2);
            }
        }
        end_destructInstanceIterator(instanceIterator);
    }
    flower_destructEndIterator(endIterator);
}

