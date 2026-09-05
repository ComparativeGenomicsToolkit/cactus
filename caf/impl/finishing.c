#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

///////////////////////////////////////////////////////////////////////////
// Convert the complete cactus graph/pinch graph into filled out set of flowers
///////////////////////////////////////////////////////////////////////////

// Functions used to record the flower ends on the pinch ends. Every pinch end's data slot holds the End (of
// some flower) it corresponds to, or NULL where that end's block has not yet been made into a flower block.

static void setPinchBlockEndsToEndsPP(stPinchBlock *pinchBlock, bool orientation, End *end) {
    stPinchEnd *pinchEnd = stPinchBlock_getEnd(pinchBlock, orientation);
    assert(pinchEnd != NULL);
    if (stPinchEnd_getData(pinchEnd) == NULL) {
        stPinchEnd_setData(pinchEnd, end);
    } else {
        assert(stPinchEnd_getData(pinchEnd) == end);
    }
}

static void setPinchBlockEndsToEndsP(stPinchSegment *pinchSegment, bool endOrientation, Cap *cap) {
    stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
    assert(pinchBlock != NULL);
    assert(cap != NULL);
    End *end = end_getPositiveOrientation(cap_getEnd(cap));
    assert(end != NULL);
    assert(!end_isBlockEnd(end));
    assert(end_getOrientation(end));
    assert(!end_getOrientation(end_getReverse(end)));
    setPinchBlockEndsToEndsPP(pinchBlock, endOrientation, end_getReverse(end));
    setPinchBlockEndsToEndsPP(pinchBlock, !endOrientation, end);
}

static void setPinchEndsToEnds(stPinchThreadSet *threadSet, Flower *parentFlower) {
    stPinchThreadSet_clearEndData(threadSet); //the slots held the cactus nodes while the graph was built
    stPinchThreadSetIt pinchThreadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *pinchThread;
    while ((pinchThread = stPinchThreadSetIt_getNext(&pinchThreadIt))) {
        Cap *cap = flower_getCap(parentFlower, stPinchThread_getName(pinchThread));
        assert(cap != NULL);
        stPinchSegment *pinchSegment = stPinchThread_getFirst(pinchThread);
        setPinchBlockEndsToEndsP(pinchSegment, stPinchSegment_getBlockOrientation(pinchSegment), cap);
        pinchSegment = stPinchThread_getLast(pinchThread);
        setPinchBlockEndsToEndsP(pinchSegment, !stPinchSegment_getBlockOrientation(pinchSegment), cap_getAdjacency(cap));
    }
}

// Functions used to build the empty flower hierarchy

static stList *makeEmptyFlowers(stCactusNode *cactusNode, Flower *flower, stPinchThreadSet *threadSet,
                                stHash *cactusNodesToFlowers, bool isRoot);

static void makeEmptyFlowers2(stCactusEdgeEnd *cactusEdgeEnd, Flower *flower, stPinchThreadSet *threadSet,
                              stHash *cactusNodesToFlowers, stList *containedStubEnds) {
    // Iterate around a chain in the flower creating the nested flowers
    cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd);
    if (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd)) { //We have a non-trivial chain
        do {
            stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd);

            // Get the cactus node for the link in the chain
            stCactusNode *cactusNode = stCactusEdgeEnd_getNode(cactusEdgeEnd);
            assert(cactusNode == stCactusEdgeEnd_getNode(linkedCactusEdgeEnd));

            // Make a new group
            Group *group = group_construct2(flower);

            //Make a nested flower
            Flower *nestedFlower = group_makeEmptyNestedFlower(group);

            // Add to hash of cactusNodes to nested flowers
            assert(stHash_search(cactusNodesToFlowers, cactusNode) == NULL);
            stHash_insert(cactusNodesToFlowers, cactusNode, nestedFlower);

            //Fill out stack
            stList *nestedStubEnds =
                    makeEmptyFlowers(cactusNode, nestedFlower, threadSet, cactusNodesToFlowers, 0);
            assert(nestedStubEnds != NULL);

            // Add all the stub ends in the nested flower
            stList_appendAll(containedStubEnds, nestedStubEnds);
            stList_destruct(nestedStubEnds);

            // Move to the next link in the chain
            cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(linkedCactusEdgeEnd);
        } while (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd));
    }
}

static int sort_pinch_ends(const void *a, const void *b) {
    End *end1 = stPinchEnd_getData((const stPinchEnd *)a);
    End *end2 = stPinchEnd_getData((const stPinchEnd *)b);
    assert(end1 != NULL && end2 != NULL);
    return (int)cactusMisc_nameCompare(end_getName(end1), end_getName(end2));
}

static stList *makeEmptyFlowers(stCactusNode *cactusNode, Flower *flower, stPinchThreadSet *threadSet,
                                stHash *cactusNodesToFlowers, bool isRoot) {
    stList *containedStubEnds = stList_construct(); // The stub ends that are contained in this flower or nested
    // versions of it

    // Iterate through each incident chain
    stCactusNodeEdgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
    stCactusEdgeEnd *cactusEdgeEnd;
    while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt))) { // For each end
        // If is a chain end and in the right orientation
        if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) {
            makeEmptyFlowers2(cactusEdgeEnd, flower, threadSet,
                              cactusNodesToFlowers, containedStubEnds); // Make the empty flowers
        }
    }

    if(!isRoot) {
        assert(flower_getEndNumber(flower) == 0);

        // Add in the stub ends that are contained in this flower
        stList *adjacencyComponents = stCactusNode_getObject(cactusNode);
        for (int64_t i = 0; i < stList_length(adjacencyComponents); i++) {
            stList *adjacencyComponent = stList_get(adjacencyComponents, i);
            for(int64_t j=0; j<stList_length(adjacencyComponent); j++) {
                stPinchEnd *pinchEnd = stList_get(adjacencyComponent, j);
                End *end = stPinchEnd_getData(pinchEnd);
                if(end != NULL) {
                    // Deal with components for dead ends of free stubs
                    if (stList_length(adjacencyComponent) == 1 && !end_getOrientation(end)) {
                        continue;
                    }
                    stList_append(containedStubEnds, pinchEnd);
                }
            }
        }

        // Sort the stub ends so they are added in order
        stList_sort(containedStubEnds, sort_pinch_ends);

        // Add the stub ends into the flower
        for (int64_t i = 0; i < stList_length(containedStubEnds); i++) {
            stPinchEnd *pinchEnd = stList_get(containedStubEnds, i);
            End *end = stPinchEnd_getData(pinchEnd);
            assert(end != NULL);
            assert(flower_getEnd(flower, end_getName(end)) == NULL);
            End *end2 = end_copyConstruct(end_getPositiveOrientation(end), flower);

            // The end now maps to its copy in this flower
            stPinchEnd_setData(pinchEnd, end2);

            // Set the group
            assert(end_getFlower(end) != NULL);
            Group *group = flower_getParentGroup(end_getFlower(end));
            if(group != NULL && group_getFlower(group) == flower) {
                end_setGroup(end2, group);
            }
        }

        return containedStubEnds;
    }
    else { // At the root of this hierarchy, so stub ends are already present
        for (int64_t i = 0; i < stList_length(containedStubEnds); i++) {
            stPinchEnd *pinchEnd = stList_get(containedStubEnds, i);
            End *end = stPinchEnd_getData(pinchEnd);
            assert(end != NULL);
            End *end2 = flower_getEnd(flower, end_getName(end));
            assert(end2 != NULL);

            // Set the group
            assert(end_getFlower(end) != NULL);
            Group *group = flower_getParentGroup(end_getFlower(end));
            if(group != NULL && group_getFlower(group) == flower) {
                end_setGroup(end2, group);
            }
        }

        stList_destruct(containedStubEnds);
        return NULL;
    }
}

//Functions for going from cactus/pinch ends to flower ends and updating flower structure as necessary

static End *convertPinchBlockEndToEnd(stPinchEnd *pinchEnd, Flower *flower) {
    End *end = stPinchEnd_getData(pinchEnd);
    if (end == NULL) { //Happens if pinch end represents end of a block in flower that has not yet been defined.
        return NULL;
    }
    End *end2 = flower_getEnd(flower, end_getName(end));
    assert(end2 != NULL);
    assert(end_getOrientation(end2));
    return end_getOrientation(end) ? end2 : end_getReverse(end2);
}

static End *convertCactusEdgeEndToEnd(stCactusEdgeEnd *cactusEdgeEnd, Flower *flower) {
    return convertPinchBlockEndToEnd(stCactusEdgeEnd_getObject(cactusEdgeEnd), flower);
}

//Functions to create blocks

static void makeBlockP(stPinchEnd *pinchEnd, End *end) {
    assert(stPinchEnd_getData(pinchEnd) == NULL);
    stPinchEnd_setData(pinchEnd, end);
}

static void makeBlock(stCactusEdgeEnd *cactusEdgeEnd, Flower *parentFlower, Flower *flower) {
    stPinchEnd *pinchEnd = stCactusEdgeEnd_getObject(cactusEdgeEnd);
    assert(pinchEnd != NULL);
    stPinchBlock *pinchBlock = stPinchEnd_getBlock(pinchEnd);
    Block *block = block_construct(stPinchBlock_getLength(pinchBlock), flower);
    stPinchSegment *pinchSegment;
    stPinchBlockIt pinchSegmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    while ((pinchSegment = stPinchBlockIt_getNext(&pinchSegmentIt))) {
        Cap *parentCap = stCaf_getThreadCap(pinchSegment, parentFlower); //The following three lines isolates the sequence associated with a segment.
        assert(parentCap != NULL);
        Sequence *parentSequence = cap_getSequence(parentCap);
        assert(parentSequence != NULL);
        Sequence *sequence = flower_getSequence(flower, sequence_getName(parentSequence));
        if (sequence == NULL) {
            sequence = cactusDisk_getSequence(flower_getCactusDisk(flower), sequence_getName(parentSequence));
            flower_addSequence(flower, sequence);
        }
        assert(sequence != NULL);
        segment_construct2(
                stPinchEnd_getOrientation(pinchEnd) ^ stPinchSegment_getBlockOrientation(pinchSegment) ? block_getReverse(block) : block,
                stPinchSegment_getStart(pinchSegment), 1, sequence);
    }
    makeBlockP(pinchEnd, block_get5End(block));
    stPinchEnd *otherPinchBlockEnd = stCactusEdgeEnd_getObject(stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd));
    makeBlockP(otherPinchBlockEnd, block_get3End(block));
}

//Functions to generate the chains of a flower

static void fillOutFlowers(stCactusNode *cactusNode, Flower *flower, bool orientation, stPinchThreadSet *threadSet,
                           Flower *parentFlower, stList *deadEndComponent,
                           stHash *cactusNodesToFlowers);

static void fillOutChain(stCactusEdgeEnd *cactusEdgeEnd, Flower *flower, bool orientation,
                         stPinchThreadSet *threadSet,  Flower *parentFlower, stList *deadEndComponent,
                         stHash *cactusNodesToFlowers, bool fillOutNestedFlowers) {
    cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd);
    if (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd)) { //We have a non-trivial chain
        Chain *chain = fillOutNestedFlowers ? chain_construct(flower) : NULL;
        do {
            stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd);
            if (convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, flower) == NULL) { //Make subsequent block
                makeBlock(linkedCactusEdgeEnd, parentFlower, flower);
            }

            if(fillOutNestedFlowers) {
                stCactusNode *cactusNode = stCactusEdgeEnd_getNode(cactusEdgeEnd);
                Flower *nestedFlower = stHash_search(cactusNodesToFlowers, cactusNode);
                assert(nestedFlower != NULL);

                assert(cactusNode == stCactusEdgeEnd_getNode(linkedCactusEdgeEnd));
                Group *group = flower_getParentGroup(nestedFlower);
                assert(group != NULL);
                End *end1 = convertCactusEdgeEndToEnd(cactusEdgeEnd, flower);
                End *end2 = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, flower);
                assert(end1 != NULL);
                assert(end2 != NULL);
                assert(end_getOrientation(end1));
                assert(end_getOrientation(end2));
                assert(!end_getSide(end1));
                assert(end_getSide(end2));
                assert(end_isBlockEnd(end1) || end_isAttached(end1));
                assert(end_isBlockEnd(end2) || end_isAttached(end2));
                if (end_getGroup(end1) == NULL) {
                    end_setGroup(end1, group);
                }
                assert(end_getGroup(end1) == group);
                if (end_getGroup(end2) == NULL) {
                    end_setGroup(end2, group);
                }
                assert(end_getGroup(end2) == group);
                link_construct(end1, end2, group, chain);

                //Add flowers to nested
                if (end_getName(end1) > end_getName(end2)) { // Swap order so added
                    End *e = end1;
                    end1 = end2;
                    end2 = e;
                }
                if (flower_getEnd(nestedFlower, end_getName(end1)) == NULL) {
                    end_copyConstruct(end1, nestedFlower);
                }
                if (flower_getEnd(nestedFlower, end_getName(end2)) == NULL) {
                    end_copyConstruct(end2, nestedFlower);
                }

                //Fill out stack
                fillOutFlowers(cactusNode, nestedFlower, orientation, threadSet,
                               parentFlower, deadEndComponent, cactusNodesToFlowers);
            }

            cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(linkedCactusEdgeEnd);
        } while (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd));
    }
}

static void fillOutChains(stCactusNode *cactusNode, Flower *flower, bool orientation,
                          stPinchThreadSet *threadSet,  Flower *parentFlower,
                          stList *deadEndComponent, stHash *cactusNodesToFlowers, bool fillOutNestedFlowers) {
    stCactusNodeEdgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
    stCactusEdgeEnd *cactusEdgeEnd;
    while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt))) {
        if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) { //We have some sort of chain
            End *end = convertCactusEdgeEndToEnd(cactusEdgeEnd, flower);
            stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd), *startCactusEdgeEnd = NULL;
            assert(linkedCactusEdgeEnd != NULL);
            bool orientation2;
            if (end != NULL) {
#ifndef NDEBUG
                End *end2;
                if ((end2 = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, flower)) != NULL) {
                    assert(end_getSide(end) != end_getSide(end2));
                }
#endif
                startCactusEdgeEnd = end_getSide(end) ? cactusEdgeEnd : linkedCactusEdgeEnd;
                orientation2 = end_getSide(end);
            } else {
                end = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, flower);
                if (end != NULL) {
                    if (end_getSide(end)) {
                        startCactusEdgeEnd = linkedCactusEdgeEnd;
                    } else {
                        makeBlock(cactusEdgeEnd, parentFlower, flower);
                        startCactusEdgeEnd = cactusEdgeEnd;
                    }
                    orientation2 = !end_getSide(end);
                } else {
                    if(orientation) {
                        makeBlock(cactusEdgeEnd, parentFlower, flower);
                        startCactusEdgeEnd = cactusEdgeEnd;
                    }
                    else {
                        makeBlock(linkedCactusEdgeEnd, parentFlower, flower);
                        startCactusEdgeEnd = linkedCactusEdgeEnd;
                    }
                    orientation2 = orientation;
                }
            }
            assert(startCactusEdgeEnd != NULL);
            fillOutChain(startCactusEdgeEnd, flower, orientation2, threadSet, parentFlower,
                         deadEndComponent, cactusNodesToFlowers, fillOutNestedFlowers);
            //fillOutChain(startCactusEdgeEnd, flower, orientation2, threadSet, parentFlower,
            //             deadEndComponent, cactusNodesToFlowers, 1);
        }
    }
}

/*
 * Adds in groups for the tangles (groups not contained as a link in a chain) in the flower.
 */
static void makeTangles(stCactusNode *cactusNode, Flower *flower, stList *deadEndComponent) {
    stList *adjacencyComponents = stCactusNode_getObject(cactusNode);
    for (int64_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        if (adjacencyComponent != deadEndComponent) {
            if (stList_length(adjacencyComponent) == 1) { //Deal with components for dead ends of free stubs
                End *end = convertPinchBlockEndToEnd(stList_get(adjacencyComponent, 0), flower);
                assert(end != NULL);
                if (!end_getOrientation(end)) {
                    continue;
                }
            }
            Group *group = group_construct2(flower);
            for (int64_t j = 0; j < stList_length(adjacencyComponent); j++) {
                End *end = convertPinchBlockEndToEnd(stList_get(adjacencyComponent, j), flower);
                assert(end != NULL);
                assert(end_getOrientation(end));
                assert(end_getGroup(end) == NULL);
                end_setGroup(end, group);
            }
        }
    }
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while((group = flower_getNextGroup(groupIt)) != NULL) {
        group_constructChainForLink(group);
    }
    flower_destructGroupIterator(groupIt);
}

/*
 * Adds in the chains and completes the groups for the flower and its nested flowers, recursively.
 */
static void fillOutFlowers(stCactusNode *cactusNode, Flower *flower, bool orientation, stPinchThreadSet *threadSet,
                           Flower *parentFlower, stList *deadEndComponent, stHash *cactusNodesToFlowers) {
    assert(flower_getAttachedStubEndNumber(flower) > 0);
    fillOutChains(cactusNode, flower, orientation, threadSet, parentFlower, deadEndComponent,
                  cactusNodesToFlowers, 0);
    fillOutChains(cactusNode, flower, orientation, threadSet, parentFlower, deadEndComponent,
                  cactusNodesToFlowers, 1); //This call is recursive
    makeTangles(cactusNode, flower, deadEndComponent);
    stCaf_addAdjacencies(flower);
    if(flower_isLeaf(flower) && flower_getBlockNumber(flower) == 0 && flower != parentFlower) { //We have a leaf with no blocks - it's effectively empty and can be removed.
        flower_destruct(flower, 1, 1); //This removes the flower completely from the database.
    }
    else {
        flower_setBuiltBlocks(flower, 1);
    }
}

//Main function

static void stCaf_convertCactusGraphToFlowers(stPinchThreadSet *threadSet, stCactusNode *startCactusNode,
                                              Flower *parentFlower, stList *deadEndComponent) {
    double t = stCaf_now();
    setPinchEndsToEnds(threadSet, parentFlower);
    stHash *cactusNodesToFlowers = stHash_construct();
    double endsHashTime = stCaf_now() - t;
    t = stCaf_now();
    makeEmptyFlowers(startCactusNode, parentFlower, threadSet, cactusNodesToFlowers, 1);
    double emptyFlowersTime = stCaf_now() - t;
    t = stCaf_now();
    fillOutFlowers(startCactusNode, parentFlower, 1, threadSet, parentFlower, deadEndComponent,
                   cactusNodesToFlowers);
    double fillOutTime = stCaf_now() - t;
    t = stCaf_now();
    stHash_destruct(cactusNodesToFlowers);
    st_logInfo("caf-timing: convert endsHash %.3fs emptyFlowers %.3fs fillOut %.3fs cleanup %.3fs\n",
               endsHashTime, emptyFlowersTime, fillOutTime, stCaf_now() - t);
}

///////////////////////////////////////////////////////////////////////////
// Functions to ensure every segment has a block
///////////////////////////////////////////////////////////////////////////

void stCaf_makeDegreeOneBlocks(stPinchThreadSet *threadSet) {
    stPinchThreadSetSegmentIt segmentIt = stPinchThreadSet_getSegmentIt(threadSet);
    stPinchSegment *segment;
    while ((segment = stPinchThreadSetSegmentIt_getNext(&segmentIt)) != NULL) {
        if (stPinchSegment_getBlock(segment) == NULL) {
            stPinchBlock_construct2(segment);
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// Functions for actually filling out cactus
///////////////////////////////////////////////////////////////////////////


void stCaf_finish(Flower *flower, stPinchThreadSet *threadSet, int64_t minLengthForChromosome,
                  double proportionOfUnalignedBasesForNewChromosome) {
    stCactusNode *startCactusNode;
    stList *deadEndComponent;
    double t = stCaf_now();
    stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 1, minLengthForChromosome,
                                                                  proportionOfUnalignedBasesForNewChromosome, 0, INT64_MAX);
    double graphTime = stCaf_now() - t;
    t = stCaf_now();

    //Convert cactus graph/pinch graph to API
    stCaf_convertCactusGraphToFlowers(threadSet, startCactusNode, flower, deadEndComponent);
    double convertTime = stCaf_now() - t;
    t = stCaf_now();

    //Cleanup
    stCaf_destructCactusGraph(cactusGraph, threadSet);
    st_logInfo("caf-timing: finish graph %.3fs convert %.3fs destruct %.3fs\n", graphTime, convertTime, stCaf_now() - t);
}
