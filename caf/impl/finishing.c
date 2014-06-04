#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

///////////////////////////////////////////////////////////////////////////
// Convert the complete cactus graph/pinch graph into filled out set of flowers
///////////////////////////////////////////////////////////////////////////

//Functions used to build hash between pinchEnds and flower ends.

static void getPinchBlockEndsToEndsHashPP(stPinchBlock *pinchBlock, bool orientation, End *end, stHash *pinchEndsToEnds) {
    stPinchEnd pinchEnd = stPinchEnd_constructStatic(pinchBlock, orientation);
    if (stHash_search(pinchEndsToEnds, &pinchEnd) == NULL) {
        stHash_insert(pinchEndsToEnds, stPinchEnd_construct(pinchBlock, orientation), end);
    } else {
        assert(stHash_search(pinchEndsToEnds, &pinchEnd) == end);
    }
}

static void getPinchBlockEndsToEndsHashP(stPinchSegment *pinchSegment, bool endOrientation, Cap *cap, stHash *pinchEndsToEnds) {
    stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
    assert(pinchBlock != NULL);
    assert(cap != NULL);
    End *end = end_getPositiveOrientation(cap_getEnd(cap));
    assert(end != NULL);
    assert(!end_isBlockEnd(end));
    assert(end_getOrientation(end));
    assert(!end_getOrientation(end_getReverse(end)));
    getPinchBlockEndsToEndsHashPP(pinchBlock, endOrientation, end_getReverse(end), pinchEndsToEnds);
    getPinchBlockEndsToEndsHashPP(pinchBlock, !endOrientation, end, pinchEndsToEnds);
}

static stHash *getPinchEndsToEndsHash(stPinchThreadSet *threadSet, Flower *parentFlower) {
    stHash *pinchEndsToEnds = stHash_construct3(stPinchEnd_hashFn, stPinchEnd_equalsFn, (void(*)(void *)) stPinchEnd_destruct, NULL);
    stPinchThreadSetIt pinchThreadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *pinchThread;
    while ((pinchThread = stPinchThreadSetIt_getNext(&pinchThreadIt))) {
        Cap *cap = flower_getCap(parentFlower, stPinchThread_getName(pinchThread));
        assert(cap != NULL);
        stPinchSegment *pinchSegment = stPinchThread_getFirst(pinchThread);
        getPinchBlockEndsToEndsHashP(pinchSegment, stPinchSegment_getBlockOrientation(pinchSegment), cap, pinchEndsToEnds);
        pinchSegment = stPinchThread_getLast(pinchThread);
        getPinchBlockEndsToEndsHashP(pinchSegment, !stPinchSegment_getBlockOrientation(pinchSegment), cap_getAdjacency(cap),
                pinchEndsToEnds);
    }
    return pinchEndsToEnds;
}

//Functions for going from cactus/pinch ends to flower ends and updating flower structure as necessary

static End *convertPinchBlockEndToEnd(stPinchEnd *pinchEnd, stHash *pinchEndsToEnds, Flower *flower) {
    End *end = stHash_search(pinchEndsToEnds, pinchEnd);
    if (end == NULL) { //Happens if pinch end represents end of a block in flower that has not yet been defined.
        return NULL;
    }
    End *end2 = flower_getEnd(flower, end_getName(end));
    if (end2 == NULL) { //Happens if is free stub end not yet defined in given flower - but must be present somewhere in the hierarchy.
        assert(end_isFree(end));
        assert(end_isStubEnd(end));
        //Copy the end down the hierarchy
        Group *parentGroup = flower_getParentGroup(flower);
        assert(parentGroup != NULL); //Can not be at the top of the hierarchy, else would be defined.
        end2 = convertPinchBlockEndToEnd(pinchEnd, pinchEndsToEnds, group_getFlower(parentGroup));
        assert(end2 != NULL);
        assert(end_getGroup(end2) == NULL); //This happens because the group of the end is only defined at the given flower.
        end_setGroup(end2, parentGroup);
        end2 = end_copyConstruct(end_getPositiveOrientation(end2), flower);
        assert(end2 != NULL);
        assert(end_getFlower(end2) == flower);
    }
    assert(end_getOrientation(end2));
    return end_getOrientation(end) ? end2 : end_getReverse(end2);
}

static End *convertCactusEdgeEndToEnd(stCactusEdgeEnd *cactusEdgeEnd, stHash *pinchEndsToEnds, Flower *flower) {
    return convertPinchBlockEndToEnd(stCactusEdgeEnd_getObject(cactusEdgeEnd), pinchEndsToEnds, flower);
}

//Functions to create blocks

static void makeBlockP(stPinchEnd *pinchEnd, End *end, stHash *pinchEndsToEnds) {
    assert(stHash_search(pinchEndsToEnds, pinchEnd) == NULL);
    stHash_insert(pinchEndsToEnds, stPinchEnd_construct(stPinchEnd_getBlock(pinchEnd), stPinchEnd_getOrientation(pinchEnd)), end);
}

static void makeBlock(stCactusEdgeEnd *cactusEdgeEnd, Flower *parentFlower, Flower *flower, stHash *pinchEndsToEnds) {
    stPinchEnd *pinchEnd = stCactusEdgeEnd_getObject(cactusEdgeEnd);
    assert(pinchEnd != NULL);
    stPinchBlock *pinchBlock = stPinchEnd_getBlock(pinchEnd);
    Block *block = block_construct(stPinchBlock_getLength(pinchBlock), flower);
    stPinchSegment *pinchSegment;
    stPinchBlockIt pinchSegmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    while ((pinchSegment = stPinchBlockIt_getNext(&pinchSegmentIt))) {
        Cap *parentCap = flower_getCap(parentFlower, stPinchSegment_getName(pinchSegment)); //The following three lines isolates the sequence associated with a segment.
        assert(parentCap != NULL);
        Sequence *parentSequence = cap_getSequence(parentCap);
        assert(parentSequence != NULL);
        Sequence *sequence = flower_getSequence(flower, sequence_getName(parentSequence));
        if (sequence == NULL) {
            sequence = sequence_construct(cactusDisk_getMetaSequence(flower_getCactusDisk(flower), sequence_getName(parentSequence)),
                    flower);
        }
        assert(sequence != NULL);
        segment_construct2(
                stPinchEnd_getOrientation(pinchEnd) ^ stPinchSegment_getBlockOrientation(pinchSegment) ? block_getReverse(block) : block,
                stPinchSegment_getStart(pinchSegment), 1, sequence);
    }
    makeBlockP(pinchEnd, block_get5End(block), pinchEndsToEnds);
    stPinchEnd *otherPinchBlockEnd = stCactusEdgeEnd_getObject(stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd));
    makeBlockP(otherPinchBlockEnd, block_get3End(block), pinchEndsToEnds);
}

//Functions to generate the chains of a flower

static void makeFlower(stCactusNode *cactusNode, Flower *flower, bool orientation,
        stPinchThreadSet *threadSet,  Flower *parentFlower,
        stList *deadEndComponent, stSet *bigFlowers, stHash *pinchEndsToEnds);

static void makeChain(stCactusEdgeEnd *cactusEdgeEnd, Flower *flower, bool orientation, bool makeSpacerFlowers,
        stPinchThreadSet *threadSet,  Flower *parentFlower,
                stList *deadEndComponent, stSet *bigFlowers, stHash *pinchEndsToEnds) {
    cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(cactusEdgeEnd);
    if (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd)) { //We have a non-trivial chain
        Chain *chain = chain_construct(flower);
        do {
            stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd);
            if (convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower) == NULL) { //Make subsequent block
                makeBlock(linkedCactusEdgeEnd, parentFlower, flower, pinchEndsToEnds);
            }
            assert(stCactusEdgeEnd_getNode(cactusEdgeEnd) == stCactusEdgeEnd_getNode(linkedCactusEdgeEnd));
            Group *group = group_construct2(flower);
            End *end1 = convertCactusEdgeEndToEnd(cactusEdgeEnd, pinchEndsToEnds, flower);
            End *end2 = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower);
            assert(end1 != NULL);
            assert(end2 != NULL);
            assert(end_getOrientation(end1));
            assert(end_getOrientation(end2));
            assert(!end_getSide(end1));
            assert(end_getSide(end2));
            assert(end_isBlockEnd(end1) || end_isAttached(end1));
            assert(end_isBlockEnd(end2) || end_isAttached(end2));
            end_setGroup(end1, group);
            end_setGroup(end2, group);
            link_construct(end1, end2, group, chain);
            //Make a nested group
            Flower *nestedFlower = group_makeEmptyNestedFlower(group);
            end_copyConstruct(end1, nestedFlower);
            end_copyConstruct(end2, nestedFlower);
            assert(flower_getGroupNumber(nestedFlower) == 0);

            //If we need the extra spacer flower
            Flower *spacerFlower;
            if(makeSpacerFlowers) {
                Group *nestedGroup = group_construct2(nestedFlower);
                End *nestedEnd1 = flower_getEnd(nestedFlower, end_getName(end1));
                End *nestedEnd2 = flower_getEnd(nestedFlower, end_getName(end2));
                end_setGroup(nestedEnd1, nestedGroup);
                end_setGroup(nestedEnd2, nestedGroup);
                Chain *nestedChain = chain_construct(nestedFlower);
                link_construct(nestedEnd1, nestedEnd2, nestedGroup, nestedChain);
                spacerFlower = nestedFlower;
                nestedFlower = group_makeEmptyNestedFlower(nestedGroup);
                end_copyConstruct(nestedEnd1, nestedFlower);
                end_copyConstruct(nestedEnd2, nestedFlower);
            }
            //Fill out stack
            makeFlower(stCactusEdgeEnd_getNode(cactusEdgeEnd), nestedFlower, orientation, threadSet, parentFlower, deadEndComponent, bigFlowers, pinchEndsToEnds);

            if(makeSpacerFlowers) { //Cleanup memory of spacer flowers
                stCaf_addAdjacencies(spacerFlower);
                flower_setBuiltBlocks(spacerFlower, 1);
                cactusDisk_addUpdateRequest(flower_getCactusDisk(spacerFlower), spacerFlower);
                flower_unload(spacerFlower);
            }

            //stList_append(stack, stCactusEdgeEnd_getNode(cactusEdgeEnd));
            //stList_append(stack, nestedFlower);
            //stList_append(stack, stIntTuple_construct1( orientation));
            cactusEdgeEnd = stCactusEdgeEnd_getOtherEdgeEnd(linkedCactusEdgeEnd);
        } while (!stCactusEdgeEnd_isChainEnd(cactusEdgeEnd));
    }
}

static void makeChains(stCactusNode *cactusNode, Flower *flower, bool orientation,
        stPinchThreadSet *threadSet,  Flower *parentFlower,
        stList *deadEndComponent, stSet *bigFlowers, stHash *pinchEndsToEnds) {
    bool makeSpacerFlowers = stSet_search(bigFlowers, cactusNode) != NULL;
    stCactusNodeEdgeEndIt cactusEdgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
    stCactusEdgeEnd *cactusEdgeEnd;
    while ((cactusEdgeEnd = stCactusNodeEdgeEndIt_getNext(&cactusEdgeEndIt))) {
        if (stCactusEdgeEnd_isChainEnd(cactusEdgeEnd) && stCactusEdgeEnd_getLinkOrientation(cactusEdgeEnd)) { //We have some sort of chain
            End *end = convertCactusEdgeEndToEnd(cactusEdgeEnd, pinchEndsToEnds, flower);
            stCactusEdgeEnd *linkedCactusEdgeEnd = stCactusEdgeEnd_getLink(cactusEdgeEnd), *startCactusEdgeEnd = NULL;
            assert(linkedCactusEdgeEnd != NULL);
            bool orientation2;
            if (end != NULL) {
#ifndef NDEBUG
                End *end2;
                if ((end2 = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower)) != NULL) {
                    assert(end_getSide(end) != end_getSide(end2));
                }
#endif
                startCactusEdgeEnd = end_getSide(end) ? cactusEdgeEnd : linkedCactusEdgeEnd;
                orientation2 = end_getSide(end);
            } else {
                end = convertCactusEdgeEndToEnd(linkedCactusEdgeEnd, pinchEndsToEnds, flower);
                if (end != NULL) {
                    if (end_getSide(end)) {
                        startCactusEdgeEnd = linkedCactusEdgeEnd;
                    } else {
                        makeBlock(cactusEdgeEnd, parentFlower, flower, pinchEndsToEnds);
                        startCactusEdgeEnd = cactusEdgeEnd;
                    }
                    orientation2 = !end_getSide(end);
                } else {
                    if(orientation) {
                        makeBlock(cactusEdgeEnd, parentFlower, flower, pinchEndsToEnds);
                        startCactusEdgeEnd = cactusEdgeEnd;
                    }
                    else {
                        makeBlock(linkedCactusEdgeEnd, parentFlower, flower, pinchEndsToEnds);
                        startCactusEdgeEnd = linkedCactusEdgeEnd;
                    }
                    orientation2 = orientation;
                }
            }
            assert(startCactusEdgeEnd != NULL);
            makeChain(startCactusEdgeEnd, flower, orientation2, makeSpacerFlowers, threadSet, parentFlower, deadEndComponent, bigFlowers, pinchEndsToEnds);
        }
    }
}

//Functions to make tangles.

static void makeTangles(stCactusNode *cactusNode, Flower *flower, stHash *pinchEndsToEnds, stList *deadEndComponent) {
    stList *adjacencyComponents = stCactusNode_getObject(cactusNode);
    for (int64_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        if (adjacencyComponent != deadEndComponent) {
            if (stList_length(adjacencyComponent) == 1) { //Deal with components for dead ends of free stubs
                End *end = convertPinchBlockEndToEnd(stList_get(adjacencyComponent, 0), pinchEndsToEnds, flower);
                assert(end != NULL);
                if (!end_getOrientation(end)) {
                    continue;
                }
            }
            Group *group = group_construct2(flower);
            for (int64_t j = 0; j < stList_length(adjacencyComponent); j++) {
                End *end = convertPinchBlockEndToEnd(stList_get(adjacencyComponent, j), pinchEndsToEnds, flower);
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

//Main function

static void makeFlower(stCactusNode *cactusNode, Flower *flower, bool orientation,
        stPinchThreadSet *threadSet,  Flower *parentFlower,
        stList *deadEndComponent, stSet *bigFlowers, stHash *pinchEndsToEnds) {
    assert(flower_getAttachedStubEndNumber(flower) > 0);
    makeChains(cactusNode, flower, orientation, threadSet, parentFlower, deadEndComponent, bigFlowers, pinchEndsToEnds); //This call is recursive
    makeTangles(cactusNode, flower, pinchEndsToEnds, deadEndComponent);
    stCaf_addAdjacencies(flower);
    if(flower_isLeaf(flower) && flower_getBlockNumber(flower) == 0 && flower != parentFlower) { //We have a leaf with no blocks - it's effectively empty and can be removed.
        flower_delete2(flower, 0); //This removes the flower completely from the database.
    }
    else {
        flower_setBuiltBlocks(flower, 1);
        //Now serialise for writing, compress and free the memory of the flower
        cactusDisk_addUpdateRequest(flower_getCactusDisk(flower), flower);
        flower_unload(flower);
    }
}

static void stCaf_convertCactusGraphToFlowers(stPinchThreadSet *threadSet, stCactusNode *startCactusNode, Flower *parentFlower,
        stList *deadEndComponent, stSet *bigFlowers) {
    stHash *pinchEndsToEnds = getPinchEndsToEndsHash(threadSet, parentFlower);
    makeFlower(startCactusNode, parentFlower, 1, threadSet, parentFlower, deadEndComponent, bigFlowers, pinchEndsToEnds);
    stHash_destruct(pinchEndsToEnds);
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

void stCaf_finish(Flower *flower, stPinchThreadSet *threadSet, int64_t chainLengthForBigFlower,
        int64_t longChain, int64_t minLengthForChromosome,
        double proportionOfUnalignedBasesForNewChromosome) {
    stCactusNode *startCactusNode;
    stList *deadEndComponent;
    stCactusGraph *cactusGraph = stCaf_getCactusGraphForThreadSet(flower, threadSet, &startCactusNode, &deadEndComponent, 1, minLengthForChromosome,
            proportionOfUnalignedBasesForNewChromosome, 0, INT64_MAX);
    //chainLengthForBigFlower = 10;
    //longChain = 10;
    stSet *bigFlowers = stCactusGraph_collapseLongChainsOfBigFlowers(cactusGraph, startCactusNode, chainLengthForBigFlower, longChain, stCaf_mergeNodeObjects, 0);
    if(stSet_size(bigFlowers) > 0) {
        printf("We are collapsing the chains of %" PRIi64 " flowers\n", stSet_size(bigFlowers));
    }
    //Convert cactus graph/pinch graph to API
    stCaf_convertCactusGraphToFlowers(threadSet, startCactusNode, flower, deadEndComponent, bigFlowers);
    //Cleanup
    stCactusGraph_destruct(cactusGraph);
    stSet_destruct(bigFlowers);
}
