#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

static void stCaf_addCapsToPinchGraphP(End *end, stPinchSegment *segment, bool orientation, stHash *endsToBlocks) {
    assert(stPinchSegment_getLength(segment) == 1);
    end = end_getPositiveOrientation(end);
    stPinchBlock *block = stHash_search(endsToBlocks, end);
    if (block == NULL) {
        stHash_insert(endsToBlocks, end, stPinchBlock_construct3(segment, orientation));
    } else {
        stPinchBlock_pinch2(block, segment, orientation);
    }
}

void stCaf_addCapsToPinchGraph(Flower *flower, stPinchThreadSet *threadSet) {
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    stHash *endsToBlocks = stHash_construct();
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        assert(!end_isBlockEnd(end));
        End_InstanceIterator *capIt = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(capIt)) != NULL) {
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if (!cap_getSide(cap)) {
                Cap *adjacentCap = cap_getAdjacency(cap);
                assert(cap_getSide(adjacentCap));
                assert(cap_getCoordinate(cap) < cap_getCoordinate(adjacentCap));
                stPinchThread *thread = stPinchThreadSet_getThread(threadSet, cap_getName(cap));
                stPinchThread_split(thread, cap_getCoordinate(cap));
                stPinchThread_split(thread, cap_getCoordinate(adjacentCap) - 1);
                stPinchSegment *_5PrimeSegment = stPinchThread_getFirst(thread);
                stPinchSegment *_3PrimeSegment = stPinchThread_getLast(thread);
                assert(stPinchSegment_getStart(_5PrimeSegment) == cap_getCoordinate(cap));
                assert(stPinchSegment_getStart(_3PrimeSegment) == cap_getCoordinate(adjacentCap));
                stCaf_addCapsToPinchGraphP(end, _5PrimeSegment, 1, endsToBlocks);
                stCaf_addCapsToPinchGraphP(cap_getEnd(adjacentCap), _3PrimeSegment, 0, endsToBlocks);
            }
        }
        end_destructInstanceIterator(capIt);
    }
    flower_destructEndIterator(endIt);
    stHash_destruct(endsToBlocks);
}

void stCaf_addThreadsToPinchGraph(Flower *flower, stPinchThreadSet *threadSet) {
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        assert(!end_isBlockEnd(end));
        End_InstanceIterator *capIt = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(capIt)) != NULL) {
            cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
            if (!cap_getSide(cap)) {
                Cap *adjacentCap = cap_getAdjacency(cap);
                assert(cap_getSide(adjacentCap));
                assert(cap_getCoordinate(cap) < cap_getCoordinate(adjacentCap));
                stPinchThreadSet_addThread(threadSet, cap_getName(cap), cap_getCoordinate(cap),
                                           cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) + 1);
            }
        }
        end_destructInstanceIterator(capIt);
    }
    flower_destructEndIterator(endIt);
}

static void initialiseFlowerForFillingOut(Flower *flower) {
    assert(!flower_builtBlocks(flower)); //We can't do this if we've already built blocks for the flower!.
    assert(flower_isTerminal(flower));
    assert(flower_getGroupNumber(flower) == 1);
    assert(group_isLeaf(flower_getFirstGroup(flower))); //this should be true by the previous assert
    //Destruct any chain
    assert(flower_getChainNumber(flower) <= 1);
    if (flower_getChainNumber(flower) == 1) {
        Chain *chain = flower_getFirstChain(flower);
        chain_destruct(chain);
    }
    group_destruct(flower_getFirstGroup(flower));
}


stPinchThreadSet *stCaf_setupForOnlineCactus(Flower *flower, stOnlineCactus **cactus) {
    //Setup the empty flower that will be filled out
    initialiseFlowerForFillingOut(flower);

    stPinchThreadSet *threadSet = stPinchThreadSet_construct();

    *cactus = stOnlineCactus_construct(
        (void *(*)(void *, bool)) stPinchBlock_getRepresentativeSegmentCap,
        (void *(*)(void *)) stPinchSegmentCap_getBlock);
    stPinchThreadSet_setAdjComponentCreationCallback(threadSet, (void (*)(void *, stConnectedComponent *)) stOnlineCactus_createNode, *cactus);
    stPinchThreadSet_setBlockCreationCallback(threadSet, (void (*)(void *, stConnectedComponent *, stConnectedComponent *, stPinchSegmentCap *, stPinchSegmentCap *, stPinchBlock *)) stOnlineCactus_addEdge, *cactus);
    stPinchThreadSet_setBlockDeletionCallback(threadSet, (void (*)(void *, stPinchSegmentCap *, stPinchSegmentCap *, stPinchBlock *)) stOnlineCactus_deleteEdge, *cactus);
    stPinchThreadSet_setAdjComponentMergeCallback(threadSet, (void (*)(void *, stConnectedComponent *, stConnectedComponent *)) stOnlineCactus_nodeMerge, *cactus);
    stPinchThreadSet_setAdjComponentCleaveCallback(threadSet, (void (*)(void *, stConnectedComponent *, stConnectedComponent *, stSet *)) stOnlineCactus_nodeCleave, *cactus);
    stPinchThreadSet_setAdjComponentDeletionCallback(threadSet, (void (*)(void *, stConnectedComponent *)) stOnlineCactus_deleteNode, *cactus);
    stOnlineCactus_setWeightFn(*cactus, (uint64_t (*)(const void *)) stPinchBlock_getLength);

    //Create empty pinch graph from flower
    stCaf_addThreadsToPinchGraph(flower, threadSet);

    return threadSet;
}

void stCaf_disableAndCleanupOnlineCactus(stPinchThreadSet *threadSet, stOnlineCactus *cactus) {
    stPinchThreadSet_setAdjComponentCreationCallback(threadSet, NULL, NULL);
    stPinchThreadSet_setBlockCreationCallback(threadSet, NULL, NULL);
    stPinchThreadSet_setBlockDeletionCallback(threadSet, NULL, NULL);
    stPinchThreadSet_setAdjComponentMergeCallback(threadSet, NULL, NULL);
    stPinchThreadSet_setAdjComponentCleaveCallback(threadSet, NULL, NULL);
    stPinchThreadSet_setAdjComponentDeletionCallback(threadSet, NULL, NULL);

    stOnlineCactus_destruct(cactus);
}

stPinchThreadSet *stCaf_setup(Flower *flower) {
    //Setup the empty flower that will be filled out
    initialiseFlowerForFillingOut(flower);

    //Create empty pinch graph from flower
    stPinchThreadSet *threadSet = stPinchThreadSet_construct();
    stCaf_addThreadsToPinchGraph(flower, threadSet);
    stCaf_addCapsToPinchGraph(flower, threadSet);

    return threadSet;
}
