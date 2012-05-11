#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

static void stCaf_constructEmptyPinchGraphP(End *end, stPinchSegment *segment, bool orientation, stHash *endsToBlocks) {
    assert(stPinchSegment_getLength(segment) == 1);
    end = end_getPositiveOrientation(end);
    stPinchBlock *block = stHash_search(endsToBlocks, end);
    if (block == NULL) {
        stHash_insert(endsToBlocks, end, stPinchBlock_construct3(segment, orientation));
    } else {
        stPinchBlock_pinch2(block, segment, orientation);
    }
}

stPinchThreadSet *stCaf_constructEmptyPinchGraph(Flower *flower) {
    stPinchThreadSet *threadSet = stPinchThreadSet_construct();
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
                stPinchThread *thread = stPinchThreadSet_addThread(threadSet, cap_getName(cap), cap_getCoordinate(cap),
                        cap_getCoordinate(adjacentCap) - cap_getCoordinate(cap) + 1);
                stPinchThread_split(thread, cap_getCoordinate(cap));
                stPinchThread_split(thread, cap_getCoordinate(adjacentCap) - 1);
                stPinchSegment *_5PrimeSegment = stPinchThread_getFirst(thread);
                stPinchSegment *_3PrimeSegment = stPinchThread_getLast(thread);
                assert(stPinchSegment_getStart(_5PrimeSegment) == cap_getCoordinate(cap));
                assert(stPinchSegment_getStart(_3PrimeSegment) == cap_getCoordinate(adjacentCap));
                stCaf_constructEmptyPinchGraphP(end, _5PrimeSegment, 1, endsToBlocks);
                stCaf_constructEmptyPinchGraphP(cap_getEnd(adjacentCap), _3PrimeSegment, 0, endsToBlocks);
            }
        }
        end_destructInstanceIterator(capIt);
    }
    flower_destructEndIterator(endIt);
    stHash_destruct(endsToBlocks);
    return threadSet;
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

stPinchThreadSet *stCaf_setup(Flower *flower) {
    //Setup the empty flower that will be filled out
    initialiseFlowerForFillingOut(flower);

    //Create empty pinch graph from flower
    stPinchThreadSet *threadSet = stCaf_constructEmptyPinchGraph(flower);

    return threadSet;
}
