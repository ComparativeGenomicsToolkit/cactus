/*
 * pinchGraphs2.h
 *
 *  Created on: 9 Mar 2012
 *      Author: benedictpaten
 */

#ifndef PINCHGRAPHS2_H_
#define PINCHGRAPHS2_H_

#include "cactus.h"

//Basic data structures

//Destruction works upwards, in the following order segment ( block ( end ( group ( net )), chain ) )

typedef struct _caSegment {
    int64_t start;
    CaSegment *pSegment;
    CaSegment *nSegment;
    CaBlock *block;
    bool blockOrientation;
    CaSegment *nBlockSegment;
} CaSegment;

typedef struct _caBlock {
    CaSegment *headSegment;
    CaSegment *tailSegment;
    CaEnd *_PEnd;
    CaEnd *_3End;
    CaChain *chain;
    CaBlock *nChainBlock;
} CaBlock;

typedef struct _caEnd {
    CaEnd *nEnd;
    CaBlock *block;
    bool blockOrientation;
    CaGroup *group;
} CaEnd;

typedef struct _caGroup {
    CaEnd *headEnd;
    CaEnd *tailEnd;
    CaNet *net;
    CaGroup *nNetGroup;
} CaGroup;

typedef struct _caNet {
    CaGroup *headGroup;
    CaGroup *tailGroup;
    CaChain *headChain;
    CaChain *tailChain;
} CaNet;

typedef struct _caChain {
    CaBlock *headBlock;
    CaBlock *tailBlock;
    CaNet *parentNet;
    CaChain *nNetChain;
} CaChain;

//Segments

CaSegment *caSegment_construct(int64_t start, int64_t end, bool attached);

CaSegment *caSegment_get3Prime(CaSegment *segment);

CaSegment *caSegment_get5Prime(CaSegment *segment);

int64_t *caSegment_getStart(CaSegment *segment);

int64_t *caSegment_getEnd(CaSegment *segment);

int64_t *caSegment_getLength(CaSegment *segment);

CaBlock *caSegment_getBlock(CaSegment *segment);

bool caSegment_getOrientation(CaSegment *segment);

PinchSegment *caSegment_split(CaSegment *segment, int32_t splitPoint);

PinchSegment *caSegment_joinIfTrivial(CaSegment *segment);

//Blocks

PinchBlock *pinchBlock_pinch(PinchBlock *pinchBlock1, PinchBlock *pinchBlock2);

void pinchBlock_cleave(PinchBlock *pinchBlock1);

PinchBlockIt pinchBlock_getSegmentIterator(PinchBlock *block);

PinchSegment *pinchBlockIt_getNext(PinchBlockIt pinchBlockIt);


//Cactus




#endif /* PINCHGRAPHS2_H_ */
