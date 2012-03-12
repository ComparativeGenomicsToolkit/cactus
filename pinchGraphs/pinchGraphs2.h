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

typedef struct _caPinchSegment {
    int64_t start;
    PinchSegment *pSegment;
    PinchSegment *nSegment;
    PinchBlock *block;
    bool blockOrientation;
    PinchSegment *nBlockSegment;
} PinchSegment;

typedef struct _pinchBlock {
    PinchSegment *headSegment;
    PinchSegment *tailSegment;
    End *_5End;
    End *_3End;
    Chain *chain;
    CactusBlock *nChainBlock;
} PinchBlock;

typedef struct _end {
    End *nEnd;
    PinchBlock *pinchBlock;
    Group *group;
} End;

typedef struct _group {
    End *headEnd;
    End *tailEnd;
    Net *net;
    Group *nNetGroup;
} Group;

typedef struct _chain {
    CactusBlock *headCactusBlock;
    CactusBlock *tailCactusBlock;
    Net *parentNet;
    Chain *nNetChain;
}

typedef struct _net {
    Group *headGroup;
    Group *tailGroup;
    Chain *headChain;
    Chain *tailChain;
} Net;


//Segments

PinchSegment *pinchSegment_construct(Sequence *sequence, int32_t start, int32_t end);

PinchSegment *pinchSegment_get3Prime(PinchSegment *segment);

PinchSegment *pinchSegment_get5Prime(PinchSegment *segment);

int64_t *pinchSegment_getStart(PinchSegment *segment);

int64_t *pinchSegment_getEnd(PinchSegment *segment);

int64_t *pinchSegment_getLength(PinchSegment *segment);

PinchBlock *pinchSegment_getBlock(PinchSegment *segment);

bool pinchSegment_getOrientation(PinchSegment *segment);

PinchSegment *pinchSegment_split(PinchSegment *segment, int32_t splitPoint);

PinchSegment *pinchSegment_joinIfTrivial(PinchSegment *segment);

//Blocks

PinchBlock *pinchBlock_pinch(PinchBlock *pinchBlock1, PinchBlock *pinchBlock2);

void pinchBlock_cleave(PinchBlock *pinchBlock1);

PinchBlockIt pinchBlock_getSegmentIterator(PinchBlock *block);

PinchSegment *pinchBlockIt_getNext(PinchBlockIt pinchBlockIt);


//Cactus




#endif /* PINCHGRAPHS2_H_ */
