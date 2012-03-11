/*
 * pinchGraphs2.h
 *
 *  Created on: 9 Mar 2012
 *      Author: benedictpaten
 */

#ifndef PINCHGRAPHS2_H_
#define PINCHGRAPHS2_H_

#include "cactus.h"

typedef struct _stPGraph {

} StPGraph;

typedef struct _pinchSegment {
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
} PinchBlock;

typedef struct _end {
        PinchBlock *pinchBlock;
        bool blockOrientation;
        Group *group;
} GroupEdge;

typedef struct _group {
        stList *groupEnds;
} Group;

typedef struct _net {

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
