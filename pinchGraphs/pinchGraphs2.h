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
} CaBlock;

typedef struct _caEnd {
    CaEnd *nEnd;
    CaBlock *block;
    bool blockOrientation;
    CaGroup *group;
    CaEnd *linkEnd;
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

void caSegment_destruct(CaSegment *segment);

CaSegment *caSegment_get3Prime(CaSegment *segment);

CaSegment *caSegment_get5Prime(CaSegment *segment);

bool caSegment_isStub(CaSegment *segment);

int64_t *caSegment_getStart(CaSegment *segment);

int64_t *caSegment_getEnd(CaSegment *segment);

int64_t *caSegment_getLength(CaSegment *segment);

CaBlock *caSegment_getBlock(CaSegment *segment);

bool caSegment_getBlockOrientation(CaSegment *segment);

CaSegment *caSegment_split(CaSegment *segment, int32_t splitPoint);

CaSegment *caSegment_joinIfTrivial(CaSegment *segment);

//Blocks

CaBlock *caBlock_construct(CaSegment *segment1, segment2);

void caBlock_destruct(CaBlock *block);

CaBlock *caBlock_pinch(CaBlock *block1, CaBlock *block2);

void caBlock_cleave(CaBlock *block);

CaBlockIt caBlock_getSegmentIterator(CaBlock *block);

CaSegment *caBlockIt_getNext(CaBlockIt caBlockIt);

CaEnd *caBlock_get5End(CaBlock *block);

CaEnd *caBlock_get3End(CaBlock *block);

//End

CaEnd *caEnd_construct(CaBlock *block, bool orientation);

CaGroup *caEnd_getGroup(CaEnd *end);

CaNet *caEnd_getEnd(CaEnd *end);






#endif /* PINCHGRAPHS2_H_ */
