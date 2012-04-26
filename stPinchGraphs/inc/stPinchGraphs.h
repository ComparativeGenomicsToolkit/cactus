/*
 * pinchGraphs2.h
 *
 *  Created on: 9 Mar 2012
 *      Author: benedictpaten
 */

#ifndef ST_PINCH_GRAPHS_H_
#define ST_PINCH_GRAPHS_H_

#include "sonLib.h"

//The exception string
extern const char *ST_PINCH_GRAPH_EXCEPTION_ID;

//Datastructures

typedef struct _stThreadSet stThreadSet;

typedef struct _stThreadIt {
    stThreadSet *threadSet;
    int32_t index;
} stThreadIt;

typedef struct _stThread stThread;

typedef struct _stSegment stSegment;

typedef struct _stThreadSetSegmentIt {
    stThreadIt threadIt;
    stSegment *segment;
} stThreadSetSegmentIt;

typedef struct _stThreadSetBlockIt {
    stThreadSetSegmentIt segmentIt;
} stThreadSetBlockIt;

typedef struct _stBlock stBlock;

typedef struct _stBlockIt {
    stSegment *segment;
} stBlockIt;

typedef struct _stEnd {
    stBlock *block;
    bool orientation;
} stEnd;

//Thread set

stThreadSet *stThreadSet_construct();

void stThreadSet_destruct(stThreadSet *threadSet);

stThread *stThreadSet_addThread(stThreadSet *threadSet, int64_t name, int64_t start, int64_t length);

stThread *stThreadSet_getThread(stThreadSet *threadSet, int64_t name);

stSegment *stThreadSet_getSegment(stThreadSet *threadSet, int64_t name, int64_t coordinate);

int32_t stThreadSet_getSize(stThreadSet *threadSet);

stThreadIt stThreadSet_getIterator(stThreadSet *threadSet);

stThread *stThreadIt_getNext(stThreadIt *);

void stThreadSet_joinTrivialBoundaries(stThreadSet *threadSet);

stSegment *stThreadSet_getSegment(stThreadSet *threadSet, int64_t name, int64_t coordinate);

stList *stThreadSet_getAdjacencyComponents(stThreadSet *threadSet);

//convenience functions

stThreadSetSegmentIt stThreadSet_getSegmentIt(stThreadSet *threadSet);

stSegment *stThreadSetSegmentIt_getNext(stThreadSetSegmentIt *segmentIt);

stThreadSetBlockIt stThreadSet_getBlockIt(stThreadSet *threadSet);

stBlock *stThreadSetBlockIt_getNext(stThreadSetBlockIt *blockIt);

//Thread

int64_t stThread_getName(stThread *stThread);

int64_t stThread_getStart(stThread *stThread);

int64_t stThread_getLength(stThread *stThread);

stSegment *stThread_getSegment(stThread *stThread, int64_t coordinate);

stSegment *stThread_getFirst(stThread *stThread);

void stThread_split(stThread *thread, int64_t leftSideOfSplitPoint);

void stThread_joinTrivialBoundaries(stThread *thread);

void stThread_pinch(stThread *thread1, stThread *thread2, int64_t start1, int64_t start2, int64_t length, bool strand2);

//Segments

int64_t stSegment_getStart(stSegment *segment);

int64_t stSegment_getLength(stSegment *segment);

int64_t stSegment_getName(stSegment *segment);

stSegment *stSegment_get3Prime(stSegment *segment);

stSegment *stSegment_get5Prime(stSegment *segment);

stThread *stSegment_getThread(stSegment *segment);

stBlock *stSegment_getBlock(stSegment *segment);

bool stSegment_getBlockOrientation(stSegment *segment);

//Blocks

stBlock *stBlock_construct2(stSegment *segment1); //Allows the specification of a block with just one element

stBlock *stBlock_construct(stSegment *segment1, bool orientation1, stSegment *segment2, bool orientation2);

stBlock *stBlock_pinch(stBlock *block1, stBlock *block2, bool orientation);

stBlock *stBlock_pinch2(stBlock *block1, stSegment *segment, bool orientation);

stBlockIt stBlock_getSegmentIterator(stBlock *block);

stSegment *stBlockIt_getNext(stBlockIt *stBlockIt);

void stBlock_destruct(stBlock *block);

int64_t stBlock_getLength(stBlock *block);

stSegment *stBlock_getFirst(stBlock *block);

uint32_t stBlock_getDegree(stBlock *block);

//Block ends

stEnd *stEnd_construct(stBlock *block, bool orientation);

void stEnd_destruct(stEnd *end);

stBlock *stEnd_getBlock(stEnd *end);

bool stEnd_getOrientation(stEnd *end);

int stEnd_equalsFn(const void *, const void *);

uint32_t stEnd_hashFn(const void *);

bool stEnd_traverse5Prime(bool endOrientation, stSegment *segment);

bool stEnd_endOrientation(bool _5PrimeTraversal, stSegment *segment);

bool stEnd_boundaryIsTrivial(stEnd end);

void stEnd_joinTrivialBoundary(stEnd end);

#endif /* ST_PINCH_GRAPHS_H_ */
