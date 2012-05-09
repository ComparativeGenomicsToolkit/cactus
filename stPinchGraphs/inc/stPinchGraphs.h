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

typedef struct _stPinchThreadSet stPinchThreadSet;

typedef struct _stPinchThreadIt {
    stPinchThreadSet *threadSet;
    int32_t index;
} stPinchThreadSetIt;

typedef struct _stPinchThread stPinchThread;

typedef struct _stPinchSegment stPinchSegment;

typedef struct _stPinchThreadSetSegmentIt {
    stPinchThreadSetIt threadIt;
    stPinchSegment *segment;
} stPinchThreadSetSegmentIt;

typedef struct _stPinchThreadSetBlockIt {
    stPinchThreadSetSegmentIt segmentIt;
} stPinchThreadSetBlockIt;

typedef struct _stPinchBlock stPinchBlock;

typedef struct _stPinchBlockIt {
    stPinchSegment *segment;
} stPinchBlockIt;

typedef struct _stPinchEnd {
    stPinchBlock *block;
    bool orientation;
} stPinchEnd;

typedef struct _stPinch {
    int64_t name1;
    int64_t name2;
    int64_t start1;
    int64_t start2;
    int64_t length;
    bool strand;
} stPinch;

typedef struct _stPinchInterval {
    int64_t name;
    int64_t start;
    int64_t length;
    void *label;
} stPinchInterval;

//Thread set

stPinchThreadSet *stPinchThreadSet_construct();

void stPinchThreadSet_destruct(stPinchThreadSet *threadSet);

stPinchThread *stPinchThreadSet_addThread(stPinchThreadSet *threadSet, int64_t name, int64_t start, int64_t length);

stPinchThread *stPinchThreadSet_getThread(stPinchThreadSet *threadSet, int64_t name);

stPinchSegment *stPinchThreadSet_getSegment(stPinchThreadSet *threadSet, int64_t name, int64_t coordinate);

int32_t stPinchThreadSet_getSize(stPinchThreadSet *threadSet);

stPinchThreadSetIt stPinchThreadSet_getIt(stPinchThreadSet *threadSet);

stPinchThread *stPinchThreadSetIt_getNext(stPinchThreadSetIt *);

void stPinchThreadSet_joinTrivialBoundaries(stPinchThreadSet *threadSet);

stPinchSegment *stPinchThreadSet_getSegment(stPinchThreadSet *threadSet, int64_t name, int64_t coordinate);

int32_t stPinchThreadSet_getTotalBlockNumber(stPinchThreadSet *threadSet);

stList *stPinchThreadSet_getAdjacencyComponents(stPinchThreadSet *threadSet);

stList *stPinchThreadSet_getAdjacencyComponents2(stPinchThreadSet *threadSet, stHash **edgeEndsToAdjacencyComponents);

stSortedSet *stPinchThreadSet_getThreadComponents(stPinchThreadSet *threadSet);

stPinchThreadSet *stPinchThreadSet_getRandomEmptyGraph();

stPinch stPinchThreadSet_getRandomPinch(stPinchThreadSet *threadSet);

stPinchThreadSet *stPinchThreadSet_getRandomGraph();

//convenience functions

stPinchThreadSetSegmentIt stPinchThreadSet_getSegmentIt(stPinchThreadSet *threadSet);

stPinchSegment *stPinchThreadSetSegmentIt_getNext(stPinchThreadSetSegmentIt *segmentIt);

stPinchThreadSetBlockIt stPinchThreadSet_getBlockIt(stPinchThreadSet *threadSet);

stPinchBlock *stPinchThreadSetBlockIt_getNext(stPinchThreadSetBlockIt *blockIt);

//Thread

int64_t stPinchThread_getName(stPinchThread *stPinchThread);

int64_t stPinchThread_getStart(stPinchThread *stPinchThread);

int64_t stPinchThread_getLength(stPinchThread *stPinchThread);

stPinchSegment *stPinchThread_getSegment(stPinchThread *stPinchThread, int64_t coordinate);

stPinchSegment *stPinchThread_getFirst(stPinchThread *stPinchThread);

stPinchSegment *stPinchThread_getLast(stPinchThread *thread);

void stPinchThread_split(stPinchThread *thread, int64_t leftSideOfSplitPoint);

void stPinchThread_joinTrivialBoundaries(stPinchThread *thread);

void stPinchThread_pinch(stPinchThread *thread1, stPinchThread *thread2, int64_t start1, int64_t start2, int64_t length, bool strand2);

//Segments

int64_t stPinchSegment_getStart(stPinchSegment *segment);

int64_t stPinchSegment_getLength(stPinchSegment *segment);

int64_t stPinchSegment_getName(stPinchSegment *segment);

stPinchSegment *stPinchSegment_get3Prime(stPinchSegment *segment);

stPinchSegment *stPinchSegment_get5Prime(stPinchSegment *segment);

stPinchThread *stPinchSegment_getThread(stPinchSegment *segment);

stPinchBlock *stPinchSegment_getBlock(stPinchSegment *segment);

bool stPinchSegment_getBlockOrientation(stPinchSegment *segment);

void stPinchSegment_split(stPinchSegment *segment, int64_t leftSideOfSplitPoint);

//Blocks

stPinchBlock *stPinchBlock_construct3(stPinchSegment *segment, bool orientation);

stPinchBlock *stPinchBlock_construct2(stPinchSegment *segment1); //Allows the specification of a block with just one element

stPinchBlock *stPinchBlock_construct(stPinchSegment *segment1, bool orientation1, stPinchSegment *segment2, bool orientation2);

stPinchBlock *stPinchBlock_pinch(stPinchBlock *block1, stPinchBlock *block2, bool orientation);

stPinchBlock *stPinchBlock_pinch2(stPinchBlock *block1, stPinchSegment *segment, bool orientation);

stPinchBlockIt stPinchBlock_getSegmentIterator(stPinchBlock *block);

stPinchSegment *stPinchBlockIt_getNext(stPinchBlockIt *stPinchBlockIt);

void stPinchBlock_destruct(stPinchBlock *block);

int64_t stPinchBlock_getLength(stPinchBlock *block);

stPinchSegment *stPinchBlock_getFirst(stPinchBlock *block);

uint32_t stPinchBlock_getDegree(stPinchBlock *block);

void stPinchBlock_trim(stPinchBlock *block, int32_t blockEndTrim);

//Block ends

void stPinchEnd_fillOut(stPinchEnd *end, stPinchBlock *block, bool orientation);

stPinchEnd *stPinchEnd_construct(stPinchBlock *block, bool orientation);

stPinchEnd stPinchEnd_constructStatic(stPinchBlock *block, bool orientation);

void stPinchEnd_destruct(stPinchEnd *end);

stPinchBlock *stPinchEnd_getBlock(stPinchEnd *end);

bool stPinchEnd_getOrientation(stPinchEnd *end);

int stPinchEnd_equalsFn(const void *, const void *);

uint32_t stPinchEnd_hashFn(const void *);

bool stPinchEnd_traverse5Prime(bool endOrientation, stPinchSegment *segment);

bool stPinchEnd_endOrientation(bool _5PrimeTraversal, stPinchSegment *segment);

bool stPinchEnd_boundaryIsTrivial(stPinchEnd end);

void stPinchEnd_joinTrivialBoundary(stPinchEnd end);

//Pinch structure

void stPinch_fillOut(stPinch *pinch, int64_t name1, int64_t name2, int64_t start1, int64_t start2, int64_t length, bool strand);

stPinch stPinch_constructStatic(int64_t name1, int64_t name2, int64_t start1, int64_t start2, int64_t length, bool strand);

stPinch *stPinch_construct(int64_t name1, int64_t name2, int64_t start1, int64_t start2, int64_t length, bool strand);

void stPinch_destruct(stPinch *pinch);

//Pinch interval structure

void stPinchInterval_fillOut(stPinchInterval *pinchInterval, int64_t name, int64_t start, int64_t length, void *label);

stPinchInterval stPinchInterval_constructStatic(int64_t name, int64_t start, int64_t length, void *label);

stPinchInterval *stPinchInterval_construct(int64_t name, int64_t start, int64_t length, void *label);

int64_t stPinchInterval_getName(stPinchInterval *pinchInterval);

int64_t stPinchInterval_getStart(stPinchInterval *pinchInterval);

int64_t stPinchInterval_getLength(stPinchInterval *pinchInterval);

void *stPinchInterval_getLabel(stPinchInterval *pinchInterval);

int stPinchInterval_compareFunction(const stPinchInterval *interval1, const stPinchInterval *interval2);

void stPinchInterval_destruct(stPinchInterval *pinchInterval);

stSortedSet *stPinchThreadSet_getLabelIntervals(stPinchThreadSet *threadSet, stHash *pinchEndsToLabels);

stPinchInterval *stPinchIntervals_getInterval(stSortedSet *pinchIntervals, int64_t name, int64_t position);

#endif /* ST_PINCH_GRAPHS_H_ */
