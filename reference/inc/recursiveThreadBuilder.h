/*
 * recursiveThreadBuilder.h
 *
 *  Created on: 21 Jun 2012
 *      Author: benedictpaten
 */

#ifndef RECURSIVETHREADBUILDER_H_
#define RECURSIVETHREADBUILDER_H_

void buildRecursiveThreads(stKVDatabase *database, stList *caps,
        char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *));

stList *buildRecursiveThreadsInList(stKVDatabase *database, stList *caps,
        char *(*segmentWriteFn)(Segment *),
        char *(*terminalAdjacencyWriteFn)(Cap *));

typedef stHash RecordHolder;

RecordHolder *recordHolder_construct();

void recordHolder_destruct(RecordHolder *rh);

void buildRecursiveThreadsNoDb(RecordHolder *rh, stList *caps, char *(*segmentWriteFn)(Segment *),
                               char *(*terminalAdjacencyWriteFn)(Cap *));

stList *buildRecursiveThreadsInListNoDb(RecordHolder *rh, stList *caps, char *(*segmentWriteFn)(Segment *),
                                        char *(*terminalAdjacencyWriteFn)(Cap *));

#endif /* RECURSIVETHREADBUILDER_H_ */
