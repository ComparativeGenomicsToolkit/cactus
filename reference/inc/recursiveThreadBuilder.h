/*
 * recursiveThreadBuilder.h
 *
 *  Created on: 21 Jun 2012
 *      Author: benedictpaten
 */

#ifndef RECURSIVETHREADBUILDER_H_
#define RECURSIVETHREADBUILDER_H_

void buildRecursiveThreads(stKVDatabase *database, stList *caps,
                           char *(*segmentWriteFn)(Segment *, void *),
                           char *(*terminalAdjacencyWriteFn)(Cap *, void *), void *extraArg);

stList *buildRecursiveThreadsInList(stKVDatabase *database, stList *caps,
        char *(*segmentWriteFn)(Segment *, void *),
        char *(*terminalAdjacencyWriteFn)(Cap *, void *), void *extraArg);

typedef stHash RecordHolder;

RecordHolder *recordHolder_construct();

void recordHolder_destruct(RecordHolder *rh);

/*
 * Removes the records from rhToAdd and puts them in rhToAddTo, leaving rhToAdd empty.
 */
void recordHolder_transferAll(RecordHolder *rhToAddTo, RecordHolder *rhToAdd);

void buildRecursiveThreadsNoDb(RecordHolder *rh, stList *caps, char *(*segmentWriteFn)(Segment *, void *),
                               char *(*terminalAdjacencyWriteFn)(Cap *, void *), void *extraArg);

stList *buildRecursiveThreadsInListNoDb(RecordHolder *rh, stList *caps, char *(*segmentWriteFn)(Segment *, void *),
                                        char *(*terminalAdjacencyWriteFn)(Cap *, void *), void *extraArg);

#endif /* RECURSIVETHREADBUILDER_H_ */
