/*
 * recursiveFileBuilder.h
 *
 *  Created on: 16 Mar 2012
 *      Author: benedictpaten
 */

#ifndef RECURSIVEFILEBUILDER_H_
#define RECURSIVEFILEBUILDER_H_

typedef struct _recursiveFileBuilder RecursiveFileBuilder;

RecursiveFileBuilder *recursiveFileBuilder_construct(const char *childDir,
        FILE *parentFileHandle, bool hasParent);

void recursiveFileBuilder_writeHead(RecursiveFileBuilder *recursiveFileBuilder, Cap *cap);

void recursiveFileBuilder_writeTail(RecursiveFileBuilder *recursiveFileBuilder,
        Cap *cap);

void recursiveFileBuilder_writeAdjacency(
        RecursiveFileBuilder *recursiveFileBuilder, Cap *cap);

void recursiveFileBuilder_writeSegment(
        RecursiveFileBuilder *recursiveFileBuilder, Segment *segment,
        void (*segmentWriteFn)(FILE *, Segment *));

#endif /* RECURSIVEFILEBUILDER_H_ */
