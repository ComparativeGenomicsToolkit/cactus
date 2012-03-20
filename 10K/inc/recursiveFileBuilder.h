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

void recursiveFileBuilder_destruct(RecursiveFileBuilder *recursiveFileBuilder);

void recursiveFileBuilder_writeThread(RecursiveFileBuilder *recursiveFileBuilder,
        Cap *cap, void (*segmentWriteFn)(FILE *, Segment *));

char *recursiveFileBuilder_getUniqueFileName(Flower *flower, const char *directory);

#endif /* RECURSIVEFILEBUILDER_H_ */
