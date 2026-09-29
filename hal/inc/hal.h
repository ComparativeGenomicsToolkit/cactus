/*
 * hal.h
 *
 *  Created on: 21 Jun 2012
 *      Author: benedictpaten
 */

#ifndef HAL_H_
#define HAL_H_

#include "sonLib.h"
#include "cactus.h"
#include "recursiveThreadBuilder.h"

void makeHalFormat(Flower *flower, stKVDatabase *database, Name referenceEventName,
                   FILE *fileHandle);

void makeHalFormatNoDb(Flower *flower, RecordHolder *rh, Name referenceEventName, FILE *fileHandle);

/*
 * As makeHalFormatNoDb at the top level (it consumes rh), writing the .c2h to fileHandle
 * (if not NULL) and/or a HAL format 3 fragment to the new directory fragmentDir (if not
 * NULL): the same threads, with each sequence's DNA, through hal's C API (hal3_c.h).
 * fragmentTree names the subproblem's ancestor and its ingroup children, e.g.
 * "(a,b)anc;"; other events (outgroups) are left out of the fragment.
 */
void makeHalOutputNoDb(Flower *flower, RecordHolder *rh, Name referenceEventName, FILE *fileHandle,
                       const char *fragmentDir, const char *fragmentTree);

void printFastaSequences(Flower *flower, FILE *fileHandle, Name referenceEventName);

#endif /* HAL_H_ */
