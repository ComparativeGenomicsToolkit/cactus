/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_SEQUENCE_H_
#define CACTUS_SEQUENCE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Meta sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a meta sequence, which contains all the essential info for a sequence.
 *
 * This function is NOT thread safe, do not try to have concurrent instances of this function!
 */
Sequence *sequence_construct(int64_t start, int64_t length, const char *string, const char *header,
		Event *event, CactusDisk *cactusDisk);

/*
 * Adds the isTrivialSequence field.
 */
Sequence *sequence_construct3(int64_t start, int64_t length, const char *string, const char *header, Event *event,
        bool isTrivialSequence, CactusDisk *cactusDisk);

/*
 * Gets the name of the sequence.
 */
Name sequence_getName(Sequence *sequence);

/*
 * Gets the start coordinate of the sequence.
 */
int64_t sequence_getStart(Sequence *sequence);

/*
 * Gets the length of the sequence.
 */
int64_t sequence_getLength(Sequence *sequence);

/*
 * Gets the associated event name.
 */
Event *sequence_getEvent(Sequence *sequence);

/*
 * Gets a string for representing a subsequence of the meta sequence.
 */
char *sequence_getString(Sequence *sequence, int64_t start, int64_t length, int64_t strand);

/*
 * Gets the header line associated with the meta sequence.
 */
const char *sequence_getHeader(Sequence *sequence);

/*
 * Returns flag indicating if sequence is trivial.
 */
bool sequence_isTrivialSequence(Sequence *sequence);

/*
 * Sets the header line associated with the meta sequence.
 */
void sequence_setHeader(Sequence *sequence, char *newHeader);

#endif
