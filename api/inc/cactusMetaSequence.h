/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_META_SEQUENCE_H_
#define CACTUS_META_SEQUENCE_H_

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
MetaSequence *metaSequence_construct(int32_t start, int32_t length, const char *string, const char *header,
		Name eventName, CactusDisk *cactusDisk);

/*
 * Gets the name of the sequence.
 */
Name metaSequence_getName(MetaSequence *metaSequence);

/*
 * Gets the start coordinate of the sequence.
 */
int32_t metaSequence_getStart(MetaSequence *metaSequence);

/*
 * Gets the length of the sequence.
 */
int32_t metaSequence_getLength(MetaSequence *metaSequence);

/*
 * Gets the associated event name.
 */
Name metaSequence_getEventName(MetaSequence *metaSequence);

/*
 * Gets a string for representing a subsequence of the meta sequence.
 */
char *metaSequence_getString(MetaSequence *metaSequence, int32_t start, int32_t length, int32_t strand);

/*
 * Gets the header line associated with the meta sequence.
 */
const char *metaSequence_getHeader(MetaSequence *metaSequence);

#endif
