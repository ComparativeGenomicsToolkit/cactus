#ifndef CACTUS_SEQUENCE_H_
#define CACTUS_SEQUENCE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic sequence functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Creates a sequence for a net, wrapping a meta sequence.
 */
Sequence *sequence_construct(MetaSequence *metaSequence, Net *net);

/*
 * Destructs the sequence.
 */
void sequence_destruct(Sequence *sequence);

/*
 * Gets the associated meta sequence.
 */
MetaSequence *sequence_getMetaSequence(Sequence *sequence);

/*
 * Gets the start coordinate of the sequence (inclusive).
 */
int32_t sequence_getStart(Sequence *sequence);

/*
 * Gets the length of the sequence.
 */
int32_t sequence_getLength(Sequence *sequence);

/*
 * Gets the name of the sequence.
 */
Name sequence_getName(Sequence *sequence);

/*
 * Gets the name associated with the sequence.
 */
Event *sequence_getEvent(Sequence *sequence);

/*
 * Gets a sub string of the the sequence, indexes must be equal to or greater than the start coordinate,
 * and less than the start coordinate plus the sequences length.
 * If the strand is negative then the sequence returned will be the reverse complement sequence, traversing
 * in the opposite direction.
 *
 * The returned string must be freed.
 */
char *sequence_getString(Sequence *sequence, int32_t start, int32_t length, int32_t strand);

/*
 * Gets the header line associated with the sequence.
 */
const char *sequence_getHeader(Sequence *sequence);

/*
 * Gets the net the sequence is associated with.
 */
Net *sequence_getNet(Sequence *sequence);

/*
 * Checks (amongst other things) the following:
 * Checks the sequence if contained in the parent has the same event.
 */
void sequence_check(Sequence *sequence);

#endif
