/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * getSequences.c
 *
 *  Created on: 24 Jun 2010
 *      Author: benedictpaten
 */

#include "adjacencySequences.h"

/*
 * Gets the raw sequence.
 */
static char *getAdjacencySequenceP(Cap *cap, int64_t maxLength) {
    Sequence *sequence = cap_getSequence(cap);
    assert(sequence != NULL);
    Cap *cap2 = cap_getAdjacency(cap);
    assert(cap2 != NULL);
    assert(!cap_getSide(cap));

    if (cap_getStrand(cap)) {
        int64_t length = cap_getCoordinate(cap2) - cap_getCoordinate(cap) - 1;
        assert(length >= 0);
        assert(maxLength >= 0);
        return sequence_getString(sequence, cap_getCoordinate(cap) + 1, length
                > maxLength ? maxLength : length, 1);
    } else {
        int64_t length = cap_getCoordinate(cap) - cap_getCoordinate(cap2) - 1;
        assert(length >= 0);
        return sequence_getString(sequence,
                length > maxLength ? cap_getCoordinate(cap) - maxLength
                        : cap_getCoordinate(cap2) + 1,
                length > maxLength ? maxLength : length, 0);
    }
}

AdjacencySequence *adjacencySequence_construct(Cap *cap, int64_t maxLength) {
    AdjacencySequence *subSequence = (AdjacencySequence *) st_malloc(
            sizeof(AdjacencySequence));
    subSequence->string = getAdjacencySequenceP(cap, maxLength);
    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    assert(!cap_getSide(cap));
    assert(cap_getSequence(cap) != NULL);
    subSequence->subsequenceIdentifier = cap_getName(cap_getStrand(cap) ? cap : adjacentCap);
    subSequence->strand = cap_getStrand(cap);
    subSequence->start = cap_getCoordinate(cap) + (cap_getStrand(cap) ? 1 : -1);
    subSequence->length = strlen(subSequence->string);
    subSequence->hasStubEnd = end_isFree(cap_getEnd(adjacentCap)) && end_isStubEnd(cap_getEnd(adjacentCap));
    return subSequence;
}

void adjacencySequence_destruct(AdjacencySequence *subSequence) {
    free(subSequence->string);
    free(subSequence);
}

