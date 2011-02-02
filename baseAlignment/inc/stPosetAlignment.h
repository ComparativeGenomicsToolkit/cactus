/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * stPosetAlignment.h
 *
 *  Created on: 21-Jun-2010
 *      Author: benedictpaten
 */

#include "sonLib.h"

#ifndef STPOSETALIGNMENT_H_
#define STPOSETALIGNMENT_H_

/*
 * The anonmyous typedef declaration of the alignment structure.
 */
typedef struct _stPosetAlignment stPosetAlignment;

/*
 * Constructs a poset alignment containing the given number of sequences
 */
stPosetAlignment *stPosetAlignment_construct(int32_t sequenceNumber);

/*
 * Destructs the poset alignment.
 */
void stPosetAlignment_destruct(stPosetAlignment *posetAlignment);

/*
 * Returns the number of sequences in the alignment.
 */
int32_t stPosetAlignment_getSequenceNumber(stPosetAlignment *posetAlignment);

/*
 * Returns non-zero if the pair can be added to the alignment, but does not actually add the pair to the alignment.
 */
bool stPosetAlignment_isPossible(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2);

/*
 * List stPosetAlignment_isPossible, but if the return value is non-zero the pair is added to the alignment.
 */
bool stPosetAlignment_add(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2);

#endif /* STPOSETALIGNMENT_H_ */
