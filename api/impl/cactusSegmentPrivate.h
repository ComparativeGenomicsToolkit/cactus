/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_ATOM_INSTANCE_PRIVATE_H_
#define CACTUS_ATOM_INSTANCE_PRIVATE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private segment functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs segment with the two caps, which must both have the same instance name.
 */
Segment *segment_construct3(Name name, Block *block,
		Cap *_5Cap, Cap *_3Cap);

/*
 * Destruct the segment, does not destruct ends.
 */
void segment_destruct(Segment *segment);

#endif
