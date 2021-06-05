/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_GLOBALS_H_
#define CACTUS_GLOBALS_H_

#include <inttypes.h>

/*
 * For the basic lib stuff
 */
#include "sonLib.h"
/*
 * For lists
 */
#include "commonC.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic data structure declarations (contents hidden)
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

#define NULL_NAME INT64_MAX
typedef int64_t Name;
typedef struct _event Event;
typedef struct _eventTree EventTree;
typedef struct _sequence Sequence;
typedef struct _end End;
typedef struct _cap Cap;
typedef struct _cap Segment;
typedef struct _end Block;
typedef struct _group Group;
typedef struct _group Link;
typedef struct _chain Chain;
typedef struct _flower Flower;
typedef struct _cactusDisk CactusDisk;
typedef stSortedSetIterator EventTree_Iterator;
typedef struct _end_instanceIterator End_InstanceIterator;
typedef struct _block_instanceIterator Block_InstanceIterator;
typedef struct _group_endIterator Group_EndIterator;
typedef stListIterator Flower_SequenceIterator;
typedef stListIterator Flower_CapIterator;
typedef stListIterator Flower_EndIterator;
typedef stListIterator Flower_GroupIterator;
typedef stListIterator Flower_ChainIterator;

#endif
