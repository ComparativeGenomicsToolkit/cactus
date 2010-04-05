#ifndef CACTUS_GLOBALS_H_
#define CACTUS_GLOBALS_H_

#include <inttypes.h>

/*
 * Includes for Tokyo Cabinet.
 */
#include <tcutil.h>
#include <tcbdb.h>
/*
 * For the sorted sets.
 */
#include "avl.h"
/*
 * For the hashtables.
 */
#include "hashTableC.h"
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

typedef int64_t Name;
#define NULL_NAME INT64_MAX

typedef struct _metaEvent MetaEvent;
typedef struct _event Event;
typedef struct _eventTree EventTree;
typedef struct _metaSequence MetaSequence;
typedef struct _sequence Sequence;
typedef struct _end End;
typedef struct _cap Cap;
typedef struct _segment Segment;
typedef struct _block Block;
typedef struct _group Group;
typedef struct _link Link;
typedef struct _chain Chain;
typedef struct _face Face;
typedef struct _net Net;
typedef struct _netDisk NetDisk;
typedef struct _pseudoChromosome PseudoChromosome;
typedef struct _pseudoAdjacency PseudoAdjacency;
typedef struct _reference Reference;
typedef struct _cactusHash Hash;

typedef struct avl_traverser EventTree_Iterator;
typedef struct _end_instanceIterator End_InstanceIterator;
typedef struct _block_instanceIterator Block_InstanceIterator;
typedef struct avl_traverser Group_EndIterator;
typedef struct avl_traverser Net_SequenceIterator;
typedef struct avl_traverser Net_CapIterator;
typedef struct avl_traverser Net_EndIterator;
typedef struct avl_traverser Net_BlockIterator;
typedef struct avl_traverser Net_GroupIterator;
typedef struct avl_traverser Net_ChainIterator;
typedef struct avl_traverser Net_FaceIterator;
typedef struct avl_traverser Net_ReferenceIterator;
typedef BDBCUR NetDisk_NetNameIterator;
typedef struct avl_traverser NetDisk_NetIterator;
typedef struct avl_traverser Reference_PseudoChromosomeIterator;
typedef struct avl_traverser PseudoChromsome_PseudoAdjacencyIterator;

#endif
