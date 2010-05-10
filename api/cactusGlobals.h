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
#include "hashTableC_itr.h"
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
typedef struct _faceEnd FaceEnd;
typedef struct _net Net;
typedef struct _netDisk NetDisk;
typedef struct _pseudoChromosome PseudoChromosome;
typedef struct _pseudoAdjacency PseudoAdjacency;
typedef struct _reference Reference;
typedef struct _cactusHash Hash;

#define SORTED_SET_ITERATOR struct avl_traverser
typedef SORTED_SET_ITERATOR EventTree_Iterator;
typedef struct _end_instanceIterator End_InstanceIterator;
typedef struct _block_instanceIterator Block_InstanceIterator;
typedef SORTED_SET_ITERATOR Group_EndIterator;
typedef SORTED_SET_ITERATOR Net_SequenceIterator;
typedef SORTED_SET_ITERATOR Net_CapIterator;
typedef SORTED_SET_ITERATOR Net_SegmentIterator;
typedef SORTED_SET_ITERATOR Net_EndIterator;
typedef SORTED_SET_ITERATOR Net_BlockIterator;
typedef SORTED_SET_ITERATOR Net_GroupIterator;
typedef SORTED_SET_ITERATOR Net_ChainIterator;
typedef SORTED_SET_ITERATOR Net_FaceIterator;
typedef SORTED_SET_ITERATOR Net_ReferenceIterator;
typedef BDBCUR NetDisk_NetNameIterator;
typedef SORTED_SET_ITERATOR NetDisk_NetIterator;
typedef SORTED_SET_ITERATOR Reference_PseudoChromosomeIterator;
typedef SORTED_SET_ITERATOR PseudoChromsome_PseudoAdjacencyIterator;
typedef struct hashtable_itr Hash_Iterator;
typedef struct _face_FaceEndIterator Face_FaceEndIterator;
typedef struct _faceEndIterator FaceEnd_BottomNodeIterator;
typedef struct avl_table SortedSet;
typedef SORTED_SET_ITERATOR SortedSet_Iterator;

#endif
