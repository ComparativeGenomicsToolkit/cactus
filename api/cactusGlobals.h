#ifndef CACTUS_GLOBALS_H_
#define CACTUS_GLOBALS_H_

#include <inttypes.h>

/*
 * Includes for Tokyo Cabinet.
 */
#include <tcutil.h>
#include <tcbdb.h>
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

typedef SortedSet_Iterator EventTree_Iterator;
typedef struct _end_instanceIterator End_InstanceIterator;
typedef struct _block_instanceIterator Block_InstanceIterator;
typedef SortedSet_Iterator Group_EndIterator;
typedef SortedSet_Iterator Net_SequenceIterator;
typedef SortedSet_Iterator Net_CapIterator;
typedef SortedSet_Iterator Net_SegmentIterator;
typedef SortedSet_Iterator Net_EndIterator;
typedef SortedSet_Iterator Net_BlockIterator;
typedef SortedSet_Iterator Net_GroupIterator;
typedef SortedSet_Iterator Net_ChainIterator;
typedef SortedSet_Iterator Net_FaceIterator;
typedef SortedSet_Iterator Net_ReferenceIterator;
typedef BDBCUR NetDisk_NetNameIterator;
typedef SortedSet_Iterator NetDisk_NetIterator;
typedef SortedSet_Iterator Reference_PseudoChromosomeIterator;
typedef SortedSet_Iterator PseudoChromsome_PseudoAdjacencyIterator;
typedef struct _face_FaceEndIterator Face_FaceEndIterator;
typedef struct _faceEndIterator FaceEnd_BottomNodeIterator;


#endif
