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

typedef stSortedSetIterator EventTree_Iterator;
typedef struct _end_instanceIterator End_InstanceIterator;
typedef struct _block_instanceIterator Block_InstanceIterator;
typedef stSortedSetIterator Group_EndIterator;
typedef stSortedSetIterator Net_SequenceIterator;
typedef stSortedSetIterator Net_CapIterator;
typedef stSortedSetIterator Net_SegmentIterator;
typedef stSortedSetIterator Net_EndIterator;
typedef stSortedSetIterator Net_BlockIterator;
typedef stSortedSetIterator Net_GroupIterator;
typedef stSortedSetIterator Net_ChainIterator;
typedef stSortedSetIterator Net_FaceIterator;
typedef stSortedSetIterator Net_ReferenceIterator;
typedef BDBCUR NetDisk_NetNameIterator;
typedef stSortedSetIterator NetDisk_NetIterator;
typedef stSortedSetIterator Reference_PseudoChromosomeIterator;
typedef stSortedSetIterator PseudoChromsome_PseudoAdjacencyIterator;
typedef struct _face_FaceEndIterator Face_FaceEndIterator;
typedef struct _faceEndIterator FaceEnd_BottomNodeIterator;


#endif
