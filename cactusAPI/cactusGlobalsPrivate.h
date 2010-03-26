#ifndef CACTUS_GLOBALS_PRIVATE_H_
#define CACTUS_GLOBALS_PRIVATE_H_

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <tcutil.h>
#include <tcbdb.h>
#include "CuTest.h"
#include "commonC.h"
#include "avl.h"
#include "hashTableC.h"

#define NAME_STRING "%llX"

#include "cactusGroup.h"
#include "cactusGroupPrivate.h"
#include "cactusBlock.h"
#include "cactusSegment.h"
#include "cactusSegmentPrivate.h"
#include "cactusBlockPrivate.h"
#include "cactusChain.h"
#include "cactusChainPrivate.h"
#include "cactusDatabase.h"
#include "cactusEnd.h"
#include "cactusCap.h"
#include "cactusCapPrivate.h"
#include "cactusEndPrivate.h"
#include "cactusEvent.h"
#include "cactusEventPrivate.h"
#include "cactusEventTree.h"
#include "cactusEventTreePrivate.h"
#include "cactusGlobals.h"
#include "cactusLink.h"
#include "cactusLinkPrivate.h"
#include "cactusMetaEvent.h"
#include "cactusMetaEventPrivate.h"
#include "cactusMetaSequence.h"
#include "cactusMetaSequencePrivate.h"
#include "cactusNet.h"
#include "cactusNetDisk.h"
#include "cactusNetDiskPrivate.h"
#include "cactusNetMisc.h"
#include "cactusNetPrivate.h"
#include "cactusPseudoChromosome.h"
#include "cactusPseudoChromosomePrivate.h"
#include "cactusPseudoAdjacency.h"
#include "cactusPseudoAdjacencyPrivate.h"
#include "cactusFace.h"
#include "cactusFacePrivate.h"
#include "cactusReference.h"
#include "cactusReferencePrivate.h"
#include "cactusSequence.h"
#include "cactusSequencePrivate.h"
#include "cactusSerialisation.h"
#include "cactusSortedSet.h"
#include "cactusTestCommon.h"



#endif
