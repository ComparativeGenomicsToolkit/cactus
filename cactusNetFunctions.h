#ifndef CACTUS_NET_FUNCTIONS_H_
#define CACTUS_NET_FUNCTIONS_H_

#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "pinchGraph.h"
#include "cactusGraph.h"
#include "net.h"

Net *constructNetFromInputs(struct CactusGraph *cactusGraph,
		struct PinchGraph *pinchGraph, struct hashtable *names, struct List *chosenAtoms,
		struct hashtable *contigStringsToSequences, struct List *contigIndexToContigStrings,
		NetDisk *netDisk, const char *(*getUniqueName)());

struct PinchGraph *constructPinchGraph(Net *net,
		struct List *contigIndexToContigStrings,
		struct IntList *contigIndexToContigStart);

#endif
