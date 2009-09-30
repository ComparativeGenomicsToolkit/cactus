#ifndef CACTUS_NET_FUNCTIONS_H_
#define CACTUS_NET_FUNCTIONS_H_

#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "pinchGraph.h"
#include "cactusGraph.h"
#include "net.h"

void fillOutNetFromInputs(Net *parentNet, struct CactusGraph *cactusGraph,
		struct PinchGraph *pinchGraph, const char *uniqueNamePrefix, struct List *chosenAtoms,
		struct List *contigIndexToContigStrings);

struct PinchGraph *constructPinchGraph(Net *net,
		struct List *contigIndexToContigStrings,
		struct IntList *contigIndexToContigStart);

void copyEndTreePhylogenies(Net *parentNet, Net *net);

#endif
