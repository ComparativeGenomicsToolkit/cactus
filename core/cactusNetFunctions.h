#ifndef CACTUS_NET_FUNCTIONS_H_
#define CACTUS_NET_FUNCTIONS_H_

#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "pinchGraph.h"
#include "cactusGraph.h"
#include "cactus.h"

void fillOutNetFromInputs(Net *parentNet, struct CactusGraph *cactusGraph,
        struct PinchGraph *pinchGraph, stSortedSet *chosenBlocks);

struct PinchGraph *constructPinchGraph(Net *net);

void copyEndTreePhylogenies(Net *parentNet, Net *net);

#endif
