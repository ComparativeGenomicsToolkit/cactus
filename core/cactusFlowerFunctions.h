#ifndef CACTUS_FLOWER_FUNCTIONS_H_
#define CACTUS_FLOWER_FUNCTIONS_H_

#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "pinchGraph.h"
#include "cactusGraph.h"
#include "cactus.h"

void fillOutFlowerFromInputs(Flower *parentFlower, struct CactusGraph *cactusGraph,
        struct PinchGraph *pinchGraph, stSortedSet *chosenBlocks);

struct PinchGraph *constructPinchGraph(Flower *flower);

void copyEndTreePhylogenies(Flower *parentFlower, Flower *flower);

#endif
