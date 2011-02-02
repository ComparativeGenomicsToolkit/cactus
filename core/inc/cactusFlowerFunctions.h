/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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
        struct PinchGraph *pinchGraph, stList *adjacencyComponents);

struct PinchGraph *constructPinchGraph(Flower *flower);

void copyEndTreePhylogenies(Flower *parentFlower, Flower *flower);

stHash *buildOrientationHash(struct List *biConnectedComponents,
        struct PinchGraph *pinchGraph, Flower *parentFlower,
        struct hashtable *endNamesHash);

#endif
