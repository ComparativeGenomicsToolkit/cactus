/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_CORE_H_
#define CACTUS_CORE_H_

#include "pairwiseAlignment.h"

typedef struct _CactusCoreInputParameters {
    /*
     * Arguments/options
     */
    bool writeDebugFiles;
    int32_t *annealingRounds;
    int32_t annealingRoundsLength;
    int32_t *deannealingRounds;
    int32_t deannealingRoundsLength;
    int32_t alignRepeatsAtRound;
    int32_t * trim;
    int32_t trimLength;

    /* Stuff for selecting chains to keep */
    float minimumTreeCoverage;
    int32_t blockTrim;
    int32_t minimumDegree;

    float requiredIngroupFraction;
    float requiredOutgroupFraction;
    float requiredAllFraction;
    int32_t requiredIngroups;
    int32_t requiredOutgroups;
    int32_t requiredAll;

    bool singleCopyIngroup;
    bool singleCopyOutgroup;
} CactusCoreInputParameters;

int32_t cactusCorePipeline(Flower *flower, CactusCoreInputParameters *cCIP,
        struct PairwiseAlignment *(*getNextAlignment)(void),
        void(*startAlignmentStack)(void), void (*cleanUpAlignment)(struct PairwiseAlignment *));

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP);

CactusCoreInputParameters *constructCactusCoreInputParameters();


#endif
