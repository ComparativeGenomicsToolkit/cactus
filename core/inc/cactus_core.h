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
    bool ignoreAllChainsLessThanMinimumTreeCoverage;
    int32_t minimumBlockLength;
    int32_t adjacencyComponentOverlap;
} CactusCoreInputParameters;

int32_t cactusCorePipeline(Flower *flower, CactusCoreInputParameters *cCIP,
        struct PairwiseAlignment *(*getNextAlignment)(),
        void(*startAlignmentStack)(), int32_t terminateRecursion);

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP);

CactusCoreInputParameters *constructCactusCoreInputParameters();

#endif
