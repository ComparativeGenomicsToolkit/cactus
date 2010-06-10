#ifndef CACTUS_CORE_H_
#define CACTUS_CORE_H_

typedef struct _CactusCoreInputParameters {
    /*
     * Arguments/options
     */
    bool writeDebugFiles;
    int32_t annealingRounds;
    int32_t alignRepeatsAtRound;
    /* Stuff for adding more homologies to graph progressively */
    int32_t trim;
    float trimChange;
    /* Stuff for selecting chains to keep */
    float minimumTreeCoverage;
    int32_t minimumBlockLength;
    float minimumBlockLengthChange;
    int32_t minimumChainLength;
    float minimumChainLengthChange;
    int32_t deannealingRounds;
} CactusCoreInputParameters;

int32_t cactusCorePipeline(Net *net, CactusCoreInputParameters *cCIP,
        struct PairwiseAlignment *(*getNextAlignment)(),
        void(*startAlignmentStack)(), int32_t terminateRecursion);

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP);

CactusCoreInputParameters *constructCactusCoreInputParameters();

#endif
