#ifndef CACTUS_CORE_H_
#define CACTUS_CORE_H_

typedef struct _CactusCoreInputParameters {
    /*
     * Arguments/options
     */
    bool writeDebugFiles;
    int32_t alignUndoLoops;
    int32_t alignRepeatsAtLoop;
    /* Stuff for adding more homologies to graph progressively */
    int32_t trim;
    int32_t trimReduction;
    /* Stuff for selecting chains to keep */
    float minimumTreeCoverage;
    int32_t minimumBlockLength;
    int32_t minimumBlockLengthIncrease;
    int32_t minimumChainLength;
    int32_t minimumChainLengthIncrease;
    int32_t minimumChainLengthCactusUndoLoopStepSize;
} CactusCoreInputParameters;

int32_t cactusCorePipeline(Net *net, CactusCoreInputParameters *cCIP,
        struct PairwiseAlignment *(*getNextAlignment)(),
        void(*startAlignmentStack)());

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP);

CactusCoreInputParameters *constructCactusCoreInputParameters();

#endif
