#ifndef CACTUS_CORE_H_
#define CACTUS_CORE_H_

typedef struct _CactusCoreInputParameters {
	/*
	 * Arguments/options
	 */
	bool writeDebugFiles;
	int32_t alignRepeats;
	int32_t alignUndoLoops;
	int32_t maxEdgeDegree;
	int32_t extensionSteps;
	int32_t extensionStepsReduction;
	int32_t trim;
	int32_t trimReduction;
	float minimumTreeCoverage;
	float minimumTreeCoverageReduction;
	int32_t minimumBlockLength;
	int32_t minimumChainLength;
	int32_t minimumChainLengthReduction;
} CactusCoreInputParameters;

int32_t cactusCorePipeline(Net *net, CactusCoreInputParameters *cCIP,
		struct PairwiseAlignment *(*getNextAlignment)(),
		void (*startAlignmentStack)());

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP);

CactusCoreInputParameters *constructCactusCoreInputParameters();

#endif
