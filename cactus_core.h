#ifndef CACTUS_CORE_H_
#define CACTUS_CORE_H_

typedef struct _CactusCoreInputParameters {
	/*
	 * Arguments/options
	 */
	int32_t extensionSteps;
	int32_t maxEdgeDegree;
	bool writeDebugFiles;
	float minimumTreeCoverage;
	float minimumTreeCoverageForAtoms;
	int32_t minimumAtomLength;
	int32_t minimumChainLength;
	int32_t trim;
	int32_t alignRepeats;
} CactusCoreInputParameters;

int32_t cactusCorePipeline(Net *net, CactusCoreInputParameters *cCIP,
		struct PairwiseAlignment *(*getNextAlignment)());

void destructCactusCoreInputParameters(CactusCoreInputParameters *cCIP);

CactusCoreInputParameters *constructCactusCoreInputParameters();

#endif
