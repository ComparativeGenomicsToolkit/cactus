#ifndef NET_H_
#define NET_H_

#include "fastCMaths.h"
#include "commonC.h"
#include "hashTableC.h"
#include "avl.h"

#include "cactusGraph.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Chains
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct Chain {
	/*
	 * A doubly linked list which can hold sub chains (to form the the net structure).
	 */
	//edge associated with the link in the chain
	struct List *segments;
	//previous link in the chain.
	struct Chain *pLink;
	//next link in the chain.
	struct Chain *nLink;
	//the sub-chains used to form the net structure.
	struct List *subChains;
	//score for chain element
	float score;
	//list of adjacency components
	struct List *adjacencyComponents;
	//list of ends, including the ends of the link and the stubs nested inside of this chain.
	struct List *ends;
};

struct Chain *constructChain(struct List *segments, struct Chain *pChain);

void destructChain(struct Chain *link);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods on chains
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void sortChains(struct List *listOfChains);

void flattenChainList(struct List *chainList, struct List *flatList);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Nets
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct Net {
	struct List *chains;
	struct List *adjacencyComponents;
	struct List *ends;
};

struct Net *constructNet(struct CactusGraph *cactusGraph);

void destructNet(struct Net *net);

void checkNet(struct CactusGraph *graph, struct Net *net);

void writeOutNet(struct Net *net, struct hashtable *names, FILE *fileHandle);

void pruneNet(struct Net *net, struct List *chosenAtoms);

struct List *identifyPseudoAdjacencyComponents(struct List *chosenAtoms, struct PinchGraph *pinchGraph);

void addAdjacencyComponents(struct Net *net, struct PinchGraph *pinchGraph, struct List *adjacencyComponents);

void removeStubAdjacencyComponents(struct Net *net);

void addEnds(struct Net *net);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Choosing which atoms in the chains to keep.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct ContigEventSet {
	struct hashtable *contigs;
	char *event;
	float branchLength;
};

struct ContigEventSet *constructContigEventSet();

void destructContigEventSet(struct ContigEventSet *contigEventSet);

void filterAtomsByTreeCoverageAndLength(struct List *chainsList,
		struct List *chosenAtoms,
		struct List *nodeSets,
		float proportionToKeep, /*Proportion of all atoms to select to keep*/
		float discardRatio, /*The proportion of an atom's chain's average atom score required to be score to be considered */
		float minimumTreeCoverage, /*Minimum tree coverage to be included */
		int32_t minimumChainLength, /* Minimum chain length to be included */
		struct PinchGraph *pinchGraph,
		struct List *contigIndexToContigStrings);

void filterAtomsByIfStubOrCap(struct List *chainsList,
		struct List *chosenAtoms, struct PinchGraph *pinchGraph);

void logTheChosenAtomSubset(struct List *allAtoms, struct List *chosenAtoms, struct PinchGraph *pinchGraph,
		struct List *nodeSets, struct List *contigIndexToContigStrings);

#endif
