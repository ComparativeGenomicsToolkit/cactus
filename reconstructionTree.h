#ifndef __RECONSTRUCTION_TREE_H__
#define __RECONSTRUCTION_TREE_H__

struct Sequence {
	char *contig;
	char *event;
	char *file;
	int32_t length;
};

struct Sequence *constructSequence(const char *contig, const char *event, const char *file, int32_t length);

void destructSequence(struct Sequence *sequence);

struct PinchGraph *constructPinchGraph(XMLNode xMainNode,
									   struct List *contigIndexToContigStrings,
									   struct IntList *contigIndexToContigStart);

struct hashtable *parseSequences(XMLNode xMainNode);

struct hashtable *createLeafEvents(XMLNode xMainNode);

int32_t createReconstructionProblem(
		const char *absolutePathPrefix, const char *relativeFilePath,
		struct List *subChains, struct List *adjacencyComponents, struct List *ends,
		struct PinchGraph *pinchGraph, struct hashtable *names,
		struct List *contigIndexToContigStrings, struct hashtable *contigNamesToSequences,
		const char *treeProgram, const char *tempFilePath);

void setUniqueNamePrefix(const char *uniqueNamePrefix);

void addBracketToAtomInstanceNamesInTree(XMLNode treeCopy, const char *name);

char *getAbsoluteFilePath(const char *absoluteFilePrefix, const char *relativeReconstructionProblem);

struct List *constructContigEventSets(XMLNode xMainNode);

//test functions

/*
 * This is the big one.
 */
void checkReconstructionTree(const char *absoluteFilePrefix, XMLNode xMainNode,
		int32_t checkChildrenRecursively, int32_t checkInternalAdjacencies);

#endif
