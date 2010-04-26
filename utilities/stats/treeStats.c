#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"


void tabulateFloatStats(struct List *unsortedValues, double *totalNumber, double *totalSum, double *min, double *max, double *avg, double *median) {
	/*
	 * Calculates basic stats from a list of float values.
	 */
	if(unsortedValues->length == 0) {
		*totalNumber = 0;
		*min = INT32_MAX;
		*max = INT32_MAX;
		*avg = INT32_MAX;
		*median = INT32_MAX;
		return;
	}
	unsortedValues = listCopy(unsortedValues); //copy the input list, to avoid altering the input.
	assert(unsortedValues->length > 0);
	qsort(unsortedValues->list, unsortedValues->length, sizeof(void *), (int (*)(const void *, const void *))floatComparator);
	*totalNumber = unsortedValues->length;
	*min = *(float *)unsortedValues->list[0];
	*max = *(float *)unsortedValues->list[unsortedValues->length-1];
	*median = *(float *)unsortedValues->list[unsortedValues->length/2];
	int32_t i;
	float j = 0;
	for(i=0; i<unsortedValues->length; i++) {
		j += *(float *)unsortedValues->list[i];
	}
	*avg = j / unsortedValues->length;
	*totalSum = j;
	unsortedValues->destructElement = NULL;
	destructList(unsortedValues);
}

void tabulateStats(struct IntList *unsortedValues, double *totalNumber, double *totalSum, double *min, double *max, double *avg, double *median) {
	/*
	 * Same as float stats, but for an intlist.
	 */
	if(unsortedValues->length == 0) {
		*totalNumber = 0;
		*min = INT32_MAX;
		*max = INT32_MAX;
		*avg = INT32_MAX;
		*median = INT32_MAX;
		return;
	}
	unsortedValues = intListCopy(unsortedValues);
	assert(unsortedValues->length > 0);
	qsort(unsortedValues->list, unsortedValues->length, sizeof(int32_t), (int (*)(const void *, const void *))intComparator_Int);
	*totalNumber = unsortedValues->length;
	*min = unsortedValues->list[0];
	*max = unsortedValues->list[unsortedValues->length-1];
	*median = unsortedValues->list[unsortedValues->length/2];
	int32_t i, j = 0;
	for(i=0; i<unsortedValues->length; i++) {
		j += unsortedValues->list[i];
	}
	*avg = (double)j / unsortedValues->length;
	*totalSum = j;
	destructIntList(unsortedValues);
}

void printOpeningTag(const char *tag, FILE *fileHandle) {
	/*
	 * Creates an opening XML tag.
	 */
	fprintf(fileHandle, "<%s>", tag);
}

void printClosingTag(const char *tag, FILE *fileHandle) {
	/*
	 * Creates a closing XML tag.
	 */
	fprintf(fileHandle, "</%s>", tag);
}

void tabulateAndPrintFloatValues(struct List *values, const char *tag, FILE *fileHandle) {
	/*
	 * Creates a node containing basic stats on the given float values and a nested "values" node containing the actual values.
	 */
	double totalNumber, totalSum, min, max, avg, median;
	tabulateFloatStats(values, &totalNumber, &totalSum, &min, &max, &avg, &median);
	fprintf(fileHandle, "<%s total=\"%f\" sum=\"%f\" min=\"%f\" max=\"%f\" avg=\"%f\" median=\"%f\">", tag, totalNumber, totalSum, min, max, avg, median);
	int32_t i;
	for(i=0; i<values->length; i++) {
		fprintf(fileHandle, "%f ", *(float *)values->list[i]);
	}
	printClosingTag(tag, fileHandle);
}

void tabulateAndPrintIntValues(struct IntList *values, const char *tag, FILE *fileHandle) {
	/*
	 * Creates a node containing basic stats on the given int values and a nested "values" node containing the actual values.
	 */
	double totalNumber, totalSum, min, max, avg, median;
	tabulateStats(values, &totalNumber, &totalSum, &min, &max, &avg, &median);
	fprintf(fileHandle, "<%s total=\"%f\" sum=\"%f\" min=\"%f\" max=\"%f\" avg=\"%f\" median=\"%f\">", tag, totalNumber, totalSum, min, max, avg, median);
	int32_t i;
	for(i=0; i<values->length; i++) {
		fprintf(fileHandle, "%i ", values->list[i]);
	}
	printClosingTag(tag, fileHandle);
}

double calculateTreeBits(Net *net, double pathBitScore) {
	/*
	 * Calculates the total number of bits to required to encode the path to every base in the net.
	 */
	double totalBitScore = 0.0;
	int32_t totalSequenceSize;
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	double followingPathBitScore = (log(net_getGroupNumber(net)) / log(2.0)) + pathBitScore;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_isTerminal(group)) {
			totalSequenceSize = group_getTotalBaseLength(group);
			totalBitScore += (totalSequenceSize > 0 ? ((log(totalSequenceSize) / log(2.0)) + followingPathBitScore) * totalSequenceSize : 0.0);
		}
		else {
			totalBitScore += calculateTreeBits(group_getNestedNet(group), followingPathBitScore);
		}
	}
	net_destructGroupIterator(groupIterator);
	Net_BlockIterator *blockIterator = net_getBlockIterator(net);
	Block *block;
	totalSequenceSize = 0.0;
	while((block = net_getNextBlock(blockIterator)) != NULL) {
		totalSequenceSize += block_getLength(block) * block_getInstanceNumber(block);
	}
	net_destructBlockIterator(blockIterator);
	return totalBitScore + (totalSequenceSize > 0 ? ((log(totalSequenceSize) / log(2.0)) + pathBitScore) * totalSequenceSize : 0.0);
}

void reportRelativeEntopyStats(Net *net, FILE *fileHandle) {
	/*
	 * Relative entropy stats. Supposed to give a metric of how balanced the tree is in how it subdivides the input sequences.
	 */
	double totalSeqSize = net_getTotalBaseLength(net);
	double totalP = calculateTreeBits(net, 0.0);
	double totalQ = (log(totalSeqSize) / log(2.0)) * totalSeqSize;
	//assert(totalP >= totalQ);
	double relativeEntropy = totalP - totalQ;
	double normalisedRelativeEntropy = relativeEntropy / totalSeqSize;

	fprintf(fileHandle, "<relative_entropy_stats totalP=\"%f\" totalQ=\"%f\" relative_entropy=\"%f\" normalised_relative_entropy=\"%f\"/>", totalP, totalQ, relativeEntropy, normalisedRelativeEntropy);
}

static void netStats(Net *net, int32_t currentDepth, struct IntList *children, struct IntList *tangleChildren, struct IntList *linkChildren, struct IntList *depths) {
	/*
	 * Calculates basic stats on nets.
	 * Children is the number of children internal nodes (those with children), have.
	 * Tangle children, like children but only including groups that are tangle groups.
	 * Link children, like children but only including groups that are link groups.
	 * Depth is the length of a path (in terms of edges/connections) from the root net to a terminal net (which are the imaginary nets that would descend from terminal groups of the tree).
	 */
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	int32_t i = 0;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_isTerminal(group)) {
			intListAppend(depths, currentDepth+1);
		}
		else {
			netStats(group_getNestedNet(group), currentDepth+1, children, tangleChildren, linkChildren, depths);
		}
		if(group_getLink(group) != NULL) {
			i++;
		}
	}
	net_destructGroupIterator(groupIterator);
	intListAppend(children, net_getGroupNumber(net));
	intListAppend(tangleChildren, net_getGroupNumber(net) - i);
	intListAppend(linkChildren, i);
}

void reportNetStats(Net *net, FILE *fileHandle) {
	/*
	 * Prints the chain stats to the XML file.
	 */
	struct IntList *children = constructEmptyIntList(0);
	struct IntList *tangleChildren = constructEmptyIntList(0);
	struct IntList *linkChildren = constructEmptyIntList(0);
	struct IntList *depths = constructEmptyIntList(0);
	netStats(net, 0, children, tangleChildren, linkChildren, depths);
	printOpeningTag("nets", fileHandle);
	tabulateAndPrintIntValues(children, "children", fileHandle);
	tabulateAndPrintIntValues(tangleChildren, "tangle_children", fileHandle);
	tabulateAndPrintIntValues(linkChildren, "link_children", fileHandle);
	tabulateAndPrintIntValues(depths, "depths", fileHandle);
	printClosingTag("nets", fileHandle);
	destructIntList(children);
	destructIntList(tangleChildren);
	destructIntList(linkChildren);
	destructIntList(depths);
}

void blockStats(Net *net, struct IntList *counts, struct IntList *lengths, struct IntList *degrees,
		struct IntList *leafDegrees, struct IntList *coverage, struct IntList *leafCoverage, int32_t minLeafDegree) {
	/*
	 * Calculates stats on the blocks.
	 * Counts is numbers of blocks per net.
	 * Lengths is lengths of blocks.
	 * Degrees is the number of segment instances in each block.
	 * Leaf degree is the number of leaf segment instances in each block.
	 * Coverage is the length * degree of each block.
	 * Leaf coverage is the length * leadf degree of each block.
	 */
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(!group_isTerminal(group)) {
			blockStats(group_getNestedNet(group), counts, lengths, degrees, leafDegrees, coverage, leafCoverage, minLeafDegree);
		}
	}
	net_destructGroupIterator(groupIterator);

	Net_BlockIterator *blockIterator = net_getBlockIterator(net);
	Block *block;
	while((block = net_getNextBlock(blockIterator)) != NULL) {
		Segment *segment;
		Block_InstanceIterator *segmentIterator = block_getInstanceIterator(block);
		int32_t i=0;
		while((segment = block_getNext(segmentIterator)) != NULL) {
			if(segment_getChildNumber(segment) == 0) {
				i++;
			}
		}
		block_destructInstanceIterator(segmentIterator);
		if(i >= minLeafDegree) {
			intListAppend(lengths, block_getLength(block));
			intListAppend(degrees, block_getInstanceNumber(block));
			intListAppend(coverage, block_getLength(block)*block_getInstanceNumber(block));
			intListAppend(leafDegrees, i);
			intListAppend(leafCoverage, block_getLength(block)*i);
		}
	}
	net_destructBlockIterator(blockIterator);
	intListAppend(counts, net_getBlockNumber(net));
}

void reportBlockStats(Net *net, FILE *fileHandle, int32_t minLeafDegree) {
	/*
	 * Prints the block stats to the XML file.
	 */
	struct IntList *counts = constructEmptyIntList(0);
	struct IntList *lengths = constructEmptyIntList(0);
	struct IntList *degrees = constructEmptyIntList(0);
	struct IntList *leafDegrees = constructEmptyIntList(0);
	struct IntList *coverage = constructEmptyIntList(0);
	struct IntList *leafCoverage = constructEmptyIntList(0);
	blockStats(net, counts, lengths, degrees, leafDegrees, coverage, leafCoverage, minLeafDegree);
	fprintf(fileHandle, "<blocks minimum_leaf_degree=\"%i\">", minLeafDegree);
	tabulateAndPrintIntValues(counts, "counts", fileHandle);
	tabulateAndPrintIntValues(lengths, "lengths", fileHandle);
	tabulateAndPrintIntValues(degrees, "degrees", fileHandle);
	tabulateAndPrintIntValues(leafDegrees, "leaf_degrees", fileHandle);
	tabulateAndPrintIntValues(coverage, "coverage", fileHandle);
	tabulateAndPrintIntValues(leafCoverage, "leaf_coverage", fileHandle);
	printClosingTag("blocks", fileHandle);
	destructIntList(counts);
	destructIntList(lengths);
	destructIntList(degrees);
	destructIntList(leafDegrees);
	destructIntList(coverage);
	destructIntList(leafCoverage);
}

static void chainStats(Net *net, struct IntList *counts, struct IntList *blockNumbers,
		struct IntList *baseBlockLengths, struct IntList *linkNumbers,
		struct IntList *avgInstanceBaseLengths,
		int32_t minNumberOfBlocksInChain) {
	/*
	 * Gets stats on the chains.
	 * Counts is numbers per net.
	 * Block number is the number of blocks per chain.
	 * Base block lengths in the number of basepairs in blocks per chain.
	 * Link numbers if the number of links per chain.
	 * Avg instance base lengths is the avg number of basepairs in an instance of a chain, per chain.
	 */
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			chainStats(group_getNestedNet(group), counts, blockNumbers, baseBlockLengths, linkNumbers, avgInstanceBaseLengths,
					minNumberOfBlocksInChain);
		}
	}
	net_destructGroupIterator(groupIterator);

	Net_ChainIterator *chainIterator = net_getChainIterator(net);
	Chain *chain;
	Block **blocks;
	int32_t i, j, k, l;
	l = 0;
	while((chain = net_getNextChain(chainIterator)) != NULL) {
		blocks = chain_getBlockChain(chain, &i);
		k = 0;
		for(j=0; j<i; j++) {
			k += block_getLength(blocks[j]);
		}

		/*Chain stats are only for those containing two or more blocks.*/
		if(i >= minNumberOfBlocksInChain) {
			intListAppend(blockNumbers, i);
			intListAppend(baseBlockLengths, k);
			intListAppend(linkNumbers, chain_getLength(chain));
			intListAppend(avgInstanceBaseLengths, chain_getAverageInstanceBaseLength(chain));
			l++;
		}
	}
	net_destructBlockIterator(chainIterator);
	intListAppend(counts, l);
}

void reportChainStats(Net *net,
					 int32_t minNumberOfBlocksInChain, FILE *fileHandle) {
	/*
	 * Prints the chain stats to the XML file.
	 */
	struct IntList *counts = constructEmptyIntList(0);
	struct IntList *blockNumbers = constructEmptyIntList(0);
	struct IntList *baseBlockLengths = constructEmptyIntList(0);
	struct IntList *linkNumbers = constructEmptyIntList(0);
	struct IntList *avgInstanceBaseLengths = constructEmptyIntList(0);
	chainStats(net, counts, blockNumbers, baseBlockLengths, linkNumbers, avgInstanceBaseLengths, minNumberOfBlocksInChain);
	fprintf(fileHandle, "<chains minimum_number_of_blocks_in_chain=\"%i\">", minNumberOfBlocksInChain);
	tabulateAndPrintIntValues(counts, "counts", fileHandle);
	tabulateAndPrintIntValues(blockNumbers, "block_numbers", fileHandle);
	tabulateAndPrintIntValues(baseBlockLengths, "base_block_lengths", fileHandle);
	tabulateAndPrintIntValues(linkNumbers, "link_numbers", fileHandle);
	tabulateAndPrintIntValues(avgInstanceBaseLengths, "avg_instance_base_length", fileHandle);
	printClosingTag("chains", fileHandle);
	destructIntList(counts);
	destructIntList(blockNumbers);
	destructIntList(baseBlockLengths);
	destructIntList(linkNumbers);
	destructIntList(avgInstanceBaseLengths);
}

void endStats(Net *net, struct IntList *counts, struct IntList *degrees,
		int32_t includeLinkGroups, int32_t includeTangleGroups,
		int32_t includeTerminalGroups, int32_t includeNonTerminalGroups) {
	/*
	 * Calculates stats on the ends in nets.
	 * Counts is the number of ends per net.
	 * Degrees is the degree of each end, where the degree of an end is the number
	 * of distinct adjacencies an end has to other ends.
	 */
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(  ((includeLinkGroups && group_getLink(group) != NULL) || (includeTangleGroups && group_getLink(group) == NULL)) &&
			 ((includeTerminalGroups && group_isTerminal(group)) || (includeNonTerminalGroups && !group_isTerminal(group))) ) {
			intListAppend(counts, group_getEndNumber(group));
			Group_EndIterator *endIterator = group_getEndIterator(group);
			End *end;
			while((end = group_getNextEnd(endIterator)) != NULL) {
				struct List *list = constructEmptyList(0, NULL);
				End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
				Cap *cap;
				while((cap = end_getNext(instanceIterator)) != NULL) {
					Cap *cap2 = cap_getAdjacency(cap);
					if(cap2 != NULL) {
						End *end =  end_getPositiveOrientation(cap_getEnd(cap2));
						if(!listContains(list, end)) {
							listAppend(list, end);
						}
					}
				}
				end_destructInstanceIterator(instanceIterator);
				intListAppend(degrees, list->length);
				destructList(list);
			}
			group_destructEndIterator(endIterator);
		}
		if(!group_isTerminal(group)) {
			endStats(group_getNestedNet(group), counts, degrees,
					includeLinkGroups, includeTangleGroups,
					includeTerminalGroups, includeNonTerminalGroups);
		}

	}
	net_destructGroupIterator(groupIterator);
}

void reportEndStats(Net *net, int32_t includeLinkGroups, int32_t includeTangleGroups,
		                      int32_t includeTerminalGroups, int32_t includeNonTerminalGroups, FILE *fileHandle) {
	/*
	 * Prints the end stats to the XML file.
	 */
	struct IntList *counts = constructEmptyIntList(0);
	struct IntList *degrees = constructEmptyIntList(0);
	endStats(net, counts, degrees, includeLinkGroups, includeTangleGroups,
			includeTerminalGroups, includeNonTerminalGroups);
	fprintf(fileHandle, "<ends include_link_groups=\"%i\" include_tangle_groups=\"%i\" include_terminal_groups=\"%i\" include_non_terminal_groups=\"%i\">",
			includeLinkGroups != 0, includeTangleGroups != 0, includeTerminalGroups != 0, includeNonTerminalGroups != 0);
	tabulateAndPrintIntValues(counts, "counts", fileHandle);
	tabulateAndPrintIntValues(degrees, "degrees", fileHandle);
	printClosingTag("ends", fileHandle);
	destructIntList(counts);
	destructIntList(degrees);
}

void terminalGroupsSizes(Net *net, struct IntList *sizes) {
	/*
	 * Reports stats on the terminal groups.
	 * Sizes = This gives the sizes of the terminal groups, i.e. the number of bases in adjacencies between ends in terminal groups.
	 * If the cactus tree has been fully decomposed then all terminal groups will contain 0 bases.
	 */
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(!group_isTerminal(group)) {
			terminalGroupsSizes(group_getNestedNet(group), sizes);
		}
		else {
			intListAppend(sizes, (int32_t)group_getTotalBaseLength(group));
		}
	}
	net_destructGroupIterator(groupIterator);
}

void reportTerminalGroupsSizes(Net *net, FILE *fileHandle) {
	/*
	 * Prints the terminal group size stats to the XML file.
	 */
	struct IntList *sizes = constructEmptyIntList(0);
	terminalGroupsSizes(net, sizes);
	tabulateAndPrintIntValues(sizes, "terminal_group_sizes", fileHandle);
	destructIntList(sizes);
}

void faceStats(Net *net, struct IntList *numberPerNet, struct IntList *cardinality,
		struct IntList *isSimple, struct IntList *isRegular, struct IntList *isCanonical,
		struct IntList *facesPerFaceAssociatedEnd,
		int32_t includeLinkGroups, int32_t includeTangleGroups,
		int32_t includeTerminalGroups, int32_t includeNonTerminalGroups) {
	/*
	 * Face stats.
	 * Number per net: faces per net.
	 * Cardinality of face.
	 * isSimple: if face is simple.
	 * isRegular: is face is regular.
	 * isCanonical: if face is canonical.
	 * facesPerFaceAssociatedEnd: the number of faces associated with each end that
	 * is associated with at least one end. Used to calculate the breakpoint reuse ratio.
	 */
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		//Call recursively..
		if(!group_isTerminal(group)) {
			faceStats(group_getNestedNet(group),
					  numberPerNet, cardinality,
					  isSimple, isRegular,
					  isCanonical, facesPerFaceAssociatedEnd,
					  includeLinkGroups, includeTangleGroups,
                      includeTerminalGroups, includeNonTerminalGroups);
		}

		//Now the actual meat.
		if(  ((includeLinkGroups && group_getLink(group) != NULL) || (includeTangleGroups && group_getLink(group) == NULL)) &&
					 ((includeTerminalGroups && group_isTerminal(group)) || (includeNonTerminalGroups && !group_isTerminal(group))) ) {
			intListAppend(numberPerNet, net_getFaceNumber(net));
			Net_FaceIterator *faceIterator = net_getFaceIterator(net);
			Face *face;
			while((face = net_getNextFace(faceIterator)) != NULL) {
				intListAppend(cardinality, face_getCardinal(face));
				intListAppend(isSimple, face_isSimple(face));
				intListAppend(isRegular, face_isRegular(face));
				intListAppend(isCanonical, face_isCanonical(face));
			}
			net_destructFaceIterator(faceIterator);

			/*
			 * Calculate the number of faces associated with each end
			 * associated with at least one face.
			 */
			End *end;
			Net_EndIterator *endIterator = net_getFaceIterator(net);
			while((end = net_getNextEnd(endIterator)) != NULL) {
				End_InstanceIterator *capIterator = end_getInstanceIterator(end);
				Cap *cap;
				Hash *faceHash = hash_construct();
				int32_t i = 0;
				while((cap = end_getNext(capIterator)) != NULL) {
					Face *face = cap_getFace(cap);
					if(face != NULL) {
						if(hash_search(faceHash, face) == NULL) {
							hash_insert(faceHash, face, face);
							i++;
						}
					}
				}
				hash_destruct(faceHash);
				end_destructInstanceIterator(capIterator);
				if(i > 0) {
					intListAppend(facesPerFaceAssociatedEnd, i);
				}
			}
			net_destructEndIterator(endIterator);
		}
	}
	net_destructGroupIterator(groupIterator);
}

void reportFaceStats(Net *net,
		int32_t includeLinkGroups, int32_t includeTangleGroups,
		int32_t includeTerminalGroups, int32_t includeNonTerminalGroups, FILE *fileHandle) {
	/*
	 * Prints the reference stats to the XML file.
	 */
	struct IntList *numberPerNet = constructEmptyIntList(0);
	struct IntList *cardinality = constructEmptyIntList(0);
	struct IntList *isSimple = constructEmptyIntList(0);
	struct IntList *isRegular = constructEmptyIntList(0);
	struct IntList *isCanonical = constructEmptyIntList(0);
	struct IntList *facesPerFaceAssociatedEnd = constructEmptyIntList(0);
	faceStats(net, numberPerNet, cardinality, isSimple, isRegular, isCanonical,
			  facesPerFaceAssociatedEnd,
			  includeLinkGroups, includeTangleGroups,
			  includeTerminalGroups, includeNonTerminalGroups);
	fprintf(fileHandle, "<faces include_link_groups=\"%i\" include_tangle_groups=\"%i\" include_terminal_groups=\"%i\" include_non_terminal_groups=\"%i\">",
				includeLinkGroups != 0, includeTangleGroups != 0, includeTerminalGroups != 0, includeNonTerminalGroups != 0);
	tabulateAndPrintIntValues(numberPerNet, "number_per_net", fileHandle);
	tabulateAndPrintIntValues(cardinality, "cardinality", fileHandle);
	tabulateAndPrintIntValues(isSimple, "is_simple", fileHandle);
	tabulateAndPrintIntValues(isRegular, "is_regular", fileHandle);
	tabulateAndPrintIntValues(isCanonical, "is_canonical", fileHandle);
	tabulateAndPrintIntValues(facesPerFaceAssociatedEnd, "faces_per_face_associated_end", fileHandle);
	printClosingTag("faces", fileHandle);
	destructIntList(numberPerNet);
	destructIntList(cardinality);
	destructIntList(isSimple);
	destructIntList(isRegular);
	destructIntList(isCanonical);
}

void referenceStats(Net *net, struct IntList *pseudoChromosomeNumber,
				    struct IntList *pseudoAdjacencyNumberPerChromosome,
				    struct IntList *truePseudoAdjacencyNumberPerChromosome,
				    struct IntList *linksPerChromosome) {
	/*
	 * Calculates stats on the reference genome structure.
	 * Stats are pretty obvious.
	 */
	//Call recursively..
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(!group_isTerminal(group)) {
			referenceStats(group_getNestedNet(group),
							pseudoChromosomeNumber,
							pseudoAdjacencyNumberPerChromosome,
							truePseudoAdjacencyNumberPerChromosome,
							linksPerChromosome);
		}
	}
	net_destructGroupIterator(groupIterator);

	//Calculate stats for first reference.
	Reference *reference = net_getFirstReference(net);
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator = reference_getPseudoChromosomeIterator(reference);
	PseudoChromosome *pseudoChromosome;
	intListAppend(pseudoChromosomeNumber, reference_getPseudoChromosomeNumber(reference));
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		intListAppend(pseudoAdjacencyNumberPerChromosome, pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome));
		PseudoChromsome_PseudoAdjacencyIterator *adjacencyIterator = pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome);
		PseudoAdjacency *pseudoAdjacency;
		int32_t i = 0, j = 0;
		while((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(adjacencyIterator)) != NULL) {
			Group *group = end_getGroup(pseudoAdjacency_get5End(pseudoAdjacency));
			assert(group != NULL);
			if(group_getLink(group) != NULL) {
				i++;
			}
			End *_5End = pseudoAdjacency_get5End(pseudoAdjacency);
			End *_3End = pseudoAdjacency_get5End(pseudoAdjacency);
			Cap *cap;
			int32_t k = 0;
			End_InstanceIterator *instanceIterator = end_getInstanceIterator(_5End);
			while((cap = end_getNext(instanceIterator)) != NULL) {
				Cap *adjacentCap = cap_getAdjacency(cap);
				if(adjacentCap != NULL) {
					if(end_getPositiveOrientation(cap_getEnd(adjacentCap)) == _3End) {
						k = 1;
					}
				}
			}
			if(k) {
				j++;
			}
		}
		intListAppend(linksPerChromosome, i);
		intListAppend(truePseudoAdjacencyNumberPerChromosome, j);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
}

void reportReferenceStats(Net *net, FILE *fileHandle) {
	/*
	 * Prints the reference stats to the XML file.
	 */
	struct IntList *pseudoChromosomeNumber = constructEmptyIntList(0);
	struct IntList *pseudoAdjacencyNumberPerChromosome = constructEmptyIntList(0);
	struct IntList *truePseudoAdjacencyNumberPerChromosome = constructEmptyIntList(0);
	struct IntList *linksPerChromosome = constructEmptyIntList(0);
	referenceStats(net, pseudoChromosomeNumber, pseudoAdjacencyNumberPerChromosome,
			truePseudoAdjacencyNumberPerChromosome, linksPerChromosome);
	printOpeningTag("reference", fileHandle);
	tabulateAndPrintIntValues(pseudoChromosomeNumber, "pseudo_chromosome_number", fileHandle);
	tabulateAndPrintIntValues(pseudoAdjacencyNumberPerChromosome, "pseudo_adjacency_number_per_pseudo_chromosome", fileHandle);
	tabulateAndPrintIntValues(truePseudoAdjacencyNumberPerChromosome, "true_pseudo_adjacency_number_per_pseudo_chromosome", fileHandle);
	tabulateAndPrintIntValues(linksPerChromosome, "links_per_chromosome", fileHandle);
	printClosingTag("reference", fileHandle);
	destructIntList(pseudoChromosomeNumber);
	destructIntList(pseudoAdjacencyNumberPerChromosome);
	destructIntList(truePseudoAdjacencyNumberPerChromosome);
	destructIntList(linksPerChromosome);
}

void reportNetDiskStats(char *netDiskName, Net *net, FILE *fileHandle) {

	double totalSeqSize = net_getTotalBaseLength(net);
	fprintf(fileHandle, "<stats net_disk=\"%s\" net_name=\"%s\" total_sequence_length=\"%f\" >", netDiskName, netMisc_nameToStringStatic(net_getName(net)), totalSeqSize);

	/*
	 * Relative entropy numbers on the balance of the tree.
	 */
	reportRelativeEntopyStats(net, fileHandle);

	/*
	 * Numbers on the structure of the tree.
	 */
	reportNetStats(net, fileHandle);

	/*
	 * Numbers on the blocks.
	 */
	reportBlockStats(net, fileHandle, 0);
	reportBlockStats(net, fileHandle, 2);

	/*
	 * Chain statistics.
	 */
	reportChainStats(net, 0, fileHandle);
	reportChainStats(net, 2, fileHandle);

	/*
	 * Stats on the ends in the problem. We look at tangle and link groups separately and at terminal and non-terminal groups seperately.
	 */
	reportEndStats(net, 1, 1, 1, 1, fileHandle);
	reportEndStats(net, 0, 1, 1, 1, fileHandle);
	reportEndStats(net, 1, 0, 1, 1, fileHandle);
	reportEndStats(net, 1, 1, 1, 0, fileHandle);
	reportEndStats(net, 0, 1, 1, 0, fileHandle);
	reportEndStats(net, 1, 0, 1, 0, fileHandle);
	reportEndStats(net, 1, 1, 0, 1, fileHandle);
	reportEndStats(net, 0, 1, 0, 1, fileHandle);
	reportEndStats(net, 1, 0, 0, 1, fileHandle);

	/*
	 * Stats on terminal groups in the tree.
	 */
	reportTerminalGroupsSizes(net, fileHandle);

	/*
	 * Stats on faces in the reconstruction..
	 */
	reportFaceStats(net, 1, 1, 1, 1, fileHandle);
	reportFaceStats(net, 0, 1, 1, 1, fileHandle);
	reportFaceStats(net, 1, 0, 1, 1, fileHandle);
	reportFaceStats(net, 1, 1, 1, 0, fileHandle);
	reportFaceStats(net, 0, 1, 1, 0, fileHandle);
	reportFaceStats(net, 1, 0, 1, 0, fileHandle);
	reportFaceStats(net, 1, 1, 0, 1, fileHandle);
	reportFaceStats(net, 0, 1, 0, 1, fileHandle);
	reportFaceStats(net, 1, 0, 0, 1, fileHandle);

	/*
	 * Stats on the reference in the reconstruction..
	 */
	reportReferenceStats(net, fileHandle);

	printClosingTag("stats", fileHandle);
}

