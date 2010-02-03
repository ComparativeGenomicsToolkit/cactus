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

void usage() {
	fprintf(stderr, "cactus_treeStats, version 0.1\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --outputFile : The file to write the stats in, XML formatted.\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

double calculateTreeBits(Net *net, double pathBitScore) {
	double totalBitScore = 0.0;
	int32_t totalSequenceSize;
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	double followingPathBitScore = (log(net_getGroupNumber(net)) / log(2.0)) + pathBitScore;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			totalBitScore += calculateTreeBits(group_getNestedNet(group), followingPathBitScore);
		}
		else {
			totalSequenceSize = group_getTotalBaseLength(group);
			totalBitScore += (totalSequenceSize > 0 ? ((log(totalSequenceSize) / log(2.0)) + followingPathBitScore) * totalSequenceSize : 0.0);
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

void tabulateFloatStats(struct List *unsortedValues, double *totalNumber, double *min, double *max, double *avg, double *median) {
	if(unsortedValues->length == 0) {
		*totalNumber = 0;
		*min = INT32_MAX;
		*max = INT32_MAX;
		*avg = INT32_MAX;
		*median = INT32_MAX;
		return;
	}
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
}

void tabulateStats(struct IntList *unsortedValues, double *totalNumber, double *min, double *max, double *avg, double *median) {
	if(unsortedValues->length == 0) {
		*totalNumber = 0;
		*min = INT32_MAX;
		*max = INT32_MAX;
		*avg = INT32_MAX;
		*median = INT32_MAX;
		return;
	}
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
}

void largestChildStatsP(Net *net, struct List *childProportions) {
	Net_GroupIterator *groupIterator;
	Group *group;
	double problemSize = net_getTotalBaseLength(net);

	int32_t internalNode = 0;
	groupIterator = net_getGroupIterator(net);
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			internalNode = 1;
		}
	}
	net_destructGroupIterator(groupIterator);

	if(internalNode) {
		double childProportion = -10.0;
		double cumProp = 0.0;
		groupIterator = net_getGroupIterator(net);
		while((group = net_getNextGroup(groupIterator)) != NULL) {
			if(group_getNestedNet(group) != NULL) {
				largestChildStatsP(group_getNestedNet(group), childProportions);
			}
			double f = group_getTotalBaseLength(group);
			assert(f >= 0.0);
			assert(f/problemSize <= 1.001);
			if(f/problemSize > childProportion) {
				childProportion = f/problemSize;
			}
			cumProp += f/problemSize;
		}
		net_destructGroupIterator(groupIterator);
		assert(childProportion != -10.0);
		assert(cumProp <= 1.001);
		listAppend(childProportions, constructFloat(childProportion));
	}
}

void largestChildStats(Net *net,
		struct List **childProportions,
		double *minProportion, double *maxProportion, double *avgProportion, double *medianProportion) {
	*childProportions = constructEmptyList(0, (void (*)(void *))destructFloat);
	largestChildStatsP(net, *childProportions);
	double f;
	tabulateFloatStats(*childProportions, &f, minProportion, maxProportion, avgProportion, medianProportion);
}

void netStatsP(Net *net, int32_t currentDepth, struct IntList *children, struct IntList *depths) {
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			netStatsP(group_getNestedNet(group), currentDepth+1, children, depths);
		}
		else {
			intListAppend(depths, currentDepth);
		}
	}
	net_destructGroupIterator(groupIterator);
	intListAppend(children, net_getGroupNumber(net));
}

void netStats(Net *net,
		struct IntList **children, struct IntList **depths,
		double *totalNetNumber,
		double *maxChildren, double *avgChildren, double *medianChildren,
		double *minDepth, double *maxDepth, double *avgDepth, double *medianDepth) {
	*children = constructEmptyIntList(0);
	*depths = constructEmptyIntList(0);
	netStatsP(net, 1, *children, *depths);
	double f;
	tabulateStats(*children, totalNetNumber, &f, maxChildren, avgChildren, medianChildren);
	tabulateStats(*depths, totalNetNumber, minDepth, maxDepth, avgDepth, medianDepth);
	*totalNetNumber = (*children)->length + (*depths)->length;
}

void blockStatsP(Net *net, struct IntList *counts, struct IntList *lengths, struct IntList *degrees,
		struct IntList *coverage) {
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			blockStatsP(group_getNestedNet(group), counts, lengths, degrees, coverage);
		}
	}
	net_destructGroupIterator(groupIterator);

	Net_BlockIterator *blockIterator = net_getBlockIterator(net);
	Block *block;
	while((block = net_getNextBlock(blockIterator)) != NULL) {
		intListAppend(lengths, block_getLength(block));
		intListAppend(degrees, block_getInstanceNumber(block));
		intListAppend(coverage, block_getLength(block)*block_getInstanceNumber(block));
	}
	net_destructBlockIterator(blockIterator);
	intListAppend(counts, net_getBlockNumber(net));
}

void blockStats(Net *net,
		struct IntList **counts, struct IntList **lengths,
		struct IntList **degrees, struct IntList **coverage,
		double *totalBlockNumber,
		double *maxNumberPerNet, double *averageNumberPerNet, double *medianNumberPerNet,
		double *maxLength, double *averageLength, double *medianLength,
		double *maxDegree, double *averageDegree, double *medianDegree,
		double *maxCoverage, double *averageCoverage, double *medianCoverage) {
	*counts = constructEmptyIntList(0);
	*lengths = constructEmptyIntList(0);
	*degrees = constructEmptyIntList(0);
	*coverage = constructEmptyIntList(0);
	blockStatsP(net, *counts, *lengths, *degrees, *coverage);
	double f;
	tabulateStats(*counts, &f, &f, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet);
	tabulateStats(*lengths, totalBlockNumber, &f, maxLength, averageLength, medianLength);
	tabulateStats(*degrees, totalBlockNumber, &f, maxDegree, averageDegree, medianDegree);
	tabulateStats(*coverage, totalBlockNumber, &f, maxCoverage, averageCoverage, medianCoverage);
}

void chainStatsP(Net *net, struct IntList *counts, struct IntList *blockLengths,  struct IntList *lengths, struct IntList *baseLengths, struct IntList *avgInstanceBaseLength, int32_t minNumberOfBlocksInChain) {
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			chainStatsP(group_getNestedNet(group), counts, blockLengths, lengths, baseLengths, avgInstanceBaseLength, minNumberOfBlocksInChain);
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
			intListAppend(baseLengths, k);
			intListAppend(blockLengths, i);
			intListAppend(lengths, chain_getLength(chain));
			intListAppend(avgInstanceBaseLength, chain_getAverageInstanceBaseLength(chain));
			l++;
		}
	}
	net_destructBlockIterator(chainIterator);
	intListAppend(counts, l);
}

void chainStats(Net *net,
		struct IntList **counts, struct IntList **blockLengths, struct IntList **lengths, struct IntList **degrees, struct IntList **avgInstanceBaseLength,
		double *totalChainNumber,
		double *maxNumberPerNet, double *averageNumberPerNet, double *medianNumberPerNet,
		double *maxBlockLength, double *averageBlockLength, double *medianBlockLength,
		double *maxLength, double *averageLength, double *medianLength,
		double *maxDegree, double *averageDegree, double *medianDegree,
		double *maxInstanceLength, double *averageInstanceLength, double *medianInstanceLength, int32_t minNumberOfBlocksInChain) {
	*counts = constructEmptyIntList(0);
	*blockLengths = constructEmptyIntList(0);
	*lengths = constructEmptyIntList(0);
	*degrees = constructEmptyIntList(0);
	*avgInstanceBaseLength = constructEmptyIntList(0);
	chainStatsP(net, *counts, *blockLengths, *lengths, *degrees, *avgInstanceBaseLength, minNumberOfBlocksInChain);
	double f;
	tabulateStats(*counts, &f, &f, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet);
	tabulateStats(*blockLengths, totalChainNumber, &f, maxBlockLength, averageBlockLength, medianBlockLength);
	tabulateStats(*lengths, totalChainNumber, &f, maxLength, averageLength, medianLength);
	tabulateStats(*degrees, totalChainNumber, &f, maxDegree, averageDegree, medianDegree);
	tabulateStats(*avgInstanceBaseLength, totalChainNumber, &f, maxInstanceLength, averageInstanceLength, medianInstanceLength);
}

void endStatsP(Net *net, struct IntList *counts, struct IntList *degrees, int32_t includeNonTrivialGroups, int32_t includeTrivialGroups) {
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			endStatsP(group_getNestedNet(group), counts, degrees, includeNonTrivialGroups, includeTrivialGroups);
		}
		if((includeNonTrivialGroups && includeTrivialGroups) ||
				(includeNonTrivialGroups && group_getLink(group) == NULL) ||
				(includeTrivialGroups && group_getLink(group) != NULL)) {
			intListAppend(counts, group_getEndNumber(group));
			Group_EndIterator *endIterator = group_getEndIterator(group);
			End *end;
			while((end = group_getNextEnd(endIterator)) != NULL) {
				struct List *list = constructEmptyList(0, NULL);
				End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
				EndInstance *endInstance;
				while((endInstance = end_getNext(instanceIterator)) != NULL) {
					End *end =  end_getPositiveOrientation(endInstance_getEnd(endInstance_getAdjacency(endInstance)));
					if(!listContains(list, end)) {
						listAppend(list, end);
					}
				}
				end_destructInstanceIterator(instanceIterator);
				intListAppend(degrees, list->length);
				destructList(list);
			}
			group_destructEndIterator(endIterator);
		}

	}
	net_destructGroupIterator(groupIterator);
}

void endStats(Net *net,
		struct IntList **counts, struct IntList **degrees,
		double *maxNumberPerNet, double *averageNumberPerNet, double *medianNumberPerNet,
		double *maxDegree, double *averageDegree, double *medianDegree,
		int32_t includeNonTrivialGroups, int32_t includeTrivialGroups) {
	*counts = constructEmptyIntList(0);
	*degrees = constructEmptyIntList(0);
	endStatsP(net, *counts, *degrees, includeNonTrivialGroups, includeTrivialGroups);
	double f;
	tabulateStats(*counts, &f, &f, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet);
	tabulateStats(*degrees, &f, &f, maxDegree, averageDegree, medianDegree);
}

void leafStatsP(Net *net, struct IntList *leafSizes) {
	/*
	 * This only works while the reconstruction is incomplete.
	 */
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			leafStatsP(group_getNestedNet(group), leafSizes);
		}
		else {
			intListAppend(leafSizes, (int32_t)group_getTotalBaseLength(group));
		}
	}
	net_destructGroupIterator(groupIterator);
}

void leafStats(Net *net,
		struct IntList **leafSizes, double *totalLeafNumber,
		double *minSeqSize, double *maxSeqSize, double *avgSeqSize,
		double *medianSeqSize) {
	*leafSizes = constructEmptyIntList(0);
	leafStatsP(net, *leafSizes);
	tabulateStats(*leafSizes, totalLeafNumber, minSeqSize, maxSeqSize, avgSeqSize, medianSeqSize);
}

void nonTrivialGroupStatsP(Net *net, struct IntList *nonTrivialGroupCounts) {
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	int32_t i = 0;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_getNestedNet(group) != NULL) {
			nonTrivialGroupStatsP(group_getNestedNet(group), nonTrivialGroupCounts);
		}
		if(group_getLink(group) == NULL) {
			i++;
		}
	}
	net_destructGroupIterator(groupIterator);
	intListAppend(nonTrivialGroupCounts, i);
}

void nonTrivialGroupStats(Net *net,
		struct IntList **nonTrivialGroupCounts,
		double *totalNonTrivialGroups,
		double *maxNonTrivialGroups,
		double *avgNonTrivialGroups,
		double *medianNonTrivialGroups) {
	*nonTrivialGroupCounts = constructEmptyIntList(0);
	nonTrivialGroupStatsP(net, *nonTrivialGroupCounts);
	double f;
	tabulateStats(*nonTrivialGroupCounts,
			totalNonTrivialGroups, &f, maxNonTrivialGroups,
			avgNonTrivialGroups, medianNonTrivialGroups);
}

void printFloatValues(struct List *values, const char *tag, FILE *fileHandle) {
	int32_t i;
	fprintf(fileHandle, "<%s>", tag);
	for(i=0; i<values->length; i++) {
		fprintf(fileHandle, "%f ", *(float *)values->list[i]);
	}
	fprintf(fileHandle, "</%s>", tag);
}

void printIntValues(struct IntList *values, const char *tag, FILE *fileHandle) {
	int32_t i;
	fprintf(fileHandle, "<%s>", tag);
	for(i=0; i<values->length; i++) {
		fprintf(fileHandle, "%i ", values->list[i]);
	}
	fprintf(fileHandle, "</%s>", tag);
}

int main(int argc, char *argv[]) {
	/*
	 * The script builds a cactus tree representation of the chains and nets.
	 * The format of the output graph is dot format.
	 */
	NetDisk *netDisk;
	Net *net;
	FILE *fileHandle;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * netDiskName = NULL;
	char * netName = "0";
	char * outputFile = NULL;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "netName", required_argument, 0, 'd' },
			{ "outputFile", required_argument, 0, 'e' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:c:d:e:h", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stringCopy(optarg);
				break;
			case 'c':
				netDiskName = stringCopy(optarg);
				break;
			case 'd':
				netName = stringCopy(optarg);
				break;
			case 'e':
				outputFile = stringCopy(optarg);
				break;
			case 'h':
				usage();
				return 0;
			default:
				usage();
				return 1;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// (0) Check the inputs.
	///////////////////////////////////////////////////////////////////////////

	assert(netDiskName != NULL);
	assert(netName != NULL);
	assert(outputFile != NULL);

	//////////////////////////////////////////////
	//Set up logging
	//////////////////////////////////////////////

	if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
		setLogLevel(LOGGING_INFO);
	}
	if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	//////////////////////////////////////////////
	//Log (some of) the inputs
	//////////////////////////////////////////////

	logInfo("Net disk name : %s\n", netDiskName);
	logInfo("Net name : %s\n", netName);
	logInfo("Output graph file : %s\n", outputFile);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
	assert(net != NULL);
	logInfo("Parsed the top level net of the cactus tree to build\n");

	///////////////////////////////////////////////////////////////////////////
	// Calculate and print to file a crap load of numbers.
	///////////////////////////////////////////////////////////////////////////

	double totalSeqSize = net_getTotalBaseLength(net);

	fileHandle = fopen(outputFile, "w");
	fprintf(fileHandle, "<stats net_disk=\"%s\" net_name=\"%s\" total_sequence_length=\"%f\" >", netDiskName, netName, totalSeqSize);

	/*
	 * Relative entropy stats. Supposed to give a metric of how balanced the tree is in how it subdivides the input sequences.
	 */
	double totalP = calculateTreeBits(net, 0.0);
	double totalQ = (log(totalSeqSize) / log(2.0)) * totalSeqSize;
	//assert(totalP >= totalQ);
	double relativeEntropy = totalP - totalQ;
	double normalisedRelativeEntropy = relativeEntropy / totalSeqSize;

	fprintf(fileHandle, "<relative_entropy_stats totalP=\"%f\" totalQ=\"%f\" relative_entropy=\"%f\" normalised_relative_entropy=\"%f\"/>", totalP, totalQ, relativeEntropy, normalisedRelativeEntropy);

	/*
	 * Largest child stats - another tree balance stat. The idea is to supplant the more complex relative entropy function with a simple stat, that provides a cross check.
	 * The largest child is the proportion of the sequence the parent covers contained in a single child. min, max, avg and median numbers are given.
	 */
	double minProportion;
	double maxProportion;
	double avgProportion;
	double medianProportion;
	struct List *childProportions;
	largestChildStats(net, &childProportions, &minProportion, &maxProportion, &avgProportion, &medianProportion);
	fprintf(fileHandle, "<largest_child min_proportion=\"%f\" max_proportion=\"%f\" avg_proportion=\"%f\" median_proportion=\"%f\">", minProportion, maxProportion, avgProportion, medianProportion);
	printFloatValues(childProportions, "child_proportions", fileHandle);
	fprintf(fileHandle, "</largest_child>");
	destructList(childProportions);

	/*
	 * Numbers on the structure of the tree. Children numbers are based on the numbers of children each internal node in the cactus tree has.
	 * The depth numbers are how deep the tree is, in terms of total nodes on the path from the root to the leaves. So the min-depth is the minimum depth to a leaf node etc..
	 */
	double totalNetNumber;
	double maxChildren;
	double avgChildren;
	double medianChildren;
	double minDepth;
	double maxDepth;
	double avgDepth;
	double medianDepth;
	struct IntList *children;
	struct IntList *depths;
	netStats(net, &children, &depths, &totalNetNumber, &maxChildren, &avgChildren, &medianChildren, &minDepth,
			&maxDepth, &avgDepth, &medianDepth);
	fprintf(fileHandle, "<nets total_net_number=\"%f\" max_children=\"%f\" avg_children=\"%f\" median_children=\"%f\" min_depth=\"%f\" max_depth=\"%f\" avg_depth=\"%f\" median_depth=\"%f\">",
			totalNetNumber, maxChildren, avgChildren, medianChildren, minDepth, maxDepth, avgDepth, medianDepth);
	printIntValues(children, "children", fileHandle);
	printIntValues(depths, "depths", fileHandle);
	fprintf(fileHandle, "</nets>");
	destructIntList(children);
	destructIntList(depths);

	/*
	 * Numbers on the blocks. Length is the number base pairs in an block (they are gapless). Degree is the number of instances in an block. Coverage is the length multiplied by the degree.
	 */
	double totalNumber;
	double maxNumberPerNet;
	double averageNumberPerNet;
	double medianNumberPerNet;
	double maxLength;
	double averageLength;
	double medianLength;
	double maxDegree;
	double averageDegree;
	double medianDegree;
	double maxCoverage;
	double averageCoverage;
	double medianCoverage;

	struct IntList *counts;
	struct IntList *lengths;
	struct IntList *degrees;
	struct IntList *coverage;
	blockStats(net,
			&counts, &lengths, &degrees, &coverage,
			&totalNumber, &maxNumberPerNet, &averageNumberPerNet, &medianNumberPerNet,
			&maxLength, &averageLength, &medianLength,
			&maxDegree, &averageDegree, &medianDegree,
			&maxCoverage, &averageCoverage, &medianCoverage);
	fprintf(fileHandle, "<blocks total_number=\"%f\" total_base_length=\"%f\" max_number_per_net=\"%f\" average_number_per_net=\"%f\" median_number_per_net=\"%f\" max_length=\"%f\" average_length=\"%f\" median_length=\"%f\" max_degree=\"%f\" average_degree=\"%f\" median_degree=\"%f\" max_coverage=\"%f\" average_coverage=\"%f\" median_coverage=\"%f\">",
			totalNumber, totalNumber*averageLength, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet, maxLength, averageLength, medianLength, maxDegree, averageDegree, medianDegree, maxCoverage, averageCoverage, medianCoverage);
	printIntValues(counts, "counts", fileHandle);
	printIntValues(lengths, "lengths", fileHandle);
	printIntValues(degrees, "degrees", fileHandle);
	printIntValues(coverage, "coverage", fileHandle);
	fprintf(fileHandle, "</blocks>");
	destructIntList(counts);
	destructIntList(lengths);
	destructIntList(degrees);
	destructIntList(coverage);

	/*
	 * Chain statistics.
	 */
	double maxInstanceLength;
	double averageInstanceLength;
	double medianInstanceLength;
	double maxBlockLength;
	double averageBlockLength;
	double medianBlockLength;
	struct IntList *blockLengths;
	struct IntList *avgInstanceBaseLength;

	int32_t i;
	int32_t minimumBlockNumber[] = { 0, 2 };
	char *chainNames[] = { "chains", "chains_with_two_or_more_blocks" };

	for(i=0; i<2; i++) {
		chainStats(net,
					&counts, &blockLengths, &lengths, &degrees, &avgInstanceBaseLength,
					&totalNumber, &maxNumberPerNet, &averageNumberPerNet, &medianNumberPerNet,
					&maxBlockLength, &averageBlockLength, &medianBlockLength,
					&maxLength, &averageLength, &medianLength,
					&maxDegree, &averageDegree, &medianDegree,
					&maxInstanceLength, &averageInstanceLength, &medianInstanceLength, minimumBlockNumber[i]);
		fprintf(fileHandle, "<%s total_number=\"%f\" max_number_per_net=\"%f\" average_number_per_net=\"%f\" median_number_per_net=\"%f\" max_block_length=\"%f\" average__block_length=\"%f\" median_block_length=\"%f\" max_length=\"%f\" average_length=\"%f\" median_length=\"%f\" max_base_length=\"%f\" average_base_length=\"%f\" median_base_length=\"%f\" max_instance_length=\"%f\" average_instance_length=\"%f\" median_instance_length=\"%f\">",
				chainNames[i], totalNumber, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet, maxBlockLength, averageBlockLength, medianBlockLength, maxLength, averageLength, medianLength, maxDegree, averageDegree, medianDegree, maxInstanceLength, averageInstanceLength, medianInstanceLength);
		printIntValues(counts, "counts", fileHandle);
		printIntValues(blockLengths, "block_lengths", fileHandle);
		printIntValues(lengths, "lengths", fileHandle);
		printIntValues(degrees, "base_lengths", fileHandle);
		printIntValues(avgInstanceBaseLength, "avg_instance_lengths", fileHandle);
		fprintf(fileHandle, "</%s>", chainNames[i]);
		destructIntList(counts);
		destructIntList(lengths);
		destructIntList(degrees);
	}

	/*
	 * Stats on the ends in the problem. An end's degree is the number of distinct ends that it has adjacencies with.
	 * We break the numbers down into those for trivial/non-trivial groups and for both together.
	 */

	int32_t includeNonTrivialGroups[] = { 1, 1, 0 };
	int32_t includeTrivialGroups[] = { 1, 0, 1};
	char *endNames[] = { "ends", "ends_non_trivial", "ends_trivial" };

	for(i=0; i<3; i++) {
		endStats(net,
				&counts, &degrees,
				&maxNumberPerNet, &averageNumberPerNet, &medianNumberPerNet,
				&maxDegree, &averageDegree, &medianDegree, includeNonTrivialGroups[i], includeTrivialGroups[i]);
		fprintf(fileHandle, "<%s max_number_per_group=\"%f\" average_number_per_group=\"%f\" median_number_per_group=\"%f\" max_degree=\"%f\" average_degree=\"%f\" median_degree=\"%f\">", endNames[i], maxNumberPerNet, averageNumberPerNet, medianNumberPerNet, maxDegree, averageDegree, medianDegree);
		printIntValues(counts, "counts", fileHandle);
		printIntValues(degrees, "degrees", fileHandle);
		fprintf(fileHandle, "</%s>", endNames[i]);
		destructIntList(counts);
		destructIntList(degrees);
	}

	/*
	 * Stats on leaf nodes in the tree. The seq size of a leaf node is the amount of sequence (in bases), that it covers.
	 */
	double totalLeafNumber;
	double minSeqSize;
	double maxSeqSize;
	double avgSeqSize;
	double medianSeqSize;
	struct IntList *leafSizes;
	leafStats(net, &leafSizes, &totalLeafNumber, &minSeqSize, &maxSeqSize, &avgSeqSize, &medianSeqSize);
	fprintf(fileHandle, "<leaves total_leaf_number=\"%f\" min_seq_size=\"%f\" max_seq_size=\"%f\" avg_seq_size=\"%f\" median_seq_size=\"%f\" >", totalLeafNumber, minSeqSize, maxSeqSize, avgSeqSize, medianSeqSize);
	printIntValues(leafSizes, "leaf_sizes", fileHandle);
	fprintf(fileHandle, "</leaves>");
	destructIntList(leafSizes);

	/*
	 * Stats on the non trivial group.
	 */
	double totalNonTrivialGroups;
	double maxNonTrivialGroupsPerNet;
	double avgNonTrivialGroupsPerNet;
	double medianNonTrivialGroupsPetNet;
	struct IntList *nonTrivialGroupCounts;
	nonTrivialGroupStats(net,
			&nonTrivialGroupCounts,
			&totalNonTrivialGroups,
			&maxNonTrivialGroupsPerNet,
			&avgNonTrivialGroupsPerNet,
			&medianNonTrivialGroupsPetNet);
	fprintf(fileHandle, "<non_trivial_groups total=\"%f\" max_per_net=\"%f\" avg_per_net=\"%f\" median_per_net=\"%f\" >", totalNonTrivialGroups,
			maxNonTrivialGroupsPerNet, avgNonTrivialGroupsPerNet, medianNonTrivialGroupsPetNet);
	printIntValues(nonTrivialGroupCounts, "trivial_groups", fileHandle);
	fprintf(fileHandle, "</non_trivial_groups>");
	destructIntList(nonTrivialGroupCounts);


	fprintf(fileHandle, "</stats>\n");
	fclose(fileHandle);
	logInfo("Finished writing out the stats.\n");

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	netDisk_destruct(netDisk);

	return 0;
}
