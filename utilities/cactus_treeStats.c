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

double calculateTotalContainedSequence(Net *net) {
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	double totalLength = 0.0;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(!end_isAtomEnd(end)) {
			End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
			EndInstance *endInstance;
			while((endInstance = end_getNext(instanceIterator)) != NULL) {
				endInstance = endInstance_getStrand(endInstance) ? endInstance : endInstance_getReverse(endInstance);
				if(!endInstance_getSide(endInstance)) {
					EndInstance *endInstance2 = endInstance_getAdjacency(endInstance);
					while(end_isAtomEnd(endInstance_getEnd(endInstance2))) {
						AtomInstance *atomInstance = endInstance_getAtomInstance(endInstance2);
						assert(atomInstance != NULL);
						assert(atomInstance_get5End(atomInstance) == endInstance2);
						endInstance2 = endInstance_getAdjacency(atomInstance_get3End(atomInstance));
						assert(endInstance_getStrand(endInstance2));
						assert(endInstance_getSide(endInstance2));
					}
					assert(endInstance_getStrand(endInstance2));
					assert(endInstance_getSide(endInstance2));
					int32_t length = endInstance_getCoordinate(endInstance2) - endInstance_getCoordinate(endInstance) - 1;
					assert(length >= 0);
					totalLength += length;
				}
			}
			end_destructInstanceIterator(instanceIterator);
		}
	}
	net_destructEndIterator(endIterator);
	return totalLength;
}

double calculateTreeBits(Net *net, double pathBitScore) {
	if(net_getAdjacencyComponentNumber(net) > 0) { //internal node
		double totalBitScore = 0.0;
		Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
		AdjacencyComponent *adjacencyComponent;
		double followingPathBitScore = (log(net_getAdjacencyComponentNumber(net)) / log(2.0)) + pathBitScore;
		while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
			totalBitScore += calculateTreeBits(adjacencyComponent_getNestedNet(adjacencyComponent), followingPathBitScore);
		}
		net_destructAdjacencyComponentIterator(adjacencyComponentIterator);
		Net_AtomIterator *atomIterator = net_getAtomIterator(net);
		Atom *atom;
		int32_t totalSequenceSize = 0.0;
		while((atom = net_getNextAtom(atomIterator)) != NULL) {
			totalSequenceSize += atom_getLength(atom) * atom_getInstanceNumber(atom);
		}
		net_destructAtomIterator(atomIterator);
		return totalBitScore + (totalSequenceSize > 0 ? ((log(totalSequenceSize) / log(2.0)) + pathBitScore) * totalSequenceSize : 0.0);
	}
	assert(net_getAtomNumber(net) == 0);
	double i = calculateTotalContainedSequence(net);
	return i > 0 ? (pathBitScore + log(i)/log(2.0)) * i : 0.0;
}

void tabulateFloatStats(struct List *unsortedValues, double *totalNumber, double *min, double *max, double *avg, double *median) {
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

double largestChildStatsP(Net *net, struct List *childProportions) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	double problemSize = calculateTotalContainedSequence(net);
	if(net_getAdjacencyComponentNumber(net) != 0  && problemSize > 0) {
		double childProportion = -10.0;
		while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
			double f = largestChildStatsP(adjacencyComponent_getNestedNet(adjacencyComponent), childProportions);
			if(f/problemSize > childProportion) {
				childProportion = f/problemSize;
			}
		}
		net_destructAdjacencyComponentIterator(adjacencyComponentIterator);
		listAppend(childProportions, constructFloat(childProportion));
	}
	return problemSize;
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
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		netStatsP(adjacencyComponent_getNestedNet(adjacencyComponent), currentDepth+1, children, depths);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	if(net_getAdjacencyComponentNumber(net) == 0) {
		intListAppend(depths, currentDepth);
	}
	else {
		intListAppend(children, net_getAdjacencyComponentNumber(net));
	}
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

void atomStatsP(Net *net, struct IntList *counts, struct IntList *lengths, struct IntList *degrees,
		struct IntList *coverage) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		atomStatsP(adjacencyComponent_getNestedNet(adjacencyComponent), counts, lengths, degrees, coverage);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	Net_AtomIterator *atomIterator = net_getAtomIterator(net);
	Atom *atom;
	while((atom = net_getNextAtom(atomIterator)) != NULL) {
		intListAppend(lengths, atom_getLength(atom));
		intListAppend(degrees, atom_getInstanceNumber(atom));
		intListAppend(coverage, atom_getLength(atom)*atom_getInstanceNumber(atom));
	}
	net_destructAtomIterator(atomIterator);
	intListAppend(counts, net_getAtomNumber(net));
}

void atomStats(Net *net,
		struct IntList **counts, struct IntList **lengths,
		struct IntList **degrees, struct IntList **coverage,
		double *totalAtomNumber,
		double *maxNumberPerNet, double *averageNumberPerNet, double *medianNumberPerNet,
		double *maxLength, double *averageLength, double *medianLength,
		double *maxDegree, double *averageDegree, double *medianDegree,
		double *maxCoverage, double *averageCoverage, double *medianCoverage) {
	*counts = constructEmptyIntList(0);
	*lengths = constructEmptyIntList(0);
	*degrees = constructEmptyIntList(0);
	*coverage = constructEmptyIntList(0);
	atomStatsP(net, *counts, *lengths, *degrees, *coverage);
	double f;
	tabulateStats(*counts, &f, &f, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet);
	tabulateStats(*lengths, totalAtomNumber, &f, maxLength, averageLength, medianLength);
	tabulateStats(*degrees, totalAtomNumber, &f, maxDegree, averageDegree, medianDegree);
	tabulateStats(*coverage, totalAtomNumber, &f, maxCoverage, averageCoverage, medianCoverage);
}

void chainStatsP(Net *net, struct IntList *counts, struct IntList *lengths, struct IntList *baseLengths) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		chainStatsP(adjacencyComponent_getNestedNet(adjacencyComponent), counts, lengths, baseLengths);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	Net_ChainIterator *chainIterator = net_getChainIterator(net);
	Chain *chain;
	Atom **atoms;
	int32_t i, j, k;
	while((chain = net_getNextChain(chainIterator)) != NULL) {
		atoms = chain_getAtomChain(chain, &i);
		k = 0;
		for(j=0; j<i; j++) {
			k += atom_getLength(atoms[j]);
		}
		intListAppend(baseLengths, k);
		intListAppend(lengths, chain_getLength(chain));
	}
	net_destructAtomIterator(chainIterator);
	intListAppend(counts, net_getChainNumber(net));
}

void chainStats(Net *net,
		struct IntList **counts, struct IntList **lengths, struct IntList **degrees,
		double *totalChainNumber,
		double *maxNumberPerNet, double *averageNumberPerNet, double *medianNumberPerNet,
		double *maxLength, double *averageLength, double *medianLength,
		double *maxDegree, double *averageDegree, double *medianDegree) {
	*counts = constructEmptyIntList(0);
	*lengths = constructEmptyIntList(0);
	*degrees = constructEmptyIntList(0);
	chainStatsP(net, *counts, *lengths, *degrees);
	double f;
	tabulateStats(*counts, &f, &f, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet);
	tabulateStats(*lengths, totalChainNumber, &f, maxLength, averageLength, medianLength);
	tabulateStats(*degrees, totalChainNumber, &f, maxDegree, averageDegree, medianDegree);
}

void endStatsP(Net *net, struct IntList *counts, struct IntList *degrees) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		endStatsP(adjacencyComponent_getNestedNet(adjacencyComponent), counts, degrees);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	intListAppend(counts, net_getEndNumber(net));
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	while((end = net_getNextEnd(endIterator)) != NULL) {
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
	net_destructEndIterator(endIterator);
}

void endStats(Net *net,
		struct IntList **counts, struct IntList **degrees,
		double *maxNumberPerNet, double *averageNumberPerNet, double *medianNumberPerNet,
		double *maxDegree, double *averageDegree, double *medianDegree) {
	*counts = constructEmptyIntList(0);
	*degrees = constructEmptyIntList(0);
	endStatsP(net, *counts, *degrees);
	double f;
	tabulateStats(*counts, &f, &f, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet);
	tabulateStats(*degrees, &f, &f, maxDegree, averageDegree, medianDegree);
}

void leafStatsP(Net *net, struct IntList *leafSizes) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		leafStatsP(adjacencyComponent_getNestedNet(adjacencyComponent), leafSizes);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	if(net_getAdjacencyComponentNumber(net) == 0) {
		intListAppend(leafSizes, (int32_t)calculateTotalContainedSequence(net));
	}
}

void leafStats(Net *net,
		struct IntList **leafSizes, double *totalLeafNumber,
		double *minSeqSize, double *maxSeqSize, double *avgSeqSize,
		double *medianSeqSize) {
	*leafSizes = constructEmptyIntList(0);
	leafStatsP(net, *leafSizes);
	tabulateStats(*leafSizes, totalLeafNumber, minSeqSize, maxSeqSize, avgSeqSize, medianSeqSize);
}

void nonTrivialAdjacencyComponentStatsP(Net *net, struct IntList *nonTrivialAdjacencyComponentCounts) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	int32_t i = 0;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		nonTrivialAdjacencyComponentStatsP(adjacencyComponent_getNestedNet(adjacencyComponent), nonTrivialAdjacencyComponentCounts);
		if(adjacencyComponent_getLink(adjacencyComponent) == NULL) {
			i++;
		}
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	if(net_getAdjacencyComponentNumber(net) != 0) {
		intListAppend(nonTrivialAdjacencyComponentCounts, i);
	}
}

void nonTrivialAdjacencyComponentStats(Net *net,
		struct IntList **nonTrivialAdjacencyComponentCounts,
		double *totalNonTrivialAdjacencyComponents,
		double *maxNonTrivialAdjacencyComponents,
		double *avgNonTrivialAdjacencyComponents,
		double *medianNonTrivialAdjacencyComponents) {
	*nonTrivialAdjacencyComponentCounts = constructEmptyIntList(0);
	nonTrivialAdjacencyComponentStatsP(net, *nonTrivialAdjacencyComponentCounts);
	double f;
	tabulateStats(*nonTrivialAdjacencyComponentCounts,
			totalNonTrivialAdjacencyComponents, &f, maxNonTrivialAdjacencyComponents,
			avgNonTrivialAdjacencyComponents, medianNonTrivialAdjacencyComponents);
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
	char * netName = NULL;
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

	assert(logLevelString == NULL || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
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

	double totalSeqSize = calculateTotalContainedSequence(net);

	fileHandle = fopen(outputFile, "w");
	fprintf(fileHandle, "<stats netDisk=\"%s\" netName=\"%s\" totalSequenceLength=\"%f\" >", netDiskName, netName, totalSeqSize);

	/*
	 * Relative entropy stats. Supposed to give a metric of how balanced the tree is in how it subdivides the input sequences.
	 */
	double totalP = calculateTreeBits(net, 0.0);
	double totalQ = (log(totalSeqSize) / log(2.0)) * totalSeqSize;
	assert(totalP >= totalQ);
	double relativeEntropy = totalP - totalQ;
	double normalisedRelativeEntropy = relativeEntropy / totalSeqSize;

	fprintf(fileHandle, "<relativeEntropyStats totalP=\"%f\" totalQ=\"%f\" relativeEntropy=\"%f\" normalisedRelativeEntropy=\"%f\"/>", totalP, totalQ, relativeEntropy, normalisedRelativeEntropy);

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
	fprintf(fileHandle, "<largestChild minProportion=\"%f\" maxProportion=\"%f\" avgProportion=\"%f\" medianProportion=\"%f\">", minProportion, maxProportion, avgProportion, medianProportion);
	printFloatValues(childProportions, "childProportions", fileHandle);
	fprintf(fileHandle, "</largestChild>");
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
	fprintf(fileHandle, "<nets totalNetNumber=\"%f\" maxChildren=\"%f\" avgChildren=\"%f\" medianChildren=\"%f\" minDepth=\"%f\" maxDepth=\"%f\" avgDepth=\"%f\" medianDepth=\"%f\">",
			totalNetNumber, maxChildren, avgChildren, medianChildren, minDepth, maxDepth, avgDepth, medianDepth);
	printIntValues(children, "children", fileHandle);
	printIntValues(depths, "depths", fileHandle);
	fprintf(fileHandle, "</nets>");
	destructIntList(children);
	destructIntList(depths);

	/*
	 * Numbers on the atoms. Length is the number base pairs in an atom (they are gapless). Degree is the number of instances in an atom. Coverage is the length multiplied by the degree.
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
	atomStats(net,
			&counts, &lengths, &degrees, &coverage,
			&totalNumber, &maxNumberPerNet, &averageNumberPerNet, &medianNumberPerNet,
			&maxLength, &averageLength, &medianLength,
			&maxDegree, &averageDegree, &medianDegree,
			&maxCoverage, &averageCoverage, &medianCoverage);
	fprintf(fileHandle, "<atoms totalNumber=\"%f\" totalBaseLength=\"%f\" maxNumberPerNet=\"%f\" averageNumberPerNet=\"%f\" medianNumberPerNet=\"%f\" maxLength=\"%f\" averageLength=\"%f\" medianLength=\"%f\" maxDegree=\"%f\" averageDegree=\"%f\" medianDegree=\"%f\" maxCoverage=\"%f\" averageCoverage=\"%f\" medianCoverage=\"%f\">",
			totalNumber, totalNumber*averageLength, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet, maxLength, averageLength, medianLength, maxDegree, averageDegree, medianDegree, maxCoverage, averageCoverage, medianCoverage);
	printIntValues(counts, "counts", fileHandle);
	printIntValues(lengths, "lengths", fileHandle);
	printIntValues(degrees, "degrees", fileHandle);
	printIntValues(coverage, "coverage", fileHandle);
	fprintf(fileHandle, "</atoms>");
	destructIntList(counts);
	destructIntList(lengths);
	destructIntList(degrees);
	destructIntList(coverage);

	/*
	 * Chain statistics.
	 */
	chainStats(net,
				&counts, &lengths, &degrees,
				&totalNumber, &maxNumberPerNet, &averageNumberPerNet, &medianNumberPerNet,
				&maxLength, &averageLength, &medianLength,
				&maxDegree, &averageDegree, &medianDegree);
	fprintf(fileHandle, "<chains totalNumber=\"%f\" maxNumberPerNet=\"%f\" averageNumberPerNet=\"%f\" medianNumberPerNet=\"%f\" maxLength=\"%f\" averageLength=\"%f\" medianLength=\"%f\" maxBaseLength=\"%f\" averageBaseLength=\"%f\" medianBaseLength=\"%f\">",
			totalNumber, maxNumberPerNet, averageNumberPerNet, medianNumberPerNet, maxLength, averageLength, medianLength, maxDegree, averageDegree, medianDegree);
	printIntValues(counts, "counts", fileHandle);
	printIntValues(lengths, "lengths", fileHandle);
	printIntValues(degrees, "degrees", fileHandle);
	fprintf(fileHandle, "</chains>");
	destructIntList(counts);
	destructIntList(lengths);
	destructIntList(degrees);

	/*
	 * Stats on the ends in the problem. An end's degree is the number of distinct ends that it has adjacencies with.
	 */
	endStats(net,
			&counts, &degrees,
			&maxNumberPerNet, &averageNumberPerNet, &medianNumberPerNet,
			&maxDegree, &averageDegree, &medianDegree);
	fprintf(fileHandle, "<ends maxNumberPerNet=\"%f\" averageNumberPerNet=\"%f\" medianNumberPerNet=\"%f\" maxDegree=\"%f\" averageDegree=\"%f\" medianDegree=\"%f\">", maxNumberPerNet, averageNumberPerNet, medianNumberPerNet, maxDegree, averageDegree, medianDegree);
	printIntValues(counts, "counts", fileHandle);
	printIntValues(degrees, "degrees", fileHandle);
	fprintf(fileHandle, "</ends>");
	destructIntList(counts);
	destructIntList(degrees);

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
	fprintf(fileHandle, "<leaves totalLeafNumber=\"%f\" minSeqSize=\"%f\" maxSeqSize=\"%f\" avgSeqSize=\"%f\" medianSeqSize=\"%f\" >", totalLeafNumber, minSeqSize, maxSeqSize, avgSeqSize, medianSeqSize);
	printIntValues(leafSizes, "leafSizes", fileHandle);
	fprintf(fileHandle, "</leaves>");
	destructIntList(leafSizes);

	/*
	 * Stats on the non trivial adjacency component.
	 */
	double totalNonTrivialAdjacencyComponents;
	double maxNonTrivialAdjacencyComponentsPerNet;
	double avgNonTrivialAdjacencyComponentsPerNet;
	double medianNonTrivialAdjacencyComponentsPetNet;
	struct IntList *nonTrivialAdjacencyComponentCounts;
	nonTrivialAdjacencyComponentStats(net,
			&nonTrivialAdjacencyComponentCounts,
			&totalNonTrivialAdjacencyComponents,
			&maxNonTrivialAdjacencyComponentsPerNet,
			&avgNonTrivialAdjacencyComponentsPerNet,
			&medianNonTrivialAdjacencyComponentsPetNet);
	fprintf(fileHandle, "<nonTrivialGroups total=\"%f\" maxPerNet=\"%f\" avgPerNet=\"%f\" medianPerNet=\"%f\" >", totalNonTrivialAdjacencyComponents,
			maxNonTrivialAdjacencyComponentsPerNet, avgNonTrivialAdjacencyComponentsPerNet, medianNonTrivialAdjacencyComponentsPetNet);
	printIntValues(nonTrivialAdjacencyComponentCounts, "trivialAdjacencyComponents", fileHandle);
	fprintf(fileHandle, "</nonTrivialGroups>");
	destructIntList(nonTrivialAdjacencyComponentCounts);


	fprintf(fileHandle, "</stats>\n");
	fclose(fileHandle);
	logInfo("Finished writing out the stats.\n");

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	netDisk_destruct(netDisk);

	return 0;
}
