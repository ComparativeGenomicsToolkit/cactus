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

void tabulateStats(struct IntList *unsortedValues, double *totalNumber, double *min, double *max, double *avg, double *median) {
	assert(unsortedValues->length > 0);
	qsort(unsortedValues->list, unsortedValues->length, sizeof(int32_t), (int (*)(const void *, const void *))intComparator_Int);
	*totalNumber = unsortedValues->length;
	*min = unsortedValues->list[0];
	*max = unsortedValues->list[unsortedValues->length-1];
	*median = unsortedValues->list[unsortedValues->length/2];
	int32_t i, j = 0;
	for(i=0; i<unsortedValues->length; i++) {
		j += unsortedValues->list[j];
	}
	*avg = j / unsortedValues->length;
}

void netStatsP(Net *net, int32_t currentDepth, struct IntList *depths) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		netStatsP(adjacencyComponent_getNestedNet(adjacencyComponent), currentDepth+1, depths);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	if(net_getAdjacencyComponentNumber(net) == 0) {
		intListAppend(depths, currentDepth);
	}
}

void netStats(Net *net, double *totalNetNumber, double *minDepth, double *maxDepth, double *avgDepth, double *medianDepth) {
	struct IntList *depths = constructEmptyIntList(0);
	netStatsP(net, 1, depths);
	tabulateStats(depths, totalNetNumber, minDepth, maxDepth, avgDepth, medianDepth);
	destructIntList(depths);
}

void atomStatsP(Net *net, int32_t currentDepth, struct IntList *depths) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		netStatsP(adjacencyComponent_getNestedNet(adjacencyComponent), currentDepth+1, depths);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	if(net_getAdjacencyComponentNumber(net) == 0) {
		intListAppend(depths, currentDepth);
	}
}

void netStats(Net *net, double *totalNetNumber, double *minDepth, double *maxDepth, double *avgDepth, double *medianDepth) {
	struct IntList *depths = constructEmptyIntList(0);
	netStatsP(net, 1, depths);
	tabulateStats(depths, totalNetNumber, minDepth, maxDepth, avgDepth, medianDepth);
	destructIntList(depths);
}

//atom number,  average atom's per node, atom avg length, average atom degree
//chain number, average chain's per node, max length, average length, median length, alignment length, max, min, average, median gap size
//max, average, median end degree


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
	// Calculate the relative entropy.
	///////////////////////////////////////////////////////////////////////////

	double totalP = calculateTreeBits(net, 0.0);
	double i = calculateTotalContainedSequence(net);
	double totalQ = (log(i) / log(2.0)) * i;
	assert(totalP >= totalQ);
	double relativeEntropy = totalP - totalQ;
	double normalisedRelativeEntropy = relativeEntropy / i;

	logInfo("The total sequence length: %f, total P: %f, Q: %f, relative entropy of the cactus tree: %f, sequence length normalised relative entropy: %f\n", i, totalP, totalQ, relativeEntropy, normalisedRelativeEntropy);
	fileHandle = fopen(outputFile, "w");
	fprintf(fileHandle, "<stats netDisk=\"%s\" netName=\"%s\" totalSequenceLength=\"%f\" totalP=\"%f\" totalQ=\"%f\" relativeEntropy=\"%f\" normalisedRelativeEntropy=\"%f\"/>\n", netDiskName, netName, i, totalP, totalQ, relativeEntropy, normalisedRelativeEntropy);
	fclose(fileHandle);

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	netDisk_destruct(netDisk);

	return 0;
}
