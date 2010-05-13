#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <getopt.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>

#include "pairwiseAlignmentModel.h"
#include "pairwiseAlignmentModelInterface.h"
#include "banding.h"
#include "algebras.h"

#include "XMLTools.h"
#include "xmlParser.h"
#include "substitutionIO.h"

#define SPANNING_TREES 2
#define MATCH_THRESHOLD 0.7
#define MAXIMUM_LENGTH 1000

extern "C" {
	#include "hashTableC.h"
	#include "bioioC.h"
	#include "commonC.h"
	#include "cactus.h"
	#include "avl.h"
	#include "pairwiseAlignment.h"
	#include "cactus_core.h"
};

void usage() {
	fprintf(stderr, "cactus_colinearAligner [net-names, ordered by order they should be processed], version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

char *getSequence(Cap *cap, int32_t maxLength) {
	Sequence *sequence = cap_getSequence(cap);
	assert(sequence != NULL);
	Cap *cap2 = cap_getAdjacency(cap);
	assert(!cap_getSide(cap));

	if(cap_getStrand(cap)) {
		int32_t length = cap_getCoordinate(cap2)
					- cap_getCoordinate(cap) - 1;
		assert(length >= 0);
		assert(maxLength >= 0);
		return sequence_getString(sequence, cap_getCoordinate(
					cap) + 1, length > maxLength ? maxLength : length, 1);
	}
	else {
		int32_t length = cap_getCoordinate(cap)
				- cap_getCoordinate(cap2) - 1;
		assert(length >= 0);
		return sequence_getString(sequence,
				length > maxLength ?
				cap_getCoordinate(cap) - maxLength:
				cap_getCoordinate(cap2) + 1,
				length > maxLength ? maxLength : length, 0);
	}
}

typedef struct _SubSequence {
	char *string;
	char *sequenceName;
	int32_t strand;
	int32_t start;
	int32_t length;
} SubSequence;

SubSequence *getSequenceAndCoordinates(Cap *cap, int32_t maxLength) {
	assert(!cap_getSide(cap));
	SubSequence *subSequence = (SubSequence *)mallocLocal(sizeof(SubSequence));
	subSequence->string = getSequence(cap, maxLength);
	Sequence *sequence = cap_getSequence(cap);
	subSequence->sequenceName = netMisc_nameToString(sequence_getName(sequence));
	subSequence->strand = cap_getStrand(cap);
	subSequence->start = cap_getCoordinate(cap) + (cap_getStrand(cap) ? 1 : -1);
	subSequence->length = strlen(subSequence->string);
	return subSequence;
}

void destructSubSequence(SubSequence *subSequence) {
	free(subSequence->string);
	free(subSequence->sequenceName);
	free(subSequence);
}

SubSequence **getForwardAndReverseSequences(Cap *cap1, int32_t maxLength) {
	assert(!cap_getSide(cap1));

	Cap *cap2 = cap_getAdjacency(cap1);

	assert(cap2 != NULL);
	assert(cap_getStrand(cap2));
	assert(cap_getSide(cap2));

	SubSequence **subSequences = (SubSequence **)mallocLocal(sizeof(void *) * 2);
	subSequences[0] = getSequenceAndCoordinates(cap1, maxLength);
	subSequences[1] = getSequenceAndCoordinates(cap_getReverse(cap2), maxLength);
	return subSequences;
}

void destructSubsequences(SubSequence **subSequences) {
	destructSubSequence(subSequences[0]);
	destructSubSequence(subSequences[1]);
	free(subSequences);
}

struct PairwiseAlignment *alignSequences(SubSequence *seqX, SubSequence *seqY) {
	vector<int> constraints;
	AlignerDPTable *fTable;
	AlignerDPTable *bTable;
	bfloat fProb;
	vector<double> gapPosteriorsX;
	vector<double> gapPosteriorsY;
	double divergence;
	struct PairwiseAlignmentInputParameters *pAIP = constructPairwiseAlignmentInputParameters();
	pAIP->matchThreshold = MATCH_THRESHOLD;
	XMLNode xMainNode = getDefaultModelParams();
	struct ModelParameters *modelParameters = readModelParameters(xMainNode);

	//Convert the sequences.
	char *convertedSeqX = convertSequence(seqX->string, convertFromDNA);
	char *convertedSeqY = convertSequence(seqY->string, convertFromDNA);

	//Do the alignment.
	vector<int> diagonals =
			computePairwiseAlignment(convertedSeqX, convertedSeqY,
					&modelParameters, pAIP, constraints,
					&fTable, &bTable, &fProb, gapPosteriorsX,
					gapPosteriorsY, &divergence);

	struct PairwiseAlignment *pairwiseAlignment =
			constructAlignment(diagonals,
						seqX->sequenceName,
						seqY->sequenceName,
						seqX->length, seqY->length,
						seqX->start, seqX->strand,
						seqY->start, seqY->strand);

	//Clean up.
	free(convertedSeqX);
	free(convertedSeqY);
	delete fTable;
	delete bTable;
	destructModelParameters(modelParameters);
	destructPairwiseAlignmentInputParameters(pAIP);

	return pairwiseAlignment;
}

void constructSpanningTree(void **items, uint32_t length, struct List *pairs) {
	/*
	 * Constructs a spanning tree linking all the nodes in items into one component.
	 */
	if(length <= 1) {
		return;
	}
	uint32_t i = length/2;
	if(length > 2) {
		constructSpanningTree(items, i, pairs);
		constructSpanningTree(items + i, length-i, pairs);
	}
	//Only add if it's a new pair.
	uint32_t j = min((uint32_t)(RANDOM() * i), i-1);
	assert(j < i);
	void *item1 = items[j];
	j = min((uint32_t)(i + (RANDOM() * (length - i))), length-1);
	assert(j >= i);
	assert(j < length);
	void *item2 = items[j];
	int32_t k=0;
	for(k=0; k<pairs->length; k+=2) {
		if((pairs->list[k] == item1 && pairs->list[k+1] == item2) ||
		   (pairs->list[k] == item2 && pairs->list[k+1] == item1)) {
			return;
		}
	}
	listAppend(pairs, item1);
	listAppend(pairs, item2);

}

static struct List *getAlignment_pairwiseAlignments;
static int32_t getAlignment_index = 0;

static struct PairwiseAlignment *getAlignments() {
	while(getAlignment_index < getAlignment_pairwiseAlignments->length) {
		return (struct PairwiseAlignment *)getAlignment_pairwiseAlignments->list[getAlignment_index++];
	}
	return NULL;
}

static void startAlignmentStack() {
	getAlignment_index = 0;
}

int main(int argc, char *argv[]) {

	char * logLevelString = NULL;
	char * netDiskName = NULL;
	int32_t i, j;

	/*
	 * Parse the options.
	 */
	while(1) {
		static struct option long_options[] = {
				{ "logLevel", required_argument, 0, 'a' },
				{ "netDisk", required_argument, 0, 'b' },
				{ "help", no_argument, 0, 'h' },
				{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:b:h", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stringCopy(optarg);
				break;
			case 'b':
				netDiskName = stringCopy(optarg);
				break;
			case 'h':
				usage();
				return 0;
			default:
				usage();
				return 1;
		}
	}

	if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
		setLogLevel(LOGGING_INFO);
	}
	if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	/*
	 * Setup the input parameters for cactus core.
	 */
	CactusCoreInputParameters *cCIP = constructCactusCoreInputParameters();
	//--maxEdgeDegree 10000000 --minimumTreeCoverage 0 --minimumTreeCoverageForBlocks 0
	//--minimumBlockLength 0 --minimumChainLength 0 --trim 0 --alignRepeats 1 --extensionSteps 0

	/*
	 * Load the netdisk
	 */
	NetDisk *netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	/*
	 * For each net.
	 */
	for (j = optind; j < argc; j++) {
		/*
		 * Read the net.
		 */
		const char *netName = argv[j];
		logInfo("Processing the net named: %s\n", netName);
		Net *net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
		assert(net != NULL);
		logInfo("Parsed the net to be aligned\n");

		/*
		 * Get the sequences
		 */
		Net_CapIterator *iterator1 = net_getCapIterator(net);
		Cap *cap1;
		struct List *subSequences = constructEmptyList(0, (void (*)(void *))destructSubsequences);
		while((cap1 = net_getNextCap(iterator1)) != NULL) {
			cap1 = cap_getStrand(cap1) ? cap1 : cap_getReverse(cap1);
			if(!cap_getSide(cap1)) {
				SubSequence **subSequences1 = getForwardAndReverseSequences(cap1, MAXIMUM_LENGTH);
				listAppend(subSequences, subSequences1);
				SubSequence *seq1 = subSequences1[0];
				SubSequence *seq1R = subSequences1[1];

				assert(seq1->length <= MAXIMUM_LENGTH);
				assert(seq1R->length <= MAXIMUM_LENGTH);
				assert((seq1->length == MAXIMUM_LENGTH && seq1R->length == MAXIMUM_LENGTH) ||
						(seq1->length != MAXIMUM_LENGTH && seq1R->length != MAXIMUM_LENGTH));
			}
		}
		net_destructCapIterator(iterator1);
		assert(subSequences->length == net_getCapNumber(net)/2);
		logInfo("Got the sequence to be aligned\n");

		/*
		 * Make a list of the pairwise comparisons to be aligned.
		 */
		struct List *pairs = constructEmptyList(0, NULL);
		for(i=0; i<SPANNING_TREES; i++) {
			arrayShuffle(subSequences->list, subSequences->length);
			constructSpanningTree(subSequences->list, subSequences->length, pairs);
		}
		//check we have the expected number of pairs.
		assert(subSequences->length == 0 ||
				(pairs->length/2 <= SPANNING_TREES*(subSequences->length-1) && (pairs->length/2 >= subSequences->length-1)));
		logInfo("Got the list of pairs to align\n");

		/*
		 * Do the alignments.
		 */
		getAlignment_pairwiseAlignments = constructEmptyList(0, (void (*)(void *))destructPairwiseAlignment);
		for(i=0; i<pairs->length; i+=2) {
			SubSequence **subSequences1 = (SubSequence **)pairs->list[i];
			SubSequence **subSequences2 = (SubSequence **)pairs->list[i+1];

			SubSequence *seq1 = subSequences1[0];
			SubSequence *seq1R = subSequences1[1];
			SubSequence *seq2 = subSequences2[0];
			SubSequence *seq2R = subSequences2[1];

			//Do forward - forward alignment.
			listAppend(getAlignment_pairwiseAlignments, alignSequences(seq1, seq2));

			//Do reverse - forward alignment.
			listAppend(getAlignment_pairwiseAlignments, alignSequences(seq1R, seq2));

			//Check if sequences are greater maximum length and do extra alignments..
			if(seq2->length == MAXIMUM_LENGTH || seq1->length == MAXIMUM_LENGTH) {
				//Do forward - reverse alignment.
				listAppend(getAlignment_pairwiseAlignments, alignSequences(seq1, seq2R));

				//Do reverse - reverse alignment.
				listAppend(getAlignment_pairwiseAlignments, alignSequences(seq1R, seq2R));
			}
		}
		logInfo("Created the alignments\n");

		/*
		 * Run the cactus core script.
		 */
		cCIP->maxEdgeDegree = subSequences->length;
		cactusCorePipeline(net, cCIP, getAlignments, startAlignmentStack);
		logInfo("Ran the cactus core script.\n");

		/*
		 * Cleanup
		 */
		destructList(subSequences);
		destructList(pairs);
		destructList(getAlignment_pairwiseAlignments);
		logInfo("Finished filling in the alignments for the net\n");
	}

	/*
	 * Write and close the netdisk.
	 */
	netDisk_write(netDisk);
	netDisk_destruct(netDisk);
	destructCactusCoreInputParameters(cCIP);
	logInfo("Finished with the net disk for this net.\n");

	return 0;
}

