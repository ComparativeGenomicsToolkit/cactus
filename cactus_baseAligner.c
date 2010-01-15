#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <getopt.h>

#include "hashTableC.h"
#include "bioioC.h"
#include "commonC.h"
#include "cactus.h"
#include "avl.h"
#include "pairwiseAlignment.h"
#include "cactus_misc.h"

void usage() {
	fprintf(stderr, "cactus_colinearAligner [net-names, ordered by order they should be processed], version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-e --tempDirRoot : The temp file root directory\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

char *getSequence(EndInstance *endInstance) {
	Sequence *sequence = endInstance_getSequence(endInstance);
	assert(sequence != NULL);
	EndInstance *endInstance2 = endInstance_getAdjacency(endInstance);
	assert(!endInstance_getSide(endInstance));
	if(endInstance_getStrand(endInstance)) {
		int32_t length = endInstance_getCoordinate(endInstance2)
					- endInstance_getCoordinate(endInstance) - 1;
		assert(length >= 0);
		return sequence_getString(sequence, endInstance_getCoordinate(
					endInstance) + 1, length, 1);
	}
	else {
		int32_t length = endInstance_getCoordinate(endInstance)
				- endInstance_getCoordinate(endInstance2) - 1;
		assert(length >= 0);
		return sequence_getString(sequence, endInstance_getCoordinate(
				endInstance2) + 1, length, 0);
	}
}

char *getTempSequenceFile(EndInstance *endInstance, char *tempDir) {
	assert(!endInstance_getSide(endInstance));
	char *string = getSequence(endInstance);
	//if(strlen(string) == 0) {
	//	free(string);
	//	return NULL;
	//}
	char *tempFile = pathJoin(tempDir, netMisc_nameToStringStatic(endInstance_getName(endInstance)));
	FILE *fileHandle = fopen(tempFile, "w");
	struct List *attributes = constructEmptyList(0, free);
	Sequence *sequence = endInstance_getSequence(endInstance);
	listAppend(attributes, netMisc_nameToString(sequence_getName(sequence)));
	listAppend(attributes, stringPrint("%i", endInstance_getStrand(endInstance)));
	listAppend(attributes, stringPrint("%i", endInstance_getCoordinate(endInstance) + (endInstance_getStrand(endInstance) ? 1 : -1)));
	char *header = fastaEncodeHeader(attributes);
	//uglyf("I am making the header: %s with endInstance coordinates: %i %i %i\n", header, endInstance_getCoordinate(endInstance), endInstance_getCoordinate(endInstance_getAdjacency(endInstance)), strlen(string));
	fprintf(fileHandle, ">%s\n%s\n", header, string);
	fclose(fileHandle);
	free(string);
	free(header);
	destructList(attributes);
	return tempFile;
}

void cleanupTempFile(char *tempFile) {
	char *command = stringPrint("rm -f %s", tempFile);
	exitOnFailure(system(command), "Failed to remove a temporary file\n");
	free(command);
	free(tempFile);
}

void alignSequences(char *tempSequenceFile1, char *tempSequenceFile2, char *tempCigarFile) {
	FILE *fileHandle = fopen(tempCigarFile, "w");
	fclose(fileHandle);
	//char *command = stringPrint("lastz --format=cigar %s[nameparse=darkspace] %s[nameparse=darkspace] --hspthresh=3000 --nogapped > %s",
	//			tempSequenceFile1, tempSequenceFile2, tempCigarFile);
	char *command = stringPrint("pecan2_pairwiseModel %s %s --cigars %s --matchThreshold 0.85", tempSequenceFile1, tempSequenceFile2, tempCigarFile);
	//uglyf("I would run the command: %s\n", command);
	exitOnFailure(system(command), "Failed to align the two sequences\n");
	//system(stringPrint("cat %s", tempSequenceFile1));
	//system(stringPrint("cat %s", tempSequenceFile2));
	//system(stringPrint("cat %s", tempCigarFile));
	free(command);
}

char **getTempSequences(EndInstance *endInstance1, char *tempDir) {
	assert(!endInstance_getSide(endInstance1));

	EndInstance *endInstance2 = endInstance_getAdjacency(endInstance1);

	assert(endInstance2 != NULL);
	assert(endInstance_getStrand(endInstance2));
	assert(endInstance_getSide(endInstance2));

	char **tempSequenceFiles = malloc(sizeof(char *) * 2);
	tempSequenceFiles[0] = getTempSequenceFile(endInstance1, tempDir);
	tempSequenceFiles[1] = getTempSequenceFile(endInstance_getReverse(endInstance2), tempDir);
	return tempSequenceFiles;
}

void destructTempSequences(char **tempSequenceFiles) {
	cleanupTempFile(tempSequenceFiles[0]);
	cleanupTempFile(tempSequenceFiles[1]);
	free(tempSequenceFiles);
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
	uint32_t j = (uint32_t)(RANDOM() * i);
	assert(j < i);
	void *item1 = items[j];
	j = (uint32_t)(i + (RANDOM() * (length - i)));
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

int main(int argc, char *argv[]) {

	char * logLevelString = NULL;
	char * netDiskName = NULL;
	char * tempDir = NULL;
	int32_t i, j;

	/*
	 * Parse the options.
	 */
	while(1) {
		static struct option long_options[] = {
				{ "logLevel", required_argument, 0, 'a' },
				{ "netDisk", required_argument, 0, 'b' },
				{ "tempDirRoot", required_argument, 0, 'e' },
				{ "help", no_argument, 0, 'h' },
				{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:b:e:h", long_options, &option_index);

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
			case 'e':
				tempDir = stringCopy(optarg);
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
	 * For each net.
	 */
	for (j = optind; j < argc; j++) {
		/*
		 * Load the netdisk
		 */

		NetDisk *netDisk = netDisk_construct(netDiskName);
		logInfo("Set up the net disk\n");

		/*
		 * Read the net.
		 */
		const char *netName = argv[j];
		logInfo("Processing the net named: %s\n", netName);
		Net *net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
		assert(net != NULL);
		logInfo("Parsed the net to be aligned\n");

		/*
		 * Make temp sequence files.
		 */
		Net_EndInstanceIterator *iterator1 = net_getEndInstanceIterator(net);
		EndInstance *endInstance1;
		struct List *tempSequenceFiles = constructEmptyList(0, (void (*)(void *))destructTempSequences);
		while((endInstance1 = net_getNextEndInstance(iterator1)) != NULL) {
			endInstance1 = endInstance_getStrand(endInstance1) ? endInstance1 : endInstance_getReverse(endInstance1);
			if(!endInstance_getSide(endInstance1)) {
				listAppend(tempSequenceFiles, getTempSequences(endInstance1, tempDir));
			}
		}
		net_destructEndInstanceIterator(iterator1);
		assert(tempSequenceFiles->length == net_getEndInstanceNumber(net)/2);
		logInfo("Got the temporary sequence files representing the sequences in to be aligned\n");

		/*
		 * Make a list of the pairwise comparisons to be aligned.
		 */
		struct List *pairs = constructEmptyList(0, NULL);
		int32_t spanningTrees = 2;
		for(i=0; i<spanningTrees; i++) {
			arrayShuffle(tempSequenceFiles->list, tempSequenceFiles->length);
			constructSpanningTree(tempSequenceFiles->list, tempSequenceFiles->length, pairs);
		}
		//check we have the expected number of pairs.
		assert(tempSequenceFiles->length == 0 ||
				(pairs->length/2 <= spanningTrees*(tempSequenceFiles->length-1) && (pairs->length/2 >= tempSequenceFiles->length-1)));
		logInfo("Got the list of pairs to align\n");

		/*
		 * Do the alignments.
		 */
		char *tempFile = pathJoin(tempDir, "temp.cigar");
		char *tempFile2 = pathJoin(tempDir, "temp2.cigar");
		FILE *fileHandle = fopen(tempFile, "w");
		for(i=0; i<pairs->length; i+=2) {
			char **tempSequenceFiles1 = pairs->list[i];
			char **tempSequenceFiles2 = pairs->list[i+1];

			char *tempSequenceFile1 = tempSequenceFiles1[0];
			char *tempSequenceFile2 = tempSequenceFiles1[1];
			char *tempSequenceFile3 = tempSequenceFiles2[0];

			//Do forward - forward alignment.
			alignSequences(tempSequenceFile1, tempSequenceFile2, tempFile2);
			convertCoordinates(tempFile2, fileHandle);
			//Do reverse - forward alignment.
			alignSequences(tempSequenceFile1, tempSequenceFile3, tempFile2);
			convertCoordinates(tempFile2, fileHandle);
		}
		fclose(fileHandle);
		logInfo("Created the alignments\n");

		/*
		 * Cleanup the temporary sequence files and other stuff.
		 */
		destructList(tempSequenceFiles);
		cleanupTempFile(tempFile2);
		destructList(pairs);
		logInfo("Cleaned up the temporary files\n");

		/*
		 * Close the netdisk.
		 */
		netDisk_destruct(netDisk);
		logInfo("Finished with the net disk for this net.\n");

		/*
		 * Run the cactus core script.
		 */
		char *command = stringPrint("cactus_core --alignments %s --netDisk %s --netName %s --tempDirRoot %s --maxEdgeDegree 10000000 --minimumTreeCoverage 0 --minimumTreeCoverageForAtoms 0 --minimumAtomLength 0 --minimumChainLength 0 --trim 0 --alignRepeats 1 --extensionSteps 0 --logLevel %s",
				tempFile, netDiskName, netName, tempDir, logLevelString);
		exitOnFailure(system(command), "Failed to run cactus_core when taking alignments from the colinear aligner");
		free(command);
		logInfo("Ran the cactus core script.");

		//assert(0);

		/*
		 * Cleanup
		 */
		cleanupTempFile(tempFile);
		logInfo("Finished filling in the alignments for the net\n");
	}

	return 0;
}

