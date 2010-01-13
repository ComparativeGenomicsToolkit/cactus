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
	if(strlen(string) == 0) {
		free(string);
		return NULL;
	}
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
	char *command = stringPrint("pecan2_pairwiseModel %s %s --cigars %s --matchThreshold 0.9", tempSequenceFile1, tempSequenceFile2, tempCigarFile);
	//uglyf("I would run the command: %s\n", command);
	exitOnFailure(system(command), "Failed to align the two sequences\n");
	//system(stringPrint("cat %s", tempSequenceFile1));
	//system(stringPrint("cat %s", tempSequenceFile2));
	//system(stringPrint("cat %s", tempCigarFile));
	free(command);
}

int main(int argc, char *argv[]) {

	char * logLevelString = NULL;
	char * netDiskName = NULL;
	char * tempDir = NULL;
	int32_t j;

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

		if(net_getEndInstanceNumber(net) > 40) {
			netDisk_destruct(netDisk);
			return 0;
		}

		/*
		 * For every pair of sequences (both forward and reverse) generate cigars using our Ninja pairwise aligner.
		 */

		char *tempFile = pathJoin(tempDir, "temp.cigar");
		char *tempFile2 = pathJoin(tempDir, "temp2.cigar");
		FILE *fileHandle = fopen(tempFile, "w");
		Net_EndInstanceIterator *iterator1 = net_getEndInstanceIterator(net);
		EndInstance *endInstance1;
		while((endInstance1 = net_getNextEndInstance(iterator1)) != NULL) {
			endInstance1 = endInstance_getStrand(endInstance1) ? endInstance1 : endInstance_getReverse(endInstance1);
			EndInstance *endInstance2 = endInstance_getAdjacency(endInstance1);
			assert(endInstance2 != NULL);
			assert(endInstance_getStrand(endInstance2));

			if(!endInstance_getSide(endInstance1)) {
				assert(endInstance_getSide(endInstance2));
				char *tempSequenceFile1 = getTempSequenceFile(endInstance1, tempDir);
				char *tempSequenceFile2 = getTempSequenceFile(endInstance_getReverse(endInstance2), tempDir);

				//uglyf("Have subsequence: %s %i %i\n", netMisc_nameToStringStatic(sequence_getName(endInstance_getSequence(endInstance1))), endInstance_getCoordinate(endInstance1), endInstance_getCoordinate(endInstance2));
				if(tempSequenceFile1 == NULL) { //returns null is zero length
					assert(tempSequenceFile2 == NULL);
					continue;
				}

				EndInstance *endInstance3;
				Net_EndInstanceIterator *iterator2 = net_copyEndInstanceIterator(iterator1);
				while((endInstance3 = net_getNextEndInstance(iterator2)) != NULL) {
					endInstance3 = endInstance_getStrand(endInstance3) ? endInstance3 : endInstance_getReverse(endInstance3);

					if(!endInstance_getSide(endInstance3)) {
						char *tempSequenceFile3 = getTempSequenceFile(endInstance3, tempDir);
						if(tempSequenceFile3 == NULL) { //returns null if zero length.
							continue;
						}

						//Do forward - forward alignment.
						alignSequences(tempSequenceFile1, tempSequenceFile3, tempFile2);
						convertCoordinates(tempFile2, fileHandle);
						//Do reverse - forward alignment.
						alignSequences(tempSequenceFile2, tempSequenceFile3, tempFile2);
						convertCoordinates(tempFile2, fileHandle);

						cleanupTempFile(tempSequenceFile3);
					}
				}
				cleanupTempFile(tempSequenceFile1);
				cleanupTempFile(tempSequenceFile2);
				net_destructEndInstanceIterator(iterator2);
			}
		}
		cleanupTempFile(tempFile2);
		net_destructEndInstanceIterator(iterator1);
		fclose(fileHandle);
		logInfo("Created the alignments\n");

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
		command = stringPrint("rm %s", tempFile);
		exitOnFailure(system(command),
					"Something went wrong cleaning up the temporary alignment file");  //Cleanup the temporary sequence files.
		free(tempFile);
		logInfo("Finished filling in the alignments for the net\n");
	}

	return 0;
}

