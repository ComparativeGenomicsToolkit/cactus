#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <inttypes.h>

#include "bioioC.h"
#include "net.h"

void usage() {
	fprintf(stderr, "cactus_setup [fastaFile]xN, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-c --netName : The name of the first net (the key in the database)\n");
	fprintf(stderr, "-d --tempDirRoot : The temp file root directory\n");
	fprintf(stderr, "-e --uniqueNamePrefix : An alpha-numeric prefix which, when appended with any alpha-numeric characters is guaranteed to produce a unique name\n");
	fprintf(stderr, "-f --speciesTree : The species tree, which will form the skeleton of the event tree\n");
	fprintf(stderr, "-f --nameMapFile : The file in to which to write an XML formatted map of the names of sequences and tree events\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

/*
 * Plenty of global variables!
 */
char * netDiskName = NULL;
char * uniqueNamePrefix = NULL;
NetDisk *netDisk;
Net *net;
EventTree *eventTree;
Event *event;
char *fastaFileName;
FILE *nameMapFileHandle;
int32_t uniqueElementIndex=0;
int32_t uniqueSequenceNumber = 0;

char *getUniqueName() {
	/*
	 * Gets a unique string.
	 */
	char *cA;
	cA = malloc(sizeof(char) * strlen(uniqueNamePrefix) + 10);
	sprintf(cA, "%s%i", uniqueNamePrefix, uniqueElementIndex++);
	return cA;
}

char *getUniqueStubName() {
	/*
	 * Gets a unique string with a '$' appended.
	 */
	char *cA;
	char *name;
	cA = getUniqueName();
	name = malloc(sizeof(char)*(strlen(name)+2));
	sprintf(name, "$%s", cA);
	free(cA);
	return name;
}

char *getSequenceFile() {
	/*
	 * Gets a sequence file (a relative name), in which to put a sequence.
	 */
	char *cA;
	cA = malloc(sizeof(char)*20);
	sprintf(cA, "%i", uniqueSequenceNumber++);
	return cA;
}

void fn(const char *fastaHeader, const char *string, int32_t length) {
	/*
	 * Processes a sequence by creating a file for it and adding it to the net.
	 */
	End *end1;
	End *end2;
	EndInstance *endInstance1;
	EndInstance *endInstance2;
	Sequence *sequence;
	FILE *fileHandle;
	char *sequenceFile;
	char *absSequenceFile;
	char *name;

	//Write into a sequence file.
	sequenceFile = getSequenceFile();
	absSequenceFile = pathJoin(netDiskName, sequenceFile);
	fileHandle = fopen(absSequenceFile, "w");
	fprintf(fileHandle, "%s", string);
	fclose(fileHandle);

	//Now put the details in a net.
	name = getUniqueName();
	sequence = sequence_construct(name, length, sequenceFile, event, net);
	end1 = end_construct(getUniqueStubName(), net);
	end2 = end_construct(getUniqueStubName(), net);
	endInstance1 = endInstance_construct2("0", end1, 1, 1, -1, sequence);
	endInstance2 = endInstance_construct2("0", end2, length, 1, 1, sequence);
	endInstance_makeAdjacent1(endInstance1, endInstance2);

	//Write header into an xml doc.
	fprintf(nameMapFileHandle, "<sequence header=\"%s\", fileName=\"%s\" name=\"%s\"/>\n",
			fastaHeader, fastaFileName, name);

	//cleanup
	free(sequenceFile);
	free(absSequenceFile);
	free(name);
}

int main(int argc, char *argv[]) {
	/*
	 * Open the database.
	 * Construct a net.
	 * Construct an event tree representing the species tree.
	 * For each sequence contruct two ends each containing an end instance.
	 * Make a file for the sequence.
	 * Link the two end instances.
	 * Finish!
	 */

	int32_t key, i;
	char *originalName;
	char *name;
	struct List *strings;
	struct List *stack;
	struct BinaryTree *binaryTree;
	FILE *fileHandle;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * netName = NULL;
	char * tempFileRootDirectory = NULL;
	char * speciesTree = NULL;
	char * nameMapFileName = NULL;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
				{ "logLevel", required_argument, 0, 'a' },
				{ "netDisk", required_argument, 0, 'b' },
				{ "netName", required_argument, 0, 'c' },
				{ "tempDirRoot", required_argument, 0, 'd' },
				{ "uniqueNamePrefix", required_argument, 0, 'e' },
				{ "speciesTree", required_argument, 0, 'f' },
				{ "nameMapFile", required_argument, 0, 'g' },
				{ "help", no_argument, 0, 'h' },
				{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		key = getopt_long(argc, argv, "a:b:c:d:e:h", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
		case 'a':
			logLevelString = optarg;
			break;
		case 'b':
			netDiskName = optarg;
			break;
		case 'c':
			netName = optarg;
			break;
		case 'd':
			tempFileRootDirectory = optarg;
			break;
		case 'e':
			uniqueNamePrefix = optarg;
			break;
		case 'f':
			speciesTree = optarg;
			break;
		case 'g':
			nameMapFileName = optarg;
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
	assert(tempFileRootDirectory != NULL);
	assert(uniqueNamePrefix != NULL);
	assert(speciesTree != NULL);
	assert(nameMapFileHandle != NULL);

	//////////////////////////////////////////////
	//Set up logging
	//////////////////////////////////////////////

	if(strcmp(logLevelString, "INFO") == 0) {
		setLogLevel(LOGGING_INFO);
	}
	if(strcmp(logLevelString, "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	//////////////////////////////////////////////
	//Log (some of) the inputs
	//////////////////////////////////////////////

	logInfo("Net disk name : %s\n", netDiskName);
	logInfo("Net name : %s\n", netName);
	logInfo("Temp file root directory : %s\n", tempFileRootDirectory);
	logInfo("Unique name prefix : %s\n", uniqueNamePrefix);

	for (i = optind; i < argc; i++) {
	   logInfo("Sequence file %s\n", argv[i]);
	}

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	//////////////////////////////////////////////
	//Construct the net
	//////////////////////////////////////////////

	net = net_construct(netName, netDisk);
	logInfo("Constructed the net\n");

	//////////////////////////////////////////////
	//Construct the event tree
	//////////////////////////////////////////////

	nameMapFileHandle = fopen(nameMapFileName, "w");
	fprintf(nameMapFileHandle, "<name_map>\n<event_map>\n");
	binaryTree = newickTreeParser(speciesTree, 0.0, &strings);
	name = getUniqueName();
	eventTree = eventTree_construct(name, net); //creates the event tree and the root even
	//now traverse the tree
	stack = constructEmptyList(0, NULL);
	listAppend(stack, eventTree_getEvent(eventTree, name));
	listAppend(stack, binaryTree);
	free(name);
	i=0;
	while(stack->length > 0) {
		event = stack->list[--stack->length];
		binaryTree = stack->list[--stack->length];
		assert(binaryTree != NULL);
		name = getUniqueName();
		event = event_construct(name, binaryTree->distance, event, eventTree);
		if(binaryTree->internal) {
			listAppend(stack, event);
			listAppend(stack, binaryTree->right);
			listAppend(stack, event);
			listAppend(stack, binaryTree->left);
		}
		else {
			assert(i < strings->length);
			originalName = strings->list[i];
			fprintf(nameMapFileHandle, "<event original_name=\"%s\" event=\"%s\"></event>\n", originalName, name);
			free(originalName);
			strings->list[i++] = name;
		}
	}
	assert(i == strings->length);
	logInfo("Constructed the event tree\n");

	//////////////////////////////////////////////
	//Process the sequences
	//////////////////////////////////////////////

	assert(strings->length == argc - optind);
	fprintf(nameMapFileHandle, "</event_map>\n<fasta_map>\n");
	for (i = optind; i < argc; i++) {
		logInfo("Sequence file %s\n", argv[i]);
		//put sequence in a special file.
		fastaFileName = argv[i];
		event = eventTree_getEvent(eventTree, strings->list[i]);
		fileHandle = fopen(fastaFileName, "r");
		fastaReadToFunction(fileHandle, fn);
		fclose(fileHandle);
	}
	fprintf(nameMapFileHandle, "</fasta_map>\n</name_map>\n");
	fclose(nameMapFileHandle);
	logInfo("Construct the event graph\n");

	///////////////////////////////////////////////////////////////////////////
	// Write the net to disk.
	///////////////////////////////////////////////////////////////////////////

	netDisk_write(netDisk);
	logInfo("Updated the net on disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Cleanup.
	///////////////////////////////////////////////////////////////////////////

	destructList(strings);

	return 0;
}
