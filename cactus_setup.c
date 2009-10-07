#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <inttypes.h>

#include "bioioC.h"
#include "cactus.h"

void usage() {
	fprintf(stderr, "cactus_setup [fastaFile]xN, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-c --netName : The name of the first net (the key in the database)\n");
	fprintf(stderr, "-e --uniqueNamePrefix : An alpha-numeric prefix which, when appended with any alpha-numeric characters is guaranteed to produce a unique name\n");
	fprintf(stderr, "-f --speciesTree : The species tree, which will form the skeleton of the event tree\n");
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

void fn(const char *fastaHeader, const char *string, int32_t length) {
	/*
	 * Processes a sequence by creating a file for it and adding it to the net.
	 */
	End *end1;
	End *end2;
	EndInstance *endInstance1;
	EndInstance *endInstance2;
	MetaSequence *metaSequence;
	Sequence *sequence;
	char *name;

	//Now put the details in a net.
	name = getUniqueName();
	metaSequence = metaSequence_construct(name, 1, length, string, fastaHeader,
			event_getName(event), netDisk);
	sequence = sequence_construct(metaSequence, net);
	end1 = end_construct(getUniqueStubName(), net);
	end2 = end_construct(getUniqueStubName(), net);
	endInstance1 = endInstance_construct2("0", end1, 1, 1, -1, sequence);
	endInstance2 = endInstance_construct2("0", end2, length, 1, 1, sequence);
	endInstance_makeAdjacent1(endInstance1, endInstance2);

	//cleanup
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
	char * speciesTree = NULL;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
				{ "logLevel", required_argument, 0, 'a' },
				{ "netDisk", required_argument, 0, 'b' },
				{ "netName", required_argument, 0, 'c' },
				{ "uniqueNamePrefix", required_argument, 0, 'e' },
				{ "speciesTree", required_argument, 0, 'f' },
				{ "help", no_argument, 0, 'h' },
				{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		key = getopt_long(argc, argv, "a:b:c:e:f:h", long_options, &option_index);

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
		case 'e':
			uniqueNamePrefix = optarg;
			break;
		case 'f':
			speciesTree = optarg;
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
	assert(uniqueNamePrefix != NULL);
	assert(speciesTree != NULL);

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

	binaryTree = newickTreeParser(speciesTree, 0.0, &strings);
	name = getUniqueName();
	eventTree = eventTree_construct(metaEvent_construct(name, "ROOT", netDisk), net); //creates the event tree and the root even
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
		if(binaryTree->internal) {
			event = event_construct(metaEvent_construct(name, NULL, netDisk), binaryTree->distance, event, eventTree);
			listAppend(stack, event);
			listAppend(stack, binaryTree->right);
			listAppend(stack, event);
			listAppend(stack, binaryTree->left);
		}
		else {
			assert(i < strings->length);
			originalName = strings->list[i];
			event = event_construct(metaEvent_construct(name, originalName, netDisk), binaryTree->distance, event, eventTree);
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
	for (i = optind; i < argc; i++) {
		logInfo("Sequence file %s\n", argv[i]);
		//put sequence in a special file.
		fastaFileName = argv[i];
		event = eventTree_getEvent(eventTree, strings->list[i]);
		fileHandle = fopen(fastaFileName, "r");
		fastaReadToFunction(fileHandle, fn);
		fclose(fileHandle);
	}
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
