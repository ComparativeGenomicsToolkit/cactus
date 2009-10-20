#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <inttypes.h>
#include <sys/types.h>
#include <dirent.h>

#include "bioioC.h"
#include "cactus.h"

void usage() {
	fprintf(stderr, "cactus_setup [fastaFile]xN, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-f --speciesTree : The species tree, which will form the skeleton of the event tree\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

/*
 * Plenty of global variables!
 */
char * netDiskName = NULL;
NetDisk *netDisk;
Net *net;
EventTree *eventTree;
Event *event;

void fn(const char *fastaHeader, const char *string, int32_t length) {
	/*
	 * Processes a sequence by adding it to the net disk.
	 */
	End *end1;
	End *end2;
	EndInstance *endInstance1;
	EndInstance *endInstance2;
	MetaSequence *metaSequence;
	Sequence *sequence;
	char *name;

	//Now put the details in a net.
	metaSequence = metaSequence_construct(1, length, string, fastaHeader,
			event_getName(event), netDisk);
	sequence = sequence_construct(metaSequence, net);
	end1 = end_construct(1, net);
	end2 = end_construct(1, net);
	endInstance1 = endInstance_construct2(end1, 1, 1, -1, sequence);
	endInstance2 = endInstance_construct2(end2, length, 1, 1, sequence);
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

	int32_t key, i, j;
	MetaEvent *metaEvent;
	struct List *strings;
	struct List *stack;
	struct BinaryTree *binaryTree;
	FILE *fileHandle;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * speciesTree = NULL;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
				{ "logLevel", required_argument, 0, 'a' },
				{ "netDisk", required_argument, 0, 'b' },
				{ "speciesTree", required_argument, 0, 'f' },
				{ "help", no_argument, 0, 'h' },
				{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		key = getopt_long(argc, argv, "a:b:f:h", long_options, &option_index);

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

	for (i = optind; i < argc; i++) {
	   logInfo("Sequence file/directory %s\n", argv[i]);
	}

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	//////////////////////////////////////////////
	//Construct the net
	//////////////////////////////////////////////

	net = net_construct(netDisk);
	logInfo("Constructed the net\n");

	//////////////////////////////////////////////
	//Construct the event tree
	//////////////////////////////////////////////

	binaryTree = newickTreeParser(speciesTree, 0.0, &strings);
	metaEvent = metaEvent_construct("ROOT", netDisk);
	eventTree = eventTree_construct(metaEvent, net); //creates the event tree and the root even
	logInfo("Constructed the basic event tree\n");

	//now traverse the tree
	stack = constructEmptyList(0, NULL);
	listAppend(stack, eventTree_getRootEvent(eventTree));
	listAppend(stack, binaryTree);
	i=0;
	j=optind;
	while(stack->length > 0) {
		binaryTree = stack->list[--stack->length];
		event = stack->list[--stack->length];
		assert(binaryTree != NULL);
		if(binaryTree->internal) {
			assert(i < strings->length);
			event = event_construct(metaEvent_construct(strings->list[i++], netDisk), binaryTree->distance, event, eventTree);
			listAppend(stack, event);
			listAppend(stack, binaryTree->right);
			listAppend(stack, event);
			listAppend(stack, binaryTree->left);
		}
		else {
			assert(i < strings->length);
			uglyf("Number of stirngs %i\n", strings->length);
			assert(j < argc);
			event = event_construct(metaEvent_construct(strings->list[i++], netDisk), binaryTree->distance, event, eventTree);
			DIR *dh;//directory handle
			struct dirent *file;//a 'directory entity' AKA file
			dh=opendir(argv[j++]);
			while((file=readdir(dh)) != NULL) {
			    printf("%s\n",file->d_name);
			}
			//fileHandle = fopen(argv[i++], "r");
			//fastaReadToFunction(fileHandle, fn);
			//fclose(fileHandle);
		}
	}
	assert(i == strings->length);
	logInfo("Constructed the initial net\n");

	///////////////////////////////////////////////////////////////////////////
	// Write the net to disk.
	///////////////////////////////////////////////////////////////////////////

	netDisk_write(netDisk);
	logInfo("Updated the net on disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Cleanup.
	///////////////////////////////////////////////////////////////////////////

	destructList(strings);
	netDisk_destruct(netDisk);

	return 0;
}
