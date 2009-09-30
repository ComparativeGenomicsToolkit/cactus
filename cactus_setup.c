#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

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
	fprintf(stderr, "-h --help : Print this help screen\n");
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
	NetDisk *netDisk;
	Net *net;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * netDiskName = NULL;
	char * netName = NULL;
	char * tempFileRootDirectory = NULL;
	char * uniqueNamePrefix = NULL;
	char * speciesTree = NULL;

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



	//////////////////////////////////////////////
	//Process the sequences
	//////////////////////////////////////////////

	//for each event..

	///////////////////////////////////////////////////////////////////////////
	// Write the net to disk.
	///////////////////////////////////////////////////////////////////////////

	netDisk_write(netDisk);
	logInfo("Updated the net on disk\n");

	return 0;
}
