/*
 * main.c
 *
 *  Created on: 15-Apr-2010
 *      Author: benedictpaten
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "treeStats.h"

/*
 * Script gathers a whole gamut of statistics about the cactus/avg datastructure and reports them in an XML formatted document.
 */

void usage() {
	fprintf(stderr, "cactus_treeStats, version 0.1\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --outputFile : The file to write the stats in, XML formatted.\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
	/*
	 * The script builds a cactus tree representation of the chains and nets.
	 * The format of the output graph is dot format.
	 */
	CactusDisk *netDisk;
	Flower *net;
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
				logLevelString = stString_copy(optarg);
				break;
			case 'c':
				netDiskName = stString_copy(optarg);
				break;
			case 'd':
				netName = stString_copy(optarg);
				break;
			case 'e':
				outputFile = stString_copy(optarg);
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
		st_setLogLevel(ST_LOGGING_INFO);
	}
	if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
		st_setLogLevel(ST_LOGGING_DEBUG);
	}

	//////////////////////////////////////////////
	//Log (some of) the inputs
	//////////////////////////////////////////////

	st_logInfo("Net disk name : %s\n", netDiskName);
	st_logInfo("Net name : %s\n", netName);
	st_logInfo("Output graph file : %s\n", outputFile);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = cactusDisk_construct(netDiskName);
	st_logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = cactusDisk_getFlower(netDisk, cactusMisc_stringToName(netName));
	assert(net != NULL);
	st_logInfo("Parsed the top level net of the cactus tree to build\n");

	///////////////////////////////////////////////////////////////////////////
	// Calculate and print to file a crap load of numbers.
	///////////////////////////////////////////////////////////////////////////

	fileHandle = fopen(outputFile, "w");
	reportNetDiskStats(netDiskName, net, fileHandle);
	st_logInfo("Finished writing out the stats.\n");
	fclose(fileHandle);

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	cactusDisk_destruct(netDisk);

	return 0;
}
