#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "pinchGraph.h"
#include "cactusGraph.h"
#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "cactus.h"
#include "pairwiseAlignment.h"
#include "cactusNetFunctions.h"
#include "cactus_core.h"

FILE *getNextAlignment_FileHandle;

struct PairwiseAlignment *getNextAlignment() {
	return cigarRead(getNextAlignment_FileHandle);
}

void usage() {
	fprintf(stderr, "cactus_core, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --alignments : The input alignments file\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-f --maxEdgeDegree : Maximum degree of aligned edges\n");
	fprintf(stderr, "-g --writeDebugFiles : Write the debug files\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
	fprintf(stderr, "-i --minimumTreeCoverage : Minimum tree coverage proportion of an atom to be included in the graph\n");
	fprintf(stderr, "-o --minimumTreeCoverageForAtoms : Minimum tree coverage for atoms, proportion of an atom to be included in the set of atoms at this level\n");
	fprintf(stderr, "-j --minimumAtomLength : The minimum length of an atom required to be included in the problem\n");
	fprintf(stderr, "-k --minimumChainLength : The minimum chain length required to be included in the problem\n");
	fprintf(stderr, "-l --trim : The length of bases to remove from the end of each alignment\n");
	fprintf(stderr, "-m --alignRepeats : Allow bases marked as repeats to be aligned (else alignments to these bases to be excluded)\n");
	fprintf(stderr, "-n --extensionSteps : The number of steps of attrition applied to edges proximal to over aligned edges.\n");
}

int main(int argc, char *argv[]) {
	/*
	 * Script for adding alignments to cactus tree.
	 */
	int32_t startTime;
	NetDisk *netDisk;
	Net *net;
	int key;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * alignmentsFile = NULL;
	char * netDiskName = NULL;
	char * netName = NULL;
	CactusCoreInputParameters *cCIP = constructCactusCoreInputParameters();

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "alignments", required_argument, 0, 'b' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "netName", required_argument, 0, 'd' },
			{ "maxEdgeDegree", required_argument, 0, 'f' },
			{ "extensionSteps", required_argument, 0, 'n' },
			{ "writeDebugFiles", no_argument, 0, 'g' },
			{ "help", no_argument, 0, 'h' },
			{ "minimumTreeCoverage", required_argument, 0, 'i' },
			{ "minimumTreeCoverageForAtoms", required_argument, 0, 'o' },
			{ "minimumAtomLength", required_argument, 0, 'j' },
			{ "minimumChainLength", required_argument, 0, 'k' },
			{ "trim", required_argument, 0, 'l' },
			{ "alignRepeats", no_argument, 0, 'm' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		key = getopt_long(argc, argv, "a:b:c:d:f:ghi:j:k:l:mn:o:", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stringCopy(optarg);
				break;
			case 'b':
				alignmentsFile = stringCopy(optarg);
				break;
			case 'c':
				netDiskName = stringCopy(optarg);
				break;
			case 'd':
				netName = stringCopy(optarg);
				break;
			case 'f':
				assert(sscanf(optarg, "%i", &cCIP->maxEdgeDegree) == 1);
				break;
			case 'g':
				cCIP->writeDebugFiles = 1;
				break;
			case 'h':
				usage();
				return 0;
			case 'i':
				assert(sscanf(optarg, "%f", &cCIP->minimumTreeCoverage) == 1);
				break;
			case 'j':
				assert(sscanf(optarg, "%i", &cCIP->minimumAtomLength) == 1);
				break;
			case 'k':
				assert(sscanf(optarg, "%i", &cCIP->minimumChainLength) == 1);
				break;
			case 'l':
				assert(sscanf(optarg, "%i", &cCIP->trim) == 1);
				break;
			case 'm':
				cCIP->alignRepeats = !cCIP->alignRepeats;
				break;
			case 'n':
				assert(sscanf(optarg, "%i", &cCIP->extensionSteps) == 1);
				break;
			case 'o':
				assert(sscanf(optarg, "%f", &cCIP->minimumTreeCoverageForAtoms) == 1);
				break;
			default:
				usage();
				return 1;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// (0) Check the inputs.
	///////////////////////////////////////////////////////////////////////////

	assert(logLevelString == NULL || strcmp(logLevelString, "CRITICAL") == 0 || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
	assert(alignmentsFile != NULL);
	assert(netDiskName != NULL);
	assert(netName != NULL);
	assert(cCIP->maxEdgeDegree > 0);
	assert(cCIP->minimumTreeCoverage >= 0.0);
	assert(cCIP->minimumTreeCoverageForAtoms >= 0.0);
	assert(cCIP->minimumAtomLength >= 0.0);
	assert(cCIP->minimumChainLength >= 0);

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

	logInfo("Pairwise alignments file : %s\n", alignmentsFile);
	logInfo("Net disk name : %s\n", netDiskName);
	logInfo("Net name : %s\n", netName);
	logInfo("Max edge degree : %i\n", cCIP->maxEdgeDegree);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
	logInfo("Parsed the net to be refined\n");

	startTime = time(NULL);

	///////////////////////////////////////////////////////////////////////////
	// Call the core program.
	///////////////////////////////////////////////////////////////////////////

	getNextAlignment_FileHandle = fopen(alignmentsFile, "r");
	exitOnFailure(cactusCorePipeline(net, cCIP, getNextAlignment), "Failed to run the cactus core pipeline\n");
	fclose(getNextAlignment_FileHandle);

	///////////////////////////////////////////////////////////////////////////
	// (9) Write the net to disk.
	///////////////////////////////////////////////////////////////////////////

	netDisk_write(netDisk);
	logInfo("Updated the net on disk\n");

	///////////////////////////////////////////////////////////////////////////
	//(10) Clean up.
	///////////////////////////////////////////////////////////////////////////

	//Destruct stuff
	startTime = time(NULL);
	netDisk_destruct(netDisk);

	logInfo("Cleaned stuff up and am finished in: %i seconds\n", time(NULL) - startTime);
	return 0;
}
