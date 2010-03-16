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

char *startAlignmentStack_fileString;
static FILE *getNextAlignment_FileHandle = NULL;

static struct PairwiseAlignment *getNextAlignment() {
	return cigarRead(getNextAlignment_FileHandle);
}

static void startAlignmentStack() {
	if(getNextAlignment_FileHandle != NULL) {
		fclose(getNextAlignment_FileHandle);
	}
	getNextAlignment_FileHandle = fopen(startAlignmentStack_fileString, "r");
}

void usage() {
	fprintf(stderr, "cactus_core, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-b --alignments : The input alignments file\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --writeDebugFiles : Write the debug files\n");
	fprintf(stderr, "-h --help : Print this help screen\n");

	fprintf(stderr, "-i --alignRepeats : Allow bases marked as repeats to be aligned (else alignments to these bases to be excluded)\n");

	fprintf(stderr, "-j --alignUndoLoops : The number of rounds of alignment, undoing of over-aligned edges and recursion into adjacency connected components (groups)\n");
	fprintf(stderr, "-k --maxEdgeDegree : Maximum degree of aligned edges\n");

	fprintf(stderr, "-l --extensionSteps : The number of steps of attrition applied to edges proximal to over aligned edges.\n");
	fprintf(stderr, "-m --extensionStepsReduction : The reduction in the extension steps after each align/undo loop (to a minimum of zero).\n");

	fprintf(stderr, "-n --trim : The length of bases to remove from the end of each alignment\n");
	fprintf(stderr, "-o --trimReduction : Trim reduction, the amount to reduce the trim after each align/undo loop (to a minimum of zero)\n");

	fprintf(stderr, "-p --minimumTreeCoverage : Minimum tree coverage proportion of an block to be included in the graph\n");
	fprintf(stderr, "-q --minimumTreeCoverageReduction : Minimum tree coverage reduction after each align/undo loop\n");

	fprintf(stderr, "-r --minimumBlockLength : The minimum length of an block required to be included in the problem\n");

	fprintf(stderr, "-s --minimumChainLength : The minimum chain length required to be included in the problem\n");
	fprintf(stderr, "-t --minimumChainLengthReduction : The minimum chain length reduction after each align/undo loop\n");
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
			{ "writeDebugFiles", no_argument, 0, 'e' },
			{ "help", no_argument, 0, 'h' },

			{ "alignRepeats", no_argument, 0, 'i' },

			{ "alignUndoLoops", required_argument, 0, 'j' },
			{ "maxEdgeDegree", required_argument, 0, 'k' },

			{ "extensionSteps", required_argument, 0, 'l' },
			{ "extensionStepsReduction", required_argument, 0, 'm' },

			{ "trim", required_argument, 0, 'n' },
			{ "trimReduction", required_argument, 0, 'o',  },

			{ "minimumTreeCoverage", required_argument, 0, 'p' },
			{ "minimumTreeCoverageReduction", required_argument, 0, 'q' },

			{ "minimumBlockLength", required_argument, 0, 'r' },

			{ "minimumChainLength", required_argument, 0, 's' },
			{ "minimumChainLengthReduction", required_argument, 0, 't',  },

			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		key = getopt_long(argc, argv, "a:b:c:d:ehij:k:l:m:n:o:p:q:r:s:t:", long_options, &option_index);

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
			case 'e':
				cCIP->writeDebugFiles = 1;
				break;
			case 'h':
				usage();
				return 0;
			case 'i':
				cCIP->alignRepeats = !cCIP->alignRepeats;
				break;
			case 'j':
				assert(sscanf(optarg, "%i", &cCIP->alignUndoLoops) == 1);
				break;
			case 'k':
				assert(sscanf(optarg, "%i", &cCIP->maxEdgeDegree) == 1);
				break;
			case 'l':
				assert(sscanf(optarg, "%i", &cCIP->extensionSteps) == 1);
				break;
			case 'm':
				assert(sscanf(optarg, "%i", &cCIP->extensionStepsReduction) == 1);
				break;
			case 'n':
				assert(sscanf(optarg, "%i", &cCIP->trim) == 1);
				break;
			case 'o':
				assert(sscanf(optarg, "%i", &cCIP->trimReduction) == 1);
				break;
			case 'p':
				assert(sscanf(optarg, "%f", &cCIP->minimumTreeCoverage) == 1);
				break;
			case 'q':
				assert(sscanf(optarg, "%f", &cCIP->minimumTreeCoverageReduction) == 1);
				break;
			case 'r':
				assert(sscanf(optarg, "%i", &cCIP->minimumBlockLength) == 1);
				break;
			case 's':
				assert(sscanf(optarg, "%i", &cCIP->minimumChainLength) == 1);
				break;
			case 't':
				assert(sscanf(optarg, "%i", &cCIP->minimumChainLengthReduction) == 1);
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
	assert(cCIP->minimumTreeCoverageReduction >= 0.0);
	assert(cCIP->minimumBlockLength >= 0.0);
	assert(cCIP->minimumChainLength >= 0);
	assert(cCIP->minimumChainLengthReduction >= 0);
	assert(cCIP->trim >= 0);
	assert(cCIP->trimReduction >= 0);
	assert(cCIP->alignUndoLoops >= 0);
	assert(cCIP->extensionSteps >= 0);
	assert(cCIP->extensionStepsReduction >= 0);

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

	startAlignmentStack_fileString = alignmentsFile;
	exitOnFailure(cactusCorePipeline(net, cCIP, getNextAlignment, startAlignmentStack), "Failed to run the cactus core pipeline\n");
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
