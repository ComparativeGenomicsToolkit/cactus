#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"

void usage() {
	fprintf(stderr, "cactus_reference [net names], version 0.1\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

void makeTopLevelPseudoChromosomes(Net *net, Reference *reference) {
	/* Get the top level attached ends, ordered according to sequence */

	/* Make pseudo chromosome */
}

void makeIntermediateLevelPseudoChromosomes(Net *net, Reference *reference) {
	/* Make intermediate leve */
}

void makePsuedoAdjacencies(Net *net, Reference *reference) {

}

void addReferenceToNet(Net *net) {
	//Function will only work if no reference has already been added.
	assert(net_getReferenceNumber(net) == 0);
	Reference *reference = reference_construct(net);

	if(net_getParentGroup(net) == NULL) {
		/*
		 * If this is the top level net then we will create the pseudo-chromosomes based
		 * upon the set of attached stubs.. (we will throw an error ?! if we don't have at least one pair of
		 * attached stubs). We do this in the order of the sequences that were passed to us.
		 */
		makeTopLevelPseudoChromosomes(net, reference);
	}
	else {
		/*
		 * Else this is not the top level net, and we must locate the pseudo-adjacencies in the parent
		 * reference, to establish the ends of pseudo-chromosomes.
		 * We do this in the order of the parent reference's pseudo-adjacencies.
		 */
		makeIntermediateLevelPseudoChromosomes(net, reference);
	}
	/*
	 * Having defined the ordered pseudo chromosomes, we fill in the pseudo adjacencies.
	 * For each pseudo-chromosome..
	 */
	makePsuedoAdjacencies(net, reference);
}

int main(int argc, char *argv[]) {
	/*
	 * Script for adding a reference genome to a net.
	 */

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * netDiskName = NULL;
	int32_t j;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:c:h", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stringCopy(optarg);
				break;
			case 'c':
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

	///////////////////////////////////////////////////////////////////////////
	// (0) Check the inputs.
	///////////////////////////////////////////////////////////////////////////

	assert(logLevelString == NULL || strcmp(logLevelString, "CRITICAL") == 0 || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
	assert(netDiskName != NULL);

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

	logInfo("Netdisk name : %s\n", netDiskName);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	NetDisk *netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Loop on the nets, doing the reference genome.
	///////////////////////////////////////////////////////////////////////////

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
		 * Now run the reference function.
		 */
		addReferenceToNet(net);
	}

	///////////////////////////////////////////////////////////////////////////
	// Write the net(s) back to disk.
	///////////////////////////////////////////////////////////////////////////

	netDisk_write(netDisk);
	logInfo("Updated the net on disk\n");

	///////////////////////////////////////////////////////////////////////////
	//Clean up.
	///////////////////////////////////////////////////////////////////////////

	//Destruct stuff
	netDisk_destruct(netDisk);

	logInfo("Cleaned stuff up and am finished\n");
	return 0;
}
