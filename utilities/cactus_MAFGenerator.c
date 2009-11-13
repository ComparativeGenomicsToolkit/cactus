#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"

/*
 * The script outputs a maf file containing all the atom in a net and its descendants.
 */

/*
 * Global variables.
 */
bool includeTreesInMafBlocks = 0;

char *formatSequenceHeader(Sequence *sequence) {
	const char *sequenceHeader = sequence_getHeader(sequence);
	if(strlen(sequenceHeader) > 0) {
		char *cA = malloc(sizeof(char) *(1 + strlen(sequenceHeader)));
		sscanf(sequenceHeader, "%s", cA);
		return cA;
	}
	else {
		return netMisc_nameToString(sequence_getName(sequence));
	}
}

void getMAFBlock(Atom *atom, FILE *fileHandle) {
	/*
	 * Outputs a MAF representation of the atom to the given file handle.
	 */
	fprintf(fileHandle, "a score=%i\n", atom_getLength(atom) *atom_getInstanceNumber(atom));
	Atom_InstanceIterator *instanceIterator = atom_getInstanceIterator(atom);
	AtomInstance *atomInstance;
	while((atomInstance = atom_getNext(instanceIterator)) != NULL) {
		Sequence *sequence = atomInstance_getSequence(atomInstance);
		if(sequence != NULL) {
			char *sequenceHeader = formatSequenceHeader(sequence);
			int32_t start;
			if(atomInstance_getStrand(atomInstance)) {
				start = atomInstance_getStart(atomInstance) - sequence_getStart(sequence);
			}
			else { //start with respect to the start of the reverse complement sequence
				start = (sequence_getStart(sequence) + sequence_getLength(sequence) - 1) - atomInstance_getStart(atomInstance);
			}
			int32_t length = atomInstance_getLength(atomInstance);
			char *strand = atomInstance_getStrand(atomInstance) ? "+" : "-";
			int32_t sequenceLength = sequence_getLength(sequence);
			char *instanceString = atomInstance_getString(atomInstance);
			fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n", sequenceHeader, start, length, strand, sequenceLength, instanceString);
			free(instanceString);
			free(sequenceHeader);
		}
	}
	atom_destructInstanceIterator(instanceIterator);
}

void getMAFs(Net *net, FILE *fileHandle) {
	/*
	 * Outputs MAF representations of all the atom sin the net and its descendants.
	 */

	//Make MAF blocks for each atom
	Net_AtomIterator *atomIterator = net_getAtomIterator(net);
	Atom *atom;
	while((atom = net_getNextAtom(atomIterator)) != NULL) {
		getMAFBlock(atom, fileHandle);
	}
	net_destructAtomIterator(atomIterator);

	//Call child nets recursively.
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		getMAFs(adjacencyComponent_getNestedNet(adjacencyComponent), fileHandle); //recursive call.
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);
}

void makeMAFHeader(Net *net, FILE *fileHandle) {
	fprintf(fileHandle, "##maf version=1 scoring=NULL\n");
	char *cA = eventTree_makeNewickString(net_getEventTree(net));
	fprintf(fileHandle, "# cactus %s\n\n", cA);
	free(cA);
}

void usage() {
	fprintf(stderr, "cactus_mafGenerator, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --outputFile : The file to write the MAFs in.\n");
	fprintf(stderr, "-f --includeTrees : Include trees for each MAF block inside of a comment line.\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
	NetDisk *netDisk;
	Net *net;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * netDiskName = NULL;
	char * netName = NULL;
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
			{ "includeTrees", no_argument, 0, 'f' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:c:d:e:fh", long_options, &option_index);

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
			case 'd':
				netName = stringCopy(optarg);
				break;
			case 'e':
				outputFile = stringCopy(optarg);
				break;
			case 'f':
				includeTreesInMafBlocks = !includeTreesInMafBlocks;
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
		setLogLevel(LOGGING_INFO);
	}
	if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	//////////////////////////////////////////////
	//Log (some of) the inputs
	//////////////////////////////////////////////

	logInfo("Net disk name : %s\n", netDiskName);
	logInfo("Net name : %s\n", netName);
	logInfo("Output MAF file : %s\n", outputFile);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
	logInfo("Parsed the top level net of the cactus tree to check\n");

	///////////////////////////////////////////////////////////////////////////
	// Recursive check the nets.
	///////////////////////////////////////////////////////////////////////////

	int64_t startTime = time(NULL);
	FILE *fileHandle = fopen(outputFile, "w");
	makeMAFHeader(net, fileHandle);
	getMAFs(net, fileHandle);
	fclose(fileHandle);
	logInfo("Got the mafs in %i seconds/\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	netDisk_destruct(netDisk);

	return 0;
}
