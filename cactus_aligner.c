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

int main(int argc, char *argv[]) {
	/*
	 * Open the database.
	 * Open the net.
	 * Open a file to write the sequences.
	 * For each adjacency construct a sequence and put it in a fasta file.
	 * Finish!
	 */
	NetDisk *netDisk;
	Net *net;
	Net_EndIterator *endIterator;
	End_InstanceIterator *instanceIterator;
	End *end;
	EndInstance *endInstance;
	EndInstance *endInstance2;
	Sequence *sequence;
	char *string;
	FILE *fileHandle;

	assert(argc == 4);
	netDisk = netDisk_construct(argv[1]);
	logInfo("Set up the net disk\n");
	net = netDisk_getNet(netDisk, argv[2]);
	logInfo("Read the net\n");
	fileHandle = fopen(argv[3], "w");
	logInfo("Opened the file %s to write the sub-sequences in\n", argv[3]);
	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		instanceIterator = end_getInstanceIterator(end);
		while((endInstance = end_getNext(instanceIterator)) != NULL) {
			if(endInstance_getStrand(endInstance)) {
				endInstance2 = endInstance_getAdjacency(endInstance);
				assert(endInstance2 != NULL);
				sequence = endInstance_getSequence(endInstance);
				assert(sequence != NULL);
				string = sequence_getString(sequence, endInstance_getCoordinate(endInstance), endInstance_getCoordinate(endInstance2)-endInstance_getCoordinate(endInstance)-1, TRUE);
				fprintf(fileHandle, ">%s|%i\n%s\n", sequence_getName(sequence), endInstance_getCoordinate(endInstance)+1, string);
				free(string);
			}
		}
	}
	logInfo("Finished writing the subsequences to the file\n");
	fclose(fileHandle);
	return 0;
}
