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
	Cap *cap;
	Cap *cap2;
	Sequence *sequence;
	char *string;
	FILE *fileHandle;
	st_setLogLevel(ST_LOGGING_DEBUG);
	assert(argc == 4);
	netDisk = netDisk_construct(argv[1]);
	st_logInfo("Set up the net disk\n");

	net = netDisk_getNet(netDisk, netMisc_stringToName(argv[2]));
	st_logInfo("Read the net\n");
	fileHandle = fopen(argv[3], "w");
	st_logInfo("Opened the file %s to write the sub-sequences in\n", argv[3]);
	endIterator = net_getEndIterator(net);
	while((end = net_getNextEnd(endIterator)) != NULL) {
		instanceIterator = end_getInstanceIterator(end);
		while((cap = end_getNext(instanceIterator)) != NULL) {
			cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
			cap2 = cap_getAdjacency(cap);
			assert(cap2 != NULL);
			assert(cap_getStrand(cap2));

			if(!cap_getSide(cap)) {
				assert(cap_getSide(cap2));
				assert(cap_getCoordinate(cap2)-cap_getCoordinate(cap)-1 >= 0);
				sequence = cap_getSequence(cap);
				assert(sequence != NULL);
				string = sequence_getString(sequence, cap_getCoordinate(cap)+1, cap_getCoordinate(cap2)-cap_getCoordinate(cap)-1, TRUE);
				fprintf(fileHandle, ">%s|1|%i\n%s\n", netMisc_nameToStringStatic(sequence_getName(sequence)), cap_getCoordinate(cap)+1, string);
				free(string);
			}
		}
	}
	st_logInfo("Finished writing the subsequences to the file\n");
	fclose(fileHandle);
	return 0;
}
