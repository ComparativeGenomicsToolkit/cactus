/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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
	 * Open the flower.
	 * Open a file to write the sequences.
	 * For each adjacency construct a sequence and put it in a fasta file.
	 * Finish!
	 */
	CactusDisk *cactusDisk;
	Flower *flower;
	Flower_EndIterator *endIterator;
	End_InstanceIterator *instanceIterator;
	End *end;
	Cap *cap;
	Cap *cap2;
	Sequence *sequence;
	char *string;
	FILE *fileHandle;
	st_setLogLevel(ST_LOGGING_DEBUG);
	assert(argc == 4);
	stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[1]);
	cactusDisk = cactusDisk_construct2(kvDatabaseConf, 0, 1);
	st_logInfo("Set up the flower disk\n");

	flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(argv[2]));
	assert(flower != NULL);
	st_logInfo("Read the flower\n");
	fileHandle = fopen(argv[3], "w");
	st_logInfo("Opened the file %s to write the sub-sequences in\n", argv[3]);
	endIterator = flower_getEndIterator(flower);
	while((end = flower_getNextEnd(endIterator)) != NULL) {
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
				fprintf(fileHandle, ">%s|1|%i\n%s\n", cactusMisc_nameToStringStatic(sequence_getName(sequence)), cap_getCoordinate(cap)+1, string);
				free(string);
			}
		}
		end_destructInstanceIterator(instanceIterator);
	}
	st_logInfo("Finished writing the subsequences to the file\n");
	fclose(fileHandle);
	flower_destructEndIterator(endIterator);
	cactusDisk_destruct(cactusDisk);
	stKVDatabaseConf_destruct(kvDatabaseConf);

	return 0;
}
