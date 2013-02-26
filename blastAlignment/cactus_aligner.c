/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
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
#include "blastAlignmentLib.h"

void makeChunksFile(Flower *flower, char *outputFile, int32_t minimumSequenceLength, int32_t chunkOverlapLength) {
    FILE *headerFileHandle = fopen(outputFile, "w");
    int64_t chunkLength = 0;
    Flower_SequenceIterator *seqIt = flower_getSequenceIterator(flower);

}

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
	st_setLogLevelFromString(argv[1]);
	assert(argc == 6);
	stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[2]);
	cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
	st_logInfo("Set up the flower disk\n");

	flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(argv[3]));
	assert(flower != NULL);
	st_logInfo("Read the flower\n");

	int32_t minimumSequenceLength;
	int i = sscanf(argv[4], "%i", &minimumSequenceLength);
	(void)i;
	assert(i == 1);
	int32_t chunkOverlapLength;
    i = sscanf(argv[4], "%i", &chunkOverlapLength);
    (void)i;
    assert(i == 1);
	makeChunksFile(flower, argv[6], minimumSequenceLength, chunkOverlapLength);
	st_logInfo("Written the chunks header file\n");

	cactusDisk_destruct(cactusDisk);
	stKVDatabaseConf_destruct(kvDatabaseConf);

	return 0;
}
