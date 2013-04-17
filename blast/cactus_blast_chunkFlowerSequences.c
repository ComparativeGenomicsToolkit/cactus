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
	assert(argc == 8);
	st_setLogLevelFromString(argv[1]);
	stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[2]);
	cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
	st_logInfo("Set up the flower disk\n");
	flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(argv[3]));
	assert(flower != NULL);
	st_logInfo("Read the flower\n");
	int64_t chunkSize, chunkOverlapSize, minimumSequenceLength;
	int64_t i = sscanf(argv[4], "%" PRIi64 "", &chunkSize);
	(void)i;
    assert(i == 1);
    i = sscanf(argv[5], "%" PRIi64 "", &chunkOverlapSize);
    assert(i == 1);
	i = sscanf(argv[6], "%" PRIi64 "", &minimumSequenceLength);
	assert(i == 1);
	setupToChunkSequences(chunkSize, chunkOverlapSize, argv[7]);
	writeFlowerSequences(flower, processSequenceToChunk, minimumSequenceLength);
	finishChunkingSequences();
	st_logInfo("Written the sequences from the flower into a file");
	cactusDisk_destruct(cactusDisk);
	stKVDatabaseConf_destruct(kvDatabaseConf);

	return 0;
}
