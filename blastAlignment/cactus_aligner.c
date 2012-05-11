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
	st_setLogLevelFromString(argv[5]);
	assert(argc == 6);
	stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[1]);
	cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
	st_logInfo("Set up the flower disk\n");

	flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(argv[2]));
	assert(flower != NULL);
	st_logInfo("Read the flower\n");

	int32_t minimumSequenceLength;
	int i = sscanf(argv[4], "%i", &minimumSequenceLength);
	(void)i;
	assert(i == 1);
	writeFlowerSequencesInFile(flower, argv[3], minimumSequenceLength);
	st_logInfo("Written the sequences from the flower into a file");

	cactusDisk_destruct(cactusDisk);
	stKVDatabaseConf_destruct(kvDatabaseConf);

	return 0;
}
