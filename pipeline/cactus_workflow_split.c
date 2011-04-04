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
#include <inttypes.h>

#include "cactus.h"

int main(int argc, char *argv[]) {
    /*
     * This code iterates through the terminal groups and returns
     * a list of the new flowers.
     */
    assert(argc >= 4);
    const char *sequenceFileNameTo = "cactusSequences";

    stKVDatabaseConf *kvDatabaseConfFrom = stKVDatabaseConf_constructFromString(argv[1]);
    stKVDatabaseConf *kvDatabaseConfTo = stKVDatabaseConf_constructFromString(argv[2]);
    Name flowerName = cactusMisc_stringToName(argv[3]);

    CactusDisk *cactusDiskFrom = cactusDisk_construct(kvDatabaseConfFrom, 0);
    CactusDisk *cactusDiskTo = cactusDisk_construct3(kvDatabaseConfTo, 0, sequenceFileNameTo);

    Flower *flower = cactusDisk_getFlower(cactusDiskFrom, flowerName);
    cactusDisk_splitFlowers(flower, cactusDiskTo);

    char *sequenceFileNameFrom = stString_copy(cactusDisk_getSequenceFileName(cactusDiskFrom));

    cactusDisk_destruct(cactusDiskFrom);
    cactusDisk_destruct(cactusDiskTo);

    //Now link the sequences file
    int32_t i = st_system("ln %s/%s %s/%s", stKVDatabaseConf_getDir(kvDatabaseConfFrom), sequenceFileNameFrom,
            stKVDatabaseConf_getDir(kvDatabaseConfTo), sequenceFileNameTo);

    assert(i == 0);
    free(sequenceFileNameFrom);

    return i;
}

