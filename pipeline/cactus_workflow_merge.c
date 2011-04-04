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
     * This code merges a flower and its descendants and puts them into a back into a cactus disk.
     */
    assert(argc >= 4);

    stKVDatabaseConf *kvDatabaseConfFrom = stKVDatabaseConf_constructFromString(argv[1]);
    stKVDatabaseConf *kvDatabaseConfTo = stKVDatabaseConf_constructFromString(argv[2]);
    Name flowerName = cactusMisc_stringToName(argv[3]);

    CactusDisk *cactusDiskFrom = cactusDisk_construct(kvDatabaseConfFrom, 0);
    CactusDisk *cactusDiskTo = cactusDisk_construct3(kvDatabaseConfTo, 0, sequenceFileNameTo);

    Flower *flower = cactusDisk_getFlower(cactusDiskFrom, flowerName);
    cactusDisk_mergeFlowers(flower, cactusDiskTo);

    return 0;
}

