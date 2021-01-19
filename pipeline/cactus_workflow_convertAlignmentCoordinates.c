/*
 * Copyright (C) 2009-2012 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "convertAlignmentCoordinates.h"

int main(int argc, char *argv[]) {
    /*
     * Gets the cigars, iterates through them converting them
     */
    assert(argc == 6);
    st_setLogLevelFromString(argv[1]);
    st_logDebug("Set up logging\n");

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[2]);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, false, true);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    convertAlignmentCoordinates(argv[3], argv[4], cactusDisk);

    cactusDisk_destruct(cactusDisk);

    return 0;
}
