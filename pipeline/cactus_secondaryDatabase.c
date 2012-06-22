/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"

int main(int argc, char *argv[]) {
    int32_t create;
    int32_t i = sscanf(argv[1], "%i", &create);
    assert(i == 1);
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[2]);
    if(create) {
         stKVDatabase *database = stKVDatabase_construct(kvDatabaseConf, 1);
         stKVDatabase_destruct(database);
    }
    else {
        stKVDatabase *database = stKVDatabase_construct(kvDatabaseConf, 0);
        stKVDatabase_deleteFromDisk(database);
    }
    return 0;
}
