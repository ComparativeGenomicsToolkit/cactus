/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"

int main(int argc, char *argv[]) {
    assert(argc == 3);
    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(argv[1]);
    int64_t create;
    int64_t i = sscanf(argv[2], "%" PRIi64 "", &create);
    (void)i;
    assert(i == 1);
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
