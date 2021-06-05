/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions shared by the test code.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

char *testCommon_getTmpTestDir(const char *testName) {
    return stFile_pathJoin("test-output/tmp", testName);
}

Name testCommon_addThreadToFlower(Flower *flower, char *header, int64_t length) {
    char *dna = stRandom_getRandomDNAString(length, true, true, true);
    EventTree *eventTree = flower_getEventTree(flower);
    assert(eventTree != NULL);
    Sequence *sequence = sequence_construct(2, length, dna, header, eventTree_getRootEvent(eventTree), flower_getCactusDisk(flower));

    End *end1 = end_construct2(0, 0, flower);
    End *end2 = end_construct2(1, 0, flower);
    Cap *cap1 = cap_construct2(end1, 1, 1, sequence);
    Cap *cap2 = cap_construct2(end2, length + 2, 1, sequence);
    cap_makeAdjacent(cap1, cap2);

    free(dna);
    return cap_getName(cap1);
}
