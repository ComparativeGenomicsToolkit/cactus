/*
 * recursiveFileBuilderTest.c
 *
 *  Created on: 19 Mar 2012
 *      Author: benedictpaten
 */

#include <stdlib.h>

#include "sonLib.h"
#include "cactus.h"
#include "CuTest.h"
#include "recursiveThreadBuilder.h"

static char *writeSegment(Segment *segment) {
    return stString_print("%" PRIi64 " %s ", segment_getStart(segment), segment_getString(segment));
}

static char *writeTerminalAdjacency(Cap *cap) {
    if(cap_getCoordinate(cap_getAdjacency(cap)) - cap_getCoordinate(cap) - 1 == 0) {
        return stString_print("");
    }
    Sequence *sequence = cap_getSequence(cap);
    assert(sequence != NULL);
    return stString_print("%" PRIi64 " %s ", cap_getCoordinate(cap), sequence_getString(sequence, cap_getCoordinate(cap)+1, cap_getCoordinate(cap_getAdjacency(cap)) - cap_getCoordinate(cap) - 1, 1));
}

static void recursiveFileBuilder_test(CuTest *testCase) {
    //Make flower with two ends and 2 blocks, and one child, one empty adjacency and two containing additional blocks.

    const char *tempDir = "recursiveFileBuilderTestTempDir";
    if(stFile_exists(tempDir)) {
        stFile_rmtree(tempDir);
    }
    stFile_mkdir(tempDir);
    stKVDatabaseConf *conf = stKVDatabaseConf_constructTokyoCabinet(
                stFile_pathJoin(tempDir, "temporaryCactusDisk"));
    CactusDisk *cactusDisk = cactusDisk_construct(conf, true, true);
    eventTree_construct2(cactusDisk);
    Flower *flower = flower_construct(cactusDisk);
    End *end1 = end_construct2(0, 1, flower);
    End *end2 = end_construct2(1, 1, flower);

    //Make event tree
    Event *referenceEvent = eventTree_getRootEvent(flower_getEventTree(flower));

    //Make sequence and thread
    MetaSequence *metaSequence1 = metaSequence_construct(1, 5, "ACGTA", "ref sequence", event_getName(referenceEvent), cactusDisk);
    Sequence *sequence1 = sequence_construct(metaSequence1, flower);
    //First reference thread
    Cap *cap1 = cap_construct2(end1, 0, 1, sequence1);
    Cap *cap2 = cap_construct2(end2, 6, 1, sequence1);
    cap_makeAdjacent(cap1, cap2);

    //Make a group
    Group *group1 = group_construct2(flower);
    end_setGroup(end1, group1);
    end_setGroup(end2, group1);

    //Make nested flower
    Flower *nestedFlower = group_makeNestedFlower(group1);

    //Now will fill in blocks at lower level
    Block *block1 = block_construct(3, nestedFlower);
    Segment *segment1 = segment_construct2(block1, 1, 1, flower_getSequence(nestedFlower, sequence_getName(sequence1)));

    //Add adjacencies at lower level
    cap_makeAdjacent(flower_getCap(nestedFlower, cap_getName(cap1)), segment_get5Cap(segment1));
    cap_makeAdjacent(segment_get3Cap(segment1), flower_getCap(nestedFlower, cap_getName(cap2)));

    Group *nestedGroup = group_construct2(nestedFlower);
    End *end;
    Flower_EndIterator *endIt = flower_getEndIterator(nestedFlower);
    while((end = flower_getNextEnd(endIt)) != NULL) {
        end_setGroup(end, nestedGroup);
    }
    flower_destructEndIterator(endIt);

    //Create the sequence database
    stKVDatabaseConf *secondaryConf = stKVDatabaseConf_constructTokyoCabinet(
                    stFile_pathJoin(tempDir, "temporaryCactusDisk2"));
    stKVDatabase *secondaryDatabase = stKVDatabase_construct(secondaryConf, 1);
    stList *caps = stList_construct();
    stList_append(caps, flower_getCap(nestedFlower, cap_getName(cap1)));
    buildRecursiveThreads(secondaryDatabase, caps, writeSegment, writeTerminalAdjacency);
    stKVDatabase_destruct(secondaryDatabase);

    //Now complete the alignment
    secondaryDatabase = stKVDatabase_construct(secondaryConf, 0);
    stList_pop(caps);
    stList_append(caps, cap1);
    stList *threadStrings = buildRecursiveThreadsInList(secondaryDatabase, caps, writeSegment, writeTerminalAdjacency);
    stKVDatabase_deleteFromDisk(secondaryDatabase);

    CuAssertIntEquals(testCase, 1, stList_length(threadStrings));
    CuAssertStrEquals(testCase, "1 ACG 3 TA ", stList_get(threadStrings, 0));

    cactusDisk_destruct(cactusDisk);
    stFile_rmtree(tempDir);
}

CuSuite* recursiveThreadBuilderTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, recursiveFileBuilder_test);
    return suite;
}
