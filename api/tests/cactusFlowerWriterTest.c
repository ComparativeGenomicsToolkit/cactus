/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static void testFlowerStream(CuTest *testCase) {
    CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
    char *tempPath = getTempFile();
    FILE *f = fopen(tempPath, "w");
    Name flowerNames[3];
    // First test flower
    Flower *flower1 = flower_construct(cactusDisk);
    flowerNames[0] = flower_getName(flower1);
    fprintf(f, "3 %" PRIi64, flowerNames[0]);
    // Second test flower
    Flower *flower2 = flower_construct(cactusDisk);
    flowerNames[1] = flower_getName(flower2);
    fprintf(f, " %" PRIi64, flowerNames[1] - flowerNames[0]);
    // Third test flower
    Flower *flower3 = flower_construct(cactusDisk);
    flowerNames[2] = flower_getName(flower3);
    fprintf(f, " %" PRIi64, flowerNames[2] - flowerNames[1]);
    fclose(f);

    // Ensure the flowers are serialized to disk, because
    // cactusDisk_getFlowers retrieves the records even if the flowers
    // are already loaded.
    cactusDisk_write(cactusDisk);
    flower_destruct(flower1, false);
    flower_destruct(flower2, false);
    flower_destruct(flower3, false);

    // Now read them back in.
    f = fopen(tempPath, "r");
    FlowerStream *flowerStream = flowerWriter_getFlowerStream(cactusDisk, f);
    CuAssertIntEquals(testCase, 3, flowerStream_size(flowerStream));
    int64_t i = 0;
    Flower *flower;
    while ((flower = flowerStream_getNext(flowerStream)) != NULL) {
        CuAssertTrue(testCase, i < 3);
        CuAssertIntEquals(testCase, flowerNames[i], flower_getName(flower));
        i++;
    }

    // Check that no flowers are loaded.
    CuAssertIntEquals(testCase, 0, stSortedSet_size(cactusDisk->flowers));
    flowerStream_destruct(flowerStream);
    fclose(f);
    removeTempFile(tempPath);
    testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
}

static void testFlowerWriter(CuTest *testCase) {
    char *tempFile = "./flowerWriterTest.txt";
    FILE *fileHandle = fopen(tempFile, "w");
    FlowerWriter *flowerWriter = flowerWriter_construct(fileHandle, 10, 5);
    flowerWriter_add(flowerWriter, 1, 5);
    flowerWriter_add(flowerWriter, 3, 5);
    flowerWriter_add(flowerWriter, 2, 5);
    flowerWriter_add(flowerWriter, 5, 4);
    flowerWriter_add(flowerWriter, 4, 5);
    flowerWriter_add(flowerWriter, 6, 1);
    flowerWriter_add(flowerWriter, 7, 9);
    flowerWriter_add(flowerWriter, 8, 1);
    flowerWriter_add(flowerWriter, 9, 1);
    flowerWriter_add(flowerWriter, 10, 1);
    flowerWriter_add(flowerWriter, 11, 1);
    flowerWriter_add(flowerWriter, 12, 7);

    flowerWriter_add(flowerWriter, 13, 1000);
    flowerWriter_destruct(flowerWriter);
    fclose(fileHandle);
    fileHandle = fopen(tempFile, "r");

    char *line = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "0 2 1 5 a 1 5 ", line);

    line = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "0 2 3 5 a 1 5 ", line);

    line = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "0 2 5 4 1 1 ", line);

    line = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "0 2 b 7 9 a 1 1 ", line);

    line = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "0 4 9 1 1 1 1 1 b 1 7 ", line);

    line = stFile_getLineFromFile(fileHandle);
    CuAssertStrEquals(testCase, "1 1 b 13 1000 ", line);

    CuAssertTrue(testCase, stFile_getLineFromFile(fileHandle) == NULL);

    stFile_rmtree(tempFile);
}

CuSuite* cactusFlowerWriterTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testFlowerStream);
    SUITE_ADD_TEST(suite, testFlowerWriter);
    return suite;
}
