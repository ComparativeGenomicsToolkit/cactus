/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

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

    stFile_rmrf(tempFile);
}

CuSuite* cactusFlowerWriterTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testFlowerWriter);
    return suite;
}
