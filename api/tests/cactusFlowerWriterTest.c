/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

static char getInt(FILE *fileHandle) {
    int32_t i;
    int32_t j = fscanf(fileHandle, "%i", &i);
    (void)j;
    assert(j == 1);
    return i;
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
    flowerWriter_add(flowerWriter, -1, 12);
    flowerWriter_destruct(flowerWriter);
    fclose(fileHandle);
    st_system("cat %s\n", tempFile);
    fileHandle = fopen(tempFile, "r");
    CuAssertIntEquals(testCase, 1, getInt(fileHandle));
    stList *flowers0 = flowerWriter_parseFlowerNames(fileHandle);
    CuAssertIntEquals(testCase, 0, getInt(fileHandle));
    stList *flowers1 = flowerWriter_parseFlowerNames(fileHandle);
    CuAssertIntEquals(testCase, 0, getInt(fileHandle));
    stList *flowers2 = flowerWriter_parseFlowerNames(fileHandle);
    CuAssertIntEquals(testCase, 0, getInt(fileHandle));
    stList *flowers3 = flowerWriter_parseFlowerNames(fileHandle);
    CuAssertIntEquals(testCase, 1, stList_length(flowers0));
    CuAssertIntEquals(testCase, -1, *(int64_t *)stList_get(flowers0, 0));
    CuAssertIntEquals(testCase, 2, stList_length(flowers1));
    CuAssertIntEquals(testCase, 1, *(int64_t *)stList_get(flowers1, 0));
    CuAssertIntEquals(testCase, 2, *(int64_t *)stList_get(flowers1, 1));
    CuAssertIntEquals(testCase, 2, stList_length(flowers2));
    CuAssertIntEquals(testCase, 3, *(int64_t *)stList_get(flowers2, 0));
    CuAssertIntEquals(testCase, 4, *(int64_t *)stList_get(flowers2, 1));
    CuAssertIntEquals(testCase, 2, stList_length(flowers3));
    CuAssertIntEquals(testCase, 5, *(int64_t *)stList_get(flowers3, 0));
    CuAssertIntEquals(testCase, 6, *(int64_t *)stList_get(flowers3, 1));
    stFile_rmrf(tempFile);
}

CuSuite* cactusFlowerWriterTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testFlowerWriter);
    return suite;
}
