/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "cactusReference.h"

char *getConsensusStringP(stList *strings, stList *outgroupStrings, int32_t blockLength);

static void testGetConsensusString(CuTest *testCase) {
    /*
     * Tests a maximum (cardinality) matching algorithm, checking that it has higher or equal
     * cardinality to the greedy algorithm.
     */
    stList *strings = stList_construct3(0, free);
    stList *outgroupStrings = stList_construct3(0, free);

    //Empty case
    char *consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "NNNNNNNNNN", consensus);
    free(consensus);
    //One sequence case
    stList_append(strings, stString_copy("ACTGNactgn"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "ACTGNactgn", consensus);
    free(consensus);
    //Two sequence case
    stList_append(strings, stString_copy("ACTGNactgn"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "ACTGNactgn", consensus);
    free(consensus);
    //Three sequence case
    stList_append(strings, stString_copy("CTGNactgnA"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "ACTGAactga", consensus);
    free(consensus);
    //five sequence case
    stList_append(strings, stString_copy("CTGNactgnA"));
    stList_append(strings, stString_copy("CTGNactgnA"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "CTGGactggA", consensus);
    free(consensus);
    stList_destruct(strings);
}

static void testGetConsensusStringWithOutgroups(CuTest *testCase) {
    /*
     * Tests a maximum (cardinality) matching algorithm, checking that it has higher or equal
     * cardinality to the greedy algorithm.
     */
    stList *strings = stList_construct3(0, free);
    stList *outgroupStrings = stList_construct3(0, free);
    stList_append(outgroupStrings, stString_copy("CTGNactgnA"));

    //Empty case
    char *consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "NNNNnnnnnN", consensus);
    free(consensus);
    //One sequence case
    stList_append(strings, stString_copy("ACTGNactgn"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "ACTGNactgN", consensus);
    free(consensus);
    //Two sequence case
    stList_append(strings, stString_copy("ACTGNactgn"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "ACTGNactgn", consensus);
    free(consensus);
    //Three sequence case
    stList_append(strings, stString_copy("CTGNactgnA"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "ACTGAactgA", consensus);
    free(consensus);
    //five sequence case
    stList_append(strings, stString_copy("CTGNactgnA")); //Second copy not needed, as outgroup breaks ties.
    //stList_append(strings, stString_copy("CTGNactgnA"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "CTGGactggA", consensus);
    free(consensus);
    stList_destruct(strings);
}

CuSuite* addReferenceCoordinatesTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testGetConsensusString);
    SUITE_ADD_TEST(suite, testGetConsensusStringWithOutgroups);
    return suite;
}
