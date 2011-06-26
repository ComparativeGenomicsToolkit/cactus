/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "cactusReference.h"

char *getConsensusStringP(stList *strings, int32_t blockLength);

static void testGetConsensusString(CuTest *testCase) {
    /*
     * Tests a maximum (cardinality) matching algorithm, checking that it has higher or equal
     * cardinality to the greedy algorithm.
     */
    stList *strings = stList_construct3(0, free);
    //Empty case
    char *consensus = getConsensusStringP(strings, 10);
    CuAssertStrEquals(testCase, "NNNNNNNNNN", consensus);
    free(consensus);
    //One sequence case
    stList_append(strings, stString_copy("ACTGNactgn"));
    consensus = getConsensusStringP(strings, 10);
    CuAssertStrEquals(testCase, "ACTGNactgn", consensus);
    free(consensus);
    //Two sequence case
    stList_append(strings, stString_copy("ACTGNactgn"));
    consensus = getConsensusStringP(strings, 10);
    CuAssertStrEquals(testCase, "ACTGNactgn", consensus);
    free(consensus);
    //Three sequence case
    stList_append(strings, stString_copy("CTGNactgnA"));
    consensus = getConsensusStringP(strings, 10);
    CuAssertStrEquals(testCase, "ACTGNactgn", consensus);
    free(consensus);
    //five sequence case
    stList_append(strings, stString_copy("CTGNactgnA"));
    stList_append(strings, stString_copy("CTGNactgnA"));
    consensus = getConsensusStringP(strings, 10);
    CuAssertStrEquals(testCase, "CTGNactgnA", consensus);
    free(consensus);
    stList_destruct(strings);
}

CuSuite* addReferenceCoordinatesTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testGetConsensusString);
    return suite;
}
