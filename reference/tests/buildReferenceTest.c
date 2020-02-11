/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "cactusReference.h"

static void constructEventTree_R(stTree *cur, EventTree *eventTree) {
    for (int64_t i = 0; i < stTree_getChildNumber(cur); i++) {
        stTree *child = stTree_getChild(cur, i);
        event_construct3(stTree_getLabel(child), stTree_getBranchLength(child),
                         eventTree_getEventByHeader(eventTree, stTree_getLabel(cur)),
                         eventTree);
        constructEventTree_R(child, eventTree);
    }
}

/*
 * Construct an event tree similar to the one cactus would build for a
 * given species tree.
 */
static EventTree *constructEventTree(const char *newick, Flower *flower) {
    stTree *tree = stTree_parseNewickString(newick);
    EventTree *eventTree = flower_getEventTree(flower);
    // Add the root of the species tree below the root of the event tree.
    event_construct3(stTree_getLabel(tree), stTree_getBranchLength(tree),
                     eventTree_getRootEvent(eventTree), eventTree);
    // Recurse on the species tree.
    constructEventTree_R(tree, eventTree);
    stTree_destruct(tree);

    return eventTree;
}

static void testEventWeighting(CuTest *testCase) {
    /*
     * Test that calculating the event weighting works correctly on a
     * simple (static) example.
     */
    CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
    eventTree_construct2(cactusDisk);
    Flower *flower = flower_construct(cactusDisk);

    EventTree *eventTree = constructEventTree("((simHuman_chr6:0.144,(simMouse_chr6:0.0845,simRat_chr6:0.0916)Anc3:0.272)Anc2:0.0206,(simCow_chr6:0.1891,simDog_chr6:0.163)Anc1:0.0329)Anc0;", flower);
    Event *referenceEvent = eventTree_getEventByHeader(eventTree, "Anc2");
    stSet *chosenEvents = stSet_construct();
    stSet_insert(chosenEvents,
                 eventTree_getEventByHeader(eventTree, "simHuman_chr6"));
    stSet_insert(chosenEvents,
                 eventTree_getEventByHeader(eventTree, "Anc3"));
    stSet_insert(chosenEvents,
                 eventTree_getEventByHeader(eventTree, "simCow_chr6"));
    stSet_insert(chosenEvents,
                 eventTree_getEventByHeader(eventTree, "simDog_chr6"));

    double phi = st_random();
    stHash *weights = getEventWeighting(referenceEvent, phi, chosenEvents);

    // Check that the ingroups' path lengths weren't adjusted.

    Event *event = eventTree_getEventByHeader(eventTree, "simHuman_chr6");
    // True path length: 0.144, adjusted path length: 0.144
    double trueScore = exp(-phi * 0.144);
    CuAssertDblEquals(testCase, trueScore, stDoubleTuple_getPosition(stHash_search(weights, event), 0), 0.001);

    event = eventTree_getEventByHeader(eventTree, "Anc3");
    // True path length: 0.272, adjusted path length: 0.272
    trueScore = exp(-phi * 0.272);
    CuAssertDblEquals(testCase, trueScore, stDoubleTuple_getPosition(stHash_search(weights, event), 0), 0.001);

    // Check that the outgroups' path lengths were adjusted to take
    // into account the multiplicity.

    event = eventTree_getEventByHeader(eventTree, "simCow_chr6");
    // True path length: 0.2426, adjusted path length: 0.0206/2 + 0.0329/2 + 0.1891 = 0.21585
    trueScore = exp(-phi * 0.2426) * 0.21585/0.2426;
    CuAssertDblEquals(testCase, trueScore, stDoubleTuple_getPosition(stHash_search(weights, event), 0), 0.001);

    testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
    stHash_destruct(weights);
    stSet_destruct(chosenEvents);
}

CuSuite* buildReferenceTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testEventWeighting);
    return suite;
}
