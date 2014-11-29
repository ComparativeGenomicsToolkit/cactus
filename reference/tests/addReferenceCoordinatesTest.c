/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>
#include "CuTest.h"
#include "sonLib.h"
#include "cactus.h"
#include "blockConsensusString.h"
#include "blockMLString.h"

char *getConsensusStringP(stList *strings, stList *outgroupStrings, int64_t blockLength);

static void testGetConsensusString(CuTest *testCase) {
    /*
     * Tests a maximum (cardinality) matching algorithm, checking that it has higher or equal
     * cardinality to the greedy algorithm.
     */
    stList *strings = stList_construct3(0, free);
    stList *outgroupStrings = stList_construct3(0, free);

    //Empty case
    char *consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "nnnnnnnnnn", consensus);
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

static void checkTree(CuTest *testCase, stTree *tree, stSet *eventsSet) {
    /*
     * Used to check the phylogenetic tree we construct for base calling.
     */
    Event *event = getEvent(tree);
    stSet_insert(eventsSet, event); //Increase the set of observed events
    stSet *connectedEvents = stSet_construct();
    if(stTree_getParent(tree) != NULL) {
        stSet_insert(connectedEvents, stTree_getParent(tree));
    }
    for(int64_t i=0; i<stTree_getChildNumber(tree); i++) {
        checkTree(testCase, stTree_getChild(tree, i), eventsSet);
        stSet_insert(connectedEvents, stTree_getChild(tree, i));
    }
    //This checks the number of connected events/nodes is equal
    CuAssertIntEquals(testCase, event_getChildNumber(event) + (event_getParent(event) != NULL), stSet_size(connectedEvents));
    //Now check that all the events connected in the tree
    for(int64_t i=0; i<event_getChildNumber(event); i++) {
        CuAssertTrue(testCase, stSet_search(connectedEvents, event_getChild(event, i)) != NULL);
    }
    if(event_getParent(event) != NULL) {
        CuAssertTrue(testCase, stSet_search(connectedEvents, event_getParent(event)) != NULL);
    }
    stSet_destruct(connectedEvents); //Cleanup loop
}

static void testMLStringRandom(CuTest *testCase) {
    for(int64_t i=0; i<100; i++) {
        CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk();
        Flower *flower = flower_construct(cactusDisk);
        //Make a random eventTree.
        stList *events = stList_construct();
        stList_append(events, eventTree_getRootEvent(flower_getEventTree(flower)));
        while(st_random() > 0.2) {
            stList_append(events, event_construct3("Boo", st_random(), st_randomChoice(events), flower_getEventTree(flower)));
        }
        //Make a random block
        Block *block = block_construct(st_randomInt(1, 100), flower);
        stList *strings = stList_construct3(0, free);
        while(st_random() > 0.2) {
            //Make a random sequence
            MetaSequence *metaSeq = metaSequence_construct(0, block_getLength(block),
                    stRandom_getRandomDNAString(block_getLength(block), 1, 0, 1),
                    "boo", event_getName(st_randomChoice(events)), cactusDisk);
            Sequence *seq = sequence_construct(metaSeq, flower);
            //Add string to set
            stList_append(strings, sequence_getString(seq, 0, block_getLength(block), 1));
            //Now make the random segment
            segment_construct2(block, 0, 1, seq);
        }
        //Make the random phylogenetic tree
        Event *refEvent = st_randomChoice(events);
        stTree *tree = getPhylogeneticTreeRootedAtGivenEvent(refEvent, generateJukesCantorMatrix);

        //Check the tree is rooted at the chosen event and is otherwise isomorphic to the random eventTree
        CuAssertTrue(testCase, getEvent(tree) == refEvent); //The returned node represents the chosen event
        CuAssertTrue(testCase, stTree_getParent(tree) == NULL); //It is at the root of the tree.
        CuAssertTrue(testCase, stTree_getNumNodes(tree) == stList_length(events)); //Check the number of nodes in the tree is equal to the number of events
        stSet *eventsSet = stSet_construct();
        checkTree(testCase, tree, eventsSet);
        CuAssertTrue(testCase, stList_length(events) == stSet_size(eventsSet));
        stSet_destruct(eventsSet);

        //Now create the ML string
        char *mlString = getMaximumLikelihoodString(tree, block);

        //Check the ML string has the right length, that each base is valid.
        CuAssertIntEquals(testCase, strlen(mlString), block_getLength(block));
        for(int64_t i=0; i<block_getLength(block); i++) {
            switch(mlString[i]) {
            case 'a':
            case 'A':
            case 'c':
            case 'C':
            case 'g':
            case 'G':
            case 't':
            case 'T':
            case 'n':
            case 'N':
                break;
            default:
                CuAssertTrue(testCase, 0);
            }
        }

        //Check the ML string is repeat masked as we would expect
        for(int64_t i=0; i<block_getLength(block); i++) {
            int64_t upperCount = 0;
            for(int64_t j=0; j<stList_length(strings); j++) {
                char c = ((char *)stList_get(strings, j))[i];
                upperCount += (c == toupper(c));
            }
            if(upperCount > stList_length(strings)/2) {
                CuAssertTrue(testCase, mlString[i] == toupper(mlString[i]));
            }
            else {
                CuAssertTrue(testCase, mlString[i] == tolower(mlString[i]));
            }
        }

        //Cleanup
        free(mlString);
        cleanupPhylogeneticTree(tree);
        stList_destruct(events);
        testCommon_deleteTemporaryCactusDisk(cactusDisk);
    }
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
    CuAssertStrEquals(testCase, "CTGNactgnA", consensus);
    free(consensus);
    //One sequence case
    stList_append(strings, stString_copy("ACTGNactgn"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "ACTGaactga", consensus);
    free(consensus);
    //Two sequence case
    stList_append(strings, stString_copy("ACTGNactgn"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "ACTGAactga", consensus);
    free(consensus);
    //Three sequence case
    stList_append(strings, stString_copy("CTGNactgnA"));
    consensus = getConsensusStringP(strings, outgroupStrings, 10);
    CuAssertStrEquals(testCase, "ACTGaactga", consensus);
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
    SUITE_ADD_TEST(suite, testMLStringRandom);

    return suite;
}
