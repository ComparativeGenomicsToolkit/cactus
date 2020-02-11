/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <ctype.h>
#include "CuTest.h"
#include "sonLib.h"
#include "cactus.h"
#include "blockMLString.h"

static void checkTree(CuTest *testCase, stTree *tree, stSet *eventsSet) {
    /*
     * Used to check the phylogenetic tree we construct for base calling.
     */
    Event *event = getEvent(tree);
    stSet_insert(eventsSet, event); //Increase the set of observed events
    stSet *connectedEvents = stSet_construct();
    if(stTree_getParent(tree) != NULL) {
        stSet_insert(connectedEvents, getEvent(stTree_getParent(tree)));
    }
    for(int64_t i=0; i<stTree_getChildNumber(tree); i++) {
        checkTree(testCase, stTree_getChild(tree, i), eventsSet);
        stSet_insert(connectedEvents, getEvent(stTree_getChild(tree, i)));
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
        CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
        eventTree_construct2(cactusDisk);
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
        //and composed of Ns where there is insufficient info.
        for(int64_t i=0; i<block_getLength(block); i++) {
            int64_t upperCount = 0;
            int64_t nCount = 0;
            for(int64_t j=0; j<stList_length(strings); j++) {
                char c = ((char *)stList_get(strings, j))[i];
                upperCount += (c == toupper(c));
                nCount += (toupper(c) == 'N');
            }
            if(upperCount > stList_length(strings)/2) {
                CuAssertTrue(testCase, mlString[i] == toupper(mlString[i]));
            }
            else {
                CuAssertTrue(testCase, mlString[i] == tolower(mlString[i]));
            }
            if(nCount == stList_length(strings)) {
                CuAssertTrue(testCase, toupper(mlString[i]) == 'N');
            }
        }

        //Cleanup
        free(mlString);
        cleanupPhylogeneticTree(tree);
        stList_destruct(events);
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
    }
}

static void testMLStringMakesScaffoldGaps(CuTest *testCase) {
    /*
     * Simply test that scaffold gaps created where the reference does
     * not have direct sequence support for an adjacency are called as Ns in
     * the ML base-calling code. We assume that all blocks that are degree-1,
     * i.e. that contain *only* a reference segment, are scaffold gaps.
     */
    for (int64_t testNum = 0; testNum < 100; testNum++) {
        CactusDisk *cactusDisk = testCommon_getTemporaryCactusDisk(testCase->name);
        eventTree_construct2(cactusDisk);
        Flower *flower = flower_construct(cactusDisk);
        Block *block = block_construct(st_randomInt64(1, 500), flower);
        Event *refEvent = eventTree_getRootEvent(flower_getEventTree(flower));
        // Add in the reference segment.
        segment_construct(block, refEvent);
        stTree *tree = getPhylogeneticTreeRootedAtGivenEvent(refEvent, generateJukesCantorMatrix);

        // Check that we call a string composed only of Ns.
        char *mlString = getMaximumLikelihoodString(tree, block);
        CuAssertIntEquals(testCase, strlen(mlString), block_getLength(block));
        for (int64_t i = 0; i < block_getLength(block); i++) {
            CuAssertTrue(testCase, toupper(mlString[i]) == 'N');
        }

        //Cleanup
        free(mlString);
        cleanupPhylogeneticTree(tree);
        testCommon_deleteTemporaryCactusDisk(testCase->name, cactusDisk);
    }
}

CuSuite* addReferenceCoordinatesTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testMLStringRandom);
    SUITE_ADD_TEST(suite, testMLStringMakesScaffoldGaps);

    return suite;
}
