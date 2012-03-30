/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * adjacencyPairsHash.c
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#include <stdlib.h>
#include "commonC.h"
#include "sonLib.h"
#include "shared.h"

const char *MATCHING_EXCEPTION = "MATCHING_EXCEPTION";

/*
 * Functions to assess matching
 */

int32_t matchingCardinality(stList *matching) {
    /*
     * Returns number of edges with weight > 0.
     */
    int32_t totalCardinality = 0;
    for (int32_t i = 0; i < stList_length(matching); i++) {
        stIntTuple *edge = stList_get(matching, i);
        totalCardinality += stIntTuple_getPosition(edge, 2) > 0 ? 1 : 0;
    }
    return totalCardinality;
}

int32_t matchingWeight(stList *matching) {
    /*
     * Returns sum of weights.
     */
    int32_t totalWeight = 0;
    for(int32_t i=0; i<stList_length(matching); i++) {
        stIntTuple *edge = stList_get(matching, i);
        totalWeight += stIntTuple_getPosition(edge, 2);
    }
    return totalWeight;
}

/*
 * Code to talk to the blossom5 maximum weight perfect matching algorithm.
 */

static void writeGraph(FILE *fileHandle, stList *edges, int32_t nodeNumber) {
    /*
     * Writes out just the adjacencies in the blossom format.
     */
    int32_t edgeNumber = stList_length(edges);
    fprintf(fileHandle, "%i %i\n", nodeNumber, edgeNumber);
    for(int32_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        int32_t from =  stIntTuple_getPosition(edge, 0);
        int32_t to = stIntTuple_getPosition(edge, 1);
        int32_t weight = stIntTuple_getPosition(edge, 2);
        //All the algorithms are minimisation algorithms, so we invert the sign.
        fprintf(fileHandle, "%i %i %i\n", from, to, weight);
    }
}

static void writeCliqueGraph(FILE *fileHandle, stList *edges, int32_t nodeNumber, bool negativeWeights) {
    /*
     * Writes out a representation of the adjacencies and ends as a graph readable by blossom.
     * Writes out additional edges so that every pair of nodes is connected.
     */
    int32_t edgeNumber = ((nodeNumber * nodeNumber)  - nodeNumber) / 2;
    fprintf(fileHandle, "%i %i\n", nodeNumber, edgeNumber);
    stSortedSet *seen = stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct);
    int32_t edgesWritten = 0;
    for(int32_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        int32_t from =  stIntTuple_getPosition(edge, 0);
        int32_t to = stIntTuple_getPosition(edge, 1);
#ifdef BEN_DEBUG
        assert(from < nodeNumber);
        assert(to < nodeNumber);
        assert(from >= 0);
        assert(to >= 0);
        assert(from != to);
#endif
        int32_t weight = stIntTuple_getPosition(edge, 2);
        //If is a minimisation algorithms we invert the sign..
        fprintf(fileHandle, "%i %i %i\n", from, to, negativeWeights ? -weight : weight);
        edgesWritten++;
        addEdgeToSet(seen, from, to);
    }
    for(int32_t i=0; i<nodeNumber; i++) {
        for(int32_t j=i+1; j<nodeNumber; j++) {
            if(!edgeInSet(seen, i, j)) {
                fprintf(fileHandle, "%i %i 0\n", i, j);
                edgesWritten++;
            }
        }
    }
    //Cleanup
    stSortedSet_destruct(seen);
    assert(edgeNumber == edgesWritten);
}

static stHash *putEdgesInHash(stList *edges) {
    stHash *intsToEdgesHash = stHash_construct3((uint32_t (*)(const void *))stIntTuple_hashKey, (int (*)(const void *, const void *))stIntTuple_equalsFn, (void (*)(void *))stIntTuple_destruct, NULL);
    for(int32_t i=0; i<stList_length(edges); i++) {
        stIntTuple *edge = stList_get(edges, i);
        stHash_insert(intsToEdgesHash, constructEdge(stIntTuple_getPosition(edge, 0), stIntTuple_getPosition(edge, 1)), edge);
    }
    return intsToEdgesHash;
}

static stList *readMatching(FILE *fileHandle, stList *originalEdges) {
    /*
     * Reads the matching created by Blossum.
     */
    stHash *originalEdgesHash = putEdgesInHash(originalEdges);
    char *line = stFile_getLineFromFile(fileHandle);
    assert(line != NULL);
    int32_t nodeNumber, edgeNumber;
    int32_t i = sscanf(line, "%i %i\n", &nodeNumber, &edgeNumber);
    assert(i == 2);
    free(line);
    stList *chosenEdges = stList_construct();
    for(int32_t j=0; j<edgeNumber; j++) {
        line = stFile_getLineFromFile(fileHandle);
        int32_t node1, node2;
        i = sscanf(line, "%i %i", &node1, &node2);
        assert(i == 2);
        free(line);
#ifdef BEN_DEBUG
        assert(node1 >= 0);
        assert(node1 < nodeNumber);
        assert(node2 >= 0);
        assert(node2 < nodeNumber);
#endif
        stIntTuple *edge = constructEdge(node1, node2);
        stIntTuple *originalEdge = stHash_search(originalEdgesHash, edge);
        if(originalEdge != NULL) {
            stList_append(chosenEdges, originalEdge);
        }
        stIntTuple_destruct(edge);
    }
    stHash_destruct(originalEdgesHash);
    return chosenEdges;
}

static stList *chooseAdjacencyPairing_externalProgram(stList *edges, int32_t nodeNumber, const char *programName) {
    /*
     * We create temp files to hold stuff.
     */
    if(nodeNumber <= 1) {
        assert(stList_length(edges) == 0);
        return stList_construct();
    }

    char *tempInputFile = getTempFile(), *tempOutputFile = getTempFile();

    /*
     * We write the graph to a temp file.
     */
    FILE *fileHandle = fopen(tempInputFile, "w");
    if(strcmp(programName, "blossom5") == 0) { //Must be all connected as
        //generates perfect matchings.
        writeCliqueGraph(fileHandle, edges, nodeNumber, 1);
    }
    else {
        writeGraph(fileHandle, edges, nodeNumber);
    }
    fclose(fileHandle);

    /*
     * We run the external program.
     */
    char *command = stString_print("%s -e %s -w %s >& /dev/null", programName, tempInputFile, tempOutputFile);
    int32_t i = st_system(command);
    if(i != 0) {
        stThrowNew(MATCHING_EXCEPTION, "Something went wrong the command %s\n", command);
    }
    free(command);

    /*
     * We get back the matching.
     */
    fileHandle = fopen(tempOutputFile, "r");
    stList *matching = readMatching(fileHandle, edges);
    fclose(fileHandle);
    st_logDebug("The adjacency matching for %i nodes with %i initial edges contains %i edges\n", nodeNumber, stList_length(edges), stList_length(matching));

    /*
     * Get rid of the temp files..
     */
    st_system("rm -rf %s %s", tempInputFile, tempOutputFile);
    free(tempInputFile);
    free(tempOutputFile);

    return matching;
}

stList *chooseMatching_blossom5(stList *edges, int32_t nodeNumber) {
    return chooseAdjacencyPairing_externalProgram(edges, nodeNumber, "blossom5");
}

stList *chooseMatching_maximumCardinalityMatching(stList *edges, int32_t nodeNumber) {
    return chooseAdjacencyPairing_externalProgram(edges, nodeNumber, "matchGraph.py -c");
}

stList *chooseMatching_maximumWeightMatching(stList *edges, int32_t nodeNumber) {
    return chooseAdjacencyPairing_externalProgram(edges, nodeNumber, "matchGraph.py");
}

/*
 * Greedy matching algorithm.
 */

int chooseMatching_greedyP(const void *a, const void *b) {
    return stIntTuple_getPosition((stIntTuple *)a, 2) - stIntTuple_getPosition((stIntTuple *)b, 2);
}

stList *chooseMatching_greedy(stList *edges, int32_t nodeNumber) {
    /*
     * Greedily picks the edge from the list such that each node has at most one edge.
     */
    //First clone the list..
    edges = stList_copy(edges, NULL);

    stSortedSet *seen = getEmptyNodeOrEdgeSetWithCleanup();
    stList *matching = stList_construct();

    //Sort the adjacency pairs..
    stList_sort(edges, chooseMatching_greedyP);

#ifdef BEN_DEBUG
    double strength = INT32_MAX;
#endif
    while (stList_length(edges) > 0) {
        stIntTuple *edge = stList_pop(edges);
#ifdef BEN_DEBUG
        double d = stIntTuple_getPosition(edge, 2);
        assert(d <= strength);
        strength = d;
#endif
        if(!nodeInSet(seen, stIntTuple_getPosition(edge, 0)) && !nodeInSet(seen, stIntTuple_getPosition(edge, 1))) {
            addNodeToSet(seen, stIntTuple_getPosition(edge, 0));
            addNodeToSet(seen, stIntTuple_getPosition(edge, 1));
            stList_append(matching,edge);
        }
    }
    assert(stList_length(edges) == 0);
    stList_destruct(edges);
    stSortedSet_destruct(seen);

    return matching;
}
