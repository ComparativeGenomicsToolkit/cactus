/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * The script builds a cactus tree representation of the chains and flowers.
 * The format of the output graph is dot format.
 */
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"

/*
 * Global variables.
 */
static bool edgeColours = 1;
static bool nameLabels = 0;

static void usage() {
    fprintf(stderr, "cactus_graphViewer, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr,
            "-e --outputFile : The file to write the dot graph file in.\n");
    fprintf(
            stderr,
            "-f --chainColours : Do not give chains distinct colours (instead of just black)\n");
    fprintf(stderr, "-g --nameLabels : Give chain and flower nodes name labels.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

void addEndNodeToGraph(End *end, FILE *fileHandle) {
    const char *nameString = cactusMisc_nameToStringStatic(end_getName(end));
    graphViz_addNodeToGraph(nameString, fileHandle, nameString, 0.5, 0.5,
            "circle", "black", 14);
}

void addEdgeToGraph(End *end1, End *end2, const char *colour,
        const char *label, double length, double weight, const char *direction,
        FILE *fileHandle) {
    char *nameString1 = cactusMisc_nameToString(end_getName(end1));
    char *nameString2 = cactusMisc_nameToString(end_getName(end2));
    graphViz_addEdgeToGraph(nameString1, nameString2, fileHandle,
            nameLabels ? label : "", colour, length, weight, direction);
    free(nameString1);
    free(nameString2);
}

void addBlockToGraph(Block *block, const char *colour, FILE *fileHandle) {
    static char label[100000];
    End *leftEnd = block_get5End(block);
    End *rightEnd = block_get3End(block);
    addEndNodeToGraph(leftEnd, fileHandle);
    addEndNodeToGraph(rightEnd, fileHandle);
    Block_InstanceIterator *iterator = block_getInstanceIterator(block);
    Segment *segment;
    while ((segment = block_getNext(iterator)) != NULL) {
        segment = segment_getStrand(segment) ? segment : segment_getReverse(
                segment);
        if (segment_getSequence(segment) != NULL) {
            sprintf(label, "%s:%i:%i", cactusMisc_nameToStringStatic(
                    sequence_getName(segment_getSequence(segment))),
                    segment_getStart(segment), segment_getStart(segment)
                            + segment_getLength(segment));
            addEdgeToGraph(cap_getEnd(segment_get5Cap(segment)), cap_getEnd(
                    segment_get3Cap(segment)), edgeColours ? colour : "black",
                    label, 1.5, 100, "forward", fileHandle);
        }
    }
    block_destructInstanceIterator(iterator);
}

void addTrivialChainsToGraph(Flower *flower, FILE *fileHandle) {
    /*
     * Add blocks not part of chain to the graph
     */
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIterator)) != NULL) {
        if (block_getChain(block) == NULL) {
            addBlockToGraph(block, "black", fileHandle);
        }
    }
    flower_destructBlockIterator(blockIterator);
}

void addChainsToGraph(Flower *flower, FILE *fileHandle) {
    /*
     * Add blocks part of a chain to the graph.
     */
    Flower_ChainIterator *chainIterator = flower_getChainIterator(flower);
    Chain *chain;
    while ((chain = flower_getNextChain(chainIterator)) != NULL) {
        int32_t i, j;
        const char *chainColour;
        while ((chainColour = graphViz_getColour()) != NULL) { //ensure the chain colours don't match the trivial block chains and the adjacencies.
            if (strcmp(chainColour, "black") != 0
                    && strcmp(chainColour, "grey") != 0) {
                break;
            }
        }
        Block **blocks = chain_getBlockChain(chain, &i);
        for (j = 0; j < i; j++) {
            addBlockToGraph(blocks[j], chainColour, fileHandle);
        }
        free(blocks);
    }
    flower_destructChainIterator(chainIterator);
}

void addAdjacencies(Group *group, FILE *fileHandle) {
    /*
     * Adds adjacency edges to the graph.
     */
    static char label[10000];
    Group_EndIterator *endIterator = group_getEndIterator(group);
    End *end;
    Flower *flower = group_getFlower(group);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
        Cap *cap;
        char *flowerName = cactusMisc_nameToString(flower_getName(flower));
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            if (cap_getSequence(cap) != NULL) {
                cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                Cap *cap2 = cap_getAdjacency(cap);
                if (!cap_getSide(cap)) {
                    assert(cap_getCoordinate(cap) < cap_getCoordinate(cap2));
                    sprintf(label, "%s:%i:%i:%s:%i",
                            cactusMisc_nameToStringStatic(sequence_getName(
                                    cap_getSequence(cap))), cap_getCoordinate(
                                    cap), cap_getCoordinate(cap2), flowerName,
                            flower_getEndNumber(flower));
                    //sprintf(label, "%s:%i",
                    //		flowerName,
                    //		flower_getEndNumber(flower));
                    addEdgeToGraph(cap_getEnd(cap), cap_getEnd(cap2), "grey",
                            label, 1.5, 1, "forward", fileHandle);
                }
            }
        }
        free(flowerName);
        end_destructInstanceIterator(instanceIterator);
    }
    group_destructEndIterator(endIterator);
}

void addStubAndCapEndsToGraph(Flower *flower, FILE *fileHandle) {
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        if (end_isStubEnd(end)) {
            addEndNodeToGraph(end, fileHandle);
        }
    }
    flower_destructEndIterator(endIterator);
}

void makeCactusGraph(Flower *flower, FILE *fileHandle) {
    if (flower_getParentGroup(flower) == NULL) {
        addStubAndCapEndsToGraph(flower, fileHandle);
    }
    addTrivialChainsToGraph(flower, fileHandle);
    addChainsToGraph(flower, fileHandle);
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        Flower *nestedFlower = group_getNestedFlower(group);
        if (nestedFlower != NULL) {
            makeCactusGraph(nestedFlower, fileHandle);
        } else { //time to add the adjacencies!
            addAdjacencies(group, fileHandle);
        }
    }
    flower_destructGroupIterator(groupIterator);
}

int main(int argc, char *argv[]) {
    Flower *flower;
    FILE *fileHandle;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;
    char * outputFile = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument, 0,
                'c' }, { "flowerName", required_argument, 0, 'd' }, {
                "outputFile", required_argument, 0, 'e' }, { "edgeColours",
                no_argument, 0, 'f' }, { "nameLabels", no_argument, 0, 'g' }, {
                "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:e:fgh", long_options,
                &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'd':
                flowerName = stString_copy(optarg);
                break;
            case 'e':
                outputFile = stString_copy(optarg);
                break;
            case 'f':
                edgeColours = !edgeColours;
                break;
            case 'g':
                nameLabels = !nameLabels;
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(cactusDiskDatabaseString != NULL);
    assert(flowerName != NULL);
    assert(outputFile != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if (logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(flowerName));
    st_logInfo("Parsed the top level flower of the cactus tree to build\n");

    ///////////////////////////////////////////////////////////////////////////
    // Build the graph.
    ///////////////////////////////////////////////////////////////////////////

    fileHandle = fopen(outputFile, "w");
    graphViz_setupGraphFile(fileHandle);
    makeCactusGraph(flower, fileHandle);
    graphViz_finishGraphFile(fileHandle);
    fclose(fileHandle);
    st_logInfo("Written the tree to file\n");

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
