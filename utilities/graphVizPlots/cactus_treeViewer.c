/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * The script builds a cactus tree representation of the chains and flowers.
 * The format of the output graph is dot format. Only works with normalised
 * cactus databases.
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
static double totalProblemSize;
static bool scaleNodeSizes = 1;
static bool nameLabels = 0;

static void usage() {
    fprintf(stderr, "cactus_treeViewer, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr,
            "-e --outputFile : The file to write the dot graph file in.\n");
    fprintf(
            stderr,
            "-f --scaleNodeSizes : Scale the node sizes according to the volume of sequence they contained.\n");
    fprintf(stderr, "-g --nameLabels : Give chain and flower nodes name labels.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

static void addNodeToGraph(const char *nodeName, FILE *graphFileHandle,
        double scalingFactor, const char *shape, const char *nodeLabel) {
    /*
     * Adds a node to the graph, scaling it's size.
     */

    double height = 1;
    double width = 0.25;
    scalingFactor *= 100.0;
    if (scaleNodeSizes && scalingFactor >= 1) {
        height = 4 * sqrt(scalingFactor);
        width = 1.0 * sqrt(scalingFactor);
    }
    graphViz_addNodeToGraph(nodeName, graphFileHandle, nameLabels ? nodeLabel
            : "", width, height, shape, "black", 14);
}

void makeCactusTree_terminalNode(Flower *flower, FILE *fileHandle,
        const char *parentNodeName, const char *parentEdgeColour) {
    char *groupNameString = cactusMisc_nameToString(flower_getName(flower));
    double scalingFactor = flower_getTotalBaseLength(flower) / totalProblemSize;
    assert(scalingFactor <= 1.001);
    assert(scalingFactor >= -0.001);
    addNodeToGraph(groupNameString, fileHandle, scalingFactor, "triangle",
            groupNameString);
    //Write in the parent edge.
    if (parentNodeName != NULL) {
        graphViz_addEdgeToGraph(parentNodeName, groupNameString, fileHandle,
                "", parentEdgeColour, 10, 1, "forward");
    }
    free(groupNameString);
}

void makeCactusTree_flower(Flower *flower, FILE *fileHandle, const char *parentNodeName,
        const char *parentEdgeColour);

void makeCactusTree_chain(Chain *chain, FILE *fileHandle,
        const char *parentNodeName, const char *parentEdgeColour) {
    //Write the flower nodes.
    char *chainNameString = cactusMisc_nameToString(chain_getName(chain));
    const char *edgeColour = graphViz_getColour();
    addNodeToGraph(chainNameString, fileHandle,
            chain_getAverageInstanceBaseLength(chain) / totalProblemSize,
            "box", chainNameString);
    //Write in the parent edge.
    if (parentNodeName != NULL) {
        graphViz_addEdgeToGraph(parentNodeName, chainNameString, fileHandle,
                "", parentEdgeColour, 10, 1, "forward");
    }
    //Create the linkers to the nested flowers.
    Link *link = chain_getFirst(chain);
    while(link != NULL) {
        Group *group = link_getGroup(link);
        assert(group != NULL);
        assert(!group_isLeaf(group));
        if (!group_isLeaf(group)) {
            makeCactusTree_flower(group_getNestedFlower(group), fileHandle,
                    chainNameString, edgeColour);
        }
        link = link_getNextLink(link);
    }
    free(chainNameString);
}

void makeCactusTree_flower(Flower *flower, FILE *fileHandle, const char *parentNodeName,
        const char *parentEdgeColour) {
    if(flower_isTerminal(flower)) {
        makeCactusTree_terminalNode(flower, fileHandle, parentNodeName, parentEdgeColour);
    }
    else {
        //Write the flower nodes.
        char *flowerNameString = cactusMisc_nameToString(flower_getName(flower));
        const char *edgeColour = graphViz_getColour();
        addNodeToGraph(flowerNameString, fileHandle, flower_getTotalBaseLength(flower)
                / totalProblemSize, "ellipse", flowerNameString);
        //Write in the parent edge.
        if (parentNodeName != NULL) {
            graphViz_addEdgeToGraph(parentNodeName, flowerNameString, fileHandle, "",
                    parentEdgeColour, 10, 1, "forward");
        }
        //Create the chains.
        Flower_ChainIterator *chainIterator = flower_getChainIterator(flower);
        Chain *chain;
        while ((chain = flower_getNextChain(chainIterator)) != NULL) {
            makeCactusTree_chain(chain, fileHandle, flowerNameString, edgeColour);
        }
        flower_destructChainIterator(chainIterator);

        //Create the diamond node
        char *diamondNodeNameString = st_malloc(sizeof(char) * (strlen(
                flowerNameString) + 2));
        sprintf(diamondNodeNameString, "z%s", flowerNameString);
        const char *diamondEdgeColour = graphViz_getColour();
        //Create all the groups linked to the diamond.
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        double size = 0.0; //get the size of the group organising node..
        int32_t nonTrivialGroupCount = 0;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            assert(!group_isLeaf(group));
            if (group_isTangle(group)) {
                size += group_getTotalBaseLength(group);
                nonTrivialGroupCount++;
            }
        }
        flower_destructGroupIterator(groupIterator);
        if(nonTrivialGroupCount == 0) {
            assert(flower_getParentGroup(flower) == 0);
        }
        else {
            //assert(nonTrivialGroupCount > 0);
            addNodeToGraph(diamondNodeNameString, fileHandle, size
                    / totalProblemSize, "diamond", "");
            graphViz_addEdgeToGraph(flowerNameString, diamondNodeNameString,
                    fileHandle, "", edgeColour, 10, 1, "forward");
            groupIterator = flower_getGroupIterator(flower);
            while ((group = flower_getNextGroup(groupIterator)) != NULL) {
                if (group_isTangle(group)) {
                    assert(!group_isLeaf(group));
                    makeCactusTree_flower(group_getNestedFlower(group), fileHandle,
                                                diamondNodeNameString, diamondEdgeColour);
                }
            }
            flower_destructGroupIterator(groupIterator);
        }

        free(flowerNameString);
        free(diamondNodeNameString);
    }
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
                "outputFile", required_argument, 0, 'e' }, { "scaleNodeSizes",
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
                scaleNodeSizes = !scaleNodeSizes;
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
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower name : %s\n", flowerName);
    st_logInfo("Output graph file : %s\n", outputFile);

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

    totalProblemSize = flower_getTotalBaseLength(flower);
    fileHandle = fopen(outputFile, "w");
    graphViz_setupGraphFile(fileHandle);
    makeCactusTree_flower(flower, fileHandle, NULL, NULL);
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
