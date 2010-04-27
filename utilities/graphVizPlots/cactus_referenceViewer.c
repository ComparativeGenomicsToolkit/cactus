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
 * The script builds a circos style plot of the reference structure of a net.
 * The format of the output graph is dot format.
 */

/*
 * Global variables.
 */
static bool chainColours = 1;
static bool tangleColours = 1;
static bool nameLabels = 0;

static void usage() {
    fprintf(stderr, "cactus_referenceViewer, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
    fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
    fprintf(stderr, "-e --outputFile : The file to write the dot graph file in.\n");
    fprintf(stderr, "-f --chainColours : Do not give chains distinct colours (instead use just black)\n");
    fprintf(stderr, "-g --tangleColours : Do not give tangles distinct colours (instead use just black)\n");
    fprintf(stderr, "-i --nameLabels : Give chain and net nodes name labels.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}


void addEdge(End *end1, End *end2, const char *colour, const char *label, double length, double weight, const char *direction, FILE *fileHandle) {
    char *nameString1 = netMisc_nameToString(end_getName(end1));
    char *nameString2 = netMisc_nameToString(end_getName(end2));
    graphViz_addEdgeToGraph(nameString1, nameString2, fileHandle, nameLabels ? label : "", colour, length, weight, direction);
    free(nameString1);
    free(nameString2);
}

void addBlockEdge(End *end, FILE *fileHandle, const char *colour) {
	assert(end_isBlockEnd(end));
	assert(end_getOtherBlockEnd(end) != NULL);
	addEdge(end_getOtherBlockEnd(end), end, colour, NULL, 1.5, 1, "forward", fileHandle);
}

void addPseudoAdjacencyEdge(End *_5End, End *_3End, FILE *fileHandle) {
	addEdge(_5End, _3End, "green", NULL, 1.5, 1, "forward", fileHandle);
}

void addAdjacencyEdge(End *_5End, End *_3End, FILE *fileHandle, const char *colour) {
	addEdge(_5End, _3End, colour, NULL, 1.5, 1, "forward", fileHandle);
}

void addEnd(End *end, FILE *fileHandle) {
	const char *nameString = netMisc_nameToStringStatic(end_getName(end));
	const char *colour = "red";
	const char *shape = "circle";
	if(end_isBlockEnd(end)) {
		colour = "black";
		shape = "circle";
	}
	else {
		assert(end_isStubEnd(end));
		assert(end_isAttached(end));
	}
	graphViz_addNodeToGraph(nameString, fileHandle, nameString, 0.5, 0.5, shape, colour, 14);
}

void addPseudoAdjacencyEnds(PseudoAdjacency *pseudoAdjacency, FILE *fileHandle) {
	End *_5End = pseudoAdjacency_get5End(pseudoAdjacency);
	End *_3End = pseudoAdjacency_get3End(pseudoAdjacency);
	addEnd(_5End, fileHandle);
	addEnd(_3End, fileHandle);
	//Now add the all important pseudo-adjacency edge..
	addPseudoAdjacencyEdge(_5End, _3End, fileHandle);
}

void addPseudoChromosomeEnds(PseudoChromosome *pseudoChromosome, FILE *fileHandle) {
	PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator = pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome);
	PseudoAdjacency *pseudoAdjacency;
	while((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(pseudoAdjacencyIterator)) != NULL) {
		addPseudoAdjacencyEnds(pseudoAdjacency, fileHandle);
	}
	pseudoChromosome_destructPseudoAdjacencyIterator(pseudoAdjacencyIterator);
}

void addReferenceEnds(Reference *reference, FILE *fileHandle) {
	Reference_PseudoChromosomeIterator *pseudoChromosomeIterator = reference_getPseudoChromosomeIterator(reference);
	PseudoChromosome *pseudoChromosome;
	while((pseudoChromosome = reference_getNextPseudoChromosome(pseudoChromosomeIterator)) != NULL) {
		addPseudoChromosomeEnds(pseudoChromosome, fileHandle);
	}
	reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
}

void addBlocks(Net *net, FILE *fileHandle) {
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	const char *colour = "black";
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(end_isBlockEnd(end)) {
			End *otherEnd = end_getOtherBlockEnd(end);
			assert(otherEnd != NULL);
			assert(end_getOrientation(end));
			assert(end_getOrientation(otherEnd));
			assert(end_getSide(end) != end_getSide(otherEnd));
			if(end_getSide(end)) {
				addBlockEdge(end, fileHandle, colour);
			}
		}
	}
	net_destructEndIterator(endIterator);
}

void addTangle(Group *group, Hash *endToPseudoAdjacencyHash, FILE *fileHandle) {
	assert(group_isTangle(group));
	Group_EndIterator *endIterator;
	End *end;
	const char *colour = "red";
	while((end = group_getNextEnd(endIterator)) != NULL) {
		if(end_isBlockEnd(end) || end_isAttached(end)) {
			PseudoAdjacency *pseudoAdjacency = hash_search(endToPseudoAdjacencyHash, end);
			assert(pseudoAdjacency != NULL);
			End *otherEnd = pseudoAdjacency_get5End(pseudoAdjacency) == end ? pseudoAdjacency_get3End(pseudoAdjacency) : pseudoAdjacency_get5End(pseudoAdjacency);
			assert(otherEnd != end);
			Cap *cap;
			End_InstanceIterator *capIterator = end_getInstanceIterator(end);
			while((cap = end_getNext(capIterator)) != NULL) {
				Cap *adjacentCap = cap_getAdjacency(cap);
				if(adjacentCap != NULL) {
					End *adjacentEnd = cap_getEnd(adjacentCap);
					if(adjacentEnd != NULL &&
					   (end_isBlockEnd(adjacentEnd) || end_isAttached(adjacentEnd)) &&
					   adjacentEnd != otherEnd) {
						addAdjacencyEdge(end, adjacentEnd, fileHandle, colour);
					}
				}
			}
			end_destructInstanceIterator(capIterator);
		}
	}
	group_destructEndIterator(endIterator);
}

void addTangles(Net *net, Reference *reference, FILE *fileHandle) {
	Hash *endToPseudoAdjacencyHash = reference_getEndToPseudoAdjacencyHash(reference);
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(group_isTangle(group)) { //the group is a tangle (not a link in a chain).
			addTangle(group, endToPseudoAdjacencyHash, fileHandle);
		}
	}
	net_destructGroupIterator(groupIterator);
	hash_destruct(endToPseudoAdjacencyHash);
}

void makeReferenceGraph(Net *net, FILE *fileHandle) {
	Reference *reference = net_getFirstReference(net);
	assert(reference != NULL);

	addReferenceEnds(reference, fileHandle);
	addBlocks(net, fileHandle);
	addTangles(net, reference, fileHandle);
}

int main(int argc, char *argv[]) {
    NetDisk *netDisk;
    Net *net;
    FILE *fileHandle;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * netDiskName = NULL;
    char * netName = NULL;
    char * outputFile = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while(1) {
        static struct option long_options[] = {
            { "logLevel", required_argument, 0, 'a' },
            { "netDisk", required_argument, 0, 'c' },
            { "netName", required_argument, 0, 'd' },
            { "outputFile", required_argument, 0, 'e' },
            { "chainColours", no_argument, 0, 'f' },
            { "tangleColours", no_argument, 0, 'g' },
            { "nameLabels", no_argument, 0, 'h' },
            { "help", no_argument, 0, 'i' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:e:fghi", long_options, &option_index);

        if(key == -1) {
            break;
        }

        switch(key) {
            case 'a':
                logLevelString = stringCopy(optarg);
                break;
            case 'c':
                netDiskName = stringCopy(optarg);
                break;
            case 'd':
                netName = stringCopy(optarg);
                break;
            case 'e':
                outputFile = stringCopy(optarg);
                break;
            case 'f':
                chainColours = !chainColours;
                break;
            case 'g':
                tangleColours = !tangleColours;
                break;
            case 'i':
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

    assert(netDiskName != NULL);
    assert(netName != NULL);
    assert(outputFile != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        setLogLevel(LOGGING_INFO);
    }
    if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        setLogLevel(LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    logInfo("Net disk name : %s\n", netDiskName);
    logInfo("Net name : %s\n", netName);
    logInfo("Output graph file : %s\n", outputFile);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    netDisk = netDisk_construct(netDiskName);
    logInfo("Set up the net disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
    logInfo("Parsed the top level net of the cactus tree to build\n");

    ///////////////////////////////////////////////////////////////////////////
    // Build the graph.
    ///////////////////////////////////////////////////////////////////////////

    fileHandle = fopen(outputFile, "w");
    graphViz_setupGraphFile(fileHandle);
    makeReferenceGraph(net, fileHandle);
    graphViz_finishGraphFile(fileHandle);
    fclose(fileHandle);
    logInfo("Written the reference graph to file\n");

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    netDisk_destruct(netDisk);

    return 0;
}
