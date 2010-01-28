#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"

static void getNets(Net *net, FILE *fileHandle, int32_t includeInternalNodes, int32_t recursive,
		int32_t extendNonZeroTrivialAdjacencyComponents) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator;
	AdjacencyComponent *adjacencyComponent;

	assert(net != NULL);
	if(includeInternalNodes) {
		fprintf(fileHandle, "%s %f\n", netMisc_nameToStringStatic(net_getName(net)), net_getTotalBaseLength(net));
	}

	if(recursive-- > 0) {
		adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
		while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
			if(adjacencyComponent_getNestedNet(adjacencyComponent) != NULL) {
				getNets(adjacencyComponent_getNestedNet(adjacencyComponent), fileHandle, includeInternalNodes,
						recursive, extendNonZeroTrivialAdjacencyComponents);
			}
			else {
				if(extendNonZeroTrivialAdjacencyComponents) {
					double size = adjacencyComponent_getTotalBaseLength(adjacencyComponent);
					if(size > 0 && (size > 400000)) {
						adjacencyComponent_makeNonTerminal(adjacencyComponent);
						fprintf(fileHandle, "%s %f\n",
								netMisc_nameToStringStatic(adjacencyComponent_getName(adjacencyComponent)),	size);
					}
				}
			}
		}
		net_destructAdjacencyComponentIterator(adjacencyComponentIterator);
	}
}

int main(int argc, char *argv[]) {
	/*
	 * This code iterates through the empty adjacency components and returns them in a list.
	 */
	NetDisk *netDisk;
	Net *net;

	assert(argc == 7);
	netDisk = netDisk_construct(argv[1]);
	logInfo("Set up the net disk\n");

	net = netDisk_getNet(netDisk, netMisc_stringToName(argv[2]));
	logInfo("Parsed the net\n");

	int32_t includeInternalNodes;
	assert(sscanf(argv[4], "%i", &includeInternalNodes) == 1);

	int32_t recursive;
	assert(sscanf(argv[5], "%i", &recursive) == 1);

	if(recursive) {
		recursive = INT32_MAX;
	}
	else {
		recursive = 1;
	}

	int32_t extendNonZeroTrivialAdjacencyComponents;
	assert(sscanf(argv[6], "%i", &extendNonZeroTrivialAdjacencyComponents) == 1);

	FILE *fileHandle = fopen(argv[3], "w");
	getNets(net, fileHandle, includeInternalNodes, recursive,
			extendNonZeroTrivialAdjacencyComponents);
	fclose(fileHandle);

	if(extendNonZeroTrivialAdjacencyComponents) {
		netDisk_write(netDisk);
	}
	logInfo("Updated the netdisk\n");

	netDisk_destruct(netDisk);

	logInfo("Am finished\n");
	return 0;
}
