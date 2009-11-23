#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"

static void getNets(Net *net, FILE *fileHandle) {
	Net_AdjacencyComponentIterator *adjacencyComponentIterator;
	AdjacencyComponent *adjacencyComponent;

	assert(net != NULL);
	if(net_getAdjacencyComponentNumber(net) == 0) {
		assert(net_getAtomNumber(net) == 0);
		fprintf(fileHandle, "%s %f\n", netMisc_nameToStringStatic(net_getName(net)), net_getTotalBaseLength(net));
	}

	adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		getNets(adjacencyComponent_getNestedNet(adjacencyComponent), fileHandle);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);
}

int main(int argc, char *argv[]) {
	/*
	 * This code iterates through the empty adjacency components and returns them in a list.
	 */
	NetDisk *netDisk;
	Net *net;

	assert(argc == 4);
	netDisk = netDisk_construct(argv[1]);
	logInfo("Set up the net disk\n");

	net = netDisk_getNet(netDisk, netMisc_stringToName(argv[2]));
	logInfo("Parsed the net\n");

	FILE *fileHandle = fopen(argv[3], "w");
	getNets(net, fileHandle);
	fclose(fileHandle);

	logInfo("Am finished\n");
	return 0;
}
