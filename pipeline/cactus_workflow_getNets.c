#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"

/*
 * This code iterates through the nets returning them recursively.
 */
static void getNets(Net *net, FILE *fileHandle) {
	Net_GroupIterator *groupIterator;
	Group *group;
	assert(net_builtBlocks(net)); //This recursion depends on the block structure having been properly defined for all nodes.
	assert(net != NULL);
	fprintf(fileHandle, "%s %lld\n", netMisc_nameToStringStatic(net_getName(net)), net_getTotalBaseLength(net));
	groupIterator = net_getGroupIterator(net);
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(!group_isTerminal(group)) {
			getNets(group_getNestedNet(group), fileHandle);
		}
	}
	net_destructGroupIterator(groupIterator);
}

int main(int argc, char *argv[]) {
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

	netDisk_destruct(netDisk);

	logInfo("Am finished\n");
	return 0;
}
