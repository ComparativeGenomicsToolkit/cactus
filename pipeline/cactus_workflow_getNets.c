#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"

int main(int argc, char *argv[]) {
	NetDisk *netDisk;
	Net *net;

	assert(argc == 4);
	netDisk = netDisk_construct(argv[1]);
	logInfo("Set up the net disk\n");

	net = netDisk_getNet(netDisk, netMisc_stringToName(argv[2]));
	assert(net != NULL);
	assert(net_builtBlocks(net)); //This recursion depends on the block structure having been properly defined for all nodes.
	logInfo("Parsed the net\n");

	FILE *fileHandle = fopen(argv[3], "w");
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		if(!group_isTerminal(group)) {
			Net *nestedNet = group_getNestedNet(group);
			assert(nestedNet != NULL);
			assert(net_builtBlocks(nestedNet)); //This recursion depends on the block structure having been properly defined for all nodes.
			fprintf(fileHandle, "%s %lld\n", netMisc_nameToStringStatic(net_getName(nestedNet)), net_getTotalBaseLength(nestedNet));
		}
	}
	net_destructGroupIterator(groupIterator);
	fclose(fileHandle);
	netDisk_destruct(netDisk);
	logInfo("Am finished\n");
	return 0;
}
