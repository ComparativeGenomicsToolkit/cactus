#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>

#include "cactus.h"

int main(int argc, char *argv[]) {
	NetDisk *netDisk;
	Net *net;

	assert(argc >= 3);
	netDisk = netDisk_construct(argv[1]);
	st_logInfo("Set up the net disk\n");

	FILE *fileHandle = fopen(argv[2], "w");

	int32_t i;
	for(i=3; i<argc; i++) {
		net = netDisk_getNet(netDisk, netMisc_stringToName(argv[i]));
		assert(net != NULL);
		assert(net_builtBlocks(net)); //This recursion depends on the block structure having been properly defined for all nodes.
		st_logInfo("Parsed the net %s\n", argv[i]);
		Net_GroupIterator *groupIterator = net_getGroupIterator(net);
		Group *group;
		while((group = net_getNextGroup(groupIterator)) != NULL) {
			if(!group_isLeaf(group)) {
				Net *nestedNet = group_getNestedNet(group);
				assert(nestedNet != NULL);
				assert(net_builtBlocks(nestedNet)); //This recursion depends on the block structure having been properly defined for all nodes.
				fprintf(fileHandle, "%s %" PRIi64 " \n", netMisc_nameToStringStatic(net_getName(nestedNet)), net_getTotalBaseLength(nestedNet));
			}
		}
		net_destructGroupIterator(groupIterator);
	}
	fclose(fileHandle);
	netDisk_destruct(netDisk);
	st_logInfo("Am finished\n");
	return 0;
}
