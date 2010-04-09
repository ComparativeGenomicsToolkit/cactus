#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"

static void getNets(Net *net, FILE *fileHandle, int32_t includeInternalNodes, int32_t recursive,
		int32_t extendNonZeroTrivialGroups, int32_t minSizeToExtend) {
	Net_GroupIterator *groupIterator;
	Group *group;

	assert(net != NULL);
	if(includeInternalNodes) {
		fprintf(fileHandle, "%s %lld\n", netMisc_nameToStringStatic(net_getName(net)), net_getTotalBaseLength(net));
	}

	if(recursive-- > 0) {
		groupIterator = net_getGroupIterator(net);
		while((group = net_getNextGroup(groupIterator)) != NULL) {
			if(!group_isTerminal(group)) {
				getNets(group_getNestedNet(group), fileHandle, includeInternalNodes,
						recursive, extendNonZeroTrivialGroups, minSizeToExtend);
			}
			else {
				if(extendNonZeroTrivialGroups) {
					int64_t size = group_getTotalBaseLength(group);
					if(size >= minSizeToExtend) {
						group_makeNonTerminal(group);
						fprintf(fileHandle, "%s %lld\n",
								netMisc_nameToStringStatic(group_getName(group)), size);
					}
				}
			}
		}
		net_destructGroupIterator(groupIterator);
	}
}

int main(int argc, char *argv[]) {
	/*
	 * This code iterates through the empty groups and returns them in a list.
	 */
	NetDisk *netDisk;
	Net *net;

	assert(argc == 8);
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

	int32_t extendNonZeroTrivialGroups;
	assert(sscanf(argv[6], "%i", &extendNonZeroTrivialGroups) == 1);

	int32_t minSizeToExtend;
	assert(sscanf(argv[7], "%i", &minSizeToExtend) == 1);

	FILE *fileHandle = fopen(argv[3], "w");
	getNets(net, fileHandle, includeInternalNodes, recursive,
			extendNonZeroTrivialGroups, minSizeToExtend);
	fclose(fileHandle);

	if(extendNonZeroTrivialGroups) {
		netDisk_write(netDisk);
	}
	logInfo("Updated the netdisk\n");

	netDisk_destruct(netDisk);

	logInfo("Am finished\n");
	return 0;
}
