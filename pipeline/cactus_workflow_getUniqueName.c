#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "cactus.h"


int main(int argc, char *argv[]) {
	/*
	 * This code simply creates a unique name and puts it in a file.
	 */
	assert(argc == 3);
	NetDisk *netDisk = netDisk_construct(argv[1]);
	Name name = netDisk_getUniqueID(netDisk);
	FILE *fileHandle = fopen(argv[2], "w");
	fprintf(fileHandle, "%s\n", netMisc_nameToStringStatic(name));
	fclose(fileHandle);
	return 0;
}
