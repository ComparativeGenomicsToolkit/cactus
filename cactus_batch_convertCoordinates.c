#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "commonC.h"
#include "cactus_misc.h"

int main(int argc, char *argv[]) {
	/*
	 * For each cigar in file, update the coordinates and write to the second file.
	 */
	assert(argc >= 2);
	FILE *fileHandle = fopen(argv[argc-1], "w");
	int32_t i;
	for(i=1; i<argc-1; i++) {
		uglyf("doing %s %s\n", argv[i], argv[argc-1]);
		convertCoordinates(argv[i], fileHandle);
	}
	fclose(fileHandle);

	//return 1;

	return 0;
}
