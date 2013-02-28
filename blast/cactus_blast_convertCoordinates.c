/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "commonC.h"
#include "hashTableC.h"
#include "bioioC.h"
#include "cactus.h"
#include "avl.h"
#include "pairwiseAlignment.h"
#include "blastAlignmentLib.h"

int main(int argc, char *argv[]) {
	/*
	 * For each cigar in file, update the coordinates and write to the second file.
	 */
	assert(argc == 3);
	FILE *fileHandleIn = fopen(argv[1], "r");
	FILE *fileHandleOut = fopen(argv[2], "w");
	int size = 100;
	char *cA = st_calloc(size+1, sizeof(char));
	int32_t i;
	do {
	    i = benLine(&cA, &size, fileHandleIn);
	    if(strlen(cA) > 0) {
	        convertCoordinates(cA, fileHandleOut);
	    }
	} while(i != -1);
	fclose(fileHandleIn);
	fclose(fileHandleOut);
	free(cA);

	//return 1;

	return 0;
}
