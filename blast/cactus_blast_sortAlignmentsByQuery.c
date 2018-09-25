/*
 * Copyright (C) 2009-2014 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "stLastzAlignments.h"

int main(int argc, char *argv[]) {
	/*
	 * Sort cigar file in descending order of query start coordinate.
	 */
	assert(argc == 4);
	st_setLogLevelFromString(argv[1]);
	stCaf_sortCigarsFileByFirstSequenceStartCoordinateInAscendingOrder(argv[2], argv[3]);
	return 0;
}
