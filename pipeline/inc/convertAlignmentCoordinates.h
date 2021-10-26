/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CONVERT_ALIGNMENT_COORDINATES_H_
#define CONVERT_ALIGNMENT_COORDINATES_H_

#include "cactus.h"

/*
 * Converts input alignments coordinates into coordinates used by cactus.
 */
void convertAlignmentCoordinates(char *inputAlignmentFile, char *outputAlignmentFile, Flower *flower);

/*
 * Strips unique identifiers from sequence IDs (which are added for leaf genomes)
 */
void stripUniqueIdsFromLeafSequences(Flower *flower);

#endif /* CONVERT_ALIGNMENT_COORDINATES_H_ */

