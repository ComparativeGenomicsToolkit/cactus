/*
 * aligner.h
 *
 *  Created on: 1 Jul 2010
 *      Author: benedictpaten
 */

#ifndef ALIGNER_H_
#define ALIGNER_H_

//#include "pairwiseAlignmentModel.h"
//#include "pairwiseAlignmentModelInterface.h"
#include "sonLib.h"

/*
 * Takes a list of DNA strings (char arrays), aligns them and returns a global alignment,
 * represented as a list of aligned stIntTuple pairs, format (score, sequence1, position1, sequence2, position2)
 * Positions and sequence indices are zero based, scores are between 1 and 1000.
 */
stList *makeAlignment(stList *sequences,
        int32_t spanningTrees, float gapGamma, bool useBanding, int32_t bandingSize);

#endif /* ALIGNER_H_ */
