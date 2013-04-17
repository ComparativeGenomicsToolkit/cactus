/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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
#include "pairwiseAligner.h"

typedef struct _column Column;
struct _column {
    int64_t seqName;
    int64_t position;
    Column *nColumn;
};

/*
 * Takes a list of DNA strings (char arrays), aligns them and returns a global alignment,
 * represented as a list of aligned stIntTuple pairs, format (score, sequence1, position1, sequence2, position2)
 * Positions and sequence indices are zero based, scores are between 1 and 1000.
 */
stList *makeAlignment(stList *sequences,
        int64_t spanningTrees, int64_t maxPairsToConsider, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters);


stList *makeAlignmentUsingAllPairs(stList *seqs, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters);

#endif /* ALIGNER_H_ */
