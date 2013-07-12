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

typedef struct _seqFrag SeqFrag;
struct _seqFrag {
    char *seq;
    int64_t length;
    bool missingLeftEnd;
    bool missingRightEnd;
};

typedef struct _column Column;
struct _column {
    int64_t seqName;
    int64_t position;
    Column *nColumn;
};

/*
 * Takes a list of DNA strings (as SeqFrags), aligns them and returns a global alignment,
 * represented as a list of aligned stIntTuple pairs, format (score, sequence1, position1, sequence2, position2)
 * Positions and sequence indices are zero based, scores are between 1 and 1000.
 */
stList *makeAlignment(stList *seqFrags,
        int64_t spanningTrees, int64_t maxPairsToConsider, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters);


stList *makeAlignmentUsingAllPairs(stList *seqFrags, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters);

SeqFrag *seqFrag_construct(const char *seq, bool missingLeftEnd, bool missingRightEnd);

void seqFrag_destruct(SeqFrag *seqFrag);

#endif /* ALIGNER_H_ */
