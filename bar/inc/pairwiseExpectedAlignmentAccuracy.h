/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * pairwiseExpectedAlignmentAccuracy.h
 *
 *  Created on: 14 Dec 2010
 *      Author: benedictpaten
 */

#ifndef PAIRWISEEXPECTEDALIGNMENTACCURACY_H_
#define PAIRWISEEXPECTEDALIGNMENTACCURACY_H_

stHash *getExpectedAlignmentAccuracyScores(stList *alignedPairs, int64_t *indelProbs1, int64_t *indelProbs2, int64_t seqLength1, int64_t seqLength2);

#endif /* PAIRWISEEXPECTEDALIGNMENTACCURACY_H_ */
