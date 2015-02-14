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
    int64_t leftEndId;
    int64_t rightEndId;
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
typedef struct {
    stSet *columns; //set of "Columns" representing alignment
    stList *alignedPairs; //set of tuples giving aligned pairs consistent with columns each of the form (similarityScore, seqX, posX, seqY, posY), positions and sequence indices are zero based.
    stList *chosenPairwiseAlignments; //set of tuples representing chosen pairwise alignments, each of the form (similarityScore, seqXIndex, seqYIndex)
} MultipleAlignment;

MultipleAlignment *makeAlignment(StateMachine *sM, stList *seqFrags,
        int64_t spanningTrees, int64_t maxPairsToConsider,
        bool useProgressiveMerging,
        float matchGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters);

MultipleAlignment *makeAlignmentUsingAllPairs(StateMachine *sM, stList *seqFrags,
        bool useProgressiveMerging, float matchGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters);

void multipleAlignment_destruct(MultipleAlignment *mA);

SeqFrag *seqFrag_construct(const char *seq, int64_t leftEndId, int64_t rightEndId);

void seqFrag_destruct(SeqFrag *seqFrag);

/*
 * This is a pairwise expected accuracy alignment function that uses the multiple alignment code, kind of odd.
 */
stList *filterPairwiseAlignmentToMakePairsOrdered(stList *alignedPairs, const char *seqX, const char *seqY, float matchGamma);

/*
 * Declarations for functions tested by unit-tests, but probably not really useful for stuff outside of this module.
 */

stSet *makeColumns(stList *seqFrags);

stList *makeAllPairwiseAlignments(StateMachine *sM, stList *seqs, PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters, stList **seqPairSimilarityScores);

stHash *makeAlignmentWeightAdjacencyLists(stSet *columns, stList *multipleAlignedPairs);

stSortedSet *makeOrderedSetOfAlignmentWeights(stHash *alignmentWeightAdjLists);

int64_t *getDistanceMatrix(stSet *columns, stList *seqs, int64_t maxPairsToConsider);

double subsPerSite(int64_t seq1, int64_t seq2, int64_t *distanceCounts, int64_t seqNo);

stSet *getMultipleSequenceAlignment(stList *seqs, stList *multipleAlignedPairs, double matchGamma);

stList *makeColumnSequences(stList *seqFrags, stSet *columns);

stSet *getMultipleSequenceAlignmentProgressive(stList *seqFrags, stList *multipleAlignedPairs, double matchGamma, stList *seqPairSimilarityScores);

stList *pairwiseAlignColumns(stList *seqXColumns, stList *seqYColumns, stHash *alignmentWeightAdjLists, stSet *columns,
        stSortedSet *alignmentWeightsOrderedByWeight, double matchGamma);

stList *filterMultipleAlignedPairs(stSet *columns, stList *multipleAlignedPairs);

#endif /* ALIGNER_H_ */
