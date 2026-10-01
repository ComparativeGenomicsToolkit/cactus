/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * endAligner.h
 *
 *  Created on: 1 Jul 2010
 *      Author: benedictpaten
 */

#ifndef ENDALIGNER_H_
#define ENDALIGNER_H_

#include "sonLib.h"
#include "cactus.h"
#include "pairwiseAligner.h"

typedef struct _AlignedPair {
    int64_t subsequenceIdentifier;
    int64_t position;
    bool strand;
    int64_t score;
    struct _AlignedPair *reverse;
} AlignedPair;

/*
 * Constructs the an aligned pair.
 */
AlignedPair *alignedPair_construct(int64_t subsequenceIdentifier1, int64_t position1, bool strand1,
        int64_t subsequenceIdentifier2, int64_t position2, bool strand2, int64_t score1, int64_t score2);

/*
 * Destruct the aligned pair.
 */
void alignedPair_destruct(AlignedPair *alignedPair);

/*
 * Compares two aligned pairs.
 */
int alignedPair_cmpFn(const AlignedPair *alignedPair1, const AlignedPair *alignedPair2);

/*
 * Creates a global alignment (as a set of aligned pairs) of the sequences from the end,
 * the pairs returned are ordered according
 * to the alignerPair comparison function.
 */
stSortedSet *makeEndAlignment(StateMachine *sM, End *end, int64_t spanningTrees, int64_t maxSequenceLength,
                              bool useProgressiveMerging, float gapGamma,
                              PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters);

/*
 * Path lengths between genomes, for choosing scoring models by how diverged two genomes are: on
 * <reference reconstructionTree> (the unscaled tree, as ancestral reconstruction uses it) when it names
 * both genomes, else on the event tree, whose branches are the ones scaled for alignment sensitivity.
 */
typedef struct _barDistances BarDistances;

BarDistances *barDistances_constructFromCactusParams(CactusParams *params);

double barDistances_get(BarDistances *distances, Event *event1, Event *event2);

/*
 * Whether the distances are on the reconstruction tree (else the event tree).
 */
bool barDistances_hasReconstructionTree(BarDistances *distances);

void barDistances_destruct(BarDistances *distances);

/*
 * Per-pair pecan models.  <pecan pecanPairModels="N"> switches them on, and its children <pairModel0> ..
 * <pairModelN-1> each carry a maxDistance and an hmm (a five-state pair HMM in the JSON hmm_jsonParse reads),
 * in ascending maxDistance.  Each pair of sequences is aligned with the first model whose maxDistance is at
 * least the distance between their genomes (see BarDistances).  Pairs from one genome (paralogs, whose
 * divergence the tree does not give) and pairs beyond the last model keep the default machine.
 */
typedef struct _pecanPairModels PecanPairModels;

/*
 * NULL when the config has no pair models.
 */
PecanPairModels *pecanPairModels_constructFromCactusParams(CactusParams *params);

void pecanPairModels_destruct(PecanPairModels *pairModels);

/*
 * As makeEndAlignment, aligning each pair of sequences with its model from pairModels (may be NULL).
 */
stSortedSet *makeEndAlignmentWithPairModels(StateMachine *sM, PecanPairModels *pairModels, End *end, int64_t spanningTrees,
                                            int64_t maxSequenceLength, bool useProgressiveMerging, float gapGamma,
                                            PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters);

/*
 * Writes an end alignment to the given file.
 */
void writeEndAlignmentToDisk(End *end, stSortedSet *endAlignment, FILE *fileHandle);

/*
 * Loads an end alignment from the given file.
 */
stSortedSet *loadEndAlignmentFromDisk(Flower *flower, FILE *fileHandle, End **end);


#endif /* ENDALIGNER_H_ */
