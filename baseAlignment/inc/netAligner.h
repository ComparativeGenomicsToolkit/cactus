/*
 * netAligner.h
 *
 *  Created on: 2 Jul 2010
 *      Author: benedictpaten
 */

#ifndef NETALIGNER_H_
#define NETALIGNER_H_

/*
 * Constructs an alignment for the net by constructing an alignment for each end
 * then filtering the alignments against each other so each position is a member of only one
 * end alignment. Spanning trees controls the number of pairwise alignments used
 * to construct the alignment, maxSequenceLength is the maximum length of a sequence to consider in the end alignment.
 * Model parameters is the parameters of the pairwise alignment model.
 */
stSortedSet *makeNetAlignment(Flower *net, int32_t spanningTrees,
        int32_t maxSequenceLength, void *modelParameters);

#endif /* NETALIGNER_H_ */
