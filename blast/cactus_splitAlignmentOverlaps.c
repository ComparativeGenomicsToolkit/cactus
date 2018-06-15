/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "stLastzAlignments.h"

uint64_t getStartCoordinate(struct PairwiseAlignment *pairwiseAlignment) {
	assert(pairwiseAlignment->strand1); // This code assumes that the alignment is reported with respect
	// to the positive strand of the first sequence
	return pairwiseAlignment->start1;
}

uint64_t getEndCoordinate(struct PairwiseAlignment *pairwiseAlignment) {
	assert(pairwiseAlignment->strand1); // This code assumes that the alignment is reported with respect
	// to the positive strand of the first sequence
	return pairwiseAlignment->end1;
}

int main(int argc, char *argv[]) {
	/*
	 * Each alignment has a unique query interval, defined by where it starts and ends on the
	 * query sequence.
	 * Two alignments partially overlap if there query intervals overlap but are not the same.
	 * This program breaks up alignments in the input file so that there are no partial overlaps between
	 * alignments, outputting the non-overlapping alignments to the output file.
	 */
	assert(argc == 4);
	st_setLogLevelFromString(argv[1]);

    FILE *fileHandleIn = fopen(argv[2], "r");
    FILE *fileHandleOut = fopen(argv[3], "w");

    // Set of alignments being progressively processed, ordered by ascending query end coordinate
    stSortedSet *activeAlignments = stSortedSet_construct();

    struct PairwiseAlignment *pairwiseAlignment;
    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {
    	// Remove overlaps in alignments up to but excluding the start of pairwiseAlignment
    	splitAlignmentOverlaps(activeAlignments, getStartCoordinate(pairwiseAlignment), fileHandleOut);

    	// Add pairwiseAlignemnt to the set of activeAlignemnts
    	stSortedSet_insert(activeAlignments, pairwiseAlignment);
    }
    // Remove remaining overlaps in alignments
    splitAlignmentOverlaps(activeAlignments, UINT64_MAX);

    fclose(fileHandleIn);
    fclose(fileHandleOut);
    return 0;
}

void splitAlignmentOverlaps(stSortedSet *activeAlignments, uint64_t splitUpto, FILE *fileHandleOut) {
	// while (minEndCoordinate = Min end coordinate in S) < splitUpto:
	uint64_t minEndCoordinate;
	while(stSortedSet_size(activeAlignments) > 0 &&
			(minEndCoordinate = getEndCoordinate(stSortedSet_getFirst(activeAlignments))) < splitUpto) {
		chopPrefixes(activeAlignments, minEndCoordinate, fileHandleOut);
	}
	if(stSortedSet_size(activeAlignments) > 0) {
		chopPrefixes(activeAlignments, splitUpto, fileHandleOut);
	}
}

void chopPrefixes(stSortedSet *activeAlignments, uint64_t prefixEnd) {
	struct PairwiseAlignment *pairwiseAlignment = stSortedSet_getFirst(activeAlignments);
	//For each alignment A in S that ends before prefix end
	while(stSortedSet_size(activeAlignments) > 0 &&
			getEndCoordinate(pairwiseAlignment = stSortedSet_getFirst(activeAlignments)) <= prefixEnd) {
		// Remove from active alignments
		stSortedSet_remove(activeAlignments, pairwiseAlignment);
		// Write out the alignment
		cigarWrite(prefixAlignment, fileHandle);
		// Delete the pairwise alignment
		pairwiseAlignment_destruct(pairwiseAlignment);
	}
	// Now split the remaining alignments, which must all end after prefixEnd
	stSortedSetIt *it = stSortedSet_getIterator(activeAlignments);
	while((pairwiseAlignment = stSortedSet_getNext(it)) != NULL) {
		assert(getEndCoordinate(pairwiseAlignment) > prefixEnd);
		// Cleave off the prefix of the alignment
		pairwiseAlignment = removeAlignmentPrefix(pairwiseAlignment, prefixEnd);
		// Write out alignment up until prefixEnd
		cigarWrite(prefixAlignment, fileHandle);
		// Delete the prefix pairwise alignment
		pairwiseAlignment_destruct(pairwiseAlignment);
	}
}

struct PairwiseAlignment *removeAlignmentPrefix(struct PairwiseAlignment *pairwiseAlignment, int64_t prefixEnd) {
	// Store the original start coordinates
	int64_t start1 = pairwiseAlignment->start1, start2 = pairwiseAlignment->start2;

	// Split the ops in the cigar string between the prefix and suffix alignments
	struct List *prefixOps = constructEmptyList(0, (void (*)(void *))destructAlignmentOperation);
	if(pairwiseAlignment->start1 < prefixEnd) {
		int64_t i = 0;
		while(pairwiseAlignment->start1 < prefixEnd) {
			struct AlignmentOperation *op = pairwiseAlignment->operationList->list[i];
			assert(op->length > 0);

			// Op is in the prefix alignment
			if(pairwiseAlignment->start1 + op->length <= prefixEnd) {
				listAppend(prefixOps, op);
				i++;
			}
			// Op spans the prefix and suffix alignments, so split it
			else {
				listAppend(prefixOps, constructAlignmentOperation(op->opType, prefixEnd-pairwiseAlignment->start1, op->score));
				op->length -= prefixEnd-pairwiseAlignment->start1;
			}

			// Update start coordinates of suffix alignments
			pairwiseAlignment->start1 += op->length;
			if(op->opType != PAIRWISE_INDEL_X) {
				pairwiseAlignment->start2 += pairwiseAlignment->strand2 ? op->length : -op->length;
			}
		}

		// Remove prefix ops from pairwise alignment
		for(int64_t j=0; i<pairwiseAlignment->operationList->length;) {
			pairwiseAlignment->operationList->list[j++] = pairwiseAlignment->operationList->list[i++];
		}
		pairwiseAlignment->operationList->length -= prefixOps->length;
	}

	// Create prefix pairwiseAlignment
	struct PairwiseAlignment *prefixAlignment = constructPairwiseAlignment(stString_copy(pairwiseAlignment->contig1),
			start1, pairwiseAlignment->start1, 1,
			stString_copy(pairwiseAlignment->contig2), start2, pairwiseAlignment->start2, pairwiseAlignment->strand2,
			pairwiseAlignment->score, prefixOps);

	return prefixAlignment;
}
