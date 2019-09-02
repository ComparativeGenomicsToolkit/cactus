/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "pairwiseAlignment.h"

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

struct PairwiseAlignment *removeAlignmentPrefix(struct PairwiseAlignment *pairwiseAlignment, int64_t prefixEnd) {
	// Store the original start coordinates
	//checkPairwiseAlignment(pairwiseAlignment);
	int64_t start1 = pairwiseAlignment->start1, start2 = pairwiseAlignment->start2;
	assert(pairwiseAlignment->end1 > prefixEnd);
	assert(pairwiseAlignment->start1 < prefixEnd);
	assert(pairwiseAlignment->strand1);

	// Split the ops in the cigar string between the prefix and suffix alignments
	struct List *prefixOps = constructEmptyList(0, (void (*)(void *))destructAlignmentOperation);
	int64_t i = 0;
	do {
		assert(i < pairwiseAlignment->operationList->length);
		struct AlignmentOperation *op = pairwiseAlignment->operationList->list[i];
		assert(op->length > 0);

		if(op->opType == PAIRWISE_INDEL_Y) { // Insert in second sequence
			listAppend(prefixOps, op);
			i++;
			pairwiseAlignment->start2 += pairwiseAlignment->strand2 ? op->length : -op->length;
		}
		else { // Not an insert in second sequence
			// Op is in the prefix alignment
			int64_t j;
			if(pairwiseAlignment->start1 + op->length <= prefixEnd) {
				listAppend(prefixOps, op);
				i++;
				j = op->length;
			}
			// Op spans the prefix and suffix alignments, so split it
			else {
				j = prefixEnd-pairwiseAlignment->start1;
				assert(j > 0);
				listAppend(prefixOps, constructAlignmentOperation(op->opType, j, op->score));
				op->length -= j;
				assert(op->length > 0);
				assert(pairwiseAlignment->start1+j == prefixEnd);
			}

			// Update start coordinates of suffix alignments
			pairwiseAlignment->start1 += j;
			if(op->opType != PAIRWISE_INDEL_X) {
				pairwiseAlignment->start2 += pairwiseAlignment->strand2 ? j : -j;
			}
		}
	} while(pairwiseAlignment->start1 < prefixEnd);

	assert(pairwiseAlignment->start1 == prefixEnd);

        if (i > 0) {
            // Remove prefix ops from pairwise alignment
            uint64_t j=0;

            while(i < pairwiseAlignment->operationList->length) {
                pairwiseAlignment->operationList->list[j++] = pairwiseAlignment->operationList->list[i++];
            }
            pairwiseAlignment->operationList->length = j;
        }
        assert(pairwiseAlignment->operationList->length > 0);

	// Create prefix pairwiseAlignment
	struct PairwiseAlignment *prefixAlignment = constructPairwiseAlignment(pairwiseAlignment->contig1,
			start1, pairwiseAlignment->start1, 1,
			pairwiseAlignment->contig2, start2, pairwiseAlignment->start2, pairwiseAlignment->strand2,
			pairwiseAlignment->score, prefixOps);

	return prefixAlignment;
}

void emitBlock(stSortedSet *activeAlignments, uint64_t from, uint64_t to, FILE *fileHandleOut) {
	/*
	 * Emits block of alignments that are all start, inclusive, at 'fromt' and end, exclusive, at 'to'.
	 */
	// Now split the remaining alignments, which must all end after prefixEnd
	struct PairwiseAlignment *pairwiseAlignment;
	while(stSortedSet_size(activeAlignments) > 0 &&
		  getStartCoordinate(pairwiseAlignment = stSortedSet_getFirst(activeAlignments)) < to) {
		assert(getStartCoordinate(pairwiseAlignment) == from);
		//assert(getEndCoordinate(pairwiseAlignment) >= to);
		
		// Remove from active alignments
		assert(stSortedSet_search(activeAlignments, pairwiseAlignment) == pairwiseAlignment);
		stSortedSet_remove(activeAlignments, pairwiseAlignment);
		assert(stSortedSet_search(activeAlignments, pairwiseAlignment) == NULL);
		
		// If the alignment needs to be split
		if(getEndCoordinate(pairwiseAlignment) > to) {
			// Cleave off the prefix of the alignment
			struct PairwiseAlignment *prefixPairwiseAlignment = removeAlignmentPrefix(pairwiseAlignment, to);
			assert(getStartCoordinate(prefixPairwiseAlignment) == from);
			assert(getEndCoordinate(prefixPairwiseAlignment) == to);
			assert(getStartCoordinate(pairwiseAlignment) == to);

			// Add back the suffix alignment to the set
			stSortedSet_insert(activeAlignments, pairwiseAlignment);

			// Set the prefix alignment to be written out and then cleaned up
			pairwiseAlignment = prefixPairwiseAlignment;
		}
		
		// Write out the prefix alignment up until end
		cigarWrite(fileHandleOut, pairwiseAlignment, 0);
		// Delete the pairwise alignment
		assert(stSortedSet_search(activeAlignments, pairwiseAlignment) == NULL);
		destructPairwiseAlignment(pairwiseAlignment);
	}
}

void splitAlignmentOverlaps(stSortedSet *activeAlignments, uint64_t splitUpto, FILE *fileHandleOut) {
	if(stSortedSet_size(activeAlignments) == 0) {
		return; // Nothing to do
	}

	// Process overlaps between alignments that precede splitUpto
	uint64_t from = getStartCoordinate(stSortedSet_getFirst(activeAlignments));
	uint64_t to;
	// while (minEndCoordinate = Min end coordinate in S) < splitUpto:
	while(stSortedSet_size(activeAlignments) > 0 &&
		  (to = getEndCoordinate(stSortedSet_getFirst(activeAlignments))) < splitUpto) {
		assert(from < to);
		emitBlock(activeAlignments, from, to, fileHandleOut);
		from = to;
	}

	// Now split at the splitUpto point
	if(stSortedSet_size(activeAlignments) > 0) {
		assert(from < to);
		emitBlock(activeAlignments, from, splitUpto, fileHandleOut);
	}
}

int comparePairwiseAlignments(const void *a, const void *b) {
	struct PairwiseAlignment *pA1 = (struct PairwiseAlignment *)a;
	struct PairwiseAlignment *pA2 = (struct PairwiseAlignment *)b;

	int i = strcmp(pA1->contig1, pA2->contig1);
	if(i == 0) {
		i = pA1->start1 > pA2->start1 ? 1 : pA1->start1 < pA2->start1 ? -1 : 0;
		if(i == 0) {
			i = pA1->end1 > pA2->end1 ? 1 : pA1->end1 < pA2->end1 ? -1 : 0;
			if(i == 0) {
				i = pA1 > pA2 ? 1 : pA1 < pA2 ? -1 : 0;
			}
		}
	}

	return i;
}

int main(int argc, char *argv[]) {
	/*
	 * Each alignment has a unique first sequence interval, defined by where it starts and ends on the
	 * first sequence.
	 * Two alignments partially overlap if their first sequence intervals overlap but are not the same.
	 * This program breaks up alignments in the input file so that there are no partial overlaps between
	 * alignments, outputting the non-partially-overlapping alignments to the output file.
	 */
	st_setLogLevelFromString(argv[1]);

	FILE *fileHandleIn;
	FILE *fileHandleOut;

	if(argc == 2) {
		fileHandleIn = stdin;
		fileHandleOut = stdout;
	}
	else {
		assert(argc == 4);
		fileHandleIn = fopen(argv[2], "r");
		fileHandleOut = fopen(argv[3], "w");
	}

    // Set of alignments being progressively processed, ordered by ascending query end coordinate
    stSortedSet *activeAlignments = stSortedSet_construct3(comparePairwiseAlignments, NULL);

    struct PairwiseAlignment *pairwiseAlignment;
    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {

    	// There are existing alignments
    	if(stSortedSet_size(activeAlignments) > 0) {
    		// If the new alignment is on the same sequence as the previous sequence
			if(strcmp(((struct PairwiseAlignment *)stSortedSet_getFirst(activeAlignments))->contig1,
					pairwiseAlignment->contig1) == 0) {
				// Remove overlaps in alignments up to but excluding the start of pairwiseAlignment
				splitAlignmentOverlaps(activeAlignments, getStartCoordinate(pairwiseAlignment), fileHandleOut);
			}
			else {
				// If pairwiseAlignment is on a new sequence
				splitAlignmentOverlaps(activeAlignments, UINT64_MAX, fileHandleOut);
				assert(stSortedSet_size(activeAlignments) == 0);
			}
    	}

    	// Add pairwiseAlignment to the set of activeAlignemnts
    	stSortedSet_insert(activeAlignments, pairwiseAlignment);
    }
    // Remove remaining overlaps in alignments
    splitAlignmentOverlaps(activeAlignments, UINT64_MAX, fileHandleOut);
    assert(stSortedSet_size(activeAlignments) == 0);

    // Cleanup
    stSortedSet_destruct(activeAlignments);
    if(argc == 4) {
    	fclose(fileHandleIn);
    	fclose(fileHandleOut);
    }

    //while(1);

    return 0;
}
