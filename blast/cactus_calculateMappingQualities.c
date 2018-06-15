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

int cmpAlignmentsFn(const void *a, const void *b) {
	struct PairwiseAlignment *pA1 = a;
	struct PairwiseAlignment *pA2 = b;
	return pA1->score < pA2->score ? -1 : (pA1->score == pA2->score ? 0 : 1);
}

int main(int argc, char *argv[]) {
	/*
	 * For each alignment in the input file copy the alignment to the output file and additionally
	 * write out the alignment with the query and target sequences reversed.
	 */
	assert(argc == 4);
	st_setLogLevelFromString(argv[1]);

    FILE *fileHandleIn = fopen(argv[2], "r");
    FILE *fileHandleOut = fopen(argv[3], "w");
    int64_t maxAlignmentsPerPosition = argv[4];
    
    // List of totally overlapping alignments
    stList *alignments = stList_construct();
    uint64_t alignmentStart = -1;
    
    struct PairwiseAlignment *pairwiseAlignment;
    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {
    
    	// If the pairwiseAlignment does not share the same interval
    	// as the other pairwise alignments
		if(stList_length(pairwiseAlignments) > 0 && 
		   	getStartCoordinate(pairwiseAlignment) != alignmentStart) {
		   
		   	// 
			updateScoresToReflectMappingQualities(alignments);
			
			// Sort by ascending score
			stList_sort(alignments, cmpAlignmentsFn);
			
			// Report the alignments
			int64_t i=0;
			while(stList_size(alignments) > 0) {
				struct *pairwiseAlignment2 = stList_pop(alignments);
				if(i++ < maxAlignmentsPerPosition) {
					// Write out modified cigar
		        	cigarWrite(fileHandleOut, pairwiseAlignment, 0);
				}

				// Cleanup
				destructPairwiseAlignment(pairwiseAlignment);
			}
		}
		else {
			alignmentStart = getStartCoordinate(pairwiseAlignment);
		}

		// Checks
		assert(stList_length(alignments) == 0);
		assert(alignmentStart == getStartCoordinate(pairwiseAlignment));

		// Adding the pairwise alignment to the set to consider
		stList_append(pairwiseAlignments, pairwiseAlignment);
    }
    
    fclose(fileHandleIn);
    fclose(fileHandleOut);
    return 0;
}

void updateScoresToReflectMappingQualities(stList *alignments) {

}
