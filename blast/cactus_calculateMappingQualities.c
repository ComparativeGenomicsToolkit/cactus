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

int cmpAlignmentsFn(const void *a, const void *b) {
	const struct PairwiseAlignment *pA1 = a;
	const struct PairwiseAlignment *pA2 = b;
	return pA1->score < pA2->score ? -1 : (pA1->score == pA2->score ? 0 : 1);
}

void updateScoresToReflectMappingQualities(stList *alignments) {
	//TODO
}

int main(int argc, char *argv[]) {
	/*
	 * For each alignment in the input file copy the alignment to the output file and additionally
	 * write out the alignment with the query and target sequences reversed.
	 */
	assert(argc == 6);
	st_setLogLevelFromString(argv[1]);

    FILE *fileHandleIn = fopen(argv[2], "r");
    FILE *fileHandleOut = fopen(argv[3], "w");
    int64_t maxAlignmentsPerSite;
    int64_t i = sscanf(argv[4], "%" PRIi64 "", &maxAlignmentsPerSite);
    assert(i == 1);
    float minimumMapQValue;
    i = sscanf(argv[5], "%f", &minimumMapQValue);
    assert(i == 1);
    
    // List of totally overlapping alignments
    stList *alignments = stList_construct();
    uint64_t alignmentStart = -1;
    
    struct PairwiseAlignment *pairwiseAlignment = NULL;
    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {
    
    	// If the pairwiseAlignment does not share the same interval
    	// as the other pairwise alignments
		if(stList_length(alignments) > 0 &&
		   	getStartCoordinate(pairwiseAlignment) != alignmentStart) {
		   
		   	// TODO
			updateScoresToReflectMappingQualities(alignments);
			
			// Sort by ascending score
			stList_sort(alignments, cmpAlignmentsFn);
			
			// Report the alignments
			for(i=0; stList_length(alignments) > 0;) {
				struct PairwiseAlignment *pairwiseAlignment2 = stList_pop(alignments);
				if(i++ < maxAlignmentsPerSite && pairwiseAlignment2->score >= minimumMapQValue) {
					// Write out modified cigar
		        	cigarWrite(fileHandleOut, pairwiseAlignment2, 0);
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
		stList_append(alignments, pairwiseAlignment);
    }
    
    fclose(fileHandleIn);
    fclose(fileHandleOut);
    return 0;
}
