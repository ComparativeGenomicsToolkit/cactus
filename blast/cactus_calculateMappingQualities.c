/*
 * Copyright (C) 2009-2018 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "sonLib.h"
#include "pairwiseAlignment.h"
#include "math.h"

uint64_t getStartCoordinate(struct PairwiseAlignment *pairwiseAlignment) {
	assert(pairwiseAlignment->strand1); // This code assumes that the alignment is reported with respect
	// to the positive strand of the first sequence
	return pairwiseAlignment->start1;
}

int cmpAlignmentsFn(const void *a, const void *b) {
	const struct PairwiseAlignment *pA1 = a;
	const struct PairwiseAlignment *pA2 = b;
	return pA1->score < pA2->score ? -1 : (pA1->score > pA2->score ? 1 : 0);
}

void updateScoresToReflectMappingQualities(stList *alignments, float alpha, uint64_t numAlignmentsToScore) {
	// Create an array of the scores
	float *alignmentScores = st_calloc(stList_length(alignments), sizeof(float));
	for(uint64_t i=0; i<stList_length(alignments); i++) {
		alignmentScores[i] = ((struct PairwiseAlignment *)stList_get(alignments, i))->score;
	}

	// Calculate mapQs for the best N alignments (N = numAlignmentsToScore).
        uint64_t start = stList_length(alignments) > numAlignmentsToScore ? stList_length(alignments) - numAlignmentsToScore : 0;
	for(uint64_t i=start; i<stList_length(alignments); i++) {
		struct PairwiseAlignment *pA = stList_get(alignments, i);

		// Cut off the calculation if clearly going to be zero
		if(alpha * (alignmentScores[i] - alignmentScores[stList_length(alignments)-1]) < -10) {
			pA->score = 0.0;
		}

		else {
			// Calculate the denominator
			double z = 0.0;
			for(uint64_t j=0; j<stList_length(alignments); j++) {
				z += pow(10, alpha * (alignmentScores[j] - alignmentScores[i]));
			}
			assert(z >= 1.0);

			if(z <= 1.000001) { // Round scores to max of 60
				pA->score = 60.0;
			}
			else {
				pA->score = -10.0 * log10(1.0 - 1.0/z);
				assert(pA->score >= 0.0);
			}
		}
	}

	// Cleanup
	free(alignmentScores);
}

void reportAlignments(stList *alignments, int64_t maxAlignmentsPerSite,
		float minimumMapQValue, float alpha, FILE **fileHandleOuts) {
	// Sort by ascending score
	stList_sort(alignments, cmpAlignmentsFn);

	// Calculate the mapping qualities
	updateScoresToReflectMappingQualities(alignments, alpha, maxAlignmentsPerSite);

	// Report the alignments
	for(int64_t i=0; stList_length(alignments) > 0;) {
		struct PairwiseAlignment *pairwiseAlignment = stList_pop(alignments);
		if(i < maxAlignmentsPerSite && pairwiseAlignment->score >= minimumMapQValue) {
			// Write out modified cigar
			cigarWrite(fileHandleOuts[i++], pairwiseAlignment, 0);
		}

		// Cleanup
		destructPairwiseAlignment(pairwiseAlignment);
	}
}

int main(int argc, char *argv[]) {
	/*
	 * For each alignment in the input file copy the alignment to the output file and additionally
	 * write out the alignment with the query and target sequences reversed.
	 */
	st_setLogLevelFromString(argv[1]);

	int64_t maxAlignmentsPerSite;
	int64_t i = sscanf(argv[2], "%" PRIi64 "", &maxAlignmentsPerSite);
	assert(i == 1);

	float minimumMapQValue;
	i = sscanf(argv[3], "%f", &minimumMapQValue);
	assert(i == 1);

	float alpha;
	i = sscanf(argv[4], "%f", &alpha);
	assert(i == 1);

	FILE **fileHandleOuts = st_malloc(sizeof(FILE *) * maxAlignmentsPerSite);
	for(i=0; i<maxAlignmentsPerSite; i++) {
		fileHandleOuts[i] = fopen(argv[i+5], "w");
	}

	FILE *fileHandleIn = stdin;
	if(argc == maxAlignmentsPerSite+6) {
		fileHandleIn = fopen(argv[maxAlignmentsPerSite+5], "r");
	}
	else {
		assert(argc == maxAlignmentsPerSite+5);
	}
    
    // List of totally overlapping alignments
    stList *alignments = stList_construct();
    
    struct PairwiseAlignment *pairwiseAlignment = NULL;
    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {

    	// If the pairwiseAlignment does not share the same interval
    	// as the previous pairwise alignments report the previous alignments
		if(stList_length(alignments) == 0 ||
			strcmp(((struct PairwiseAlignment *)stList_peek(alignments))->contig1, pairwiseAlignment->contig1) != 0 ||
		   	getStartCoordinate(stList_peek(alignments)) != getStartCoordinate(pairwiseAlignment)) {
		   
			reportAlignments(alignments, maxAlignmentsPerSite, minimumMapQValue, alpha, fileHandleOuts);
		}

		// Adding the pairwise alignment to the set to consider
		stList_append(alignments, pairwiseAlignment);
    }
    
    reportAlignments(alignments, maxAlignmentsPerSite, minimumMapQValue, alpha, fileHandleOuts);

    assert(stList_length(alignments) == 0);
    // Cleanup
    stList_destruct(alignments);
    for(i=0; i<maxAlignmentsPerSite; i++) {
    	fclose(fileHandleOuts[i]);
    }
    if(argc == maxAlignmentsPerSite+6) {
    	fclose(fileHandleIn);
    }

    //while(1);

    return 0;
}
