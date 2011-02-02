#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "commonC.h"
#include "hashTableC.h"
#include "bioioC.h"
#include "commonC.h"
#include "cactus.h"
#include "avl.h"
#include "pairwiseAlignment.h"

static void convertCoordinatesP(char *header, char **contig1, int32_t *start, int32_t *strand) {
	struct List *attributes = fastaDecodeHeader(header);
	//uglyf(" strings: %s\n", stringsJoin(" $ ", attributes->list, attributes->length));
	assert(attributes->length >= 3);
	assert(sscanf((const char *)attributes->list[attributes->length-1], "%i", start) == 1);
	assert(sscanf((const char *)attributes->list[attributes->length-2], "%i", strand) == 1);
	assert(*strand == 0 || *strand == 1);
	free(attributes->list[--attributes->length]);
	free(attributes->list[--attributes->length]);
	*contig1 = fastaEncodeHeader(attributes);
	destructList(attributes);
}

void convertCoordinates(char *tempCigarFile, FILE *fileHandle) {
	FILE *fileHandle2 = fopen(tempCigarFile, "r");
	struct PairwiseAlignment *pairwiseAlignment;

	while((pairwiseAlignment = cigarRead(fileHandle2)) != NULL) {
		//uglyf(" I have got %s %s %i %i %i %i %i %i\n", pairwiseAlignment->contig1, pairwiseAlignment->contig2, pairwiseAlignment->start1, pairwiseAlignment->start2, pairwiseAlignment->end1, pairwiseAlignment->end2, pairwiseAlignment->strand1, pairwiseAlignment->strand2);
		//Correct coordinates.
		char *contig;
		int32_t start, strand;

		checkPairwiseAlignment(pairwiseAlignment);

		convertCoordinatesP(pairwiseAlignment->contig1, &contig, &start, &strand);
		free(pairwiseAlignment->contig1);
		pairwiseAlignment->contig1 = contig;
		pairwiseAlignment->strand1 = strand ? pairwiseAlignment->strand1 : !pairwiseAlignment->strand1;
		pairwiseAlignment->start1 = strand ? pairwiseAlignment->start1 + start : start - pairwiseAlignment->start1 + 1;
		pairwiseAlignment->end1 = strand ? pairwiseAlignment->end1 + start : start - pairwiseAlignment->end1 + 1;

		convertCoordinatesP(pairwiseAlignment->contig2, &contig, &start, &strand);
		free(pairwiseAlignment->contig2);
		pairwiseAlignment->contig2 = contig;
		pairwiseAlignment->strand2 = strand ? pairwiseAlignment->strand2 : !pairwiseAlignment->strand2;
		pairwiseAlignment->start2 = strand ? pairwiseAlignment->start2 + start : start - pairwiseAlignment->start2 + 1;
		pairwiseAlignment->end2 = strand ? pairwiseAlignment->end2 + start : start - pairwiseAlignment->end2 + 1;

		checkPairwiseAlignment(pairwiseAlignment);

		cigarWrite(fileHandle, pairwiseAlignment, 0);

		//uglyf(" I have got2 %s %s %i %i %i %i %i %i\n", pairwiseAlignment->contig1, pairwiseAlignment->contig2, pairwiseAlignment->start1, pairwiseAlignment->start2, pairwiseAlignment->end1, pairwiseAlignment->end2, pairwiseAlignment->strand1, pairwiseAlignment->strand2);

		destructPairwiseAlignment(pairwiseAlignment);
	}
	fclose(fileHandle2);
}

int main(int argc, char *argv[]) {
	/*
	 * For each cigar in file, update the coordinates and write to the second file.
	 */
	assert(argc == 3);
	FILE *fileHandleIn = fopen(argv[1], "r");
	FILE *fileHandleOut = fopen(argv[2], "w");
	int size = 100;
	char *cA = st_calloc(size+1, sizeof(char));
	int32_t i;
	do {
	    i = benLine(&cA, &size, fileHandleIn);
	    if(strlen(cA) > 0) {
	        convertCoordinates(cA, fileHandleOut);
	    }
	} while(i != -1);
	fclose(fileHandleIn);
	fclose(fileHandleOut);
	free(cA);

	//return 1;

	return 0;
}
