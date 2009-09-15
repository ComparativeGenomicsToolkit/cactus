#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "bioioC.h"
#include "commonC.h"

/*
 * Assistance routine for cactus_setup.py.
 * Routine reads in fasta file progressively, writes each one two a separate sequence
 * file. Returns a list of sequence headers, lengths and file names.
 */

char *dirForFiles;
int32_t counter;
FILE *fileHandle;

void fn(const char *fastaHeader, const char *sequence, int32_t length) {
	FILE *fileHandle2;
	static char cA[100000];
	assert(strlen(dirForFiles) < 900000);
	sprintf(cA, "%s/%i.seq", dirForFiles, counter);
	fileHandle2 = fopen(cA, "w");
	fprintf(fileHandle2, sequence);
	fclose(fileHandle2);
	fprintf(fileHandle, "%s\n", fastaHeader);
	fprintf(fileHandle, "%i.seq\n", counter++);
	fprintf(fileHandle, "%i\n", length);
}

int main(int argc, char *argv[]) {
	assert(argc == 5);
	FILE *fileHandle2;
	dirForFiles = argv[1];
	fileHandle = fopen(argv[2], "w");
	fileHandle2 = fopen(argv[3], "r");
	assert(sscanf(argv[4], "%i", &counter) == 1);
	fastaReadToFunction(fileHandle2, fn);
	fclose(fileHandle);
	fclose(fileHandle2);
	return 0;
}
