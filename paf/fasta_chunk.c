/*
 * fasta_chunk: Break a set of sequence files into a series of overlapping chunks
 * suitable for parallel computation.
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <getopt.h>
#include <time.h>
#include "bioioC.h"
#include "commonC.h"
#include "sonLib.h"

static FILE *chunkFileHandle = NULL;
static const char *chunksDir = "./temp_fastas";
static int64_t chunkSize = 10000000;
static int64_t chunkOverlapSize = 100000;
static int64_t chunkNo = 0;
static int64_t chunkRemaining; // must be initialized
static char *tempChunkFile = NULL;

void usage() {
    fprintf(stderr, "fasta_chunk [fasta_file]xN [options], version 0.1\n");
    fprintf(stderr, "Breaks up a set of fasta sequences into a series of overlapping chunks, "
                    "printing the names of each chunk file to standard out.\n"
                    "To encode the chunking information each faster header is appended |x,"
                    "where x is the start coordinate (0-based) of the chunked sequence in the original sequence");
    fprintf(stderr, "-c --chunkSize : The chunk size, by default: %" PRIi64 "\n", chunkSize);
    fprintf(stderr, "-o --overlap : The chunk overlap size, by default: %" PRIi64 "\n", chunkOverlapSize);
    fprintf(stderr, "-d --dir : An empty directory to place the chunk files in, by default: %s\n", chunksDir);
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

void finishChunkingSequences() {
    if (chunkFileHandle != NULL) {
        fclose(chunkFileHandle);
        fprintf(stdout, "%s\n", tempChunkFile); // Print the output file to stdout
        free(tempChunkFile);
        tempChunkFile = NULL;
        chunkFileHandle = NULL;
    }
}

static void updateChunkRemaining(int64_t seqLength) {
    //Update remaining portion of the chunk.
    assert(seqLength >= 0);
    chunkRemaining -= seqLength;
    if (chunkRemaining <= 0) {
        finishChunkingSequences();
        chunkRemaining = chunkSize;
    }
}

static int64_t processSubsequenceChunk(char *fastaHeader, int64_t start, char *sequence, int64_t seqLength, int64_t lengthOfChunkRemaining) {
    if (chunkFileHandle == NULL) {
        tempChunkFile = stString_print("%s/%" PRIi64 "", chunksDir, chunkNo++);
        chunkFileHandle = fopen(tempChunkFile, "w");
    }

    int64_t i = 0;
    fastaHeader = stString_copy(fastaHeader);
    while (fastaHeader[i] != '\0') {
        if (fastaHeader[i] == ' ' || fastaHeader[i] == '\t') {
            fastaHeader[i] = '\0';
            break;
        }
        i++;
    }
    char *chunkHeader = stString_print("%s|%" PRIi64 "|%" PRIi64 "", fastaHeader, seqLength, start);
    free(fastaHeader);
    assert(lengthOfChunkRemaining <= chunkSize);
    assert(start >= 0);
    int64_t lengthOfSubsequence = lengthOfChunkRemaining;
    if (start + lengthOfChunkRemaining > seqLength) {
        lengthOfSubsequence = seqLength - start;
    }
    assert(lengthOfSubsequence > 0);
    char c = sequence[start + lengthOfSubsequence];
    sequence[start + lengthOfSubsequence] = '\0';
    fastaWrite(&sequence[start], chunkHeader, chunkFileHandle);
    //fprintf(chunkFileHandle, "%s\n", &sequence[start]);
    free(chunkHeader);
    sequence[start + lengthOfSubsequence] = c;

    updateChunkRemaining(lengthOfSubsequence);
    return lengthOfSubsequence;
}

void processSequenceToChunk(void* dest, const char *fastaHeader, const char *sequence, int64_t sequenceLength) {
    if (sequenceLength > 0) {
        int64_t lengthOfSubsequence = processSubsequenceChunk((char *) fastaHeader, 0, (char *) sequence, sequenceLength, chunkRemaining);
        while (sequenceLength - lengthOfSubsequence > 0) {
            //Make the non overlap file
            int64_t lengthOfFollowingSubsequence = processSubsequenceChunk((char *) fastaHeader, lengthOfSubsequence, (char *) sequence, sequenceLength, chunkRemaining);

            //Make the overlap file
            if(chunkOverlapSize > 0) {
                int64_t i = lengthOfSubsequence - chunkOverlapSize / 2;
                if (i < 0) {
                    i = 0;
                }
                processSubsequenceChunk((char *) fastaHeader, i, (char *) sequence, sequenceLength, chunkOverlapSize);
            }
            lengthOfSubsequence += lengthOfFollowingSubsequence;
        }
    }
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "chunkSize", required_argument, 0, 'c' },
                                                { "overlap", required_argument, 0, 'o' },
                                                { "dir", required_argument, 0, 'd' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:c:o:d:h", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'c':
                chunkSize = atoi(optarg);
                break;
            case 'o':
                chunkOverlapSize = atoi(optarg);
                break;
            case 'd':
                chunksDir = optarg;
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Chunks output directory : %s\n", chunksDir);
    st_logInfo("Chunk size : %" PRIi64 "\n", chunkSize);
    st_logInfo("Chunk overlap size : %" PRIi64 "\n", chunkOverlapSize);

    //////////////////////////////////////////////
    // Make the output directory
    //////////////////////////////////////////////

    if(stFile_exists(chunksDir)) {
        if(!stFile_isDir(chunksDir)) { // Directory does not yet exist
            fprintf(stderr, "Output directory is not a directory: %s", chunksDir);
            exit(1);
        }
        stList *files = stFile_getFileNamesInDirectory(chunksDir);
        if(stList_length(files) != 0) {
            fprintf(stderr, "Output directory is not empty, please specify an empty directory ");
            exit(1);
        }
        stList_destruct(files);
    }
    else {
        st_logCritical("Output directory does not exist, trying to create it: %s\n", chunksDir);
        stFile_mkdir(chunksDir); // Make output directory
    }

    //////////////////////////////////////////////
    // Chunk the sequences
    //////////////////////////////////////////////

    chunkRemaining = chunkSize;
    while(optind < argc) {
        char *seq_file = argv[optind++];
        st_logInfo("Chunking sequence file : %s\n", seq_file);
        FILE *fileHandle2 = fopen(seq_file, "r");
        fastaReadToFunction(fileHandle2, NULL, processSequenceToChunk);
        fclose(fileHandle2);
    }
    finishChunkingSequences();

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    st_logInfo("Fasta chunk is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);
    return 0;
}
