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
                    "where x is the start coordinate (0-based) of the chunked sequence in the original sequence\n");
    fprintf(stderr, "-c --chunkSize : The chunk size, by default: %" PRIi64 "\n", chunkSize);
    fprintf(stderr, "-o --overlap : The chunk overlap size (each chunked file will have length chunkSize + overlap), by default: %" PRIi64 "\n", chunkOverlapSize);
    fprintf(stderr, "-d --dir : An empty directory to place the chunk files in, by default: %s\n", chunksDir);
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

void startChunkingSequences() {
    if(chunkFileHandle == NULL) {
        assert(tempChunkFile == NULL);
        tempChunkFile = stString_print("%s/%"
        PRIi64
        "", chunksDir, chunkNo++);
        chunkFileHandle = fopen(tempChunkFile, "w");
        chunkRemaining = chunkSize;
        st_logDebug("Starting chunk %s\n", tempChunkFile);
    }
}

void finishChunkingSequences() {
    if (chunkFileHandle != NULL) {
        fclose(chunkFileHandle);
        st_logDebug("Finishing chunk %s\n", tempChunkFile);
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
    }
}

void processSequenceToChunk(void* dest, const char *fastaHeader, const char *sequence, int64_t sequenceLength) {
    // For each chunk of these sequence
    assert(chunkSize > chunkOverlapSize);
    for (int64_t i = 0; i < sequenceLength; i += chunkSize) {
        // be ready to print more sequence
        startChunkingSequences();

        // print the header to the file
        fprintf(chunkFileHandle, ">%s|%" PRIi64 "|%" PRIi64 "\n", fastaHeader, sequenceLength, i);

        // Get end of chunk, including the extra overlap
        int64_t j = (i + chunkSize + chunkOverlapSize) <= sequenceLength ? (i + chunkSize + chunkOverlapSize) : sequenceLength;

        // Get chunk sequence
        char *seq_chunk = stString_getSubString(sequence, i, j-i);
        assert(strlen(seq_chunk) == j - i);

        // print the sequence to the file
        fprintf(chunkFileHandle, "%s\n", seq_chunk);
        free(seq_chunk); // cleanup the fragment

        // update the remaining chunk
        updateChunkRemaining(j - i);
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
                chunkSize = atol(optarg);
                break;
            case 'o':
                chunkOverlapSize = atol(optarg);
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
