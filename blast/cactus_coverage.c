#include <getopt.h>
#include <errno.h>
#include "sonLib.h"
#include "bioioC.h"
#include "pairwiseAlignment.h"

// For blocks on the same contig.
struct block {
    int64_t start;
    int64_t end;
    int64_t value;
};

// For calculating total % coverage, and % coverage per seq.
static stHash *sequenceLengths;
static stList *sequenceNames;
static stHash *sequenceCoverage;

static void addSequenceLength(const char *name, const char *seq, int64_t len)
{
    char *identifier = stString_copy(name);
    // lastz only takes the first token of a fasta header as the seq ID.
    // not thread-safe
    identifier = strtok(identifier, " ");
    stList_append(sequenceNames, identifier);
    int64_t *heapLen = malloc(sizeof(int64_t));
    *heapLen = len;
    if(stHash_search(sequenceLengths, identifier) != NULL) {
        fprintf(stderr, "Duplicate sequence identifier %s found: make sure "
                "the first tokens in the headers are unique\n", identifier);
        exit(1);
    }

    // extra copy in case the hash is deleted before the list, or vice
    // versa.
    stHash_insert(sequenceLengths, stString_copy(identifier), heapLen);
}

static void usage(void)
{
    fprintf(stderr, "cactus_coverage <--fasta or --cactusDisk and --genome> alignmentsFile\n");
}

static void printCoverage(char *name, uint16_t *array, int64_t length) {
    int64_t i, regionStart = 0;
    uint16_t prevCoverage = 0;
    for(i = 0; i < length; i++) {
        if(array[i] != prevCoverage) {
            if(prevCoverage != 0) {
                printf("%s\t%" PRIi64 "\t%" PRIi64 "\t\t%u\n", name, 
                       regionStart, i, prevCoverage);
            }
            regionStart = i;
        }
        prevCoverage = array[i];
    }
}

static void fillCoverage(struct PairwiseAlignment *pA, int contigNum,
                         uint16_t *coverageArray)
{
    int strand = contigNum == 1 ? pA->strand1 : pA->strand2;
    int64_t startPos = contigNum == 1 ? pA->start1 : pA->start2;
    int64_t endPos = contigNum == 1 ? pA->end1 : pA->end2;
    (void)endPos; //Avoid compile warning of unused variable.
    int64_t i, j;
    int64_t *lenPtr = stHash_search(sequenceLengths, contigNum == 1 ? pA->contig1 : pA->contig2);
    assert(lenPtr != NULL);
    int64_t len = *lenPtr;
    if(endPos > len) {
        fprintf(stderr, "Error: alignment on %s:%" PRIi64 "-%" PRIi64 " is past chr end\n", contigNum == 1 ? pA->contig1 : pA->contig2, startPos, endPos);
        exit(1);
    }
    int64_t curAlignmentPos = startPos;
    for(i = 0; i < pA->operationList->length; i++) {
        struct AlignmentOperation *op = pA->operationList->list[i];
        switch(op->opType) {
        case PAIRWISE_INDEL_Y:
            if(contigNum == 2) {
                if(strand) {
                    curAlignmentPos += op->length;
                } else {
                    curAlignmentPos -= op->length;
                }
            }
            break;
        case PAIRWISE_INDEL_X:
            if(contigNum == 1) {
                if(strand) {
                    curAlignmentPos += op->length;
                } else {
                    curAlignmentPos -= op->length;
                }
            }
            break;
        case PAIRWISE_MATCH:
            if(strand) {
                for(j=curAlignmentPos; j < curAlignmentPos + op->length; j++) {
                    if(coverageArray[j] == 65535) {
                        fprintf(stderr, "ERROR: Coverage overflow on contig: "
                                "%s pos: %" PRIi64 "\n",
                                contigNum == 1 ? pA->contig1 : pA->contig2, j);
                        exit(1);
                    }
                    coverageArray[j]++;
                }
                curAlignmentPos += op->length;
                assert(curAlignmentPos <= endPos);
            } else {
                for(j=curAlignmentPos-1; j >= curAlignmentPos - op->length; j--) {
                    if(coverageArray[j] == 65535) {
                        fprintf(stderr, "ERROR: Coverage overflow on contig: "
                                "%s pos: %" PRIi64 "\n",
                                contigNum == 1 ? pA->contig1 : pA->contig2, j);
                        exit(1);
                    }
                    coverageArray[j]++;
                }
                curAlignmentPos -= op->length;
                assert(curAlignmentPos >= endPos);
            }
        }
    }
}

int main(int argc, char *argv[])
{
    char *fastaPath = NULL;
    struct option opts[] = { {"onlyContig1", no_argument, NULL, '1'},
                             {"onlyContig2", no_argument, NULL, '2'},
                             {0, 0, 0, 0} };
    int outputOnContig1 = TRUE, outputOnContig2 = TRUE;
    int64_t flag, i ;
    while((flag = getopt_long(argc, argv, "", opts, NULL)) != -1) {
        switch(flag) {
        case '1':
            outputOnContig2 = FALSE;
            break;
        case '2':
            outputOnContig1 = FALSE;
            break;
        case '?':
        default:
            usage();
            return 1;
        }
    }
    if(!(outputOnContig1 || outputOnContig2)) {
        fprintf(stderr, "--onlyContig1 and --onlyContig2 options are "
                "mutually exclusive\n");
        return 1;
    }

    sequenceLengths = stHash_construct3(stHash_stringKey,
                                        stHash_stringEqualKey, free, free);
    sequenceCoverage = stHash_construct3(stHash_stringKey,
                                         stHash_stringEqualKey, free, free);
    sequenceNames = stList_construct3(0, free);

    if (optind >= argc - 1) {
        fprintf(stderr, "fasta file for sequence and alignments file (in "
                "cigar format) must be provided\n");
        return 1;
    }
    fastaPath = argv[optind];
    // Fill in sequence lengths from fasta.
    FILE *fastaHandle = fopen(fastaPath, "r");
    if (!fastaHandle) {
    	st_errAbort("Could not open fasta file %s", fastaPath);
    }
    fastaReadToFunction(fastaHandle, addSequenceLength);
    fclose(fastaHandle);

    // Fill coverage arrays with the alignments
    FILE *alignmentsHandle = fopen(argv[optind + 1], "r");
    for(;;) {
        int64_t *lengthPtr;
        struct PairwiseAlignment *pA = cigarRead(alignmentsHandle);
        if(pA == NULL) {
            // Reached end of alignment file
            break;
        }
        if(outputOnContig1 && (lengthPtr = stHash_search(sequenceLengths, pA->contig1))) {
            // contig 1 is present in the fasta
            uint16_t *array;
            if((array = stHash_search(sequenceCoverage, pA->contig1)) == NULL) {
                // Coverage array for this seq doesn't exist yet, so
                // create it
                int64_t length = *lengthPtr;
                array = st_malloc(length * sizeof(uint16_t));
                memset(array, 0, length * sizeof(uint16_t));
                stHash_insert(sequenceCoverage, stString_copy(pA->contig1),
                              array);
            }
            fillCoverage(pA, 1, array);
        } else if(outputOnContig2 && (lengthPtr = stHash_search(sequenceLengths, pA->contig2))) {
            //contig 2 is present in the fasta
            uint16_t *array;
            if((array = stHash_search(sequenceCoverage, pA->contig2)) == NULL) {
                // Coverage array for this seq doesn't exist yet, so
                // create it
                int64_t length = *lengthPtr;
                array = st_malloc(length * sizeof(uint16_t));
                memset(array, 0, length * sizeof(uint16_t));
                stHash_insert(sequenceCoverage, stString_copy(pA->contig2),
                    array);
            }
            fillCoverage(pA, 2, array);
        }
        destructPairwiseAlignment(pA);
    }
    fclose(alignmentsHandle);

    // Print results as BED
    for(i = 0; i < stList_length(sequenceNames); i++) {
        uint16_t *array;
        char *name = stList_get(sequenceNames, i);
        int64_t *lengthPtr = stHash_search(sequenceLengths, name);
        assert(lengthPtr != NULL);
        int64_t length = *lengthPtr;
        if((array = stHash_search(sequenceCoverage, name))) {
            printCoverage(name, array, length);
        }
    }

    // Cleanup
    stList_destruct(sequenceNames);
    stHash_destruct(sequenceCoverage);
    stHash_destruct(sequenceLengths);
}
