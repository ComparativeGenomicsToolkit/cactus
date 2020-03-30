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

// For calculating coverage on the target genome
static stHash *sequenceLengths = NULL;
static stList *sequenceNames = NULL;
static stHash *sequenceCoverage = NULL;
// For determining if a sequence belongs to the "query" genome
// (although there is no relation to the query contig in the cigar):
// i.e. the genome specified in --from, if any
static stSet *otherGenomeSequences = NULL;
// For splitting sequence coverage arrays by ID, if we're using the
// --depthByID option.
static stHash *IDToSequenceCoverage;

// Add a sequence from the genome to sequenceLength and sequenceNames
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

static void addOtherGenomeSequence(const char *name, const char *seq,
                                   int64_t len)
{
    // Only use the first space-separated token to identify seq.
    // not thread-safe
    char *identifier = stString_copy(name);
    identifier = strtok(identifier, " ");
    stSet_insert(otherGenomeSequences, identifier);
}

static void usage(void)
{
    fprintf(stderr, "cactus_coverage fastaFile alignmentsFile\n");
    fprintf(stderr, "Prints a bed file representing coverage from a CIGAR file "
            "on the sequences provided in the fasta file.\n");
    fprintf(stderr, "Format: seq\tregionStart\tregionStop\tcoverageDepth");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "--onlyContig1: Only print coverage that occurs when a "
            "sequence from the fasta is in the contig1 field of the CIGAR.\n");
    fprintf(stderr, "--onlyContig2: Only print coverage that occurs when a "
            "sequence from the fasta is in the contig2 field of the CIGAR.\n");
    fprintf(stderr, "--depthById: Assume that headers have an 'id=N|' prefix, "
            "where N is an integer. Score coverage depth by the number of "
            "different prefixes that align to a region, rather than the total "
            "number of alignments. Uses much more memory than the standard mode."
            "\n");
    fprintf(stderr, "--from <fromFastaFile>: Only consider alignments for which one sequence is in fastaFile and the other is in fromFastaFile (multiple allowed).\n");
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
    if(prevCoverage != 0) {
        printf("%s\t%" PRIi64 "\t%" PRIi64 "\t\t%u\n", name,
               regionStart, i, prevCoverage);
    }
}

// Fill in the part of a coverage array that is covered by a
// particular pairwise alignment. contigNum is which contig this
// coverage array corresponds to in the CIGAR.
static void fillCoverage(struct PairwiseAlignment *pA, int contigNum,
                         uint16_t *coverageArray)
{
    int strand = contigNum == 1 ? pA->strand1 : pA->strand2;
    int64_t startPos = contigNum == 1 ? pA->start1 : pA->start2;
    int64_t endPos = contigNum == 1 ? pA->end1 : pA->end2;
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
                        fprintf(stderr, "WARNING: Coverage hit cap (65535) on contig: "
                                "%s pos: %" PRIi64 "\n",
                                contigNum == 1 ? pA->contig1 : pA->contig2, j);
                    } else {
                        coverageArray[j]++;
                    }
                }
                curAlignmentPos += op->length;
                assert(curAlignmentPos <= endPos);
            } else {
                for(j=curAlignmentPos-1; j >= curAlignmentPos - op->length; j--) {
                    if(coverageArray[j] == 65535) {
                        fprintf(stderr, "WARNING: Coverage hit cap (65535) on contig: "
                                "%s pos: %" PRIi64 "\n",
                                contigNum == 1 ? pA->contig1 : pA->contig2, j);
                    } else {
                        coverageArray[j]++;
                    }
                }
                curAlignmentPos -= op->length;
                assert(curAlignmentPos >= endPos);
            }
        }
    }
}

// Get the proper coverage array to fill in, given the "on" header
// (i.e. a header in the fasta provided in the arguments to this
// program), and the "from" header (the other header in the CIGAR file,
// which may or may not be in that fasta). Initialize the array if necessary.
static uint16_t *getCoverageArray(char *onHeader, char *fromHeader,
                                  int depthById) {
    stHash *coverageHash;
    if (depthById) {
        // We're splitting arrays by "id=N|" of the "from" header.
        stList *attributes = fastaDecodeHeader(fromHeader);
        char *id = stList_get(attributes, 0);
        if (strncmp(id, "id=", 3)) {
            st_errAbort("Using --depthById mode, but header %s does not have an "
                        "'id=N|' prefix", fromHeader);
        }
        stHash *sequenceSubCoverage = stHash_search(IDToSequenceCoverage, id);
        if (sequenceSubCoverage == NULL) {
            // Initialize coverage sub-hash.
            sequenceSubCoverage = stHash_construct3(stHash_stringKey,
                                                    stHash_stringEqualKey, free,
                                                    free);
            stHash_insert(IDToSequenceCoverage, id, sequenceSubCoverage);
        }
        coverageHash = sequenceSubCoverage;
    } else {
        coverageHash = sequenceCoverage;
    }
    int64_t *lengthPtr = stHash_search(sequenceLengths, onHeader);
    assert(lengthPtr != NULL);
    uint16_t *array;
    if((array = stHash_search(coverageHash, onHeader)) == NULL) {
        // Coverage array for this seq doesn't exist yet, so
        // create it
        int64_t length = *lengthPtr;
        array = st_malloc(length * sizeof(uint16_t));
        memset(array, 0, length * sizeof(uint16_t));
        stHash_insert(coverageHash, stString_copy(onHeader),
                      array);
    }
    return array;
}

int main(int argc, char *argv[])
{
    char *fastaPath = NULL;
    stList* otherGenomeFastaPaths = stList_construct3(0, free);
    struct option opts[] = { {"onlyContig1", no_argument, NULL, '1'},
                             {"onlyContig2", no_argument, NULL, '2'},
                             {"depthById", no_argument, NULL, 'i'},
                             {"from", required_argument, NULL, 'f'},
                             {0, 0, 0, 0} };
    int outputOnContig1 = TRUE, outputOnContig2 = TRUE, depthById = FALSE;
    int64_t flag, i;
    while((flag = getopt_long(argc, argv, "", opts, NULL)) != -1) {
        switch(flag) {
        case '1':
            outputOnContig2 = FALSE;
            break;
        case '2':
            outputOnContig1 = FALSE;
            break;
        case 'i':
            depthById = TRUE;
            break;
        case 'f':
            stList_append(otherGenomeFastaPaths, stString_copy(optarg));
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

    if(stList_length(otherGenomeFastaPaths) > 0) {
        otherGenomeSequences = stSet_construct3(stHash_stringKey,
                                                stHash_stringEqualKey, free);

        stListIterator* it = stList_getIterator(otherGenomeFastaPaths);
        char* otherGenomeFastaPath;
        while ((otherGenomeFastaPath = stList_getNext(it)) != NULL) {
            // Load the "query" genome (we will only look for alignments
            // from/to the main genome and this other genome).
            FILE *otherGenomeFastaHandle = fopen(otherGenomeFastaPath, "r");
            if(otherGenomeFastaHandle == NULL) {
                st_errAbort("Could not open fasta file %s", fastaPath);
            }
            fastaReadToFunction(otherGenomeFastaHandle, addOtherGenomeSequence);
            fclose(otherGenomeFastaHandle);
        }
        stList_destructIterator(it);
    }
    stList_destruct(otherGenomeFastaPaths);
    otherGenomeFastaPaths = NULL;

    sequenceLengths = stHash_construct3(stHash_stringKey,
                                        stHash_stringEqualKey, free, free);
    sequenceCoverage = stHash_construct3(stHash_stringKey,
                                         stHash_stringEqualKey, free, free);
    sequenceNames = stList_construct3(0, free);
    IDToSequenceCoverage = stHash_construct3(stHash_stringKey,
                                             stHash_stringEqualKey,
                                             free,
                                             (void (*)(void *)) stHash_destruct);

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
        if((outputOnContig1 && (lengthPtr = stHash_search(sequenceLengths, pA->contig1))) && ((otherGenomeSequences == NULL) || stSet_search(otherGenomeSequences, pA->contig2))) {
            // contig 1 is present in the fasta and contig 2 is in the
            // "from" genome if it exists
            uint16_t *array = getCoverageArray(pA->contig1, pA->contig2,
                                               depthById);
            fillCoverage(pA, 1, array);
        }
        if((outputOnContig2 && (lengthPtr = stHash_search(sequenceLengths, pA->contig2))) && ((otherGenomeSequences == NULL) || stSet_search(otherGenomeSequences, pA->contig1))) {
            // contig 2 is present in the fasta and contig 1 is in the
            // "from" genome if it exists
            uint16_t *array = getCoverageArray(pA->contig2, pA->contig1,
                                               depthById);
            fillCoverage(pA, 2, array);
        }
        destructPairwiseAlignment(pA);
    }
    fclose(alignmentsHandle);

    if (depthById) {
        // Have to merge all coverage arrays that are divided by
        // source ID into the main sequenceCoverage hash.
        stHashIterator *idIt = stHash_getIterator(IDToSequenceCoverage);
        char *id;
        while ((id = stHash_getNext(idIt)) != NULL) {
            stHash *subHash = stHash_search(IDToSequenceCoverage, id);
            assert(subHash != NULL);
            stHashIterator *sequenceIt = stHash_getIterator(subHash);
            char *sequence;
            while ((sequence = stHash_getNext(sequenceIt)) != NULL) {
                uint16_t *srcArray = stHash_search(subHash, sequence);
                assert(srcArray != NULL);
                uint16_t *destArray = getCoverageArray(sequence, NULL, FALSE);
                int64_t *length = stHash_search(sequenceLengths, sequence);
                assert(length != NULL && *length >= 0);
                for (int64_t i = 0; i < *length; i++) {
                    if (srcArray[i] > 0) {
                        destArray[i]++;
                    }
                }
            }
        }
        stHash_destruct(IDToSequenceCoverage);
    }

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
    if(otherGenomeSequences) {
//        stSet_destruct(otherGenomeSequences);
    }
}
