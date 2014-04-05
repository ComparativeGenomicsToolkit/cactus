#include <getopt.h>
#include <error.h>
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

static void addSequenceLength(const char *name, const char *seq, int64_t len)
{
    stList_append(sequenceNames, stString_copy(name));
    int64_t *heapLen = malloc(sizeof(int64_t));
    *heapLen = len;
    stHash_insert(sequenceLengths, stString_copy(name), heapLen);
}

static void usage(void)
{
    fprintf(stderr, "cactus_coverage <--fasta or --cactusDisk and --genome> alignmentsFile\n");
}

// Comparison function for sorting pairwise alignments by increasing
// start position on contig 1.
static int sortOnContig1(const void *val1, const void *val2)
{
    struct PairwiseAlignment *pA1 = (struct PairwiseAlignment *) val1;
    struct PairwiseAlignment *pA2 = (struct PairwiseAlignment *) val2;
    int contigsCmp;
    int64_t start1, start2;
    if((contigsCmp = strcmp(pA1->contig1, pA2->contig2))) {
        return contigsCmp;
    }
    if(pA1->strand1) {
        start1 = pA1->start1;
    } else {
        start1 = pA1->end1;
    }
    if(pA2->strand1) {
        start2 = pA2->start1;
    } else {
        start2 = pA2->end1;
    }
    if(start1 > start2) {
        return 1;
    } else if(start1 == start2) {
        return 0;
    } else {
        return -1;
    }
}

// Comparison function for sorting pairwise alignments by increasing
// start position on contig 1.
static int sortOnContig2(const void *val1, const void *val2)
{
    struct PairwiseAlignment *pA1 = (struct PairwiseAlignment *) val1;
    struct PairwiseAlignment *pA2 = (struct PairwiseAlignment *) val2;
    int contigsCmp;
    int64_t start1, start2;
    if ((contigsCmp = strcmp(pA1->contig2, pA2->contig2))) {
        return contigsCmp;
    }
    if(pA1->strand2) {
        start1 = pA1->start2;
    } else {
        start1 = pA1->end2;
    }
    if(pA2->strand2) {
        start2 = pA2->start2;
    } else {
        start2 = pA2->end2;
    }
    if(start1 > start2) {
        return 1;
    } else if(start1 == start2) {
        return 0;
    } else {
        return -1;
    }
}

static int sortBlocksByStart(const void *val1, const void *val2)
{
    struct block *block1 = (struct block *) val1;
    struct block *block2 = (struct block *) val2;
    if(block1->start > block2->start) {
        return 1;
    } else if (block1->start == block2->start) {
        return 0;
    } else {
        return -1;
    }
}

// Combine blocks from a sorted list of blocks.
static void combineBlocks(stList *blocks)
{
    int64_t i;
    for(i = 1; i < stList_length(blocks); i++) {
        struct block *block = stList_get(blocks, i);
        struct block *prevBlock = stList_get(blocks, i - 1);
        struct block *newBlock;
        assert(block->start >= prevBlock->start);
        if(block->start >= prevBlock->end) {
            continue;
        }
        if(block->end > prevBlock->end) {
            if(prevBlock->start != block->start) {
                newBlock = malloc(sizeof(struct block));
                newBlock->start = prevBlock->start;
                newBlock->end = block->start;
                newBlock->value = prevBlock->value;
                stList_append(blocks, newBlock);
            }
            prevBlock->value += block->value;
            prevBlock->start = block->start;
            block->start = prevBlock->end;
            stList_sort(blocks, sortBlocksByStart);
            i = 1;
        } else if(block->end == prevBlock->end) {
            assert(prevBlock->start == block->start);
            prevBlock->value += block->value;
            free(stList_remove(blocks, i));
            // To keep the iteration correct
            i -= 1;
        } else {
            // prevBlock overlaps this block completely
            assert(block->end < prevBlock->end);
            if (prevBlock->start != block->start) {
                newBlock = malloc(sizeof(struct block));
                newBlock->start = prevBlock->start;
                newBlock->end = block->start;
                newBlock->value = prevBlock->value;
                stList_append(blocks, newBlock);
            }
            block->value += prevBlock->value;
            prevBlock->start = block->end;
            stList_sort(blocks, sortBlocksByStart);
            i = 1;
        }
    }
}

// Print and remove all blocks up to (but not including) start
// position upToPos.
static void printAndFlushBlocks(stList *blocks, char *contig, int64_t upToPos)
{
    while(stList_length(blocks) != 0) {
        struct block *block = stList_get(blocks, 0);
        if(block->end >= upToPos) {
            break;
        }
        printf("%s\t%" PRIi64 "\t%" PRIi64 "\t\t%" PRIi64 "\n", contig,
               block->start, block->end, block->value);
        free(stList_remove(blocks, 0));
    }
}

static void printCoverageBed(stList *alignments, int contigNum)
{
    int64_t i;
    stList *curBlocks = stList_construct();
    char *prevContig = NULL;
    for(i = 0; i < stList_length(alignments); i++) {
        struct PairwiseAlignment *pA = stList_get(alignments, i);
        int64_t j, startPos, endPos, curAlignmentPos;
        int strand;
        char *contig;
        if ((contigNum == 1 && !stHash_search(sequenceLengths, pA->contig1)) ||
            (contigNum == 2 && !stHash_search(sequenceLengths, pA->contig2))) {
            // Contig isn't relevant to our query
            continue;
        }
        if(contigNum == 1) {
            startPos = pA->strand1 ? pA->start1 : pA->end1;
            endPos = pA->strand1 ? pA->end1 : pA->start1;
            contig = pA->contig1;
            strand = pA->strand1;
        } else {
            startPos = pA->strand2 ? pA->start2 : pA->end2;
            endPos = pA->strand2 ? pA->end2 : pA->start2;
            contig = pA->contig2;
            strand = pA->strand2;
        }
        if (prevContig && strcmp(prevContig, contig)) {
            // Print all blocks in the list
            printAndFlushBlocks(curBlocks, prevContig, INT64_MAX);
        } else {
            // Print all blocks before this alignment
            printAndFlushBlocks(curBlocks, contig, startPos);
        }
        // add alignment blocks to the list and combine them
        curAlignmentPos = strand ? startPos : endPos;
        for(j = 0; j < pA->operationList->length; j++) {
            struct AlignmentOperation *op = pA->operationList->list[j];
            struct block *newBlock;
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
                newBlock = malloc(sizeof(struct block));
                if(strand) {
                    newBlock->start = curAlignmentPos;
                    newBlock->end = curAlignmentPos + op->length;
                    curAlignmentPos += op->length;
                } else {
                    newBlock->start = curAlignmentPos - op->length;
                    newBlock->end = curAlignmentPos;
                    curAlignmentPos -= op->length;
                }
                newBlock->value = 1;
                stList_append(curBlocks, newBlock);
                break;
            }
        }
        stList_sort(curBlocks, sortBlocksByStart);
        combineBlocks(curBlocks);
        prevContig = contig;
    }
    printAndFlushBlocks(curBlocks, prevContig, INT64_MAX);
}

int main(int argc, char *argv[])
{
    char *fastaPath = NULL;
    struct option opts[] = { {"onlyContig1", no_argument, NULL, '1'},
                             {"onlyContig2", no_argument, NULL, '2'},
                             {0, 0, 0, 0} };
    int outputOnContig1 = TRUE, outputOnContig2 = TRUE;
    int64_t flag;
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
    sequenceNames = stList_construct();

    if (optind >= argc - 1) {
        fprintf(stderr, "fasta file for sequence and alignments file (in "
                "cigar format) must be provided\n");
        return 1;
    }
    fastaPath = argv[optind];
    // Fill in sequence lengths from fasta.
    FILE *fastaHandle = fopen(fastaPath, "r");
    if (!fastaHandle) {
        error(1, errno, "Could not open fasta file %s", fastaPath);
    }
    fastaReadToFunction(fastaHandle, addSequenceLength);
    fclose(fastaHandle);
    FILE *alignmentsHandle = fopen(argv[optind + 1], "r");

    // output bed of covered regions
    stList *alignments = stList_construct();
    for(;;) {
        // load alignments into memory
        struct PairwiseAlignment *pA = cigarRead(alignmentsHandle);
        if(pA == NULL) {
            break;
        }
        stList_append(alignments, pA);
    }
    if (outputOnContig1) {
        stList_sort(alignments, sortOnContig1);
        printCoverageBed(alignments, 1);
    }
    if (outputOnContig2) {
        stList_sort(alignments, sortOnContig2);
        printCoverageBed(alignments, 2);
    }
}
