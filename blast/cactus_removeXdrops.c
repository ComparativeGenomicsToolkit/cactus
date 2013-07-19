/*
 * Script implements algorithm in "Post-processing long pairwise alignments", Zhang, Berman, Wiehe and Miller, 1999, Bioinformatics.
 *
 * Using language of paper, script processes cigar alignments, using lastz scoring matrix, and reports all
 * full normal-rises containing normal-drops scoring at most X, where X is an input parameter.
 *
 */

/*
 * Data structure that represents the "useful tree" in the paper.
 */
typedef struct _segment Segment;

struct _segment {
   int64_t score; //Called sigma in paper
   int64_t minScore; //Call sigma_* in the paper
   struct PairwiseAlignment *pA;
   Segment *child1, child2, child3;
};

void segment_destruct(Segment *segment) {
    if(segment != NULL) {
        segment_destruct(segment->child1);
        segment_destruct(segment->child2);
        segment_destruct(segment->child3);
        if(segment->pA != NULL) {
            destructPairwiseAlignment(segment->pA);
        }
        free(segment);
    }
}

void segment_lineariseUsefulTree(Segment *segment, stList *l) {
    /*
     * Puts the leaves of the subtree into a pre-order sequence.
     */
    if(segment->child1 != NULL) {
        segment_lineariseUsefulTree(segment->child1, l);
        segment_lineariseUsefulTree(segment->child2, l);
        segment_lineariseUsefulTree(segment->child3, l);
    }
    else {
        stList_append(l, segment);
    }
}

struct PairwiseAlignment *segment_getSubAlignment(Segment *segment) {
    /*
     * Gets a cigar subalignment representing the alignment covered by the cigar.
     */
    stList *l = stList_construct();
    segment_lineariseUsefulTree(segment, l); //Get a list of the alignment operations to make into a cigar.
    //These two segments are used to set the coordinates
    Segment *leftMost = stList_get(l, 0);
    Segment *rightMost = stList_peek(l);
    //Now we must convert the sequence of segments to alignment operations
    struct List *opList = constructEmptyList(0, (void (*)(void *))destructAlignmentOperation);
    for(int64_t i=0; i<stList_length(l); i++) {
        Segment *segment2 = stList_get(l, i);
        for(int64_t j=0; j<segment2->pA->operationList->length; j++) {
            List_append(opList, segment2->pA->operationList->list[j]);
        }
    }
    //Make the alignment operations
    pA = constructPairwiseAlignment(segment->pA->contig1, leftMost->pA->start1, rightMost->pA->end1, segment->pA->strand1,
                                    segment->pA->contig2, leftMost->pA->start2, rightMost->pA->end2, segment->pA->strand2,
                                    segment->score, opList);
    checkPairwiseAlignment(pA); //This does some sanity checking
    //Cleanup
    stList_destruct(l);
    return pA;
}

void segment_reportXFullSubalignments(Segment *segment, struct PairwiseAlignment *pA, int64_t X, int64_t Y, FILE *output) {
    /*
     * Prints the X-full sub alignments in the useful tree. These are by definition non-overlapping.
     */
    if(segment != NULL) {
        if(segment->minScore >= X) {
            if(segment->score >= Y) {
                struct PairwiseAlignment *pA2 = segment_getSubAlignment(segment, pA);
                cigarWrite(output, pA2, 0);
                destructPairwiseAlignment(pA2);
            }
        }
        else {
            segment_reportXFullSubalignments(segment->child1, pA, X, Y);
            segment_reportXFullSubalignments(segment->child2, pA, X, Y);
            segment_reportXFullSubalignments(segment->child3, pA, X, Y);
        }
    }
}

/*
 * Following gives lastz scoring values.
 */

int64_t lastZGapOpenPenalty = -400;
int64_t lastZGapExtendPenalty = -30;
int64_t index(char a) {
    assert(a != '-');
    switch(toupper(a)) {
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return 4; //treated as wildcard
    }
}

int64_t lastZScoreFn(char a, char b) {
    /*
     * Scores a pair of residues.
     */
    int64_t *scoringMatrix = {
            /*A-A*/ 91, /*A-C*/, -114, /*A-G*/ -31, /*A-T*/ -123, /*A - N*/ -44,
            /*C-A*/ -114, /*C-C*/, 100, /*C-G*/ -125, /*C-T*/ -31, /*C - N*/ -44,
            /*G-A*/ -31, /*G-C*/, -125, /*G-G*/ 100, /*G-T*/ -114, /*G - N*/ -44,
            /*T-A*/ -123, /*T-C*/, -31, /*T-G*/ -114, /*T-T*/ -91, /*T - N*/ -44,
            /*N-A*/ -44, /*N-C*/, -44, /*N-G*/ -44, /*N-T*/ -44, /*N - N*/ -44 };
    return scoringMatrix[index(a)*5 + index(b)];
}

int64_t splitPoint(const char *seqX,
        const char *seqY, int64_t length, int64_t (*matchScoreFn)(char, char), int64_t *score) {
    *score = 0;
    for(int64_t i=0; i<length; i++) {
        int64_t j = matchScoreFn(seqX[i], seqY[i]);
        if((score < 0 && j > 0) || (score > 0 && j < 0)) {
            return i;
        }
        *score += j;
    }
    return length;
}

Segment *split(struct PairwiseAlignment *pA,
        int64_t (*matchScoreFn)(char, char), int64_t gapOpenPenalty,
        int64_t gapExtendPenalty, const char *seqX, const char *seqY) {
    /*
     * Returns a segment representing a maximal prefix of the pA with either positive or negative score,
     * removing the prefix from pA in the process.
     */
    assert(pA->strand1);
    assert(pA->strand2);
    Segment *segment = st_calloc(1, sizeof(Segment));
    struct List *opList = constructEmptyList(0, (void (*)(void))destructAlignmentOperation);
    int64_t lX = 0, lY = 0;
    for(int64_t i=0; i<pA->operationList->length; i++) {
        struct AlignmentOperation *op = pA->operationList->list[i];
        if(op->opType == PAIRWISE_MATCH) {
            int64_t k;
            int64_t m = splitPoint(seqX, pA->strand1 ? lX + pA->start1 : pA->start1 - lX, pA->strand1.
                    seqY, pA->strand2 ? lY + pA->start2 : pA->start2 - lY, pA->strand2,
                    op->length, matchScoreFn, &k);
            assert(k != 0);
            if((segment->score < 0 && k > 0) || (segment->score > 0 && k < 0)) { //If the match has the opposite sign to the existing match then we're good.
                break;
            }
            segment->score += k;
            lX += m;
            lY += m;
            if(m < op->length) {
                op->length -= m;
                List_append(opList, constructAlignmentOperation(PAIRWISE_MATCH, m, 0));
                break;
            }
            else {
                List_append(opList, List_remove(pA, op)); //This could be improved it if proves slow.
            }
        }
        else {
            if(segment->score > 0) {
                break;
            }
            segment->score += gapOpenPenalty + (segment->length-1) * gapExtendPenalty;
            List_append(opList, List_remove(pA, op)); //This could be improved it if proves slow.
            if(op->opType == PAIRWISE_INDEL_X) {
                lX += op->length;
            }
            else {
                assert(op->opType == PAIRWISE_INDEL_Y);
                lY += op->length;
            }
        }
    }
    int64_t end1 = pA->start1 + (pA->strand1 ? lX : -lX);
    int64_t end2 = pA->start2 + (pA->strand2 ? lY : -lY);
    segment->pA = constructPairwiseAlignment(pA->contig1, pA->start1, end1, pA->strand1,
            pA->contig2, pA->start2, end2, pA->strand2, 0, opList);
    pA->start1 = end1;
    pA->start2 = end2;
    segment->minScore = segment->score;
    return segment;
}

stList *partitionCigarIntoUsefulTrees(struct PairwiseAlignment *pA,
                                      int64_t (*scoreFn)(char, char),
                                      int64_t gapOpenPenalty, int64_t gapExtendPenalty,
                                      const char *seqX, const char *seqY) {
    /*
     * Builds the useful tree for a given alignment.
     */
    //Partitions an alignment into a sequence of "segments" with alternating signed scores.
    stList *segments = stList_construct();
    while(pA->operationList->length > 0) {
        stList_append(segments, split(pA, seqX, seqY, matchScoreFn, gapOpenPenalty, gapExtendPenalty));
    }

    //Now perform algorithm described in paper to construct useful trees.
}

static stHash *hashAdditionalSequences_sequences;
static void addSeq(const char *header, const char *sequence, int64_t length) {
    stHash_insert(hashAdditionalSequences_sequences, header, sequence);
}

void hashAdditionalSequences(const char *seqFile, stHash *sequences) {
    /*
     * Gets sequences from fasta and hashes them by name, like hashSequences, but add to the "sequences" hash.
     */
    FILE *fH = fopen(seqfile, "r");
    hashAdditionalSequences_sequences = sequences;
    fastaReadToFunction(fH, addSeq);
}

stHash *hashSequences(const char *seqFile) {
    /*
     * Gets sequences from fasta and hashes them by name.
     */
    stHash *sequences = stHash_construct3((uint64_t (*)(const void *))stHash_stringKey,
            (int (*)(const void *, const void *))stHash_stringEqualKey, free, free);
    hashAdditionalSequences(seqFile, sequences);
    return sequences;
}



int main(int argc, char *argv[]) {
    //arguments are log-string, X, Y, inputCigarFile, outputCigarFile, seqFilesxN (where N >= 1)
    assert(argc >= 7);

    //set log level
    st_setLogLevelFromString(argv[1]);

    //Set X
    int64_t X;
    int64_t i = sscanf(argv[2], "%" PRIi64 "", &X);
    assert(i == 1);

    //Set Y
    int64_t X;
    i = sscanf(argv[3], "%" PRIi64 "", &Y);
    assert(i == 1);

    //Get the file-handles for the cigars
    FILE *inputCigarFileHandle = fopen(argv[4], "r"), *outputCigarFileHandle = fopen(argv[5], "w");

    //Get sequences hashed by fasta header
    stHash *sequences = hashSequences(argv[6]);
    for(i=6; i<argc; i++) {
        hashAdditionalSequences(argv[i], sequences);
    }
    struct PairwiseAlignment *pA;
    while((pA = cigarRead(inputCigarFileHandle)) != NULL) {
        stList *usefulTrees = partitionCigarIntoUsefulTrees(pA, lastZScoreFn, stHash_search(sequences, pA->contig1), stHash_search(sequences, pA->contig2));
        for(int64_t i=0; i<stList_length(usefulTrees); i++) {
            segment_reportXFullSubalignments(stList_get(usefulTrees, i), pA, X, outputCigarFileHandle);
        }
        stList_destruct(usefulTrees);
        destructPairwiseAlignment(pA);
    }

    //Cleanup
    fclose(inputCigarFileHandle);
    fclose(outputCigarFileHandle);
    stHash_destruct(sequences);

    return 0;
}
