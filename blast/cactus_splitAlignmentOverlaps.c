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

uint64_t getEndCoordinate(struct PairwiseAlignment *pairwiseAlignment) {
    assert(pairwiseAlignment->strand1); // This code assumes that the alignment is reported with respect
    // to the positive strand of the first sequence
    return pairwiseAlignment->end1;
}

static int64_t max(int64_t i, int64_t j) {
    if (i < j) {
        return j;
    }
    return i;
}

static int64_t min(int64_t i, int64_t j) {
    if (i < j) {
        return i;
    }
    return j;
}

int comparePairwiseAlignmentsByScore(const void *a, const void *b) {
    struct PairwiseAlignment *pA1 = (struct PairwiseAlignment *)a;
    struct PairwiseAlignment *pA2 = (struct PairwiseAlignment *)b;

    int i = pA1->score < pA2->score ? 1 : pA1->score > pA2->score ? -1 : 0;
    if(i == 0) {
        i = pA1 > pA2 ? 1 : pA1 < pA2 ? -1 : 0;
    }

    return i;
}

int comparePairwiseAlignmentsByStart(const void *a, const void *b) {
    struct PairwiseAlignment *pA1 = (struct PairwiseAlignment *)a;
    struct PairwiseAlignment *pA2 = (struct PairwiseAlignment *)b;

    int i = strcmp(pA1->contig1, pA2->contig1);
    if(i == 0) {
        i = pA1->start1 > pA2->start1 ? 1 : pA1->start1 < pA2->start1 ? -1 : 0;
        if(i == 0) {
            i = pA1->end1 > pA2->end1 ? 1 : pA1->end1 < pA2->end1 ? -1 : 0;
            if(i == 0) {
                i = pA1 > pA2 ? 1 : pA1 < pA2 ? -1 : 0;
            }
        }
    }

    return i;
}

int comparePairwiseAlignmentsByEnd(const void *a, const void *b) {
    struct PairwiseAlignment *pA1 = (struct PairwiseAlignment *)a;
    struct PairwiseAlignment *pA2 = (struct PairwiseAlignment *)b;

    int i = strcmp(pA1->contig1, pA2->contig1);
    if(i == 0) {
        i = pA1->end1 > pA2->end1 ? 1 : pA1->end1 < pA2->end1 ? -1 : 0;
        if(i == 0) {
            i = pA1 > pA2 ? 1 : pA1 < pA2 ? -1 : 0;
        }
    }

    return i;
}

struct PairwiseAlignment *removeAlignmentPrefix(struct PairwiseAlignment *pairwiseAlignment, int64_t prefixEnd) {
    // Store the original start coordinates
    int64_t start1 = pairwiseAlignment->start1, start2 = pairwiseAlignment->start2;
    assert(pairwiseAlignment->end1 > prefixEnd);
    assert(pairwiseAlignment->start1 < prefixEnd);
    assert(pairwiseAlignment->strand1);

    // Split the ops in the cigar string between the prefix and suffix alignments
    struct List *prefixOps = constructEmptyList(0, (void (*)(void *))destructAlignmentOperation);
    int64_t i = 0;
    do {
        assert(i < pairwiseAlignment->operationList->length);
        struct AlignmentOperation *op = pairwiseAlignment->operationList->list[i];
        assert(op->length > 0);

        if(op->opType == PAIRWISE_INDEL_Y) { // Insert in second sequence
            listAppend(prefixOps, op);
            i++;
            pairwiseAlignment->start2 += pairwiseAlignment->strand2 ? op->length : -op->length;
        }
        else { // Not an insert in second sequence
            // Op is in the prefix alignment
            int64_t j;
            if(pairwiseAlignment->start1 + op->length <= prefixEnd) {
                listAppend(prefixOps, op);
                i++;
                j = op->length;
            }
            // Op spans the prefix and suffix alignments, so split it
            else {
                j = prefixEnd-pairwiseAlignment->start1;
                assert(j > 0);
                listAppend(prefixOps, constructAlignmentOperation(op->opType, j, op->score));
                op->length -= j;
                assert(op->length > 0);
                assert(pairwiseAlignment->start1+j == prefixEnd);
            }

            // Update start coordinates of suffix alignments
            pairwiseAlignment->start1 += j;
            if(op->opType != PAIRWISE_INDEL_X) {
                pairwiseAlignment->start2 += pairwiseAlignment->strand2 ? j : -j;
            }
        }
    } while(pairwiseAlignment->start1 < prefixEnd);

    assert(pairwiseAlignment->start1 == prefixEnd);

    if (i > 0) {
        // Remove prefix ops from pairwise alignment
        uint64_t j=0;

        while(i < pairwiseAlignment->operationList->length) {
            pairwiseAlignment->operationList->list[j++] = pairwiseAlignment->operationList->list[i++];
        }
        pairwiseAlignment->operationList->length = j;
    }
    assert(pairwiseAlignment->operationList->length > 0);

    // Create prefix pairwiseAlignment
    struct PairwiseAlignment *prefixAlignment = constructPairwiseAlignment(pairwiseAlignment->contig1,
                                                                           start1, pairwiseAlignment->start1, 1,
                                                                           pairwiseAlignment->contig2, start2, pairwiseAlignment->start2, pairwiseAlignment->strand2,
                                                                           pairwiseAlignment->score, prefixOps);

    return prefixAlignment;
}

bool covers(struct PairwiseAlignment *pA, int64_t from, int64_t to, int64_t *coverStart, int64_t *coverEnd) {
    if (from == to) {
        st_errAbort("here");
    }
    if (pA->start1 == pA->end1) {
        return false;
    }
    if (pA->start1 >= from && pA->end1 <= to) {
        // The idea here is, if the entire alignment is within the
        // boundaries, we probably don't need to worry about keeping
        // accurate track of the depth within its indels. This is so
        // that we don't split *every* alignment into only its match
        // pieces.
        *coverStart = pA->start1;
        *coverEnd = pA->end1;
        return true;
    }
    int64_t start = pA->start1;
    for (int64_t i = 0; i < pA->operationList->length; i++) {
        struct AlignmentOperation *op = pA->operationList->list[i];
        if (op->opType == PAIRWISE_MATCH) {
            if (from >= start && from < start + op->length) {
                st_logDebug("bar");
                *coverStart = max(start, from);
                *coverEnd = min(start + op->length, to);
                return true;
            }
            start += op->length;
        } else if (op->opType == PAIRWISE_INDEL_X) {
            start += op->length;
        }
        if (start >= to) {
            break;
        }
    }
    return false;
}

void getNonOverlappingSubset(stSortedSet *activeAlignments, uint64_t from, uint64_t to, uint64_t maxDepth, struct PairwiseAlignment *startAlignment, stList *output) {
    struct PairwiseAlignment *pairwiseAlignment;
    stSortedSetIterator *it;
    if (startAlignment == NULL) {
        it = stSortedSet_getIterator(activeAlignments);
    } else {
        st_logDebug("in set: %p", stSortedSet_search(activeAlignments, startAlignment));
        it = stSortedSet_getIteratorFrom(activeAlignments, startAlignment);
    }
    st_logDebug("get overlapping pairwise alignments from %" PRIi64 " to %" PRIi64 "\n", from, to);
    while (maxDepth > 0 && (pairwiseAlignment = stSortedSet_getNext(it)) != NULL) {
        int64_t coverageStart, coverageEnd;
        if (covers(pairwiseAlignment, from, to, &coverageStart, &coverageEnd)) {
            st_logDebug("found pairwise alignment covering from %" PRIi64 " to %" PRIi64 "\n", coverageStart, coverageEnd);
            assert(coverageStart >= from);
            assert(coverageEnd <= to);
            if (pairwiseAlignment->start1 < from) {
                // Discard the section of the alignment before we're interested in.
                destructPairwiseAlignment(removeAlignmentPrefix(pairwiseAlignment, from));
            }
            if (coverageStart == from && coverageEnd == to) {
                // Easy case: the alignment covers the entire from..to range.
                if (pairwiseAlignment->end1 > to) {
                    // Save the suffix of the alignment after we're
                    // interested in, and restrict ourselves to the part
                    // that covers our region of interest.
                    pairwiseAlignment = removeAlignmentPrefix(pairwiseAlignment, to);
                }
                stList_append(output, pairwiseAlignment);
                maxDepth--;
            } else {
                // Annoying case. Recurse.
                struct PairwiseAlignment *nextPairwiseAlignment = stSortedSet_getNext(it);
                st_logDebug("in set: %p", stSortedSet_search(activeAlignments, nextPairwiseAlignment));
                do {
                    st_logDebug("handling partial coverage from %" PRIi64 " to %" PRIi64 "\n", coverageStart, coverageEnd);
                    stList_append(output, removeAlignmentPrefix(pairwiseAlignment, coverageEnd));
                    if (nextPairwiseAlignment != NULL) {
                        st_logDebug("in set a: %p", stSortedSet_search(activeAlignments, nextPairwiseAlignment));
                        getNonOverlappingSubset(activeAlignments, coverageStart, coverageEnd, maxDepth - 1, nextPairwiseAlignment, output);
                        if (from != coverageStart) {
                            st_logDebug("in set b: %p", stSortedSet_search(activeAlignments, nextPairwiseAlignment));
                            getNonOverlappingSubset(activeAlignments, from, coverageStart, maxDepth, nextPairwiseAlignment, output);
                        }
                    }
                    from = coverageEnd;
                } while (from != to && covers(pairwiseAlignment, from, to, &coverageStart, &coverageEnd));
                if (nextPairwiseAlignment != NULL && to != coverageEnd) {
                    st_logDebug("in set c: %p", stSortedSet_search(activeAlignments, nextPairwiseAlignment));
                    getNonOverlappingSubset(activeAlignments, coverageEnd, to, maxDepth, nextPairwiseAlignment, output);
                }

                break;
            }
        }
    }
    stSortedSet_destructIterator(it);
}

void splitAlignmentOverlaps(stSortedSet *activeAlignments, uint64_t splitFrom, uint64_t splitUpTo, FILE *fileHandleOut, uint64_t maxDepth) {
    if (stSortedSet_size(activeAlignments) == 0 || splitFrom == splitUpTo) {
        return; // Nothing to do
    }

    stList *alignmentsToOutput = stList_construct();
    getNonOverlappingSubset(activeAlignments, splitFrom, splitUpTo, maxDepth, NULL, alignmentsToOutput);
    stList_sort(alignmentsToOutput, comparePairwiseAlignmentsByStart);
    st_logDebug("outputting\n");
    for (int64_t i = 0; i < stList_length(alignmentsToOutput); i++) {
        struct PairwiseAlignment *pA = stList_get(alignmentsToOutput, i);
        cigarWrite(fileHandleOut, pA, false);
    }
}

void cleanupOldAlignments(char *curContig, int64_t curEnd, stSortedSet *activeAlignments, stSortedSet *alignmentsByEnd) {
    if (curContig == NULL) {
        return;
    }
    while (stSortedSet_size(alignmentsByEnd) > 0) {
        assert(stSortedSet_size(activeAlignments) > 0);
        struct PairwiseAlignment *pA = stSortedSet_getFirst(alignmentsByEnd);
        if (strcmp(curContig, pA->contig1) != 0 || pA->end1 > curEnd) {
            break;
        }
        stSortedSet_remove(alignmentsByEnd, pA);
        stSortedSet_remove(activeAlignments, pA);
        destructPairwiseAlignment(pA);
    }
}

int main(int argc, char *argv[]) {
    /*
     * Each alignment has a unique first sequence interval, defined by where it starts and ends on the
     * first sequence.
     * Two alignments partially overlap if their first sequence intervals overlap but are not the same.
     * This program breaks up alignments in the input file so that there are no partial overlaps between
     * alignments, outputting the non-partially-overlapping alignments to the output file.
     */

    FILE *fileHandleIn;
    FILE *fileHandleOut;

    uint64_t maxDepth;

    if (!(argc == 3 || argc == 5)) {
        st_errAbort("Usage: %s logLevel maxDepth [inFile] [outFile]", argv[0]);
    }

    st_setLogLevelFromString(argv[1]);
    if (sscanf(argv[2], "%" PRIu64, &maxDepth) != 1) {
        st_errAbort("Couldn't understand maxDepth");
    }

    if(argc == 3) {
        fileHandleIn = stdin;
        fileHandleOut = stdout;
    }
    else {
        assert(argc == 4);
        fileHandleIn = fopen(argv[3], "r");
        fileHandleOut = fopen(argv[4], "w");
    }

    // Set of alignments being progressively processed, ordered by descending score
    stSortedSet *activeAlignments = stSortedSet_construct3(comparePairwiseAlignmentsByScore, NULL);

    struct PairwiseAlignment *pairwiseAlignment;
    // Alignments being processed, sorted by end position. This is so
    // that we can destruct irrelevant alignments later on.
    stSortedSet *alignmentsByEnd = stSortedSet_construct3(comparePairwiseAlignmentsByEnd, NULL);
    // Contig and end position of non-overlapping alignments already output.
    char *curContig = NULL;
    int64_t curEnd = 0;
    while ((pairwiseAlignment = cigarRead(fileHandleIn)) != NULL) {
        st_logDebug("handling new incoming alignment\n");
        cleanupOldAlignments(curContig, curEnd, activeAlignments, alignmentsByEnd);
    	// There are existing alignments
    	if (stSortedSet_size(activeAlignments) > 0) {
            // If the new alignment is not on the same sequence as the previous alignment
            if(curContig != NULL && strcmp(curContig,
                                           pairwiseAlignment->contig1) != 0) {
                // If pairwiseAlignment is on a new sequence
                splitAlignmentOverlaps(activeAlignments, curEnd, INT64_MAX, fileHandleOut, maxDepth);
                assert(stSortedSet_size(activeAlignments) == 0);
                curContig = stString_copy(pairwiseAlignment->contig1);
                curEnd = 0;
            } else {
                // Remove overlaps in alignments up to but excluding the start of pairwiseAlignment
                int64_t startCoordinate = getStartCoordinate(pairwiseAlignment);
                splitAlignmentOverlaps(activeAlignments, curEnd, startCoordinate, fileHandleOut, maxDepth);
                curEnd = startCoordinate;
            }
    	}

    	// Add pairwiseAlignment to the set of activeAlignemnts
    	stSortedSet_insert(activeAlignments, pairwiseAlignment);
        stSortedSet_insert(alignmentsByEnd, pairwiseAlignment);
    }
    // Remove remaining overlaps in alignments
    splitAlignmentOverlaps(activeAlignments, curEnd, INT64_MAX, fileHandleOut, maxDepth);
    assert(stSortedSet_size(activeAlignments) == 0);

    // Cleanup
    stSortedSet_destruct(activeAlignments);
    free(curContig);
    if(argc == 4) {
    	fclose(fileHandleIn);
    	fclose(fileHandleOut);
    }

    return 0;
}
