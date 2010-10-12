#include "endAligner.h"
#include "cactus.h"
#include "sonLib.h"
#include "adjacencySequences.h"

/*
 * Gets an ordered list of pairs from the end alignment for the given adjacency sequence.
 */
stList *getInducedAlignment(stSortedSet *endAlignment, AdjacencySequence *adjacencySequence) {
    stList *inducedAlignment = stList_construct();
    if(adjacencySequence->strand) {
        AlignedPair *alignedPair = alignedPair_construct(adjacencySequence->sequenceName,
                    adjacencySequence->start-1, adjacencySequence->strand, 0, 0, 0, 0);
        AlignedPair *alignedPair2 = stSortedSet_searchGreaterThan(endAlignment, alignedPair);
        if (alignedPair2 != NULL) {
           stSortedSetIterator *it = stSortedSet_getIteratorFrom(endAlignment, alignedPair2);
           while ((alignedPair2 = stSortedSet_getNext(it)) != NULL) {
               if (alignedPair2->sequence == adjacencySequence->sequenceName) {
                   if (alignedPair2->position >= adjacencySequence->start + adjacencySequence->length) {
                       break;
                   }
                   if (alignedPair2->strand && alignedPair2->position >= adjacencySequence->start) {
                       stList_append(inducedAlignment, alignedPair2);
                   }
               }
               else {
                   break;
               }
           }
           stSortedSet_destructIterator(it);
        }
        alignedPair_destruct(alignedPair->reverse);
        alignedPair_destruct(alignedPair);
    }
    else {
        AlignedPair *alignedPair = alignedPair_construct(adjacencySequence->sequenceName,
                            adjacencySequence->start+1, adjacencySequence->strand, 0, 0, 0, 0);
        AlignedPair *alignedPair2 = stSortedSet_searchLessThan(endAlignment, alignedPair);
        if (alignedPair2 != NULL) {
            stSortedSetIterator *it = stSortedSet_getIteratorFrom(endAlignment, alignedPair2);
            stSortedSet_getNext(it);
            stSortedSet_getNext(it); //shift it up..
            while ((alignedPair2 = stSortedSet_getPrevious(it)) != NULL) {
                if (alignedPair2->sequence == adjacencySequence->sequenceName) {
                    if (alignedPair2->position <= adjacencySequence->start - adjacencySequence->length) {
                        break;
                    }
                    if (!alignedPair2->strand && alignedPair2->position <= adjacencySequence->start) {
                        stList_append(inducedAlignment, alignedPair2);
                    }
                }
                else {
                    break;
                }
            }
            stSortedSet_destructIterator(it);
        }
        alignedPair_destruct(alignedPair->reverse);
        alignedPair_destruct(alignedPair);
    }
    return inducedAlignment;
}

/*
 * Runs along and cumulate the score of the pairs, traversing forward through the induced alignment.
 */
static int64_t *cumulateScoreForward(stList *inducedAlignment1) {
    int64_t *iA = st_malloc(sizeof(int64_t) * stList_length(inducedAlignment1));
    int64_t totalScore = 0;
    for(int32_t i=0; i<stList_length(inducedAlignment1); i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment1, i);
        totalScore += alignedPair->score;
        iA[i] = totalScore;
    }
    return iA;
}

/*
 * Runs along and cumulate the score of the pairs, traversing backward through the induced alignment.
 */
static int64_t *cumulateScoreBackward(stList *inducedAlignment1) {
    int64_t *iA = st_malloc(sizeof(int64_t) * stList_length(inducedAlignment1));
    int64_t totalScore = 0;
    for(int32_t i=stList_length(inducedAlignment1)-1; i>=0; i--) {
        AlignedPair *alignedPair = stList_get(inducedAlignment1, i);
        totalScore += alignedPair->score;
        iA[i] = totalScore;
    }
    return iA;
}

/*
 * Chooses a point along the adjacency sequence at which to filter the two alignments,
 * then filters the aligned pairs by this point.
 */
static void pruneAlignments2(stList *inducedAlignment1, stList *inducedAlignment2,
        stSortedSet *endAlignment1, stSortedSet *endAlignment2) {
    int64_t *cScore1 = cumulateScoreForward(inducedAlignment1);
    int64_t *cScore2 = cumulateScoreBackward(inducedAlignment2);

#ifdef BEN_DEBUG //Check the score arrays for sanity..
    for(int32_t i=1; i<stList_length(inducedAlignment1); i++) {
        assert(cScore1[i-1] < cScore1[i]);
    }
    for(int32_t i=1; i<stList_length(inducedAlignment2); i++) {
        assert(cScore2[i-1] > cScore2[i]);
    }
#endif

    //Find the score threshold by iterating along each alignment.
    int32_t cutOff1 = 0, cutOff2 = 0;
    int64_t maxScore = -1;
    if(stList_length(inducedAlignment2) > 0) {
        maxScore = cScore2[0];
    }
    int32_t j=0;
    for(int32_t i=0; i<stList_length(inducedAlignment1); i++) {
        AlignedPair *alignedPair1 = stList_get(inducedAlignment1, i);
        assert(alignedPair1->strand);
        if(j<stList_length(inducedAlignment2)) {
            do {
                AlignedPair *alignedPair2 = stList_get(inducedAlignment2, j);
#ifdef BEN_DEBUG
                assert(!alignedPair2->strand);
#endif
                if(alignedPair1->position < alignedPair2->position) {
                    if(cScore1[i] + cScore2[j] >= maxScore) {
                        maxScore = cScore1[i] + cScore2[j];
                        cutOff1 = i+1;
                        cutOff2 = j;
                    }
                    break;
                }
                else {
                    j++;
                }
            } while(j<stList_length(inducedAlignment2));
        }
        else {
            if(cScore1[i] >= maxScore) {
                cutOff1 = stList_length(inducedAlignment1);
                cutOff2 = j;
                break;
            }
        }
    }

    //Now do the actual filtering of the alignments.
    for(int32_t i=cutOff1; i<stList_length(inducedAlignment1); i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment1, i);
#ifdef BEN_DEBUG
        assert(stSortedSet_search(endAlignment1, alignedPair) != NULL);
        assert(stSortedSet_search(endAlignment1, alignedPair->reverse) != NULL);
#endif
        stSortedSet_remove(endAlignment1, alignedPair);
        stSortedSet_remove(endAlignment1, alignedPair->reverse);
        alignedPair_destruct(alignedPair->reverse);
        alignedPair_destruct(alignedPair);
    }

    for(int32_t i=0; i<cutOff2; i++) {
        AlignedPair *alignedPair = stList_get(inducedAlignment2, i);
#ifdef BEN_DEBUG
        assert(stSortedSet_search(endAlignment2, alignedPair) != NULL);
        assert(stSortedSet_search(endAlignment2, alignedPair->reverse) != NULL);
#endif
        stSortedSet_remove(endAlignment2, alignedPair);
        stSortedSet_remove(endAlignment2, alignedPair->reverse);
        alignedPair_destruct(alignedPair->reverse);
        alignedPair_destruct(alignedPair);
    }

    //Cleanup
    free(cScore1);
    free(cScore2);
}

/*
 * Removes pairs in the alignment so that each position in the adjacency sequence is a member of one or other alignment.
 */
static void pruneAlignments(Cap *cap, stSortedSet *endAlignment1, Cap *adjacentCap, stSortedSet *endAlignment2) {
    if (endAlignment1 == endAlignment2) {
        return; //We ignore self loops, which have been dealt with already by the
        //poset aligner.
    }
    AdjacencySequence *adjacencySequence1 = adjacencySequence_construct(cap, INT32_MAX);
    AdjacencySequence *adjacencySequence2 = adjacencySequence_construct(adjacentCap, INT32_MAX);
#ifdef BEN_DEBUG
    assert(adjacencySequence1->length == adjacencySequence2->length);
    assert(adjacencySequence1->sequenceName == adjacencySequence2->sequenceName);
    assert(adjacencySequence1->strand == !adjacencySequence2->strand);
    assert(adjacencySequence2->start == adjacencySequence1->start + adjacencySequence1->length - 1);
#endif

    stList *inducedAlignment1 = getInducedAlignment(endAlignment1, adjacencySequence1);
    stList *inducedAlignment2 = getInducedAlignment(endAlignment2, adjacencySequence2);
    stList_reverse(inducedAlignment2);

    pruneAlignments2(inducedAlignment1, inducedAlignment2, endAlignment1, endAlignment2);

    //Cleanup.
    adjacencySequence_destruct(adjacencySequence1);
    adjacencySequence_destruct(adjacencySequence2);
    stList_destruct(inducedAlignment1);
    stList_destruct(inducedAlignment2);
}

stSortedSet *makeFlowerAlignment(Flower *flower, int32_t spanningTrees,
        int32_t maxSequenceLength, float gapGamma, bool useBanding,
        int32_t bandingSize, float bandingThreshold) {
    //Make the end alignments, representing each as an adjacency alignment.
    End *end;
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    stHash *endAlignments = stHash_construct2(NULL,
            (void(*)(void *)) stSortedSet_destruct);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        stSortedSet *endAlignment = makeEndAlignment(end, spanningTrees, maxSequenceLength, gapGamma, useBanding,
                bandingSize, bandingThreshold);
        stHash_insert(endAlignments, end, endAlignment);
    }
    flower_destructEndIterator(endIterator);

    //Prune end alignments
    endIterator = flower_getEndIterator(flower);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        stSortedSet *endAlignment1 = stHash_search(endAlignments, end);
        assert(endAlignment1 != NULL);
        Cap *cap;
        End_InstanceIterator *capIterator = end_getInstanceIterator(end);
        while ((cap = end_getNext(capIterator)) != NULL) {
            if(cap_getSide(cap)) {
                cap = cap_getReverse(cap);
            }
            if(cap_getStrand(cap)) {
                Cap *adjacentCap = cap_getAdjacency(cap);
#ifdef BEN_DEBUG
                assert(adjacentCap != NULL);
                assert(cap_getSide(adjacentCap));
                assert(cap_getStrand(adjacentCap));
#endif
                adjacentCap = cap_getReverse(adjacentCap);
                stSortedSet *endAlignment2 = stHash_search(endAlignments, cap_getEnd(adjacentCap));
                if(endAlignment2 == NULL) {
                    endAlignment2 = stHash_search(endAlignments, cap_getEnd(cap_getReverse(adjacentCap)));
                }
#ifdef BEN_DEBUG
                assert(endAlignment2 != NULL);
#endif
                //Now traverse the alignments, pruning stuff
                pruneAlignments(cap, endAlignment1, adjacentCap, endAlignment2);
            }
        }
        end_destructInstanceIterator(capIterator);
    }
    flower_destructEndIterator(endIterator);

    //Now convert to set of final aligned pairs to return.
    stSortedSet *sortedAlignment = stSortedSet_construct3((int (*)(const void *, const void *))alignedPair_cmpFn,
                (void (*)(void *))alignedPair_destruct);
    stList *endAlignmentsList = stHash_getValues(endAlignments);
    while (stList_length(endAlignmentsList) > 0) {
        stSortedSet *endAlignment = stList_pop(endAlignmentsList);
        while(stSortedSet_size(endAlignment) > 0) {
            AlignedPair *alignedPair = stSortedSet_getFirst(endAlignment);
            stSortedSet_remove(endAlignment, alignedPair);
#ifdef BEN_DEBUG
            assert(stSortedSet_search(endAlignment, alignedPair->reverse) != NULL);
#endif
            stSortedSet_remove(endAlignment, alignedPair->reverse);
#ifdef BEN_DEBUG
            assert(stSortedSet_search(sortedAlignment, alignedPair) == NULL);
            assert(stSortedSet_search(sortedAlignment, alignedPair->reverse) == NULL);
#endif
            stSortedSet_insert(sortedAlignment, alignedPair);
            stSortedSet_insert(sortedAlignment, alignedPair->reverse);
        }
    }
    stList_destruct(endAlignmentsList);
    stHash_destruct(endAlignments);

    return sortedAlignment;
}

