#include "multipleAligner.h"
#include "sonLib.h"
#include "stPosetAlignment.h"
#include "pairwiseAligner.h"

/*
 * Functions to align a bunch of sequences, creating a global alignment.
 */

/*
 * Constructs a random spanning tree linking all the nodes in items into one component.
 */
void constructSpanningTree(int32_t numberOfSequences,
        stSortedSet *pairwiseAlignments) {
    stList *list = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    for (int32_t i = 0; i < numberOfSequences; i++) {
        stList_append(list, stIntTuple_construct(1, i));
    }
    while (stList_length(list) > 1) {
        stIntTuple *i = st_randomChoice(list);
        stList_removeItem(list, i);
        int32_t j = stIntTuple_getPosition(i, 0);
        int32_t k = stIntTuple_getPosition(st_randomChoice(list), 0);
        assert(j != k);

        stIntTuple *l = j < k ? stIntTuple_construct(2, j, k) : stIntTuple_construct(2, k, j);
        stIntTuple_destruct(i);
        if (!stSortedSet_search(pairwiseAlignments, l)) {
            stSortedSet_insert(pairwiseAlignments, l);
        } else {
            stIntTuple_destruct(l);
        }
    }
    stList_destruct(list);
}

stList *makeAlignment(stList *sequences,
        int32_t spanningTrees, void *modelParameters) {
    //Get the set of pairwise alignments (by constructing spanning trees)
    stSortedSet *pairwiseAlignments = stSortedSet_construct3((int(*)(
            const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    //for (int i = 0; i < spanningTrees; i++) {
    //    constructSpanningTree(stList_length(sequences), pairwiseAlignments);
    //}

    for(int32_t i=0; i<stList_length(sequences); i++) {
        for(int32_t j=i+1; j<stList_length(sequences); j++) {
            stSortedSet_insert(pairwiseAlignments, stIntTuple_construct(2, i, j));
        }
    }

    //Construct the alignments
    //and sort them by weight
    stSortedSetIterator *pairwiseAlignmentsIterator = stSortedSet_getIterator(
            pairwiseAlignments);
    stIntTuple *pairwiseAlignment;
    stList *alignedPairs = stList_construct();
    while ((pairwiseAlignment = (stIntTuple *)stSortedSet_getNext(pairwiseAlignmentsIterator))
            != NULL) {
        int32_t sequence1 = stIntTuple_getPosition(pairwiseAlignment, 0);
        int32_t sequence2 = stIntTuple_getPosition(pairwiseAlignment, 1);
        stList *alignedPairs2 = getAlignedPairs(stList_get(sequences, sequence1), stList_get(sequences, sequence2), modelParameters);
        while(stList_length(alignedPairs2) > 0) {
            stIntTuple *alignedPair = (stIntTuple *)stList_pop(alignedPairs2);
#ifdef BEN_DEBUG
            assert(stIntTuple_length(alignedPair) == 3);
#endif
            stList_append(alignedPairs, stIntTuple_construct(5,
                    /* score */ stIntTuple_getPosition(alignedPair, 0),
                    /*seq 1 */ sequence1, stIntTuple_getPosition(alignedPair, 1),
                    /*seq 2 */ sequence2, stIntTuple_getPosition(alignedPair, 2)));
            stIntTuple_destruct(alignedPair);
        }
        stList_destruct(alignedPairs2);
    }
    stSortedSet_destructIterator(pairwiseAlignmentsIterator);

    //Sort the pairs (by weight)
    stList_sort(alignedPairs, (int (*)(const void *, const void *))stIntTuple_cmpFn);

    //Greedily construct poset and filter pairs..
    stPosetAlignment *posetAlignment = stPosetAlignment_construct(stList_length(sequences));
    stList *acceptedAlignedPairs = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    int32_t pScore = INT32_MAX;
    while (stList_length(alignedPairs) > 0) {
        stIntTuple *alignedPair = stList_pop(alignedPairs);
        int32_t score = stIntTuple_getPosition(alignedPair, 0);
#ifdef BEN_DEBUG
        assert(score > 0);
        assert(score <= pScore);
#endif
        pScore = score;
        int32_t sequence1 = stIntTuple_getPosition(alignedPair, 1);
        int32_t position1 = stIntTuple_getPosition(alignedPair, 2);
        int32_t sequence2 = stIntTuple_getPosition(alignedPair, 3);
        int32_t position2 = stIntTuple_getPosition(alignedPair, 4);
        if (stPosetAlignment_isPossible(posetAlignment, sequence1, position1, sequence2, position2)) {
            stPosetAlignment_add(posetAlignment, sequence1, position1, sequence2, position2);
            //Add a converted version to the accepted aligned pairs.
            stList_append(acceptedAlignedPairs, alignedPair);
        }
        else {
            stIntTuple_destruct(alignedPair);
        }
    }
    stList_destruct(alignedPairs);
    stPosetAlignment_destruct(posetAlignment);
    stSortedSet_destruct(pairwiseAlignments);

    //Return the accepted pairs
    return acceptedAlignedPairs;
}
