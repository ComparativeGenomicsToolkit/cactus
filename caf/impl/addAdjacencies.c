#include "cactus.h"
#include "sonLib.h"

static int addAdjacenciesPP(Cap *cap1, Cap *cap2) {
    assert(cap_getStrand(cap1) && cap_getStrand(cap2));
    Sequence *sequence1 = cap_getSequence(cap1);
    Sequence *sequence2 = cap_getSequence(cap2);
    int64_t i = cactusMisc_nameCompare(sequence_getName(sequence1), sequence_getName(sequence2));
    if (i == 0) {
        int64_t j = cap_getCoordinate(cap1);
        int64_t k = cap_getCoordinate(cap2);
        i = j > k ? 1 : (j < k ? -1 : 0);
        if (i == 0) {
            assert(cap_getSegment(cap1) == cap_getSegment(cap2));
            j = cap_getSide(cap1);
            k = cap_getSide(cap2);
            assert((j && !k) || (!j && k));
            i = j ? -1 : 1;
        }
    }
    return i;
}

void stCaf_addAdjacencies(Flower *flower) {
    //Build a list of caps.
    stList *list = stList_construct();
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            if (!cap_getStrand(cap)) {
                cap = cap_getReverse(cap);
            }
            stList_append(list, cap);
        }
        end_destructInstanceIterator(instanceIterator);
    }
    flower_destructEndIterator(endIterator);
    assert(stList_length(list) % 2 == 0);
    //Sort the list of caps.
    stList_sort(list, (int(*)(const void *, const void *)) addAdjacenciesPP);
    //Now make the adjacencies.
    for (int64_t i = 1; i < stList_length(list); i += 2) {
        Cap *cap = stList_get(list, i - 1);
        Cap *cap2 = stList_get(list, i);
        cap_makeAdjacent(cap, cap2);
    }
    //Clean up.
    stList_destruct(list);
}
