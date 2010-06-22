
#include "sonLib.h"
#include "stPosetAlignment.h"

struct _stPosetAlignment {
    stList sequences;
    stSortedSet *pairAlignments;
    stHash *constraintLists;
};

stPosetAlignment *stPosetAlignment_construct(void) {
    stPosetAlignment *posetAlignment = st_malloc(sizeof(stPosetAlignment));
    posetAlignment->sequences = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    posetAlignment->pairAlignments =
        stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn,
            (void (*)(void *))stIntTuple_destruct);
    posetAlignment->constraintLists = stHash_construct3((uint32_t (*)(const void *))stIntTuple_hashKey,
            (int (*)(const void *, const void *))stIntTuple_equalsFn, (void (*)(void *))stIntTuple_destruct,
            (void (*)(void *))stList_destruct);
    return posetAlignment;
}

void stPosetAlignment_destruct(stPosetAlignment *posetAlignment) {
    stList_destruct(posetAlignment->sequences);
    stSortedSet_destruct(posetAlignment->pairAlignments);
    stHash_destruct(posetAlignment->constraintLists);
    free(posetAlignment);
}

void stPosetAlignment_addSequence(stPosetAlignment *posetAlignment, int32_t length) {
   stList_append(posetAlignment->sequences, stIntTuple_construct(1, length);
}

int32_t stPosetAlignment_getSequenceNumber(stPosetAlignment *posetAlignment) {
    return stList_length(posetAlignment->sequences);
}

int32_t stPosetAlignment_getSequenceLength(stPosetAlignment *posetAlignment, int32_t sequence) {
    stIntTuple *i = stList_get(posetAlignment->sequences, sequence);
    return stIntTuple_getPosition(i, 1);
}

void stPosetAlignment_addPairwiseAlignment(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t sequence2) {
    if(sequence1 > sequence2) {
        int32_t i = sequence1;
        sequence1 = sequence2;
        sequence2 = sequence1;
    }
    stIntTuple *pair = stIntTuple_construct(2, sequence1, sequence2);
    if(stSortedSet_search(i) == NULL) {
        stSortedSet_insert(posetAlignment->pairAlignments, pair);
        stHash_insert(posetAlignment->constraintLists, stIntTuple_construct(2, sequence1, sequence2),
        stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct));
        stHash_insert(posetAlignment->constraintLists, stIntTuple_construct(2, sequence2, sequence1),
        stSortedSet_construct3((int (*)(const void *, const void *))stIntTuple_cmpFn, (void (*)(void *))stIntTuple_destruct));
    }
    else {
        stIntTuple_destruct(i);
    }
}

int32_t stPosetAlignment_getNumberPairwiseAlignments(stPosetAlignment *posetAlignment) {
    return stSortedSet_size(posetAlignment->pairAlignments);
}

stSortedSet *stPosetAlignment_getPairwiseAlignments(stPosetAlignment *posetAlignment) {
    return stSortedSet_copyConstruct(posetAlignment->pairAlignments);
}

static int32_t getConstraint(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2) {
    stIntTuple *pair = stIntTuple_construct(2, sequence1, sequence2);
    assert(stSortedSet_search(posetAlignment->pairAlignments, pair) != NULL);
    stSortedSet *constraintList = stHash_search(posetAlignment->constraintLists, pair);
    stIntTuple_destruct(pair);
    stIntTuple *pos = stIntTuple_construct(2, position1, INT32_MAX);
    //Get less than or equal
    return NULL;
}

bool stPosetAlignment_isPossible(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {
    return position2 <= getConstraint(posetAlignment, sequence1, sequence2, position1) && position2 <= getConstraint(posetAlignment, sequence2, sequence1, position2);
}

bool stPosetAlignment_add(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {
    if(stPosetAlignment_isPossible(posetAlignment, sequence1, position1, sequence2, position2)) {
        //for all pairs do check..
        return 1;
    }
    return 0;
}
