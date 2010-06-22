
#include "sonLib.h"
#include "stPosetAlignment.h"

struct _stPosetAlignment {
    stList sequences
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
            (int (*)(const void *, const void *))stIntTuple_equalsFn, NULL, (void (*)(void *))stSortedSet_destruct);
    return posetAlignment;
}

void stPosetAlignment_destruct(stPosetAlignment *posetAlignment) {
    stSortedSet_destruct(posetAlignment->pairAlignments);
    stHash_destruct(posetAlignment->constraintLists);
    free(posetAlignment);
}

void stPosetAlignment_addSequence(stPosetAlignment *posetAlignment, int32_t length) {
   stList_append(posetAlignment->sequences, stIntTuple_construct(1, length);
}

int32_t stPosetAlignment_getSequenceNumber(stPosetAlignment *posetAlignment) {

}

int32_t stPosetAlignment_getSequenceLength(stPosetAlignment *posetAlignment, int32_t sequence) {

}

void stPosetAlignment_addPairwiseAlignment(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t sequence2) {

}

int32_t stPosetAlignment_getNumberPairwiseAlignments(stPosetAlignment *posetAlignment) {

}

stSortedSet *stPosetAlignment_getPairwiseAlignments(stPosetAlignment *posetAlignment) {

}

bool stPosetAlignment_isPossible(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {

}

bool stPosetAlignment_add(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {

}
