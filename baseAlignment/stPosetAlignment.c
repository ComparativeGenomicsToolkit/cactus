
#include "sonLib.h"
#include "stPosetAlignment.h"
#include "stdlib.h"

struct _stPosetAlignment {
    int32_t sequenceNumber;
    stSortedSet **constraintLists;
};

static int comparePositions(stIntTuple *position1, stIntTuple *position2) {
    if(stIntTuple_getPosition(position1, 0) == INT32_MAX || stIntTuple_getPosition(position2, 0) == INT32_MAX) { //Indicates we should ignore the first position and compare the second.
        assert(stIntTuple_getPosition(position1, 1) != INT32_MAX);
        assert(stIntTuple_getPosition(position2, 1) != INT32_MAX);
        return stIntTuple_getPosition(position1, 1) - stIntTuple_getPosition(position1, 1);
    }
    return stIntTuple_getPosition(position1, 0) - stIntTuple_getPosition(position1, 0);
}

stPosetAlignment *stPosetAlignment_construct(int32_t sequenceNumber) {
    stPosetAlignment *posetAlignment = st_malloc(sizeof(stPosetAlignment));
    posetAlignment->sequenceNumber = sequenceNumber;
    posetAlignment->constraintLists = st_malloc(sizeof(stSortedSet *) * sequenceNumber * sequenceNumber);
    for(int32_t i=0; i<posetAlignment->sequenceNumber; i++) {
        for(int32_t j=0; j<posetAlignment->sequenceNumber; j++) {
            if(i != j) {
                posetAlignment->constraintLists[i*posetAlignment->sequenceNumber + j] =
                    stSortedSet_construct3((int (*)(const void *, const void *))comparePositions,
                            (void (*)(void *))stIntTuple_destruct);
            }
        }
    }
    return posetAlignment;
}

void stPosetAlignment_destruct(stPosetAlignment *posetAlignment) {
    for(int32_t i=0; i<posetAlignment->sequenceNumber; i++) {
        for(int32_t j=0; j<posetAlignment->sequenceNumber; j++) {
            if(i != j) {
                stSortedSet_destruct(posetAlignment->constraintLists[i*posetAlignment->sequenceNumber + j]);
            }
        }
    }
    free(posetAlignment->constraintLists);
    free(posetAlignment);
}

int32_t stPosetAlignment_getSequenceNumber(stPosetAlignment *posetAlignment) {
    return posetAlignment->sequenceNumber;
}

static stSortedSet *getConstraintList(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t sequence2) {
    assert(sequence1 >= 0 && sequence1 < posetAlignment->sequenceNumber);
    assert(sequence2 >= 0 && sequence2 < posetAlignment->sequenceNumber);
    assert(sequence1 != sequence2);
    return posetAlignment->constraintLists[sequence1 * posetAlignment->sequenceNumber + sequence2];
}

/*
 * Gets the position in sequence2 that the position in sequence2 must be less or equal to in the alignment
 */
static int32_t getConstraint_lessThanOrEquals(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2) {
    stIntTuple *pos = stIntTuple_construct(2, position1, INT32_MAX);
    //Get less than or equal
    stIntTuple *constraint = stSortedSet_searchGreaterThanOrEqual(getConstraintList(posetAlignment, sequence1, sequence2), pos);
    int32_t i = INT32_MAX;
    if(constraint != NULL) {
        i = stIntTuple_getPosition(constraint, 1) - (stIntTuple_getPosition(constraint, 0) == position1 ? 0 : 1);
    }
    stIntTuple_destruct(pos);
    return i;
}

/*
 * Gets the position in sequence2 that the position in sequence1 must be greater than or equal to in the alignment
 */
static int32_t getConstraint_greaterThanOrEquals(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2) {
    stIntTuple *pos = stIntTuple_construct(2, INT32_MAX, position1);
    //Get less than or equal
    stIntTuple *constraint = stSortedSet_searchLessThanOrEqual(getConstraintList(posetAlignment, sequence2, sequence1), pos);
    int32_t i = INT32_MAX;
    if(constraint != NULL) {
        i = stIntTuple_getPosition(constraint, 0) + (stIntTuple_getPosition(constraint, 1) == position1 ? 0 : 1);
    }
    stIntTuple_destruct(pos);
    return i;
}

/*
 * Adds a less than or equal constraint to the list of prime constraints, removing any redundant constraints in the process.
 */
void addConstraint_lessThanOrEquals(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {
    stSortedSet *constraintList = getConstraintList(posetAlignment, sequence1, sequence2);
    stIntTuple *pos = stIntTuple_construct(2, position1, position2);
    stIntTuple *pos2;
    while((pos2 = stSortedSet_searchLessThanOrEqual(constraintList, pos)) != NULL) {
        if(stIntTuple_getPosition(pos2, 0) <= position1 && stIntTuple_getPosition(pos2, 1) >= position2) {
            stSortedSet_remove(constraintList, pos2);
            stIntTuple_destruct(pos2);
        }
    }
    stSortedSet_insert(constraintList, pos);
}

bool stPosetAlignment_isPossible(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {
    return position2 <= getConstraint_lessThanOrEquals(posetAlignment, sequence1, sequence2, position1) && position2 <= getConstraint_lessThanOrEquals(posetAlignment, sequence2, sequence1, position2);
}

static void stPosetAlignment_addP2(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t sequence3, int32_t position3, int32_t sequence2, int32_t position2) {
    for(int32_t sequence4=0; sequence4<posetAlignment->sequenceNumber; sequence4++) {
        if(sequence4 != sequence1 && sequence4 != sequence2 && sequence4 != sequence3) {
            int32_t position4 = getConstraint_lessThanOrEquals(posetAlignment, sequence2, position2, sequence4);
            if(position4 < getConstraint_lessThanOrEquals(posetAlignment, sequence3, position3, sequence4)) { //We have a new transitive constraint..
                addConstraint_lessThanOrEquals(posetAlignment, sequence3, position3, sequence4, position4);
            }
        }
    }
}

static void stPosetAlignment_addP(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {
    //for all pairs do check..
    if(position2 < getConstraint_lessThanOrEquals(posetAlignment, sequence1, position1, sequence2)) {
        addConstraint_lessThanOrEquals(posetAlignment, sequence1, position1, sequence2, position2);
        for(int32_t sequence3=0; sequence3<posetAlignment->sequenceNumber; sequence3++) {
            if(sequence3 != sequence2) {
                if(sequence3 != sequence1) {
                    int32_t position3 = getConstraint_greaterThanOrEquals(posetAlignment, sequence1, position1, sequence3);
                    if(position2 < getConstraint_lessThanOrEquals(posetAlignment, sequence3, position3, sequence2)) { //new constraint found, so add it to the set..
                        addConstraint_lessThanOrEquals(posetAlignment, sequence3, position3, sequence2, position2);
                        stPosetAlignment_addP2(posetAlignment, sequence1, sequence3, position3, sequence2, position2);
                    }
                }
                else {
                    stPosetAlignment_addP2(posetAlignment, INT32_MAX, sequence1, position1, sequence2, position2);
                }
            }
        }
    }
}

bool stPosetAlignment_add(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {
    if(stPosetAlignment_isPossible(posetAlignment, sequence1, position1, sequence2, position2)) {
        stPosetAlignment_addP(posetAlignment, sequence1, position1, sequence2, position2);
        stPosetAlignment_addP(posetAlignment, sequence2, position2, sequence1, position1);
        return 1;
    }
    return 0;
}
