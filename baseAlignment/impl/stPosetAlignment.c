
#include "sonLib.h"
#include "stPosetAlignment.h"
#include "stdlib.h"

struct _stPosetAlignment {
    int32_t sequenceNumber;
    stSortedSet **constraintLists;
};

static int comparePositions(stIntTuple *position1, stIntTuple *position2) {
    if(stIntTuple_getPosition(position1, 0) == INT32_MAX || stIntTuple_getPosition(position2, 0) == INT32_MAX) { //Indicates we should ignore the first position and compare the second.
#ifdef BEN_DEBUG
        assert(stIntTuple_getPosition(position1, 1) != INT32_MAX);
        assert(stIntTuple_getPosition(position2, 1) != INT32_MAX);
#endif
        return stIntTuple_getPosition(position1, 1) - stIntTuple_getPosition(position2, 1);
    }
    return stIntTuple_getPosition(position1, 0) - stIntTuple_getPosition(position2, 0);
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
#ifdef BEN_DEBUG
    assert(sequence1 >= 0 && sequence1 < posetAlignment->sequenceNumber);
    assert(sequence2 >= 0 && sequence2 < posetAlignment->sequenceNumber);
    assert(sequence1 != sequence2);
#endif
    return posetAlignment->constraintLists[sequence1 * posetAlignment->sequenceNumber + sequence2];
}

/*
 * Gets the position in sequence2 that the position in sequence2 must be less or equal to in the alignment
 */
static stIntTuple *getConstraint_lessThan(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2) {
    stIntTuple *pos = stIntTuple_construct(2, position1, INT32_MAX);
    //Get less than or equal
    stIntTuple *constraint = stSortedSet_searchGreaterThanOrEqual(getConstraintList(posetAlignment, sequence1, sequence2), pos);
    stIntTuple_destruct(pos);
    return constraint;
}

/*
 * Gets the position in sequence2 that the position in sequence1 must be greater than or equal to in the alignment
 */
static stIntTuple *getConstraint_greaterThan(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2) {
    stIntTuple *pos = stIntTuple_construct(2, INT32_MAX, position1);
    //Get less than or equal
    stIntTuple *constraint = stSortedSet_searchLessThanOrEqual(getConstraintList(posetAlignment, sequence2, sequence1), pos);
    stIntTuple_destruct(pos);
    return constraint;
}

/*
 * Returns non-zero iff the constraint is prime. The lessThanOrEquals argument, if non-zero specifies the constraint is less than equals.
 */
static bool lessThanConstraintIsPrime(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2, int32_t lessThanOrEquals) {
    stIntTuple *constraint = getConstraint_lessThan(posetAlignment, sequence1, position1, sequence2);
    if(constraint == NULL) {
        return 1;
    }
    if(stIntTuple_getPosition(constraint, 2)) { //less than or equals
        return position2 < stIntTuple_getPosition(constraint, 1) + (lessThanOrEquals ? 0 : 1);
    }
    return position2 < stIntTuple_getPosition(constraint, 1); //just less than
}

/*
 * Adds a prime less than (or equals) constraint to the list of prime constraints, removing any redundant constraints in the process.
 * The or equals is specified by making the lessThanOrEquals argument non-zero.
 */
void addConstraint_lessThan(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2, int32_t lessThanOrEquals) {
    stSortedSet *constraintList = getConstraintList(posetAlignment, sequence1, sequence2);
#ifdef BEN_DEBUG
    assert(position1 != INT32_MAX);
    assert(position2 != INT32_MAX);
#endif
    stIntTuple *constraint1 = stIntTuple_construct(3, position1, position2, lessThanOrEquals);
    stIntTuple *constraint2;
    while((constraint2 = stSortedSet_searchLessThanOrEqual(constraintList, constraint1)) != NULL) {
        assert(stIntTuple_getPosition(constraint2, 0) <= position1);
        if(stIntTuple_getPosition(constraint2, 1) >= position2) {
#ifdef BEN_DEBUG
            if(stIntTuple_getPosition(constraint2, 1) == position2) { //Check we are not removing an equivalent or more severe constraint.
                assert((!lessThanOrEquals && stIntTuple_getPosition(constraint2, 2)) || stIntTuple_getPosition(constraint2, 0) < position1);
            }
#endif
            stSortedSet_remove(constraintList, constraint2);
            stIntTuple_destruct(constraint2);
        }
        else {
#ifdef BEN_DEBUG
            assert(stIntTuple_getPosition(constraint2, 0) < position1); //Check the constraint does not overshadow our proposed constraint.
#endif
            break;
        }
    }
    stSortedSet_insert(constraintList, constraint1);
}

bool stPosetAlignment_isPossibleP(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {
    stIntTuple *constraint = getConstraint_lessThan(posetAlignment, sequence1, position1, sequence2);
    if(constraint == NULL) {
        return 1;
    }
    if(stIntTuple_getPosition(constraint, 2) && stIntTuple_getPosition(constraint, 0) == position1) { //less than or equals
        return position2 <= stIntTuple_getPosition(constraint, 1);
    }
    else {
        return position2 < stIntTuple_getPosition(constraint, 1);
    }
}

bool stPosetAlignment_isPossible(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {
    return stPosetAlignment_isPossibleP(posetAlignment, sequence1, position1, sequence2, position2) &&
    stPosetAlignment_isPossibleP(posetAlignment, sequence2, position2, sequence1, position1);
}

static void stPosetAlignment_addP2(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t sequence3, int32_t position3, int32_t sequence2, int32_t position2, int32_t lessThanOrEqual) {
    for(int32_t sequence4=0; sequence4<posetAlignment->sequenceNumber; sequence4++) {
        if(sequence4 != sequence1 && sequence4 != sequence2 && sequence4 != sequence3) {
            stIntTuple *constraint = getConstraint_lessThan(posetAlignment, sequence2, position2, sequence4);
            if(constraint != NULL) {
                int32_t position4 = stIntTuple_getPosition(constraint, 1);
                int32_t transLessThanOrEqual = lessThanOrEqual && stIntTuple_getPosition(constraint, 2) && stIntTuple_getPosition(constraint, 0) == position2; //stuff which maintains the less than or equals
                if(lessThanConstraintIsPrime(posetAlignment, sequence3, position3, sequence4, position4, transLessThanOrEqual)) {//We have a new transitive constraint..
                    addConstraint_lessThan(posetAlignment, sequence3, position3, sequence4, position4, transLessThanOrEqual);
                }
            }
        }
    }
}

static void stPosetAlignment_addP(stPosetAlignment *posetAlignment, int32_t sequence1, int32_t position1, int32_t sequence2, int32_t position2) {
    //for all pairs do check..
    if(lessThanConstraintIsPrime(posetAlignment, sequence1, position1, sequence2, position2, 1)) {
        addConstraint_lessThan(posetAlignment, sequence1, position1, sequence2, position2, 1);
        for(int32_t sequence3=0; sequence3<posetAlignment->sequenceNumber; sequence3++) {
            if(sequence3 != sequence2) {
                if(sequence3 != sequence1) {
                    stIntTuple *constraint = getConstraint_greaterThan(posetAlignment, sequence1, position1, sequence3);
                    if(constraint != NULL) {
                        int32_t position3 = stIntTuple_getPosition(constraint, 0); //its reversed
                        int32_t lessThanOrEqual = stIntTuple_getPosition(constraint, 2) && stIntTuple_getPosition(constraint, 1) == position1;
                        if(lessThanConstraintIsPrime(posetAlignment, sequence3, position3, sequence2, position2, lessThanOrEqual)) { //new constraint found, so add it to the set..
                            addConstraint_lessThan(posetAlignment, sequence3, position3, sequence2, position2, lessThanOrEqual);
                            stPosetAlignment_addP2(posetAlignment, sequence1, sequence3, position3, sequence2, position2, lessThanOrEqual);
                        }
                    }
                }
                else {
                    stPosetAlignment_addP2(posetAlignment, INT32_MAX, sequence1, position1, sequence2, position2, 1);
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
