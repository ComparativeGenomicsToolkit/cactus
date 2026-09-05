#include <stdlib.h>
#include "cactus.h"
#include "sonLib.h"

/*
 * A cap with its sort key pulled out once: the comparison below is run O(n log n) times per flower and
 * the cap accessors it replaced are not cheap, so the key is computed once per cap instead.
 */
typedef struct _capWithKey {
    Name sequenceName;
    int64_t coordinate;
    bool side;
    Cap *cap;
} CapWithKey;

static int addAdjacenciesPP(const void *a, const void *b) {
    const CapWithKey *cap1 = a, *cap2 = b;
    int64_t i = cactusMisc_nameCompare(cap1->sequenceName, cap2->sequenceName);
    if (i == 0) {
        int64_t j = cap1->coordinate;
        int64_t k = cap2->coordinate;
        i = j > k ? 1 : (j < k ? -1 : 0);
        if (i == 0) {
            assert(cap_getSegment(cap1->cap) == cap_getSegment(cap2->cap));
            j = cap1->side;
            k = cap2->side;
            assert((j && !k) || (!j && k));
            i = j ? -1 : 1;
        }
    }
    return i;
}

void stCaf_addAdjacencies(Flower *flower) {
    //Build an array of caps with their keys.
    int64_t capNumber = 0, capCapacity = 64;
    CapWithKey *caps = st_malloc(capCapacity * sizeof(CapWithKey));
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
        Cap *cap;
        while ((cap = end_getNext(instanceIterator)) != NULL) {
            if (!cap_getStrand(cap)) {
                cap = cap_getReverse(cap);
            }
            if (capNumber == capCapacity) {
                capCapacity *= 2;
                caps = st_realloc(caps, capCapacity * sizeof(CapWithKey));
            }
            caps[capNumber].sequenceName = sequence_getName(cap_getSequence(cap));
            caps[capNumber].coordinate = cap_getCoordinate(cap);
            caps[capNumber].side = cap_getSide(cap);
            caps[capNumber].cap = cap;
            capNumber++;
        }
        end_destructInstanceIterator(instanceIterator);
    }
    flower_destructEndIterator(endIterator);
    assert(capNumber % 2 == 0);
    //Sort the caps (the key is a total order over the caps of a flower, so the sort algorithm does not matter).
    qsort(caps, capNumber, sizeof(CapWithKey), addAdjacenciesPP);
    //Now make the adjacencies.
    for (int64_t i = 1; i < capNumber; i += 2) {
        cap_makeAdjacent(caps[i - 1].cap, caps[i].cap);
    }
    //Clean up.
    free(caps);
}
