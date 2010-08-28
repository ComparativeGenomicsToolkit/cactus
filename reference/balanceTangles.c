
#include "cactus.h"
#include "sonLib.h"

/*
 * balanceTangles.c
 *
 *  Created on: 10 Aug 2010
 *      Author: benedictpaten
 */

/*
 * Code to ensure all flowers contain an even number of stub ends, the condition
 * necessary to construct a reference genome. It does this by inventing blocks, of zero
 * length and with zero elements to balance the groups.
 */

static void balanceTanglesP(Group *oddTangle, stList *ends) {
    int32_t i = group_getAttachedStubEndNumber(oddTangle)
            + group_getBlockEndNumber(oddTangle);
#ifdef BEN_DEBUG
    assert(i > 0);
    assert(i % 2 != 0);
    assert(stList_length(ends) == 1 || stList_length(ends) == 3);
#endif
    Flower *flower = group_getFlower(oddTangle);
    if (i == 1 && stList_length(ends) == 1) { //We add an extra block in
        Block *block = block_construct(1, flower);
        stList_append(ends, block_get5End(block));
        stList_append(ends, block_get3End(block));
    }
    for (int32_t i = 0; i < stList_length(ends); i++) {
        End *end = stList_get(ends, i);
        assert(end_getFlower(end) == group_getFlower(oddTangle));
        end_setGroup(end, oddTangle);
    }
#ifdef BEN_DEBUG
    i = group_getAttachedStubEndNumber(oddTangle)
                + group_getBlockEndNumber(oddTangle);
    assert(i > 0);
    assert(i % 2 == 0);
    assert(i != 2);
#endif
    if (!group_isLeaf(oddTangle)) {
        Flower_GroupIterator *groupIt = flower_getGroupIterator(
                group_getNestedFlower(oddTangle));
        Group *group;
        while ((group = flower_getNextGroup(groupIt)) != NULL) {
            int32_t i = group_getAttachedStubEndNumber(group)
                    + group_getBlockEndNumber(group);
            assert(i > 0);
            if (i % 2 != 0) {
                //Update the ends..
                for(int32_t j=0; j<stList_length(ends); j++) {
                    stList_set(ends, j, end_copyConstruct(stList_get(ends, j), group_getFlower(group)));
                }
                balanceTanglesP(group, ends);
                flower_destructGroupIterator(groupIt);
                return;
            }
        }
        assert(0);
        st_errAbort("We are exiting because we could not find an odd tangle\n");
    }
    //ensure we have normalised it (but after the recursion on odd tangles)
    flower_makeTerminalNormal(flower);
}

void balanceTangles(Flower *flower) {
    assert(flower_getAttachedStubEndNumber(flower) > 0);
    assert(flower_getAttachedStubEndNumber(flower) % 2 == 0);
    assert(flower_getBlockEndNumber(flower) % 2 == 0);

    stList *oddTangles = stList_construct();
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        int32_t i = group_getAttachedStubEndNumber(group)
                + group_getBlockEndNumber(group);
        assert(i > 0);
        if (i % 2 != 0) {
            stList_append(oddTangles, group);
        }
    }
    flower_destructGroupIterator(groupIt);
    assert(stList_length(oddTangles) % 2 == 0);

    //Create the link blocks
    stList *linkEnds = stList_construct();
    for (int32_t i = stList_length(oddTangles) / 2; i > 0; i--) {
        Block *block = block_construct(1, flower);
        stList_append(linkEnds, block_get5End(block));
        stList_append(linkEnds, block_get3End(block));
    }

    //Assign the link ends to the groups.
    while (stList_length(oddTangles) > 0) {
        Group *oddTangle = stList_pop(oddTangles);
        stList *ends = stList_construct();
        stList_append(ends, stList_pop(linkEnds));
        balanceTanglesP(oddTangle, ends);
        stList_destruct(ends);
    }
    assert(stList_length(linkEnds) == 0);
    stList_destruct(linkEnds);
    stList_destruct(oddTangles);
}
