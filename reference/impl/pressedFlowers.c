
#include "cactus.h"
#include "sonLib.h"

stList *getListOfPressedFlowers(Flower *flower) {
    stList *pressedFlowers = stList_construct();
    stList *stack = stList_construct();
    stList_append(stack, flower);
    while (stList_length(stack) > 0) {
        flower = stList_pop(stack);
        assert(flower != NULL);
        if (flower_isTerminal(flower)) {
            stList_append(pressedFlowers, flower);
        } else {
            Group *group;
            Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
            while ((group = flower_getNextGroup(groupIt)) != NULL) {
#ifdef BEN_DEBUG
                assert(!group_isLeaf(group));
#endif
                if (group_isTangle(group)) {
                    stList_append(stack, group_getNestedFlower(group));
                }
            }
            flower_destructGroupIterator(groupIt);
        }
    }
    stList_destruct(stack);
    return pressedFlowers;
}
