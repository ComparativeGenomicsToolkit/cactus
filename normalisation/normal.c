#include "normal.h"

void makeTerminalNormal(Flower *flower) {
    if (!flower_isTerminal(flower)) {
        Flower_GroupIterator *groupIterator;
        Group *group;
        groupIterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (group_isLeaf(group)) {
                //assert(group_getTotalBaseLength(group) == 0);
                group_makeNestedFlower(group);
                flower_setBuiltBlocks(group_getNestedFlower(group), 1);
            }
        }
        flower_destructGroupIterator(groupIterator);
    }
}
