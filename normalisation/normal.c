#include "normal.h"

void makeTerminalNormal(Flower *net) {
    if (!flower_isTerminal(net)) {
        Flower_GroupIterator *groupIterator;
        Group *group;
        groupIterator = flower_getGroupIterator(net);
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (group_isLeaf(group)) {
                //assert(group_getTotalBaseLength(group) == 0);
                group_makeNestedNet(group);
                flower_setBuiltBlocks(group_getNestedNet(group), 1);
            }
        }
        flower_destructGroupIterator(groupIterator);
    }
}
