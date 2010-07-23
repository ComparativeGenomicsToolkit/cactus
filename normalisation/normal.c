#include "normal.h"

void makeTerminalNormal(Net *net) {
    if (!net_isTerminal(net)) {
        Net_GroupIterator *groupIterator;
        Group *group;
        groupIterator = net_getGroupIterator(net);
        while ((group = net_getNextGroup(groupIterator)) != NULL) {
            if (group_isLeaf(group)) {
                //assert(group_getTotalBaseLength(group) == 0);
                group_makeNestedNet(group);
                net_setBuiltBlocks(group_getNestedNet(group), 1);
            }
        }
        net_destructGroupIterator(groupIterator);
    }
}
