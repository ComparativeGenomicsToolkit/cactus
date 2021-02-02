#include "traverseFlowers.h"

void extendFlowers(Flower *flower, stList *extendedFlowers, int64_t minFlowerSize) {
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    assert(flower_builtBlocks(
            flower)); //This recursion depends on the block structure having been properly defined for all nodes.
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        if (group_isLeaf(group)) { // Has no nested flower
            int64_t size = group_getTotalBaseLength(group);
            assert(size >= 0);
            if (size >= minFlowerSize) {
                Flower *nestedFlower = group_makeNestedFlower(group);
                stList_append(extendedFlowers, nestedFlower);
            }
        }
        else { // Recursively search for more nested flowers to exnted
            Flower *nestedFlower = group_getNestedFlower(group);
            assert(nestedFlower != NULL);
            extendFlowers(nestedFlower, extendedFlowers, minFlowerSize);
        }
    }
    flower_destructGroupIterator(groupIterator);
}

/*
 * Get all the child flowers of a given flower.
 */
void getChildFlowers(Flower *flower, stList *children) {
    if (!flower_isLeaf(flower)) {
        assert(flower_builtBlocks(
                flower)); //This recursion depends on the block structure having been properly defined for all nodes.
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if (!group_isLeaf(group)) {
                Flower *nestedFlower = group_getNestedFlower(group);
                assert(nestedFlower != NULL);
                stList_append(children, nestedFlower);
            }
        }
        flower_destructGroupIterator(groupIterator);
    }
}

stList *getFlowerHierarchyInLayers(Flower *rootFlower) {
    stList *flowers = stList_construct();
    stList_append(flowers, rootFlower);
    stList *flowerLayers = stList_construct3(0, (void (*)(void *)) stList_destruct);
    while (stList_length(flowers) > 0) {
        stList *childFlowers = stList_construct();
        for(int64_t i=0; i<stList_length(flowers); i++) {
            getChildFlowers(stList_get(flowers, i), childFlowers);
        }
        stList_append(flowerLayers, flowers);
        flowers = childFlowers;
    }
    stList_destruct(flowers);
    return flowerLayers;
}
