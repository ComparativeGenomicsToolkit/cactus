#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
static Flower *flower;
static Flower *nestedFlower1;
static Flower *nestedFlower2;
static End *end1;
static Block *block;
static End *end2;
static Group *group1;
static Group *group2;
static Chain *chain;
static Link *link1;
static Link *link2;

static void cactusChainsSharedTestTeardown() {
    if (cactusDisk != NULL) {
        cactusDisk_destruct(cactusDisk);
        testCommon_deleteTemporaryCactusDisk();
        cactusDisk = NULL;
    }
}

static void cactusChainsSharedTestSetup() {
    cactusChainsSharedTestTeardown();
    cactusDisk = cactusDisk_construct(testCommon_getTemporaryCactusDisk());
    flower = flower_construct(cactusDisk);
    nestedFlower1 = flower_construct(cactusDisk);
    nestedFlower2 = flower_construct(cactusDisk);
    end1 = end_construct2(0, 0, flower);
    end2 = end_construct(0, flower);
    block = block_construct(2, flower);
    end_copyConstruct(end1, nestedFlower1);
    end_copyConstruct(block_get5End(block), nestedFlower1);
    end_copyConstruct(block_get3End(block), nestedFlower2);
    end_copyConstruct(end2, nestedFlower2);
    group1 = group_construct(flower, nestedFlower1);
    group2 = group_construct(flower, nestedFlower2);
    chain = chain_construct(flower);
    link1 = link_construct(end1, block_get5End(block), group1, chain);
    link2 = link_construct(block_get3End(block), end2, group2, chain);
}
