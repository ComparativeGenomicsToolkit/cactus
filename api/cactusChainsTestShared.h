#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
static Flower *flower;
static Flower *nestedNet1;
static Flower *nestedNet2;
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
        testCommon_deleteTemporaryNetDisk();
        cactusDisk = NULL;
    }
}

static void cactusChainsSharedTestSetup() {
    cactusChainsSharedTestTeardown();
    cactusDisk = cactusDisk_construct(testCommon_getTemporaryNetDisk());
    flower = flower_construct(cactusDisk);
    nestedNet1 = flower_construct(cactusDisk);
    nestedNet2 = flower_construct(cactusDisk);
    end1 = end_construct2(0, 0, flower);
    end2 = end_construct(0, flower);
    block = block_construct(2, flower);
    end_copyConstruct(end1, nestedNet1);
    end_copyConstruct(block_get5End(block), nestedNet1);
    end_copyConstruct(block_get3End(block), nestedNet2);
    end_copyConstruct(end2, nestedNet2);
    group1 = group_construct(flower, nestedNet1);
    group2 = group_construct(flower, nestedNet2);
    chain = chain_construct(flower);
    link1 = link_construct(end1, block_get5End(block), group1, chain);
    link2 = link_construct(block_get3End(block), end2, group2, chain);
}
