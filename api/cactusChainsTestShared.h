#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk;
static Net *net;
static Net *nestedNet1;
static Net *nestedNet2;
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
    net = net_construct(cactusDisk);
    nestedNet1 = net_construct(cactusDisk);
    nestedNet2 = net_construct(cactusDisk);
    end1 = end_construct2(0, 0, net);
    end2 = end_construct(0, net);
    block = block_construct(2, net);
    end_copyConstruct(end1, nestedNet1);
    end_copyConstruct(block_get5End(block), nestedNet1);
    end_copyConstruct(block_get3End(block), nestedNet2);
    end_copyConstruct(end2, nestedNet2);
    group1 = group_construct(net, nestedNet1);
    group2 = group_construct(net, nestedNet2);
    chain = chain_construct(net);
    link1 = link_construct(end1, block_get5End(block), group1, chain);
    link2 = link_construct(block_get3End(block), end2, group2, chain);
}
