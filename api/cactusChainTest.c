#include "cactusChainsTestShared.h"

static bool nestedTest = 0;

static void cactusChainTestTeardown() {
    if (!nestedTest) {
        cactusChainsSharedTestTeardown();
    }
}

static void cactusChainTestSetup() {
    if (!nestedTest) {
        cactusChainsSharedTestSetup();
    }
}

void testChain_construct(CuTest* testCase) {
    nestedTest = 0;
    cactusChainTestSetup();
    CuAssertTrue(testCase, chain != NULL);
    CuAssertTrue(testCase, link1 != NULL);
    CuAssertTrue(testCase, link2 != NULL);
    cactusChainTestTeardown();
}

void testChain_getLink(CuTest* testCase) {
    cactusChainTestSetup();
    CuAssertTrue(testCase, chain_getLink(chain, 0) == link1);
    CuAssertTrue(testCase, chain_getLink(chain, 1) == link2);
    //CuAssertTrue(testCase, 0);
    cactusChainTestTeardown();
}

void testChain_getLength(CuTest* testCase) {
    cactusChainTestSetup();
    CuAssertTrue(testCase, chain_getLength(chain) == 2);
    cactusChainTestTeardown();
}

void testChain_getBlockChain(CuTest* testCase) {
    cactusChainTestSetup();
    int32_t i;
    Block **blockChain = chain_getBlockChain(chain, &i);
    CuAssertTrue(testCase, i == 1);
    CuAssertTrue(testCase, blockChain[0] == block);
    free(blockChain);
    cactusChainTestTeardown();
}

void testChain_getName(CuTest* testCase) {
    cactusChainTestSetup();
    CuAssertTrue(testCase, chain_getName(chain) != NULL_NAME);
    CuAssertTrue(testCase, flower_getChain(flower, chain_getName(chain)) == chain);
    cactusChainTestTeardown();
}

void testChain_getNet(CuTest* testCase) {
    cactusChainTestSetup();
    CuAssertTrue(testCase, chain_getNet(chain) == flower);
    cactusChainTestTeardown();
}

void testChain_rename(CuTest* testCase) {
    cactusChainTestSetup();
    Name newName = cactusDisk_getUniqueID(cactusDisk);
    Name oldName = chain_getName(chain);
    chain_rename(chain, newName);
    assert(chain_getName(chain) == newName);
    assert(flower_getChain(flower, newName) == chain);
    assert(flower_getChain(flower, oldName) == NULL);
    cactusChainTestTeardown();
}

void testChain_serialisation(CuTest* testCase) {
    cactusChainTestSetup();
    int32_t i;
    void
            *vA =
                    binaryRepresentation_makeBinaryRepresentation(chain,
                            (void(*)(void *, void(*)(const void *, size_t,
                                    size_t))) chain_writeBinaryRepresentation,
                            &i);
    CuAssertTrue(testCase, i> 0);
    chain_destruct(chain);
    void *vA2 = vA;
    chain = chain_loadFromBinaryRepresentation(&vA2, flower);
    CuAssertTrue(testCase, chain_getLength(chain) == 2);
    link1 = chain_getLink(chain, 0);
    link2 = chain_getLink(chain, 1);
    nestedTest = 1;
    testChain_getLink(testCase);
    testChain_getLength(testCase);
    testChain_getName(testCase);
    testChain_getNet(testCase);
    nestedTest = 0;
    cactusChainTestTeardown();
}

CuSuite* cactusChainTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, testChain_getLink);
    SUITE_ADD_TEST(suite, testChain_getLength);
    SUITE_ADD_TEST(suite, testChain_getBlockChain);
    SUITE_ADD_TEST(suite, testChain_getName);
    SUITE_ADD_TEST(suite, testChain_getNet);
    SUITE_ADD_TEST(suite, testChain_serialisation);
    SUITE_ADD_TEST(suite, testChain_rename);
    SUITE_ADD_TEST(suite, testChain_construct);
    return suite;
}
