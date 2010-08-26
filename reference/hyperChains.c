
#include "sonLib.h"
#include "cactus.h"
#include "pressedFlowers.h"

static void constructHyperChainsP2(Group *group, Name endName, Flower *parentFlower, stHash *hyperChains) {
    if(group == NULL) { //The problem was at the root, so no parent.
        return;
    }
    Flower *flower = group_getFlower(group);
    End *end = flower_getEnd(flower, endName);
    assert(end != NULL);
    if(end_isBlockEnd(end)) { //We have found the hyper link, so we use it.
        End *end2 = end_getOtherBlockEnd(end);
        assert(end_getOrientation(end2));
        Link *link = group_getLink(end_getGroup(end2));
        if(link != NULL) {
            Chain *chain = link_getChain(link);
            End *end3;
            if(end_getSide(end)) {
                assert(end2 == link_get5End(chain_getLink(chain, chain_getLength(chain)-1)));
                end3 = link_get3End(chain_getLink(chain, 0));
            }
            else {
                assert(end2 == link_get3End(chain_getLink(chain, 0)));
                end3 = link_get5End(chain_getLink(chain, chain_getLength(chain)-1));
            }
            assert(end3 != end);
            assert(end_getOrientation(end3));
            assert(end_isBlockEnd(end3));
            End *end4 = end_getOtherBlockEnd(end3);
            assert(end_getOrientation(end4));
            assert(end4 != end);
            stHash_insert(hyperChains, end, end4);
        }
        else { //Is a trivial chain
            stHash_insert(hyperChains, end, end2);
        }
        return;
    }
    if(flower == parentFlower) { //We can not go up further than the parent flower.
        return;
    }
    constructHyperChainsP2(flower_getParentGroup(flower), endName, parentFlower, hyperChains);
}

static void constructHyperChainsP(Flower *flower, stHash *hyperChains, Flower *parentFlower) {
    End *end;
    assert(flower_isTerminal(flower));
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    while((end = flower_getNextEnd(endIt)) != NULL) {
        assert(!end_isBlockEnd(end));
        if(end_isAttached(end)) {
            if(flower != parentFlower) { //We can not go above the parent to define links..
                constructHyperChainsP2(flower_getParentGroup(flower), end_getName(end), parentFlower, hyperChains);
            }
        }
    }
    flower_destructEndIterator(endIt);
}

stHash *constructHyperChains(Flower *flower) {
    stList *pressedFlowers = getListOfPressedFlowers(flower);
    stHash *hyperChains = stHash_construct();
    while(stList_length(pressedFlowers)) {
        constructHyperChainsP(stList_pop(pressedFlowers), hyperChains, flower);
    }
    stList_destruct(pressedFlowers);
    return hyperChains;
}
