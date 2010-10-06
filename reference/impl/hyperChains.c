
#include "sonLib.h"
#include "cactus.h"
#include "pressedFlowers.h"

static End *getTerminalEnd(End *end) {
    /*
     * Gets the terminal version of the end.
     */
    return flower_isTerminal(end_getFlower(end)) ? end : getTerminalEnd(flower_getEnd(group_getNestedFlower(end_getGroup(end)), end_getName(end)));
}

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
            if(!end_getSide(end)) {
#ifdef BEN_DEBUG
                assert(end2 == link_get5End(chain_getLink(chain, chain_getLength(chain)-1)));
#endif
                end3 = link_get3End(chain_getLink(chain, 0));
            }
            else {
#ifdef BEN_DEBUG
                assert(end2 == link_get3End(chain_getLink(chain, 0)));
#endif
                end3 = link_get5End(chain_getLink(chain, chain_getLength(chain)-1));
            }
#ifdef BEN_DEBUG
            assert(end3 != end);
            assert(end_getOrientation(end3));
#endif
            if(end_isBlockEnd(end3)) {
                End *end4 = end_getOtherBlockEnd(end3);
#ifdef BEN_DEBUG
                assert(end_getOrientation(end4));
                assert(end4 != end);
#endif
                stHash_insert(hyperChains, getTerminalEnd(end), getTerminalEnd(end4));
            }
#ifdef BEN_DEBUG
            else { //We must be in the root flower
                assert(flower_getParentGroup(flower) == NULL);
                assert(flower == parentFlower);
            }
#endif
        }
        else { //Is a trivial chain
            stHash_insert(hyperChains, getTerminalEnd(end), getTerminalEnd(end2));
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
    assert(flower != NULL);
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
    while(stList_length(pressedFlowers) > 0) {
        constructHyperChainsP(stList_pop(pressedFlowers), hyperChains, flower);
    }
    stList_destruct(pressedFlowers);

#ifdef BEN_DEBUG
    assert(stHash_size(hyperChains) % 2 == 0);
    stHashIterator *endIt = stHash_getIterator(hyperChains);
    End *end;
    while((end = stHash_getNext(endIt)) != NULL) {
        End *end2 = stHash_search(hyperChains, end);
        assert(end2 != NULL);
        assert(stHash_search(hyperChains, end2) == end); //The link should be reciprocal
        assert(flower_isTerminal(end_getFlower(end))); //the ends must be the terminal versions
        assert(flower_isTerminal(end_getFlower(end2))); //the ends must be the terminal versions
    }
    stHash_destructIterator(endIt);
#endif

    return hyperChains;
}
