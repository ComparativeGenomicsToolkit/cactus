#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions used to normalise the structure of the cactus graph.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static void promoteBlock(Block *block, Flower *flower, Flower *parentFlower) {
    /*
     * Redirects all the pointers in the block to the higher level.
     */
    Segment *segment;
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    while ((segment = block_getNext(it)) != NULL) {
        flower_removeSegment(flower, segment);
        flower_addSegment(parentFlower, segment);
    }
    block_destructInstanceIterator(it);
    flower_removeBlock(flower, block);
    flower_addBlock(parentFlower, block);
    block->blockContents->flower = parentFlower;
}

static void mergeStubEnd(End *end, Flower *flower, Flower *parentFlower) {
    /*
     * Removes the lower level stub end, reconnects the adjacencies of the higher level end to that of the lower level stub end.
     */
    assert(end_isStubEnd(end));
    assert(group_getLink(end_getGroup(end)) != NULL);
    End *parentEnd = flower_getEnd(parentFlower, end_getName(end));
    assert(parentEnd != NULL); //end is already present, so we will get rid of the lower level end, reconnecting the adjacencies..
    assert(end_getInstanceNumber(parentEnd) == end_getInstanceNumber(end));

    Group *group = end_getGroup(end);
    End_InstanceIterator *it = end_getInstanceIterator(end); //We ensure that the adjacencies match the higher level flower..
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        Cap *parentCap = end_getInstance(parentEnd, cap_getName(cap));
        assert(parentCap != NULL);
        Cap *adjacentCap = cap_getAdjacency(cap);
        if (cap_getAdjacency(cap) != NULL) {
            Cap *parentAdjacentCap = flower_getCap(parentFlower, cap_getName(
                    adjacentCap));
            assert(parentAdjacentCap != NULL);
            assert(cactusDisk_getFlower(flower_getCactusDisk(flower), flower_getName(end_getFlower(cap_getEnd(adjacentCap)))) == end_getFlower(cap_getEnd(adjacentCap)));
            cap_breakAdjacency(cap); //We must remove the old adjacency so it doesn't hang around.
            cap_makeAdjacent(parentCap, parentAdjacentCap);
        }
    }
    end_destruct(end); //Destruct the old end
    assert(group_getFlower(group) == parentFlower); //the group should already have been promoted
    end_setGroup(parentEnd, group); //ensures the parent end is in the group of the lower level flower..
}

static void promoteBlockEnd(End *end, Flower *flower, Flower *parentFlower) {
    /*
     * Redirects all the pointers in the block end to the higher level flower.
     */
    assert(end_isBlockEnd(end));
    assert(flower_getEnd(parentFlower, end_getName(end)) == NULL);
    Cap *cap;
    End_InstanceIterator *it = end_getInstanceIterator(end);
    EventTree *eventTree = flower_getEventTree(parentFlower);
    while ((cap = end_getNext(it)) != NULL) {
        Event *event = eventTree_getEvent(eventTree, event_getName(
                cap_getEvent(cap)));
        assert(event != NULL);
        cap->capContents->event = event;
        if (cap_getSequence(cap) != NULL) {
            Sequence *sequence = flower_getSequence(parentFlower, sequence_getName(
                    cap_getSequence(cap)));
            assert(sequence != NULL);
            cap->capContents->sequence = sequence;
        }
        flower_removeCap(flower, cap);
        flower_addCap(parentFlower, cap);
    }
    end_destructInstanceIterator(it);
    flower_removeEnd(flower, end);
    flower_addEnd(parentFlower, end);
    end->endContents->flower = parentFlower;
}

static void promoteEndsBlocksAndGroups(Chain *chain, Flower *flower, Flower *parentFlower) {
    /*
     * Promotes all the ends and blocks in the chain..
     */
    for (int32_t i = 0; i < chain_getLength(chain); i++) {
        Link *link = chain_getLink(chain, i);
        End *_3End = link_get3End(link);
        End *_5End = link_get5End(link);
        if (end_isBlockEnd(_3End)) {
            promoteBlock(end_getBlock(_3End), flower, parentFlower);
            promoteBlockEnd(_3End, flower, parentFlower);
        } else {
            assert(i == 0);
        }
        if (end_isBlockEnd(_5End)) {
            if (i + 1 == chain_getLength(chain)) {
                promoteBlock(end_getBlock(_5End), flower, parentFlower);
            }
            promoteBlockEnd(_5End, flower, parentFlower);
        } else {
            assert(i == chain_getLength(chain)-1);
        }
        Group *group = link_getGroup(link);
        //Get the list of free stub ends to promote (we make a list as
        //when we promote them we disrupt the ordering of the list)
        stList *freeStubEndsToPromote = stList_construct();
        Group_EndIterator *groupEndIt = group_getEndIterator(group);
        End *end;
        while ((end = group_getNextEnd(groupEndIt)) != NULL) {
            if (end_isStubEnd(end) && end_isFree(end)) {
                stList_append(freeStubEndsToPromote, end);
            }
        }
        group_destructEndIterator(groupEndIt);
        //Promote group
        flower_removeGroup(flower, group);
        flower_addGroup(parentFlower, group);
        group->flower = parentFlower;
        Flower *nestedFlower = group_getNestedFlower(group);
        if (nestedFlower != NULL) {
            nestedFlower->parentFlowerName = flower_getName(parentFlower);
        }
        //Promote any free stub ends..
        while (stList_length(freeStubEndsToPromote) > 0) {
            mergeStubEnd(stList_pop(freeStubEndsToPromote), flower, parentFlower);
        }
        stList_destruct(freeStubEndsToPromote);
    }
    //Now do the stub ends of the chain..
    End *_3End = link_get3End(chain_getLink(chain, 0));
    if (end_isStubEnd(_3End)) {
        assert(end_isAttached(_3End));
        assert(end_getFlower(_3End) == flower);
        mergeStubEnd(_3End, flower, parentFlower);
    }
    End *_5End = link_get5End(chain_getLink(chain, chain_getLength(chain) - 1));
    if (end_isStubEnd(_5End)) {
        assert(end_isAttached(_5End));
        assert(end_getFlower(_5End) == flower);
        mergeStubEnd(_5End, flower, parentFlower);
    }
}

static Cap *promoteChainEndP(Cap *cap, Flower *parentFlower) {
    End *end = cap_getEnd(cap);
    if (end_isBlockEnd(end)) {
        assert(flower_getEnd(parentFlower, end_getName(end)) == NULL);
        Cap *otherCap = cap_getOtherSegmentCap(cap);
        assert(otherCap != NULL);
        assert(cap != otherCap);
        assert(cap_getEnd(otherCap) != cap_getEnd(cap));
        Cap *adjacentCap = cap_getAdjacency(otherCap);
        assert(adjacentCap != NULL);
        return promoteChainEndP(adjacentCap, parentFlower);
    }
    Cap *higherCap = flower_getCap(parentFlower, cap_getName(cap));
    assert(higherCap != NULL);
    assert(higherCap != cap);
    assert(cap_getName(higherCap) == cap_getName(cap));
    return higherCap;
}

static void promoteChainEnd(End *end, Flower *flower, Flower *parentFlower) {
    /*
     * Promotes the block end of a chain to the parent flower, creating a stub end
     * in the lower level flower.
     */
    assert(end_isBlockEnd(end));
    assert(end_getOrientation(end));
    assert(end_getFlower(end) == flower);
    //promote the block end
    promoteBlockEnd(end, flower, parentFlower);
    //Sort out the group...
    Group *childGroup = end_getGroup(end);
    end_setGroup(end, flower_getParentGroup(flower));

    //copy construct them into the children
    End *childEnd = end_copyConstruct(end, flower);
    end_setGroup(childEnd, childGroup);

    //copy the adjacencies in the lower level
    assert(end_getInstanceNumber(end) == end_getInstanceNumber(childEnd));
    Cap *cap;
    End_InstanceIterator *capIt = end_getInstanceIterator(end);
    while ((cap = end_getNext(capIt)) != NULL) {
        Cap *childCap = end_getInstance(childEnd, cap_getName(cap));
        Cap *childAdjacentCap = cap_getAdjacency(cap);
        assert(childCap != NULL);
        assert(childAdjacentCap != NULL);
        assert(end_getFlower(cap_getEnd(childCap)) == flower);
        assert(end_getFlower(cap_getEnd(childAdjacentCap)) == flower);
        assert(end_getGroup(cap_getEnd(childCap)) == end_getGroup(cap_getEnd(childAdjacentCap)));
        cap_makeAdjacent(childCap, childAdjacentCap); //reconnect the adjacency.
    }
    end_destructInstanceIterator(capIt);

    //create the attachments in the parent
    capIt = end_getInstanceIterator(childEnd);
    while ((cap = end_getNext(capIt)) != NULL) {
        Cap *parentCap = end_getInstance(end, cap_getName(cap));
        assert(parentCap != NULL);
        Cap *adjacentCap = cap_getAdjacency(cap);
        assert(adjacentCap != NULL);
        assert(end_getFlower(cap_getEnd(adjacentCap)) == flower);
        //Search for the higher level adjacency..
        Cap *adjacentParentCap = promoteChainEndP(adjacentCap, parentFlower);
        assert(adjacentParentCap != NULL);
        assert(end_getGroup(cap_getEnd(parentCap)) == end_getGroup(cap_getEnd(adjacentParentCap)));
        cap_makeAdjacent(parentCap, adjacentParentCap);
    }
    end_destructInstanceIterator(capIt);
}

/////////
//Functions which calculate the chain structure after the promotion of a chain.
////////

static void addToChainList(Chain *chain, int32_t start, int32_t end,
        stList *chainList) {
    /*
     * Adds the links in the chain between the specified start (inc) and end (exc) indices to the list
     * chainList.
     */
    for (int32_t i = start; i < end; i++) {
        Link *link = chain_getLink(chain, i);
        stList_append(chainList, cactusMisc_nameToString(end_getName(link_get3End(
                link))));
        stList_append(chainList, cactusMisc_nameToString(end_getName(link_get5End(
                link))));
    }
}

static void getMaximalChain_between(Link *parentLink, Chain *chain,
        stList *chainList) {
    Chain *parentChain = link_getChain(parentLink);
    addToChainList(parentChain, 0, link_getIndex(parentLink), chainList);
    //Get the ends..
    End *parent3End = link_get3End(parentLink);
    End *parent5End = link_get5End(parentLink);
    End *_3End = link_get3End(chain_getLink(chain, 0));
    End *_5End = link_get5End(chain_getLink(chain, chain_getLength(chain) - 1));
    Name parent3EndName = end_getName(parent3End);
    Name parent5EndName = end_getName(parent5End);
    Name _3EndName = end_getName(_3End);
    Name _5EndName = end_getName(_5End);
    assert(parent3EndName != _5EndName);
    if (parent3EndName == _3EndName) {
        addToChainList(chain, 0, chain_getLength(chain), chainList);
        if (parent5EndName != _5EndName) { //We need to introduce a link between the other block end and this link..
            assert(end_isBlockEnd(_5End));
            stList_append(chainList, cactusMisc_nameToString(end_getName(
                    end_getOtherBlockEnd(_5End))));
            stList_append(chainList, cactusMisc_nameToString(end_getName(
                    parent5End)));
        }
    } else { //We need to introduce a new link..
        stList_append(chainList, cactusMisc_nameToString(end_getName(parent3End)));
        assert(parent5EndName == _5EndName);
        assert(end_isBlockEnd(_3End));
        stList_append(chainList, cactusMisc_nameToString(end_getName(
                end_getOtherBlockEnd(_3End))));
        addToChainList(chain, 0, chain_getLength(chain), chainList);
    }
    addToChainList(parentChain, link_getIndex(parentLink) + 1, chain_getLength(
            parentChain), chainList);
}

void getMaximalChain_extension(End *parentEnd, End *end, stList *chainList,
        int32_t orientation) {
    if (parentEnd != NULL) {
        if (end_isBlockEnd(parentEnd)) {
            End *parentOtherBlockEnd = end_getOtherBlockEnd(parentEnd);
            Link *parentOtherLink = group_getLink(end_getGroup(
                    parentOtherBlockEnd));
            if (parentOtherLink != NULL) { //We can add this chain to the final chain
                Chain *parentChain = link_getChain(parentOtherLink);
#ifdef BEN_DEBUG
                assert(link_get3End(parentOtherLink) == parentOtherBlockEnd || link_get5End(parentOtherLink) == parentOtherBlockEnd);
                if (link_get3End(parentOtherLink) == parentOtherBlockEnd) {
                    assert(link_getIndex(parentOtherLink) == 0);
                } else {
                    assert(link_getIndex(parentOtherLink) == chain_getLength(parentChain)-1);
                }
                if (orientation) {
                    assert(link_get5End(parentOtherLink) == parentOtherBlockEnd);
                } else {
                    assert(link_get3End(parentOtherLink) == parentOtherBlockEnd);
                }
#endif
                addToChainList(parentChain, 0, chain_getLength(parentChain),
                        chainList);
            }
        } else {
            assert(end_isStubEnd(parentEnd) && end_isAttached(parentEnd));
        }
    } else {
        assert(end_isBlockEnd(end));
    }
}

static stList *getMaximalChain(Chain *chain, Flower *flower, Flower *parentFlower) {
    /*
     * Calculates the structure of the promoted chain in relation to the existing parent chains.
     * This may result in the merging of the child chain with a parent chain.
     */
    //This is the list we are adding to.
    stList *chainList = stList_construct3(0, free);

    //Get the ends of the chain..
    Link *_5Link = chain_getLink(chain, 0);
    Link *_3Link = chain_getLink(chain, chain_getLength(chain) - 1);
    End *_3End = link_get3End(_5Link);
    End *_5End = link_get5End(_3Link);
    End *parent3End = flower_getEnd(parentFlower, end_getName(_3End));
    End *parent5End = flower_getEnd(parentFlower, end_getName(_5End));
    Link *parentLink;
    if (parent3End != NULL && (parentLink = group_getLink(end_getGroup(
            parent3End))) != NULL) { //The chain is within an existing chain
        getMaximalChain_between(parentLink, chain, chainList);
    } else if (parent5End != NULL && (parentLink = group_getLink(end_getGroup(
            parent5End))) != NULL) {
        getMaximalChain_between(parentLink, chain, chainList);
    } else { //The chain lies within a flower, and may extend one chain or join two seperate chains..
        getMaximalChain_extension(parent3End, _3End, chainList, 1);
        addToChainList(chain, 0, chain_getLength(chain), chainList);
        getMaximalChain_extension(parent5End, _5End, chainList, 0);
    }
    assert(stList_length(chainList) % 2 == 0);
    return chainList;
}

void chain_promote(Chain *chain) {
    /*
     * Pushes the chain into the higher level flower.
     */
    assert(chain_getLength(chain)> 0);
    Flower *flower = chain_getFlower(chain);
    Group *parentGroup = flower_getParentGroup(flower);
    assert(parentGroup != NULL);
    Flower *parentFlower = group_getFlower(parentGroup);

#ifdef BEN_DEBUG
    if (group_getLink(parentGroup) != NULL) { //Check we will be merging it into a higher level chain..
        End *_3End = link_get3End(chain_getLink(chain, 0));
        End *_5End = link_get5End(chain_getLink(chain, chain_getLength(chain)
                - 1));
        assert(end_isStubEnd(_3End) || end_isStubEnd(_5End));
    }
    flower_check(flower);
    flower_check(parentFlower);
#endif

    //Calculate the final chain structure..
    stList *finalChainList = getMaximalChain(chain, flower, parentFlower);

    //Calculate the chains that are involved in the final chain list.
    stSortedSet *chainsToExpunge = stSortedSet_construct2(
            (void(*)(void *)) chain_destruct);
    stListIterator *endIt = stList_getIterator(finalChainList);
    char *endName;
    while ((endName = stList_getNext(endIt)) != NULL) { //Get chains in the parent which we extend..
        End *end = flower_getEnd(parentFlower, cactusMisc_stringToName(endName));
        if (end != NULL) {
            Link *link = group_getLink(end_getGroup(end));
            if (link != NULL) {
                stSortedSet_insert(chainsToExpunge, link_getChain(link));
            }
        }
    }
    stList_destructIterator(endIt);
    stSortedSet_insert(chainsToExpunge, chain);

    //Handling the ends of the chain
    End *_3End = link_get3End(chain_getLink(chain, 0));
    if (end_isBlockEnd(_3End)) {
        promoteChainEnd(end_getOtherBlockEnd(_3End), flower, parentFlower);
    }
    End *_5End = link_get5End(chain_getLink(chain, chain_getLength(chain) - 1));
    if (end_isBlockEnd(_5End)) {
        promoteChainEnd(end_getOtherBlockEnd(_5End), flower, parentFlower);
    }
    //Redirect and reorganise all the ends, blocks and groups in the chain
    promoteEndsBlocksAndGroups(chain, flower, parentFlower);

    //Get rid of the old chains (at this point we have everything barring the final chain structure).
    stSortedSet_destruct(chainsToExpunge);

    //Make the final chain..
    Chain *newChain = chain_construct(parentFlower);
    for (int32_t i = 0; i < stList_length(finalChainList); i += 2) {
        Name _3EndName = cactusMisc_stringToName(stList_get(finalChainList, i));
        Name _5EndName =
                cactusMisc_stringToName(stList_get(finalChainList, i + 1));
        End *_3End = flower_getEnd(parentFlower, _3EndName);
        assert(_3End != NULL);
        End *_5End = flower_getEnd(parentFlower, _5EndName);
        assert(_5End != NULL);
        Group *group = end_getGroup(_3End);
#ifdef BEN_DEBUG
        assert(group == end_getGroup(_5End));
        assert(group_getLink(group) == NULL);
        End *end;
        Group_EndIterator *endIt = group_getEndIterator(group);
        int32_t endNumber = 0;
        while ((end = group_getNextEnd(endIt)) != NULL) {
            assert(end_getGroup(end) == group);
            if (end_isBlockEnd(end) || end_isAttached(end)) {
                endNumber++;
                assert(end == _3End || end == _5End);
            }
        }
        assert(endNumber == 2);
        group_destructEndIterator(endIt);
#endif
        link_construct(_3End, _5End, group, newChain);
    }
    stList_destruct(finalChainList);

    if (flower_getAttachedStubEndNumber(flower) == 2 && group_isTangle(parentGroup)) {
        //We've inadvertantly created a length one chain involving just the ends of the flower
        Chain *littleChain = chain_construct(parentFlower);
        End *_3End = NULL, *_5End = NULL;
        End *end;
        Group_EndIterator *endIt = group_getEndIterator(parentGroup);
        while ((end = group_getNextEnd(endIt)) != NULL) {
            if (end_isBlockEnd(end) || end_isAttached(end)) {
                assert(flower_getEnd(flower, end_getName(end)) != NULL);
                if (_3End == NULL) {
                    _3End = end;
                } else {
                    assert(_5End == NULL);
                    _5End = end;
                }
            }
        }
        group_destructEndIterator(endIt);
        link_construct(_3End, _5End, parentGroup, littleChain);
    }

#ifdef BEN_DEBUG
    assert(flower_getParentGroup(flower) == parentGroup);
    //assert(group_getLink(parentGroup) == NULL);
    if (flower_getEndNumber(flower) == 0) { //Check the properties of the flower if we've gutted it
        assert(flower_getBlockNumber(flower) == 0);
        assert(flower_getChainNumber(flower) == 0);
        assert(flower_getGroupNumber(flower) == 0);
    } else {
        assert(flower_getGroupNumber(flower)> 0);
    }
    flower_check(flower);
    flower_check(parentFlower);
#endif
}

void chain_promoteChainsThatExtendHigherLevelChains(Flower *flower) {
    Group *parentGroup = flower_getParentGroup(flower);
    if (parentGroup != NULL) {
        //Find links which need to be promoted
        Chain *chain;
        Flower_ChainIterator *chainIt = flower_getChainIterator(flower);
        stList *chains = stList_construct();
        while ((chain = flower_getNextChain(chainIt)) != NULL) {
            assert(chain_getLength(chain) >= 1);
            assert(chain_getFlower(chain) == flower);
            stList_append(chains, chain);
        }
        flower_destructChainIterator(chainIt);
        while (stList_length(chains) > 0) {
            chain = stList_pop(chains);
            End *_3End = link_get3End(chain_getLink(chain, 0));
            End *_5End = link_get5End(chain_getLink(chain, chain_getLength(
                    chain) - 1));
            if (end_isStubEnd(_3End) || end_isStubEnd(_5End)) { //Is part of higher chain..
                chain_promote(chain);
            }
        }
    }
}
