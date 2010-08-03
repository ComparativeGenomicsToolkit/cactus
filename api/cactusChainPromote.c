#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions used to normalise the structure of the cactus graph.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static void promoteBlock(Block *block, Net *net, Net *parentNet) {
    /*
     * Redirects all the pointers in the block to the higher level.
     */
    Segment *segment;
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    while ((segment = block_getNext(it)) != NULL) {
        net_removeSegment(net, segment);
        net_addSegment(parentNet, segment);
    }
    block_destructInstanceIterator(it);
    net_removeBlock(net, block);
    net_addBlock(parentNet, block);
    block->blockContents->net = parentNet;
}

static void mergeStubEnd(End *end, Net *net, Net *parentNet) {
    /*
     * Removes the lower level stub end, reconnects the adjacencies of the higher level end to that of the lower level stub end.
     */
    assert(end_isStubEnd(end));
    assert(group_getLink(end_getGroup(end)) != NULL);
    End *parentEnd = net_getEnd(parentNet, end_getName(end));
    assert(parentEnd != NULL); //end is already present, so we will get rid of the lower level end, reconnecting the adjacencies..
    assert(end_getInstanceNumber(parentEnd) == end_getInstanceNumber(end));

    Group *group = end_getGroup(end);
    End_InstanceIterator *it = end_getInstanceIterator(end); //We ensure that the adjacencies match the higher level net..
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        Cap *parentCap = end_getInstance(parentEnd, cap_getName(cap));
        assert(parentCap != NULL);
        Cap *adjacentCap = cap_getAdjacency(cap);
        if (cap_getAdjacency(cap) != NULL) {
            Cap *parentAdjacentCap = net_getCap(parentNet, cap_getName(
                    adjacentCap));
            assert(parentAdjacentCap != NULL);
            assert(netDisk_getNet(net_getNetDisk(net), net_getName(end_getNet(cap_getEnd(adjacentCap)))) == end_getNet(cap_getEnd(adjacentCap)));
            cap_breakAdjacency(cap); //We must remove the old adjacency so it doesn't hang around.
            cap_makeAdjacent(parentCap, parentAdjacentCap);
        }
    }
    end_destruct(end); //Destruct the old end
    assert(group_getNet(group) == parentNet); //the group should already have been promoted
    end_setGroup(parentEnd, group); //ensures the parent end is in the group of the lower level net..
}

static void promoteBlockEnd(End *end, Net *net, Net *parentNet) {
    /*
     * Redirects all the pointers in the block end to the higher level net.
     */
    assert(end_isBlockEnd(end));
    assert(net_getEnd(parentNet, end_getName(end)) == NULL);
    Cap *cap;
    End_InstanceIterator *it = end_getInstanceIterator(end);
    EventTree *eventTree = net_getEventTree(parentNet);
    while ((cap = end_getNext(it)) != NULL) {
        Event *event = eventTree_getEvent(eventTree, event_getName(
                cap_getEvent(cap)));
        assert(event != NULL);
        cap->capContents->event = event;
        if (cap_getSequence(cap) != NULL) {
            Sequence *sequence = net_getSequence(parentNet, sequence_getName(
                    cap_getSequence(cap)));
            assert(sequence != NULL);
            cap->capContents->sequence = sequence;
        }
        net_removeCap(net, cap);
        net_addCap(parentNet, cap);
    }
    end_destructInstanceIterator(it);
    net_removeEnd(net, end);
    net_addEnd(parentNet, end);
    end->endContents->net = parentNet;
}

static void promoteEndsBlocksAndGroups(Chain *chain, Net *net, Net *parentNet) {
    /*
     * Promotes all the ends and blocks in the chain..
     */
    for (int32_t i = 0; i < chain_getLength(chain); i++) {
        Link *link = chain_getLink(chain, i);
        End *_5End = link_get5End(link);
        End *_3End = link_get3End(link);
        if (end_isBlockEnd(_5End)) {
            promoteBlock(end_getBlock(_5End), net, parentNet);
            promoteBlockEnd(_5End, net, parentNet);
        } else {
            assert(i == 0);
        }
        if (end_isBlockEnd(_3End)) {
            if (i + 1 == chain_getLength(chain)) {
                promoteBlock(end_getBlock(_3End), net, parentNet);
            }
            promoteBlockEnd(_3End, net, parentNet);
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
        net_removeGroup(net, group);
        net_addGroup(parentNet, group);
        group->net = parentNet;
        Net *nestedNet = group_getNestedNet(group);
        if (nestedNet != NULL) {
            nestedNet->parentNetName = net_getName(parentNet);
        }
        //Promote any free stub ends..
        while (stList_length(freeStubEndsToPromote) > 0) {
            mergeStubEnd(stList_pop(freeStubEndsToPromote), net, parentNet);
        }
        stList_destruct(freeStubEndsToPromote);
    }
    //Now do the stub ends of the chain..
    End *_5End = link_get5End(chain_getLink(chain, 0));
    if (end_isStubEnd(_5End)) {
        assert(end_isAttached(_5End));
        assert(end_getNet(_5End) == net);
        mergeStubEnd(_5End, net, parentNet);
    }
    End *_3End = link_get3End(chain_getLink(chain, chain_getLength(chain) - 1));
    if (end_isStubEnd(_3End)) {
        assert(end_isAttached(_3End));
        assert(end_getNet(_3End) == net);
        mergeStubEnd(_3End, net, parentNet);
    }
}

static Cap *promoteChainEndP(Cap *cap, Net *parentNet) {
    End *end = cap_getEnd(cap);
    if (end_isBlockEnd(end)) {
        assert(net_getEnd(parentNet, end_getName(end)) == NULL);
        Cap *otherCap = cap_getOtherSegmentCap(cap);
        assert(otherCap != NULL);
        assert(cap != otherCap);
        assert(cap_getEnd(otherCap) != cap_getEnd(cap));
        Cap *adjacentCap = cap_getAdjacency(otherCap);
        assert(adjacentCap != NULL);
        return promoteChainEndP(adjacentCap, parentNet);
    }
    Cap *higherCap = net_getCap(parentNet, cap_getName(cap));
    assert(higherCap != NULL);
    assert(higherCap != cap);
    assert(cap_getName(higherCap) == cap_getName(cap));
    return higherCap;
}

static void promoteChainEnd(End *end, Net *net, Net *parentNet) {
    /*
     * Promotes the block end of a chain to the parent net, creating a stub end
     * in the lower level net.
     */
    assert(end_isBlockEnd(end));
    assert(end_getOrientation(end));
    assert(end_getNet(end) == net);
    //promote the block end
    promoteBlockEnd(end, net, parentNet);
    //Sort out the group...
    Group *childGroup = end_getGroup(end);
    end_setGroup(end, net_getParentGroup(net));

    //copy construct them into the children
    End *childEnd = end_copyConstruct(end, net);
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
        assert(end_getNet(cap_getEnd(childCap)) == net);
        assert(end_getNet(cap_getEnd(childAdjacentCap)) == net);
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
        assert(end_getNet(cap_getEnd(adjacentCap)) == net);
        //Search for the higher level adjacency..
        Cap *adjacentParentCap = promoteChainEndP(adjacentCap, parentNet);
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
        int32_t orientation, stList *chainList) {
    /*
     * Adds the links in the chain between the specified start (inc) and end (exc) indices to the list
     * chainList.
     */
    if (orientation) {
        for (int32_t i = start; i < end; i++) {
            Link *link = chain_getLink(chain, i);
            stList_append(chainList, netMisc_nameToString(end_getName(
                    link_get5End(link))));
            stList_append(chainList, netMisc_nameToString(end_getName(
                    link_get3End(link))));
        }
    } else {
        for (int32_t i = end - 1; i >= start; i--) {
            Link *link = chain_getLink(chain, i);
            stList_append(chainList, netMisc_nameToString(end_getName(
                    link_get3End(link))));
            stList_append(chainList, netMisc_nameToString(end_getName(
                    link_get5End(link))));
        }
    }
}

static void getMaximalChain_between(Link *parentLink, Chain *chain,
        stList *chainList) {
    Chain *parentChain = link_getChain(parentLink);
    addToChainList(parentChain, 0, link_getIndex(parentLink), 1, chainList);
    //Get the ends..
    End *parent5End = link_get5End(parentLink);
    End *parent3End = link_get3End(parentLink);
    End *_5End = link_get5End(chain_getLink(chain, 0));
    End *_3End = link_get3End(chain_getLink(chain, chain_getLength(chain) - 1));
    Name parent5EndName = end_getName(parent5End);
    Name parent3EndName = end_getName(parent3End);
    Name _5EndName = end_getName(_5End);
    Name _3EndName = end_getName(_3End);
    if (parent5EndName == _5EndName) {
        addToChainList(chain, 0, chain_getLength(chain), 1, chainList);
        if (parent3EndName != _3EndName) { //We need to introduce a link between the other block end and this link..
            assert(end_isBlockEnd(_3End));
            stList_append(chainList, netMisc_nameToString(end_getName(
                    end_getOtherBlockEnd(_3End))));
            stList_append(chainList, netMisc_nameToString(end_getName(
                    parent3End)));
        }
    } else if (parent5EndName == _3EndName) {
        addToChainList(chain, 0, chain_getLength(chain), 0, chainList);
        if (parent3EndName != _5EndName) { //We need to introduce a link between the other block end and this link..
            assert(end_isBlockEnd(_5End));
            stList_append(chainList, netMisc_nameToString(end_getName(
                    end_getOtherBlockEnd(_5End))));
            stList_append(chainList, netMisc_nameToString(end_getName(
                    parent3End)));
        }
    } else { //We need to introduce a new link..
        stList_append(chainList, netMisc_nameToString(end_getName(parent5End)));
        if (parent3EndName == _3EndName) {
            assert(end_isBlockEnd(_5End));
            stList_append(chainList, netMisc_nameToString(end_getName(
                    end_getOtherBlockEnd(_5End))));
            addToChainList(chain, 0, chain_getLength(chain), 1, chainList);
        } else {
            assert(parent3EndName == _5EndName);
            assert(end_isBlockEnd(_3End));
            stList_append(chainList, netMisc_nameToString(end_getName(
                    end_getOtherBlockEnd(_3End))));
            addToChainList(chain, 0, chain_getLength(chain), 0, chainList);
        }
    }
    addToChainList(parentChain, link_getIndex(parentLink) + 1, chain_getLength(
            parentChain), 1, chainList);
}

void getMaximalChain_extension(End *parentEnd, End *end, stList *chainList, int32_t orientation) {
    if (parentEnd != NULL) {
        if (end_isBlockEnd(parentEnd)) {
            End *parentOtherBlockEnd = end_getOtherBlockEnd(parentEnd);
            Link *parentOtherLink = group_getLink(end_getGroup(
                    parentOtherBlockEnd));
            if (parentOtherLink != NULL) { //We can add this chain to the final chain
                Chain *parentChain = link_getChain(parentOtherLink);
#ifdef BEN_DEBUG
                assert(link_get5End(parentOtherLink) == parentOtherBlockEnd || link_get3End(parentOtherLink) == parentOtherBlockEnd);
                if (link_get5End(parentOtherLink) == parentOtherBlockEnd) {
                    assert(link_getIndex(parentOtherLink) == 0);
                } else {
                    assert(link_getIndex(parentOtherLink) == chain_getLength(parentChain)-1);
                }
#endif
                addToChainList(parentChain, 0, chain_getLength(parentChain),
                        (orientation ? link_get3End(parentOtherLink) == parentOtherBlockEnd : link_get5End(parentOtherLink) == parentOtherBlockEnd),
                        chainList);
            }
        } else {
            assert(end_isStubEnd(parentEnd) && end_isAttached(parentEnd));
        }
    } else {
        assert(end_isBlockEnd(end));
    }
}

static stList *getMaximalChain(Chain *chain, Net *net, Net *parentNet) {
    /*
     * Calculates the structure of the promoted chain in relation to the existing parent chains.
     * This may result in the merging of the child chain with a parent chain.
     */
    //This is the list we are adding to.
    stList *chainList = stList_construct3(0, free);

    //Get the ends of the chain..
    Link *_5Link = chain_getLink(chain, 0);
    Link *_3Link = chain_getLink(chain, chain_getLength(chain) - 1);
    End *_5End = link_get5End(_5Link);
    End *_3End = link_get3End(_3Link);
    End *parent5End = net_getEnd(parentNet, end_getName(_5End));
    End *parent3End = net_getEnd(parentNet, end_getName(_3End));
    Link *parentLink;
    if (parent5End != NULL && (parentLink = group_getLink(end_getGroup(
            parent5End))) != NULL) { //The chain is within an existing chain
        getMaximalChain_between(parentLink, chain, chainList);
    } else if (parent3End != NULL && (parentLink = group_getLink(end_getGroup(
            parent3End))) != NULL) {
        getMaximalChain_between(parentLink, chain, chainList);
    } else { //The chain lies within a net, and may extend one chain or join two seperate chains..
        getMaximalChain_extension(parent5End, _5End, chainList, 1);
        addToChainList(chain, 0, chain_getLength(chain), 1, chainList);
        getMaximalChain_extension(parent3End, _3End, chainList, 0);
    }
    assert(stList_length(chainList) % 2 == 0);
    return chainList;
}

void chain_promote(Chain *chain) {
    /*
     * Pushes the chain into the higher level net.
     */
    assert(chain_getLength(chain)> 0);
    Net *net = chain_getNet(chain);
    Group *parentGroup = net_getParentGroup(net);
    assert(parentGroup != NULL);
    Net *parentNet = group_getNet(parentGroup);

#ifdef BEN_DEBUG
    if (group_getLink(parentGroup) != NULL) { //Check we will be merging it into a higher level chain..
        End *_5End = link_get5End(chain_getLink(chain, 0));
        End *_3End = link_get3End(chain_getLink(chain, chain_getLength(chain)
                - 1));
        assert(end_isStubEnd(_5End) || end_isStubEnd(_3End));
    }
    net_check(net);
    net_check(parentNet);
#endif

    //Calculate the final chain structure..
    stList *finalChainList = getMaximalChain(chain, net, parentNet);

    //Calculate the chains that are involved in the final chain list.
    stSortedSet *chainsToExpunge = stSortedSet_construct2(
            (void(*)(void *)) chain_destruct);
    stListIterator *endIt = stList_getIterator(finalChainList);
    char *endName;
    while ((endName = stList_getNext(endIt)) != NULL) { //Get chains in the parent which we extend..
        End *end = net_getEnd(parentNet, netMisc_stringToName(endName));
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
    End *_5End = link_get5End(chain_getLink(chain, 0));
    if (end_isBlockEnd(_5End)) {
        promoteChainEnd(end_getOtherBlockEnd(_5End), net, parentNet);
    }
    End *_3End = link_get3End(chain_getLink(chain, chain_getLength(chain) - 1));
    if (end_isBlockEnd(_3End)) {
        promoteChainEnd(end_getOtherBlockEnd(_3End), net, parentNet);
    }
    //Redirect and reorganise all the ends, blocks and groups in the chain
    promoteEndsBlocksAndGroups(chain, net, parentNet);

    //Get rid of the old chains (at this point we have everything barring the final chain structure).
    stSortedSet_destruct(chainsToExpunge);

    //Make the final chain..
    Chain *newChain = chain_construct(parentNet);
    for (int32_t i = 0; i < stList_length(finalChainList); i += 2) {
        Name _5EndName = netMisc_stringToName(stList_get(finalChainList, i));
        Name _3EndName =
                netMisc_stringToName(stList_get(finalChainList, i + 1));
        End *_5End = net_getEnd(parentNet, _5EndName);
        assert(_5End != NULL);
        End *_3End = net_getEnd(parentNet, _3EndName);
        assert(_3End != NULL);
        Group *group = end_getGroup(_5End);
#ifdef BEN_DEBUG
        assert(group == end_getGroup(_3End));
        assert(group_getLink(group) == NULL);
        End *end;
        Group_EndIterator *endIt = group_getEndIterator(group);
        int32_t endNumber = 0;
        while((end = group_getNextEnd(endIt)) != NULL) {
            assert(end_getGroup(end) == group);
            if(end_isBlockEnd(end) || end_isAttached(end)) {
                endNumber++;
                assert(end == _5End || end == _3End);
            }
        }
        assert(endNumber == 2);
        group_destructEndIterator(endIt);
#endif
        link_construct(_5End, _3End, group, newChain);
    }
    stList_destruct(finalChainList);

    if(net_getAttachedStubEndNumber(net) == 2 && group_isTangle(parentGroup)) {
        //We've inadvertantly created a length one chain involving just the ends of the net
        Chain *littleChain = chain_construct(parentNet);
        End *_5End = NULL, *_3End = NULL;
        End *end;
        Group_EndIterator *endIt = group_getEndIterator(parentGroup);
        while((end = group_getNextEnd(endIt)) != NULL) {
            if(end_isBlockEnd(end) || end_isStubEnd(end)) {
                assert(net_getEnd(net, end_getName(end)) != NULL);
                if(_5End == NULL) {
                    _5End = end;
                }
                else {
                    assert(_3End == NULL);
                    _3End = end;
                }
            }
        }
        group_destructEndIterator(endIt);
        link_construct(_5End, _3End, parentGroup, littleChain);
    }

#ifdef BEN_DEBUG
    assert(net_getParentGroup(net) == parentGroup);
    //assert(group_getLink(parentGroup) == NULL);
    if (net_getEndNumber(net) == 0) { //Check the properties of the net if we've gutted it
        assert(net_getBlockNumber(net) == 0);
        assert(net_getChainNumber(net) == 0);
        assert(net_getGroupNumber(net) == 0);
    } else {
        assert(net_getGroupNumber(net)> 0);
    }
    net_check(net);
    net_check(parentNet);
#endif
}

void chain_promoteChainsThatExtendHigherLevelChains(Net *net) {
    Group *parentGroup = net_getParentGroup(net);
    if (parentGroup != NULL) {
        //Find links which need to be promoted
        Chain *chain;
        Net_ChainIterator *chainIt = net_getChainIterator(net);
        stList *chains = stList_construct();
        while ((chain = net_getNextChain(chainIt)) != NULL) {
            assert(chain_getLength(chain) >= 1);
            assert(chain_getNet(chain) == net);
            stList_append(chains, chain);
        }
        net_destructChainIterator(chainIt);
        while (stList_length(chains) > 0) {
            chain = stList_pop(chains);
            End *_5End = link_get5End(chain_getLink(chain, 0));
            End *_3End = link_get3End(chain_getLink(chain, chain_getLength(
                    chain) - 1));
            if (end_isStubEnd(_5End) || end_isStubEnd(_3End)) { //Is part of higher chain..
                chain_promote(chain);
            }
        }
    }
}

