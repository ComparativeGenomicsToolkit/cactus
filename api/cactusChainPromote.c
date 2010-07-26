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
    block->blockContents->net = net;
}

static void promoteEnd(End *end, Net *net, Net *parentNet) {
    /*
     * Redirects all the pointers in the end to the higher level net.
     */
    End *parentEnd = net_getEnd(parentNet, end_getName(end));
    if (parentEnd != NULL) { //end is already present, so we will get rid of the lower level end, reconnecting the adjacencies..
        assert(end_getInstanceNumber(parentEnd) == end_getInstanceNumber(end));
        assert(end_isStubEnd(end));
        End_InstanceIterator *it = end_getInstanceIterator(end); //We ensure that the adjacencies match the higher level net..
        Cap *cap;
        while ((cap = end_getNext(it)) != NULL) {
            Cap *parentCap = end_getInstance(parentEnd, cap_getName(cap));
            assert(parentCap != NULL);
            assert(cap_getAdjacency(cap) != NULL);
            //cap_makeAdjacent(parentCap, cap_getAdjacency(cap));
            Cap *adjacentCap = cap->capContents->adjacency; //reconnect the adjacency.
            parentCap->capContents->adjacency = adjacentCap;
            adjacentCap->capContents->adjacency = parentCap;
        }
        end_setGroup(parentEnd, end_getGroup(end)); //ensures the parent end is in the group of the lower level net..
    } else { //end is not present, so we add it to the parent net..
        Cap *cap;
        End_InstanceIterator *it = end_getInstanceIterator(end);
        EventTree *eventTree = net_getEventTree(parentNet);
        while ((cap = end_getNext(it)) != NULL) {
            Event *event = eventTree_getEvent(eventTree, event_getName(
                    cap_getEvent(cap)));
            assert(event != NULL);
            cap->capContents->event = event;
            if (cap_getSequence(cap) != NULL) {
                Sequence *sequence = net_getSequence(net, sequence_getName(
                        cap_getSequence(cap)));
                assert(sequence != NULL);
                cap->capContents->sequence = sequence;
            }
            net_removeCap(net, cap);
            net_addCap(parentNet, cap);
        }
        end_destructInstanceIterator(it);
        net_removeEnd(net, end);
        net_addEnd(net, end);
        end->endContents->net = parentNet;
    }
}

static void promoteEndsBlocksAndGroups(Chain *chain, Net *net, Net *parentNet) {
    /*
     * Promotes all the ends and blocks in the chain..
     */
    Link *link;
    for (int32_t i = 0; i < chain_getLength(chain); i++) {
        link = chain_getLink(chain, i);
        End *_5End = link_get5End(link);
        if (end_isBlockEnd(_5End)) {
            promoteBlock(end_getBlock(_5End), net, parentNet);
        }
        promoteEnd(_5End, net, parentNet);
        End *_3End = link_get3End(link);
        if (i + 1 == chain_getLength(chain) && end_isBlockEnd(_3End)) {
            promoteBlock(end_getBlock(_3End), net, parentNet);
        }
        promoteEnd(link_get3End(link), net, parentNet);
        //Promote group
        Group *group = link_getGroup(link);
        net_removeGroup(net, group);
        net_addGroup(parentNet, group);
        group->net = parentNet;
    }
}

static Cap *promoteChainEndP(Cap *cap, Net *parentNet) {
    End *end = cap_getEnd(cap);
    if (end_isBlockEnd(end)) {
        Cap *otherCap = cap_getOtherSegmentCap(cap);
        assert(otherCap != NULL);
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
    assert(end_isBlockEnd(end));
    //promote the block end
    promoteEnd(end, net, parentNet);
    //Sort out the group...
    Group *childGroup = end_getGroup(end);
    end_setGroup(end, net_getParentGroup(net));

    //copy construct them into the children
    End *childEnd = end_copyConstruct(end, net);
    end_setGroup(end, childGroup);

    //copy the adjacencies in the lower level
    assert(end_getInstanceNumber(end) == end_getInstanceNumber(childEnd));
    Cap *cap;
    End_InstanceIterator *capIt = end_getInstanceIterator(end);
    while ((cap = end_getNext(capIt)) != NULL) {
        Cap *childCap = end_getInstance(childEnd, cap_getName(cap));
        Cap *childAdjacentCap = cap_getAdjacency(cap);
        assert(childCap != NULL);
        assert(childAdjacentCap != NULL);
        cap_makeAdjacent(childCap, childAdjacentCap);
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
        cap_makeAdjacent(parentCap, promoteChainEndP(adjacentCap, parentNet));
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

static void getMaximalChain_between(End *parentEnd, Link *parentLink,
        Chain *chain, stList *chainList) {
    Chain *parentChain = link_getChain(parentLink);
    addToChainList(parentChain, 0, link_getIndex(parentLink), 1, chainList);
    assert(link_get5End(parentLink) == parentEnd || link_get3End(parentLink) == parentEnd);
    bool parentChainOrientation = link_get5End(parentLink) == parentEnd;
    if (parentChainOrientation) {
        addToChainList(chain, 0, chain_getLength(chain), 1, chainList);
    } else {
        addToChainList(chain, 0, chain_getLength(chain), 0, chainList);
    }
    addToChainList(parentChain, link_getIndex(parentLink) + 1, chain_getLength(parentChain), 1, chainList);
    //chain_destruct(parentChain); //We destroy the old parent chain at this point..
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

    if (parent5End != NULL && group_getLink(end_getGroup(parent5End)) != NULL) { //The chain is within an existing chain
        Link *parent5Link = group_getLink(end_getGroup(parent5End));
        assert(parent3End != NULL);
        Link *parent3Link = group_getLink(end_getGroup(parent3End));
        assert(parent3Link != NULL);
        assert(link_getChain(parent5Link) == link_getChain(parent3Link));
        //Now to the actual work..
        getMaximalChain_between(parent5End, parent5Link, chain, chainList);
    } else { //The chain lies within a net, and may extend one chain or join two seperate chains..
        if (parent5End != NULL) {
            if(end_isBlockEnd(parent5End)) {
                End *parentOtherBlockEnd = end_getOtherBlockEnd(parent5End);
                Link *parentOtherLink = group_getLink(end_getGroup(
                        parentOtherBlockEnd));
                if (parentOtherLink != NULL) { //We can add this chain to the final chain
                    Chain *parentChain = link_getChain(parentOtherLink);
                    assert(link_get5End(parentOtherLink) == parentOtherBlockEnd || link_get3End(parentOtherLink) == parentOtherBlockEnd);
                    addToChainList(parentChain, 0, chain_getLength(parentChain),
                            link_get5End(parentOtherLink) != parentOtherBlockEnd,
                            chainList);
                }
            }
            else {
                assert(end_isStubEnd(parent5End) && end_isAttached(parent5End));
            }
        }
        else {
            assert(end_isBlockEnd(_5End));
        }
        addToChainList(chain, 0, chain_getLength(chain), 0, chainList);
        if (parent3End != NULL) {
            assert(group_getLink(end_getGroup(parent3End)) == NULL); //Check it is not in a link..
            if(end_isBlockEnd(parent3End)) {
                End *parentOtherBlockEnd = end_getOtherBlockEnd(parent3End);
                Link *parentOtherLink = group_getLink(end_getGroup(
                        parentOtherBlockEnd));
                if (parentOtherLink != NULL) {
                    Chain *parentChain = link_getChain(parentOtherLink);
                    assert(link_get5End(parentOtherLink) == parentOtherBlockEnd || link_get3End(parentOtherLink) == parentOtherBlockEnd);
                    addToChainList(parentChain, 0, chain_getLength(parentChain),
                            link_get5End(parentOtherLink) == parentOtherBlockEnd,
                            chainList);
                }
            }
            else {
                assert(end_isStubEnd(parent3End) && end_isAttached(parent3End));
            }
        }
        else {
            assert(end_isBlockEnd(_3End));
        }
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

    //Calculate the final chain structure..
    stList *finalChainList = getMaximalChain(chain, net, parentNet);

    return;

    //Redirect and reorganise all the ends, blocks and groups in the chain
    promoteEndsBlocksAndGroups(chain, net, parentNet);
    //Handling the ends of the chain
    End *_5End = link_get5End(chain_getLink(chain, 0));
    if (end_isBlockEnd(_5End)) {
        promoteChainEnd(end_getOtherBlockEnd(_5End), net, parentNet);
    }
    End *_3End = link_get3End(chain_getLink(chain, chain_getLength(chain) - 1));
    if (end_isBlockEnd(_3End)) {
        promoteChainEnd(end_getOtherBlockEnd(_3End), net, parentNet);
    }
    //Get rid of the lower level chain..
    chain_destruct(chain);

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
        assert(end_getGroup(_5End) == end_getGroup(_3End));
        link_construct(_5End, _3End, end_getGroup(_5End), newChain);
    }
    stList_destruct(finalChainList);

    //Destroy the net if it is now empty..
    if (net_getEndNumber(net) == 0) {
        assert(net_getGroupNumber(net) == 0);
        netDisk_deleteNetFromDisk(net_getNetDisk(net), net_getName(net));
        net_destruct(net, 0);
    }
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
            stList_append(chains, chain);
        }
        net_destructChainIterator(chainIt);
        while (stList_length(chains) > 0) {
            chain = stList_pop(chains);
            End *_5End = link_get5End(chain_getLink(chain, 0));
            End *_3End = link_get3End(chain_getLink(chain, chain_getLength(
                    chain) - 1));
            if (end_isStubEnd(_5End) || end_isStubEnd(_3End)) {
                chain_promote(chain);
            }
        }
    }
}
