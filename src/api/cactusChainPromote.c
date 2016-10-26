/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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
    assert(end_getOrientation(end));
    assert(group_getLink(end_getGroup(end)) != NULL);
    End *parentEnd = flower_getEnd(parentFlower, end_getName(end));
    assert(parentEnd != NULL); //end is already present, so we will get rid of the lower level end, reconnecting the adjacencies..
    assert(end_getInstanceNumber(parentEnd) == end_getInstanceNumber(end));

    Group *group = end_getGroup(end);
    if(end_isAttached(end)) {
        Link *link = group_getLink(group);
        if(link != NULL) {
            if(end_getSide(end)) {
                assert(link_get5End(link) == end);
                link->_5End = parentEnd;
            }
            else {
                assert(link_get3End(link) == end);
                link->_3End = parentEnd;
            }
        }
    }

    End_InstanceIterator *it = end_getInstanceIterator(end); //We ensure that the adjacencies match the higher level flower..
    Cap *cap;
    while ((cap = end_getNext(it)) != NULL) {
        Cap *parentCap = end_getInstance(parentEnd, cap_getName(cap));
        assert(parentCap != NULL);
        Cap *adjacentCap = cap_getAdjacency(cap);
        if (cap_getAdjacency(cap) != NULL) {
            Cap *parentAdjacentCap = flower_getCap(parentFlower, cap_getName(adjacentCap));
            assert(parentAdjacentCap != NULL);
            assert(
                    cactusDisk_getFlower(flower_getCactusDisk(flower),
                            flower_getName(end_getFlower(cap_getEnd(adjacentCap)))) == end_getFlower(
                            cap_getEnd(adjacentCap)));
            cap_breakAdjacency(cap); //We must remove the old adjacency so it doesn't hang around.
            cap_makeAdjacent(parentCap, parentAdjacentCap);
        }
    }
    end_destructInstanceIterator(it);
    end_destruct(end); //Destruct the old end
    assert(group_getFlower(group) == parentFlower); //the group should already have been promoted
    end_setGroup(parentEnd, group); //ensures the parent end is in the group of the lower level flower..
}

static void promoteBlockEnd(End *end, Flower *flower, Flower *parentFlower) {
    /*
     * Redirects all the pointers in the block end to the higher level flower.
     */
    assert(end_isBlockEnd(end));
    assert(end_getFlower(end) == flower);
    assert(flower != parentFlower);
    assert(flower_getEnd(parentFlower, end_getName(end)) == NULL);
    Cap *cap;
    End_InstanceIterator *it = end_getInstanceIterator(end);
    EventTree *eventTree = flower_getEventTree(parentFlower);
    while ((cap = end_getNext(it)) != NULL) {
        Event *event = eventTree_getEvent(eventTree, event_getName(cap_getEvent(cap)));
        assert(event != NULL);
        cap->capContents->event = event;
        if (cap_getSequence(cap) != NULL) {
            Sequence *sequence = flower_getSequence(parentFlower, sequence_getName(cap_getSequence(cap)));
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
    assert(chain_getFlower(chain) == flower);
    bool circular = 0;
    Link *link = chain_getFirst(chain);
    while (link != NULL) {
        End *_3End = link_get3End(link);
        End *_5End = link_get5End(link);
        if (end_isBlockEnd(_3End)) {
            if (chain_getLength(chain) == 1 && end_getOtherBlockEnd(_3End) == _5End) { //Decide if circular or not
                circular = 1;
            }
            promoteBlock(end_getBlock(_3End), flower, parentFlower);
            promoteBlockEnd(_3End, flower, parentFlower);
        } else {
            assert(link == chain_getFirst(chain));
        }
        if (end_isBlockEnd(_5End)) {
            if (link == chain_getLast(chain) && !circular) {
                promoteBlock(end_getBlock(_5End), flower, parentFlower);
            }
            promoteBlockEnd(_5End, flower, parentFlower);
        } else {
            assert(link == chain_getLast(chain));
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
        link = link_getNextLink(link);
    }
    //Now do the stub ends of the chain..
    End *_3End = link_get3End(chain_getFirst(chain));
    if (end_isStubEnd(_3End)) {
        assert(end_isAttached(_3End));
        assert(end_getFlower(_3End) == flower);
        mergeStubEnd(_3End, flower, parentFlower);
    }
    End *_5End = link_get5End(chain_getLast(chain));
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

static void promoteChainEnds(stList *ends, Flower *flower, Flower *parentFlower) {
    /*
     * Promotes the block end of a chain to the parent flower, creating a stub end
     * in the lower level flower.
     */
    stListIterator *endIt = stList_getIterator(ends);
    End *end;
    stHash *childEnds = stHash_construct();
    while ((end = stList_getNext(endIt)) != NULL) {
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
        stHash_insert(childEnds, end, childEnd);
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
            childAdjacentCap = flower_getCap(flower, cap_getName(childAdjacentCap)); //This is necessary in the case that there was a self loop!
            assert(childAdjacentCap != NULL);
            assert(end_getFlower(cap_getEnd(childCap)) == flower);
            assert(end_getFlower(cap_getEnd(childAdjacentCap)) == flower);
            assert(end_getGroup(cap_getEnd(childCap)) == end_getGroup(cap_getEnd(childAdjacentCap)));
            cap_makeAdjacent(childCap, childAdjacentCap); //reconnect the adjacency.
        }
        end_destructInstanceIterator(capIt);
    }
    stList_destructIterator(endIt);

    endIt = stList_getIterator(ends);
    while ((end = stList_getNext(endIt)) != NULL) {
        //create the attachments in the parent
        End *childEnd = stHash_search(childEnds, end);
        Cap *cap;
        End_InstanceIterator *capIt = end_getInstanceIterator(childEnd);
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
    stList_destructIterator(endIt);
    stHash_destruct(childEnds);
}

/////////
//Functions which calculate the chain structure after the promotion of a chain.
////////

static void addToChainList(Chain *chain, Link *linkFirst, Link *linkEnd, stList *chainList) {
    /*
     * Adds the links in the chain between the specified start (inc) and end (exc) indices to the list
     * chainList.
     */
    Link *link = linkFirst;
    while (link != linkEnd) {
        stList_append(chainList, stIntTuple_construct1( end_getName(link_get3End(link))));
        stList_append(chainList, stIntTuple_construct1( end_getName(link_get5End(link))));
        link = link_getNextLink(link);
    }
}

static void getMaximalChain_between(Link *parentLink, Chain *chain, stList *chainList) {
    Chain *parentChain = link_getChain(parentLink);
    addToChainList(parentChain, chain_getFirst(parentChain), parentLink, chainList);
    //Get the ends..
    End *parent3End = link_get3End(parentLink);
    End *parent5End = link_get5End(parentLink);
    End *_3End = link_get3End(chain_getFirst(chain));
    End *_5End = link_get5End(chain_getLast(chain));
    Name parent3EndName = end_getName(parent3End);
    Name parent5EndName = end_getName(parent5End);
    Name _3EndName = end_getName(_3End);
    Name _5EndName = end_getName(_5End);
    assert(parent3EndName != _5EndName);
    if (parent3EndName == _3EndName) {
        addToChainList(chain, chain_getFirst(chain), NULL, chainList);
        if (parent5EndName != _5EndName) { //We need to introduce a link between the other block end and this link..
            assert(end_isBlockEnd(_5End));
            stList_append(chainList, stIntTuple_construct1( end_getName(end_getOtherBlockEnd(_5End))));
            stList_append(chainList, stIntTuple_construct1( end_getName(parent5End)));
        }
    } else { //We need to introduce a new link..
        stList_append(chainList, stIntTuple_construct1( end_getName(parent3End)));
        assert(parent5EndName == _5EndName);
        assert(end_isBlockEnd(_3End));
        stList_append(chainList, stIntTuple_construct1( end_getName(end_getOtherBlockEnd(_3End))));
        addToChainList(chain, chain_getFirst(chain), NULL, chainList);
    }
    addToChainList(parentChain, link_getNextLink(parentLink), NULL, chainList);
}

void getMaximalChain_extension(End *parentEnd, End *end, stList *chainList, int64_t orientation) {
    if (parentEnd != NULL) {
        if (end_isBlockEnd(parentEnd)) {
            End *parentOtherBlockEnd = end_getOtherBlockEnd(parentEnd);
            Link *parentOtherLink = group_getLink(end_getGroup(parentOtherBlockEnd));
            if (parentOtherLink != NULL) { //We can add this chain to the final chain
                Chain *parentChain = link_getChain(parentOtherLink);
#ifndef NDEBUG
                assert(link_get3End(parentOtherLink) == parentOtherBlockEnd || link_get5End(parentOtherLink) == parentOtherBlockEnd);
                if (link_get3End(parentOtherLink) == parentOtherBlockEnd) {
                    assert(chain_getFirst(parentChain) == parentOtherLink);
                } else {
                    assert(chain_getLast(parentChain) == parentOtherLink);
                }
                if (orientation) {
                    assert(link_get5End(parentOtherLink) == parentOtherBlockEnd);
                } else {
                    assert(link_get3End(parentOtherLink) == parentOtherBlockEnd);
                }
#endif
                addToChainList(parentChain, chain_getFirst(parentChain), NULL, chainList);
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
    Link *_5Link = chain_getFirst(chain);
    Link *_3Link = chain_getLast(chain);
    End *_3End = link_get3End(_5Link);
    End *_5End = link_get5End(_3Link);
    End *parent3End = flower_getEnd(parentFlower, end_getName(_3End));
    End *parent5End = flower_getEnd(parentFlower, end_getName(_5End));
    Link *parentLink;
    if (parent3End != NULL && (parentLink = group_getLink(end_getGroup(parent3End))) != NULL) { //The chain is within an existing chain
        getMaximalChain_between(parentLink, chain, chainList);
    } else if (parent5End != NULL && (parentLink = group_getLink(end_getGroup(parent5End))) != NULL) {
        getMaximalChain_between(parentLink, chain, chainList);
    } else { //The chain lies within a flower, and may extend one chain or join two separate chains..
        getMaximalChain_extension(parent3End, _3End, chainList, 1);
        addToChainList(chain, chain_getFirst(chain), NULL, chainList);
        getMaximalChain_extension(parent5End, _5End, chainList, 0);
    }
    assert(stList_length(chainList) % 2 == 0);
    return chainList;
}

void block_promote(Block *block) {
    Flower *flower = block_getFlower(block);
    Group *parentGroup = flower_getParentGroup(flower);
    assert(parentGroup != NULL);
    Flower *parentFlower = group_getFlower(parentGroup);
    assert(group_getLink(parentGroup) == NULL);
    assert(group_getLink(end_getGroup(block_get5End(block))) == NULL);
    assert(group_getLink(end_getGroup(block_get3End(block))) == NULL);

    //Handling the ends of the chain..
    stList *chainEnds = stList_construct();
    End *_3End = block_get3End(block);
    if (end_isBlockEnd(_3End)) {
        stList_append(chainEnds, end_getOtherBlockEnd(_3End));
    }
    End *_5End = block_get5End(block);
    if (end_isBlockEnd(_5End)) {
        stList_append(chainEnds, end_getOtherBlockEnd(_5End));
    }
    promoteChainEnds(chainEnds, flower, parentFlower);
    stList_destruct(chainEnds);

    //Now finally promote the block.
    promoteBlock(block, flower, parentFlower);

    //We've inadvertantly created a length one chain involving just the ends of the flower
    group_constructChainForLink(parentGroup);

    //We have not removed any ends.. so we're done.
    assert(flower_getParentGroup(flower) == parentGroup);
}

void chain_promote(Chain *chain) {
    /*
     * Pushes the chain into the higher level flower.
     */

    Flower *flower = chain_getFlower(chain);
    Group *parentGroup = flower_getParentGroup(flower);
    assert(chain_getLength(chain)> 0);
    assert(parentGroup != NULL);
    Flower *parentFlower = group_getFlower(parentGroup);

#ifndef NDEBUG
    if (group_getLink(parentGroup) != NULL) { //Check we will be merging it into a higher level chain..
        End *_3End = link_get3End(chain_getFirst(chain));
        End *_5End = link_get5End(chain_getLast(chain));
        assert(end_isStubEnd(_3End) || end_isStubEnd(_5End));
    }
#endif

    //if(chain_getLength(chain) == 1) {
    Link *_5Link = chain_getFirst(chain);
    Link *_3Link = chain_getLast(chain);
    if (end_isStubEnd(link_get3End(_5Link)) && end_isStubEnd(link_get5End(_3Link))) {
        End *_5ChildEnd = link_get3End(_5Link);
        End *_3ChildEnd = link_get5End(_3Link);
        End *_5End = flower_getEnd(parentFlower, end_getName(_5ChildEnd)); //Remember the reversal of sides here
        End *_3End = flower_getEnd(parentFlower, end_getName(_3ChildEnd));

        //Various completely paranoid checks
        assert(_5End != NULL);
        assert(_3End != NULL);
        assert(end_getName(_5ChildEnd) == end_getName(_5End));
        assert(end_getName(_3ChildEnd) == end_getName(_3End));
        assert(end_getSide(_5ChildEnd) == end_getSide(_5End));
        assert(end_getSide(_3ChildEnd) == end_getSide(_3End));
        assert(end_getSide(_5End) != end_getSide(_3End));
        assert(end_getSide(_3ChildEnd));
        assert(!end_getSide(_5ChildEnd));
        assert(end_getSide(_3End));
        assert(!end_getSide(_5End));

        if (end_getGroup(_5End) == end_getGroup(_3End) && group_isLink(end_getGroup(_5End))) { //We have a trivial promotion, should cover most cases
            assert(end_getGroup(_5End) == parentGroup);
            assert(end_getGroup(_3End) == parentGroup);
            Link *parentLink = group_getLink(end_getGroup(_5End));
            assert(group_getLink(end_getGroup(_3End)) == parentLink);

            assert(parentLink != NULL);
            assert(link_get3End(parentLink) == _5End);
            assert(link_get5End(parentLink) == _3End);
            Chain *parentChain = link_getChain(parentLink);
            Link *pLink = link_getPreviousLink(parentLink);
            Link *nLink = link_getNextLink(parentLink);
            //Excise the parent link
            int64_t finalLength = chain_getLength(chain) + chain_getLength(parentChain) - 1;
            (void)finalLength;
            if (pLink != NULL) {
                pLink->nLink = nLink;
                if (nLink != NULL) { //Fix the reverse link
                    nLink->pLink = pLink;
                }
            } else {
                parentChain->link = nLink;
                if (nLink != NULL) {
                    nLink->pLink = NULL;
                }
            }
            parentChain->linkNumber--;
            assert(parentChain->linkNumber >= 0);

            //Cleanup old parent link
            free(parentLink); //Avoiding destructor because of hierarchy of side effects.
            //Now cleanup parent group
            parentGroup->link = NULL;

            //Promote the chain and insert it into the existing chain
            promoteEndsBlocksAndGroups(chain, flower, parentFlower);
            assert(group_getEndNumber(parentGroup) == 0);

            assert(_5Link->pLink == NULL);
            while (_5Link != NULL) {
                Link *linkToPromote = _5Link;
                _5Link = _5Link->nLink; //Setup next loop
                //Fix the linked list
                if (pLink != NULL) {
                    pLink->nLink = linkToPromote;
                } else {
                    parentChain->link = linkToPromote;
                }
                linkToPromote->pLink = pLink;
                //Fix the parent chain
                linkToPromote->chain = parentChain;
                parentChain->linkNumber++;
                pLink = linkToPromote;
                chain->linkNumber--;
            }
            assert(pLink != NULL);
            //Attach the other end.
            pLink->nLink = nLink;
            if (nLink != NULL) {
                nLink->pLink = pLink;
            } else { //Else put pLink at the end of the chain
                parentChain->endLink = pLink;
            }

            //Cleanup old chain
            chain->link = NULL;
            assert(chain->linkNumber == 0);
            chain_destruct(chain);

            assert(chain_getLength(parentChain) == finalLength);

            return;
        }
    }
    //}

    //Calculate the final chain structure..
    stList *finalChainList = getMaximalChain(chain, flower, parentFlower);

    //Calculate the chains that are involved in the final chain list.
    stSortedSet *chainsToExpunge = stSortedSet_construct2((void(*)(void *)) chain_destruct);
    stListIterator *endIt = stList_getIterator(finalChainList);
    stIntTuple *endName;
    while ((endName = stList_getNext(endIt)) != NULL) { //Get chains in the parent which we extend..
        End *end = flower_getEnd(parentFlower, stIntTuple_get(endName, 0));
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
    End *_3End = link_get3End(chain_getFirst(chain));
    End *_5End = link_get5End(chain_getLast(chain));
    if (chain_getLength(chain) == 1 && end_isBlockEnd(_3End) && end_getOtherBlockEnd(_3End) == _5End) { //Is a circle
        assert(end_isBlockEnd(_5End));
        assert(end_getOtherBlockEnd(_5End) == _3End);
    } else {
        stList *chainEnds = stList_construct();
        if (end_isBlockEnd(_3End)) { //Is not a circle
            assert(end_getOtherBlockEnd(_3End) != _5End);
            stList_append(chainEnds, end_getOtherBlockEnd(_3End));
        }
        if (end_isBlockEnd(_5End)) { //Is not a circle
            assert(end_getOtherBlockEnd(_5End) != _3End);
            stList_append(chainEnds, end_getOtherBlockEnd(_5End));
        }
        promoteChainEnds(chainEnds, flower, parentFlower);
        stList_destruct(chainEnds);
    }

    //Redirect and reorganise all the ends, blocks and groups in the chain
    promoteEndsBlocksAndGroups(chain, flower, parentFlower);

    //Get rid of the old chains (at this point we have everything barring the final chain structure).
    stSortedSet_destruct(chainsToExpunge);

    //Make the final chain..
    Chain *newChain = chain_construct(parentFlower);
    for (int64_t i = 0; i < stList_length(finalChainList); i += 2) {
        Name _3EndName = stIntTuple_get(stList_get(finalChainList, i), 0);
        Name _5EndName = stIntTuple_get(stList_get(finalChainList, i + 1), 0);
        End *_3End = flower_getEnd(parentFlower, _3EndName);
        assert(_3End != NULL);
        End *_5End = flower_getEnd(parentFlower, _5EndName);
        assert(_5End != NULL);
        Group *group = end_getGroup(_3End);
#ifndef NDEBUG
        assert(group == end_getGroup(_5End));
        assert(group_getLink(group) == NULL);
        End *end;
        Group_EndIterator *endIt = group_getEndIterator(group);
        int64_t endNumber = 0;
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

    //We've inadvertantly created a length one chain involving just the ends of the flower
    group_constructChainForLink(parentGroup);

#ifndef NDEBUG
    assert(flower_getParentGroup(flower) == parentGroup);
    //assert(group_getLink(parentGroup) == NULL);
    if (flower_getEndNumber(flower) == 0) { //Check the properties of the flower if we've gutted it
        assert(flower_getBlockNumber(flower) == 0);
        assert(flower_getChainNumber(flower) == 0);
        assert(flower_getGroupNumber(flower) == 0);
    } else {
        assert(flower_getGroupNumber(flower)> 0);
    }
#endif
}

