#include "cactus.h"
#include "sonLib.h"
#include "stCheckEdges.h"
#include "stPerfectMatching.h"
#include "stCheckEdges.h"
#include "stMatchingAlgorithms.h"
#include "stReferenceProblem2.h"
#include <math.h>

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#endif

const char *REFERENCE_BUILDING_EXCEPTION = "REFERENCE_BUILDING_EXCEPTION";

////////////////////////////////////
////////////////////////////////////
//Get the reference event
////////////////////////////////////
////////////////////////////////////

static Event *getReferenceEvent(Flower *flower, const char *referenceEventHeader) {
    /*
     * Create the reference event, with the given header, for the flower.
     */
    EventTree *eventTree = flower_getEventTree(flower);
    Event *referenceEvent = eventTree_getEventByHeader(eventTree, referenceEventHeader);
    if (referenceEvent == NULL) {
        st_errAbort("Couldn't find reference event %s", referenceEventHeader);
    }
    return referenceEvent;
}

////////////////////////////////////
////////////////////////////////////
//Get ends from parent
////////////////////////////////////
////////////////////////////////////

static stList *getExtraAttachedStubsFromParent(Flower *flower) {
    /*
     * Copy any stubs not present in the flower from the parent group.
     * At the end there will be an even number of attached stubs in the problem.
     * Returns the list of new ends, which need to be assigned a group.
     */
    Group *parentGroup = flower_getParentGroup(flower);
    stList *endsToAdd = stList_construct();
    if (parentGroup != NULL) {
        Group_EndIterator *parentEndIt = group_getEndIterator(parentGroup);
        End *parentEnd;
        while ((parentEnd = group_getNextEnd(parentEndIt)) != NULL) {
            if (end_isAttached(parentEnd) || end_isBlockEnd(parentEnd)) {
                End *end = flower_getEnd(flower, end_getName(parentEnd));
                if (end == NULL) { //We have found an end that needs to be pushed into the child.
                    stList_append(endsToAdd, parentEnd); // end_copyConstruct(parentEnd, flower)); //At this point it has no associated group;
                }
            }
        }
        group_destructEndIterator(parentEndIt);
    }

    stList *newEnds = end_bulkCopyConstruct(endsToAdd, flower); // Now add ends in bulk
    stList_destruct(endsToAdd); // Cleanup

    assert(flower_getAttachedStubEndNumber(flower) % 2 == 0);
    if (flower_getBlockEndNumber(flower) > 0) {
        assert(flower_getAttachedStubEndNumber(flower) > 0);
    }
    return newEnds;
}

////////////////////////////////////
////////////////////////////////////
//Adjacency edge scoring
////////////////////////////////////
////////////////////////////////////

static stList *calculateZP(Cap *cap, stHash *endsToNodes) {
    /*
     * Get the list of caps that represent the ends of the chains and stubs within a sequence.
     */
    assert(!cap_getSide(cap));
    assert(end_isStubEnd(end_getPositiveOrientation(cap_getEnd(cap))));
    stList *caps = stList_construct();
    bool b = 0;
    while (1) {
        End *end = end_getPositiveOrientation(cap_getEnd(cap));
        if (stHash_search(endsToNodes, end) != NULL) {
            assert(!cap_getSide(cap));
            if (stList_length(caps) > 0) {
                assert(b);
            }
            b = 0;
            stList_append(caps, cap);
        }
        cap = cap_getAdjacency(cap);
        assert(cap != NULL);
        end = end_getPositiveOrientation(cap_getEnd(cap));
        if (stHash_search(endsToNodes, end) != NULL) {
            assert(cap_getSide(cap));
            if (stList_length(caps) > 0) {
                assert(!b);
            }
            b = 1;
            stList_append(caps, cap);
        }
        if (end_isStubEnd(end)) {
            return caps;
        }
        assert(cap != cap_getOtherSegmentCap(cap));
        cap = cap_getOtherSegmentCap(cap);
        assert(cap != NULL);
    }
    return NULL;
}

static Cap *calculateZP4(Cap *cap, stHash *endsToNodes) {
    if (cap_getOtherSegmentCap(cap) == NULL) {
        return NULL;
    }
    while (1) {
        cap = cap_getOtherSegmentCap(cap);
        assert(cap != NULL);
        End *end = end_getPositiveOrientation(cap_getEnd(cap));
        if (stHash_search(endsToNodes, end) != NULL) {
            return cap;
        }
        cap = cap_getAdjacency(cap);
        assert(cap != NULL);
        end = end_getPositiveOrientation(cap_getEnd(cap));
        assert(stHash_search(endsToNodes, end) == NULL);
        if (end_isStubEnd(end)) {
            //if(!end_isFree(end)) { //Not true is normalisation is disabled
            //    assert(!flower_hasParentGroup(end_getFlower(end)));
            //}
            return NULL;
        }
    }
    return NULL;
}

static int64_t calculateZP2(Cap *cap, stHash *endsToNodes) {
    /*
     * Calculate the length of a segment that can be traversed from a cap,
     * before hitting the end of the sequence of one of the other ends in the set
     * endsToNodes.
     */
    assert(cap_getStrand(cap));
    Sequence *sequence = cap_getSequence(cap);
    assert(sequence != NULL);
    Cap *otherCap = calculateZP4(cap, endsToNodes);
    int64_t capLength;
    if (otherCap == NULL) {
        //capLength = 1000000000; //make the length really long if attached, so that we don't bias toward one or the other end.
        capLength =
                cap_getSide(cap) ?
                        sequence_getLength(sequence) + sequence_getStart(sequence) - cap_getCoordinate(cap) :
                        cap_getCoordinate(cap) - sequence_getStart(sequence) + 1;
    } else {
        capLength =
                cap_getSide(cap) ?
                        cap_getCoordinate(otherCap) - cap_getCoordinate(cap) + 1 : cap_getCoordinate(cap) - cap_getCoordinate(otherCap) + 1;
    }
    if (capLength == 0) {
        capLength = 1;
    }
    assert(capLength > 0);
    return capLength;
}

static int64_t getBranchMultiplicitiesP(Event *pEvent, Event *event,
        stHash *branchesToMultiplicity, stSet *chosenEvents) {
    /*
     * See getBranchMultiplicities.
     */
    int64_t multiplicity = 0;
    for(int64_t i=0; i<event_getChildNumber(event); i++) {
        Event *nEvent = event_getChild(event, i);
        assert(nEvent != NULL);
        if(nEvent != pEvent) { //Don't go backwards towards the reference event.
            multiplicity += getBranchMultiplicitiesP(event, nEvent, branchesToMultiplicity, chosenEvents);
        }
    }
    if(event_getParent(event) != pEvent && event_getParent(event) != NULL) { //Case we're traversing up the tree from the reference event.
        multiplicity += getBranchMultiplicitiesP(event, event_getParent(event), branchesToMultiplicity, chosenEvents);
    }
    if(stSet_search(chosenEvents, event) != NULL) {
        multiplicity++;
    }
    if(pEvent != NULL) {
        stHash_insert(branchesToMultiplicity, event, stIntTuple_construct1(multiplicity));
    }
    return multiplicity;
}

stHash *getBranchMultiplicities(Event *referenceEvent, stSet *chosenEvents) {
    /*
     * The multiplicity of a branch is the number of simple paths that traverse it from chosen events to the given reference event.
     * Returns a map of events to multiplicities, where each event represents the branch incident with it on the path to the reference event.
     */
    stHash *branchesToMultiplicity = stHash_construct2(NULL, (void (*)(void *))stIntTuple_destruct);
    getBranchMultiplicitiesP(NULL, referenceEvent, branchesToMultiplicity, chosenEvents);
    assert(stHash_size(branchesToMultiplicity) == eventTree_getEventNumber(event_getEventTree(referenceEvent))-1); //Excludes the reference event.
    return branchesToMultiplicity;
}

static void getEventWeightingP(Event *pEvent, Event *event,
        double pathLength, double adjustedPathLength,
        stHash *branchesToMultiplicity, stHash *eventToWeights, double phi, stSet *chosenEvents) {
    /*
     * See getEventWeighting.
     */
    for(int64_t i=0; i<event_getChildNumber(event); i++) {
        Event *nEvent = event_getChild(event, i);
        if (nEvent != pEvent) {
            // Don't go backwards towards the reference event
            assert(stHash_search(branchesToMultiplicity, nEvent) != NULL);
            int64_t multiplicity = stIntTuple_get(stHash_search(branchesToMultiplicity, nEvent), 0);

            if(multiplicity > 0) { // Don't traverse paths not leading to interesting events
                getEventWeightingP(event, nEvent,
                                   pathLength + event_getBranchLength(nEvent),
                                   adjustedPathLength + event_getBranchLength(nEvent)/multiplicity,
                                   branchesToMultiplicity, eventToWeights, phi,
                                   chosenEvents);
            }
        }
    }
    Event *nEvent;
    if((nEvent = event_getParent(event)) != pEvent && nEvent != NULL) { //Case we're traversing up the tree from the reference event to a node with another child lineage.
        assert(stHash_search(branchesToMultiplicity, nEvent) != NULL);
        int64_t multiplicity = stIntTuple_get(stHash_search(branchesToMultiplicity, nEvent), 0);
        if(multiplicity > 0) {
            getEventWeightingP(event, nEvent,
                    pathLength + event_getBranchLength(event),
                    adjustedPathLength + event_getBranchLength(event)/multiplicity,
                    branchesToMultiplicity, eventToWeights, phi, chosenEvents);
        }
    }
    if(stSet_search(chosenEvents, event) != NULL) { //Defines a leaf event that is not the reference, which we need to give a score for
        assert(adjustedPathLength <= pathLength);
        double score = exp(-phi * pathLength);
        assert(score <= 1.0);
        assert(score >= 0.0);
        if (pathLength > 0.0) {
            score *= adjustedPathLength/pathLength;
            assert(score <= 1.0);
            assert(score >= 0.0);
        }
        // fprintf(stdout, "Chose weight %lf for event %s (multiplicity %" PRIi64 "). Adj path length %lf, path length %lf.\n", score, event_getHeader(event), stIntTuple_get(stHash_search(branchesToMultiplicity, event), 0), adjustedPathLength, pathLength);
        stHash_insert(eventToWeights, event, stDoubleTuple_construct(1, score));
    }
}

stHash *getEventWeighting(Event *referenceEvent, double phi, stSet *chosenEvents) {
    /*
     * Weights events by how informative they are for inferring the reference event.
     * Accounts for both distance and the sharing of branches. Returns a hash of chosen events to weights.
     *
     * Let R be the chosen reference event, S the set of chosen event events and A a member of S. Let b_1, b_2, ..., b_n be the branches on
     * the simple path from R to A. Let d(b_i) be the length of the branch b_i and s(b_i) the number of simple paths from R to a member of S
     * that pass through b_i, termed multiplicity.
     *
     * The independence weight a_{R,A} is \sum_{i} d(b_i)/s(b_i) / \sum_{i} d(b_i).
     *
     * We calculate the total weight for A as e^(-phi * \sum_{i} d(b_i)) * a_{R,A}
     */
    stHash *eventToWeightHash = stHash_construct2(NULL, (void (*)(void *))stDoubleTuple_destruct);

    //Calculate the multiplicity (s function) of branches.
    stHash *branchesToMultiplicity = getBranchMultiplicities(referenceEvent, chosenEvents);

    //Calculate total weights
    getEventWeightingP(NULL, referenceEvent, 0.0, 0.0, branchesToMultiplicity, eventToWeightHash, phi, chosenEvents);
    stHash_destruct(branchesToMultiplicity);
    assert(stSet_size(chosenEvents) == stHash_size(eventToWeightHash));

    return eventToWeightHash;
}

static stSet *getEventsWithSequences(Flower *flower) {
    /*
     * Returns all the events in the event tree that have sequences associated with them.
     */
    stSet *seqSet = stSet_construct();
    Flower_SequenceIterator *seqIt = flower_getSequenceIterator(flower);
    Sequence *seq;
    while((seq = flower_getNextSequence(seqIt)) != NULL) {
        stSet_insert(seqSet, sequence_getEvent(seq));
    }
    flower_destructSequenceIterator(seqIt);
    return seqSet;
}

static double calculateZScoreWeightedAdapterFn(Cap *_5Cap, int64_t length5Segment, int64_t length3Segment, int64_t gap, void *extraArgs) {
    double theta = *((double *)((void **) extraArgs)[0]);
    assert(theta >= 0.0);
    assert(cap_getEvent(_5Cap) != NULL);
    //fprintf(stderr, " hello event: %" PRIi64 "\n", cap_getEvent(_5Cap));
    //fprintf(stderr, " hello seq: %" PRIi64 "\n", cap_getSequence(_5Cap));
    //fprintf(stderr, " hello event name: %" PRIi64 " \n", event_getName(cap_getEvent(_5Cap)));
    assert(stHash_search(((void **) extraArgs)[1], cap_getEvent(_5Cap)) != NULL);
    assert(stDoubleTuple_length(stHash_search(((void **) extraArgs)[1], cap_getEvent(_5Cap))) == 1);
    double weight = stDoubleTuple_getPosition(stHash_search(((void **) extraArgs)[1], cap_getEvent(_5Cap)), 0);
    return calculateZScore(length5Segment, length3Segment, gap, theta) * weight;
}

static double countAdapterFn(Cap *_5Cap, int64_t length5Segment, int64_t length3Segment, int64_t gap, void *extraArgs) {
    return 1;
}

refAdjList *calculateZ(Flower *flower, stHash *endsToNodes, int64_t nodeNumber, int64_t maxWalkForCalculatingZ,
bool ignoreUnalignedGaps, double (*zScoreFn)(Cap *, int64_t, int64_t, int64_t, void *), void *zScoreExtraArgs) {
    /*
     * Calculate the zScores between all ends.
     */
    refAdjList *aL = refAdjList_construct(nodeNumber);
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end)) {
            End_InstanceIterator *capIt = end_getInstanceIterator(end);
            Cap *cap;
            while ((cap = end_getNext(capIt)) != NULL) {
                cap = cap_getStrand(cap) ? cap : cap_getReverse(cap);
                if (!cap_getSide(cap) && cap_getSequence(cap) != NULL) {
                    stList *caps = calculateZP(cap, endsToNodes);

                    /*
                     * Calculate the lengths of the sequences following the 3 caps, for efficiency.
                     */
                    int64_t *capSizes = st_malloc(sizeof(int64_t) * stList_length(caps));
                    for (int64_t i = 0; i < stList_length(caps); i++) {
                        Cap *cap = stList_get(caps, i);
                        capSizes[i] = calculateZP2(cap, endsToNodes);
                    }

                    /*
                     * Iterate through all pairs of 5' and 3' caps to calculate additions to scores.
                     */
                    for (int64_t i = (stList_length(caps) > 0 && cap_getSide(stList_get(caps, 0))) ? 1 : 0; i < stList_length(caps); i += 2) {
                        Cap *_3Cap = stList_get(caps, i);
                        assert(!cap_getSide(_3Cap));
                        int64_t _3CapSize = capSizes[i];
                        int64_t _3Node = stIntTuple_get(stHash_search(endsToNodes, end_getPositiveOrientation(cap_getEnd(_3Cap))), 0);
                        int64_t unaligned = 0;
                        for (int64_t k = 0; k < maxWalkForCalculatingZ; k++) {
                            int64_t j = k * 2 + i + 1;
                            if (j >= stList_length(caps)) {
                                break;
                            }
                            Cap *_5Cap = stList_get(caps, j);
                            assert(cap_getSide(_5Cap));
                            assert(cap_getAdjacency(_5Cap) != NULL);
                            if (ignoreUnalignedGaps) {
                                assert(cap_getCoordinate(_5Cap) - cap_getCoordinate(cap_getAdjacency(_5Cap)) - 1 >= 0);
                                unaligned += cap_getCoordinate(_5Cap) - cap_getCoordinate(cap_getAdjacency(_5Cap)) - 1;
                            }
                            int64_t _5Node = stIntTuple_get(stHash_search(endsToNodes, end_getPositiveOrientation(cap_getEnd(_5Cap))), 0);
                            int64_t _5CapSize = capSizes[j];
                            assert(cap_getCoordinate(_5Cap) - cap_getCoordinate(_3Cap) > 0);
                            int64_t diff = cap_getCoordinate(_5Cap) - cap_getCoordinate(_3Cap) - unaligned;
                            assert(diff >= 1);
                            if (zScoreFn(_5Cap, 1, 1, diff, zScoreExtraArgs) < 0.0000000001) { //no point walking when score gets too small, should be effective for theta >= 0.000001
                                break;
                            }
                            double score = zScoreFn(_5Cap, _5CapSize, _3CapSize, diff, zScoreExtraArgs);
                            assert(score >= -0.0001);
                            if (score <= 0.0) {
                                score = 1e-10; //Make slightly non-zero.
                            }
                            assert(score > 0.0);
                            refAdjList_addToWeight(aL, _3Node, _5Node, score);
                            assert(refAdjList_getWeight(aL, _3Node, _5Node) == refAdjList_getWeight(aL, _5Node, _3Node));
                            assert(refAdjList_getWeight(aL, _3Node, _5Node) >= 0.0);
                        }
                    }
                    stList_destruct(caps);
                    free(capSizes);
                }
            }
            end_destructInstanceIterator(capIt);
        }
    }
    flower_destructEndIterator(endIt);

    return aL;
}

////////////////////////////////////
////////////////////////////////////
//Chain edges
////////////////////////////////////
////////////////////////////////////

static End *getEndFromNode(stHash *nodesToEnds, int64_t node) {
    /*
     * Get the end for the given node.
     */
    stIntTuple *i = stIntTuple_construct1(node);
    End *end = stHash_search(nodesToEnds, i);
    assert(end != NULL);
    assert(end_getOrientation(end));
    stIntTuple_destruct(i);
    return end;
}

static void addToNodes(End *end, int64_t n, stHash *endsToNodes) {
    stHash_insert(endsToNodes, end_getPositiveOrientation(end), stIntTuple_construct1(n));
}

static void getNonTrivialChainNodes(Flower *flower, stHash *endsToNodes, int64_t *nodeCounter) {
    /*
     * Get the non-trivial chain edges for the flower.
     */
    Flower_ChainIterator *chainIt = flower_getChainIterator(flower);
    Chain *chain;
    while ((chain = flower_getNextChain(chainIt)) != NULL) {
        End *_5End = link_get3End(chain_getFirst(chain));
        End *_3End = link_get5End(chain_getLast(chain));
        if (end_isBlockEnd(_5End) && end_isBlockEnd(_3End)) {
            End *end1 = end_getOtherBlockEnd(_5End);
            End *end2 = end_getOtherBlockEnd(_3End);

            assert(end_getGroup(end1) != NULL);
            assert(end_getGroup(end2) != NULL);
            assert(group_isTangle(end_getGroup(end1)));
            assert(group_isTangle(end_getGroup(end2)));

            addToNodes(end1, *nodeCounter, endsToNodes);
            addToNodes(end2, -(*nodeCounter), endsToNodes);
            (*nodeCounter)++;
        }
    }
    flower_destructChainIterator(chainIt);
}

static void getTrivialChainNodes(Flower *flower, stHash *endsToNodes, int64_t *nodeCounter) {
    /*
     * Get the trivial chain edges for the flower.
     */
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if(end_partOfBlock(end) && end_left(end)) {
            Block *block = end_getBlock(end);
            End *_5End = block_get5End(block);
            End *_3End = block_get3End(block);
            assert(end_getGroup(_5End) != NULL);
            assert(end_getGroup(_3End) != NULL);
            if (group_isTangle(end_getGroup(_5End)) && group_isTangle(end_getGroup(_3End))) {
                addToNodes(_5End, *nodeCounter, endsToNodes);
                addToNodes(_3End, -(*nodeCounter), endsToNodes);
                (*nodeCounter)++;
            }
        }
    }
    flower_destructEndIterator(endIt);
}

static stHash *getChainNodes(Flower *flower) {
    /*
     * Get the trivial chain edges for the flower.
     */
    int64_t nodeCounter = 1;
    stHash *endsToNodes = stHash_construct2(NULL, NULL); //We clean up the tuple memory when we invert the hash to go nodesToEnds
    getNonTrivialChainNodes(flower, endsToNodes, &nodeCounter);
    getTrivialChainNodes(flower, endsToNodes, &nodeCounter);
    return endsToNodes;
}

////////////////////////////////////
////////////////////////////////////
//Stub edges
////////////////////////////////////
////////////////////////////////////

static stList *getTangleStubEnds(Flower *flower, stHash *endsToNodes) {
    /*
     * Iterates over attached stub ends not in chain and adds them to nodes list (from the parent, without a group, and those ends in tangles).
     */
    assert(stHash_size(endsToNodes) % 2 == 0);
    int64_t endCount = stHash_size(endsToNodes) / 2 + 1;
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    stList *stubEnds = stList_construct();
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        Group *group = end_getGroup(end);
        if (group != NULL) {
            if (group_isTangle(group)) {
                if ((end_isAttached(end) || end_isBlockEnd(end)) && stHash_search(endsToNodes, end) == NULL) {
                    stHash_insert(endsToNodes, end, stIntTuple_construct1(endCount++));
                    stList_append(stubEnds, end);
                }
            } else {
                assert(group_isLink(group));
            }
        } else {
            assert(end_isStubEnd(end));
            assert(end_isAttached(end));
            assert(stHash_search(endsToNodes, end) == NULL);
            stHash_insert(endsToNodes, end, stIntTuple_construct1(endCount++));
            stList_append(stubEnds, end);
        }
    }
    flower_destructEndIterator(endIt);
    return stubEnds;
}

static Cap *getCapWithEvent(End *end, Name eventName) {
    /*
     * Get the first encountered cap in the end with the given event.
     */
    Cap *cap;
    End_InstanceIterator *instanceIt = end_getInstanceIterator(end);
    while ((cap = end_getNext(instanceIt)) != NULL) {
        if (event_getName(cap_getEvent(cap)) == eventName) {
            end_destructInstanceIterator(instanceIt);
            return cap;
        }
    }
    end_destructInstanceIterator(instanceIt);
    return NULL;
}

static End *getAdjacentEndFromParent(End *end, Event *referenceEvent) {
    /*
     * Get the adjacent stub end by looking at the reference adjacency in the parent.
     */
    Flower *flower = end_getFlower(end);
    Group *parentGroup = flower_getParentGroup(flower);
    assert(parentGroup != NULL);
    End *parentEnd = group_getEnd(parentGroup, end_getName(end));
    assert(parentEnd != NULL);
    //Now trace the parent end..
    Cap *cap = getCapWithEvent(parentEnd, event_getName(referenceEvent));
    assert(cap != NULL);
    Cap *adjacentCap = cap_getAdjacency(cap);
    assert(adjacentCap != NULL);
    End *adjacentParentEnd = cap_getEnd(adjacentCap);
    End *adjacentEnd = flower_getEnd(flower, end_getName(adjacentParentEnd));
    assert(adjacentEnd != NULL);
    return adjacentEnd;
}

static End *getStubEdgesFromParent2(End *end) {
    assert(end_getOrientation(end));
    assert(end_getGroup(end) != NULL);
    Link *link = group_getLink(end_getGroup(end));
    assert(link != NULL);
    Chain *chain = link_getChain(link);
    assert(chain_getFirst(chain) == link || chain_getLast(chain) == link);
    assert(link_get3End(link) == end || link_get5End(link) == end);
    if (chain_getFirst(chain) == chain_getLast(chain)) {
        return link_get3End(link) == end ? link_get5End(link) : link_get3End(link);
    }
    if (chain_getFirst(chain) == link) {
        assert(link_get3End(link) == end);
        return link_get5End(chain_getLast(chain));
    } else {
        assert(link_get5End(link) == end);
        return link_get3End(chain_getFirst(chain));
    }
}

static void getStubEdgesFromParent(refOrdering *ref, Flower *flower, Event *referenceEvent, stHash *endsToNodes, stList *stubEnds) {
    /*
     * For each attached stub in the flower, get the end in the parent group, and add its
     * adjacency to the group.
     */
    stSortedSet *endsSeen = stSortedSet_construct();
    for (int64_t i = 0; i < stList_length(stubEnds); i++) {
        End *end = stList_get(stubEnds, i);
        assert(end_getOrientation(end));
        if (stSortedSet_search(endsSeen, end) == NULL) {
            End *end2 = end_isBlockEnd(end) ? getStubEdgesFromParent2(end_getOtherBlockEnd(end)) : end;
            assert(end_isStubEnd(end2));
            End *adjacentEnd = getAdjacentEndFromParent(end2, referenceEvent);
            if (stHash_search(endsToNodes, adjacentEnd) == NULL) {
                adjacentEnd = getStubEdgesFromParent2(adjacentEnd);
                assert(end_isBlockEnd(adjacentEnd));
                adjacentEnd = end_getOtherBlockEnd(adjacentEnd);
                assert(stHash_search(endsToNodes, adjacentEnd) != NULL);
            }
            assert(end_getOrientation(adjacentEnd));
            stSortedSet_insert(endsSeen, end);
            stSortedSet_insert(endsSeen, adjacentEnd);
            reference_makeNewInterval(ref, -stIntTuple_get(stHash_search(endsToNodes, end), 0),
                    stIntTuple_get(stHash_search(endsToNodes, adjacentEnd), 0));
        }
    }
    stSortedSet_destruct(endsSeen);
}

stHash *makeStubEdgesToNodesHash(stList *stubEnds, stHash *endsToNodes) {
    stHash *stubEndsToNodes = stHash_construct();
    for (int64_t i = 0; i < stList_length(stubEnds); i++) {
        End *end = stList_get(stubEnds, i);
        stHash_insert(stubEndsToNodes, end, stHash_search(endsToNodes, end));
    }
    return stubEndsToNodes;
}

static void getStubEdgesInTopLevelFlower(refOrdering *ref, Flower *flower, stHash *endsToNodes, int64_t nodeNumber, Event *referenceEvent,
        stList *(*matchingAlgorithm)(stList *edges, int64_t nodeNumber), stList *stubEnds, double phi) {
    /*
     * Create a matching for the parent stub edges.
     */
    stHash *stubEndsToNodes = makeStubEdgesToNodesHash(stubEnds, endsToNodes);
    double theta = 0.0;
    stSet *chosenEvents = getEventsWithSequences(flower);
    stHash *eventWeighting = getEventWeighting(referenceEvent, phi, chosenEvents);
    stSet_destruct(chosenEvents);
    void *zArgs[2] = { &theta, eventWeighting };
    refAdjList *stubAL = calculateZ(flower, stubEndsToNodes, nodeNumber,
    INT64_MAX, 1, calculateZScoreWeightedAdapterFn, zArgs);
    stHash_destruct(eventWeighting);
    st_logInfo(
            "Building a matching for %" PRIi64 " stub nodes in the top level problem from %" PRIi64 " total stubs of which %"
            PRIi64 " attached , %" PRIi64 " total ends, %" PRIi64 " chains, %" PRIi64 " blocks %" PRIi64 " groups and %" PRIi64 " sequences\n",
            stList_length(stubEnds), flower_getStubEndNumber(flower), flower_getAttachedStubEndNumber(flower), flower_getEndNumber(flower),
            flower_getChainNumber(flower), flower_getBlockNumber(flower), flower_getGroupNumber(flower), flower_getSequenceNumber(flower));

    stList *adjacencyEdges = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < stList_length(stubEnds); i++) {
        int64_t node1 = stIntTuple_get(stHash_search(endsToNodes, stList_get(stubEnds, i)), 0);
        for (int64_t j = i + 1; j < stList_length(stubEnds); j++) {
            int64_t node2 = stIntTuple_get(stHash_search(endsToNodes, stList_get(stubEnds, j)), 0);
            double score = refAdjList_getWeight(stubAL, node1, node2);
            assert(score >= 0);
            int64_t score2 = score > INT64_MAX ? INT64_MAX : score;
            assert(score2 >= 0);
            stList_append(adjacencyEdges, constructWeightedEdge(node1, node2, score2));
        }
    }
    stList *stubNodes = stHash_getValues(stubEndsToNodes);
    stSortedSet *stubNodesSet = stList_getSortedSet(stubNodes, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
    stList_destruct(stubNodes);
    checkEdges(adjacencyEdges, stubNodesSet, 1, 0);

    stList *chosenAdjacencyEdges = getPerfectMatching(stubNodesSet, adjacencyEdges, matchingAlgorithm);

    for (int64_t i = 0; i < stList_length(chosenAdjacencyEdges); i++) {
        stIntTuple *adjacencyEdge = stList_get(chosenAdjacencyEdges, i);
        assert(stIntTuple_length(adjacencyEdge) == 3);
        reference_makeNewInterval(ref, -stIntTuple_get(adjacencyEdge, 0), stIntTuple_get(adjacencyEdge, 1));
    }

    stHash_destruct(stubEndsToNodes);
    stList_destruct(chosenAdjacencyEdges);
    stList_destruct(adjacencyEdges);
    stSortedSet_destruct(stubNodesSet);
    refAdjList_destruct(stubAL);
}

static refOrdering *getEmptyReference(Flower *flower, stHash *endsToNodes, int64_t nodeNumber, Event *referenceEvent,
        stList *(*matchingAlgorithm)(stList *edges, int64_t nodeNumber), stList *stubEnds, double phi) {
    refOrdering *ref = reference_construct(nodeNumber);
    if (flower_getParentGroup(flower) != NULL) {
        getStubEdgesFromParent(ref, flower, referenceEvent, endsToNodes, stubEnds);
    } else {
        getStubEdgesInTopLevelFlower(ref, flower, endsToNodes, nodeNumber, referenceEvent, matchingAlgorithm, stubEnds, phi);
    }
    return ref;
}

////////////////////////////////////
////////////////////////////////////
//Functions to add adjacencies and segments, given the chosen edges
////////////////////////////////////
////////////////////////////////////

static Cap *makeStubCap(End *end, Event *referenceEvent) {
    //assert(end_isAttached(end));
    assert(end_isStubEnd(end));
    Cap *cap = getCapWithEvent(end, event_getName(referenceEvent));
    if (cap != NULL) {
        return cap_getStrand(cap) ? cap : cap_getReverse(cap);
    }
    Group *parentGroup = flower_getParentGroup(end_getFlower(end));
    End *parentEnd;
    if (parentGroup != NULL && (parentEnd = group_getEnd(parentGroup, end_getName(end))) != NULL) {
        Cap *parentCap = getCapWithEvent(parentEnd, event_getName(referenceEvent));
        assert(parentCap != NULL);
        parentCap = cap_getStrand(parentCap) ? parentCap : cap_getReverse(parentCap);
        if (end_getSide(end) != cap_getSide(parentCap)) {
            end = end_getReverse(end);
        }
        cap = cap_copyConstruct(end, parentCap);
        assert(cap_getSide(cap) == cap_getSide(parentCap));
        assert(cap_getStrand(cap));
        assert(cap_getSide(cap) == end_getSide(end));
    } else {
        cap = cap_construct(end, referenceEvent);
        cap_setCoordinates(cap, INT64_MAX, 1, NULL);
        assert(end_getSide(end) == cap_getSide(cap));
        assert(cap_getStrand(cap));
    }
    return cap;
}

static void makeAdjacent(Cap *cap, Cap *cap2) {
    assert(cap_getStrand(cap));
    assert(cap_getStrand(cap2));
    assert(cap_getSide(cap) != cap_getSide(cap2));
    cap_makeAdjacent(cap, cap2);
}

static void makeThread(End *end, stHash *endsToEnds, Flower *flower, Event *referenceEvent) {
    Cap *cap = makeStubCap(end, referenceEvent);
    assert(cap_getStrand(cap));
    end = cap_getEnd(cap);
    assert(cap_getSide(cap) == end_getSide(end));
    while (1) {
        End *end2 = stHash_search(endsToEnds, end_getPositiveOrientation(end));
        assert(end2 != NULL);
        if (end_getSide(end) == end_getSide(end2)) {
            end2 = end_getReverse(end2);
        }
        assert(end_getSide(end) != end_getSide(end2));
        if (end_isStubEnd(end2)) {
            Cap *cap2 = makeStubCap(end2, referenceEvent);
            makeAdjacent(cap, cap2);
            break;
        }
        if (getCapWithEvent(end2, event_getName(referenceEvent)) != NULL) { //Already connected, case where thread is already generated.
            assert(cap_getAdjacency(getCapWithEvent(end2, event_getName(referenceEvent))) == cap);
            break;
        }
        Block *block = end_getBlock(end2);
        Segment *segment = segment_construct(block, referenceEvent);
        Cap *cap2 = segment_get5Cap(segment);
        assert(cap_getSide(cap2));
        Cap *cap3 = segment_get3Cap(segment);
        assert(!cap_getSide(cap3));
        cap_setCoordinates(cap2, INT64_MAX, 1, NULL);
        cap_setCoordinates(cap3, INT64_MAX, 1, NULL);
        if (end_getSide(end2)) {
            assert(block_get5End(block) == end2);
            makeAdjacent(cap, cap2);
            cap = cap3;
        } else {
            assert(block_get3End(block) == end2);
            makeAdjacent(cap, cap3);
            cap = cap2;
        }
        end = cap_getEnd(cap);
    }
}

static void makeThreads(stHash *endsToEnds, Flower *flower, Event *referenceEvent) {
    Flower_EndIterator *endIt = flower_getEndIterator(flower);
    End *end;
    stList *stack = stList_construct();
    //Group *parentGroup = flower_getParentGroup(flower);
    while ((end = flower_getNextEnd(endIt)) != NULL) {
        if (end_isStubEnd(end)) {
            if (end_isAttached(end)) {
                makeThread(end, endsToEnds, flower, referenceEvent); //We first start from stub ends that already have
                //ends in parent flower, because they define the 5' to 3' direction of traversal.
            } else if (stHash_search(endsToEnds, end) != NULL) { //In the reference
                stList_append(stack, end);
            }
        }
    }
    flower_destructEndIterator(endIt);
    while (stList_length(stack) > 0) {
        makeThread(stList_pop(stack), endsToEnds, flower, referenceEvent);
    }
    stList_destruct(stack);
}

static void mapEnds(stHash *endsToEnds, End *end1, End *end2) {
    assert(end1 != end2);
    assert(end_getGroup(end1) == end_getGroup(end2) || end_getGroup(end1) == NULL || end_getGroup(end2) == NULL);
    assert(end_getOrientation(end1));
    assert(end_getOrientation(end2));
    assert(stHash_search(endsToEnds, end1) == NULL);
    assert(stHash_search(endsToEnds, end2) == NULL);
    stHash_insert(endsToEnds, end1, end2);
    stHash_insert(endsToEnds, end2, end1);
}

static stHash *getEndsToEnds(Flower *flower, stList *chosenAdjacencyEdges, stHash *nodesToEnds, int64_t numberOfNsForScaffoldGap) {
    /*
     * Get a hash of matched ends.
     */
    stHash *endsToEnds = stHash_construct();
    for (int64_t i = 0; i < stList_length(chosenAdjacencyEdges); i++) {
        stIntTuple *edge = stList_get(chosenAdjacencyEdges, i);
        End *end1 = getEndFromNode(nodesToEnds, stIntTuple_get(edge, 0));
        End *end2 = getEndFromNode(nodesToEnds, stIntTuple_get(edge, 1));
        assert(end1 != NULL);
        assert(end2 != NULL);
        assert(end1 != end2);

        if (end_getGroup(end1) != NULL && end_getGroup(end2) != NULL && end_getGroup(end1) != end_getGroup(end2)) {
            if(group_getNestedFlower(end_getGroup(end1)) == NULL && !group_isLink(end_getGroup(end1)) &&
               group_getNestedFlower(end_getGroup(end2)) == NULL && !group_isLink(end_getGroup(end2))) {
                /*
                 * Merge the groups together so we can join them together, providing they are part of terminal groups
                 * and not part of chains
                 */
                group_mergeTerminalGroups(end_getGroup(end1), end_getGroup(end2));
                assert(end_getGroup(end1) == end_getGroup(end2));
                mapEnds(endsToEnds, end1, end2);
            }
            else {
                /*
                 * We build a 'scaffolding block' between ends
                 *
                 * We require the ends that are being scaffolded to already have groups, else they must themselves the
                 * ends of scaffolding blocks, so we don't add additional scaffold gaps.
                 */
                Block *block = block_construct(numberOfNsForScaffoldGap, flower);
                end_setGroup(block_get5End(block), end_getGroup(end1));
                assert(group_isTangle(end_getGroup(end1)));
                mapEnds(endsToEnds, end1, block_get5End(block));
                end_setGroup(block_get3End(block), end_getGroup(end2));
                assert(group_isTangle(end_getGroup(end2)));
                mapEnds(endsToEnds, end2, block_get3End(block));
            }
        } else {
            mapEnds(endsToEnds, end1, end2);
        }
    }
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIt)) != NULL) {
        if (group_isLink(group)) {
            Link *link = group_getLink(group);
            mapEnds(endsToEnds, link_get3End(link), link_get5End(link));
        }
    }
    flower_destructGroupIterator(groupIt);
    assert(flower_getEndNumber(flower) - flower_getFreeStubEndNumber(flower) <= stHash_size(endsToEnds));
    return endsToEnds;
}

static void makeReferenceThreads(Flower *flower, stList *chosenAdjacencyEdges, stHash *nodesToEnds, Event *referenceEvent,
        int64_t numberOfNsForScaffoldGap) {
    stHash *endsToEnds = getEndsToEnds(flower, chosenAdjacencyEdges, nodesToEnds, numberOfNsForScaffoldGap);
    makeThreads(endsToEnds, flower, referenceEvent);
    stHash_destruct(endsToEnds);
}

////////////////////////////////////
////////////////////////////////////
//Put the new ends into groups.
////////////////////////////////////
////////////////////////////////////

static void assignGroups(stList *newEnds, Flower *flower, Event *referenceEvent) {
    /*
     * Put ends into groups.
     */
    for (int64_t i = 0; i < stList_length(newEnds); i++) {
        End *end = stList_get(newEnds, i);
        if (end_getGroup(end) == NULL) { //It may already have been assigned.
            Cap *cap = getCapWithEvent(end, event_getName(referenceEvent));
            assert(cap != NULL);
            Cap *adjacentCap = cap_getAdjacency(cap);
            assert(adjacentCap != NULL);
            End *adjacentEnd = cap_getEnd(adjacentCap);
            Group *group = end_getGroup(adjacentEnd);
            if (group == NULL) {
                /*
                 * Two new ends link to one another, hence we need to find an existing tangle group.
                 */
                Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
                while ((group = flower_getNextGroup(groupIt)) != NULL) {
                    if (group_isTangle(group)) {
                        break;
                    }
                }
                assert(group_isTangle(group));
                flower_destructGroupIterator(groupIt);
                end_setGroup(adjacentEnd, group);
            }
            end_setGroup(end, group);
        }
    }
}

static stList *convertReferenceToAdjacencyEdges2(refOrdering *ref) {
    stList *edges = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t i = 0; i < reference_getIntervalNumber(ref); i++) {
        int64_t n = -reference_getFirstOfInterval(ref, i);
        while (reference_getNext(ref, n) != INT64_MAX) {
            stList_append(edges, constructEdge(n, reference_getNext(ref, n)));
            n = -reference_getNext(ref, n);
        }
    }
    return edges;
}

/*
 * Function to add additional ends to graph, representing breaks in the reference.
 */
static void addAdditionalStubEnds(stList *extraStubNodes, Flower *flower, stHash *nodesToEnds, stList *newEnds) {
    for (int64_t i = 0; i < stList_length(extraStubNodes); i++) {
        stIntTuple *node = stList_get(extraStubNodes, i);
        End *end = end_construct(0, flower); //Ensure we get the right orientation
        stHash_insert(nodesToEnds, stIntTuple_construct1(stIntTuple_get(node, 0)), end);
        stList_append(newEnds, end);
    }
}


bool referenceSplitFn(int64_t pNode, refOrdering *ref, void *extraArgs) {
    assert(reference_getNext(ref, pNode) != INT64_MAX);

    // Unpack extraArgs
    stHash *nodesToEnds = ((void **) extraArgs)[0];
    refAdjList *dAL = ((void **) extraArgs)[1];
    int64_t minNumberOfSequencesToSupportAdjacency = *((int64_t *) ((void **) extraArgs)[2]);
    bool canBreakAdjacency = *((bool *) ((void **) extraArgs)[3]);

    if(!canBreakAdjacency) {  // If we've determined we can not break the adjacency because of surrounding structure
        return false;
    }

    // First, check for sufficient direct adjacency support. If it exists, we never split.
    if (refAdjList_getWeight(dAL, -pNode, reference_getNext(ref, pNode)) >= minNumberOfSequencesToSupportAdjacency) {
        return false;
    }

    // An adjacency is only broken if it is between ends in a leaf group.
    stIntTuple *i = stIntTuple_construct1(-pNode);
    End *end = stHash_search(nodesToEnds, i);
    stIntTuple_destruct(i);
    assert(end != NULL);
    Group *group = end_getGroup(end);
    if (group == NULL) {
        i = stIntTuple_construct1(reference_getNext(ref, pNode));
        End *adjacentEnd = stHash_search(nodesToEnds, i);
        stIntTuple_destruct(i);
        assert(adjacentEnd != NULL);
        group = end_getGroup(adjacentEnd);
    }
    if (group == NULL || !group_isLeaf(group)) {
        // This can happen for new ends from scaffolding, or if it's not a tangle. Don't split.
        return false;
    }
    return true;
}

stList *getReferenceIntervalsToPreserve(refOrdering *ref, refAdjList *dAL, int64_t minNumberOfSequencesToSupportAdjacency) {
    stList *referenceIntervalsToPreserve = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    for (int64_t interval = 0; interval < reference_getIntervalNumber(ref); interval++) {
        int64_t firstNode = reference_getFirstOfInterval(ref, interval);
        int64_t lastNode = reference_getLast(ref, firstNode);
        assert(reference_getNext(ref, firstNode) == lastNode);
        if (refAdjList_getWeight(dAL, -firstNode, lastNode) >= minNumberOfSequencesToSupportAdjacency) { //Decide if we want to preserve the interval
            stList_append(referenceIntervalsToPreserve, stIntTuple_construct2(firstNode, lastNode));
        }
    }
    return referenceIntervalsToPreserve;
}

bool canBreakAdjacencies(Flower *flower,
                         int64_t minimumNestedBasesToBreakAdjacency,
                         int64_t maximumChainBasesToBreakAdjacency) {
    /*
     * The following is the logic for deciding if we can split an adjacency in a given flower
     */

    // Look in parents to decide if nested in large chain
    Group *ancestor_group;
    while((ancestor_group = flower_getParentGroup(flower)) != NULL) {
        // If flower contains more than x bases, regardless of if in large chain, determine we can break
        if(flower_getTotalBaseLength(flower) > minimumNestedBasesToBreakAdjacency) {
            return true;
        }

        // If parent group is part of chain, then determine if chain is large enough to
        // decide not to break adajcencies
        if(group_isLink(ancestor_group)) {
            Chain *parent_chain = link_getChain(ancestor_group);
            // Count the bases in the chain
            Link *l = chain_getFirst(parent_chain);
            int64_t i = 0;
            while(l != NULL && i < maximumChainBasesToBreakAdjacency) {
                End *end = link_get5End(l);
                if(end_partOfBlock(end)) {
                    i += block_getLength(end_getBlock(end));
                }
                l = link_getNextLink(l);
            }
            if(i >= maximumChainBasesToBreakAdjacency) {
                return false;
            }
        }
        flower = group_getFlower(ancestor_group);
    }
    return true;
}

////////////////////////////////////
////////////////////////////////////
//Main function
////////////////////////////////////
////////////////////////////////////

void buildReferenceTopDown(Flower *flower, const char *referenceEventHeader, int64_t permutations,
        stList *(*matchingAlgorithm)(stList *edges, int64_t nodeNumber), double (*temperature)(double),
        double theta, double phi, int64_t maxWalkForCalculatingZ,
        bool ignoreUnalignedGaps, double wiggle, int64_t numberOfNsForScaffoldGap, int64_t minNumberOfSequencesToSupportAdjacency, bool makeScaffolds,
        int64_t minimumNestedBasesToBreakAdjacency, int64_t maximumChainBasesToBreakAdjacency) {
    /**
     * @brief Constructs the reference sequence paths for a given flower.
     *
     * This is the main function for building the reference. It orchestrates a multi-step process
     * to determine the optimal ordering and connections of sequence segments (blocks) within a flower,
     * effectively creating a reference genome layout for that part of the pangenome graph.
     *
     * The process involves:
     * 1.  Identifying all relevant endpoints (chain ends and unchained "stub" ends) which are
     *     the nodes in the adjacency problem.
     * 2.  Calculating weighted adjacency scores (Z-scores) between all pairs of nodes. These scores
     *     are informed by phylogenetic distance (phi) and sequence continuity (theta).
     * 3.  For the top-level flower, it solves an initial matching problem for stub ends that connect
     *     to the parent, or uses a matching algorithm if there is no parent.
     * 4.  Building an initial linear ordering of nodes using a greedy algorithm (`makeReferenceGreedily2`).
     * 5.  Refining this ordering through iterative sampling (`updateReferenceGreedily`) and local
     *     rearrangements (`nudgeGreedily`) to improve the overall score.
     * 6.  Breaking the reference at adjacencies that lack sufficient sequence support. An adjacency is broken if:
     *     (1) it has fewer than `minNumberOfSequencesToSupportAdjacency` direct adjacencies, AND
     *     (2) the flower contains at least `minimumNestedBasesToBreakAdjacency` bases, AND
     *     (3) the flower is either the top-level flower OR it is a nested flower ("snarl") within a
     *         chain that is no longer than `maximumChainBasesToBreakAdjacency` links long.
     * 7.  Optionally, creating new "scaffold" blocks to bridge gaps between ends that are inferred
     *     to be adjacent but are not connected by sequence.
     * 8.  Finally, creating the actual reference threads (sequences and adjacencies) in the flower,
     *     modifying it in-place.
     *
     * @param flower The flower for which to build the reference. It will be modified in-place.
     * @param referenceEventHeader The header of the event to be treated as the reference, defining the
     *        coordinate space and phylogenetic root for weighting.
     * @param permutations The number of greedy permutation rounds to perform for refining the reference ordering.
     * @param matchingAlgorithm A function pointer to the algorithm used for perfect matching of stub ends
     *        (e.g., greedy, max-weight).
     * @param temperature A function pointer for a simulated annealing schedule (currently unused, but available).
     * @param theta A parameter for the Z-score calculation, influencing the penalty for gaps.
     * @param phi A parameter for weighting events by phylogenetic distance. A higher phi gives more weight
     *        to closer relatives of the reference event.
     * @param maxWalkForCalculatingZ The maximum number of segments to traverse when searching for indirect adjacencies.
     * @param ignoreUnalignedGaps If true, unaligned regions between adjacent segments are not penalized in Z-score calculation.
     * @param wiggle A factor allowing the greedy algorithm to choose a sub-optimal adjacency if it is close to the best score,
     *        adding stochasticity.
     * @param numberOfNsForScaffoldGap The number of 'N' characters to insert for a synthetic scaffold gap.
     * @param minNumberOfSequencesToSupportAdjacency The minimum number of direct adjacency observations required to
     *        prevent an adjacency from being broken.
     * @param makeScaffolds If true, enables the creation of scaffold blocks to bridge gaps.
     * @param maximumChainBasesToBreakAdjacency The maximum length, in terms of bases, of a containing chain for an
     *        adjacency in a nested flower (snarl) to be considered for breaking.
     * @param minimumNestedBasesToBreakAdjacency The minimum number of bases a flower must contain for its
     *        adjacencies to be considered for breaking.
     * @note Links within chains are considered pre-determined adjacencies and are not part of the ordering problem solved
     *       here. If a flower consists only of one or more complete chains, this function may not create any new
     *       reference intervals, as all adjacencies are already resolved.
     */

    /*
     * Get any extra ends to balance the group from the parent problem.
     */
    stList *newEnds = getExtraAttachedStubsFromParent(flower);

    /*
     * Get the chain edges.
     */
    stHash *endsToNodes = getChainNodes(flower);
    assert(stHash_size(endsToNodes) % 2 == 0);
    int64_t chainNumber = stHash_size(endsToNodes) / 2; // This number includes "trivial" chains composed of isolated blocks
    // but excludes chains that form a larger chain with parent stubs


    /*
     * Create the stub nodes.
     */
    stList *stubTangleEnds = getTangleStubEnds(flower, endsToNodes); // These are all the attached stubs from the parent,
    // excluding any that form a chain
    assert(stList_length(stubTangleEnds) % 2 == 0);
    int64_t nodeNumber = chainNumber + stList_length(stubTangleEnds);

    /*
     * Get the reference event
     */
    Event *referenceEvent = getReferenceEvent(flower, referenceEventHeader);
    st_logDebug("Chose reference event %" PRIi64 ": %s\n", event_getName(referenceEvent), event_getHeader(referenceEvent));

    /*
     * Log info about the flower
     */
    st_logDebug("For flower: %" PRIi64 " we have %" PRIi64 " nodes in reference problem for: %" PRIi64 " chains (including trivial chains) and %" PRIi64 " reference intervals %" PRIi64 "\n",
                flower_getName(flower), nodeNumber, chainNumber, stList_length(stubTangleEnds)/2);
    st_logDebug("In flower: %" PRIi64 " we have %" PRIi64 " ends, %" PRIi64 " chains (including those in chains with higher level stubs), %" PRIi64 " stubs and %" PRIi64 " blocks\n",
                flower_getName(flower), flower_getEndNumber(flower), flower_getChainNumber(flower), stList_length(stubTangleEnds), flower_getBlockNumber(flower));

    /*
     * Get the reference with chosen stub matched intervals
     */
    refOrdering *ref = getEmptyReference(flower, endsToNodes, nodeNumber, referenceEvent, matchingAlgorithm, stubTangleEnds, phi);
    assert(reference_getIntervalNumber(ref) == stList_length(stubTangleEnds) / 2);

    /*
     * Invert the hash from ends to nodes to nodes to ends.
     */
    stHash *nodesToEnds = stHash_invert(endsToNodes, (uint64_t (*)(const void *)) stIntTuple_hashKey,
            (int (*)(const void *, const void *)) stIntTuple_equalsFn, (void (*)(void *)) stIntTuple_destruct, NULL);

    /*
     * Determine which adjacencies between stubs must be preserved (i.e. scaffolded if necessary)
     */
    stList *referenceIntervalsToPreserve = NULL;
    if (makeScaffolds) {
        stHash *stubEndsToNodes = makeStubEdgesToNodesHash(stubTangleEnds, endsToNodes);
        refAdjList *stubDAL = calculateZ(flower, stubEndsToNodes, nodeNumber, 1, 1, countAdapterFn, NULL); //Gets set of adjacencies between stub ends.
        stHash_destruct(stubEndsToNodes);
        referenceIntervalsToPreserve = getReferenceIntervalsToPreserve(ref, stubDAL, minNumberOfSequencesToSupportAdjacency); //List of int-tuple pairs identifying the matchings between ends that should be preserved.
        refAdjList_destruct(stubDAL);
    }

    /*
     * Calculate z functions, using phylogenetic weighting.
     */
    stSet *chosenEvents = getEventsWithSequences(flower);
    stHash *eventWeighting = getEventWeighting(referenceEvent, phi, chosenEvents);
    stSet_destruct(chosenEvents);
    void *zArgs[2] = { &theta, eventWeighting };
    refAdjList *aL = calculateZ(flower, endsToNodes, nodeNumber, maxWalkForCalculatingZ, ignoreUnalignedGaps, calculateZScoreWeightedAdapterFn, zArgs);
    int64_t directTheta = 0.0;
    zArgs[0] = &directTheta;
    refAdjList *dAL = calculateZ(flower, endsToNodes, nodeNumber, 1, ignoreUnalignedGaps, calculateZScoreWeightedAdapterFn, zArgs); //Gets set of direct of direct adjacencies
    stHash_destruct(eventWeighting);

    /*
     * Check the edges and nodes before starting to calculate the matching.
     */
    st_logDebug("Starting to build the reference for flower %lli, with %" PRIi64 " reference intervals and %" PRIi64 " chains and %" PRIi64 " nodes in the flowers tangle\n",
                flower_getName(flower), reference_getIntervalNumber(ref), chainNumber, nodeNumber);

    double maxPossibleScore = refAdjList_getMaxPossibleScore(aL);
    makeReferenceGreedily2(aL, dAL, ref, wiggle);
    int64_t badAdjacenciesAfterGreedy = getBadAdjacencyCount(dAL, ref);
    double totalScoreAfterGreedy = getReferenceScore(aL, ref);
    st_logDebug("The score of the initial solution is %f/%" PRIi64 " out of a max possible %f\n", totalScoreAfterGreedy, badAdjacenciesAfterGreedy,
                maxPossibleScore);

    updateReferenceGreedily(aL, dAL, ref, permutations);
    int64_t badAdjacenciesAfterGreedySampling = getBadAdjacencyCount(dAL, ref);
    double totalScoreAfterGreedySampling = getReferenceScore(aL, ref);
    st_logDebug("The score of the solution after permutation sampling is %f/%" PRIi64 " after %" PRIi64 " rounds of greedy permutation out of a max possible %f\n",
                totalScoreAfterGreedySampling, badAdjacenciesAfterGreedySampling, permutations, maxPossibleScore);

    //reorderReferenceToAvoidBreakpoints(dAL2, ref);
    //int64_t badAdjacenciesAfterTopologicalReordering = getBadAdjacencyCount(dAL, ref);
    //double totalScoreAfterTopologicalReordering = getReferenceScore(aL, ref);
    //log_fn("The score of the solution after topological reordering is %f/%" PRIi64 " after %" PRIi64 " rounds of greedy permutation out of a max possible %f\n",
    //        totalScoreAfterTopologicalReordering, badAdjacenciesAfterTopologicalReordering, permutations, maxPossibleScore);

    int64_t maxNudge = 100;
    int64_t nudgePermutations = 100;
    nudgeGreedily(dAL, aL, ref, nudgePermutations, maxNudge);
    int64_t badAdjacenciesAfterNudging = getBadAdjacencyCount(dAL, ref);
    double totalScoreAfterNudging = getReferenceScore(aL, ref);

    st_logDebug("The score of the final reference solution is %f/%" PRIi64 " after %" PRIi64 " rounds of greedy nudging out of a max possible %f\n",
                totalScoreAfterNudging, badAdjacenciesAfterNudging, nudgePermutations, maxPossibleScore);

    //The aL and dAL arrays are no longer valid as we've added additional nodes to the reference, let's clean up the arrays explicitly.
    refAdjList_destruct(aL);
    refAdjList_destruct(dAL);

    /*
     * Split reference intervals where the ordering of adjacent nodes
     * is not supported by sufficient direct adjacencies, subject to conditions on the flower's
     * topology.
     */
    // The canBreakAdjacencies function is used to determined if okay to break adjacencies in the given flower
    // the logic is controlled by the number of bases in the flower and if the flower is nested within a large chain
    bool okay_to_break_adjacencies = canBreakAdjacencies(flower, minimumNestedBasesToBreakAdjacency, maximumChainBasesToBreakAdjacency);
    refAdjList *countDAL = calculateZ(flower, endsToNodes, nodeNumber, 1, 1, countAdapterFn, NULL); //Gets set of adjacencies between stub ends.
    void *extraArgs[4] = { nodesToEnds, countDAL, &minNumberOfSequencesToSupportAdjacency, &okay_to_break_adjacencies };
    stList *extraStubNodes = splitReferenceAtIndicatedLocations(ref, referenceSplitFn, extraArgs);
    assert(stList_length(extraStubNodes)%2 == 0); // Sanity check
    refAdjList_destruct(countDAL);
    stHash_destruct(endsToNodes); //Note this does not destroy the associated memory.

    /*
     * Now re-join together pairs that need to be scaffolded together.
     */
    stList *prunedExtraStubNodes;
    if (makeScaffolds) {
        prunedExtraStubNodes = remakeReferenceIntervals(ref, referenceIntervalsToPreserve, extraStubNodes);
        stList_destruct(referenceIntervalsToPreserve); //Clean this up.
    } else {
        prunedExtraStubNodes = stList_copy(extraStubNodes, NULL);
    }
    assert(stList_length(prunedExtraStubNodes)%2 == 0); // Sanity check
    st_logDebug("Making %" PRIi64 " extra stub nodes after breaking adjacencies\n", stList_length(prunedExtraStubNodes));

    if(stList_length(prunedExtraStubNodes) > 0) {
        // Add a bunch of logging information whenever we break adjacencies in the ancestor
        st_logInfo("For flower: %" PRIi64 " we have %" PRIi64 " nodes in reference problem (including ends) for: %" PRIi64 " chains (including trivial chains) and %" PRIi64 " reference intervals\n"
                   "\t\twe have %" PRIi64 " ends, %" PRIi64 " chains, %" PRIi64 " stubs and %" PRIi64 " blocks\n"
                   "\t\tThe score of the initial solution is %f/%" PRIi64 " out of a max possible %f\n"
                   "\t\tThe score of the solution after permutation sampling is %f/%" PRIi64 " after %" PRIi64 " rounds of greedy permutation out of a max possible %f\n"
                   "\t\tThe score of the final reference solution is %f/%" PRIi64 " after %" PRIi64 " rounds of greedy nudging out of a max possible %f\n"
                   "\t\tBreaking %" PRIi64 " adjacencies within %" PRIi64 " reference intervals, after rescuing from %" PRIi64 " proposed adjacency breaks\n",
           flower_getName(flower), nodeNumber, chainNumber, stList_length(stubTangleEnds)/2,
           flower_getEndNumber(flower), flower_getChainNumber(flower), stList_length(stubTangleEnds), flower_getBlockNumber(flower),
           totalScoreAfterGreedy, badAdjacenciesAfterGreedy, maxPossibleScore,
           totalScoreAfterGreedySampling, badAdjacenciesAfterGreedySampling, permutations, maxPossibleScore,
           totalScoreAfterNudging, badAdjacenciesAfterNudging, nudgePermutations, maxPossibleScore,
           stList_length(prunedExtraStubNodes)/2, stList_length(stubTangleEnds)/2, stList_length(extraStubNodes)/2);
    }

    /*
     * Convert the additional stub nodes into new stub ends, updating the endsToNodes and nodesToEnds sets.
     */
    addAdditionalStubEnds(prunedExtraStubNodes, flower, nodesToEnds, newEnds);
    stList_destruct(prunedExtraStubNodes);

    /*
     * Convert the reference into a list of adjacency edges.
     */
    stList *chosenEdges = convertReferenceToAdjacencyEdges2(ref);

    /*
     * Check the matching we have.
     */
    assert(stList_length(chosenEdges) * 2 == stHash_size(nodesToEnds));
    /*stList *nodes = stHash_getValues(nodesToEnds);
    stSortedSet *nodesSet = stList_getSortedSet(nodes, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
    assert(stHash_size(nodesToEnds) == stSortedSet_size(nodesSet));
    stSortedSet_destruct(nodesSet);
    stList_destruct(nodes);*/

    /*
     * Add the reference genome into flower
     */
    makeReferenceThreads(flower, chosenEdges, nodesToEnds, referenceEvent, numberOfNsForScaffoldGap);

    /*
     * Ensure the newly created ends have a group.
     */
    assignGroups(newEnds, flower, referenceEvent);

    /*
     * Cleanup
     */
    stList_destruct(newEnds);
    stHash_destruct(nodesToEnds);
    stList_destruct(chosenEdges);
    reference_destruct(ref);
    stList_destruct(stubTangleEnds);
    stList_destruct(extraStubNodes);
}

////////////////////////////////////
////////////////////////////////////
//Overall coordination function
////////////////////////////////////
////////////////////////////////////

void cactus_make_reference(stList *flowers, char *referenceEventString,
                           CactusDisk *cactusDisk, CactusParams *params) {
    ///////////////////////////////////////////////////////////////////////////
    // Build the reference
    ///////////////////////////////////////////////////////////////////////////

    int64_t permutations = cactusParams_get_int(params, 2, "reference", "permutations");
    double theta = cactusParams_get_float(params, 2, "reference", "theta");
    double phi = cactusParams_get_float(params, 2, "reference", "phi");
    bool useSimulatedAnnealing = cactusParams_get_int(params, 2, "reference", "useSimulatedAnnealing");
    int64_t maxWalkForCalculatingZ = cactusParams_get_int(params, 2, "reference", "maxWalkForCalculatingZ");
    bool ignoreUnalignedGaps = cactusParams_get_int(params, 2, "reference", "ignoreUnalignedGaps");
    double wiggle = cactusParams_get_float(params, 2, "reference", "wiggle");
    int64_t numberOfNsForScaffoldGap = cactusParams_get_int(params, 2, "reference", "numberOfNs");
    int64_t minNumberOfSequencesToSupportAdjacency = cactusParams_get_int(params, 2, "reference", "minNumberOfSequencesToSupportAdjacency");
    int64_t minimumNestedBasesToBreakAdjacency = cactusParams_get_int(params, 2, "reference", "minimumNestedBasesToBreakAdjacency");
    int64_t maximumChainBasesToBreakAdjacency = cactusParams_get_int(params, 2, "reference", "maximumChainBasesToBreakAdjacency");
    bool makeScaffolds = cactusParams_get_int(params, 2, "reference", "makeScaffolds");

    stList *(*matchingAlgorithm)(stList *edges, int64_t nodeNumber) = chooseMatching_greedy;
    char *matchAlgorithmString = cactusParams_get_string(params, 2, "reference", "matchingAlgorithm");
    if (strcmp("greedy", matchAlgorithmString) == 0) {
        matchingAlgorithm = chooseMatching_greedy;
    } else if (strcmp("maxCardinality", matchAlgorithmString) == 0) {
        matchingAlgorithm = chooseMatching_maximumCardinalityMatching;
    } else if (strcmp("maxWeight", matchAlgorithmString) == 0) {
        matchingAlgorithm = chooseMatching_maximumWeightMatching;
    } else if (strcmp("blossom5", matchAlgorithmString) == 0) {
        matchingAlgorithm = chooseMatching_blossom5;
    } else {
        stThrowNew(REFERENCE_BUILDING_EXCEPTION, "Input error: unrecognized matching algorithm: %s", matchAlgorithmString);
    }
    free(matchAlgorithmString);

    /*st_logDebug("The reference event string: %s\n", referenceEventString);
    st_logDebug("The theta parameter has been set to %lf\n", theta);
    st_logDebug("The ignore unaligned gaps parameter is %i\n", ignoreUnalignedGaps);
    st_logDebug("The number of permutations is %" PRIi64 "\n", permutations);
    st_logDebug("Simulated annealing is %" PRIi64 "\n", useSimulatedAnnealing);
    st_logDebug("Max number of segments in thread to calculate z-score between is %" PRIi64 "\n",
            maxWalkForCalculatingZ);
    st_logDebug("Wiggle is %f\n", wiggle);
    st_logDebug("Max number of Ns for a scaffold gap is: %" PRIi64 "\n", numberOfNsForScaffoldGap);
    st_logDebug("Min number of sequences to required to support an adjacency is: %" PRIi64 "\n",
            minNumberOfSequencesToSupportAdjacency);
    st_logDebug("Make scaffolds is: %i\n", makeScaffolds);*/

    double (*temperatureFn)(double) = useSimulatedAnnealing ? exponentiallyDecreasingTemperatureFn : constantTemperatureFn;

#pragma omp parallel for schedule(dynamic, 1)
    for(int64_t i=0; i<stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        st_logDebug("Processing flower %" PRIi64 "\n", flower_getName(flower));
        buildReferenceTopDown(flower, referenceEventString, permutations, matchingAlgorithm, temperatureFn, theta,
                              phi, maxWalkForCalculatingZ, ignoreUnalignedGaps, wiggle, numberOfNsForScaffoldGap,
                              minNumberOfSequencesToSupportAdjacency, makeScaffolds,
                              minimumNestedBasesToBreakAdjacency, maximumChainBasesToBreakAdjacency);
    }
}
