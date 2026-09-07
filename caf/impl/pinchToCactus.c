#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stCaf.h"

///////////////////////////////////////////////////////////////////////////
// Construct dead end component
///////////////////////////////////////////////////////////////////////////

static void attachPinchBlockEndToAnotherComponent(stPinchBlock *pinchBlock, bool orientation, stList *anotherComponent) {
    assert(pinchBlock != NULL);
    assert(stPinchBlock_getLength(pinchBlock) == 1);
    stPinchEnd *pinchEnd = stPinchBlock_getEnd(pinchBlock, orientation);
    assert(pinchEnd != NULL);
    stList *component = stPinchEnd_getComponent(pinchEnd);
    assert(component != NULL);
    assert(stList_length(component) == 1);
    assert(stList_get(component, 0) == pinchEnd);
    // The emptied list stays where it is in the list of components so that the order the remaining
    // components are visited in does not change; stCaf_getCactusGraphForThreadSet frees the empties.
    stList_pop(component);
    stList_append(anotherComponent, pinchEnd);
    stPinchEnd_setComponent(pinchEnd, anotherComponent);
}

static stList *stCaf_constructDeadEndComponent(Flower *flower, stPinchThreadSet *threadSet) {
    /*
     * Locates the ends of all the attached ends and merges together their 'dead end' components to create a single
     * 'dead end' component, as described in the JCB cactus paper.
     */
    //For each block end at the end of a thread, attach to dead end component if associated end is attached
    stList *deadEndAdjacencyComponent = stList_construct(); //the ends belong to the thread set's block end records
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *pinchThread;
    while ((pinchThread = stPinchThreadSetIt_getNext(&threadIt))) {
        Cap *cap = flower_getCap(flower, stPinchThread_getName(pinchThread));
        assert(cap != NULL);
        End *end1 = cap_getEnd(cap), *end2 = cap_getEnd(cap_getAdjacency(cap));
        assert(end1 != NULL && end2 != NULL);
        if (end_isAttached(end1)) {
            stPinchSegment *pinchSegment = stPinchThread_getFirst(pinchThread);
            stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
            if (stPinchBlock_getFirst(pinchBlock) == pinchSegment) { //We only want to do this once
                attachPinchBlockEndToAnotherComponent(pinchBlock, stPinchSegment_getBlockOrientation(pinchSegment),
                        deadEndAdjacencyComponent);
            }
        }
        if (end_isAttached(end2)) {
            stPinchSegment *pinchSegment = stPinchThread_getLast(pinchThread);
            stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
            if (stPinchBlock_getFirst(pinchBlock) == pinchSegment) { //And only once for the other end
                attachPinchBlockEndToAnotherComponent(pinchBlock, !stPinchSegment_getBlockOrientation(pinchSegment),
                        deadEndAdjacencyComponent);
            }
        }
    }
    return deadEndAdjacencyComponent;
}

///////////////////////////////////////////////////////////////////////////
// Attach unatttached thread components
///////////////////////////////////////////////////////////////////////////

static bool threadIsAttachedToDeadEndComponent5Prime(stPinchThread *thread, stList *deadEndComponent) {
    stPinchSegment *pinchSegment = stPinchThread_getFirst(thread);
    stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
    assert(pinchBlock != NULL);
    stPinchEnd *pinchEnd = stPinchBlock_getEnd(pinchBlock, stPinchSegment_getBlockOrientation(pinchSegment));
    assert(pinchEnd != NULL && stPinchEnd_getComponent(pinchEnd) != NULL);
    return stPinchEnd_getComponent(pinchEnd) == deadEndComponent;
}

static bool threadIsAttachedToDeadEndComponent3Prime(stPinchThread *thread, stList *deadEndComponent) {
    stPinchSegment *pinchSegment = stPinchThread_getLast(thread);
    stPinchBlock *pinchBlock = stPinchSegment_getBlock(pinchSegment);
    assert(pinchBlock != NULL);
    stPinchEnd *pinchEnd = stPinchBlock_getEnd(pinchBlock, !stPinchSegment_getBlockOrientation(pinchSegment));
    assert(pinchEnd != NULL && stPinchEnd_getComponent(pinchEnd) != NULL);
    return stPinchEnd_getComponent(pinchEnd) == deadEndComponent;
}

static bool threadIsAttachedToDeadEndComponent(stPinchThread *thread, stList *deadEndComponent) {
    return threadIsAttachedToDeadEndComponent5Prime(thread, deadEndComponent)
            || threadIsAttachedToDeadEndComponent3Prime(thread, deadEndComponent);
}

static void attachThreadToDeadEndComponent(stPinchThread *thread, stList *deadEndAdjacencyComponent,
        bool markEndsAttached, Flower *flower) {
    stPinchSegment *segment = stPinchThread_getFirst(thread);
    attachPinchBlockEndToAnotherComponent(stPinchSegment_getBlock(segment), stPinchSegment_getBlockOrientation(segment),
            deadEndAdjacencyComponent);
    segment = stPinchThread_getLast(thread);
    attachPinchBlockEndToAnotherComponent(stPinchSegment_getBlock(segment), !stPinchSegment_getBlockOrientation(segment),
            deadEndAdjacencyComponent);
    if (markEndsAttached) { //Get the ends and attach them
        Cap *cap = flower_getCap(flower, stPinchSegment_getName(stPinchThread_getFirst(thread))); //The following three lines isolates the sequence associated with a segment.
        assert(cap != NULL);
        end_makeAttached(cap_getEnd(cap));
        end_makeAttached(cap_getEnd(cap_getAdjacency(cap)));
    }
}

static int comparePinchThreadsByName(const void *a, const void *b) {
    int64_t i = stPinchThread_getName((stPinchThread *) a);
    int64_t j = stPinchThread_getName((stPinchThread *) b);
    return i > j ? 1 : (i < j ? -1 : 0);
}

int comparePinchThreadsByLength(const void *a, const void *b) {
    int64_t i = stPinchThread_getLength((stPinchThread *) a);
    int64_t j = stPinchThread_getLength((stPinchThread *) b);
    if (i != j) {
        return i > j ? 1 : -1;
    }
    return comparePinchThreadsByName(a, b); // Break ties by name: which of two equally long
    // threads gets attached first has to be decided the same way on every run
}

static int compareThreadComponentsByFirstThread(const void *a, const void *b) {
    // The components are disjoint and thread names are unique, so comparing their
    // lowest-named thread is a total order over the components
    return comparePinchThreadsByName(stList_get((stList *) a, 0), stList_get((stList *) b, 0));
}

static void attachThreadComponentToDeadEndComponent(stList *threadComponent, stList *deadEndComponent,
        bool markEndsAttached, int64_t minLengthForChromosome,
        double proportionOfUnalignedBasesForNewChromosome, Flower *flower, int64_t *basesAligned) {
    /*
     * Algorithm walks the threads in the connected component, in descending order of length and
     * for each thread attaches it to the dead-end component if its length of bases in blocks not contained in chromosome
     * length fragments.
     *
     * basesAligned is indexed by thread index and shared between the components: the blocks of a component only
     * ever hold that component's threads, so each entry is touched by one component and can start at zero.
     * Blocks already counted are marked in the data slot of their first end, which the caller clears afterwards.
     */

    //First get threads in order, already attached first then unattached, sorted by descending length
    bool first = 1; //This is used to ensure at least one thread is attached in component.
    stList *l = stList_construct();
    stList *l2 = stList_construct();
    for (int64_t i = 0; i < stList_length(threadComponent); i++) {
        stPinchThread *pinchThread = stList_get(threadComponent, i);
        if (threadIsAttachedToDeadEndComponent(pinchThread, deadEndComponent)) {
            first = 0;
            stList_append(l2, pinchThread);
        } else {
            stList_append(l, pinchThread);
        }
    }
    stList_sort(l, comparePinchThreadsByLength);
    stList_appendAll(l, l2); // Add any already attached threads to the end (so they are considered first)
    stList_destruct(l2);

    //Now iterate on threads, attaching as needed.
    while (stList_length(l) > 0) {
        stPinchThread *pinchThread = stList_pop(l);
        if(markEndsAttached) {
            Cap *cap = flower_getCap(flower, stPinchSegment_getName(stPinchThread_getFirst(pinchThread))); //The following three lines isolates the sequence associated with a segment.
            assert(cap != NULL);
            Sequence *sequence = cap_getSequence(cap);
            assert(sequence != NULL);
            st_logDebug("Considering attaching the following sequence to the cactus root %"
            PRIi64 ", header %s with length %" PRIi64
            ", is already attached: %s, have already attached something: %s\n",
                    sequence_getName(sequence), sequence_getHeader(sequence), sequence_getLength(sequence),
                    threadIsAttachedToDeadEndComponent(pinchThread, deadEndComponent) ?
                    "True" : "False", first ? "False" : "True");
        }
        if (stPinchThread_getLength(pinchThread) < minLengthForChromosome && !first) { // If too short and nothing in the component is attached
            continue;
        }
        int64_t basesAlignedToChromosomeThreads = basesAligned[stPinchThread_getIndex(pinchThread)]; //This is the number of bases already aligned in threads with attached ends (chromosomes);
        assert(basesAlignedToChromosomeThreads >= 0);
        //Walk the thread
        stPinchSegment *segment = stPinchThread_getFirst(pinchThread);
        assert(segment != NULL);
        assert(stPinchSegment_get5Prime(segment) == NULL);
        do {
            stPinchBlock *block;
            if ((block = stPinchSegment_getBlock(segment)) != NULL) {
                stPinchEnd *seenMark = stPinchBlock_getEnd(block, 0);
                assert(seenMark != NULL);
                if (stPinchEnd_getData(seenMark) == NULL) {
                    stPinchEnd_setData(seenMark, block); //any non-NULL value marks the block as counted
                    stPinchBlockIt segIt = stPinchBlock_getSegmentIterator(block);
                    stPinchSegment *segment2;
                    while ((segment2 = stPinchBlockIt_getNext(&segIt)) != NULL) {
                        basesAligned[stPinchThread_getIndex(stPinchSegment_getThread(segment2))] += stPinchBlock_getLength(block);
                    }
                }
            }
            segment = stPinchSegment_get3Prime(segment);
        } while (segment != NULL);
        if (threadIsAttachedToDeadEndComponent(pinchThread, deadEndComponent)) { //If this is already attached we can stop at this point
            continue;
        }
        int64_t totalBasesAligned = basesAligned[stPinchThread_getIndex(pinchThread)]; //This is the number of bases already aligned in chromosomes;
        assert(totalBasesAligned >= basesAlignedToChromosomeThreads);
        if((totalBasesAligned - basesAlignedToChromosomeThreads) >= proportionOfUnalignedBasesForNewChromosome * totalBasesAligned || first) { // Attach if sufficiently distinct or nothing is yet attached
            attachThreadToDeadEndComponent(pinchThread, deadEndComponent, markEndsAttached, flower);
            first = 0; // We have officially attached an end in the component
            if(markEndsAttached) {
                Cap *cap = flower_getCap(flower, stPinchSegment_getName(stPinchThread_getFirst(pinchThread))); //The following three lines isolates the sequence associated with a segment.
                assert(cap != NULL);
                Sequence *sequence = cap_getSequence(cap);
                assert(sequence != NULL);
                if(stPinchThread_getLength(pinchThread) > minLengthForChromosome) { // Just log the longer sequences
                    st_logInfo("Attaching the sequence to the cactus root %" PRIi64 ", header %s with length %" PRIi64 " and %" PRIi64 " total bases aligned and %" PRIi64 " bases aligned to other chromosome threads\n",
                                sequence_getName(sequence), sequence_getHeader(sequence), sequence_getLength(sequence), totalBasesAligned, basesAlignedToChromosomeThreads);
                }
            }
        }
    }
    stList_destruct(l);
}

static void stCaf_attachUnattachedThreadComponents(Flower *flower, stPinchThreadSet *threadSet, stList *deadEndComponent,
        bool markEndsAttached, int64_t minLengthForChromosome,
        double proportionOfUnalignedBasesForNewChromosome) {
    /*
     * Locates threads components which have no dead ends part of the dead end component, and then
     * connects them, picking the longest thread to attach them.
     */
    if (flower_getName(flower) != 0) {
        return; //We don't need to attach ends if flower is not at top level, as everything is anchored to attached components in lower level flowers.
    }
    stSortedSet *threadComponents = stPinchThreadSet_getThreadComponents(threadSet);
    assert(stSortedSet_size(threadComponents) > 0);
    /*
     * The set holds the components in address order, and each component holds its
     * threads in whatever order the union-find produced them.  Which threads get
     * attached to the root depends on the order they are visited in, so put both
     * into name order first, which is the same on every run.
     */
    stList *threadComponents2 = stSortedSet_getList(threadComponents);
    for (int64_t i = 0; i < stList_length(threadComponents2); i++) {
        stList_sort(stList_get(threadComponents2, i), comparePinchThreadsByName);
    }
    stList_sort(threadComponents2, compareThreadComponentsByFirstThread);
    int64_t *basesAligned = st_calloc(stPinchThreadSet_getSize(threadSet), sizeof(int64_t));
    for (int64_t i = 0; i < stList_length(threadComponents2); i++) {
        attachThreadComponentToDeadEndComponent(stList_get(threadComponents2, i), deadEndComponent,
                markEndsAttached,
                minLengthForChromosome, proportionOfUnalignedBasesForNewChromosome, flower, basesAligned);
    }
    free(basesAligned);
    stPinchThreadSet_clearEndData(threadSet); //the blocks-seen marks; the cactus graph construction needs the slots empty
    stList_destruct(threadComponents2);
    stSortedSet_destruct(threadComponents);
}

///////////////////////////////////////////////////////////////////////////
// Create a cactus graph from a pinch graph
///////////////////////////////////////////////////////////////////////////

/*
 * During construction the data slot of every end holds the cactus node of the end's adjacency component.
 */
static void setCactusNodeForComponent(stList *adjacencyComponent, stCactusNode *cactusNode) {
    for (int64_t i = 0; i < stList_length(adjacencyComponent); i++) {
        stPinchEnd_setData(stList_get(adjacencyComponent, i), cactusNode);
    }
}

static stCactusNode *getCactusNode(stPinchEnd *pinchEnd) {
    stCactusNode *cactusNode = stPinchEnd_getData(pinchEnd);
    assert(cactusNode != NULL);
    return cactusNode;
}

void *stCaf_mergeNodeObjects(void *a, void *b) {
    stList *adjacencyComponents1 = a;
    stList *adjacencyComponents2 = b;
    assert(adjacencyComponents1 != adjacencyComponents2);
    stList_appendAll(adjacencyComponents1, adjacencyComponents2);
    stList_setDestructor(adjacencyComponents2, NULL);
    stList_destruct(adjacencyComponents2);
    return adjacencyComponents1;
}

static void *makeNodeObject(stList *adjacencyComponent) {
    stList *adjacencyComponents = stList_construct3(0, (void(*)(void *)) stList_destruct);
    stList_append(adjacencyComponents, adjacencyComponent);
    return adjacencyComponents;
}

static bool isDeadEndStubComponent(stList *adjacencyComponent, stPinchEnd *pinchEnd) {
    if (stList_length(adjacencyComponent) != 1) {
        return 0;
    }
    stPinchSegment *pinchSegment = stPinchBlock_getFirst(stPinchEnd_getBlock(pinchEnd));
    return (stPinchEnd_traverse5Prime(stPinchEnd_getOrientation(pinchEnd), pinchSegment) ? stPinchSegment_get5Prime(pinchSegment)
            : stPinchSegment_get3Prime(pinchSegment)) == NULL;
}

/*
 * stCactusGraph_breakChainsByEndsNotInChains rescans a chain from its start after every merge it makes,
 * so on a long chain the predicate is asked the same question many times over. The answer depends only
 * on the pinch graph, which does not change during the pass, and on which block the end is currently
 * linked to, so it is remembered in the end's data slot as that block's pointer tagged with the answer
 * (blocks are 16 byte aligned, so the low bits are free; the tag bit tells a memo from anything else
 * the slot may have held). The counters go into the caf-timing line.
 */
static int64_t reversalAtEndCalls = 0, reversalAtEndCached = 0;

static bool stCaf_reversalAtEnd(stCactusEdgeEnd *cactusEdgeEnd, void *extraArg) {
    assert(extraArg == NULL);
    assert(stCactusEdgeEnd_getObject(cactusEdgeEnd) != NULL);
    assert(stCactusEdgeEnd_getLink(cactusEdgeEnd) != NULL);
    assert(stCactusEdgeEnd_getObject(stCactusEdgeEnd_getLink(cactusEdgeEnd)) != NULL);
    stPinchEnd *end = stCactusEdgeEnd_getObject(cactusEdgeEnd);
    stPinchBlock *otherBlock = stPinchEnd_getBlock(stCactusEdgeEnd_getObject(stCactusEdgeEnd_getLink(cactusEdgeEnd)));
    assert(((uintptr_t) otherBlock & (uintptr_t) 3) == 0);
    uintptr_t memo = (uintptr_t) stPinchEnd_getData(end);
    reversalAtEndCalls++;
    if ((memo & (uintptr_t) 2) && (memo & ~(uintptr_t) 3) == (uintptr_t) otherBlock) {
        reversalAtEndCached++;
        return memo & (uintptr_t) 1;
    }
    bool reversal = stPinchEnd_hasSelfLoopWithRespectToOtherBlock(end, otherBlock);
    stPinchEnd_setData(end, (void *) ((uintptr_t) otherBlock | (uintptr_t) 2 | (reversal ? (uintptr_t) 1 : (uintptr_t) 0)));
    return reversal;
}

static bool stCaf_breakAtTooGreatEnoughMedianSeparation(stCactusEdgeEnd *cactusEdgeEnd, void *extraArg) {
    int64_t maximumMedianSpacingBetweenLinkedEnds = *(int64_t *)extraArg;
    assert(maximumMedianSpacingBetweenLinkedEnds >= 0 && maximumMedianSpacingBetweenLinkedEnds < INT64_MAX);
    int64_t median = stPinchEnd_getMedianSubSequenceLengthConnectingEnds(stCactusEdgeEnd_getObject(cactusEdgeEnd), stCactusEdgeEnd_getObject(stCactusEdgeEnd_getLink(cactusEdgeEnd)));
    return median >= 0 && median > maximumMedianSpacingBetweenLinkedEnds;
}

/*
 * Per-phase wall-clock times of one pinch graph -> cactus graph build, for the caf-timing log line.
 */
typedef struct _stCafBuildTiming {
    double adjacency, deadEnd, attach, build, collapse, bridges, tandems, median;
    int64_t ends, nodesBeforeCollapse, nodesAfterCollapse, tandemCalls, tandemCached;
} stCafBuildTiming;

/*
 * Debugging aid: with CACTUS_CAF_DUMP_CACTUS=<dir> set, every cactus graph built is written to <dir>/cactus-<n>.txt in
 * iteration order, so that a change which alters iteration order (which feeds the names in the output) shows up as a diff
 * even before it changes the c2h.
 */
static void stCaf_dumpCactusGraph(stCactusGraph *cactusGraph, stCactusNode *startCactusNode) {
    const char *dir = getenv("CACTUS_CAF_DUMP_CACTUS");
    if (dir == NULL) {
        return;
    }
    static int64_t buildCount = 0;
    char *path = stString_print("%s/cactus-%" PRIi64 ".txt", dir, buildCount++);
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        st_errAbort("Could not open %s for the cactus graph dump", path);
    }
    stCactusGraphNodeIt *nodeIt = stCactusGraphNodeIterator_construct(cactusGraph);
    stCactusNode *cactusNode;
    int64_t nodeIndex = 0;
    while ((cactusNode = stCactusGraphNodeIterator_getNext(nodeIt)) != NULL) {
        fprintf(f, "node %" PRIi64 "%s\n", nodeIndex++, cactusNode == startCactusNode ? " start" : "");
        stCactusNodeEdgeEndIt edgeEndIt = stCactusNode_getEdgeEndIt(cactusNode);
        stCactusEdgeEnd *edgeEnd;
        while ((edgeEnd = stCactusNodeEdgeEndIt_getNext(&edgeEndIt)) != NULL) {
            stPinchEnd *end = stCactusEdgeEnd_getObject(edgeEnd);
            stPinchSegment *segment = stPinchBlock_getFirst(stPinchEnd_getBlock(end));
            fprintf(f, " %" PRIi64 ":%" PRIi64 ":%" PRIi64 ":%i lo=%i ce=%i", stPinchSegment_getName(segment), stPinchSegment_getStart(segment),
                    stPinchBlock_getLength(stPinchEnd_getBlock(end)), stPinchEnd_getOrientation(end),
                    stCactusEdgeEnd_getLinkOrientation(edgeEnd), stCactusEdgeEnd_isChainEnd(edgeEnd));
            stCactusEdgeEnd *link = stCactusEdgeEnd_getLink(edgeEnd);
            if (link != NULL) {
                stPinchEnd *linkEnd = stCactusEdgeEnd_getObject(link);
                stPinchSegment *linkSegment = stPinchBlock_getFirst(stPinchEnd_getBlock(linkEnd));
                fprintf(f, " link=%" PRIi64 ":%" PRIi64 ":%i", stPinchSegment_getName(linkSegment), stPinchSegment_getStart(linkSegment),
                        stPinchEnd_getOrientation(linkEnd));
            }
            fprintf(f, "\n");
        }
    }
    stCactusGraphNodeIterator_destruct(nodeIt);
    fclose(f);
    free(path);
}

static void makeCactusNodeForEnd(stCactusGraph *cactusGraph, stPinchEnd *pinchEnd) {
    if (stPinchEnd_getData(pinchEnd) != NULL) { //The end's component already has a node
        return;
    }
    stList *adjacencyComponent = stPinchEnd_getComponent(pinchEnd);
    assert(adjacencyComponent != NULL);
    if (isDeadEndStubComponent(adjacencyComponent, pinchEnd)) { //Going to be a bridge to nowhere, so we join it - this ensures
        //that all dead end nodes of free stubs end up in the same node as their non-dead end counterparts.
        assert(pinchEnd == stList_get(adjacencyComponent, 0));
        stPinchEnd *otherPinchEnd = stPinchEnd_getOtherEnd(pinchEnd);
        stList *otherAdjacencyComponent = stPinchEnd_getComponent(otherPinchEnd);
        assert(otherAdjacencyComponent != NULL && adjacencyComponent != otherAdjacencyComponent);
        stCactusNode *cactusNode = stPinchEnd_getData(otherPinchEnd);
        if (cactusNode == NULL) {
            cactusNode = stCactusNode_construct(cactusGraph, makeNodeObject(otherAdjacencyComponent));
            setCactusNodeForComponent(otherAdjacencyComponent, cactusNode);
        }
        stList *adjacencyComponents = stCactusNode_getObject(cactusNode);
        stList_append(adjacencyComponents, adjacencyComponent);
        setCactusNodeForComponent(adjacencyComponent, cactusNode);
    } else {
        setCactusNodeForComponent(adjacencyComponent, stCactusNode_construct(cactusGraph, makeNodeObject(adjacencyComponent)));
    }
}

static void makeCactusEdgeForEnd(stCactusGraph *cactusGraph, stPinchEnd *pinchEnd) {
    if (stPinchEnd_getOrientation(pinchEnd)) { //Assure we make the edge only once
        assert(stPinchEnd_getBlock(pinchEnd) != NULL);
        stPinchEnd *pinchEnd2 = stPinchEnd_getOtherEnd(pinchEnd);
        assert(pinchEnd != pinchEnd2);
        stCactusEdgeEnd_construct(cactusGraph, getCactusNode(pinchEnd), getCactusNode(pinchEnd2), pinchEnd, pinchEnd2);
    }
}

static stCactusGraph *stCaf_constructCactusGraph(stPinchThreadSet *threadSet, stList *deadEndComponent, stList *adjacencyComponents,
        stCactusNode **startCactusNode, bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds,
        stCafBuildTiming *timing) {
    /*
     * Constructs a cactus graph from a set of pinch graph components, including the dead end component. Returns a cactus
     * graph, and assigns 'startCactusNode' to the cactus node containing the dead end component.
     *
     * The order nodes and edges are made in decides the order the graph is later traversed in, and with it the names
     * in the output, so it has to stay what it was when the ends were kept in an insertion-ordered hash: components
     * in discovery order, each end in the order it was discovered, and the ends moved into the dead end component
     * last, in the order they were moved.
     */
    double t = stCaf_now();
    stCactusGraph *cactusGraph = stCactusGraph_construct2((void(*)(void *)) stList_destruct, NULL);

    //Make the nodes
    *startCactusNode = stCactusNode_construct(cactusGraph, makeNodeObject(deadEndComponent));
    setCactusNodeForComponent(deadEndComponent, *startCactusNode);
    int64_t ends = stList_length(deadEndComponent);
    for (int64_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        ends += stList_length(adjacencyComponent);
        for (int64_t j = 0; j < stList_length(adjacencyComponent); j++) {
            makeCactusNodeForEnd(cactusGraph, stList_get(adjacencyComponent, j));
        }
    }
    for (int64_t j = 0; j < stList_length(deadEndComponent); j++) {
        makeCactusNodeForEnd(cactusGraph, stList_get(deadEndComponent, j));
    }

    //Make the edges
    for (int64_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        for (int64_t j = 0; j < stList_length(adjacencyComponent); j++) {
            makeCactusEdgeForEnd(cactusGraph, stList_get(adjacencyComponent, j));
        }
    }
    for (int64_t j = 0; j < stList_length(deadEndComponent); j++) {
        makeCactusEdgeForEnd(cactusGraph, stList_get(deadEndComponent, j));
    }
    timing->ends = ends;
    timing->nodesBeforeCollapse = stCactusGraph_getNodeNumber(cactusGraph);
    timing->build = stCaf_now() - t;
    t = stCaf_now();

    //Run the cactus-ifying functions
    stCactusGraph_collapseToCactus(cactusGraph, stCaf_mergeNodeObjects, *startCactusNode);
    timing->collapse = stCaf_now() - t;
    t = stCaf_now();
    stCactusGraph_collapseBridges(cactusGraph, *startCactusNode, stCaf_mergeNodeObjects);
    timing->bridges = stCaf_now() - t;
    t = stCaf_now();
    reversalAtEndCalls = reversalAtEndCached = 0;
    if(breakChainsAtReverseTandems) {
        stPinchThreadSet_clearEndData(threadSet); //the slots held the cactus nodes, no longer needed; the predicate memoises in them
        *startCactusNode = stCactusGraph_breakChainsByEndsNotInChains(cactusGraph, *startCactusNode, stCaf_mergeNodeObjects, stCaf_reversalAtEnd, NULL);
    }
    timing->tandems = stCaf_now() - t;
    timing->tandemCalls = reversalAtEndCalls;
    timing->tandemCached = reversalAtEndCached;
    t = stCaf_now();
    if(maximumMedianSpacingBetweenLinkedEnds < INT64_MAX) {
        *startCactusNode = stCactusGraph_breakChainsByEndsNotInChains(cactusGraph, *startCactusNode, stCaf_mergeNodeObjects, stCaf_breakAtTooGreatEnoughMedianSeparation, &maximumMedianSpacingBetweenLinkedEnds);
    }
    timing->median = stCaf_now() - t;
    timing->nodesAfterCollapse = stCactusGraph_getNodeNumber(cactusGraph);

    return cactusGraph;
}

///////////////////////////////////////////////////////////////////////////
// Function that draws together above functions to generate a cactus graph from a pinch graph.
///////////////////////////////////////////////////////////////////////////

stCactusGraph *stCaf_getCactusGraphForThreadSet(Flower *flower, stPinchThreadSet *threadSet, stCactusNode **startCactusNode,
        stList **deadEndComponent, bool attachEndsInFlower, int64_t minLengthForChromosome,
        double proportionOfUnalignedBasesForNewChromosome,
        bool breakChainsAtReverseTandems, int64_t maximumMedianSpacingBetweenLinkedEnds) {
    stCafBuildTiming timing;
    double t = stCaf_now();

    //Get adjacency components, recorded on the block ends. The records live until stCaf_destructCactusGraph.
    stPinchThreadSet_attachEnds(threadSet);
    stList *adjacencyComponents = stPinchThreadSet_getAdjacencyComponents(threadSet);
    timing.adjacency = stCaf_now() - t;
    t = stCaf_now();

    //Merge together dead end component
    *deadEndComponent = stCaf_constructDeadEndComponent(flower, threadSet);
    timing.deadEnd = stCaf_now() - t;
    t = stCaf_now();

    //Join unattached components of graph by dead ends to dead end component, and make other ends 'attached' if necessary
    stCaf_attachUnattachedThreadComponents(flower, threadSet, *deadEndComponent, attachEndsInFlower,
            minLengthForChromosome, proportionOfUnalignedBasesForNewChromosome);
    timing.attach = stCaf_now() - t;

    //Create cactus
    stCactusGraph *cactusGraph = stCaf_constructCactusGraph(threadSet, *deadEndComponent, adjacencyComponents, startCactusNode,
            breakChainsAtReverseTandems, maximumMedianSpacingBetweenLinkedEnds, &timing);

    //Cleanup: the non-empty components now belong to the cactus nodes; the ones emptied by the dead end merging do not
    for (int64_t i = 0; i < stList_length(adjacencyComponents); i++) {
        stList *adjacencyComponent = stList_get(adjacencyComponents, i);
        if (stList_length(adjacencyComponent) == 0) {
            stList_destruct(adjacencyComponent);
        }
    }
    stList_setDestructor(adjacencyComponents, NULL);
    stList_destruct(adjacencyComponents);

    st_logInfo("caf-timing: cactus-graph ends=%" PRIi64 " adjacency %.3fs deadend %.3fs attach %.3fs build %.3fs collapse %.3fs bridges %.3fs tandems %.3fs median %.3fs nodes %" PRIi64 "->%" PRIi64 " tandem-calls %" PRIi64 " cached %" PRIi64 "\n",
               timing.ends, timing.adjacency, timing.deadEnd, timing.attach, timing.build, timing.collapse, timing.bridges, timing.tandems, timing.median,
               timing.nodesBeforeCollapse, timing.nodesAfterCollapse, timing.tandemCalls, timing.tandemCached);
    stCaf_dumpCactusGraph(cactusGraph, *startCactusNode);

    return cactusGraph;
}

void stCaf_destructCactusGraph(stCactusGraph *cactusGraph, stPinchThreadSet *threadSet) {
    stCactusGraph_destruct(cactusGraph);
    stPinchThreadSet_detachEnds(threadSet);
}
