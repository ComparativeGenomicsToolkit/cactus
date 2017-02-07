#include "cactus.h"
#include "sonLib.h"
#include "commonC.h"
#include "stCaf.h"

// This global is a bit gross but needs to be used to work around the
// fact that the pinch filter functions don't have an "extra data"
// parameter.
static Flower *flower;

void stCaf_setFlowerForAlignmentFiltering(Flower *input) {
    flower = input;
}

/*
 * Functions used for prefiltering the alignments.
 */

Event *stCaf_getEvent(stPinchSegment *segment, Flower *flower) {
    Event *event = cap_getEvent(flower_getCap(flower, stPinchSegment_getName(segment)));
    assert(event != NULL);
    return event;
}

/*
 * Filtering by presence of outgroup. This code is efficient and scales linearly with depth.
 */

static bool containsOutgroupSegment(stPinchBlock *block, Flower *flower) {
    stPinchBlockIt it = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&it)) != NULL) {
        if (event_isOutgroup(stCaf_getEvent(segment, flower))) {
            stPinchSegment_putSegmentFirstInBlock(segment);
            assert(stPinchBlock_getFirst(block) == segment);
            return 1;
        }
    }
    return 0;
}

static bool isOutgroupSegment(stPinchSegment *segment, Flower *flower) {
    if (event_isOutgroup(stCaf_getEvent(segment, flower))) {
        return 1;
    } else {
        return 0;
    }
}

bool stCaf_filterByOutgroup(stPinchSegment *segment1,
                            stPinchSegment *segment2) {
    stPinchBlock *block1, *block2;
    if ((block1 = stPinchSegment_getBlock(segment1)) != NULL) {
        if ((block2 = stPinchSegment_getBlock(segment2)) != NULL) {
            if (block1 == block2) {
                return stPinchBlock_getLength(block1) == 1 ? 0 : containsOutgroupSegment(block1, flower);
            }
            if (stPinchBlock_getDegree(block1) < stPinchBlock_getDegree(block2)) {
                return containsOutgroupSegment(block1, flower) && containsOutgroupSegment(block2, flower);
            }
            return containsOutgroupSegment(block2, flower) && containsOutgroupSegment(block1, flower);
        }
        return isOutgroupSegment(segment2, flower) && containsOutgroupSegment(block1, flower);
    }
    if ((block2 = stPinchSegment_getBlock(segment2)) != NULL) {
        return isOutgroupSegment(segment1, flower) && containsOutgroupSegment(block2, flower);
    }
    return isOutgroupSegment(segment1, flower) && isOutgroupSegment(segment2, flower);
}

bool stCaf_relaxedFilterByOutgroup(stPinchSegment *segment1,
                                   stPinchSegment *segment2) {
    stPinchBlock *block1, *block2;
    if ((block1 = stPinchSegment_getBlock(segment1)) != NULL) {
        if ((block2 = stPinchSegment_getBlock(segment2)) != NULL) {
            if (block1 == block2) {
                return stPinchBlock_getLength(block1) == 1 ? 0 : containsOutgroupSegment(block1, flower);
            }
            if (stPinchBlock_getDegree(block1) < stPinchBlock_getDegree(block2)) {
                return containsOutgroupSegment(block1, flower) && containsOutgroupSegment(block2, flower);
            }
            return containsOutgroupSegment(block2, flower) && containsOutgroupSegment(block1, flower);
        }
    }
    // If we get here, we are just adding a segment to a block, not
    // pinching two blocks together.
    return false;
}

/*
 * Filtering by presence of repeat species in block. This code is inefficient and does not scale.
 */

static bool checkIntersection(stSortedSet *names1, stSortedSet *names2) {
    stSortedSet *n12 = stSortedSet_getIntersection(names1, names2);
    bool b = stSortedSet_size(n12) > 0;
    stSortedSet_destruct(names1);
    stSortedSet_destruct(names2);
    stSortedSet_destruct(n12);
    return b;
}

static stSortedSet *getEvents(stPinchSegment *segment, Flower *flower) {
    stSortedSet *events = stSortedSet_construct();
    if (stPinchSegment_getBlock(segment) != NULL) {
        stPinchBlock *block = stPinchSegment_getBlock(segment);
        stPinchBlockIt it = stPinchBlock_getSegmentIterator(block);
        while ((segment = stPinchBlockIt_getNext(&it)) != NULL) {
            stSortedSet_insert(events, stCaf_getEvent(segment, flower));
        }
    } else {
        stSortedSet_insert(events, stCaf_getEvent(segment, flower));
    }
    return events;
}

bool stCaf_filterByRepeatSpecies(stPinchSegment *segment1,
                                 stPinchSegment *segment2) {
    return checkIntersection(getEvents(segment1, flower), getEvents(segment2, flower));
}

bool stCaf_relaxedFilterByRepeatSpecies(stPinchSegment *segment1,
                                        stPinchSegment *segment2) {
    return stPinchSegment_getBlock(segment1) != NULL
        && stPinchSegment_getBlock(segment2) != NULL
        && checkIntersection(getEvents(segment1, flower), getEvents(segment2, flower));
}

static stSortedSet *getChrNames(stPinchSegment *segment, Flower *flower) {
    stSortedSet *names = stSortedSet_construct();
    if (stPinchSegment_getBlock(segment) != NULL) {
        stPinchBlock *block = stPinchSegment_getBlock(segment);
        stPinchBlockIt it = stPinchBlock_getSegmentIterator(block);
        while ((segment = stPinchBlockIt_getNext(&it)) != NULL) {
            Cap *cap = flower_getCap(flower, stPinchSegment_getName(segment));
            Sequence *sequence = cap_getSequence(cap);
            stSortedSet_insert(names, (void *) sequence_getName(sequence));
        }
    } else {
        Cap *cap = flower_getCap(flower, stPinchSegment_getName(segment));
        Sequence *sequence = cap_getSequence(cap);
        stSortedSet_insert(names, (void *) sequence_getName(sequence));
    }
    return names;
}

bool stCaf_singleCopyChr(stPinchSegment *segment1,
                         stPinchSegment *segment2) {
    return checkIntersection(getChrNames(segment1, flower), getChrNames(segment2, flower));
}

static stSortedSet *getIngroupEvents(stPinchSegment *segment, Flower *flower) {
    stSortedSet *events = stSortedSet_construct();
    if (stPinchSegment_getBlock(segment) != NULL) {
        stPinchBlock *block = stPinchSegment_getBlock(segment);
        stPinchBlockIt it = stPinchBlock_getSegmentIterator(block);
        while ((segment = stPinchBlockIt_getNext(&it)) != NULL) {
            Event *event = stCaf_getEvent(segment, flower);
            if (!event_isOutgroup(event)) {
                stSortedSet_insert(events, event);
            }
        }
    } else {
        Event *event = stCaf_getEvent(segment, flower);
        if (!event_isOutgroup(event)) {
            stSortedSet_insert(events, event);
        }
    }
    return events;
}

bool stCaf_singleCopyIngroup(stPinchSegment *segment1,
                             stPinchSegment *segment2) {
    return checkIntersection(getIngroupEvents(segment1, flower), getIngroupEvents(segment2, flower));
}

bool stCaf_relaxedSingleCopyIngroup(stPinchSegment *segment1,
                                    stPinchSegment *segment2) {
    return stPinchSegment_getBlock(segment1) != NULL
        && stPinchSegment_getBlock(segment2) != NULL
        && checkIntersection(getIngroupEvents(segment1, flower), getIngroupEvents(segment2, flower));
}

/*
 * Special filtering for a draft HGVM, where every chromosome should
 * be in its own component and should have no within-component cycles.
 */

static stHash *threadToComponent;

void stCaf_setThreadsToBeCycleFreeIsolatedComponents(stSet *input) {
    threadToComponent = stHash_construct();
    stSetIterator *it = stSet_getIterator(input);
    stPinchThread *thread;
    while ((thread = stSet_getNext(it)) != NULL) {
        stHash_insert(threadToComponent, thread, thread);
    }
}

// Assign a segment's thread (and all threads in its block) to a component.
static stPinchThread *assignComponent(stPinchSegment *segment, stPinchThread *component) {
    stPinchBlock *block = stPinchSegment_getBlock(segment);
    if (block != NULL) {
        stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
        while ((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
            stPinchThread *thread = stPinchSegment_getThread(segment);
            assert(stHash_search(threadToComponent, thread) == NULL);
            stHash_insert(threadToComponent, thread, component);
        }
    } else {
        stPinchThread *thread = stPinchSegment_getThread(segment);
        assert(stHash_search(threadToComponent, thread) == NULL);
        stHash_insert(threadToComponent, thread, component);
    }
    return false;
}

bool stCaf_filterToEnsureCycleFreeIsolatedComponents(stPinchSegment *segment1,
                                                     stPinchSegment *segment2) {
    stPinchThread *component1 = stHash_search(threadToComponent, stPinchSegment_getThread(segment1));
    stPinchThread *component2 = stHash_search(threadToComponent, stPinchSegment_getThread(segment2));

    // Check if this alignment will bridge two special components, or
    // create a cycle within one. If so, reject it.
    bool bridgesTwoComponents = false;
    if (component1 != NULL
        && component2 != NULL
        && component1 != component2) {
        // We reject this alignment.
        bridgesTwoComponents = true;
    } else {
        // This alignment will be applied, so update the components.
        if (component1 != component2 && component1 != NULL) {
            assignComponent(segment2, component1);
        } else if (component1 != component2 && component2 != NULL) {
            assignComponent(segment1, component2);
        }
    }
    return bridgesTwoComponents;
}

/*
 * Functions used for filtering blocks/chains on certain criteria.
 */

bool stCaf_chainHasUnequalNumberOfIngroupCopies(stCactusEdgeEnd *chainEnd,
                                                Flower *flower) {
    stPinchEnd *end = stCactusEdgeEnd_getObject(chainEnd);
    stPinchBlockIt it = stPinchBlock_getSegmentIterator(end->block);
    stPinchSegment *segment;
    stHash *ingroupToNumCopies = stHash_construct2(NULL, free);
    while ((segment = stPinchBlockIt_getNext(&it)) != NULL) {
        Cap *cap = flower_getCap(flower, stPinchSegment_getName(segment));
        Event *event = cap_getEvent(cap);
        if (!event_isOutgroup(event)) {
            if (stHash_search(ingroupToNumCopies, event) == NULL) {
                stHash_insert(ingroupToNumCopies, event, calloc(1, sizeof(uint64_t)));
            }
            uint64_t *numCopies = stHash_search(ingroupToNumCopies, event);
            (*numCopies)++;
        }
    }

    EventTree *eventTree = flower_getEventTree(flower);
    EventTree_Iterator *eventIt = eventTree_getIterator(eventTree);
    bool equalNumIngroupCopies = true;
    Event *event;
    uint64_t prevCount = 0;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        if (event_isOutgroup(event) || event_getChildNumber(event) != 0) {
            continue;
        }
        uint64_t *count = stHash_search(ingroupToNumCopies, event);
        if (count == NULL) {
            equalNumIngroupCopies = false;
            break;
        }
        if (prevCount == 0) {
            prevCount = *count;
        } else if (prevCount != *count) {
            equalNumIngroupCopies = false;
            break;
        }
    }

    eventTree_destructIterator(eventIt);
    stHash_destruct(ingroupToNumCopies);
    return !equalNumIngroupCopies;
}

bool stCaf_chainHasUnequalNumberOfIngroupCopiesOrNoOutgroup(stCactusEdgeEnd *chainEnd,
                                                            Flower *flower) {
    EventTree *eventTree = flower_getEventTree(flower);
    EventTree_Iterator *eventIt = eventTree_getIterator(eventTree);
    bool equalNumIngroupCopies = !stCaf_chainHasUnequalNumberOfIngroupCopies(chainEnd, flower);
    Event *event;
    uint64_t numOutgroups = 0;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        if (event_isOutgroup(event)) {
            numOutgroups++;
        }
    }

    uint64_t numOutgroupCopies = 0;
    stPinchEnd *end = stCactusEdgeEnd_getObject(chainEnd);
    stPinchBlockIt it = stPinchBlock_getSegmentIterator(end->block);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&it)) != NULL) {
        Cap *cap = flower_getCap(flower, stPinchSegment_getName(segment));
        Event *event = cap_getEvent(cap);
        if (event_isOutgroup(event)) {
            numOutgroupCopies++;
        }
    }

    eventTree_destructIterator(eventIt);
    return !equalNumIngroupCopies
        || (numOutgroups > 0 && numOutgroupCopies == 0);
}

bool stCaf_containsRequiredSpecies(stPinchBlock *pinchBlock,
                                   Flower *flower,
                                   int64_t minimumIngroupDegree,
                                   int64_t minimumOutgroupDegree,
                                   int64_t minimumDegree,
                                   int64_t minimumNumberOfSpecies) {
    stSet *seenEvents = stSet_construct();
    int64_t numberOfSpecies = 0;
    int64_t outgroupSequences = 0;
    int64_t ingroupSequences = 0;
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        Event *event = stCaf_getEvent(segment, flower);
        if (!stSet_search(seenEvents, event)) {
            stSet_insert(seenEvents, event);
            numberOfSpecies++;
        }
        if (event_isOutgroup(event)) {
            outgroupSequences++;
        } else {
            ingroupSequences++;
        }
    }
    stSet_destruct(seenEvents);
    return ingroupSequences >= minimumIngroupDegree &&
        outgroupSequences >= minimumOutgroupDegree &&
        outgroupSequences + ingroupSequences >= minimumDegree &&
        numberOfSpecies >= minimumNumberOfSpecies;
}

bool stCaf_treeCoverage(stPinchBlock *pinchBlock, Flower *flower) {
    EventTree *eventTree = flower_getEventTree(flower);
    Event *commonAncestorEvent = NULL;
    stPinchSegment *segment;
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    while ((segment = stPinchBlockIt_getNext(&segmentIt))) {
        Event *event = stCaf_getEvent(segment, flower);
        commonAncestorEvent = commonAncestorEvent == NULL ? event : eventTree_getCommonAncestor(event, commonAncestorEvent);
    }
    assert(commonAncestorEvent != NULL);
    float treeCoverage = 0.0;
    stHash *hash = stHash_construct();

    segmentIt = stPinchBlock_getSegmentIterator(pinchBlock);
    while ((segment = stPinchBlockIt_getNext(&segmentIt))) {
        Event *event = stCaf_getEvent(segment, flower);
        while (event != commonAncestorEvent && stHash_search(hash, event) == NULL) {
            treeCoverage += event_getBranchLength(event);
            stHash_insert(hash, event, event);
            event = event_getParent(event);
        }
    }

    float wholeTreeCoverage = event_getSubTreeBranchLength(event_getChild(eventTree_getRootEvent(eventTree), 0));
    assert(wholeTreeCoverage >= 0.0);
    if (wholeTreeCoverage <= 0.0) { //deal with case all leaf branches are not empty.
        return 0.0;
    }
    treeCoverage /= wholeTreeCoverage;
    assert(treeCoverage >= -0.001);
    assert(treeCoverage <= 1.0001);
    return treeCoverage;
}
