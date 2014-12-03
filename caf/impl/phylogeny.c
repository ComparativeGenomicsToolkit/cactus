/*
 * phlogeny.c
 *
 *  Created on: Jun 2, 2014
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "cactus.h"
#include "stPinchGraphs.h"
#include "stCactusGraphs.h"
#include "stPinchPhylogeny.h"
#include "stCaf.h"
#include "stCafPhylogeny.h"

// Doesn't have to be exported since nothing outside this file should
// really care about split branches.
typedef struct {
    stTree *child; // Child of the branch.
    stPinchBlock *block; // Block the tree refers to (can and should
                         // refer to more than the child subtree).
    double support; // Bootstrap support for this branch.
} stCaf_SplitBranch;

// gross, but useful for debug prints of what the scale of the
// single-degree segment problem is
// FIXME: get rid of these once done
int64_t numSingleDegreeSegmentsDropped = 0;
int64_t numBasesDroppedFromSingleDegreeSegments = 0;

int64_t numSimpleBlocksSkipped = 0;
int64_t numSingleCopyBlocksSkipped = 0;

stHash *stCaf_getThreadStrings(Flower *flower, stPinchThreadSet *threadSet) {
    stHash *threadStrings = stHash_construct2(NULL, free);
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *thread;
    while((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        assert(cap != NULL);
        Sequence *sequence = cap_getSequence(cap);
        assert(sequence != NULL);
        assert(stPinchThread_getLength(thread)-2 >= 0);
        char *string = sequence_getString(sequence, stPinchThread_getStart(thread)+1, stPinchThread_getLength(thread)-2, 1); //Gets the sequence excluding the empty positions representing the caps.
        char *paddedString = stString_print("N%sN", string); //Add in positions to represent the flanking bases
        stHash_insert(threadStrings, thread, paddedString);
        free(string);
    }
    return threadStrings;
}

stSet *stCaf_getOutgroupThreads(Flower *flower, stPinchThreadSet *threadSet) {
    stSet *outgroupThreads = stSet_construct();
    stPinchThreadSetIt threadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *thread;
    while ((thread = stPinchThreadSetIt_getNext(&threadIt)) != NULL) {
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        assert(cap != NULL);
        Event *event = cap_getEvent(cap);
        if(event_isOutgroup(event)) {
            stSet_insert(outgroupThreads, thread);
        }
    }
    return outgroupThreads;
}


/*
 * Gets a list of the segments in the block that are part of outgroup threads.
 * The list contains stIntTuples, each of length 1, representing the index of a particular segment in
 * the block.
 */
static stList *getOutgroupThreads(stPinchBlock *block, stSet *outgroupThreads) {
    stList *outgroups = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    int64_t i=0;
    while((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        if(stSet_search(outgroupThreads, stPinchSegment_getThread(segment)) != NULL) {
            stList_append(outgroups, stIntTuple_construct1(i));
        }
        i++;
    }
    assert(i == stPinchBlock_getDegree(block));
    return outgroups;
}

/*
 * Splits the block using the given partition into a set of new blocks.
 */
void splitBlock(stPinchBlock *block, stList *partitions, bool allowSingleDegreeBlocks) {
    assert(stList_length(partitions) > 0);
    if(stList_length(partitions) == 1) {
        return; //Nothing to do.
    }
    //Build a mapping of indices of the segments in the block to the segments
    int64_t blockDegree = stPinchBlock_getDegree(block);
    stPinchSegment **segments = st_calloc(blockDegree, sizeof(stPinchSegment *));
    bool *orientations = st_calloc(blockDegree, sizeof(bool));
    stPinchBlockIt segmentIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    int64_t i=0;
    while((segment = stPinchBlockIt_getNext(&segmentIt)) != NULL) {
        segments[i] = segment;
        assert(segments[i] != NULL);
        orientations[i++] = stPinchSegment_getBlockOrientation(segment);
    }
    assert(i == stPinchBlock_getDegree(block));
    //Destruct old block, as we build new blocks now.
    stPinchBlock_destruct(block);
    //Now build the new blocks.
    for(int64_t i=0; i<stList_length(partitions); i++) {
        stList *partition = stList_get(partitions, i);
        assert(stList_length(partition) > 0);
        int64_t k = stIntTuple_get(stList_get(partition, 0), 0);
        assert(segments[k] != NULL);
        assert(stPinchSegment_getBlock(segments[k]) == NULL);

        if (!allowSingleDegreeBlocks && stList_length(partition) == 1) {
            // We need to avoid assigning this single-degree block
            numBasesDroppedFromSingleDegreeSegments += stPinchSegment_getLength(segments[k]);
            numSingleDegreeSegmentsDropped++;
            segments[k] = NULL;
            continue;
        }

        block = stPinchBlock_construct3(segments[k], orientations[k]);
        assert(stPinchSegment_getBlock(segments[k]) == block);
        assert(stPinchSegment_getBlockOrientation(segments[k]) == orientations[k]);
        segments[k] = NULL; //Defensive, and used for debugging.
        for(int64_t j=1; j<stList_length(partition); j++) {
            k = stIntTuple_get(stList_get(partition, j), 0);
            assert(segments[k] != NULL);
            assert(stPinchSegment_getBlock(segments[k]) == NULL);
            stPinchBlock_pinch2(block, segments[k], orientations[k]);
            assert(stPinchSegment_getBlock(segments[k]) == block);
            assert(stPinchSegment_getBlockOrientation(segments[k]) == orientations[k]);
            segments[k] = NULL; //Defensive, and used for debugging.
        }
    }
    //Now check the segments have all been used - this is just debugging.
    for(int64_t i=0; i<blockDegree; i++) {
        assert(segments[i] == NULL);
    }
    //Cleanup
    free(segments);
    free(orientations);
}

/*
 * For logging purposes gets the total number of similarities and differences in the matrix.
 */
static void getTotalSimilarityAndDifferenceCounts(stMatrix *matrix, double *similarities, double *differences) {
    *similarities = 0.0;
    *differences = 0.0;
    for(int64_t i=0; i<stMatrix_n(matrix); i++) {
        for(int64_t j=i+1; j<stMatrix_n(matrix); j++) {
            *similarities += *stMatrix_getCell(matrix, i, j);
            *differences += *stMatrix_getCell(matrix, j, i);
        }
    }
}

// If the tree contains any zero branch lengths (i.e. there were
// negative branch lengths when neighbor-joining), fudge the branch
// lengths so that both children have non-zero branch lengths, but are
// still the same distance apart. When both children have zero branch
// lengths, give them both a small branch length. This makes
// likelihood methods usable.

// Only works on binary trees.
static void fudgeZeroBranchLengths(stTree *tree, double fudgeFactor, double smallNonZeroBranchLength) {
    assert(stTree_getChildNumber(tree) == 2 || stTree_getChildNumber(tree) == 0);
    assert(fudgeFactor < 1.0 && fudgeFactor > 0.0);
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree *child = stTree_getChild(tree, i);
        fudgeZeroBranchLengths(child, fudgeFactor, smallNonZeroBranchLength);
        if (stTree_getBranchLength(child) == 0.0) {
            stTree *otherChild = stTree_getChild(tree, !i);
            if (stTree_getBranchLength(otherChild) == 0.0) {
                // Both children have zero branch lengths, set them
                // both to some very small but non-zero branch length
                // so that probabilistic methods can actually work
                stTree_setBranchLength(child, smallNonZeroBranchLength);
                stTree_setBranchLength(otherChild, smallNonZeroBranchLength);
            } else {
                // Keep the distance between the children equal, but
                // move it by fudgeFactor so that no branch length is
                // zero.
                stTree_setBranchLength(child, fudgeFactor * stTree_getBranchLength(otherChild));
                stTree_setBranchLength(otherChild, (1 - fudgeFactor) * stTree_getBranchLength(otherChild));
            }
        }
    }
}
/*
 * Get an Event pointer -> species tree node mapping.
 */
static stHash *getEventToSpeciesNode(EventTree *eventTree,
                                     stTree *speciesTree) {
    stHash *ret = stHash_construct();
    EventTree_Iterator *eventIt = eventTree_getIterator(eventTree);
    Event *event;
    while ((event = eventTree_getNext(eventIt)) != NULL) {
        char *speciesLabel = stString_print("%" PRIi64, event_getName(event));
        stTree *species = stTree_findChild(speciesTree, speciesLabel);
        if (species != NULL) {
            stHash_insert(ret, event, species);
        } else {
            // Every node should be in the species tree besides the root event.
            assert(event == eventTree_getRootEvent(eventTree));
        }
        free(speciesLabel);
    }
    eventTree_destructIterator(eventIt);
    return ret;
}

/*
 * Get a gene node->species node mapping from a gene tree, a species
 * tree, and the pinch block.
 */

static stHash *getLeafToSpecies(stTree *geneTree, stPinchBlock *block,
                                Flower *flower, stHash *eventToSpeciesNode) {
    stHash *leafToSpecies = stHash_construct();
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    int64_t i = 0; // Current segment index in block.
    while((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        stPinchThread *thread = stPinchSegment_getThread(segment);
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        Event *event = cap_getEvent(cap);
        stTree *species = stHash_search(eventToSpeciesNode, event);
        assert(species != NULL);
        stTree *gene = stPhylogeny_getLeafByIndex(geneTree, i);
        assert(gene != NULL);
        stHash_insert(leafToSpecies, gene, species);
        i++;
    }
    return leafToSpecies;
}

/*
 * Get a mapping from matrix index -> join cost index for use in
 * neighbor-joining guided by a species tree.
 */

static stHash *getMatrixIndexToJoinCostIndex(stPinchBlock *block, Flower *flower, stHash* eventToSpeciesNode, stHash *speciesToJoinCostIndex) {
    stHash *matrixIndexToJoinCostIndex = stHash_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey, (int (*)(const void *, const void *)) stIntTuple_equalsFn, (void (*)(void *)) stIntTuple_destruct, (void (*)(void *)) stIntTuple_destruct);
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    int64_t i = 0; // Current segment index in block.
    while((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        stPinchThread *thread = stPinchSegment_getThread(segment);
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        Event *event = cap_getEvent(cap);
        stTree *species = stHash_search(eventToSpeciesNode, event);
        assert(species != NULL);

        stIntTuple *joinCostIndex = stHash_search(speciesToJoinCostIndex, species);
        assert(joinCostIndex != NULL);

        stIntTuple *matrixIndex = stIntTuple_construct1(i);
        stHash_insert(matrixIndexToJoinCostIndex, matrixIndex,
                      // Copy the join cost index so it has the same
                      // lifetime as the hash
                      stIntTuple_construct1(stIntTuple_get(joinCostIndex, 0)));
        i++;
    }
    return matrixIndexToJoinCostIndex;
}

static double scoreTree(stTree *tree, enum stCaf_ScoringMethod scoringMethod, stTree *speciesStTree, stPinchBlock *block, Flower *flower, stList *featureColumns, stHash *eventToSpeciesNode) {
    double ret = 0.0;
    if (scoringMethod == RECON_COST) {
        stHash *leafToSpecies = getLeafToSpecies(tree, block, flower,
                                                 eventToSpeciesNode);
        int64_t dups = 0, losses = 0;
        stPhylogeny_reconciliationCostAtMostBinary(tree, &dups, &losses);
        ret = -dups - losses;

        stHash_destruct(leafToSpecies);
    } else if (scoringMethod == NUCLEOTIDE_LIKELIHOOD) {
        ret = stPinchPhylogeny_likelihood(tree, featureColumns);
    } else if (scoringMethod == RECON_LIKELIHOOD) {
        // FIXME: hardcoding dup-rate parameter for now
        ret = stPinchPhylogeny_reconciliationLikelihood(tree, speciesStTree, 1.0);
    } else if (scoringMethod == COMBINED_LIKELIHOOD) {
        // FIXME: hardcoding dup-rate parameter for now
        ret = stPinchPhylogeny_reconciliationLikelihood(tree, speciesStTree, 1.0);
        ret += stPinchPhylogeny_likelihood(tree, featureColumns);
    }
    return ret;
}

// Build a tree from a set of feature columns and root it according to
// the rooting method.
static stTree *buildTree(stList *featureColumns,
                         stCaf_PhylogenyParameters *params,
                         bool bootstrap,
                         stList *outgroups, stPinchBlock *block,
                         Flower *flower, stTree *speciesStTree,
                         stMatrix *joinCosts,
                         stHash *speciesToJoinCostIndex,
                         int64_t **speciesMRCAMatrix,
                         stHash *eventToSpeciesNode,
                         stMatrixDiffs *snpDiffs,
                         stMatrixDiffs *breakpointDiffs) {
    // Make substitution matrix
    stMatrix *substitutionMatrix = stPinchPhylogeny_constructMatrixFromDiffs(snpDiffs, bootstrap);
    assert(stMatrix_n(substitutionMatrix) == stPinchBlock_getDegree(block));
    assert(stMatrix_m(substitutionMatrix) == stPinchBlock_getDegree(block));
    //Make breakpoint matrix
    stMatrix *breakpointMatrix = stPinchPhylogeny_constructMatrixFromDiffs(breakpointDiffs, bootstrap);
    
    //Combine the matrices into distance matrices
    stMatrix_scale(breakpointMatrix, params->breakpointScalingFactor, 0.0);
    stMatrix *combinedMatrix = stMatrix_add(substitutionMatrix, breakpointMatrix);
    stMatrix *distanceMatrix = stPinchPhylogeny_getSymmetricDistanceMatrix(combinedMatrix);

    stTree *tree = NULL;
    if (params->rootingMethod == OUTGROUP_BRANCH) {
        if (params->treeBuildingMethod == NEIGHBOR_JOINING) {
            tree = stPhylogeny_neighborJoin(distanceMatrix, outgroups);
        } else {
            assert(params->treeBuildingMethod == GUIDED_NEIGHBOR_JOINING);
            st_errAbort("Longest-outgroup-branch rooting not supported with guided neighbor joining");
        }
    } else if (params->rootingMethod == LONGEST_BRANCH) {
        if (params->treeBuildingMethod == NEIGHBOR_JOINING) {
            tree = stPhylogeny_neighborJoin(distanceMatrix, NULL);
        } else {
            assert(params->treeBuildingMethod == GUIDED_NEIGHBOR_JOINING);
            st_errAbort("Longest-branch rooting not supported with guided neighbor joining");
        }
    } else if (params->rootingMethod == BEST_RECON) {
        if (params->treeBuildingMethod == NEIGHBOR_JOINING) {
            tree = stPhylogeny_neighborJoin(distanceMatrix, NULL);
        } else {
            // FIXME: Could move this out of the function as
            // well. It's the same for each tree generated for the
            // block.
            stHash *matrixIndexToJoinCostIndex = getMatrixIndexToJoinCostIndex(block, flower, eventToSpeciesNode,
                                                                               speciesToJoinCostIndex);
            tree = stPhylogeny_guidedNeighborJoining(combinedMatrix, joinCosts, matrixIndexToJoinCostIndex, speciesToJoinCostIndex, speciesMRCAMatrix, speciesStTree);
            stHash_destruct(matrixIndexToJoinCostIndex);
        }
        stHash *leafToSpecies = getLeafToSpecies(tree, block, flower,
                                                 eventToSpeciesNode);

        stTree *newTree = stPhylogeny_rootByReconciliationAtMostBinary(tree, leafToSpecies);
        stPhylogeny_addStIndexedTreeInfo(newTree);

        stPhylogenyInfo_destructOnTree(tree);
        stTree_destruct(tree);
        stHash_destruct(leafToSpecies);
        tree = newTree;
    }

    // Need to get leaf->species mapping again even when the tree is
    // already reconciled, as the nodes have changed.
    stHash *leafToSpecies = getLeafToSpecies(tree, block, flower,
                                             eventToSpeciesNode);
    // add stReconciliationInfo.
    stPhylogeny_reconcileAtMostBinary(tree, leafToSpecies, false);
    stHash_destruct(leafToSpecies);

    // Needed for likelihood methods not to have 0/100% probabilities
    // overly often (normally, almost every other leaf has a branch
    // length of 0)
    fudgeZeroBranchLengths(tree, 0.02, 0.0001);

    stMatrix_destruct(substitutionMatrix);
    stMatrix_destruct(breakpointMatrix);
    stMatrix_destruct(combinedMatrix);
    stMatrix_destruct(distanceMatrix);
    return tree;
}

// Check if the block's phylogeny is simple:
// - the block has only one event, or
// - the block has < 3 segments.
static bool hasSimplePhylogeny(stPinchBlock *block,
                               Flower *flower) {
    if(stPinchBlock_getDegree(block) <= 2) {
        return true;
    }
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment = NULL;
    bool found2Events = 0;
    Event *currentEvent = NULL;
    while((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        stPinchThread *thread = stPinchSegment_getThread(segment);
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        assert(cap != NULL);
        Event *event = cap_getEvent(cap);
        if(currentEvent == NULL) {
            currentEvent = event;
        } else if(currentEvent != event) {
            found2Events = 1;
        }
    }
    return !found2Events;
}

// Check if the block contains as many species as segments
static bool isSingleCopyBlock(stPinchBlock *block, Flower *flower) {
    stSet *seenEvents = stSet_construct();
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment = NULL;
    while((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        stPinchThread *thread = stPinchSegment_getThread(segment);
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        assert(cap != NULL);
        Event *event = cap_getEvent(cap);
        if(stSet_search(seenEvents, event) != NULL) {
            stSet_destruct(seenEvents);
            return false;
        }
        stSet_insert(seenEvents, event);
    }
    stSet_destruct(seenEvents);
    return true;
}

// relabel a tree so it's useful for debug output
static void relabelMatrixIndexedTree(stTree *tree, stHash *matrixIndexToName) {
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        relabelMatrixIndexedTree(stTree_getChild(tree, i), matrixIndexToName);
    }
    if (stTree_getChildNumber(tree) == 0) {
        stPhylogenyInfo *info = stTree_getClientData(tree);
        assert(info != NULL);
        assert(info->index->matrixIndex != -1);
        stIntTuple *query = stIntTuple_construct1(info->index->matrixIndex);
        char *header = stHash_search(matrixIndexToName, query);
        assert(header != NULL);
        stTree_setLabel(tree, stString_copy(header));
        stIntTuple_destruct(query);
    }
}

// Print the debug info for blocks that are normally not printed
// (those that cannot be partitioned). The debug info contains just
// the "partition", i.e. the sequences and positions within the block
// in a list of lists.
static void printSimpleBlockDebugInfo(Flower *flower, stPinchBlock *block, FILE *outFile) {
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment = NULL;
    fprintf(outFile, "[[");
    int64_t i = 0;
    while ((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        stPinchThread *thread = stPinchSegment_getThread(segment);
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        const char *seqHeader = sequence_getHeader(cap_getSequence(cap));
        Event *event = cap_getEvent(cap);
        const char *eventHeader = event_getHeader(event);
        char *segmentHeader = stString_print("%s.%s|%" PRIi64 "-%" PRIi64, eventHeader, seqHeader, stPinchSegment_getStart(segment), stPinchSegment_getStart(segment) + stPinchSegment_getLength(segment));
        
        if (i != 0) {
            fprintf(outFile, ",");
        }
        fprintf(outFile, "\"%s\"", segmentHeader);
        free(segmentHeader);
        i++;
    }
    assert(i == stPinchBlock_getDegree(block));
    fprintf(outFile, "]]\n");
}

// print debug info: "tree\tpartition\n" to the file
static void printTreeBuildingDebugInfo(Flower *flower, stPinchBlock *block, stTree *bestTree, stList *partition, stMatrix *matrix, double score, FILE *outFile) {
    // First get a map from matrix indices to names
    // The format we will use for leaf names is "genome.seq|posStart-posEnd"
    int64_t blockDegree = stPinchBlock_getDegree(block);
    stHash *matrixIndexToName = stHash_construct3((uint64_t (*)(const void *)) stIntTuple_hashKey,
                                                  (int (*)(const void *, const void *)) stIntTuple_equalsFn,
                                                  (void (*)(void *)) stIntTuple_destruct, free);
    int64_t i = 0;
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment = NULL;
    while ((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        stPinchThread *thread = stPinchSegment_getThread(segment);
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        const char *seqHeader = sequence_getHeader(cap_getSequence(cap));
        Event *event = cap_getEvent(cap);
        const char *eventHeader = event_getHeader(event);
        char *segmentHeader = stString_print("%s.%s|%" PRIi64 "-%" PRIi64, eventHeader, seqHeader, stPinchSegment_getStart(segment), stPinchSegment_getStart(segment) + stPinchSegment_getLength(segment));
        stHash_insert(matrixIndexToName, stIntTuple_construct1(i), segmentHeader);
        i++;
    }
    assert(i == blockDegree);

    // Relabel (our copy of) the best tree.
    stTree *treeCopy = stTree_clone(bestTree);
    relabelMatrixIndexedTree(treeCopy, matrixIndexToName);
    char *newick = stTree_getNewickTreeString(treeCopy);

    fprintf(outFile, "%s\t", newick);

    // Print the partition
    fprintf(outFile, "[");
    for (i = 0; i < stList_length(partition); i++) {
        if (i != 0) {
            fprintf(outFile, ",");
        }
        stList *subList = stList_get(partition, i);
        fprintf(outFile, "[");
        for (int64_t j = 0; j < stList_length(subList); j++) {
            if (j != 0) {
                fprintf(outFile, ",");
            }
            stIntTuple *index = stList_get(subList, j);
            assert(stIntTuple_get(index, 0) < blockDegree);
            char *header = stHash_search(matrixIndexToName, index);
            assert(header != NULL);
            fprintf(outFile, "\"%s\"", header);
        }
        fprintf(outFile, "]");
    }
    fprintf(outFile, "]\t");

    if (matrix != NULL) {
        // print the matrix
        fprintf(outFile, "[");
        for (i = 0; i < stMatrix_m(matrix); i++) {
            if (i != 0) {
                fprintf(outFile, ",");
            }
            fprintf(outFile, "[");
            for (int64_t j = 0; j < stMatrix_n(matrix); j++) {
                if (j != 0) {
                    fprintf(outFile, ",");
                }
                fprintf(outFile, "%lf", *stMatrix_getCell(matrix, i, j));
            }
            fprintf(outFile, "]");
        }
        fprintf(outFile, "]\t");
    }
    // print the sequences corresponding to the matrix indices
    fprintf(outFile, "[");
    for (int64_t i = 0; i < blockDegree; i++) {
        if (i != 0) {
            fprintf(outFile, ",");
        }
        stIntTuple *query = stIntTuple_construct1(i);
        char *header = stHash_search(matrixIndexToName, query);
        assert(header != NULL);
        fprintf(outFile, "\"%s\"", header);
        stIntTuple_destruct(query);
    }
    fprintf(outFile, "]\t");

    // print the score
    fprintf(outFile, "%lf\n", score);
    stTree_destruct(treeCopy);
    free(newick);
    stHash_destruct(matrixIndexToName);
}

static int64_t countBasesBetweenSingleDegreeBlocks(stPinchThreadSet *threadSet) {
    stPinchThreadSetIt pinchThreadIt = stPinchThreadSet_getIt(threadSet);
    stPinchThread *thread;
    int64_t numBases = 0;
    int64_t numBasesInSingleCopyBlocks = 0;
    while ((thread = stPinchThreadSetIt_getNext(&pinchThreadIt)) != NULL) {
        stPinchSegment *segment = stPinchThread_getFirst(thread);
        if (segment == NULL) {
            // No segments on this thread.
            continue;
        }
        bool wasInSingleDegreeBlock = stPinchBlock_getDegree(stPinchSegment_getBlock(segment)) == 1;
        stPinchSegment *oldSegment = NULL;
        while ((segment = stPinchSegment_get3Prime(segment)) != NULL) {
            stPinchBlock *block = stPinchSegment_getBlock(segment);
            if (block == NULL) {
                // Segment without a block.
                continue;
            }
            bool isInSingleDegreeBlock = stPinchBlock_getDegree(block) == 1;
            if (isInSingleDegreeBlock) {
                numBasesInSingleCopyBlocks += stPinchBlock_getLength(block);
            }
            int64_t numBasesBetweenSegments = 0;
            if (oldSegment != NULL) {
                numBasesBetweenSegments = stPinchSegment_getStart(segment) - (stPinchSegment_getStart(oldSegment) + stPinchSegment_getLength(oldSegment));
            }
            assert(numBasesBetweenSegments >= 0); // could be 0 if the
                                                  // blocks aren't
                                                  // identical
            if (wasInSingleDegreeBlock && isInSingleDegreeBlock) {
                numBases += numBasesBetweenSegments;
            }
            oldSegment = segment;
            wasInSingleDegreeBlock = isInSingleDegreeBlock;
        }
    }
    // FIXME: tmp
    fprintf(stdout, "There were %" PRIi64 " bases in single degree blocks.\n", numBasesInSingleCopyBlocks);
    return numBases;
}

// Compare two bootstrap scores of split branches. Use the pointer
// value of the branches as a tiebreaker since we are using a sorted
// set and don't want to merge together all branches with the same
// support value.
int compareSplitBranches(stCaf_SplitBranch *branch1,
                         stCaf_SplitBranch *branch2) {
    if (branch1->support > branch2->support) {
        return 2;
    } else if (branch1->support == branch2->support) {
        if (branch1->child == branch2->child) {
            return 0;
        } else if  (branch1->child > branch2->child) {
            return 1;
        } else {
            return -1;
        }
    } else {
        return -2;
    }
}

stCaf_SplitBranch *stCaf_SplitBranch_construct(stTree *child,
                                               stPinchBlock *block,
                                               double support) {
    stCaf_SplitBranch *ret = calloc(1, sizeof(stCaf_SplitBranch));
    ret->support = support;
    ret->child = child;
    ret->block = block;
    return ret;
}

// Find new split branches from the block and add them to the sorted set.
// speciesToSplitOn is just the species that are on the path from the
// reference node to the root.
void findSplitBranches(stPinchBlock *block, stTree *tree,
                       stSortedSet *splitBranches,
                       stSet *speciesToSplitOn) {
    stTree *parent = stTree_getParent(tree);
    if (parent != NULL) {
        stPhylogenyInfo *parentInfo = stTree_getClientData(parent);
        assert(parentInfo != NULL);
        stReconciliationInfo *parentReconInfo = parentInfo->recon;
        assert(parentReconInfo != NULL);
        if (parentReconInfo->event == DUPLICATION && stSet_search(speciesToSplitOn, parentReconInfo->species)) {
            // Found a split branch.
            stPhylogenyInfo *info = stTree_getClientData(tree);
            stCaf_SplitBranch *splitBranch = stCaf_SplitBranch_construct(tree, block, info->index->bootstrapSupport);
            stSortedSet_insert(splitBranches, splitBranch);
        } else {
            // Since the reconciliation must follow the order of the
            // species tree, any child of this node cannot be
            // reconciled to the speciesToSplitOnSet.
            return;
        }
    }

    // Recurse down the tree as far as makes sense.
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        findSplitBranches(block, stTree_getChild(tree, i), splitBranches,
                          speciesToSplitOn);
    }
}

// Get species that are a candidate to split on if a duplication node
// is reconciled to them (i.e. everything at or above the reference
// event node).
static void getSpeciesToSplitOn(stTree *speciesTree, EventTree *eventTree,
                                const char *referenceEventHeader,
                                stSet *speciesToSplitOn) {
    if (referenceEventHeader == NULL) {
        st_errAbort("unrecoverable error: cannot find species to split on "
                    "without knowing reference event header.");
    }
    // Look through the entire species tree (which is labeled by event
    // Name as a string) and find the event whose header matches the
    // reference event.
    stTree *referenceSpecies = NULL;
    stList *bfQueue = stList_construct();
    stList_append(bfQueue, speciesTree);
    while (stList_length(bfQueue) > 0) {
        stTree *node = stList_pop(bfQueue);
        for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
            stList_append(bfQueue, stTree_getChild(node, i));
        }

        Name name;
        int ret = sscanf(stTree_getLabel(node), "%" PRIi64, &name);
        (void) ret;
        assert(ret == 1);
        Event *event = eventTree_getEvent(eventTree, name);
        const char *header = event_getHeader(event);
        if (header == NULL) {
            continue;
        }
        if (strcmp(header, referenceEventHeader) == 0) {
            referenceSpecies = node;
        }
    }
    stList_destruct(bfQueue);

    if (referenceSpecies == NULL) {
        st_errAbort("Could not find reference event header %s in species "
                    "tree.\n", referenceEventHeader);
    }

    stTree *aSpeciesToSplitOn = referenceSpecies;
    do {
        stSet_insert(speciesToSplitOn, aSpeciesToSplitOn);
        aSpeciesToSplitOn = stTree_getParent(aSpeciesToSplitOn);
    } while (aSpeciesToSplitOn != NULL);
}

// O(n).
static stPinchSegment *getSegmentByBlockIndex(stPinchBlock *block,
                                              int64_t index) {
    assert(index >= 0);
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment = NULL;
    int64_t i = 0;
    while((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        if (i == index) {
            return segment;
        }
        i++;
    }
    return NULL;
}

// Add any stPinchBlocks close enough to the given block to be
// affected by its breakpoint information to the given set.
static void addContextualBlocksToSet(stPinchBlock *block,
                                     int64_t maxBaseDistance,
                                     int64_t maxBlockDistance,
                                     bool ignoreUnalignedBases,
                                     stSet *contextualBlocks) {
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    while ((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        // Go toward the 5' end of this thread adding blocks, until we
        // reach the end or maxBaseDistance or maxBlockDistance.
        stPinchSegment *curSegment = stPinchSegment_get5Prime(segment);
        int64_t curBlockDistance = 0;
        int64_t curBaseDistance = stPinchSegment_getLength(segment) / 2;
        while ((curSegment != NULL) && (curBlockDistance < maxBlockDistance)
               && (curBaseDistance < maxBaseDistance)) {
            stPinchBlock *curBlock = stPinchSegment_getBlock(curSegment);
            if (curBlock != NULL) {
                stSet_insert(contextualBlocks, curBlock);
                curBaseDistance += stPinchSegment_getLength(segment);
                curBlockDistance++;
            } else if (!ignoreUnalignedBases) {
                curBaseDistance += stPinchSegment_getLength(segment);
                curBlockDistance++;
            }
            curSegment = stPinchSegment_get5Prime(curSegment);
        }

        // Do the same for the 3' side.
        curSegment = stPinchSegment_get3Prime(segment);
        curBlockDistance = 0;
        curBaseDistance = stPinchSegment_getLength(segment) / 2;
        while ((curSegment != NULL) && (curBlockDistance < maxBlockDistance)
               && (curBaseDistance < maxBaseDistance)) {
            stPinchBlock *curBlock = stPinchSegment_getBlock(curSegment);
            if (curBlock != NULL) {
                stSet_insert(contextualBlocks, curBlock);
                curBaseDistance += stPinchSegment_getLength(segment);
                curBlockDistance++;
            } else if (!ignoreUnalignedBases) {
                curBaseDistance += stPinchSegment_getLength(segment);
                curBlockDistance++;
            }
            curSegment = stPinchSegment_get3Prime(curSegment);
        }
    }
}

static void removeOldSplitBranches(stPinchBlock *block, stTree *tree,
                                   stSet *speciesToSplitOn,
                                   stSortedSet *splitBranches) {
    if (block == NULL || tree == NULL) {
        return;
    }
    stSortedSet *splitBranchesToDelete = stSortedSet_construct3((int (*)(const void *, const void *)) compareSplitBranches, free);
    findSplitBranches(block, tree, splitBranchesToDelete,
                      speciesToSplitOn);
    // Could set splitBranches = splitBranches \ splitBranchesToDelete.
    // But this is probably faster.
    stSortedSetIterator *splitBranchToDeleteIt = stSortedSet_getIterator(splitBranchesToDelete);
    stCaf_SplitBranch *splitBranchToDelete;
    while ((splitBranchToDelete = stSortedSet_getNext(splitBranchToDeleteIt)) != NULL) {
        stSortedSet_remove(splitBranches, splitBranchToDelete);
    }
    stSortedSet_destructIterator(splitBranchToDeleteIt);
    stSortedSet_destruct(splitBranchesToDelete);
}

// tmp for debug
static double getAvgSupportValue(stSortedSet *splitBranches) {
    double totalSupport = 0.0;
    stSortedSetIterator *splitBranchIt = stSortedSet_getIterator(splitBranches);
    stCaf_SplitBranch *splitBranch;
    while ((splitBranch = stSortedSet_getNext(splitBranchIt)) != NULL) {
        totalSupport += splitBranch->support;
    }
    stSortedSet_destructIterator(splitBranchIt);
    if (stSortedSet_size(splitBranches) == 0) {
        return 0.0;
    }
    return totalSupport/stSortedSet_size(splitBranches);
}

// Struct of constant things that gets passed around. Since these are
// only set once in a run, they could be global variables, but this is
// just in case we ever need to run in parallel on sub-flowers or
// something weird.
typedef struct {
    stHash *threadStrings;
    stSet *outgroupThreads;
    Flower *flower;
    stCaf_PhylogenyParameters *params;
    stMatrix *joinCosts;
    stHash *speciesToJoinCostIndex;
    int64_t **speciesMRCAMatrix;
    stHash *eventToSpeciesNode;
    stTree *speciesStTree;
} TreeBuildingConstants;

// Gets passed to buildTreeForBlock.
typedef struct {
    stPinchBlock *block;
    TreeBuildingConstants *constants;
    stHash *blocksToTrees;
} TreeBuildingInput;

// Gets returned from buildTreeForBlock and passed into
// addBlockTreeToHash.
typedef struct {
    stHash *blocksToTrees;
    stTree *tree;
    stPinchBlock *block;
    TreeBuildingConstants *constants; // Temp -- just so we can get
                                      // access to the flower for the
                                      // single-copy debug print.
} TreeBuildingResult;

// Small wrapper function to tell the pool to build, reconcile, and
// bootstrap a tree for a block.
static void pushBlockToPool(stPinchBlock *block,
                            TreeBuildingConstants *constants,
                            stHash *blocksToTrees,
                            stThreadPool *threadPool) {
    TreeBuildingInput *input = st_malloc(sizeof(TreeBuildingInput));
    input->constants = constants;
    input->block = block;
    input->blocksToTrees = blocksToTrees;
    stThreadPool_push(threadPool, input);
}

// Gets run as a worker in a thread.
static TreeBuildingResult *buildTreeForBlock(TreeBuildingInput *input) {
    stPinchBlock *block = input->block;
    stCaf_PhylogenyParameters *params = input->constants->params;

    TreeBuildingResult *ret = st_calloc(1, sizeof(TreeBuildingResult));
    ret->blocksToTrees = input->blocksToTrees;
    ret->block = block;
    ret->constants = input->constants;

    if (hasSimplePhylogeny(block, input->constants->flower)) {
        // No point trying to build a phylogeny for certain blocks.
        free(input);
        return ret;
    }
    if (isSingleCopyBlock(block, input->constants->flower)
        && params->skipSingleCopyBlocks) {
        free(input);
        return ret;
    }

    // Get the feature blocks
    stList *featureBlocks = stFeatureBlock_getContextualFeatureBlocks(
        block, params->maxBaseDistance,
        params->maxBlockDistance,
        params->ignoreUnalignedBases, params->onlyIncludeCompleteFeatureBlocks,
        input->constants->threadStrings);

    // Make feature columns
    stList *featureColumns = stFeatureColumn_getFeatureColumns(featureBlocks,
                                                               block);

    // Get the outgroup threads
    stList *outgroups = getOutgroupThreads(block,
                                           input->constants->outgroupThreads);

    // Get the matrix diffs.
    stMatrixDiffs *snpDiffs = stPinchPhylogeny_getMatrixDiffsFromSubstitutions(featureColumns, block, NULL);
    stMatrixDiffs *breakpointDiffs = stPinchPhylogeny_getMatrixDiffsFromBreakpoints(featureColumns, block, NULL);

    // Build the canonical tree.
    stTree *blockTree = buildTree(featureColumns, params,
                                  0, outgroups, block, input->constants->flower,
                                  input->constants->speciesStTree,
                                  input->constants->joinCosts,
                                  input->constants->speciesToJoinCostIndex,
                                  input->constants->speciesMRCAMatrix,
                                  input->constants->eventToSpeciesNode,
                                  snpDiffs, breakpointDiffs);

    // Sample the rest of the trees.
    stList *trees = stList_construct();
    stList_append(trees, blockTree);
    for (int64_t i = 0; i < params->numTrees - 1; i++) {
        stTree *tree = buildTree(featureColumns, params,
                                 1, outgroups, block, input->constants->flower,
                                 input->constants->speciesStTree,
                                 input->constants->joinCosts,
                                 input->constants->speciesToJoinCostIndex,
                                 input->constants->speciesMRCAMatrix,
                                 input->constants->eventToSpeciesNode,
                                 snpDiffs, breakpointDiffs);
        stList_append(trees, tree);
    }

    stMatrixDiffs_destruct(snpDiffs);
    stMatrixDiffs_destruct(breakpointDiffs);

    // Get the best-scoring tree.
    double maxScore = -INFINITY;
    stTree *bestTree = NULL;
    for (int64_t i = 0; i < stList_length(trees); i++) {
        stTree *tree = stList_get(trees, i);
        double score = scoreTree(tree, params->scoringMethod,
                                 input->constants->speciesStTree, block,
                                 input->constants->flower, featureColumns,
                                 input->constants->eventToSpeciesNode);
        if (score > maxScore) {
            maxScore = score;
            bestTree = tree;
        }
    }

    if(bestTree == NULL) {
        // Can happen if/when the nucleotide likelihood score
        // is used and a block is all N's. Just use the
        // canonical NJ tree in that case.
        bestTree = blockTree;
    }

    assert(bestTree != NULL);

    // Update the bootstrap support for each branch.
    stTree *bootstrapped = stPhylogeny_scoreReconciliationFromBootstraps(bestTree, trees);

    // Cleanup
    for (int64_t i = 0; i < stList_length(trees); i++) {
        stTree *tree = stList_get(trees, i);
        stPhylogenyInfo_destructOnTree(tree);
        stTree_destruct(tree);
    }
    stList_destruct(trees);
    stList_destruct(featureColumns);
    stList_destruct(featureBlocks);
    stList_destruct(outgroups);
    free(input);

    ret->tree = bootstrapped;

    return ret;
}

// Gets run as a "finisher" in the thread pool, so it's run in series
// and we don't have to lock the hash.
static void addBlockTreeToHash(TreeBuildingResult *result) {
    if (stHash_search(result->blocksToTrees, result->block)) {
        stHash_remove(result->blocksToTrees, result->block);
    }
    if (result->tree != NULL) {
        stHash_insert(result->blocksToTrees, result->block, result->tree);
    } else {
        // Need to take this section out as soon as the debug prints
        // are no longer helpful, as this is in a critical section!
        if (hasSimplePhylogeny(result->block, result->constants->flower)) {
            numSimpleBlocksSkipped++;
        } else if (isSingleCopyBlock(result->block, result->constants->flower)) {
            numSingleCopyBlocksSkipped++;
        }
    }
    free(result); // Sucks to have to do this in a critical section,
                  // but shouldn't matter too much.
}

void splitBlockOnSplitBranch(stPinchBlock *block,
                             stCaf_SplitBranch *splitBranch,
                             stSortedSet *splitBranches,
                             TreeBuildingConstants *constants,
                             stSet *speciesToSplitOn,
                             stThreadPool *treeBuildingPool,
                             stHash *blocksToTrees,
                             stSet *blocksToUpdate) {
    block = splitBranch->block;
    // Do the partition on the block.
    stList *partition = stList_construct();
    // Create a leaf set with all leaves below this split branch.
    stList *leafSet = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    int64_t segmentBelowBranchIndex = -1; // Arbitrary index of
                                          // segment below the
                                          // branch so we can
                                          // recover the blocks we
                                          // split into later.
    stList *bfQueue = stList_construct();
    stList_append(bfQueue, splitBranch->child);
    while (stList_length(bfQueue) != 0) {
        stTree *node = stList_pop(bfQueue);
        for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
            stList_append(bfQueue, stTree_getChild(node, i));
        }

        if (stTree_getChildNumber(node) == 0) {
            stPhylogenyInfo *info = stTree_getClientData(node);
            assert(info->index->matrixIndex != -1);
            stList_append(leafSet, stIntTuple_construct1(info->index->matrixIndex));
            segmentBelowBranchIndex = info->index->matrixIndex;
        }
    }
    stList_append(partition, leafSet);
    // Create a leaf set with all leaves that aren't below the
    // split branch.
    leafSet = stList_construct3(0, (void (*)(void *))stIntTuple_destruct);
    stTree *root = splitBranch->child;
    while (stTree_getParent(root) != NULL) {
        root = stTree_getParent(root);
    }
    int64_t segmentNotBelowBranchIndex = -1; // Arbitrary index of
                                             // segment not below the
                                             // branch so we can
                                             // recover the blocks
                                             // we split into
                                             // later.
    stList_append(bfQueue, root);
    while (stList_length(bfQueue) != 0) {
        stTree *node = stList_pop(bfQueue);
        if (node != splitBranch->child) {
            for (int64_t i = 0; i < stTree_getChildNumber(node); i++) {
                stList_append(bfQueue, stTree_getChild(node, i));
            }

            if (stTree_getChildNumber(node) == 0) {
                stPhylogenyInfo *info = stTree_getClientData(node);
                assert(info->index->matrixIndex != -1);
                stList_append(leafSet, stIntTuple_construct1(info->index->matrixIndex));
                segmentNotBelowBranchIndex = info->index->matrixIndex;
            }
        }
    }
    stList_append(partition, leafSet);
    // Get arbitrary segments from the blocks we will split
    // into. This is so that we can recover the blocks after the
    // split.
    stPinchSegment *segmentBelowBranch = getSegmentByBlockIndex(block, segmentBelowBranchIndex);
    stPinchSegment *segmentNotBelowBranch = getSegmentByBlockIndex(block, segmentNotBelowBranchIndex);

    // Remove all split branch entries for this block tree.
    removeOldSplitBranches(block, root, speciesToSplitOn,
                           splitBranches);

    // Destruct the block tree.
    // This also invalidates "splitBranch" from here on!
    stPhylogenyInfo_destructOnTree(root);
    stTree_destruct(root);
    stHash_remove(blocksToTrees, block);

    // Actually perform the split according to the partition.
    splitBlock(block, partition, constants->params->keepSingleDegreeBlocks);
    // Recover the blocks.
    stPinchBlock *blockBelowBranch = stPinchSegment_getBlock(segmentBelowBranch);
    stPinchBlock *blockNotBelowBranch = stPinchSegment_getBlock(segmentNotBelowBranch);

    // Make a new tree for both of the blocks we split into, if
    // they are not too simple to make a tree.
    if (blockBelowBranch != NULL) {
        pushBlockToPool(blockBelowBranch, constants, blocksToTrees,
                        treeBuildingPool);
    }
    if (blockNotBelowBranch != NULL) {
        pushBlockToPool(blockNotBelowBranch, constants, blocksToTrees,
                        treeBuildingPool);
    }
    stThreadPool_wait(treeBuildingPool);

    stTree *treeNotBelowBranch = stHash_search(blocksToTrees,
                                               blockNotBelowBranch);
    stTree *treeBelowBranch = stHash_search(blocksToTrees, blockBelowBranch);
    if (treeBelowBranch != NULL) {
        findSplitBranches(blockBelowBranch, treeBelowBranch,
                          splitBranches, speciesToSplitOn);
    }
    if (treeNotBelowBranch != NULL) {
        findSplitBranches(blockNotBelowBranch, treeNotBelowBranch,
                          splitBranches, speciesToSplitOn);
    }

    // Finally, update the trees for all blocks close enough to
    // either of the new blocks to be affected by this breakpoint
    // change.
    if (blockBelowBranch != NULL) {
        addContextualBlocksToSet(blockBelowBranch,
                                 constants->params->maxBaseDistance,
                                 constants->params->maxBlockDistance,
                                 constants->params->ignoreUnalignedBases,
                                 blocksToUpdate);
    }
    if (blockNotBelowBranch != NULL) {
        addContextualBlocksToSet(blockNotBelowBranch,
                                 constants->params->maxBaseDistance,
                                 constants->params->maxBlockDistance,
                                 constants->params->ignoreUnalignedBases,
                                 blocksToUpdate);
    }
}

static void destructBlockTree(stTree *tree) {
    stPhylogenyInfo_destructOnTree(tree);
    stTree_destruct(tree);
}

void stCaf_buildTreesToRemoveAncientHomologies(stPinchThreadSet *threadSet,
                                               stHash *threadStrings,
                                               stSet *outgroupThreads,
                                               Flower *flower,
                                               stCaf_PhylogenyParameters *params,
                                               FILE *debugFile,
                                               const char *referenceEventHeader) {
    // Functions we aren't using right now but should stick around anyway.
    (void) printSimpleBlockDebugInfo;
    (void) getTotalSimilarityAndDifferenceCounts;
    (void) printTreeBuildingDebugInfo;

    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block;

    //Hash in which we store a map of blocks to their trees
    stHash *blocksToTrees = stHash_construct2(NULL, (void (*)(void *)) destructBlockTree);

    //Get species tree as an stTree
    EventTree *eventTree = flower_getEventTree(flower);
    stTree *speciesStTree = eventTree_getStTree(eventTree);

    // Get info for reconciliation.
    stHash *eventToSpeciesNode = getEventToSpeciesNode(eventTree, speciesStTree);

    // Get info for guided neighbor-joining
    stHash *speciesToJoinCostIndex = stHash_construct2(NULL, (void (*)(void *)) stIntTuple_destruct);
    stMatrix *joinCosts = stPhylogeny_computeJoinCosts(
        speciesStTree, speciesToJoinCostIndex,
        params->costPerDupPerBase * 2 * params->maxBaseDistance,
        params->costPerLossPerBase * 2 * params->maxBaseDistance);
    int64_t **speciesMRCAMatrix = stPhylogeny_getMRCAMatrix(speciesStTree, speciesToJoinCostIndex);

    stSortedSet *splitBranches = stSortedSet_construct3((int (*)(const void *, const void *)) compareSplitBranches, free);
    stSet *speciesToSplitOn = stSet_construct();
    getSpeciesToSplitOn(speciesStTree, eventTree, referenceEventHeader,
                        speciesToSplitOn);

    // Temp debug print
    printf("Chose events in the species tree \"%s\" to split on:", eventTree_makeNewickString(eventTree));
    stSetIterator *speciesToSplitOnIt = stSet_getIterator(speciesToSplitOn);
    stTree *speciesNodeToSplitOn;
    while ((speciesNodeToSplitOn = stSet_getNext(speciesToSplitOnIt)) != NULL) {
        Name name;
        sscanf(stTree_getLabel(speciesNodeToSplitOn), "%" PRIi64, &name);
        Event *event = eventTree_getEvent(eventTree, name);
        printf(" %s", event_getHeader(event));
    }
    printf("\n");
    stSet_destructIterator(speciesToSplitOnIt);

    stThreadPool *treeBuildingPool = stThreadPool_construct(
        params->numTreeBuildingThreads,
        (void *(*)(void *)) buildTreeForBlock,
        (void (*)(void *)) addBlockTreeToHash);
    TreeBuildingConstants constants;
    constants.threadStrings = threadStrings;
    constants.outgroupThreads = outgroupThreads;
    constants.flower = flower;
    constants.params = params;
    constants.joinCosts = joinCosts;
    constants.speciesToJoinCostIndex = speciesToJoinCostIndex;
    constants.speciesMRCAMatrix = speciesMRCAMatrix;
    constants.eventToSpeciesNode = eventToSpeciesNode;
    constants.speciesStTree = speciesStTree;

    // The loop to build a tree for each block
    while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
        pushBlockToPool(block, &constants, blocksToTrees, treeBuildingPool);
    }

    // We need the trees to be done before we can continue.
    stThreadPool_wait(treeBuildingPool);

    fprintf(stdout, "First tree-building round done. Skipped %" PRIi64
            " simple blocks (those with 1 event or with < 3 segments, and %"
            PRIi64 " single-copy blocks (those with (# events) = (# segments))."
            "\n", numSimpleBlocksSkipped, numSingleCopyBlocksSkipped);

    // All the blocks have their trees computed. Find the split
    // branches in those trees.
    stHashIterator *blocksToTreesIt = stHash_getIterator(blocksToTrees);
    while ((block = stHash_getNext(blocksToTreesIt)) != NULL) {
        stTree *tree = stHash_search(blocksToTrees, block);
        assert(tree != NULL);
        findSplitBranches(block, tree, splitBranches, speciesToSplitOn);
    }

    fprintf(stdout, "Before partitioning, there were %" PRIi64 " bases lost in between single-degree blocks\n", countBasesBetweenSingleDegreeBlocks(threadSet));
    fprintf(stdout, "Found %" PRIi64 " split branches initially in %" PRIi64
            " blocks (%" PRIi64 " of which have trees), with an "
            "average support value of %lf\n", stSortedSet_size(splitBranches),
            stPinchThreadSet_getTotalBlockNumber(threadSet),
            stHash_size(blocksToTrees),
            getAvgSupportValue(splitBranches));
    // Now walk through the split branches, doing the most confident
    // splits first, and updating the blocks whose breakpoint
    // information is modified.
    int64_t numberOfSplitsMade = 0;
    double totalSupport = 0;
    int64_t totalNumberOfBlocksRecomputed = 0;
    stCaf_SplitBranch *splitBranch = stSortedSet_getLast(splitBranches);
    if (params->doSplitsWithSupportHigherThanThisAllAtOnce <= 1.0) {
        // Save a lot of time by doing the large fraction of highly
        // confident splits first, then updating the other affected
        // blocks in one go.
        stSet *blocksToUpdate = stSet_construct();
        stSet *blocksSplit = stSet_construct();
        while (splitBranch != NULL && splitBranch->support >= params->doSplitsWithSupportHigherThanThisAllAtOnce) {
            totalSupport += splitBranch->support;
            block = splitBranch->block;
            stSet_insert(blocksSplit, block);
            splitBlockOnSplitBranch(block, splitBranch, splitBranches,
                                    &constants, speciesToSplitOn,
                                    treeBuildingPool, blocksToTrees,
                                    blocksToUpdate);
            // Get the next split branch.
            splitBranch = stSortedSet_getLast(splitBranches);
            numberOfSplitsMade++;
        }

        // Update the blocks that were just affected by all the
        // perfectly supported split branches.
        // But first ensure that the blocks weren't already split!
        stSetIterator *blocksToUpdateIt = stSet_getIterator(blocksToUpdate);
        stPinchBlock *blockToUpdate;
        while ((blockToUpdate = stSet_getNext(blocksToUpdateIt)) != NULL) {
            if (stSet_search(blocksSplit, blockToUpdate)) {
                // This block is invalid!
                continue;
            }
            totalNumberOfBlocksRecomputed++;
            stTree *oldTree = stHash_search(blocksToTrees, blockToUpdate);
            removeOldSplitBranches(blockToUpdate, oldTree, speciesToSplitOn,
                                   splitBranches);
            pushBlockToPool(blockToUpdate, &constants, blocksToTrees,
                            treeBuildingPool);
        }
        stSet_destructIterator(blocksToUpdateIt);
        // Wait for the trees to be done.
        stThreadPool_wait(treeBuildingPool);
        blocksToUpdateIt = stSet_getIterator(blocksToUpdate);
        while ((blockToUpdate = stSet_getNext(blocksToUpdateIt)) != NULL) {
            stTree *tree = stHash_search(blocksToTrees, blockToUpdate);
            if (tree != NULL) {
                findSplitBranches(blockToUpdate, tree,
                                  splitBranches, speciesToSplitOn);
            }
        }
        stSet_destructIterator(blocksToUpdateIt);
        stSet_destruct(blocksToUpdate);
    }
    splitBranch = stSortedSet_getLast(splitBranches);
    while (splitBranch != NULL) {
        totalSupport += splitBranch->support;
        block = splitBranch->block;
        stSet *blocksToUpdate = stSet_construct();
        splitBlockOnSplitBranch(block, splitBranch, splitBranches,
                                &constants, speciesToSplitOn, treeBuildingPool,
                                blocksToTrees, blocksToUpdate);
        stSetIterator *blocksToUpdateIt = stSet_getIterator(blocksToUpdate);
        stPinchBlock *blockToUpdate;
        while ((blockToUpdate = stSet_getNext(blocksToUpdateIt)) != NULL) {
            totalNumberOfBlocksRecomputed++;
            stTree *oldTree = stHash_search(blocksToTrees, blockToUpdate);
            removeOldSplitBranches(blockToUpdate, oldTree, speciesToSplitOn,
                                   splitBranches);
            pushBlockToPool(blockToUpdate, &constants, blocksToTrees,
                            treeBuildingPool);
        }
        stSet_destructIterator(blocksToUpdateIt);
        stThreadPool_wait(treeBuildingPool);
        blocksToUpdateIt = stSet_getIterator(blocksToUpdate);
        while ((blockToUpdate = stSet_getNext(blocksToUpdateIt)) != NULL) {
            stTree *tree = stHash_search(blocksToTrees, blockToUpdate);
            if (tree != NULL) {
                findSplitBranches(blockToUpdate, tree,
                                  splitBranches, speciesToSplitOn);
            }
        }
        stSet_destructIterator(blocksToUpdateIt);
        stSet_destruct(blocksToUpdate);
        numberOfSplitsMade++;
        splitBranch = stSortedSet_getLast(splitBranches);
    }

    st_logDebug("Finished partitioning the blocks\n");
    fprintf(stdout, "There were %" PRIi64 " splits made overall in the end.\n",
            numberOfSplitsMade);
    fprintf(stdout, "The split branches that we actually used had an average "
            "support of %lf.\n",
            numberOfSplitsMade != 0 ? totalSupport/numberOfSplitsMade : 0.0);
    fprintf(stdout, "We recomputed the trees for an average of %" PRIi64
            " blocks after every split, and the pinch graph had a total of "
            "%" PRIi64 " blocks.\n",
            numberOfSplitsMade != 0 ? totalNumberOfBlocksRecomputed/numberOfSplitsMade : 0,
            stPinchThreadSet_getTotalBlockNumber(threadSet));
    fprintf(stdout, "After partitioning, there were %" PRIi64 " bases lost in between single-degree blocks\n", countBasesBetweenSingleDegreeBlocks(threadSet));
    fprintf(stdout, "We stopped %" PRIi64 " single-degree segments from becoming"
            " blocks (avg %lf per block) for a total of %" PRIi64 " bases\n",
            numSingleDegreeSegmentsDropped,
            ((float)numSingleDegreeSegmentsDropped)/stPinchThreadSet_getTotalBlockNumber(threadSet),
            numBasesDroppedFromSingleDegreeSegments);

    //Cleanup
    for (int64_t i = 0; i < stTree_getNumNodes(speciesStTree); i++) {
        free(speciesMRCAMatrix[i]);
    }
    free(speciesMRCAMatrix);
    stTree_destruct(speciesStTree);
    stThreadPool_destruct(treeBuildingPool);
    stHash_destruct(blocksToTrees);
}
