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
static void splitBlock(stPinchBlock *block, stList *partitions) {
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

/*
 * Get a gene node->species node mapping from a gene tree, a species
 * tree, and the pinch block.
 */

static stHash *getLeafToSpecies(stTree *geneTree, stTree *speciesTree,
                                stPinchBlock *block, Flower *flower) {
    stHash *leafToSpecies = stHash_construct();
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment;
    int64_t i = 0; // Current segment index in block.
    while((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        stPinchThread *thread = stPinchSegment_getThread(segment);
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        Event *event = cap_getEvent(cap);
        char *eventNameString = stString_print("%" PRIi64, event_getName(event));
        stTree *species = stTree_findChild(speciesTree, eventNameString);
        free(eventNameString);
        assert(species != NULL);
        stTree *gene = stPhylogeny_getLeafByIndex(geneTree, i);
        assert(gene != NULL);
        stHash_insert(leafToSpecies, gene, species);
        i++;
    }
    return leafToSpecies;
}

static stTree *eventTreeToStTree_R(Event *event) {
    stTree *ret = stTree_construct();
    stTree_setLabel(ret, stString_print("%" PRIi64, event_getName(event)));
    stTree_setBranchLength(ret, event_getBranchLength(event));
    for(int64_t i = 0; i < event_getChildNumber(event); i++) {
        Event *child = event_getChild(event, i);
        stTree *childStTree = eventTreeToStTree_R(child);
        stTree_setParent(childStTree, ret);
    }
    return ret;
}

// Get species tree from event tree (labeled by the event Names),
// which requires ignoring the root event.
static stTree *eventTreeToStTree(EventTree *eventTree) {
    Event *rootEvent = eventTree_getRootEvent(eventTree);
    // Need to skip the root event, since it is added onto the real
    // species tree.
    assert(event_getChildNumber(rootEvent) == 1);
    Event *speciesRoot = event_getChild(rootEvent, 0);
    return eventTreeToStTree_R(speciesRoot);
}

static double scoreTree(stTree *tree, enum stCaf_ScoringMethod scoringMethod, stTree *speciesStTree, stPinchBlock *block, Flower *flower, stList *featureColumns) {
    double ret = 0.0;
    if(scoringMethod == RECON_COST) {
        stHash *leafToSpecies = getLeafToSpecies(tree,
                                                 speciesStTree,
                                                 block, flower);
        int64_t dups, losses;
        stPinchPhylogeny_reconciliationCostBinary(tree, speciesStTree,
                                                  leafToSpecies, &dups,
                                                  &losses);
        ret = -dups - losses;

        stHash_destruct(leafToSpecies);
    } else if(scoringMethod == NUCLEOTIDE_LIKELIHOOD) {
        ret = stPinchPhylogeny_likelihood(tree, featureColumns);
    }
    return ret;
}

// Build a tree from a set of feature columns and root it according to
// the rooting method.
static stTree *buildTree(stList *featureColumns,
                         enum stCaf_RootingMethod rootingMethod,
                         double breakPointScalingFactor,
                         bool bootstrap,
                         stList *outgroups, stPinchBlock *block,
                         Flower *flower, stTree *speciesStTree) {
    // Make substitution matrix
    stMatrix *substitutionMatrix = stPinchPhylogeny_getMatrixFromSubstitutions(featureColumns, block, NULL, bootstrap);
    assert(stMatrix_n(substitutionMatrix) == stPinchBlock_getDegree(block));
    assert(stMatrix_m(substitutionMatrix) == stPinchBlock_getDegree(block));
    //Make breakpoint matrix
    stMatrix *breakpointMatrix = stPinchPhylogeny_getMatrixFromBreakpoints(featureColumns, block, NULL, bootstrap);
    
    //Combine the matrices into distance matrices
    stMatrix_scale(breakpointMatrix, breakPointScalingFactor, 0.0);
    stMatrix *combinedMatrix = stMatrix_add(substitutionMatrix, breakpointMatrix);
    stMatrix *distanceMatrix = stPinchPhylogeny_getSymmetricDistanceMatrix(combinedMatrix);
    
    stTree *tree = NULL;
    if(rootingMethod == OUTGROUP_BRANCH) {
        tree = stPhylogeny_neighborJoin(distanceMatrix, outgroups);
    } else if(rootingMethod == LONGEST_BRANCH) {
        tree = stPhylogeny_neighborJoin(distanceMatrix, NULL);
    } else if(rootingMethod == BEST_RECON) {
        tree = stPhylogeny_neighborJoin(distanceMatrix, NULL);
        stHash *leafToSpecies = getLeafToSpecies(tree,
                                                 speciesStTree,
                                                 block, flower);
        stTree *newTree = stPinchPhylogeny_reconcileBinary(tree, speciesStTree, leafToSpecies);

        stPhylogenyInfo_destructOnTree(tree);
        stTree_destruct(tree);
        stHash_destruct(leafToSpecies);
        tree = newTree;
    }

    return tree;
}

// Check if the block's phylogeny is simple:
// - the block has only one event, or
// - the block has < 3 segments, or
// - the block does not contain any segments that are part of an
//   outgroup thread.
static bool hasSimplePhylogeny(stPinchBlock *block,
                               stSet *outgroupThreads,
                               Flower *flower) {
    if(stPinchBlock_getDegree(block) <= 2) {
        return true;
    }
    stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
    stPinchSegment *segment = NULL;
    bool foundOutgroup = 0, found2Events = 0;
    Event *currentEvent = NULL;
    while((segment = stPinchBlockIt_getNext(&blockIt)) != NULL) {
        stPinchThread *thread = stPinchSegment_getThread(segment);
        if(stSet_search(outgroupThreads, thread) != NULL) {
            foundOutgroup = 1;
        }
        Cap *cap = flower_getCap(flower, stPinchThread_getName(thread));
        assert(cap != NULL);
        Event *event = cap_getEvent(cap);
        if(currentEvent == NULL) {
            currentEvent = event;
        } else if(currentEvent != event) {
            found2Events = 1;
        }
    }
    return !(foundOutgroup && found2Events);
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

void stCaf_buildTreesToRemoveAncientHomologies(stPinchThreadSet *threadSet, stHash *threadStrings, stSet *outgroupThreads, Flower *flower, int64_t numTrees, enum stCaf_RootingMethod rootingMethod, enum stCaf_ScoringMethod scoringMethod, double breakPointScalingFactor, bool skipSingleCopyBlocks) {
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    stPinchBlock *block;

    //Hash in which we store a map of blocks to the partitions
    stHash *blocksToPartitions = stHash_construct2(NULL, NULL);

    //Get species tree as an stTree
    EventTree *eventTree = flower_getEventTree(flower);
    stTree *speciesStTree = eventTreeToStTree(eventTree);

    //Count of the total number of blocks partitioned by an ancient homology
    int64_t totalBlocksSplit = 0;
    int64_t totalSingleCopyBlocksSplit = 0;
    double totalSubstitutionSimilarities = 0.0;
    double totalSubstitutionDifferences = 0.0;
    double totalBreakpointSimilarities = 0.0;
    double totalBreakpointDifferences = 0.0;
    double totalBestTreeScore = 0.0;
    double totalTreeScore = 0.0;
    int64_t sampledTreeWasBetterCount = 0;
    int64_t totalOutgroupThreads = 0;
    int64_t totalSimpleBlocks = 0;
    int64_t totalSingleCopyBlocks = 0;

    //The loop to build a tree for each block
    while ((block = stPinchThreadSetBlockIt_getNext(&blockIt)) != NULL) {
        if(!hasSimplePhylogeny(block, outgroupThreads, flower)) { //No point trying to build a phylogeny for certain blocks.
            bool singleCopy = 0;
            if(isSingleCopyBlock(block, flower)) {
                singleCopy = 1;
                totalSingleCopyBlocks++;
                if(skipSingleCopyBlocks) {
                    continue;
                }
            }
            //Parameters.
            int64_t maxBaseDistance = 1000;
            int64_t maxBlockDistance = 100;
            bool ignoreUnalignedBases = 1;
            bool onlyIncludeCompleteFeatureBlocks = 0;

            //Get the feature blocks
            stList *featureBlocks = stFeatureBlock_getContextualFeatureBlocks(block, maxBaseDistance, maxBlockDistance,
                    ignoreUnalignedBases, onlyIncludeCompleteFeatureBlocks, threadStrings);

            //Make feature columns
            stList *featureColumns = stFeatureColumn_getFeatureColumns(featureBlocks, block);

            //Make substitution matrix
            stMatrix *substitutionMatrix = stPinchPhylogeny_getMatrixFromSubstitutions(featureColumns, block, NULL, 0);
            assert(stMatrix_n(substitutionMatrix) == stPinchBlock_getDegree(block));
            assert(stMatrix_m(substitutionMatrix) == stPinchBlock_getDegree(block));
            double similarities, differences;
            getTotalSimilarityAndDifferenceCounts(substitutionMatrix, &similarities, &differences);
            totalSubstitutionSimilarities += similarities;
            totalSubstitutionDifferences += differences;

            //Make breakpoint matrix
            stMatrix *breakpointMatrix = stPinchPhylogeny_getMatrixFromBreakpoints(featureColumns, block, NULL, 0);
            getTotalSimilarityAndDifferenceCounts(breakpointMatrix, &similarities, &differences);
            totalBreakpointSimilarities += similarities;
            totalBreakpointDifferences += differences;

            //Get the outgroup threads
            stList *outgroups = getOutgroupThreads(block, outgroupThreads);
            totalOutgroupThreads += stList_length(outgroups);

            //Build the canonical tree.
            stTree *blockTree = buildTree(featureColumns, rootingMethod,
                                          breakPointScalingFactor,
                                          0, outgroups, block, flower,
                                          speciesStTree);

            // Sample the rest of the trees.
            stList *trees = stList_construct();
            stList_append(trees, blockTree);
            for(int64_t i = 0; i < numTrees - 1; i++) {
                stTree *tree = buildTree(featureColumns, rootingMethod,
                                         breakPointScalingFactor,
                                         1, outgroups, block, flower,
                                         speciesStTree);
                stList_append(trees, tree);
            }

            // Get the best-scoring tree.
            double maxScore = -INFINITY;
            stTree *bestTree = NULL;
            for(int64_t i = 0; i < stList_length(trees); i++) {
                stTree *tree = stList_get(trees, i);
                double score = scoreTree(tree, scoringMethod,
                                         speciesStTree, block, flower,
                                         featureColumns);
                if(score != -INFINITY) { // avoid counting impossible
                                         // trees
                    totalTreeScore += score;
                }
                if(score > maxScore) {
                    sampledTreeWasBetterCount += (i == 0) ? 0 : 1;
                    maxScore = score;
                    bestTree = tree;
                }
            }
            assert(bestTree != NULL);

            totalBestTreeScore += maxScore;
            //Get the partitions
            stList *partition = stPinchPhylogeny_splitTreeOnOutgroups(bestTree, outgroups);
            if(stList_length(partition) > 1) {
                if(singleCopy) {
                    totalSingleCopyBlocksSplit++;
                }
                totalBlocksSplit++;
            }
            stHash_insert(blocksToPartitions, block, partition);

            //Cleanup
            stMatrix_destruct(substitutionMatrix);
            stMatrix_destruct(breakpointMatrix);
            for(int64_t i = 0; i < stList_length(trees); i++) {
                stTree *tree = stList_get(trees, i);
                stPhylogenyInfo_destructOnTree(tree);
                stTree_destruct(tree);
            }
            stList_destruct(featureColumns);
            stList_destruct(featureBlocks);
            stList_destruct(outgroups);
        } else {
            totalSimpleBlocks++;
        }
    }
    // Block count including skipped blocks.
    int64_t totalBlockCount = stPinchThreadSet_getTotalBlockNumber(threadSet);
    // Number of blocks that were actually considered while partitioning.
    int64_t blockCount = totalBlockCount - totalSimpleBlocks;
    fprintf(stdout, "Using phylogeny building, of %" PRIi64 " blocks considered (%" PRIi64 " total), %" PRIi64 " blocks were partitioned\n", blockCount, totalBlockCount, totalBlocksSplit);
    fprintf(stdout, "There were %" PRIi64 " outgroup threads seen total over all blocks\n", totalOutgroupThreads);
    fprintf(stdout, "In phylogeny building there were %f avg. substitution similarities %f avg. substitution differences\n", totalSubstitutionSimilarities/blockCount, totalSubstitutionDifferences/blockCount);
    fprintf(stdout, "In phylogeny building there were %f avg. breakpoint similarities %f avg. breakpoint differences\n", totalBreakpointSimilarities/blockCount, totalBreakpointDifferences/blockCount);
    fprintf(stdout, "In phylogeny building we saw an average score of %f for the best tree in each block, an average score of %f overall, and %" PRIi64 " trees total.\n", totalBestTreeScore/blockCount, totalTreeScore/(numTrees*blockCount), numTrees*blockCount);
    fprintf(stdout, "In phylogeny building we used a sampled tree instead of the canonical tree %" PRIi64 " times.\n", sampledTreeWasBetterCount);
    fprintf(stdout, "In phylogeny building we skipped %" PRIi64 " simple blocks.\n", totalSimpleBlocks);
    fprintf(stdout, "In phylogeny building there were %" PRIi64 " single copy blocks considered, of which %" PRIi64 " were partitioned.\n", totalSingleCopyBlocks, totalSingleCopyBlocksSplit);

    st_logDebug("Got homology partition for each block\n");

    //Now walk through the blocks and do the actual splits, must be done after the fact using the blocks
    //in the original hash, as we are now disrupting and changing the original graph.
    stHashIterator *blockIt2 = stHash_getIterator(blocksToPartitions);
    while ((block = stHash_getNext(blockIt2)) != NULL) {
        stList *partition = stHash_search(blocksToPartitions, block);
        assert(partition != NULL);
        splitBlock(block, partition);
        stList_destruct(partition);
    }
    stHash_destructIterator(blockIt2);

    st_logDebug("Finished partitioning the blocks\n");

    //Cleanup
    stHash_destruct(blocksToPartitions);
    stTree_destruct(speciesStTree);
}
