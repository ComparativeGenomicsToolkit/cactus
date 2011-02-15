/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"

/*
 * Stats for a cactus tree that passes cactus_check.
 */

void tabulateFloatStats(struct List *unsortedValues, double *totalNumber,
        double *totalSum, double *min, double *max, double *avg, double *median) {
    /*
     * Calculates basic stats from a list of float values.
     */
    if (unsortedValues->length == 0) {
        *totalNumber = 0;
        *totalSum = 0;
        *min = INT32_MAX;
        *max = INT32_MAX;
        *avg = INT32_MAX;
        *median = INT32_MAX;
        return;
    }
    unsortedValues = listCopy(unsortedValues); //copy the input list, to avoid altering the input.
    assert(unsortedValues->length > 0);
    qsort(unsortedValues->list, unsortedValues->length, sizeof(void *),
            (int(*)(const void *, const void *)) floatComparator);
    *totalNumber = unsortedValues->length;
    *min = *(float *) unsortedValues->list[0];
    *max = *(float *) unsortedValues->list[unsortedValues->length - 1];
    *median = *(float *) unsortedValues->list[unsortedValues->length / 2];
    int32_t i;
    float j = 0;
    for (i = 0; i < unsortedValues->length; i++) {
        j += *(float *) unsortedValues->list[i];
    }
    *avg = j / unsortedValues->length;
    *totalSum = j;
    unsortedValues->destructElement = NULL;
    destructList(unsortedValues);
}

void tabulateStats(struct IntList *unsortedValues, double *totalNumber,
        double *totalSum, double *min, double *max, double *avg, double *median) {
    /*
     * Same as float stats, but for an intlist.
     */
    if (unsortedValues->length == 0) {
        *totalNumber = 0;
        *totalSum = 0;
        *min = INT32_MAX;
        *max = INT32_MAX;
        *avg = INT32_MAX;
        *median = INT32_MAX;
        return;
    }
    unsortedValues = intListCopy(unsortedValues);
    assert(unsortedValues->length > 0);
    qsort(unsortedValues->list, unsortedValues->length, sizeof(int32_t),
            (int(*)(const void *, const void *)) intComparator_Int);
    *totalNumber = unsortedValues->length;
    *min = unsortedValues->list[0];
    *max = unsortedValues->list[unsortedValues->length - 1];
    *median = unsortedValues->list[unsortedValues->length / 2];
    int32_t i, j = 0;
    for (i = 0; i < unsortedValues->length; i++) {
        j += unsortedValues->list[i];
    }
    *avg = (double) j / unsortedValues->length;
    *totalSum = j;
    destructIntList(unsortedValues);
}

void printOpeningTag(const char *tag, FILE *fileHandle) {
    /*
     * Creates an opening XML tag.
     */
    fprintf(fileHandle, "<%s>", tag);
}

void printClosingTag(const char *tag, FILE *fileHandle) {
    /*
     * Creates a closing XML tag.
     */
    fprintf(fileHandle, "</%s>", tag);
}

void tabulateAndPrintFloatValues(struct List *values, const char *tag,
        FILE *fileHandle) {
    /*
     * Creates a node containing basic stats on the given float values and a nested "values" node containing the actual values.
     */
    double totalNumber, totalSum, min, max, avg, median;
    tabulateFloatStats(values, &totalNumber, &totalSum, &min, &max, &avg,
            &median);
    fprintf(
            fileHandle,
            "<%s total=\"%f\" sum=\"%f\" min=\"%f\" max=\"%f\" avg=\"%f\" median=\"%f\">",
            tag, totalNumber, totalSum, min, max, avg, median);
    int32_t i;
    for (i = 0; i < values->length; i++) {
        fprintf(fileHandle, "%f ", *(float *) values->list[i]);
    }
    printClosingTag(tag, fileHandle);
}

void tabulateAndPrintIntValues(struct IntList *values, const char *tag,
        FILE *fileHandle) {
    /*
     * Creates a node containing basic stats on the given int values and a nested "values" node containing the actual values.
     */
    double totalNumber, totalSum, min, max, avg, median;
    tabulateStats(values, &totalNumber, &totalSum, &min, &max, &avg, &median);
    fprintf(
            fileHandle,
            "<%s total=\"%f\" sum=\"%f\" min=\"%f\" max=\"%f\" avg=\"%f\" median=\"%f\">",
            tag, totalNumber, totalSum, min, max, avg, median);
    int32_t i;
    for (i = 0; i < values->length; i++) {
        fprintf(fileHandle, "%i ", values->list[i]);
    }
    printClosingTag(tag, fileHandle);
}

struct IntList *convertToIntList(stList *list) {
    struct IntList *tempList = constructEmptyIntList(0);
    for (int32_t i = 0; i < stList_length(list); i++) {
        intListAppend(tempList, stIntTuple_getPosition(stList_get(list, i), 0));
    }
    return tempList;
}

void tabulateAndPrintIntTupleValues(stList *list, const char *name,
        FILE *fileHandle) {
    struct IntList *tempList = convertToIntList(list);
    tabulateAndPrintIntValues(tempList, name, fileHandle);
    destructIntList(tempList);
}

/////
//Now on to the actual stats
/////

double calculateTreeBits(Flower *flower, double pathBitScore) {
    /*
     * Calculates the total number of bits to required to encode the path to every base in the flower.
     */
    double totalBitScore = 0.0;
    int32_t totalSequenceSize;
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    double followingPathBitScore = (log(flower_getGroupNumber(flower)) / log(
            2.0)) + pathBitScore;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        if (group_isLeaf(group)) {
            totalSequenceSize = group_getTotalBaseLength(group);
            totalBitScore += (totalSequenceSize > 0 ? ((log(totalSequenceSize)
                    / log(2.0)) + followingPathBitScore) * totalSequenceSize
                    : 0.0);
        } else {
            totalBitScore += calculateTreeBits(group_getNestedFlower(group),
                    followingPathBitScore);
        }
    }
    flower_destructGroupIterator(groupIterator);
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    totalSequenceSize = 0.0;
    while ((block = flower_getNextBlock(blockIterator)) != NULL) {
        totalSequenceSize += block_getLength(block) * block_getInstanceNumber(
                block);
    }
    flower_destructBlockIterator(blockIterator);
    return totalBitScore + (totalSequenceSize > 0 ? ((log(totalSequenceSize)
            / log(2.0)) + pathBitScore) * totalSequenceSize : 0.0);
}

void reportRelativeEntopyStats(Flower *flower, FILE *fileHandle) {
    /*
     * Relative entropy stats. Supposed to give a metric of how balanced the tree is in how it subdivides the input sequences.
     */
    double totalSeqSize = flower_getTotalBaseLength(flower);
    double totalP = calculateTreeBits(flower, 0.0);
    double totalQ = (log(totalSeqSize) / log(2.0)) * totalSeqSize;
    //assert(totalP >= totalQ);
    double relativeEntropy = totalP - totalQ;
    double normalisedRelativeEntropy = relativeEntropy / totalSeqSize;

    fprintf(
            fileHandle,
            "<relative_entropy_stats totalP=\"%f\" totalQ=\"%f\" relative_entropy=\"%f\" normalised_relative_entropy=\"%f\"/>",
            totalP, totalQ, relativeEntropy, normalisedRelativeEntropy);
}

static void flowerStats(Flower *flower, int32_t currentDepth,
        struct IntList *children, struct IntList *tangleChildren,
        struct IntList *linkChildren, struct IntList *depths) {
    /*
     * Calculates basic stats on flowers.
     * Children is the number of children internal nodes (those with children), have.
     * Tangle children, like children but only including groups that are tangle groups.
     * Link children, like children but only including groups that are link groups.
     * Depth is the length of a path (in terms of edges/connections) from the root flower to a terminal flower (which are the leaf flowers of the tree, if terminally normalised).
     */
    if (!flower_isTerminal(flower)) {
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        int32_t i = 0;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            assert(!group_isLeaf(group));
            flowerStats(group_getNestedFlower(group), currentDepth + 1,
                    children, tangleChildren, linkChildren, depths);
            if (group_getLink(group) != NULL) {
                i++;
            }
        }
        flower_destructGroupIterator(groupIterator);
        intListAppend(children, flower_getGroupNumber(flower));
        intListAppend(tangleChildren, flower_getGroupNumber(flower) - i);
        intListAppend(linkChildren, i);
    } else {
        intListAppend(depths, currentDepth);
    }
}

void reportFlowerStats(Flower *flower, FILE *fileHandle) {
    /*
     * Prints the chain stats to the XML file.
     */
    struct IntList *children = constructEmptyIntList(0);
    struct IntList *tangleChildren = constructEmptyIntList(0);
    struct IntList *linkChildren = constructEmptyIntList(0);
    struct IntList *depths = constructEmptyIntList(0);
    flowerStats(flower, 0, children, tangleChildren, linkChildren, depths);
    printOpeningTag("flowers", fileHandle);
    tabulateAndPrintIntValues(children, "children", fileHandle);
    tabulateAndPrintIntValues(tangleChildren, "tangle_children", fileHandle);
    tabulateAndPrintIntValues(linkChildren, "link_children", fileHandle);
    tabulateAndPrintIntValues(depths, "depths", fileHandle);
    printClosingTag("flowers", fileHandle);
    destructIntList(children);
    destructIntList(tangleChildren);
    destructIntList(linkChildren);
    destructIntList(depths);
}

void blockStats(Flower *flower, struct IntList *counts,
        struct IntList *lengths, struct IntList *degrees,
        struct IntList *leafDegrees, struct IntList *coverage,
        struct IntList *leafCoverage, bool(*includeBlock)(Block *),
        struct IntList *columnDegrees, struct IntList *columnLeafDegrees, bool perColumnStats) {
    /*
     * Calculates stats on the blocks outside of terminal flowers.
     * Counts is numbers of blocks per non-terminal flower.
     * Lengths is lengths of blocks.
     * Degrees is the number of segment instances in each block.
     * Leaf degree is the number of leaf segment instances in each block.
     * Coverage is the length * degree of each block.
     * Leaf coverage is the length * leadf degree of each block.
     */
    if (!flower_isTerminal(flower)) {
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            assert(!group_isLeaf(group));
            blockStats(group_getNestedFlower(group), counts, lengths, degrees,
                    leafDegrees, coverage, leafCoverage, includeBlock,
                    columnDegrees, columnLeafDegrees, perColumnStats);
        }
        flower_destructGroupIterator(groupIterator);
        Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
        Block *block;
        while ((block = flower_getNextBlock(blockIterator)) != NULL) {
            if (includeBlock(block)) {
                intListAppend(lengths, block_getLength(block));
                intListAppend(degrees, block_getInstanceNumber(block));
                intListAppend(coverage, block_getLength(block)
                        * block_getInstanceNumber(block));
                Segment *segment;
                Block_InstanceIterator *segmentIterator =
                        block_getInstanceIterator(block);
                int32_t i = 0;
                while ((segment = block_getNext(segmentIterator)) != NULL) {
                    if (segment_getChildNumber(segment) == 0) {
                        i++;
                    }
                }
                block_destructInstanceIterator(segmentIterator);
                intListAppend(leafDegrees, i);
                intListAppend(leafCoverage, block_getLength(block) * i);
                if(perColumnStats) {
                    for (int32_t j = 0; j < block_getLength(block); j++) {
                        intListAppend(columnDegrees, block_getInstanceNumber(block));
                        intListAppend(columnLeafDegrees, i);
                    }
                }
            }
        }
        flower_destructBlockIterator(blockIterator);
        intListAppend(counts, flower_getBlockNumber(flower));
    }
}

void reportBlockStatsP(Flower *flower, FILE *fileHandle, bool(*includeBlock)(
        Block *), const char *attribString, bool perColumnStats) {
    /*
     * Prints the block stats to the XML file.
     */
    struct IntList *counts = constructEmptyIntList(0);
    struct IntList *lengths = constructEmptyIntList(0);
    struct IntList *degrees = constructEmptyIntList(0);
    struct IntList *leafDegrees = constructEmptyIntList(0);
    struct IntList *coverage = constructEmptyIntList(0);
    struct IntList *leafCoverage = constructEmptyIntList(0);
    struct IntList *columnDegrees = constructEmptyIntList(0);
    struct IntList *columnLeafDegrees = constructEmptyIntList(0);
    blockStats(flower, counts, lengths, degrees, leafDegrees, coverage,
            leafCoverage, includeBlock, columnDegrees, columnLeafDegrees, perColumnStats);
    fprintf(fileHandle, "<blocks %s>", attribString);
    tabulateAndPrintIntValues(counts, "counts", fileHandle);
    tabulateAndPrintIntValues(lengths, "lengths", fileHandle);
    tabulateAndPrintIntValues(degrees, "degrees", fileHandle);
    tabulateAndPrintIntValues(leafDegrees, "leaf_degrees", fileHandle);
    tabulateAndPrintIntValues(coverage, "coverage", fileHandle);
    tabulateAndPrintIntValues(leafCoverage, "leaf_coverage", fileHandle);
    if (perColumnStats) {
        tabulateAndPrintIntValues(columnDegrees, "column_degrees", fileHandle);
        tabulateAndPrintIntValues(columnLeafDegrees, "column_leaf_degrees",
                fileHandle);
    }
    printClosingTag("blocks", fileHandle);
    destructIntList(counts);
    destructIntList(lengths);
    destructIntList(degrees);
    destructIntList(leafDegrees);
    destructIntList(coverage);
    destructIntList(leafCoverage);
    destructIntList(columnDegrees);
    destructIntList(columnLeafDegrees);
}

int32_t reportBlockStats_minBlockDegree;
bool reportBlockStatsP2(Block *block) {
    Segment *segment;
    Block_InstanceIterator *segmentIterator = block_getInstanceIterator(block);
    int32_t i = 0;
    while ((segment = block_getNext(segmentIterator)) != NULL) {
        if (segment_getChildNumber(segment) == 0) {
            i++;
        }
    }
    block_destructInstanceIterator(segmentIterator);
    return i >= reportBlockStats_minBlockDegree;
}

static void reportBlockStats(Flower *flower, FILE *fileHandle, int32_t minBlockDegree, bool perColumnStats) {
    reportBlockStats_minBlockDegree = minBlockDegree;
    char *cA = stString_print("minimum_leaf_degree=\"%i\"", minBlockDegree);
    reportBlockStatsP(flower, fileHandle, reportBlockStatsP2, cA,
            perColumnStats);
}

static void chainStats(Flower *flower, struct IntList *counts,
        struct IntList *blockNumbers, struct IntList *baseBlockLengths,
        struct IntList *linkNumbers, struct IntList *avgInstanceBaseLengths,
        int32_t minNumberOfBlocksInChain) {
    /*
     * Gets stats on the chains.
     * Counts is numbers per non-terminal flower.
     * Block number is the number of blocks per chain.
     * Base block lengths in the number of basepairs in blocks per chain.
     * Link numbers if the number of links per chain.
     * Avg instance base lengths is the avg number of basepairs in an instance of a chain, per chain.
     */
    if (!flower_isTerminal(flower)) {
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            assert(group_getNestedFlower(group) != NULL);
            chainStats(group_getNestedFlower(group), counts, blockNumbers,
                    baseBlockLengths, linkNumbers, avgInstanceBaseLengths,
                    minNumberOfBlocksInChain);
        }
        flower_destructGroupIterator(groupIterator);

        Flower_ChainIterator *chainIterator = flower_getChainIterator(flower);
        Chain *chain;
        Block **blocks;
        int32_t i, j, k, l;
        l = 0;
        while ((chain = flower_getNextChain(chainIterator)) != NULL) {
            blocks = chain_getBlockChain(chain, &i);
            k = 0;
            for (j = 0; j < i; j++) {
                k += block_getLength(blocks[j]);
            }

            /*Chain stats are only for those containing two or more blocks.*/
            if (i >= minNumberOfBlocksInChain) {
                intListAppend(blockNumbers, i);
                intListAppend(baseBlockLengths, k);
                intListAppend(linkNumbers, chain_getLength(chain));
                intListAppend(avgInstanceBaseLengths,
                        chain_getAverageInstanceBaseLength(chain));
                l++;
            }
        }
        flower_destructBlockIterator(chainIterator);
        intListAppend(counts, l);
    }
}

static void reportChainStats(Flower *flower, int32_t minNumberOfBlocksInChain,
        FILE *fileHandle) {
    /*
     * Prints the chain stats to the XML file.
     */
    struct IntList *counts = constructEmptyIntList(0);
    struct IntList *blockNumbers = constructEmptyIntList(0);
    struct IntList *baseBlockLengths = constructEmptyIntList(0);
    struct IntList *linkNumbers = constructEmptyIntList(0);
    struct IntList *avgInstanceBaseLengths = constructEmptyIntList(0);
    chainStats(flower, counts, blockNumbers, baseBlockLengths, linkNumbers,
            avgInstanceBaseLengths, minNumberOfBlocksInChain);
    fprintf(fileHandle, "<chains minimum_number_of_blocks_in_chain=\"%i\">",
            minNumberOfBlocksInChain);
    tabulateAndPrintIntValues(counts, "counts", fileHandle);
    tabulateAndPrintIntValues(blockNumbers, "block_numbers", fileHandle);
    tabulateAndPrintIntValues(baseBlockLengths, "base_block_lengths",
            fileHandle);
    tabulateAndPrintIntValues(linkNumbers, "link_numbers", fileHandle);
    tabulateAndPrintIntValues(avgInstanceBaseLengths,
            "avg_instance_base_length", fileHandle);
    printClosingTag("chains", fileHandle);
    destructIntList(counts);
    destructIntList(blockNumbers);
    destructIntList(baseBlockLengths);
    destructIntList(linkNumbers);
    destructIntList(avgInstanceBaseLengths);
}

void terminalFlowerSizes(Flower *flower, struct IntList *sizes) {
    /*
     * Reports stats on the size of terminal flowers..
     * Sizes = This gives the sizes of the terminal flowers, i.e. the number of bases in adjacencies between ends in terminal flowers.
     * If the cactus tree has been fully decomposed then all terminal flowers will contain 0 bases.
     */
    if (flower_isTerminal(flower)) {
        intListAppend(sizes, flower_getTotalBaseLength(flower));
    } else {
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            assert(!group_isLeaf(group));
            terminalFlowerSizes(group_getNestedFlower(group), sizes);
        }
        flower_destructGroupIterator(groupIterator);
    }
}

void reportTerminalFlowerSizes(Flower *flower, FILE *fileHandle) {
    /*
     * Prints the terminal group size stats to the XML file.
     */
    struct IntList *sizes = constructEmptyIntList(0);
    terminalFlowerSizes(flower, sizes);
    tabulateAndPrintIntValues(sizes, "terminal_group_sizes", fileHandle);
    destructIntList(sizes);
}

static int32_t endDegree(End *end) {
    /*
     * Returns the number of distint ends and end is connected to.
     */
    struct List *list = constructEmptyList(0, NULL);
    End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
    Cap *cap;
    while ((cap = end_getNext(instanceIterator)) != NULL) {
        Cap *cap2 = cap_getAdjacency(cap);
        if (cap2 != NULL) {
            End *end = end_getPositiveOrientation(cap_getEnd(cap2));
            if (!listContains(list, end)) {
                listAppend(list, end);
            }
        }
    }
    end_destructInstanceIterator(instanceIterator);
    int32_t i = list->length;
    destructList(list);
    return i;
}

int32_t netStats(Flower *flower, stList *totalEndNumbersPerTerminalGroup,
        stList *totalNonFreeStubEndNumbersPerTerminalGroup,
        struct List *endDegrees, stList *totalGroupsPerNet) {
    /*
     * Calculates stats on the nets which contain tangle groups, so called 'tangle nets'
     * Reports ends per tangle net, non-free stub ends per tangle net, avg number of distinct
     * end an end is connected to in a tangle net and the number of tangle groups per net.
     */
    if (flower_isTerminal(flower)) {
        stList_append(totalEndNumbersPerTerminalGroup, stIntTuple_construct(1,
                flower_getEndNumber(flower)));
        stList_append(totalNonFreeStubEndNumbersPerTerminalGroup,
                stIntTuple_construct(1, flower_getEndNumber(flower)
                        - flower_getFreeStubEndNumber(flower)));

        End *end;
        Flower_EndIterator *flowerEndIt = flower_getEndIterator(flower);
        int32_t endConnectivity = 0;
        while ((end = flower_getNextEnd(flowerEndIt))) {
            endConnectivity += endDegree(end);
        }
        flower_destructEndIterator(flowerEndIt);
        listAppend(endDegrees, constructFloat((0.0 + endConnectivity)
                / flower_getEndNumber(flower)));
        return 1;
    } else {
        int32_t totalGroups = 0;
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            assert(group_getNestedFlower(group) != NULL);
            totalGroups += netStats(group_getNestedFlower(group),
                    totalEndNumbersPerTerminalGroup,
                    totalNonFreeStubEndNumbersPerTerminalGroup, endDegrees,
                    totalGroupsPerNet);
        }
        flower_destructGroupIterator(groupIterator);
        if (flower_getParentGroup(flower) != NULL) {
            Group *parentGroup = flower_getParentGroup(flower);
            if (group_isTangle(parentGroup)) {
                return totalGroups;
            }
        }
        stList_append(totalGroupsPerNet, stIntTuple_construct(1, totalGroups));
        return 0;
    }
}

void reportNetStats(Flower *flower, FILE *fileHandle) {
    /*
     * Prints the end stats to the XML file.
     */
    stList *totalEndNumbersPerTerminalGroup = stList_construct3(0, (void(*)(
            void *)) stIntTuple_destruct);
    stList *totalNonFreeStubEndNumbersPerTerminalGroup = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    struct List *endDegrees = constructEmptyList(0, free);
    stList *totalGroupsPerNet = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);

    netStats(flower, totalEndNumbersPerTerminalGroup,
            totalNonFreeStubEndNumbersPerTerminalGroup, endDegrees,
            totalGroupsPerNet);

    fprintf(fileHandle, "<nets>");
    tabulateAndPrintIntTupleValues(totalEndNumbersPerTerminalGroup,
            "total_end_numbers_per_terminal_group", fileHandle);
    tabulateAndPrintIntTupleValues(totalNonFreeStubEndNumbersPerTerminalGroup,
            "total_non_free_stub_end_numbers_per_terminal_group", fileHandle);
    tabulateAndPrintIntTupleValues(totalGroupsPerNet, "total_groups_per_net",
            fileHandle);
    tabulateAndPrintFloatValues(endDegrees, "end_degrees_per_terminal_group",
            fileHandle);
    printClosingTag("nets", fileHandle);

    stList_destruct(totalEndNumbersPerTerminalGroup);
    stList_destruct(totalNonFreeStubEndNumbersPerTerminalGroup);
    destructList(endDegrees);
    stList_destruct(totalGroupsPerNet);
}

void faceStats(Flower *flower, struct IntList *numberPerGroup,
        struct IntList *cardinality, struct IntList *isSimple,
        struct IntList *isRegular, struct IntList *isCanonical,
        struct IntList *facesPerFaceAssociatedEnd, int32_t includeLinkGroups,
        int32_t includeTangleGroups) {
    /*
     * Face stats for the terminal AVGs.
     * Number per group: faces per group.
     * Cardinality of face.
     * isSimple: if face is simple.
     * isRegular: is face is regular.
     * isCanonical: if face is canonical.
     * facesPerFaceAssociatedEnd: the number of faces associated with each end that
     * is associated with at least one end. Used to calculate the breakpoint reuse ratio.
     */
    if (flower_isTerminal(flower)) {
        Group *group = flower_getParentGroup(flower);
        if (group != NULL) { //Only works when parent is not empty.
            if ((includeLinkGroups && group_getLink(group) != NULL)
                    || (includeTangleGroups && group_getLink(group) == NULL)) {
                Flower_FaceIterator *faceIterator = flower_getFaceIterator(
                        flower);
                Face *face;
                while ((face = flower_getNextFace(faceIterator)) != NULL) {
                    intListAppend(cardinality, face_getCardinal(face));
                    intListAppend(isSimple, face_isSimple(face));
                    intListAppend(isRegular, face_isRegular(face));
                    intListAppend(isCanonical, face_isCanonical(face));
                }
                flower_destructFaceIterator(faceIterator);
            }
        }
    } else {
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            //Call recursively..
            assert(!group_isLeaf(group));
            faceStats(group_getNestedFlower(group), numberPerGroup,
                    cardinality, isSimple, isRegular, isCanonical,
                    facesPerFaceAssociatedEnd, includeLinkGroups,
                    includeTangleGroups);
        }
        flower_destructGroupIterator(groupIterator);
    }
}

void reportFaceStats(Flower *flower, int32_t includeLinkGroups,
        int32_t includeTangleGroups, FILE *fileHandle) {
    /*
     * Prints the reference stats to the XML file.
     */
    struct IntList *numberPerGroup = constructEmptyIntList(0);
    struct IntList *cardinality = constructEmptyIntList(0);
    struct IntList *isSimple = constructEmptyIntList(0);
    struct IntList *isRegular = constructEmptyIntList(0);
    struct IntList *isCanonical = constructEmptyIntList(0);
    struct IntList *facesPerFaceAssociatedEnd = constructEmptyIntList(0);
    faceStats(flower, numberPerGroup, cardinality, isSimple, isRegular,
            isCanonical, facesPerFaceAssociatedEnd, includeLinkGroups,
            includeTangleGroups);
    fprintf(fileHandle,
            "<faces include_link_groups=\"%i\" include_tangle_groups=\"%i\">",
            includeLinkGroups != 0, includeTangleGroups != 0);
    tabulateAndPrintIntValues(numberPerGroup, "number_per_group", fileHandle);
    tabulateAndPrintIntValues(cardinality, "cardinality", fileHandle);
    tabulateAndPrintIntValues(isSimple, "is_simple", fileHandle);
    tabulateAndPrintIntValues(isRegular, "is_regular", fileHandle);
    tabulateAndPrintIntValues(isCanonical, "is_canonical", fileHandle);
    tabulateAndPrintIntValues(facesPerFaceAssociatedEnd,
            "faces_per_face_associated_end", fileHandle);
    printClosingTag("faces", fileHandle);
    destructIntList(numberPerGroup);
    destructIntList(cardinality);
    destructIntList(isSimple);
    destructIntList(isRegular);
    destructIntList(isCanonical);
}

bool isTruePseudoAdjacency(PseudoAdjacency *pseudoAdjacency) {
    /*
     * Returns 1 iff the pseudoAdjacency is a true pseudo-adjacency (not linked by and adjacency).
     */
    End *_5End = pseudoAdjacency_get5End(pseudoAdjacency);
    End *_3End = pseudoAdjacency_get3End(pseudoAdjacency);
    Cap *cap;
    bool k = 1;
    End_InstanceIterator *instanceIterator = end_getInstanceIterator(_5End);
    while ((cap = end_getNext(instanceIterator)) != NULL) {
        Cap *adjacentCap = cap_getAdjacency(cap);
        if (adjacentCap != NULL) {
            assert(end_getOrientation(_3End));
            if (end_getPositiveOrientation(cap_getEnd(adjacentCap)) == _3End) {
                k = 0;
            }
        }
    }
    end_destructInstanceIterator(instanceIterator);
    return k;
}

void referenceStats(Flower *flower,
        stList *truePseudoAdjacencyPerTerminalGroup,
        stList *pseudoAdjacencyPerTerminalGroup) {
    /*
     * Calculates stats on the reference genome structure.
     * Stats are pretty obvious.
     */
    //Call recursively..
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        if (!group_isLeaf(group)) {
            referenceStats(group_getNestedFlower(group),
                    truePseudoAdjacencyPerTerminalGroup,
                    pseudoAdjacencyPerTerminalGroup);
        }
    }
    flower_destructGroupIterator(groupIterator);

    //Calculate stats for first reference.
    if (flower_isTerminal(flower)) { //By counting only the terminal problems we
        //Calculate the number of pseudo ends
        Reference_PseudoChromosomeIterator *pseudoChromosomeIt =
                reference_getPseudoChromosomeIterator(flower_getReference(
                        flower));
        PseudoChromosome *pseudoChromosome;
        int32_t i = 0, j = 0;
        while ((pseudoChromosome = reference_getNextPseudoChromosome(
                pseudoChromosomeIt)) != NULL) {
            PseudoChromsome_PseudoAdjacencyIterator *adjacencyIt =
                    pseudoChromosome_getPseudoAdjacencyIterator(
                            pseudoChromosome);
            PseudoAdjacency *pseudoAdjacency;

            while ((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(
                    adjacencyIt)) != NULL) {
                if (isTruePseudoAdjacency(pseudoAdjacency)) {
                    i++;
                }
                j++;
            }
            pseudoChromosome_destructPseudoAdjacencyIterator(adjacencyIt);
        }
        reference_destructPseudoChromosomeIterator(pseudoChromosomeIt);
        stList_append(truePseudoAdjacencyPerTerminalGroup,
                stIntTuple_construct(1, i));
        stList_append(pseudoAdjacencyPerTerminalGroup, stIntTuple_construct(1,
                j));
    }
}

void reportReferenceStats(Flower *flower, FILE *fileHandle) {
    /*
     * Prints the reference stats to the XML file.
     */
    if (flower_getReference(flower) != NULL) {
        stList *truePseudoAdjacencyPerTerminalGroup = stList_construct3(0,
                (void(*)(void *)) stIntTuple_destruct);
        stList *pseudoAdjacencyPerTerminalGroup = stList_construct3(0,
                (void(*)(void *)) stIntTuple_destruct);
        referenceStats(flower, truePseudoAdjacencyPerTerminalGroup,
                pseudoAdjacencyPerTerminalGroup);
        fprintf(fileHandle, "<reference method=\"default\">");
        tabulateAndPrintIntTupleValues(truePseudoAdjacencyPerTerminalGroup,
                "true_pseudo_adjacencies_per_terminal_group", fileHandle);
        tabulateAndPrintIntTupleValues(pseudoAdjacencyPerTerminalGroup,
                "pseudo_adjacencies_per_terminal_group", fileHandle);
        printClosingTag("reference", fileHandle);
        stList_destruct(truePseudoAdjacencyPerTerminalGroup);
        stList_destruct(pseudoAdjacencyPerTerminalGroup);
    }
}

int64_t reportReferenceStats2P(PseudoAdjacency *pseudoAdjacency,
        int64_t length, int64_t *blockNumber,
        stList *locationsOfTruePseudoAdjacencies) {
    End *_5End = pseudoAdjacency_get5End(pseudoAdjacency);
    End *_3End = pseudoAdjacency_get3End(pseudoAdjacency);
    Group *group = end_getGroup(_5End);
    assert(end_getGroup(_3End) == group);
    Flower *nestedFlower;
    if ((nestedFlower = group_getNestedFlower(group)) != NULL) {
        End *nested5End = flower_getEnd(nestedFlower, end_getName(_5End));
        End *nested3End = flower_getEnd(nestedFlower, end_getName(_3End));
        assert(end_isAttached(nested5End));
        assert(end_isAttached(nested3End));
        assert(nested5End != NULL);
        assert(nested3End != NULL);
        PseudoAdjacency *nestedPseudoAdjacency = end_getPseudoAdjacency(
                nested5End);
        assert(nestedPseudoAdjacency != NULL);
        PseudoChromosome *nestedPseudoChromosome =
                pseudoAdjacency_getPseudoChromosome(nestedPseudoAdjacency);
        assert(nestedPseudoChromosome != NULL);
        if (pseudoAdjacency_get5End(nestedPseudoAdjacency) == nested5End) {
            assert(pseudoAdjacency_getIndex(nestedPseudoAdjacency) == 0);
            assert(pseudoChromosome_get5End(nestedPseudoChromosome)
                    == nested5End);
            assert(pseudoChromosome_get3End(nestedPseudoChromosome)
                    == nested3End);
            for (int32_t i = 0; i < pseudoChromosome_getPseudoAdjacencyNumber(
                    nestedPseudoChromosome); i++) {
                nestedPseudoAdjacency
                        = pseudoChromosome_getPseudoAdjacencyByIndex(
                                nestedPseudoChromosome, i);
                length = reportReferenceStats2P(nestedPseudoAdjacency, length,
                        blockNumber, locationsOfTruePseudoAdjacencies);
                //Add the length of the intervening block
                End *otherEnd = pseudoAdjacency_get3End(nestedPseudoAdjacency);
                if (end_isBlockEnd(otherEnd)) {
                    if (block_getInstanceNumber(end_getBlock(otherEnd)) > 0) {
                        length += block_getLength(end_getBlock(otherEnd));
                    } else {
                        assert(block_getLength(end_getBlock(otherEnd)) == 1);
                    }
                    assert(i != pseudoChromosome_getPseudoAdjacencyNumber(
                            nestedPseudoChromosome) - 1);
                } else {
                    assert(i == pseudoChromosome_getPseudoAdjacencyNumber(
                            nestedPseudoChromosome) - 1);
                }
            }
        } else {
            assert(pseudoAdjacency_get3End(nestedPseudoAdjacency) == nested5End);
            assert(pseudoChromosome_get5End(nestedPseudoChromosome)
                    == nested3End);
            assert(pseudoAdjacency_getIndex(nestedPseudoAdjacency)
                    == pseudoChromosome_getPseudoAdjacencyNumber(
                            nestedPseudoChromosome) - 1);
            for (int32_t i = pseudoChromosome_getPseudoAdjacencyNumber(
                    nestedPseudoChromosome) - 1; i >= 0; i--) {
                nestedPseudoAdjacency
                        = pseudoChromosome_getPseudoAdjacencyByIndex(
                                nestedPseudoChromosome, i);
                length = reportReferenceStats2P(nestedPseudoAdjacency, length,
                        blockNumber, locationsOfTruePseudoAdjacencies);
                //Add the length of the intervening block
                End *otherEnd = pseudoAdjacency_get5End(nestedPseudoAdjacency);
                if (end_isBlockEnd(otherEnd)) {
                    if (block_getInstanceNumber(end_getBlock(otherEnd)) > 0) {
                        length += block_getLength(end_getBlock(otherEnd));
                    } else {
                        assert(block_getLength(end_getBlock(otherEnd)) == 1);
                    }
                    assert(i != 0);
                } else {
                    assert(i == 0);
                }
            }
        }
    } else {
        if (isTruePseudoAdjacency(pseudoAdjacency)) {
            stList_append(locationsOfTruePseudoAdjacencies,
                    stIntTuple_construct(1, length));
        }
    }
    return length;
}

int32_t calculateN50(stList *list) {
    stList *sortedList = stList_copy(list, NULL);
    stList_sort(sortedList,
            (int(*)(const void *, const void *)) stIntTuple_cmpFn);
    int64_t j = 0;
    for (int32_t i = 0; i < stList_length(sortedList); i++) {
        j += stIntTuple_getPosition(stList_get(sortedList, i), 0);
    }
    int64_t k = 0;
    int32_t n50 = 0;
    for (int32_t i = 0; i < stList_length(sortedList); i++) {
        k += stIntTuple_getPosition(stList_get(sortedList, i), 0);
        if (k >= j / 2) {
            n50 = stIntTuple_getPosition(stList_get(sortedList, i), 0);
            break;
        }
    }
    stList_destruct(sortedList);
    return n50;
}

stList *filterZeroValues(stList *values) {
    stList *filteredValues = stList_construct();
    for (int32_t i = 0; i < stList_length(values); i++) {
        stIntTuple *intTuple = stList_get(values, i);
        if (stIntTuple_getPosition(intTuple, 0) > 0) {
            stList_append(filteredValues, intTuple);
        }
    }
    return filteredValues;
}

void reportReferenceStats2(Flower *flower, FILE *fileHandle) {
    if (flower_getReference(flower) != NULL) {
        stList *contigLengths = stList_construct3(0,
                (void(*)(void *)) stIntTuple_destruct);
        Reference *reference = flower_getReference(flower);
        Reference_PseudoChromosomeIterator *pseudoChromosomeIterator =
                reference_getPseudoChromosomeIterator(reference);
        PseudoChromosome *pseudoChromosome;
        stList *topLevelPseudoChromosomeLengths = stList_construct3(0,
                (void(*)(void *)) stIntTuple_destruct);
        stList *topLevelPseudoChromosomeBlockNumbers = stList_construct3(0,
                (void(*)(void *)) stIntTuple_destruct);
        while ((pseudoChromosome = reference_getNextPseudoChromosome(
                pseudoChromosomeIterator)) != NULL) {
            //Calculate the lengths of each pseudo chromosome
            stList *locationsOfTruePseudoAdjacencies = stList_construct3(0,
                    (void(*)(void *)) stIntTuple_destruct);
            int64_t length = 0;
            //Collate the lengths of the stuff in the lower level nets.
            int64_t blockNumber = 0;
            for (int32_t i = 0; i < pseudoChromosome_getPseudoAdjacencyNumber(
                    pseudoChromosome); i++) {
                PseudoAdjacency *pseudoAdjacency =
                        pseudoChromosome_getPseudoAdjacencyByIndex(
                                pseudoChromosome, i);
                length = reportReferenceStats2P(pseudoAdjacency, length,
                        &blockNumber, locationsOfTruePseudoAdjacencies);
                //Add the length of the intervening block
                End *otherEnd = pseudoAdjacency_get3End(pseudoAdjacency);
                if (end_isBlockEnd(otherEnd)) {
                    if (block_getInstanceNumber(end_getBlock(otherEnd)) > 0) {
                        blockNumber++;
                        length += block_getLength(end_getBlock(otherEnd));
                    } else {
                        assert(block_getLength(end_getBlock(otherEnd)) == 1);
                    }
                    assert(i < pseudoChromosome_getPseudoAdjacencyNumber(
                            pseudoChromosome) - 1);
                } else {
                    assert(i == pseudoChromosome_getPseudoAdjacencyNumber(
                            pseudoChromosome) - 1);
                }
            }
            stList_append(topLevelPseudoChromosomeLengths,
                    stIntTuple_construct(1, length));
            stList_append(topLevelPseudoChromosomeBlockNumbers,
                    stIntTuple_construct(1, blockNumber));
            //Now convert locations of true pseudo adjacencies into lengths of contigs
            int32_t k = 0;
            for (int32_t i = 0; i < stList_length(
                    locationsOfTruePseudoAdjacencies); i++) {
                int32_t j = stIntTuple_getPosition(stList_get(
                        locationsOfTruePseudoAdjacencies, i), 0);
                stList_append(contigLengths, stIntTuple_construct(1, j - k));
                k = j;
            }
            assert(k <= length);
            stList_append(contigLengths, stIntTuple_construct(1, length - k)); //Account for the last contig
            stList_destruct(locationsOfTruePseudoAdjacencies);
        }
        reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);

        //Get values filtered of zero length connections
        stList *contigLengthsFiltered = filterZeroValues(contigLengths);

        //Now report the results


        //Calculate N50
        int32_t n50 = calculateN50(contigLengths);
        int32_t n50Filtered = calculateN50(contigLengthsFiltered);

        fprintf(
                fileHandle,
                "<reference2 method=\"default\" n50=\"%i\" n50Filtered=\"%i\">",
                n50, n50Filtered);

        tabulateAndPrintIntTupleValues(topLevelPseudoChromosomeLengths,
                "top_level_pseudo_chromosome_lengths", fileHandle);
        tabulateAndPrintIntTupleValues(topLevelPseudoChromosomeBlockNumbers,
                "top_level_pseudo_chromosome_block_numbers", fileHandle);
        tabulateAndPrintIntTupleValues(contigLengths, "contig_lengths",
                fileHandle);
        tabulateAndPrintIntTupleValues(contigLengthsFiltered,
                "contig_lengths_filtered", fileHandle);

        printClosingTag("reference2", fileHandle);
        //Clean up
        stList_destruct(topLevelPseudoChromosomeLengths);
        stList_destruct(topLevelPseudoChromosomeBlockNumbers);
        stList_destruct(contigLengths);
        stList_destruct(contigLengthsFiltered);
    }
}

void reportCactusDiskStats(char *cactusDiskName, Flower *flower,
        FILE *fileHandle, bool perColumnStats,
        stSortedSet *includeSpecies, stSortedSet *excludeSpecies) {

    double totalSeqSize = flower_getTotalBaseLength(flower);
    fprintf(
            fileHandle,
            "<stats flower_disk=\"%s\" flower_name=\"%s\" total_sequence_length=\"%f\" >",
            cactusDiskName, cactusMisc_nameToStringStatic(
                    flower_getName(flower)), totalSeqSize);

    /*
     * Relative entropy numbers on the balance of the tree.
     */
    reportRelativeEntopyStats(flower, fileHandle);

    /*
     * Numbers on the structure of the tree.
     */
    reportFlowerStats(flower, fileHandle);

    /*
     * Numbers on the blocks.
     */
    reportBlockStats(flower, fileHandle, 0, perColumnStats);
    reportBlockStats(flower, fileHandle, 2, perColumnStats);

    /*
     * Chain statistics.
     */
    reportChainStats(flower, 0, fileHandle);
    reportChainStats(flower, 2, fileHandle);

    /*
     * Stats on terminal flowers in the tree.
     */
    reportTerminalFlowerSizes(flower, fileHandle);

    /*
     * Stats on the ends in the problem. Currently just the numbers of ends in each net.
     */
    reportNetStats(flower, fileHandle);

    /*
     * Stats on faces in the reconstruction..
     */
    reportFaceStats(flower, 1, 1, fileHandle);
    reportFaceStats(flower, 0, 1, fileHandle);
    reportFaceStats(flower, 1, 0, fileHandle);

    /*
     * Stats on the reference in the reconstruction..
     */
    reportReferenceStats(flower, fileHandle);
    reportReferenceStats2(flower, fileHandle);

    printClosingTag("stats", fileHandle);
}

