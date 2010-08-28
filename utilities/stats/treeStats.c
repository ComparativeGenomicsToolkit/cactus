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
        *min = INT32_MAX;
        *max = INT32_MAX;
        *avg = INT32_MAX;
        *median = INT32_MAX;
        return;
    }
    unsortedValues = listCopy(unsortedValues); //copy the input list, to avoid altering the input.
    assert(unsortedValues->length> 0);
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
        *min = INT32_MAX;
        *max = INT32_MAX;
        *avg = INT32_MAX;
        *median = INT32_MAX;
        return;
    }
    unsortedValues = intListCopy(unsortedValues);
    assert(unsortedValues->length> 0);
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
    double followingPathBitScore = (log(flower_getGroupNumber(flower)) / log(2.0))
            + pathBitScore;
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

static void flowerStats(Flower *flower, int32_t currentDepth, struct IntList *children,
        struct IntList *tangleChildren, struct IntList *linkChildren,
        struct IntList *depths) {
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
            flowerStats(group_getNestedFlower(group), currentDepth + 1, children,
                    tangleChildren, linkChildren, depths);
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

void blockStats(Flower *flower, struct IntList *counts, struct IntList *lengths,
        struct IntList *degrees, struct IntList *leafDegrees,
        struct IntList *coverage, struct IntList *leafCoverage,
        int32_t minLeafDegree) {
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
                    leafDegrees, coverage, leafCoverage, minLeafDegree);
        }
        flower_destructGroupIterator(groupIterator);
        Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
        Block *block;
        while ((block = flower_getNextBlock(blockIterator)) != NULL) {
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
            if (i >= minLeafDegree) {
                intListAppend(lengths, block_getLength(block));
                intListAppend(degrees, block_getInstanceNumber(block));
                intListAppend(coverage, block_getLength(block)
                        * block_getInstanceNumber(block));
                intListAppend(leafDegrees, i);
                intListAppend(leafCoverage, block_getLength(block) * i);
            }
        }
        flower_destructBlockIterator(blockIterator);
        intListAppend(counts, flower_getBlockNumber(flower));
    }
}

void reportBlockStats(Flower *flower, FILE *fileHandle, int32_t minLeafDegree) {
    /*
     * Prints the block stats to the XML file.
     */
    struct IntList *counts = constructEmptyIntList(0);
    struct IntList *lengths = constructEmptyIntList(0);
    struct IntList *degrees = constructEmptyIntList(0);
    struct IntList *leafDegrees = constructEmptyIntList(0);
    struct IntList *coverage = constructEmptyIntList(0);
    struct IntList *leafCoverage = constructEmptyIntList(0);
    blockStats(flower, counts, lengths, degrees, leafDegrees, coverage,
            leafCoverage, minLeafDegree);
    fprintf(fileHandle, "<blocks minimum_leaf_degree=\"%i\">", minLeafDegree);
    tabulateAndPrintIntValues(counts, "counts", fileHandle);
    tabulateAndPrintIntValues(lengths, "lengths", fileHandle);
    tabulateAndPrintIntValues(degrees, "degrees", fileHandle);
    tabulateAndPrintIntValues(leafDegrees, "leaf_degrees", fileHandle);
    tabulateAndPrintIntValues(coverage, "coverage", fileHandle);
    tabulateAndPrintIntValues(leafCoverage, "leaf_coverage", fileHandle);
    printClosingTag("blocks", fileHandle);
    destructIntList(counts);
    destructIntList(lengths);
    destructIntList(degrees);
    destructIntList(leafDegrees);
    destructIntList(coverage);
    destructIntList(leafCoverage);
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

void reportChainStats(Flower *flower, int32_t minNumberOfBlocksInChain,
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

void netStats(Flower *flower, struct IntList *totalEndNumbersPerNet,
        struct IntList *totalNonFreeStubEndNumbersPerNet,
        struct List *endDegrees,
        struct IntList *totalGroupsPerNet) {
    /*
     * Calculates stats on the nets which contain tangle groups, so called 'tangle nets'
     * Reports ends per tangle net, non-free stub ends per tangle net, avg number of distinct
     * end an end is connected to in a tangle net and the number of tangle groups per net.
     */
    if (!flower_isTerminal(flower)) { //Do not double count terminal groups when doing the math.
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        int32_t tangleGroupEndSum = 0;
        int32_t tangleGroupEndSumExcludingFreeStubs = 0;
        int32_t endConnectivity = 0;
        int32_t tangleGroups = 0;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            if(group_isTangle(group)) {
                tangleGroups++;
                tangleGroupEndSum += group_getEndNumber(group);
                tangleGroupEndSumExcludingFreeStubs += group_getAttachedStubEndNumber(group) + group_getBlockEndNumber(group);
                End *end;
                Group_EndIterator *groupEndIt = group_getEndIterator(group);
                while((end = group_getNextEnd(groupEndIt))) {
                    endConnectivity += endDegree(end);
                }
                group_destructEndIterator(groupEndIt);
            }
            assert(!group_isLeaf(group));
            netStats(group_getNestedFlower(group), totalEndNumbersPerNet, totalNonFreeStubEndNumbersPerNet, endDegrees, totalGroupsPerNet);
        }
        flower_destructGroupIterator(groupIterator);
        if(tangleGroups > 0) {
            intListAppend(totalEndNumbersPerNet, tangleGroupEndSum);
            intListAppend(totalNonFreeStubEndNumbersPerNet, tangleGroupEndSumExcludingFreeStubs);
            listAppend(endDegrees, constructFloat((0.0 + endConnectivity)/tangleGroupEndSum));
            intListAppend(totalGroupsPerNet, tangleGroups);
        }
    }
}

void reportNetStats(Flower *flower, FILE *fileHandle) {
    /*
     * Prints the end stats to the XML file.
     */
    struct IntList *totalEndNumbersPerNet = constructEmptyIntList(0);
    struct IntList *totalNonFreeStubEndNumbersPerNet = constructEmptyIntList(0);
    struct List *endDegrees = constructEmptyList(0, (void (*)(void *))destructFloat);
    struct IntList *totalGroupsPerNet = constructEmptyIntList(0);
    netStats(flower, totalEndNumbersPerNet, totalNonFreeStubEndNumbersPerNet, endDegrees, totalGroupsPerNet);
    fprintf(
            fileHandle,
            "<nets>");
    tabulateAndPrintIntValues(totalEndNumbersPerNet, "end_numbers_per_net", fileHandle);
    tabulateAndPrintIntValues(totalNonFreeStubEndNumbersPerNet, "non_free_stub_end_numbers_per_net", fileHandle);
    tabulateAndPrintFloatValues(endDegrees, "end_degrees_per_net", fileHandle);
    tabulateAndPrintIntValues(totalGroupsPerNet, "total_groups_per_net", fileHandle);
    printClosingTag("nets", fileHandle);
    destructIntList(totalEndNumbersPerNet);
    destructIntList(totalNonFreeStubEndNumbersPerNet);
    destructList(endDegrees);
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
        assert(group != NULL);
        if ((includeLinkGroups && group_getLink(group) != NULL)
                || (includeTangleGroups && group_getLink(group) == NULL)) {
            Flower_FaceIterator *faceIterator = flower_getFaceIterator(flower);
            Face *face;
            while ((face = flower_getNextFace(faceIterator)) != NULL) {
                intListAppend(cardinality, face_getCardinal(face));
                intListAppend(isSimple, face_isSimple(face));
                intListAppend(isRegular, face_isRegular(face));
                intListAppend(isCanonical, face_isCanonical(face));
            }
            flower_destructFaceIterator(faceIterator);
        }
    } else {
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            //Call recursively..
            assert(!group_isLeaf(group));
            faceStats(group_getNestedFlower(group), numberPerGroup, cardinality,
                    isSimple, isRegular, isCanonical,
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

void referenceStats(Flower *flower, struct IntList *pseudoChromosomeNumber,
        struct IntList *pseudoAdjacencyNumberPerChromosome,
        struct IntList *truePseudoAdjacencyNumberPerChromosome,
        struct IntList *truePseudoEndsPerFlower,
        struct IntList *linksPerChromosome) {
    /*
     * Calculates stats on the reference genome structure.
     * Stats are pretty obvious.
     */
    //Call recursively..
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        if (!group_isLeaf(group)) {
            referenceStats(group_getNestedFlower(group), pseudoChromosomeNumber,
                    pseudoAdjacencyNumberPerChromosome,
                    truePseudoAdjacencyNumberPerChromosome, linksPerChromosome, truePseudoEndsPerFlower);
        }
    }
    flower_destructGroupIterator(groupIterator);

    //Calculate stats for first reference.
    if (flower_isTerminal(flower)) { //By counting only the terminal problems we
        //include each end only once in our statistics
        Reference *reference = flower_getReference(flower);
        assert(reference != NULL);
        Reference_PseudoChromosomeIterator *pseudoChromosomeIterator =
                reference_getPseudoChromosomeIterator(reference);
        PseudoChromosome *pseudoChromosome;
        intListAppend(pseudoChromosomeNumber,
                reference_getPseudoChromosomeNumber(reference));
        while ((pseudoChromosome = reference_getNextPseudoChromosome(
                pseudoChromosomeIterator)) != NULL) {
            intListAppend(pseudoAdjacencyNumberPerChromosome,
                    pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome));
            PseudoChromsome_PseudoAdjacencyIterator *adjacencyIterator =
                    pseudoChromosome_getPseudoAdjacencyIterator(
                            pseudoChromosome);
            PseudoAdjacency *pseudoAdjacency;
            int32_t i = 0, j = 0;
            while ((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(
                    adjacencyIterator)) != NULL) {
                Group *group = end_getGroup(pseudoAdjacency_get5End(
                        pseudoAdjacency));
                assert(group != NULL);
                if (group_getLink(group) != NULL) {
                    i++;
                }
                End *_5End = pseudoAdjacency_get5End(pseudoAdjacency);
                End *_3End = pseudoAdjacency_get3End(pseudoAdjacency);
                Cap *cap;
                int32_t k = 1;
                End_InstanceIterator *instanceIterator =
                        end_getInstanceIterator(_5End);
                while ((cap = end_getNext(instanceIterator)) != NULL) {
                    Cap *adjacentCap = cap_getAdjacency(cap);
                    if (adjacentCap != NULL) {
                        assert(end_getOrientation(_3End));
                        if (end_getPositiveOrientation(cap_getEnd(adjacentCap))
                                == _3End) {
                            k = 0;
                        }
                    }
                }
                end_destructInstanceIterator(instanceIterator);
                if (k) {
                    j++;
                }
            }
            intListAppend(linksPerChromosome, i);
            intListAppend(truePseudoAdjacencyNumberPerChromosome, j);
        }
        reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
        //Calculate the number of pseudo ends
        Flower_EndIterator *endIt = flower_getEndIterator(flower);
        End *end;
        int32_t i = 0;
        while((end = flower_getNextEnd(endIt)) != NULL) {
            if(end_getInstanceNumber(end) == 0) {
                i++;
            }
        }
        flower_destructEndIterator(endIt);
        intListAppend(truePseudoEndsPerFlower, i);
    }
}

void reportReferenceStats(Flower *flower, FILE *fileHandle) {
    /*
     * Prints the reference stats to the XML file.
     */
    if (flower_getReference(flower) != NULL) {
        struct IntList *pseudoChromosomeNumber = constructEmptyIntList(0);
        struct IntList *pseudoAdjacencyNumberPerChromosome =
                constructEmptyIntList(0);
        struct IntList *truePseudoAdjacencyNumberPerChromosome =
                constructEmptyIntList(0);
        struct IntList *linksPerChromosome = constructEmptyIntList(0);
        struct IntList *truePseudoEndsPerFlower = constructEmptyIntList(0);
        referenceStats(flower, pseudoChromosomeNumber,
                pseudoAdjacencyNumberPerChromosome,
                truePseudoAdjacencyNumberPerChromosome, linksPerChromosome, truePseudoEndsPerFlower);
        fprintf(fileHandle, "<reference method=\"default\">");
        tabulateAndPrintIntValues(pseudoChromosomeNumber,
                "pseudo_chromosome_number", fileHandle);
        tabulateAndPrintIntValues(pseudoAdjacencyNumberPerChromosome,
                "pseudo_adjacency_number_per_pseudo_chromosome", fileHandle);
        tabulateAndPrintIntValues(truePseudoAdjacencyNumberPerChromosome,
                "true_pseudo_adjacency_number_per_pseudo_chromosome",
                fileHandle);
        tabulateAndPrintIntValues(linksPerChromosome, "links_per_chromosome",
                fileHandle);
        tabulateAndPrintIntValues(truePseudoEndsPerFlower, "true_pseudo_ends_per_flower",
                        fileHandle);
        printClosingTag("reference", fileHandle);
        destructIntList(pseudoChromosomeNumber);
        destructIntList(pseudoAdjacencyNumberPerChromosome);
        destructIntList(truePseudoAdjacencyNumberPerChromosome);
        destructIntList(linksPerChromosome);
        destructIntList(truePseudoEndsPerFlower);
    }
}

void reportCactusDiskStats(char *cactusDiskName, Flower *flower, FILE *fileHandle) {

    double totalSeqSize = flower_getTotalBaseLength(flower);
    fprintf(
            fileHandle,
            "<stats flower_disk=\"%s\" flower_name=\"%s\" total_sequence_length=\"%f\" >",
            cactusDiskName, cactusMisc_nameToStringStatic(flower_getName(flower)),
            totalSeqSize);

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
    reportBlockStats(flower, fileHandle, 0);
    reportBlockStats(flower, fileHandle, 2);

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

    printClosingTag("stats", fileHandle);
}

