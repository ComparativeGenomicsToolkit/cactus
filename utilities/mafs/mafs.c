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
//#include "cactus_addReferenceSeq.h"


/*
 * Library for generating mafs from cactus.
 */

static char *formatSequenceHeader(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) * (1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

static void getMAFBlockP2(Segment *segment, FILE *fileHandle) {
    assert(segment != NULL);
    Sequence *sequence = segment_getSequence(segment);
    if (sequence != NULL) {
        char *sequenceHeader = formatSequenceHeader(sequence);
        int32_t start;
        if (segment_getStrand(segment)) {
            start = segment_getStart(segment) - sequence_getStart(sequence);
        } else { //start with respect to the start of the reverse complement sequence
            start = (sequence_getStart(sequence) + sequence_getLength(sequence)
                    - 1) - segment_getStart(segment);
        }
        int32_t length = segment_getLength(segment);
        char *strand = segment_getStrand(segment) ? "+" : "-";
        int32_t sequenceLength = sequence_getLength(sequence);
        char *instanceString = segment_getString(segment);
        fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n", sequenceHeader,
                start, length, strand, sequenceLength, instanceString);
        free(instanceString);
        free(sequenceHeader);
    }
}

static void getMAFBlockP(Segment *segment, FILE *fileHandle) {
    int32_t i;
    for (i = 0; i < segment_getChildNumber(segment); i++) {
        getMAFBlockP(segment_getChild(segment, i), fileHandle);
    }
    getMAFBlockP2(segment, fileHandle);
}

static int32_t getNumberOnPositiveStrand(Block *block) {
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment;
    int32_t i = 0;
    while ((segment = block_getNext(it)) != NULL) {
        if (segment_getChildNumber(segment) == 0) {
            if (segment_getStrand(segment)) {
                i++;
            }
        }
    }
    block_destructInstanceIterator(it);
    return i;
}

void getMAFBlock(Block *block, FILE *fileHandle) {
//void getMAFBlock(Block *block, FILE *fileHandle, ReferenceSequence *referenceSequence) {
    /*
     * Outputs a MAF representation of the block to the given file handle.
     */
    //Correct the orientation..
    if (getNumberOnPositiveStrand(block) == 0) {
        block = block_getReverse(block);
    }
    if (block_getInstanceNumber(block) > 0) {
        //Add in the header
        if (block_getRootInstance(block) != NULL) {
            /* Get newick tree string with internal labels and no unary events */
            char *newickTreeString = block_makeNewickString(block, 1, 0);
            assert(newickTreeString != NULL);
            fprintf(fileHandle, "a score=%i tree='%s'\n",
                    block_getLength(block) * block_getInstanceNumber(block),
                    newickTreeString);
            free(newickTreeString);
        } else {
            fprintf(fileHandle, "a score=%i\n", block_getLength(block)
                    * block_getInstanceNumber(block));
        }
        //Now for the reference segment
        /*if (referenceSequence != NULL) {
            char *instanceString = getConsensusString(block);
            fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n",
                    referenceSequence->header, referenceSequence->index,
                    block_getLength(block), "+", referenceSequence->length,
                    instanceString);
            free(instanceString);
            referenceSequence->index += block_getLength(block);
        }*/
        //Now add the blocks in
        if (block_getRootInstance(block) != NULL) {
            assert(block_getRootInstance(block) != NULL);
            getMAFBlockP(block_getRootInstance(block), fileHandle);
            fprintf(fileHandle, "\n");
        } else {
            Block_InstanceIterator *iterator = block_getInstanceIterator(block);
            Segment *segment;
            while ((segment = block_getNext(iterator)) != NULL) {
                getMAFBlockP2(segment, fileHandle);
            }
            block_destructInstanceIterator(iterator);
            fprintf(fileHandle, "\n");
        }
    }
}

//void getMAFSReferenceOrdered_walkDown(End *end, FILE *fileHandle, ReferenceSequence *referenceSequence);
static void getMAFSReferenceOrdered_walkDown(End *end, FILE *fileHandle, void (*getMafBlockFn)(Block *, FILE *));

//void getMAFSReferenceOrdered_walkUp(End *end, FILE *fileHandle, ReferenceSequence *referenceSequence) {
static void getMAFSReferenceOrdered_walkUp(End *end, FILE *fileHandle, void (*getMafBlockFn)(Block *, FILE *)) {
    assert(end != NULL);
    if (end_isBlockEnd(end)) {
        getMafBlockFn(end_getBlock(end), fileHandle);
        getMAFSReferenceOrdered_walkDown(end_getOtherBlockEnd(end), fileHandle, getMafBlockFn);
        //getMAFBlock(end_getBlock(end), fileHandle, referenceSequence);
        //getMAFSReferenceOrdered_walkDown(end_getOtherBlockEnd(end), fileHandle, referenceSequence);
    } else {
        assert(end_isAttached(end));
        Group *parentGroup = flower_getParentGroup(end_getFlower(end));
        if (parentGroup != NULL) {
            //getMAFSReferenceOrdered_walkUp(group_getEnd(parentGroup,
            //        end_getName(end)), fileHandle, referenceSequence);
            getMAFSReferenceOrdered_walkUp(group_getEnd(parentGroup,
                    end_getName(end)), fileHandle, getMafBlockFn);
        } else { //We reached the end of a pseudo-chromosome!
            assert(pseudoChromosome_get3End(pseudoAdjacency_getPseudoChromosome(end_getPseudoAdjacency(end))) == end);
        }
    }
}

//void getMAFSReferenceOrdered_walkDown(End *end, FILE *fileHandle, ReferenceSequence *referenceSequence) {
static void getMAFSReferenceOrdered_walkDown(End *end, FILE *fileHandle, void (*getMafBlockFn)(Block *, FILE *)) {
    assert(end != NULL);
    //assert(end_isAttached(end));
    Group *group = end_getGroup(end);
    if (group_isLeaf(group)) { //Walk across
        PseudoAdjacency *pseudoAdjacency = end_getPseudoAdjacency(end);
        assert(pseudoAdjacency != NULL);
        assert(pseudoAdjacency_get3End(pseudoAdjacency) == end || pseudoAdjacency_get5End(pseudoAdjacency) == end);
        if (pseudoAdjacency_get3End(pseudoAdjacency) == end) {
            end = pseudoAdjacency_get5End(pseudoAdjacency);
        } else {
            end = pseudoAdjacency_get3End(pseudoAdjacency);
        }
        //Now walk up
        //getMAFSReferenceOrdered_walkUp(end, fileHandle, referenceSequence);
        getMAFSReferenceOrdered_walkUp(end, fileHandle, getMafBlockFn);
    } else { //Walk down
        //getMAFSReferenceOrdered_walkDown(flower_getEnd(group_getNestedFlower(
        //        group), end_getName(end)), fileHandle, referenceSequence);
        getMAFSReferenceOrdered_walkDown(flower_getEnd(group_getNestedFlower(
                group), end_getName(end)), fileHandle, getMafBlockFn);
    }
}

//void getMAFsReferenceOrdered(Flower *flower, FILE *fileHandle, ReferenceSequence *referenceSequence) {
void getMAFsReferenceOrdered(Flower *flower, FILE *fileHandle, void (*getMafBlockFn)(Block *, FILE *)) {
    /*
     * Outputs MAF representations of all the block in the flower and its descendants, ordered
     * according to the reference ordering.
     */
    Reference *reference = flower_getReference(flower);
    assert(reference != NULL);
    Reference_PseudoChromosomeIterator *it =
            reference_getPseudoChromosomeIterator(reference);
    PseudoChromosome *pseudoChromosome;
    while ((pseudoChromosome = reference_getNextPseudoChromosome(it)) != NULL) {
        End *end = pseudoChromosome_get5End(pseudoChromosome);
        assert(!end_isBlockEnd(end));
        //getMAFSReferenceOrdered_walkDown(end, fileHandle, referenceSequence);
        getMAFSReferenceOrdered_walkDown(end, fileHandle, getMafBlockFn);
    }
    reference_destructPseudoChromosomeIterator(it);
}

void getMAFs(Flower *flower, FILE *fileHandle, void (*getMafBlock)(Block *, FILE *)) {
    /*
     * Outputs MAF representations of all the block sin the flower and its descendants.
     */

    //Make MAF blocks for each block
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIterator)) != NULL) {
        getMAFBlock(block, fileHandle);
        //getMAFBlock(block, fileHandle, NULL);
    }
    flower_destructBlockIterator(blockIterator);

    //Call child flowers recursively.
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        if (!group_isLeaf(group)) {
            getMAFs(group_getNestedFlower(group), fileHandle, getMafBlock); //recursive call.
        }
    }
    flower_destructGroupIterator(groupIterator);
}

void makeMAFHeader(Flower *flower, FILE *fileHandle) {
    fprintf(fileHandle, "##maf version=1 scoring=N/A\n");
    char *cA = eventTree_makeNewickString(flower_getEventTree(flower));
    fprintf(fileHandle, "# cactus %s\n\n", cA);
    free(cA);
}
