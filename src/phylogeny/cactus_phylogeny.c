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
#include <sys/stat.h>
#include <time.h>
#include <getopt.h>

#include "hashTableC.h"
#include "bioioC.h"
#include "commonC.h"
#include "cactus.h"
#include "avl.h"
#include "phylogeny.h"

#include "treelib.h"

#define closeEnough 0.001

typedef struct _chainAlignment {
    /*
     * Structure to represent a concatenated list of blocks as a single 2d alignment.
     */
    //Matrix of alignment, matrix[column][row]
    Segment ***matrix; //NULL values are okay, where an instance of an block is missing from an instance of the chain.
    Block **blocks; //this list of blocks, in order.
    int64_t columnNumber; //the number of blocks.
    int64_t rowNumber; //the number of rows of the alignment, each row containing an instance of blocks in the chain.
    int64_t totalAlignmentLength; //the length in base pairs of the alignment.
} ChainAlignment;

char **chainAlignment_getAlignment(ChainAlignment *chainAlignment) {
    /*
     * Constructs a concatenated 2d matrix of chars referenced by [block][instance], representing the
     * base pair alignment of the chain alignment.
     */
    char **alignment;
    int64_t i, j, k, l;
    Segment *segment;
    char *cA;

    //alloc the memory for the char alignment.
    alignment = st_malloc(sizeof(void *) * chainAlignment->rowNumber);
    for (i = 0; i < chainAlignment->rowNumber; i++) {
        alignment[i] = st_malloc(sizeof(char)
                * (chainAlignment->totalAlignmentLength + 1));
    }

    /*
     * Fills in the alignment.
     */
    for (j = 0; j < chainAlignment->rowNumber; j++) {
        l = 0;
        for (i = 0; i < chainAlignment->columnNumber; i++) {
            segment = chainAlignment->matrix[i][j];
            if (segment == NULL) {
                for (k = 0; k < block_getLength(chainAlignment->blocks[i]); k++) {
                    alignment[j][l++] = 'N';
                }
            } else {
                cA = segment_getString(segment);
                for (k = 0; k < segment_getLength(segment); k++) {
                    alignment[j][l++] = cA[k];
                }
                free(cA);
            }
        }
        alignment[j][l] = '\0';
        assert(l == chainAlignment->totalAlignmentLength);
    }

    return alignment;
}

Segment *chainAlignment_getFirstNonNullSegment(ChainAlignment *chainAlignment,
        int64_t row, bool increasing) {
    /*
     * Gets the first instance on an row in the chain alignment which is non null. If increasing is false, gets the last.
     */
    Segment *segment;
    int64_t j = 0, k = chainAlignment->columnNumber, l = 1;
    if (!increasing) {
        j = chainAlignment->columnNumber - 1;
        k = -1;
        l = -1;
    }
    for (; j != k; j += l) {
        segment = chainAlignment->matrix[j][row];
        if (segment != NULL) {
            return segment;
        }
    }
    assert(0);
    return NULL;
}

Name *chainAlignment_getEndNames(ChainAlignment *chainAlignment, bool _5End) {
    /*
     * Gets the names of the ends associated with the instances at the ends of a chain alignment.
     */
    Name *names;
    int64_t i;
    Segment *segment;

    names = st_malloc(sizeof(Name) * chainAlignment->rowNumber);
    for (i = 0; i < chainAlignment->rowNumber; i++) {
        segment = chainAlignment_getFirstNonNullSegment(chainAlignment, i,
                _5End);
        names[i] = end_getName(cap_getEnd(cap_getAdjacency(
                _5End ? segment_get5Cap(segment) : segment_get3Cap(segment))));
    }
    return names;
}

Name *chainAlignment_getLeafEventNames(ChainAlignment *chainAlignment) {
    /*
     * Gets the names of the leaf events of the rows (instances), in the alignment.
     */
    Name *names;
    int64_t i;
    Segment *segment;

    names = st_malloc(sizeof(Name) * chainAlignment->rowNumber);
    for (i = 0; i < chainAlignment->rowNumber; i++) {
        segment = chainAlignment_getFirstNonNullSegment(chainAlignment, i, 1);
        names[i] = event_getName(segment_getEvent(segment));
    }
    return names;
}

int64_t *chainAlignment_getBlockBoundaries(ChainAlignment *chainAlignment) {
    /*
     * Gets the boundaries of block in the chain alignment.
     */
    int64_t *blockBoundaries;
    int64_t i, j;

    blockBoundaries = st_malloc(sizeof(int64_t) * chainAlignment->columnNumber);
    j = 0;
    for (i = 0; i < chainAlignment->columnNumber; i++) {
        blockBoundaries[i] = j + block_getLength(chainAlignment->blocks[i]);
        j = blockBoundaries[i];
    }

    return blockBoundaries;
}

void chomp(const char *s) {
    char *p;
    while (NULL != s && NULL != (p = strrchr(s, '\n'))) {
        *p = '\0';
    }
}

void buildChainTrees_Bernard(int64_t blockNumber, char ***concatenatedBlocks,
        Name **_5Ends, Name **_3Ends, Name **leafEventLabels,
        int64_t **blockBoundaries, char *eventTreeString,
        ChainAlignment **chainAlignments, char **modifiedEventTreeString,
        char ****blockTreeStrings, int64_t ***refinedBlockBoundaries,
        int64_t **refinedBlockNumbers) {
    /*
     * Here's the function you need to fill in, I haven't defined the outputs yet - you get it working and then we can discuss.
     *
     * Arrays are all indexed first by the chain/concatenated blocks.
     * Alignments are column/row indexed with rows as chain instances and columns as aligned bases.
     * So for example: concatenatedBlocks[i][j][k] is the ith chain/concatenated block, jth row (chain instance), kth column (position in concatenated block)
     *
     * Names can be converted to strings with: flowerMisc_nameToString() and flowerMisc_nameToStringStatic() (See the API).
     *
     * The output arrays (last three arguments) must be initialised and are *pointers to* the point containing the array/string. Which must be
     * created.
     */

    st_logInfo("Started Bernard's function\n");

    int64_t i = 0;
    int64_t j = 0;

    int64_t rowNumber = 0;
    int64_t colNumber = 0;

    char ***blockTreeArray = NULL;
    blockTreeArray = st_malloc(sizeof(void *) * blockNumber);
    for (i = 0; i < blockNumber; i++) {
        colNumber = chainAlignments[i]->columnNumber;
        blockTreeArray[i] = st_malloc(sizeof(void *) * colNumber);
    }

    char buffer[128];

    for (i = 0; i < blockNumber; i++) {

        sprintf(buffer, "\tBN start: %" PRIi64 " of %" PRIi64 "\n", i, blockNumber);
        st_logInfo(buffer);

        rowNumber = chainAlignments[i]->rowNumber;
        colNumber = chainAlignments[i]->columnNumber;

        char *treestring = NULL;
        treestring = msa2tree(concatenatedBlocks[i], rowNumber);

        for (j = 0; j < colNumber; j++) {
            blockTreeArray[i][j] = treestring;
        }
        st_logInfo("\tEnd the loop\n");
    }
    *blockTreeStrings = blockTreeArray;

    /* Clone Block boundaries to refined Block Boundaries*/
    int64_t **newBlockBoundaries = NULL;
    int64_t *newBlockNumbers = NULL;

    newBlockBoundaries = st_malloc(sizeof(void *) * blockNumber);
    newBlockNumbers = st_malloc(sizeof(int64_t) * blockNumber);
    for (i = 0; i < blockNumber; i++) {
        colNumber = chainAlignments[i]->columnNumber;
        newBlockBoundaries[i] = st_malloc(sizeof(int64_t) * colNumber);
        newBlockNumbers[i] = colNumber;
        for (j = 0; j < colNumber; j++) {
            newBlockBoundaries[i][j] = blockBoundaries[i][j];
        }
    }
    *refinedBlockBoundaries = newBlockBoundaries;
    *refinedBlockNumbers = newBlockNumbers;

    st_logInfo("Ended Bernard's function\n");

    return;
}

Segment *getSegment(struct BinaryTree *binaryTree, Segment **segments) {
    /*
     * Gets associated segment from integer index of leaf name.
     */
    int64_t i;
    int j = sscanf(binaryTree->label, "%" PRIi64 "", &i);
    (void)j;
    assert(j == 1);
    assert(i >= 0);
    return segments[i];
}

Event *getEvent(struct BinaryTree *binaryTree, Segment **segments) {
    Segment *segment = getSegment(binaryTree, segments);
    return segment != NULL ? segment_getEvent(segment) : NULL;
}

Segment *buildChainTrees3P(Block *block, Segment **segments,
        int64_t blockNumber, struct BinaryTree *binaryTree) {
    /*
     * Recursive partner to buildChainTree3 function, recurses on the binary tree constructing the block tree.
     * The labels of the leaves are indexes into the segments array, the internal node's labels are events in the event tree.
     */
    if (binaryTree->internal) { //deal with an internal node of the block tree.
        Segment *leftInstance = buildChainTrees3P(block, segments, blockNumber,
                binaryTree->left);
        Segment *rightInstance = buildChainTrees3P(block, segments,
                blockNumber, binaryTree->right);
        if (leftInstance != NULL) {
            if (rightInstance != NULL) {
                Event *event = eventTree_getEvent(flower_getEventTree(
                        block_getFlower(block)), cactusMisc_stringToName(
                        binaryTree->label));
                assert(event != NULL); //check event is present in the event tree.
                //Check that this does not create a cycle with respect to the event tree.
                assert(event_isAncestor(segment_getEvent(leftInstance), event));
                assert(event_isAncestor(segment_getEvent(rightInstance), event));

                Segment *segment = segment_construct(block, event);
                segment_makeParentAndChild(segment, leftInstance);
                segment_makeParentAndChild(segment, rightInstance);
                return segment;
            }
            return leftInstance;
        }
        return rightInstance;
    } else { //a leaf, so find the leaf instance in the list of segments (this may be null if missing data).
        return getSegment(binaryTree, segments);
    }
}

void buildChainTrees3(Block *block, Segment **segments, int64_t blockNumber,
        struct BinaryTree *binaryTree) {
    /*
     * Constructs a block tree for the block.
     */
    EventTree *eventTree = flower_getEventTree(block_getFlower(block));
    //Now do the reconcilliation of events to events in the event tree.
    reconcile(binaryTree, eventTree,
            (Event *(*)(struct BinaryTree *, void *)) getEvent, segments);
    //Now build the tree
    Segment *mostAncestralEvent = buildChainTrees3P(block, segments,
            blockNumber, binaryTree);
    assert(block_getInstanceNumber(block)> 0);
    assert(mostAncestralEvent != NULL);
    //Make a root event.
    Event *rootEvent = eventTree_getRootEvent(eventTree);
    Segment *rootSegment = segment_construct(block, rootEvent);
    assert(event_isAncestor(segment_getEvent(mostAncestralEvent), rootEvent));
    segment_makeParentAndChild(rootSegment, mostAncestralEvent);
    block_setRootInstance(block, rootSegment);

#ifndef NDEBUG //Now go through all events checking they have a parent.
    Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
    Segment *segment;
    while ((segment = block_getNext(instanceIterator)) != NULL) {
        Segment *parent = segment_getParent(segment);
        if (parent == NULL) {
            assert(segment == rootSegment);
        }
    }
    block_destructInstanceIterator(instanceIterator);
#endif
}

void buildChainTrees2(ChainAlignment *chainAlignment,
        struct BinaryTree **refinedBlockTrees, int64_t *refinedBlockBoundaries,
        int64_t refinedBlockNumber) {
    /*
     * Iterates through a chain alignment, constructing the block trees and splitting blocks as needed.
     */
    assert(chainAlignment->columnNumber <= refinedBlockNumber);
    int64_t i, j, k;
    Block *block, *leftBlock, *rightBlock;

    j = 0;
    k = 0;
    for (i = 0; i < chainAlignment->columnNumber; i++) {
        block = chainAlignment->blocks[i];
        k += block_getLength(block);
        //walk along the refined block boundaries, splitting the considered as needed.
        assert(j < refinedBlockNumber && refinedBlockBoundaries[j] <= k);
        do {
            if (refinedBlockBoundaries[j] < k) { //we need to split the block
                assert(refinedBlockBoundaries[j] >= k - block_getLength(block)); //boundary must break block so that left block is at least one base pair long.
                assert(refinedBlockBoundaries[j] < k); //boundary must break block so that right block is at least one base pair long.
                block_split(block, block_getLength(block) - (k
                        - refinedBlockBoundaries[j]), &leftBlock, &rightBlock);
                buildChainTrees3(leftBlock, chainAlignment->matrix[i],
                        chainAlignment->rowNumber, refinedBlockTrees[j]);
                block = rightBlock;
                assert(k - block_getLength(block) == refinedBlockBoundaries[j]); //check the split did what we expect
                assert(j+1 < refinedBlockNumber && refinedBlockBoundaries[j+1] <= k); //check that we have another block tree to deal with the right side of the split.
            } else {
                buildChainTrees3(block, chainAlignment->matrix[i],
                        chainAlignment->rowNumber, refinedBlockTrees[j]);
            }
            j++;
        } while (j < refinedBlockNumber && refinedBlockBoundaries[j] <= k);
    }
    assert(j == refinedBlockNumber);
}

void buildChainTrees(ChainAlignment **chainAlignments,
        int64_t chainAlignmentNumber, EventTree *eventTree) {
    /*
     * This function builds a load of inputs which are then passed to Bernard's code.
     */
    int64_t i, j;

    //Make the block alignments.
    char ***concatenatedBlocks = st_malloc(sizeof(void *)
            * chainAlignmentNumber); //this is the list of 2d alignments.
    Name **_5Ends = st_malloc(sizeof(void *) * chainAlignmentNumber); //these are the lists of ends associated with each end.
    Name **_3Ends = st_malloc(sizeof(void *) * chainAlignmentNumber);
    Name **leafEventLabels = st_malloc(sizeof(void *) * chainAlignmentNumber); //this is the list of leaf event labels.
    int64_t **blockBoundaries =
            st_malloc(sizeof(void *) * chainAlignmentNumber); //each chain alignment has a list of block lengths to demark where the block boundaries are.

    //now fill out the the various arrays.
    for (i = 0; i < chainAlignmentNumber; i++) {
        concatenatedBlocks[i] = chainAlignment_getAlignment(chainAlignments[i]);
        _5Ends[i] = chainAlignment_getEndNames(chainAlignments[i], 1);
        _3Ends[i] = chainAlignment_getEndNames(chainAlignments[i], 0);
        leafEventLabels[i] = chainAlignment_getLeafEventNames(
                chainAlignments[i]);
        blockBoundaries[i] = chainAlignment_getBlockBoundaries(
                chainAlignments[i]);
    }
    //Event tree string
    char *eventTreeString = eventTree_makeNewickString(eventTree);

    //call to Bernard's code
    char *augmentedEventTreeString; //pointer to string holding the augmented event tree.
    char ***blockTreeStrings; //array of string pointers, for holding the constructed block trees.
    int64_t **refinedBlockBoundaries; //like the block boundaries, but revised by allowing for splits in the existing blocks.
    int64_t *refinedBlockNumbers; //the lengths of the block boundary arrays.
    buildChainTrees_Bernard(chainAlignmentNumber, concatenatedBlocks, _5Ends,
            _3Ends, leafEventLabels, blockBoundaries, eventTreeString,
            chainAlignments, &augmentedEventTreeString, &blockTreeStrings,
            &refinedBlockBoundaries, &refinedBlockNumbers);
    st_logDebug("Ran Bernard's code apparently okay\n");

    /*
     * Now process and reconcile each new block tree.
     */
    st_logDebug("Starting process and reconcile of new block tree\n");
    for (i = 0; i < chainAlignmentNumber; i++) {
        struct BinaryTree **blockTrees = st_malloc(sizeof(void *)
                * refinedBlockNumbers[i]);
        for (j = 0; j < refinedBlockNumbers[i]; j++) {
            blockTrees[j] = newickTreeParser(blockTreeStrings[i][j], 0.0, 0);
        }
        buildChainTrees2(chainAlignments[i], blockTrees,
                refinedBlockBoundaries[i], refinedBlockNumbers[i]);
        for (j = 0; j < refinedBlockNumbers[i]; j++) {
            destructBinaryTree(blockTrees[j]);
        }
        free(blockTrees);
    }
    st_logDebug("Processed the new block trees\n");

    /*
     * Cleanup the inputs.
     */
    for (i = 0; i < chainAlignmentNumber; i++) {
        free(concatenatedBlocks[i]);
        free(_5Ends[i]);
        free(_3Ends[i]);
        free(leafEventLabels[i]);
        free(blockBoundaries[i]);
    }
    free(concatenatedBlocks);
    free(_5Ends);
    free(_3Ends);
    free(leafEventLabels);
    free(blockBoundaries);
    free(eventTreeString);
    st_logDebug("Cleaned up the inputs\n");

    //done!
}

int chainAlignment_cmpFn(ChainAlignment **cA1, ChainAlignment **cA2) {
    /*
     * Compares chain alignments by their total base pair alignment length, sorting them in descending order.
     */
    return (*cA2)->totalAlignmentLength - (*cA1)->totalAlignmentLength;
}

static int oComparator(const void *o1, const void *o2, void *a) {
    /*
     * Compares the objects by there address.
     */
    assert(a == NULL);
    return o1 > o2 ? 1 : o1 < o2 ? -1 : 0;
}

ChainAlignment *chainAlignment_construct(Block **blocks, int64_t blocksLength) {
    /*
     * Constructs a chain alignment structure from a chain of blocks.
     */
    int64_t i, j, k;
    Block *block;
    Segment *segment;
    Segment *segment2;
    Cap *cap;
    ChainAlignment *chainAlignment;
    struct hashtable *hash;
    struct List *list;
    struct List *list2;
    struct avl_table *avlTable = avl_create(oComparator, NULL, NULL);
    avl_destroy(avlTable, NULL);

    /*
     * First iterate through all the segments in the chain, in order, to construct instances of the chain.
     */
    hash = create_hashtable(1, hashtable_key, hashtable_equalKey, NULL, NULL); //to keep track of the instances included in a chain.
    list = constructEmptyList(0, (void(*)(void *)) destructList); //the list of chain instances.
    k = 0;
    assert(blocksLength> 0);
    for (i = 0; i < blocksLength; i++) {
        block = blocks[i];
        Block_InstanceIterator *instanceIterator = block_getInstanceIterator(
                block);
        while ((segment = block_getNext(instanceIterator)) != NULL) {
            k++;
            assert(segment_getOrientation(segment));
            if (hashtable_search(hash, segment) == NULL) { //not yet in a chain instance
                list2 = constructEmptyList(blocksLength, NULL);
                for (j = 0; j < blocksLength; j++) { //this list will contain one instance for each block, or NULL, of missing.
                    list2->list[j] = NULL;
                }
                listAppend(list, list2);
                j = i; //start from the block we're up to.
                while (1) {
                    assert(segment_getOrientation(segment));
                    hashtable_insert(hash, segment, segment); //put in the hash to say we've seen it.
                    list2->list[j++] = segment; //put in the list of chains list.
                    if (j == blocksLength) { //end of chain
                        break;
                    }
                    cap = cap_getAdjacency(segment_get3Cap(segment));
                    //assert(!end_isStub(cap_getEnd(cap))); //can not be inherited.
                    if (end_isStubEnd(cap_getEnd(cap))) { //terminates with missing information.
                        break;
                    }
                    segment2 = cap_getSegment(cap);
                    assert(segment != NULL); //must be connected to a segment (not a cap).
                    assert(segment_getBlock(segment2) == blocks[j] || segment_getBlock(segment2) == block_getReverse(blocks[j-1])); //must be connected to reverse of itself or the next block
                    if (segment_getBlock(segment2) != blocks[j]) {
                        break;
                    }
                    segment = segment2;
                }
            } else {
                assert(i != 0);
            }
        }
        block_destructInstanceIterator(instanceIterator);
    }
    assert(k == (int64_t)hashtable_count(hash));
    assert(list->length > 0);

    /*
     * Now convert the chain instances in the list into the desired format for the chain alignment.
     */

    //alloc the chain alignment and the matrix of instances.
    chainAlignment = st_malloc(sizeof(ChainAlignment));
    chainAlignment->matrix = st_malloc(sizeof(Segment **) * blocksLength);
    for (i = 0; i < blocksLength; i++) {
        chainAlignment->matrix[i] = st_malloc(sizeof(Segment *) * list->length);
    }
    //fill out the fields of the chain alignment, including the matrix.
    chainAlignment->columnNumber = blocksLength;
    chainAlignment->rowNumber = list->length;
    for (i = 0; i < list->length; i++) {
        list2 = list->list[i];
        assert(list2->length == blocksLength);
        for (j = 0; j < blocksLength; j++) {
            chainAlignment->matrix[j][i] = list2->list[j];
        }
    }
    //Calculate the total length.
    chainAlignment->totalAlignmentLength = 0;
    for (i = 0; i < blocksLength; i++) {
        chainAlignment->totalAlignmentLength += block_getLength(blocks[i]);
    }
    assert(chainAlignment->totalAlignmentLength > 0);
    //Add chain of blocks.
    chainAlignment->blocks = st_malloc(sizeof(void *) * blocksLength);
    for (i = 0; i < blocksLength; i++) {
        chainAlignment->blocks[i] = blocks[i];
    }

    //Cleanup
    destructList(list);
    hashtable_destroy(hash, FALSE, FALSE);

    return chainAlignment;
}

void chainAlignment_destruct(ChainAlignment *chainAlignment) {
    /*
     * Destructs a chain alignment.
     */
    int64_t i;
    for (i = 0; i < chainAlignment->columnNumber; i++) {
        free(chainAlignment->matrix[i]);
    }
    free(chainAlignment->matrix);
    free(chainAlignment);
}

Event *copyConstructUnaryEvent(Event *event, EventTree *eventTree2) {
    /*
     * Adds the unary event to the event tree, allowing for the possibility that other unary events in the tree
     * of the event are not yet present in eventTree2.
     */
    assert(event_getChildNumber(event) == 1);

    double branchLength = event_getBranchLength(event);

    //Search for the first ancestor of event which is also in eventTree2, adding to the branch length as we go.
    Event *parentEvent = event_getParent(event);
    while (eventTree_getEvent(eventTree2, event_getName(parentEvent)) == NULL) {
        branchLength += event_getBranchLength(parentEvent);
        parentEvent = event_getParent(parentEvent);
        assert(parentEvent != NULL);
    } //at this point branch length is equal to branch length in eventTree2 from new event to common ancestor in both trees.
    parentEvent = eventTree_getEvent(eventTree2, event_getName(parentEvent)); //now get the event in the other tree.
    assert(parentEvent != NULL);

    assert(eventTree_getEvent(eventTree2, event_getName(event)) == NULL); //check it isn't already in the tree.
    //Search for the first descendant of the event which is also in eventTree2
    Event *childEvent = event_getChild(event, 0);
    assert(childEvent != NULL);
    while (eventTree_getEvent(eventTree2, event_getName(childEvent)) == NULL) {
        assert(event_getChildNumber(childEvent) == 1);
        childEvent = event_getChild(childEvent, 0);
        assert(childEvent != NULL);
    } //at this point the child event is present in both event trees.
    childEvent = eventTree_getEvent(eventTree2, event_getName(childEvent)); //get the child event.
    assert(childEvent != NULL);
    assert(event_getParent(childEvent) != NULL);
    assert(event_getParent(childEvent) == parentEvent);
    //assert(event_getBranchLength(childEvent) >= branchLength);

    return event_construct2(event_getName(event), event_getHeader(event), branchLength,
            parentEvent, childEvent, eventTree2);
}

void usage() {
    fprintf(
            stderr,
            "cactus_tree [flower-names, ordered by order they should be processed], version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    /*
     * The script builds trees.
     *
     * (1) It reads the flower.
     *
     * (2) It builds trees for the blocks.
     *
     * (3) It augments the event tree.
     *
     * (4) It copies the relevant block end trees into the ends of its descendant flowers.
     *
     * (5) It copies the relevant augmented events into the event trees of its descendants.
     *
     */
    CactusDisk *cactusDisk;
    int64_t startTime;
    int64_t i, j;
    Chain *chain;
    Block *block;
    struct List *sortedChainAlignments;
    Group *group;
    Flower *flower2;
    End *end;
    End *end2;
    Cap *cap;
    Cap *cap2;
    Cap *cap3;
    Event *event;
    Flower_EndIterator *endIterator;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument, 0,
                'c' }, { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:h", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    //assert(logLevelString == NULL || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
    assert(cactusDiskDatabaseString != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    //////////////////////////////////////////////
    //For each flower do tree building..
    //////////////////////////////////////////////

    stList *flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    preCacheNestedFlowers(cactusDisk, flowers);
    for(j = 0; j < stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);
        st_logInfo("Processing a flower\n");

        ///////////////////////////////////////////////////////////////////////////
        // Do nothing if we have already built the trees.
        ///////////////////////////////////////////////////////////////////////////

        if (flower_builtTrees(flower)) {
            st_logInfo("We have already built trees for flower %s\n", cactusMisc_nameToStringStatic(flower_getName(flower)));
            continue;
        }

        ///////////////////////////////////////////////////////////////////////////
        //Sets up the 'trees' for the caps for the top level problem.
        //In other words, it ensure ends from stubs have there root instance set - only needs doing for
        //the top level flower, after which the stub roots will be set recursively.
        ///////////////////////////////////////////////////////////////////////////

        if (flower_getName(flower) == 0) {
            endIterator = flower_getEndIterator(flower);
            while ((end = flower_getNextEnd(endIterator)) != NULL) {
                if (end_isStubEnd(end)) {
                    assert(end_getInstanceNumber(end) == 1);
                    Event *rootEvent = eventTree_getRootEvent(flower_getEventTree(
                            flower));
                    Cap *childCap = end_getFirst(end);
                    Cap *rootCap = cap_construct(end, rootEvent);
                    assert(event_isAncestor(cap_getEvent(childCap), rootEvent));
                    cap_makeParentAndChild(rootCap, childCap);
                    end_setRootInstance(end, rootCap);
                } else {
                    assert(!end_isStubEnd(end));
                    assert(end_isBlockEnd(end));
                }
            }
            flower_destructEndIterator(endIterator);
        }

#ifndef NDEBUG
        endIterator = flower_getEndIterator(flower);
        while ((end = flower_getNextEnd(endIterator)) != NULL) {
            if (!end_isBlockEnd(end)) {
                assert(end_getRootInstance(end) != NULL);
            }
        }
        flower_destructEndIterator(endIterator);
#endif

        ///////////////////////////////////////////////////////////////////////////
        //Construct the chain alignments for all the non-trivial chains.
        ///////////////////////////////////////////////////////////////////////////

        startTime = time(NULL);
        sortedChainAlignments = constructEmptyList(0,
                (void(*)(void *)) chainAlignment_destruct);
        Flower_ChainIterator *chainIterator = flower_getChainIterator(flower);
        while ((chain = flower_getNextChain(chainIterator)) != NULL) {
            Block **blockChain = chain_getBlockChain(chain, &i);
            if (i > 0) {
                listAppend(sortedChainAlignments, chainAlignment_construct(
                        blockChain, i));
            }
            free(blockChain);
        }
        flower_destructChainIterator(chainIterator);
        st_logInfo(
                "Constructed the block trees for the non-trivial chains in the flower in: %" PRIi64 " seconds\n",
                time(NULL) - startTime);

        ///////////////////////////////////////////////////////////////////////////
        //Construct the chain alignment for each trivial chain.
        ///////////////////////////////////////////////////////////////////////////

        startTime = time(NULL);
        Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
        while ((block = flower_getNextBlock(blockIterator)) != NULL) {
            if (block_getChain(block) == NULL) {
                listAppend(sortedChainAlignments, chainAlignment_construct(
                        &block, 1));
            }
        }
        flower_destructBlockIterator(blockIterator);

        qsort(sortedChainAlignments->list, sortedChainAlignments->length,
                sizeof(void *),
                (int(*)(const void *, const void *)) chainAlignment_cmpFn);
#ifndef NDEBUG
        for (i = 1; i < sortedChainAlignments->length; i++) {
            assert(((ChainAlignment *)sortedChainAlignments->list[i-1])->totalAlignmentLength >=
                    ((ChainAlignment *)sortedChainAlignments->list[i])->totalAlignmentLength);
        }
#endif

        st_logInfo(
                "Constructed the block trees for the trivial chains in the flower in: %" PRIi64 " seconds\n",
                time(NULL) - startTime);

        ///////////////////////////////////////////////////////////////////////////
        //For each chain call the tree construction function.
        ///////////////////////////////////////////////////////////////////////////

        startTime = time(NULL);
        buildChainTrees((ChainAlignment **) sortedChainAlignments->list,
                sortedChainAlignments->length, flower_getEventTree(flower));
        destructList(sortedChainAlignments);

        st_logInfo("Augmented the block trees in: %" PRIi64 " seconds\n", time(NULL)
                - startTime);

#ifndef NDEBUG
        endIterator = flower_getEndIterator(flower);
        while ((end = flower_getNextEnd(endIterator)) != NULL) {
            assert(end_getRootInstance(end) != NULL);
        }
        flower_destructEndIterator(endIterator);
#endif

        ///////////////////////////////////////////////////////////////////////////
        //Pass the end trees and augmented events to the child flowers.
        ///////////////////////////////////////////////////////////////////////////

        startTime = time(NULL);
        Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
        while ((group = flower_getNextGroup(groupIterator)) != NULL) {
            flower2 = group_getNestedFlower(group);
            if (flower2 != NULL) {
                //add in the end trees and augment the event trees.
                Group_EndIterator *endIterator = group_getEndIterator(group);
                while ((end = group_getNextEnd(endIterator)) != NULL) {
                    end2 = flower_getEnd(flower2, end_getName(end));
                    assert(end2 != NULL);
                    //copy the caps.
                    End_InstanceIterator *instanceIterator =
                            end_getInstanceIterator(end);
                    while ((cap = end_getNext(instanceIterator)) != NULL) {
                        if (end_getInstance(end2, cap_getName(cap)) == NULL) {
                            assert(cap_getChildNumber(cap)> 0); //can not be a leaf
                            //make sure the augmented event is in there.
                            event = cap_getEvent(cap);
                            if (eventTree_getEvent(flower_getEventTree(flower2),
                                    event_getName(event)) == NULL) {
                                assert(event_getChildNumber(event) == 1); //must be a unary event
                                copyConstructUnaryEvent(event,
                                        flower_getEventTree(flower2));
                            }
                            event = eventTree_getEvent(flower_getEventTree(flower2),
                                    event_getName(event));
                            assert(event != NULL);
                            cap_copyConstruct(end2, cap);
                        }
                    }
                    //now copy the parent links.
                    while ((cap = end_getPrevious(instanceIterator)) != NULL) {
                        cap2 = end_getInstance(end2, cap_getName(cap));
                        assert(cap2 != NULL);
                        if (cap_getParent(cap) != NULL) {
                            cap3 = end_getInstance(end2, cap_getName(
                                    cap_getParent(cap)));
                            assert(cap3 != NULL);
                            cap_makeParentAndChild(cap3, cap2);
                        } else {
                            assert(end_getRootInstance(end) != NULL);
                            assert(cap == end_getRootInstance(end));

                            assert(end_getRootInstance(end2) == NULL);
                            end_setRootInstance(end2, cap2);
                        }
                    }
                    end_destructInstanceIterator(instanceIterator);
                }
                group_destructEndIterator(endIterator);
            }
        }
        flower_destructGroupIterator(groupIterator);
        st_logInfo(
                "Filled in end trees and augmented the event trees for the child flowers in: %" PRIi64 " seconds\n",
                time(NULL) - startTime);

        ///////////////////////////////////////////////////////////////////////////
        //Set the trees in the flower to 'built' status.
        ///////////////////////////////////////////////////////////////////////////

        assert(!flower_builtTrees(flower));
        flower_setBuiltTrees(flower, 1);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Unload the parent flowers
    ///////////////////////////////////////////////////////////////////////////

    for(j = 0; j < stList_length(flowers); j++) {
        Flower *flower = stList_get(flowers, j);
        assert(flower != NULL);
        flower_unloadParent(flower); //We have this line just in case we are loading the parent..
    }
    stList_destruct(flowers);

    ///////////////////////////////////////////////////////////////////////////
    // (9) Write the flower to disk.
    ///////////////////////////////////////////////////////////////////////////

    startTime = time(NULL);
    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk in: %" PRIi64 " seconds\n", time(NULL)
            - startTime);

    ///////////////////////////////////////////////////////////////////////////
    //(15) Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);

    return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.

    //Destruct stuff
    startTime = time(NULL);
    stKVDatabaseConf_destruct(kvDatabaseConf);
    if(logLevelString != NULL) {
        free(logLevelString);
    }
    free(cactusDiskDatabaseString);

    st_logInfo("Cleaned stuff up and am finished in: %" PRIi64 " seconds\n", time(NULL)
            - startTime);
    return 0;
}
