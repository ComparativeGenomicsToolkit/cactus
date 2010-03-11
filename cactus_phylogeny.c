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

#define closeEnough 0.001

typedef struct _chainAlignment {
	/*
	 * Structure to represent a concatenated list of blocks as a single 2d alignment.
	 */
	//Matrix of alignment, matrix[column][row]
	Segment ***matrix; //NULL values are okay, where an instance of an block is missing from an instance of the chain.
	Block **blocks; //this list of blocks, in order.
	int32_t columnNumber; //the number of blocks.
	int32_t rowNumber; //the number of rows of the alignment, each row containing an instance of blocks in the chain.
	int32_t totalAlignmentLength; //the length in base pairs of the alignment.
} ChainAlignment;

char **chainAlignment_getAlignment(ChainAlignment *chainAlignment) {
	/*
	 * Constructs a concatenated 2d matrix of chars referenced by [block][instance], representing the
	 * base pair alignment of the chain alignment.
	 */
	char **alignment;
	int32_t i, j, k, l;
	Segment *segment;
	char *cA;

	//alloc the memory for the char alignment.
	alignment = malloc(sizeof(void *) * chainAlignment->rowNumber);
	for(i=0; i<chainAlignment->rowNumber; i++) {
		alignment[i] = malloc(sizeof(char) * (chainAlignment->totalAlignmentLength + 1));
	}

	/*
	 * Fills in the alignment.
	 */
	for(j=0; j<chainAlignment->rowNumber; j++) {
		l = 0;
		for(i=0; i<chainAlignment->columnNumber; i++) {
			segment = chainAlignment->matrix[i][j];
			if(segment == NULL) {
				for(k=0; k<block_getLength(chainAlignment->blocks[i]); k++) {
					alignment[j][l++]  = 'N';
				}
			}
			else {
				cA = segment_getString(segment);
				for(k=0; k<segment_getLength(segment); k++) {
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

Segment *chainAlignment_getFirstNonNullSegment(ChainAlignment *chainAlignment, int32_t row, bool increasing) {
	/*
	 * Gets the first instance on an row in the chain alignment which is non null. If increasing is false, gets the last.
	 */
	Segment *segment;
	int32_t j = 0, k = chainAlignment->columnNumber, l = 1;
	if(!increasing) {
		j = chainAlignment->columnNumber-1;
		k = -1;
		l = -1;
	}
	for(; j!=k; j += l) {
		segment = chainAlignment->matrix[j][row];
		if(segment != NULL) {
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
	int32_t i;
	Segment *segment;

	names = malloc(sizeof(Name) * chainAlignment->rowNumber);
	for(i=0; i<chainAlignment->rowNumber; i++) {
		segment = chainAlignment_getFirstNonNullSegment(chainAlignment, i, _5End);
		names[i] = end_getName(cap_getEnd(cap_getAdjacency(_5End ? segment_get5Cap(segment) : segment_get3Cap(segment))));
	}
	return names;
}

Name *chainAlignment_getLeafEventNames(ChainAlignment *chainAlignment) {
	/*
	 * Gets the names of the leaf events of the rows (instances), in the alignment.
	 */
	Name *names;
	int32_t i;
	Segment *segment;

	names = malloc(sizeof(Name) * chainAlignment->rowNumber);
	for(i=0; i<chainAlignment->rowNumber; i++) {
		segment = chainAlignment_getFirstNonNullSegment(chainAlignment, i, 1);
		names[i] = event_getName(segment_getEvent(segment));
	}
	return names;
}

int32_t *chainAlignment_getBlockBoundaries(ChainAlignment *chainAlignment) {
	/*
	 * Gets the boundaries of block in the chain alignment.
	 */
	int32_t *blockBoundaries;
	int32_t i, j;

	blockBoundaries = malloc(sizeof(int32_t) * chainAlignment->columnNumber);
	j = 0;
	for(i=0; i<chainAlignment->columnNumber; i++) {
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

void buildChainTrees_Bernard(int32_t blockNumber, char ***concatenatedBlocks, Name **_5Ends, Name **_3Ends, Name **leafEventLabels,
							int32_t **blockBoundaries, char *eventTreeString, const char *tempDir, ChainAlignment **chainAlignments,
							char **modifiedEventTreeString, char ****blockTreeStrings, int32_t ***refinedBlockBoundaries, int32_t **refinedBlockNumbers) {
	/*
	 * Here's the function you need to fill in, I haven't defined the outputs yet - you get it working and then we can discuss.
	 *
	 * Arrays are all indexed first by the chain/concatenated blocks.
	 * Alignments are column/row indexed with rows as chain instances and columns as aligned bases.
	 * So for example: concatenatedBlocks[i][j][k] is the ith chain/concatenated block, jth row (chain instance), kth column (position in concatenated block)
	 *
	 * Names can be converted to strings with: netMisc_nameToString() and netMisc_nameToStringStatic() (See the API).
	 *
	 * The output arrays (last three arguments) must be initialised and are *pointers to* the point containing the array/string. Which must be
	 * created.
	 */

	int32_t i = 0;
	int32_t j = 0;

	char *randomDirName = NULL;

	int32_t rowNumber = 0;
	int32_t colNumber = 0;

	char eventTreeFileName[] = "pre.event.tree";
	char blockFileName[] = "annot";
	char mafDirFileName[] = "mafs";
	char augTreeFileName[] = "remap.augtree";
	char dupTreeFileName[] = "event.duptree";

	const int32_t TMP_BUFFER_SIZE = 256; // Assume that 256 chars is enough for any temp filename strings
	char tmpStringBuffer[TMP_BUFFER_SIZE];
	FILE *fp = NULL;

	const int32_t LINE_BUFF_SIZE = 4096; // Assume that 4096 chars is enough for reading in a single line of a file
	char lineBuffer[LINE_BUFF_SIZE];
	char *tmpString = NULL;

	int32_t **newBlockBoundaries = NULL;
	int32_t *newBlockNumbers = NULL;

	char ***blockTreeArray = NULL;

	char *tempDirName;
	tempDirName = (char *) tempDir;
//	char debugDir[] = "/Users/bsuh/TESTDIR";
//	tempDirName = (char *) &debugDir;

	randomDirName = strrchr(tempDirName, '/');
	if (randomDirName == NULL) {
		; // Assume this doesn't happen
	} else {
		randomDirName += 1;
	}

	/* Output the eventTree to file */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, eventTreeFileName);
	fp = fopen(tmpStringBuffer, "w");
	fprintf(fp, "%s\n", eventTreeString);
	fclose(fp);

	/* Output the block definition file */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, blockFileName);
	fp = fopen(tmpStringBuffer, "w");
	for (i=0; i<blockNumber; i++) {
		rowNumber = chainAlignments[i]->rowNumber;

		/* NOTE: Block Numbers start at 1 */
		fprintf(fp, ">%d %d %d\n", i+1, chainAlignments[i]->totalAlignmentLength, rowNumber);
		for (j=0; j<rowNumber; j++) {
			fprintf(fp, "%s.chr0:1-%d + %s %s\n", netMisc_nameToString(leafEventLabels[i][j]), chainAlignments[i]->totalAlignmentLength, netMisc_nameToString(_5Ends[i][j]), netMisc_nameToString(_3Ends[i][j]));
		}
	}
	fclose(fp);

	/* Create the block map file */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, "block.map");
	fp = fopen(tmpStringBuffer, "w");
	for (i=0; i<blockNumber; i++) {
		/* NOTE: Block Numbers start at 1 */
		fprintf(fp, "%d\t%d\n", i+1, i+1);
	}
	fclose(fp);

	/* Create the maf files */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, mafDirFileName);
	mkdir(tmpStringBuffer, 0777);

	for (i=0; i<blockNumber; i++) {
		/* NOTE: Block Numbers start at 1 */
		snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s/%d", tempDirName, mafDirFileName, i+1);

		fp = fopen(tmpStringBuffer, "w");

		rowNumber = chainAlignments[i]->rowNumber;
		for (j=0; j<rowNumber; j++) {
			fprintf(fp, ">%s\n", netMisc_nameToString(leafEventLabels[i][j]));
			fprintf(fp, "%s\n", concatenatedBlocks[i][j]);
		}

		fclose(fp);
	}

	/* Run the tree pipeline */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "conTrees_PhyloBuilder.py %s", tempDirName);
	exitOnFailure(system(tmpStringBuffer), "conTrees_PhyloBuilder.py failed\n");
	printf("Completed running tree pipeline, now onto parsing\n");

	/* Clone Block boundaries to refined Block Boundaries*/
	newBlockBoundaries = malloc(sizeof(void *) * blockNumber);
	newBlockNumbers = malloc(sizeof(int32_t) * blockNumber);
	for (i=0; i<blockNumber; i++) {
		colNumber = chainAlignments[i]->columnNumber;
		newBlockBoundaries[i] = malloc(sizeof(int32_t) * colNumber);
		newBlockNumbers[i] = colNumber;
		for (j=0; j<colNumber; j++) {
			newBlockBoundaries[i][j] = blockBoundaries[i][j];
		}
	}
	*refinedBlockBoundaries = newBlockBoundaries;
	*refinedBlockNumbers = newBlockNumbers;

	/* Read in tree pipeline output for the event tree */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, dupTreeFileName);
	fp = fopen(tmpStringBuffer, "r");
	if (fp != NULL) {
		if (fgets(lineBuffer, LINE_BUFF_SIZE, fp) != NULL) {
			chomp(lineBuffer);
			tmpString = stringCopy(lineBuffer);
			*modifiedEventTreeString = tmpString;
		} else {
			perror("Failed to read new eventTree file");
		}
	} else {
		perror("Failed to open new eventTree file");
	}
	fclose(fp);

	/* Read in tree pipeline output for the aug tree */
	blockTreeArray = malloc(sizeof(void *) * blockNumber);
	for (i=0; i<blockNumber; i++) {
		colNumber = chainAlignments[i]->columnNumber;
		blockTreeArray[i] = malloc(sizeof(void *) * colNumber);
	}
	i = 0; // Variable to store blockNumber
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, augTreeFileName);
	fp = fopen(tmpStringBuffer, "r");
	if (fp != NULL) {
		while(fgets(lineBuffer, LINE_BUFF_SIZE, fp) != NULL) {
			chomp(lineBuffer);
			if (lineBuffer[0] == '>') {
				sscanf(lineBuffer, ">%d", &i);
				i -= 1; // Block number ids start at 1
			} else {
				colNumber = chainAlignments[i]->columnNumber;
				for (j=0; j<colNumber; j++) {
					tmpString = stringCopy(lineBuffer);
					blockTreeArray[i][j] = tmpString;
				}
			}
		}
	} else {
		perror("Failed to open augTree file");
	}
	fclose(fp);
	*blockTreeStrings = blockTreeArray;

	return;
}

bool isNewEvent(const char *eventName) {
	return strlen(eventName) >= 8 && eventName[0] == 'N' && eventName[1] == 'E' && eventName[2] == 'W';
}

Event *augmentEventTree(struct BinaryTree *augmentedEventTree,
		EventTree *eventTree, struct hashtable *newEventNameMap) {
	/*
	 * Function takes an augmented event tree and adds in the extra (unary) events to the original
	 * event tree.
	 */
	Event *event;

	if(augmentedEventTree->internal) {
		Event *childEvent = augmentEventTree(augmentedEventTree->left, eventTree, newEventNameMap);
		if(augmentedEventTree->right != NULL) { //is a speciation node.
			Event *childEvent2 = augmentEventTree(augmentedEventTree->right, eventTree, newEventNameMap);
			assert(!isNewEvent(augmentedEventTree->label));
			Name eventName = netMisc_stringToName(augmentedEventTree->label);
			event = eventTree_getEvent(eventTree, eventName);
#ifdef BEN_DEBUG
			assert(event != NULL);
			assert(event_getChildNumber(event) == 2);
			assert(event_getParent(childEvent) == event);
			assert(event_getParent(childEvent2) == event);
			assert(floatValuesClose(event_getBranchLength(childEvent), augmentedEventTree->left->distance, closeEnough));
			assert(floatValuesClose(event_getBranchLength(childEvent2), augmentedEventTree->right->distance, closeEnough));
#endif
			return event;
		}
		else { //a unary event
			if(isNewEvent(augmentedEventTree->label)) {
				MetaEvent *metaEvent = metaEvent_construct("", net_getNetDisk(eventTree_getNet(eventTree)));
				//We set the branch length so that of the child branch is correct.
				assert(augmentedEventTree->left->distance <= event_getBranchLength(childEvent));
				event = event_construct2(metaEvent, event_getBranchLength(childEvent) - augmentedEventTree->left->distance, event_getParent(childEvent), childEvent, eventTree);
				hashtable_insert(newEventNameMap,  //add to the map of new event names.
								 stringCopy(augmentedEventTree->label),
								 netMisc_nameToString(event_getName(event)));
			}
			else {
				Name eventName = netMisc_stringToName(augmentedEventTree->label);
				event = eventTree_getEvent(eventTree, eventName);
				assert(event != NULL);
			}
#ifdef BEN_DEBUG
			assert(event_getParent(childEvent) == event);
			assert(floatValuesClose(event_getBranchLength(childEvent), augmentedEventTree->left->distance, closeEnough));
#endif
			return event;
		}
	}
	else { //is a leaf
		assert(!isNewEvent(augmentedEventTree->label));
		Name eventName = netMisc_stringToName(augmentedEventTree->label);
		event = eventTree_getEvent(eventTree, eventName);
#ifdef BEN_DEBUG
		assert(event != NULL);
		assert(event_getChildNumber(event) == 0);
#endif
		return event;
	}
}


Segment *buildChainTrees3P(Block *block, Segment **segments, int32_t blockNumber,
		struct BinaryTree *binaryTree, struct hashtable *newEventNameMap) {
	/*
	 * Recursive partner to buildChainTree3 function, recurses on the binary tree constructing the block tree.
	 * The labels of the leaves are indexes into the segments array, the internal node's labels are events in the event tree.
	 */
	if(binaryTree->internal) { //deal with an internal node of the block tree.
		Segment *leftInstance = buildChainTrees3P(block, segments, blockNumber, binaryTree->left, newEventNameMap);
		Segment *rightInstance = buildChainTrees3P(block, segments, blockNumber, binaryTree->right, newEventNameMap);
		if(leftInstance != NULL) {
			if(rightInstance != NULL) {

				Event *event;
				if(isNewEvent(binaryTree->label)) {
					const char *cA = hashtable_search(newEventNameMap, binaryTree->label);
					assert(cA != NULL);
					event = eventTree_getEvent(net_getEventTree(block_getNet(block)), netMisc_stringToName(cA));
				}
				else {
					event = eventTree_getEvent(net_getEventTree(block_getNet(block)), netMisc_stringToName(binaryTree->label));
					if(event != NULL) {
						// printBinaryTree(stderr, binaryTree);
					}
				}

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
	}
	else { //a leaf, so find the leaf instance in the list of segments (this may be null if missing data).
		int32_t i;
		assert(sscanf(binaryTree->label, "%i", &i) == 1);
		assert(i < blockNumber);
		assert(i >= 0);
		return segments[i];
	}
}

void buildChainTrees3(Block *block, Segment **segments, int32_t blockNumber, struct BinaryTree *binaryTree, struct hashtable *newEventNameMap) {
	/*
	 * Constructs an block tree for the block.
	 */
	Segment *mostAncestralEvent = buildChainTrees3P(block, segments, blockNumber, binaryTree, newEventNameMap);
	assert(block_getInstanceNumber(block) > 0);
	assert(mostAncestralEvent != NULL);
	//Make a root event.
	Event *rootEvent = eventTree_getRootEvent(net_getEventTree(block_getNet(block)));
	Segment *rootSegment = segment_construct(block, rootEvent);
	assert(event_isAncestor(segment_getEvent(mostAncestralEvent), rootEvent));
	segment_makeParentAndChild(rootSegment, mostAncestralEvent);
	block_setRootInstance(block, rootSegment);

#ifdef BEN_DEBUG //Now go through all events checking they have a parent.
	Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
	Segment *segment;
	while((segment = block_getNext(instanceIterator)) != NULL) {
		Segment *parent = segment_getParent(segment);
		if(parent == NULL) {
			assert(segment == rootSegment);
		}
	}
	block_destructInstanceIterator(instanceIterator);
#endif
}

void buildChainTrees2(ChainAlignment *chainAlignment,
					  struct BinaryTree **refinedBlockTrees,
					  int32_t *refinedBlockBoundaries,
					  int32_t refinedBlockNumber,
					  struct hashtable *newEventNameMap) {
	/*
	 * Iterates through a chain alignment, constructing the block trees and splitting blocks as needed.
	 */
	assert(chainAlignment->columnNumber <= refinedBlockNumber);
	int32_t i, j, k;
	Block *block, *leftBlock, *rightBlock;

	j=0;
	k=0;
	for(i=0; i < chainAlignment->columnNumber; i++) {
		block = chainAlignment->blocks[i];
		k += block_getLength(block);
		//walk along the refined block boundaries, splitting the considered as needed.
		assert(j < refinedBlockNumber && refinedBlockBoundaries[j] <= k);
		do {
			if(refinedBlockBoundaries[j] < k) { //we need to split the block
				assert(refinedBlockBoundaries[j] >= k - block_getLength(block)); //boundary must break block so that left block is at least one base pair long.
				assert(refinedBlockBoundaries[j] < k); //boundary must break block so that right block is at least one base pair long.
				block_split(block, block_getLength(block) - (k - refinedBlockBoundaries[j]), &leftBlock, &rightBlock);
				buildChainTrees3(leftBlock, chainAlignment->matrix[i], chainAlignment->rowNumber, refinedBlockTrees[j], newEventNameMap);
				block = rightBlock;
				assert(k - block_getLength(block) == refinedBlockBoundaries[j]); //check the split did what we expect
				assert(j+1 < refinedBlockNumber && refinedBlockBoundaries[j+1] <= k); //check that we have another block tree to deal with the right side of the split.
			}
			else {
				buildChainTrees3(block, chainAlignment->matrix[i], chainAlignment->rowNumber, refinedBlockTrees[j], newEventNameMap);
			}
			j++;
		} while(j < refinedBlockNumber && refinedBlockBoundaries[j] <= k);
	}
	assert(j == refinedBlockNumber);
}

void buildChainTrees(ChainAlignment **chainAlignments, int32_t chainAlignmentNumber, EventTree *eventTree,
		const char *tempFilePath) {
	/*
	 * This function builds a load of inputs which are then passed to Bernard's code.
	 */
	int32_t i, j;

	//Make the block alignments.
	char ***concatenatedBlocks = malloc(sizeof(void *) * chainAlignmentNumber); //this is the list of 2d alignments.
	Name **_5Ends = malloc(sizeof(void *) * chainAlignmentNumber); //these are the lists of ends associated with each end.
	Name **_3Ends = malloc(sizeof(void *) * chainAlignmentNumber);
	Name **leafEventLabels = malloc(sizeof(void *) * chainAlignmentNumber); //this is the list of leaf event labels.
	int32_t **blockBoundaries = malloc(sizeof(void *) * chainAlignmentNumber); //each chain alignment has a list of block lengths to demark where the block boundaries are.

	//now fill out the the various arrays.
	for(i=0; i<chainAlignmentNumber; i++) {
		concatenatedBlocks[i] = chainAlignment_getAlignment(chainAlignments[i]);
		_5Ends[i] = chainAlignment_getEndNames(chainAlignments[i], 1);
		_3Ends[i] = chainAlignment_getEndNames(chainAlignments[i], 0);
		leafEventLabels[i] = chainAlignment_getLeafEventNames(chainAlignments[i]);
		blockBoundaries[i] = chainAlignment_getBlockBoundaries(chainAlignments[i]);
	}
	//Event tree string
	char *eventTreeString = eventTree_makeNewickString(eventTree);

	//Construct the random dir.
	char *randomDir;
	exitOnFailure(constructRandomDir(tempFilePath, &randomDir), "Tried to make a recursive directory of temp files but failed\n");

	//call to Bernard's code
	char *augmentedEventTreeString; //pointer to string holding the augmented event tree.
	char ***blockTreeStrings; //array of string pointers, for holding the constructed block trees.
	int32_t **refinedBlockBoundaries; //like the block boundaries, but revised by allowing for splits in the existing blocks.
	int32_t *refinedBlockNumbers; //the lengths of the block boundary arrays.
	buildChainTrees_Bernard(chainAlignmentNumber, concatenatedBlocks, _5Ends, _3Ends, leafEventLabels,
							blockBoundaries, eventTreeString, randomDir, chainAlignments,
							&augmentedEventTreeString, &blockTreeStrings, &refinedBlockBoundaries, &refinedBlockNumbers);
	logDebug("Ran Bernard's code apparently okay\n");

	/*
	 * Clean up the temporary files directory.
	 */
	exitOnFailure(destructRandomDir(randomDir), "Tried to destroy a recursive directory of temp files but failed\n");

	/*
	 * Augment the event tree with the new events.
	 */
	struct BinaryTree *modifiedEventTree = newickTreeParser(augmentedEventTreeString, 0.0, 1);
	struct hashtable *newEventNameMap = create_hashtable(1, hashtable_stringHashKey,
			hashtable_stringEqualKey, free, free);
	augmentEventTree(modifiedEventTree, eventTree, newEventNameMap);
	logDebug("Augmented the event tree\n");

	/*
	 * Now process each new block tree.
	 */
	for(i=0; i<chainAlignmentNumber; i++) {
		struct BinaryTree **blockTrees = malloc(sizeof(void *)* refinedBlockNumbers[i]);
		for(j=0; j<refinedBlockNumbers[i]; j++) {
			blockTrees[j] = newickTreeParser(blockTreeStrings[i][j], 0.0, 0);
		}
		buildChainTrees2(chainAlignments[i], blockTrees, refinedBlockBoundaries[i], refinedBlockNumbers[i], newEventNameMap);
		for(j=0; j<refinedBlockNumbers[i]; j++) {
			destructBinaryTree(blockTrees[j]);
		}
		free(blockTrees);
	}
	logDebug("Processed the new block trees\n");

	/*
	 * Cleanup the inputs.
	 */
	for(i=0; i<chainAlignmentNumber; i++) {
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
	hashtable_destroy(newEventNameMap, 1, 1);
	logDebug("Cleaned up the inputs\n");

	//done!
}

int chainAlignment_cmpFn(ChainAlignment **cA1, ChainAlignment **cA2) {
	/*
	 * Compares chain alignments by their total base pair alignment length, sorting them in descending order.
	 */
	return (*cA2)->totalAlignmentLength - (*cA1)->totalAlignmentLength;
}

static int32_t oComparator(const void *o1, const void *o2, void *a) {
	/*
	 * Compares the objects by there address.
	 */
	assert(a == NULL);
	return o1 > o2 ? 1 : o1 < o2 ? -1 : 0;
}

ChainAlignment *chainAlignment_construct(Block **blocks, int32_t blocksLength) {
	/*
	 * Constructs a chain alignment structure from a chain of blocks.
	 */
	int32_t i, j, k;
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
	list = constructEmptyList(0, (void (*)(void *))destructList); //the list of chain instances.
	k = 0;
	assert(blocksLength > 0);
	for(i=0; i<blocksLength; i++) {
		block = blocks[i];
		Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
		while((segment = block_getNext(instanceIterator)) != NULL) {
			k++;
			assert(segment_getOrientation(segment));
			if(hashtable_search(hash, segment) == NULL) { //not yet in a chain instance
				list2 = constructEmptyList(blocksLength, NULL);
				for(j=0; j<blocksLength; j++) { //this list will contain one instance for each block, or NULL, of missing.
					list2->list[j] = NULL;
				}
				listAppend(list, list2);
				j = i; //start from the block we're up to.
				while(1) {
					assert(segment_getOrientation(segment));
					hashtable_insert(hash, segment, segment); //put in the hash to say we've seen it.
					list2->list[j++] = segment; //put in the list of chains list.
					if(j == blocksLength) { //end of chain
						break;
					}
					cap = cap_getAdjacency(segment_get3Cap(segment));
					assert(!end_isCap(cap_getEnd(cap))); //can not be a cap.
					if(end_isStub(cap_getEnd(cap))) { //terminates with missing information.
						break;
					}
					segment2 = cap_getSegment(cap);
					assert(segment != NULL); //must be connected to an segment (not a cap).
					assert(segment_getBlock(segment2) == blocks[j] || segment_getBlock(segment2) == block_getReverse(blocks[j-1])); //must be connected to reverse of itself or the next block
					if(segment_getBlock(segment2) != blocks[j]) {
						break;
					}
					segment = segment2;
				}
			}
			else {
				assert(i != 0);
			}
		}
		block_destructInstanceIterator(instanceIterator);
	}
	assert(k == (int32_t)hashtable_count(hash));
	assert(list->length > 0);

	/*
	 * Now convert the chain instances in the list into the desired format for the chain alignment.
	 */

	//alloc the chain alignment and the matrix of instances.
	chainAlignment = malloc(sizeof(ChainAlignment));
	chainAlignment->matrix = malloc(sizeof(Segment **) * blocksLength);
	for(i=0; i<blocksLength; i++) {
		chainAlignment->matrix[i] = malloc(sizeof(Segment *) * list->length);
	}
	//fill out the fields of the chain alignment, including the matrix.
	chainAlignment->columnNumber = blocksLength;
	chainAlignment->rowNumber = list->length;
	for(i=0; i<list->length; i++) {
		list2 = list->list[i];
		assert(list2->length == blocksLength);
		for(j=0; j<blocksLength; j++) {
			chainAlignment->matrix[j][i] = list2->list[j];
		}
	}
	//Calculate the total length.
	chainAlignment->totalAlignmentLength = 0;
	for(i=0;i<blocksLength; i++) {
		chainAlignment->totalAlignmentLength += block_getLength(blocks[i]);
	}
	assert(chainAlignment->totalAlignmentLength > 0);
	//Add chain of blocks.
	chainAlignment->blocks = malloc(sizeof(void *)*blocksLength);
	for(i=0; i<blocksLength; i++) {
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
	int32_t i;
	for(i=0; i<chainAlignment->columnNumber; i++) {
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
	while(eventTree_getEvent(eventTree2, event_getName(parentEvent)) == NULL) {
		branchLength += event_getBranchLength(parentEvent);
		parentEvent = event_getParent(parentEvent);
		assert(parentEvent != NULL);
	} //at this point branch length is equal to branch length in eventTree2 from new event to common ancestor in both trees.
	parentEvent = eventTree_getEvent(eventTree2, event_getName(parentEvent)); //now get the event in the other tree.
	assert(parentEvent != NULL);

	//Search for the first descendant of the event which is also in eventTree2, adding to the branch length as we go.
	Event *childEvent = event_getChild(event, 0);
	while(eventTree_getEvent(eventTree2, event_getName(childEvent)) == NULL) {
		assert(event_getChildNumber(childEvent) == 1);
		childEvent = event_getChild(parentEvent, 0);
		assert(childEvent != NULL);
	} //at this point the child event is present in both event trees.
	childEvent = eventTree_getEvent(eventTree2, event_getName(childEvent)); //get the child event.
	assert(childEvent != NULL);
	//assert(event_getBranchLength(childEvent) >= branchLength);

	return event_construct2(event_getMetaEvent(event), branchLength, parentEvent, childEvent, eventTree2);
}

void usage() {
	fprintf(stderr, "cactus_tree [net-names, ordered by order they should be processed], version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-e --tempDirRoot : The temp file root directory\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
	/*
	 * The script builds trees.
	 *
	 * (1) It reads the net.
	 *
	 * (2) It builds trees for the blocks.
	 *
	 * (3) It augments the event tree.
	 *
	 * (4) It copies the relevant block end trees into the ends of its descendant nets.
	 *
	 * (5) It copies the relevant augmented events into the event trees of its descendants.
	 *
	 */
	NetDisk *netDisk;
	Net *net;
	int32_t startTime;
	int32_t i, j;
	Chain *chain;
	Block *block;
	struct List *sortedChainAlignments;
	Group *group;
	Net *net2;
	End *end;
	End *end2;
	Cap *cap;
	Cap *cap2;
	Cap *cap3;
	Event *event;
	Net_EndIterator *endIterator;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * netDiskName = NULL;
	char * tempFileRootDirectory = NULL;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "tempDirRoot", required_argument, 0, 'e' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:c:e:h", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stringCopy(optarg);
				break;
			case 'c':
				netDiskName = stringCopy(optarg);
				break;
			case 'e':
				tempFileRootDirectory = stringCopy(optarg);
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

	assert(logLevelString == NULL || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
	assert(netDiskName != NULL);
	assert(tempFileRootDirectory != NULL);

	//////////////////////////////////////////////
	//Set up logging
	//////////////////////////////////////////////

	if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
		setLogLevel(LOGGING_INFO);
	}
	if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
		setLogLevel(LOGGING_DEBUG);
	}

	//////////////////////////////////////////////
	//Log (some of) the inputs
	//////////////////////////////////////////////

	logInfo("Net disk name : %s\n", netDiskName);
	logInfo("Temp file root directory : %s\n", tempFileRootDirectory);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	//////////////////////////////////////////////
	//For each net do tree building..
	//////////////////////////////////////////////

	for (j = optind; j < argc; j++) {
		const char *netName = argv[j];
		logInfo("Processing the net named: %s", netName);

		///////////////////////////////////////////////////////////////////////////
		// Parse the basic reconstruction problem
		///////////////////////////////////////////////////////////////////////////

		net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
		logInfo("Parsed the net to be refined\n");

		///////////////////////////////////////////////////////////////////////////
		//Setups the 'trees' for the caps for the top level problem.
		//In other words, it ensure ends from stubs have there root instance set - only needs doing for
		//the top level net, after which the stub roots will be set recursively.
		///////////////////////////////////////////////////////////////////////////

		if(net_getName(net) == 0) {
			endIterator = net_getEndIterator(net);
			while((end = net_getNextEnd(endIterator)) != NULL) {
				if(end_isStub(end)) {
					assert(end_getInstanceNumber(end) == 1);
					Event *rootEvent = eventTree_getRootEvent(net_getEventTree(net));
					Cap *rootCap = cap_construct(end, rootEvent);
					assert(event_isAncestor(cap_getEvent(end_getFirst(end)), rootEvent));
					cap_makeParentAndChild(rootCap, end_getFirst(end));
					end_setRootInstance(end, rootCap);
				}
				else {
					assert(!end_isCap(end));
					assert(end_isBlockEnd(end));
				}
			}
			net_destructEndIterator(endIterator);
		}

#ifdef BEN_DEBUG
		endIterator = net_getEndIterator(net);
		while((end = net_getNextEnd(endIterator)) != NULL) {
			if(!end_isBlockEnd(end)) {
				assert(end_getRootInstance(end) != NULL);
			}
		}
		net_destructEndIterator(endIterator);
#endif

		///////////////////////////////////////////////////////////////////////////
		//Construct the chain alignments for all the non-trivial chains.
		///////////////////////////////////////////////////////////////////////////

		startTime = time(NULL);
		sortedChainAlignments = constructEmptyList(0, (void (*)(void *))chainAlignment_destruct);
		Net_ChainIterator *chainIterator = net_getChainIterator(net);
		while((chain = net_getNextChain(chainIterator)) != NULL) {
			Block **blockChain = chain_getBlockChain(chain, &i);
			if(i > 0) {
				listAppend(sortedChainAlignments, chainAlignment_construct(blockChain, i));
			}
			free(blockChain);
		}
		net_destructChainIterator(chainIterator);
		logInfo("Constructed the block trees for the non-trivial chains in the net in: %i seconds\n", time(NULL) - startTime);

		///////////////////////////////////////////////////////////////////////////
		//Construct the chain alignment for each trivial chain.
		///////////////////////////////////////////////////////////////////////////

		startTime = time(NULL);
		Net_BlockIterator *blockIterator = net_getBlockIterator(net);
		while((block = net_getNextBlock(blockIterator)) != NULL) {
			if(block_getChain(block) == NULL) {
				listAppend(sortedChainAlignments, chainAlignment_construct(&block, 1));
			}
		}
		net_destructBlockIterator(blockIterator);

		qsort(sortedChainAlignments->list, sortedChainAlignments->length, sizeof(void *), (int (*)(const void *, const void *))chainAlignment_cmpFn);
#ifdef BEN_DEBUG
		for(i=1; i<sortedChainAlignments->length; i++) {
			assert(((ChainAlignment *)sortedChainAlignments->list[i-1])->totalAlignmentLength >=
					((ChainAlignment *)sortedChainAlignments->list[i])->totalAlignmentLength);
		}
#endif

		logInfo("Constructed the block trees for the trivial chains in the net in: %i seconds\n", time(NULL) - startTime);

		///////////////////////////////////////////////////////////////////////////
		//For each chain call the tree construction function.
		///////////////////////////////////////////////////////////////////////////

		startTime = time(NULL);
		buildChainTrees((ChainAlignment **)sortedChainAlignments->list, sortedChainAlignments->length, net_getEventTree(net), tempFileRootDirectory);
		destructList(sortedChainAlignments);

		logInfo("Augmented the block trees in: %i seconds\n", time(NULL) - startTime);

#ifdef BEN_DEBUG
		endIterator = net_getEndIterator(net);
		while((end = net_getNextEnd(endIterator)) != NULL) {
			assert(end_getRootInstance(end) != NULL);
		}
		net_destructEndIterator(endIterator);
#endif

		///////////////////////////////////////////////////////////////////////////
		//Pass the end trees and augmented events to the child nets.
		///////////////////////////////////////////////////////////////////////////

		startTime = time(NULL);
		Net_GroupIterator *groupIterator = net_getGroupIterator(net);
		while((group = net_getNextGroup(groupIterator)) != NULL) {
			net2 = group_getNestedNet(group);
			if(net2 != NULL) {
				//add in the end trees and augment the event trees.
				Group_EndIterator *endIterator = group_getEndIterator(group);
				while((end = group_getNextEnd(endIterator)) != NULL) {
					end2 = net_getEnd(net2, end_getName(end));
					assert(end2 != NULL);
					//copy the caps.
					End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
					while((cap = end_getNext(instanceIterator)) != NULL) {
						if(end_getInstance(end2, cap_getName(cap)) == NULL) {
							assert(cap_getChildNumber(cap) > 0); //can not be a leaf
							//make sure the augmented event is in there.
							event = cap_getEvent(cap);
							if(eventTree_getEvent(net_getEventTree(net2), event_getName(event)) == NULL) {
								assert(event_getChildNumber(event) == 1); //must be a unary event
								copyConstructUnaryEvent(event, net_getEventTree(net2));
							}
							event = eventTree_getEvent(net_getEventTree(net2), event_getName(event));
							assert(event != NULL);
							cap_copyConstruct(end2, cap);
						}
					}
					//now copy the parent links.
					while((cap = end_getPrevious(instanceIterator)) != NULL) {
						cap2 = end_getInstance(end2, cap_getName(cap));
						assert(cap2 != NULL);
						if(cap_getParent(cap) != NULL) {
							cap3 = end_getInstance(end2, cap_getName(cap_getParent(cap)));
							assert(cap3 != NULL);
							cap_makeParentAndChild(
									cap3,
									cap2);
						}
						else {
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
		net_destructGroupIterator(groupIterator);
		logInfo("Filled in end trees and augmented the event trees for the child nets in: %i seconds\n", time(NULL) - startTime);
	}

	///////////////////////////////////////////////////////////////////////////
	// (9) Write the net to disk.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	netDisk_write(netDisk);
	logInfo("Updated the net on disk in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	//(15) Clean up.
	///////////////////////////////////////////////////////////////////////////

	//Destruct stuff
	startTime = time(NULL);
	netDisk_destruct(netDisk);

	logInfo("Cleaned stuff up and am finished in: %i seconds\n", time(NULL) - startTime);
	return 0;
}
