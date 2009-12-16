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
	 * Structure to represent a concatenated list of atoms as a single 2d alignment.
	 */
	//Matrix of alignment, matrix[column][row]
	AtomInstance ***matrix; //NULL values are okay, where an instance of an atom is missing from an instance of the chain.
	Atom **atoms; //this list of atoms, in order.
	int32_t columnNumber; //the number of atoms.
	int32_t rowNumber; //the number of rows of the alignment, each row containing an instance of atoms in the chain.
	int32_t totalAlignmentLength; //the length in base pairs of the alignment.
} ChainAlignment;

char **chainAlignment_getAlignment(ChainAlignment *chainAlignment) {
	/*
	 * Constructs a concatenated 2d matrix of chars referenced by [atom][instance], representing the
	 * base pair alignment of the chain alignment.
	 */
	char **alignment;
	int32_t i, j, k, l;
	AtomInstance *atomInstance;
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
			atomInstance = chainAlignment->matrix[i][j];
			if(atomInstance == NULL) {
				for(k=0; k<atom_getLength(chainAlignment->atoms[i]); k++) {
					alignment[j][l++]  = 'N';
				}
			}
			else {
				cA = atomInstance_getString(atomInstance);
				for(k=0; k<atomInstance_getLength(atomInstance); k++) {
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

AtomInstance *chainAlignment_getFirstNonNullAtomInstance(ChainAlignment *chainAlignment, int32_t row, bool increasing) {
	/*
	 * Gets the first instance on an row in the chain alignment which is non null. If increasing is false, gets the last.
	 */
	AtomInstance *atomInstance;
	int32_t j = 0, k = chainAlignment->columnNumber, l = 1;
	if(!increasing) {
		j = chainAlignment->columnNumber-1;
		k = -1;
		l = -1;
	}
	for(; j!=k; j += l) {
		atomInstance = chainAlignment->matrix[j][row];
		if(atomInstance != NULL) {
			return atomInstance;
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
	AtomInstance *atomInstance;

	names = malloc(sizeof(Name) * chainAlignment->rowNumber);
	for(i=0; i<chainAlignment->rowNumber; i++) {
		atomInstance = chainAlignment_getFirstNonNullAtomInstance(chainAlignment, i, _5End);
		names[i] = end_getName(endInstance_getEnd(endInstance_getAdjacency(_5End ? atomInstance_get5End(atomInstance) : atomInstance_get3End(atomInstance))));
	}
	return names;
}

Name *chainAlignment_getLeafEventNames(ChainAlignment *chainAlignment) {
	/*
	 * Gets the names of the leaf events of the rows (instances), in the alignment.
	 */
	Name *names;
	int32_t i;
	AtomInstance *atomInstance;

	names = malloc(sizeof(Name) * chainAlignment->rowNumber);
	for(i=0; i<chainAlignment->rowNumber; i++) {
		atomInstance = chainAlignment_getFirstNonNullAtomInstance(chainAlignment, i, 1);
		names[i] = event_getName(atomInstance_getEvent(atomInstance));
	}
	return names;
}

int32_t *chainAlignment_getAtomBoundaries(ChainAlignment *chainAlignment) {
	/*
	 * Gets the boundaries of atom in the chain alignment.
	 */
	int32_t *atomBoundaries;
	int32_t i, j;

	atomBoundaries = malloc(sizeof(int32_t) * chainAlignment->columnNumber);
	j = 0;
	for(i=0; i<chainAlignment->columnNumber; i++) {
		atomBoundaries[i] = j + atom_getLength(chainAlignment->atoms[i]);
		j = atomBoundaries[i];
	}

	return atomBoundaries;
}

void chomp(const char *s) {
	char *p;
	while (NULL != s && NULL != (p = strrchr(s, '\n'))) {
		*p = '\0';
	}
}

void buildChainTrees_Bernard(int32_t atomNumber, char ***concatenatedAtoms, Name **_5Ends, Name **_3Ends, Name **leafEventLabels,
							int32_t **atomBoundaries, char *eventTreeString, const char *tempDir, ChainAlignment **chainAlignments,
							char **modifiedEventTreeString, char ****atomTreeStrings, int32_t ***refinedAtomBoundaries, int32_t **refinedAtomNumbers) {
	/*
	 * Here's the function you need to fill in, I haven't defined the outputs yet - you get it working and then we can discuss.
	 *
	 * Arrays are all indexed first by the chain/concatenated atoms.
	 * Alignments are column/row indexed with rows as chain instances and columns as aligned bases.
	 * So for example: concatenatedAtoms[i][j][k] is the ith chain/concatenated atom, jth row (chain instance), kth column (position in concatenated atom)
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
	char atomFileName[] = "annot";
	char mafDirFileName[] = "mafs";
	char augTreeFileName[] = "remap.augtree";
	char dupTreeFileName[] = "event.duptree";

	const int32_t TMP_BUFFER_SIZE = 256; // Assume that 256 chars is enough for any temp filename strings
	char tmpStringBuffer[TMP_BUFFER_SIZE];
	FILE *fp = NULL;

	const int32_t LINE_BUFF_SIZE = 4096; // Assume that 4096 chars is enough for reading in a single line of a file
	char lineBuffer[LINE_BUFF_SIZE];
	char *tmpString = NULL;

	int32_t **newAtomBoundaries = NULL;
	int32_t *newAtomNumbers = NULL;

	char ***atomTreeArray = NULL;

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

	/* Output the atom definition file */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, atomFileName);
	fp = fopen(tmpStringBuffer, "w");
	for (i=0; i<atomNumber; i++) {
		rowNumber = chainAlignments[i]->rowNumber;

		/* NOTE: Atom Numbers start at 1 */
		fprintf(fp, ">%d %d %d\n", i+1, chainAlignments[i]->totalAlignmentLength, rowNumber);
		for (j=0; j<rowNumber; j++) {
			fprintf(fp, "%s.chr0:1-%d + %s %s\n", netMisc_nameToString(leafEventLabels[i][j]), chainAlignments[i]->totalAlignmentLength, netMisc_nameToString(_5Ends[i][j]), netMisc_nameToString(_3Ends[i][j]));
		}
	}
	fclose(fp);

	/* Create the atom map file */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, "atom.map");
	fp = fopen(tmpStringBuffer, "w");
	for (i=0; i<atomNumber; i++) {
		/* NOTE: Atom Numbers start at 1 */
		fprintf(fp, "%d\t%d\n", i+1, i+1);
	}
	fclose(fp);

	/* Create the maf files */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, mafDirFileName);
	mkdir(tmpStringBuffer, 0777);

	for (i=0; i<atomNumber; i++) {
		/* NOTE: Atom Numbers start at 1 */
		snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s/%d", tempDirName, mafDirFileName, i+1);

		fp = fopen(tmpStringBuffer, "w");

		rowNumber = chainAlignments[i]->rowNumber;
		for (j=0; j<rowNumber; j++) {
			fprintf(fp, ">%s\n", netMisc_nameToString(leafEventLabels[i][j]));
			fprintf(fp, "%s\n", concatenatedAtoms[i][j]);
		}

		fclose(fp);
	}

	/* Run the tree pipeline */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "conTrees_PhyloBuilder.py %s", tempDirName);
	exitOnFailure(system(tmpStringBuffer), "conTrees_PhyloBuilder.py failed\n");
	printf("Completed running tree pipeline, now onto parsing\n");

	/* Clone Atom boundaries to refined Atom Boundaries*/
	newAtomBoundaries = malloc(sizeof(void *) * atomNumber);
	newAtomNumbers = malloc(sizeof(int32_t) * atomNumber);
	for (i=0; i<atomNumber; i++) {
		colNumber = chainAlignments[i]->columnNumber;
		newAtomBoundaries[i] = malloc(sizeof(int32_t) * colNumber);
		newAtomNumbers[i] = colNumber;
		for (j=0; j<colNumber; j++) {
			newAtomBoundaries[i][j] = atomBoundaries[i][j];
		}
	}
	*refinedAtomBoundaries = newAtomBoundaries;
	*refinedAtomNumbers = newAtomNumbers;

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
	atomTreeArray = malloc(sizeof(void *) * atomNumber);
	for (i=0; i<atomNumber; i++) {
		colNumber = chainAlignments[i]->columnNumber;
		atomTreeArray[i] = malloc(sizeof(void *) * colNumber);
	}
	i = 0; // Variable to store atomNumber
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirName, augTreeFileName);
	fp = fopen(tmpStringBuffer, "r");
	if (fp != NULL) {
		while(fgets(lineBuffer, LINE_BUFF_SIZE, fp) != NULL) {
			chomp(lineBuffer);
			if (lineBuffer[0] == '>') {
				sscanf(lineBuffer, ">%d", &i);
				i -= 1; // Atom number ids start at 1
			} else {
				colNumber = chainAlignments[i]->columnNumber;
				for (j=0; j<colNumber; j++) {
					tmpString = stringCopy(lineBuffer);
					atomTreeArray[i][j] = tmpString;
				}
			}
		}
	} else {
		perror("Failed to open augTree file");
	}
	fclose(fp);
	*atomTreeStrings = atomTreeArray;

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


AtomInstance *buildChainTrees3P(Atom *atom, AtomInstance **atomInstances, int32_t atomNumber,
		struct BinaryTree *binaryTree, struct hashtable *newEventNameMap) {
	/*
	 * Recursive partner to buildChainTree3 function, recurses on the binary tree constructing the atom tree.
	 * The labels of the leaves are indexes into the atom instances array, the internal node's labels are events in the event tree.
	 */
	if(binaryTree->internal) { //deal with an internal node of the atom tree.
		AtomInstance *leftInstance = buildChainTrees3P(atom, atomInstances, atomNumber, binaryTree->left, newEventNameMap);
		AtomInstance *rightInstance = buildChainTrees3P(atom, atomInstances, atomNumber, binaryTree->right, newEventNameMap);
		if(leftInstance != NULL) {
			if(rightInstance != NULL) {

				Event *event;
				if(isNewEvent(binaryTree->label)) {
					const char *cA = hashtable_search(newEventNameMap, binaryTree->label);
					assert(cA != NULL);
					event = eventTree_getEvent(net_getEventTree(atom_getNet(atom)), netMisc_stringToName(cA));
				}
				else {
					event = eventTree_getEvent(net_getEventTree(atom_getNet(atom)), netMisc_stringToName(binaryTree->label));
					if(event != NULL) {
						printBinaryTree(stderr, binaryTree);
					}
				}

				assert(event != NULL); //check event is present in the event tree.
				//Check that this does not create a cycle with respect to the event tree.
				assert(event_isAncestor(atomInstance_getEvent(leftInstance), event));
				assert(event_isAncestor(atomInstance_getEvent(rightInstance), event));

				AtomInstance *atomInstance = atomInstance_construct(atom, event);
				atomInstance_makeParentAndChild(atomInstance, leftInstance);
				atomInstance_makeParentAndChild(atomInstance, rightInstance);
				return atomInstance;
			}
			return leftInstance;
		}
		return rightInstance;
	}
	else { //a leaf, so find the leaf instance in the list of atom instances (this may be null if missing data).
		int32_t i;
		assert(sscanf(binaryTree->label, "%i", &i) == 1);
		assert(i < atomNumber);
		assert(i >= 0);
		return atomInstances[i];
	}
}

void buildChainTrees3(Atom *atom, AtomInstance **atomInstances, int32_t atomNumber, struct BinaryTree *binaryTree, struct hashtable *newEventNameMap) {
	/*
	 * Constructs an atom tree for the atom.
	 */
	AtomInstance *mostAncestralEvent = buildChainTrees3P(atom, atomInstances, atomNumber, binaryTree, newEventNameMap);
	assert(atom_getInstanceNumber(atom) > 0);
	assert(mostAncestralEvent != NULL);
	//Make a root event.
	Event *rootEvent = eventTree_getRootEvent(net_getEventTree(atom_getNet(atom)));
	AtomInstance *rootAtomInstance = atomInstance_construct(atom, rootEvent);
	assert(event_isAncestor(atomInstance_getEvent(mostAncestralEvent), rootEvent));
	atomInstance_makeParentAndChild(rootAtomInstance, mostAncestralEvent);
	atom_setRootInstance(atom, rootAtomInstance);

#ifdef BEN_DEBUG //Now go through all events checking they have a parent.
	Atom_InstanceIterator *instanceIterator = atom_getInstanceIterator(atom);
	AtomInstance *atomInstance;
	while((atomInstance = atom_getNext(instanceIterator)) != NULL) {
		AtomInstance *parent = atomInstance_getParent(atomInstance);
		if(parent == NULL) {
			assert(atomInstance == rootAtomInstance);
		}
	}
	atom_destructInstanceIterator(instanceIterator);
#endif
}

void buildChainTrees2(ChainAlignment *chainAlignment,
					  struct BinaryTree **refinedAtomTrees,
					  int32_t *refinedAtomBoundaries,
					  int32_t refinedAtomNumber,
					  struct hashtable *newEventNameMap) {
	/*
	 * Iterates through a chain alignment, constructing the atom trees and splitting atoms as needed.
	 */
	assert(chainAlignment->columnNumber <= refinedAtomNumber);
	int32_t i, j, k;
	Atom *atom, *leftAtom, *rightAtom;

	j=0;
	k=0;
	for(i=0; i < chainAlignment->columnNumber; i++) {
		atom = chainAlignment->atoms[i];
		k += atom_getLength(atom);
		//walk along the refined atom boundaries, splitting the considered as needed.
		assert(j < refinedAtomNumber && refinedAtomBoundaries[j] <= k);
		do {
			if(refinedAtomBoundaries[j] < k) { //we need to split the atom
				assert(refinedAtomBoundaries[j] >= k - atom_getLength(atom)); //boundary must break atom so that left atom is at least one base pair long.
				assert(refinedAtomBoundaries[j] < k); //boundary must break atom so that right atom is at least one base pair long.
				atom_split(atom, atom_getLength(atom) - (k - refinedAtomBoundaries[j]), &leftAtom, &rightAtom);
				buildChainTrees3(leftAtom, chainAlignment->matrix[i], chainAlignment->rowNumber, refinedAtomTrees[j], newEventNameMap);
				atom = rightAtom;
				assert(k - atom_getLength(atom) == refinedAtomBoundaries[j]); //check the split did what we expect
				assert(j+1 < refinedAtomNumber && refinedAtomBoundaries[j+1] <= k); //check that we have another atom tree to deal with the right side of the split.
			}
			else {
				buildChainTrees3(atom, chainAlignment->matrix[i], chainAlignment->rowNumber, refinedAtomTrees[j], newEventNameMap);
			}
			j++;
		} while(j < refinedAtomNumber && refinedAtomBoundaries[j] <= k);
	}
	assert(j == refinedAtomNumber);
}

void buildChainTrees(ChainAlignment **chainAlignments, int32_t chainAlignmentNumber, EventTree *eventTree,
		const char *tempFilePath) {
	/*
	 * This function builds a load of inputs which are then passed to Bernard's code.
	 */
	int32_t i, j;

	//Make the atom alignments.
	char ***concatenatedAtoms = malloc(sizeof(void *) * chainAlignmentNumber); //this is the list of 2d alignments.
	Name **_5Ends = malloc(sizeof(void *) * chainAlignmentNumber); //these are the lists of ends associated with each end.
	Name **_3Ends = malloc(sizeof(void *) * chainAlignmentNumber);
	Name **leafEventLabels = malloc(sizeof(void *) * chainAlignmentNumber); //this is the list of leaf event labels.
	int32_t **atomBoundaries = malloc(sizeof(void *) * chainAlignmentNumber); //each chain alignment has a list of atom lengths to demark where the atom boundaries are.

	//now fill out the the various arrays.
	for(i=0; i<chainAlignmentNumber; i++) {
		concatenatedAtoms[i] = chainAlignment_getAlignment(chainAlignments[i]);
		_5Ends[i] = chainAlignment_getEndNames(chainAlignments[i], 1);
		_3Ends[i] = chainAlignment_getEndNames(chainAlignments[i], 0);
		leafEventLabels[i] = chainAlignment_getLeafEventNames(chainAlignments[i]);
		atomBoundaries[i] = chainAlignment_getAtomBoundaries(chainAlignments[i]);
	}
	//Event tree string
	char *eventTreeString = eventTree_makeNewickString(eventTree);

	//Construct the random dir.
	char *randomDir;
	exitOnFailure(constructRandomDir(tempFilePath, &randomDir), "Tried to make a recursive directory of temp files but failed\n");

	//call to Bernard's code
	char *augmentedEventTreeString; //pointer to string holding the augmented event tree.
	char ***atomTreeStrings; //array of string pointers, for holding the constructed atom trees.
	int32_t **refinedAtomBoundaries; //like the atom boundaries, but revised by allowing for splits in the existing atoms.
	int32_t *refinedAtomNumbers; //the lengths of the atom boundary arrays.
	buildChainTrees_Bernard(chainAlignmentNumber, concatenatedAtoms, _5Ends, _3Ends, leafEventLabels,
							atomBoundaries, eventTreeString, randomDir, chainAlignments,
							&augmentedEventTreeString, &atomTreeStrings, &refinedAtomBoundaries, &refinedAtomNumbers);
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
	 * Now process each new atom tree.
	 */
	for(i=0; i<chainAlignmentNumber; i++) {
		struct BinaryTree **atomTrees = malloc(sizeof(void *)* refinedAtomNumbers[i]);
		for(j=0; j<refinedAtomNumbers[i]; j++) {
			atomTrees[j] = newickTreeParser(atomTreeStrings[i][j], 0.0, 0);
		}
		buildChainTrees2(chainAlignments[i], atomTrees, refinedAtomBoundaries[i], refinedAtomNumbers[i], newEventNameMap);
		for(j=0; j<refinedAtomNumbers[i]; j++) {
			destructBinaryTree(atomTrees[j]);
		}
		free(atomTrees);
	}
	logDebug("Processed the new atom trees\n");

	/*
	 * Cleanup the inputs.
	 */
	for(i=0; i<chainAlignmentNumber; i++) {
		free(concatenatedAtoms[i]);
		free(_5Ends[i]);
		free(_3Ends[i]);
		free(leafEventLabels[i]);
		free(atomBoundaries[i]);
	}
	free(concatenatedAtoms);
	free(_5Ends);
	free(_3Ends);
	free(leafEventLabels);
	free(atomBoundaries);
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

ChainAlignment *chainAlignment_construct(Atom **atoms, int32_t atomsLength) {
	/*
	 * Constructs a chain alignment structure from a chain of atoms.
	 */
	int32_t i, j, k;
	Atom *atom;
	AtomInstance *atomInstance;
	AtomInstance *atomInstance2;
	EndInstance *endInstance;
	ChainAlignment *chainAlignment;
	struct hashtable *hash;
	struct List *list;
	struct List *list2;
	struct avl_table *avlTable = avl_create(NULL, NULL, NULL);
	avl_destroy(avlTable, NULL);

	/*
	 * First iterate through all the atom instances in the chain, in order, to construct instances of the chain.
	 */
	hash = create_hashtable(1, hashtable_key, hashtable_equalKey, NULL, NULL); //to keep track of the instances included in a chain.
	list = constructEmptyList(0, (void (*)(void *))destructList); //the list of chain instances.
	k = 0;
	assert(atomsLength > 0);
	for(i=0; i<atomsLength; i++) {
		atom = atoms[i];
		Atom_InstanceIterator *instanceIterator = atom_getInstanceIterator(atom);
		while((atomInstance = atom_getNext(instanceIterator)) != NULL) {
			k++;
			assert(atomInstance_getOrientation(atomInstance));
			if(hashtable_search(hash, atomInstance) == NULL) { //not yet in a chain instance
				list2 = constructEmptyList(atomsLength, NULL);
				for(j=0; j<atomsLength; j++) { //this list will contain one instance for each atom, or NULL, of missing.
					list2->list[j] = NULL;
				}
				listAppend(list, list2);
				j = i; //start from the atom we're up to.
				while(1) {
					assert(atomInstance_getOrientation(atomInstance));
					hashtable_insert(hash, atomInstance, atomInstance); //put in the hash to say we've seen it.
					list2->list[j++] = atomInstance; //put in the list of chains list.
					if(j == atomsLength) { //end of chain
						break;
					}
					endInstance = endInstance_getAdjacency(atomInstance_get3End(atomInstance));
					assert(!end_isCap(endInstance_getEnd(endInstance))); //can not be a cap.
					if(end_isStub(endInstance_getEnd(endInstance))) { //terminates with missing information.
						break;
					}
					atomInstance2 = endInstance_getAtomInstance(endInstance);
					assert(atomInstance != NULL); //must be connected to an atom instance (not a cap).
					assert(atomInstance_getAtom(atomInstance2) == atoms[j] || atomInstance_getAtom(atomInstance2) == atom_getReverse(atoms[j-1])); //must be connected to reverse of itself or the next atom
					if(atomInstance_getAtom(atomInstance2) != atoms[j]) {
						break;
					}
					atomInstance = atomInstance2;
				}
			}
			else {
				assert(i != 0);
			}
		}
		atom_destructInstanceIterator(instanceIterator);
	}
	assert(k == (int32_t)hashtable_count(hash));
	assert(list->length > 0);

	/*
	 * Now convert the chain instances in the list into the desired format for the chain alignment.
	 */

	//alloc the chain alignment and the matrix of instances.
	chainAlignment = malloc(sizeof(ChainAlignment));
	chainAlignment->matrix = malloc(sizeof(AtomInstance **) * atomsLength);
	for(i=0; i<atomsLength; i++) {
		chainAlignment->matrix[i] = malloc(sizeof(AtomInstance *) * list->length);
	}
	//fill out the fields of the chain alignment, including the matrix.
	chainAlignment->columnNumber = atomsLength;
	chainAlignment->rowNumber = list->length;
	for(i=0; i<list->length; i++) {
		list2 = list->list[i];
		assert(list2->length == atomsLength);
		for(j=0; j<atomsLength; j++) {
			chainAlignment->matrix[j][i] = list2->list[j];
		}
	}
	//Calculate the total length.
	chainAlignment->totalAlignmentLength = 0;
	for(i=0;i<atomsLength; i++) {
		chainAlignment->totalAlignmentLength += atom_getLength(atoms[i]);
	}
	assert(chainAlignment->totalAlignmentLength > 0);
	//Add chain of atoms.
	chainAlignment->atoms = malloc(sizeof(void *)*atomsLength);
	for(i=0; i<atomsLength; i++) {
		chainAlignment->atoms[i] = atoms[i];
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
	assert(event_getBranchLength(childEvent) > branchLength);

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
	 * (2) It builds trees for the atoms.
	 *
	 * (3) It augments the event tree.
	 *
	 * (4) It copies the relevant atom end trees into the ends of its descendant nets.
	 *
	 * (5) It copies the relevant augmented events into the event trees of its descendants.
	 *
	 */
	NetDisk *netDisk;
	Net *net;
	int32_t startTime;
	int32_t i, j;
	Chain *chain;
	Atom *atom;
	struct List *sortedChainAlignments;
	AdjacencyComponent *adjacencyComponent;
	Net *net2;
	End *end;
	End *end2;
	EndInstance *endInstance;
	EndInstance *endInstance2;
	EndInstance *endInstance3;
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
					EndInstance *rootEndInstance = endInstance_construct(end, rootEvent);
					assert(event_isAncestor(endInstance_getEvent(end_getFirst(end)), rootEvent));
					endInstance_makeParentAndChild(rootEndInstance, end_getFirst(end));
					end_setRootInstance(end, rootEndInstance);
				}
				else {
					assert(!end_isCap(end));
					assert(end_isAtomEnd(end));
				}
			}
			net_destructEndIterator(endIterator);
		}

#ifdef BEN_DEBUG
		endIterator = net_getEndIterator(net);
		while((end = net_getNextEnd(endIterator)) != NULL) {
			if(!end_isAtomEnd(end)) {
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
			Atom **atomChain = chain_getAtomChain(chain, &i);
			if(i > 0) {
				listAppend(sortedChainAlignments, chainAlignment_construct(atomChain, i));
			}
			free(atomChain);
		}
		net_destructChainIterator(chainIterator);
		logInfo("Constructed the atom trees for the non-trivial chains in the net in: %i seconds\n", time(NULL) - startTime);

		///////////////////////////////////////////////////////////////////////////
		//Construct the chain alignment for each trivial chain.
		///////////////////////////////////////////////////////////////////////////

		startTime = time(NULL);
		Net_AtomIterator *atomIterator = net_getAtomIterator(net);
		while((atom = net_getNextAtom(atomIterator)) != NULL) {
			if(atom_getChain(atom) == NULL) {
				listAppend(sortedChainAlignments, chainAlignment_construct(&atom, 1));
			}
		}
		net_destructAtomIterator(atomIterator);

		qsort(sortedChainAlignments->list, sortedChainAlignments->length, sizeof(void *), (int (*)(const void *, const void *))chainAlignment_cmpFn);
#ifdef BEN_DEBUG
		for(i=1; i<sortedChainAlignments->length; i++) {
			assert(((ChainAlignment *)sortedChainAlignments->list[i-1])->totalAlignmentLength >=
					((ChainAlignment *)sortedChainAlignments->list[i])->totalAlignmentLength);
		}
#endif

		logInfo("Constructed the atom trees for the trivial chains in the net in: %i seconds\n", time(NULL) - startTime);

		///////////////////////////////////////////////////////////////////////////
		//For each chain call the tree construction function.
		///////////////////////////////////////////////////////////////////////////

		startTime = time(NULL);
		buildChainTrees((ChainAlignment **)sortedChainAlignments->list, sortedChainAlignments->length, net_getEventTree(net), tempFileRootDirectory);
		destructList(sortedChainAlignments);

		logInfo("Augmented the atom trees in: %i seconds\n", time(NULL) - startTime);

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
		Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
		while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
			net2 = adjacencyComponent_getNestedNet(adjacencyComponent);
			//add in the end trees and augment the event trees.
			AdjacencyComponent_EndIterator *endIterator = adjacencyComponent_getEndIterator(adjacencyComponent);
			while((end = adjacencyComponent_getNextEnd(endIterator)) != NULL) {
				end2 = net_getEnd(net2, end_getName(end));
				assert(end2 != NULL);
				//copy the end instances.
				End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
				while((endInstance = end_getNext(instanceIterator)) != NULL) {
					if(end_getInstance(end2, endInstance_getName(endInstance)) == NULL) {
						assert(endInstance_getChildNumber(endInstance) > 0); //can not be a leaf
						//make sure the augmented event is in there.
						event = endInstance_getEvent(endInstance);
						if(eventTree_getEvent(net_getEventTree(net2), event_getName(event)) == NULL) {
							assert(event_getChildNumber(event) == 1); //must be a unary event
							copyConstructUnaryEvent(event, net_getEventTree(net2));
						}
						event = eventTree_getEvent(net_getEventTree(net2), event_getName(event));
						assert(event != NULL);
						endInstance_copyConstruct(end2, endInstance);
					}
				}
				//now copy the parent links.
				while((endInstance = end_getPrevious(instanceIterator)) != NULL) {
					endInstance2 = end_getInstance(end2, endInstance_getName(endInstance));
					assert(endInstance2 != NULL);
					if(endInstance_getParent(endInstance) != NULL) {
						endInstance3 = end_getInstance(end2, endInstance_getName(endInstance_getParent(endInstance)));
						assert(endInstance3 != NULL);
						endInstance_makeParentAndChild(
								endInstance3,
								endInstance2);
					}
					else {
						assert(end_getRootInstance(end) != NULL);
						assert(endInstance == end_getRootInstance(end));

						assert(end_getRootInstance(end2) == NULL);
						end_setRootInstance(end2, endInstance2);
					}
				}
				end_destructInstanceIterator(instanceIterator);
			}
			adjacencyComponent_destructEndIterator(endIterator);
		}
		net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

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
