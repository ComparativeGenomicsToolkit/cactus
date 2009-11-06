#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"
#include "bioioC.h"

typedef struct _chainAlignment {
	//Matrix of alignment, matrix[column][row]
	AtomInstance ***matrix;
	Atom **atoms;
	int32_t columnNumber;
	int32_t rowNumber;
	int32_t totalAlignmentLength;
} ChainAlignment;

char **chainAlignment_getAlignment(ChainAlignment *chainAlignment) {
	char **alignment;
	int32_t i, j, k, l;
	AtomInstance *atomInstance;
	char *cA;

	alignment = malloc(sizeof(void *) * chainAlignment->rowNumber);
	for(i=0; i<chainAlignment->rowNumber; i++) {
		alignment[i] = malloc(sizeof(char) * (chainAlignment->totalAlignmentLength + 1));
	}

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
					alignment[j][l++]  = cA[k];
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

void buildChainTrees_Bernard(int32_t atomNumber, char ***concatenatedAtoms, Name **_5Ends, Name **_3Ends, Name **leafEventLabels,
							int32_t **atomBoundaries, char *eventTreeString, const char *tempDir, ChainAlignment **chainAlignments,
							char **modifiedEventTreeString, char ****atomTreeStrings, int32_t ***refinedAtomBoundaries, int32_t **refinedAtomNumbers) {
	/*
	 * Here's the function you need to fill in, I haven't defined the outputs yet - you get it working and then we can discuss.
	 *
	 * Arrays are all indexed first by the chain/concatenated atoms.
	 * Alignments are column/row indexed with rows as chain instances and columns as aligned bases.
	 * So for example: concatenatedAtoms[i][j][k] is the ith chain/concatenated atom, jth column (position in concatenated atom), kth row (chain instance).
	 *
	 * Names can be converted to strings with: netMisc_nameToString() and netMisc_nameToStringStatic() (See the API).
	 */

	int32_t i = 0;
	int32_t j = 0;
	int32_t atomIdx = 0;

	char *randomDirName = NULL;

	int32_t rowNumber = 0;

	char eventTreeFileName[] = "pre.event.tree";
	char atomFileName[] = "annot";
	char mafDirFileName[] = "mafs";
	char augTreeFileName[] = "remap.augtree";
	char dupTreeFileName[] = "R.duptree";

	const int32_t TMP_BUFFER_SIZE = 256; // Assume that 256 chars is enough for any temp filename strings
	char tmpStringBuffer[TMP_BUFFER_SIZE];
	FILE *fp = NULL;

	randomDirName = strrchr(tempDir, '/');
	if (randomDirName == NULL) {
		; // Assume this doesn't happen
	} else {
		randomDirName += 1;
	}

//	char tempDirDir[] = "/Users/bsuh/TOKYO";

	/* Output the eventTree to file */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDir, eventTreeFileName);
//	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirDir, eventTreeFileName);
	fp = fopen(tmpStringBuffer, "w");
	fprintf(fp, "%s\n", eventTreeString);
	fclose(fp);

	/* Output the atom definition file */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDir, atomFileName);
//	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirDir, atomFileName);
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

	/* Atom boundaries */
	for (i=0; i<atomNumber; i++) {
		rowNumber = chainAlignments[i]->rowNumber;
		for (j=0; j<rowNumber; j++) {
			atomIdx = atomBoundaries[i][j];
		}
	}

	/* Atom map */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDir, "atom.map");
//	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirDir, "atom.map");
	fp = fopen(tmpStringBuffer, "w");
	for (i=0; i<atomNumber; i++) {
		/* NOTE: Atom Numbers start at 1 */
		fprintf(fp, "%d\t%d\n", i+1, i+1);
	}
	fclose(fp);

	/* Create the maf files */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDir, mafDirFileName);
//	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDirDir, mafDirFileName);
	mkdir(tmpStringBuffer, 0777);

	for (i=0; i<atomNumber; i++) {
		/* NOTE: Atom Numbers start at 1 */
		snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s/%d", tempDir, mafDirFileName, i+1);
//		snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s/%d", tempDirDir, mafDirFileName, i+1);

		fp = fopen(tmpStringBuffer, "w");

		rowNumber = chainAlignments[i]->rowNumber;
		for (j=0; j<rowNumber; j++) {
			fprintf(fp, ">%s\n", netMisc_nameToString(leafEventLabels[i][j]));
			fprintf(fp, "%s\n", concatenatedAtoms[i][j]);
		}

		fclose(fp);
	}

	/* Run the tree pipeline */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "conTrees_PhyloBuilder.py %s", tempDir);
	exitOnFailure(system(tmpStringBuffer), "conTrees_PhyloBuilder.py failed\n");

	/* Readin tree pipeline output */
	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDir, dupTreeFileName);
	fp = fopen(tmpStringBuffer, "r");
	fclose(fp);

	snprintf(tmpStringBuffer, TMP_BUFFER_SIZE, "%s/%s", tempDir, augTreeFileName);
	fp = fopen(tmpStringBuffer, "r");
	fclose(fp);
}

void augmentEventTree(struct BinaryTree *modifiedEventTree, EventTree *eventTree) {

}

void buildAtomTree(struct BinaryTree *atomTree, Net *net, Atom *atom) {

}

void buildChainTrees(ChainAlignment **chainAlignments, int32_t chainAlignmentNumber, EventTree *eventTree,
		const char *tempFilePath) {
	/*
	 * This is the function to do the tree construction in.
	 */
	int32_t i, j;

	char ***concatenatedAtoms;
	Name **_5Ends;
	Name **_3Ends;
	Name **leafEventLabels;
	int32_t **atomBoundaries;
	char *eventTreeString;
	char *randomDir;

	//Make the atom alignments.
	concatenatedAtoms = malloc(sizeof(void *) * chainAlignmentNumber);
	_5Ends = malloc(sizeof(void *) * chainAlignmentNumber);
	_3Ends = malloc(sizeof(void *) * chainAlignmentNumber);
	leafEventLabels = malloc(sizeof(void *) * chainAlignmentNumber);
	atomBoundaries = malloc(sizeof(void *) * chainAlignmentNumber);
	for(i=0; i<chainAlignmentNumber; i++) {
		concatenatedAtoms[i] = chainAlignment_getAlignment(chainAlignments[i]);
		_5Ends[i] = chainAlignment_getEndNames(chainAlignments[i], 1);
		_3Ends[i] = chainAlignment_getEndNames(chainAlignments[i], 0);
		leafEventLabels[i] = chainAlignment_getLeafEventNames(chainAlignments[i]);
		atomBoundaries[i] = chainAlignment_getAtomBoundaries(chainAlignments[i]);
	}
	//Event tree string
	eventTreeString = eventTree_makeNewickString(eventTree);

	exitOnFailure(constructRandomDir(tempFilePath, &randomDir), "Tried to make a recursive directory of temp files but failed\n");

	//call to Bernard's code
	char *modifiedEventTreeString;
	char ***atomTreeStrings;
	int32_t **refinedAtomBoundaries;
	int32_t *refinedAtomNumbers;
	buildChainTrees_Bernard(chainAlignmentNumber, concatenatedAtoms, _5Ends, _3Ends, leafEventLabels,
							atomBoundaries, eventTreeString, randomDir, chainAlignments,
							&modifiedEventTreeString, &atomTreeStrings, &refinedAtomBoundaries, &refinedAtomNumbers);

	/*
	 * Augment the event tree with the new events.
	 */
	struct BinaryTree *modifiedEventTree = newickTreeParser(modifiedEventTreeString, 0.0);
	augmentEventTree(modifiedEventTree, eventTree);

	/*
	 * Now process each new atom tree.
	 */
	for(i=0; i<chainAlignmentNumber; i++) {
		for(j=0; j<refinedAtomNumbers[i]; j++) {
			newickTreeParser(atomTreeStrings[i][j], 0.0);
		}
	}

	exitOnFailure(destructRandomDir(randomDir), "Tried to destroy a recursive directory of temp files but failed\n");

	//build trees, augmented event trees.

	//clean up.
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

	//done!
}

int32_t chainAlignment_cmpFn(ChainAlignment *cA1, ChainAlignment *cA2, const void *a) {
	assert(a == NULL);
	return cA2->totalAlignmentLength - cA1->totalAlignmentLength; //sort descending
}

int chainAlignment_constructP(AtomInstance **atomInstance1, AtomInstance **atomInstance2) {
	assert(atomInstance_getStrand(*atomInstance1) && atomInstance_getStrand(*atomInstance2));
	Sequence *sequence1 = atomInstance_getSequence(*atomInstance1);
	Sequence *sequence2 = atomInstance_getSequence(*atomInstance2);
	int32_t i = netMisc_nameCompare(sequence_getName(sequence1), sequence_getName(sequence2));
	if(i == 0) {
		int32_t j = atomInstance_getStart(*atomInstance1);
		int32_t k = atomInstance_getStart(*atomInstance2);
		i = j - k;
		assert(i != 0);
	}
	return i;
}

ChainAlignment *chainAlignment_construct(Atom **atoms, int32_t atomsLength) {
	int32_t i, j, k;
	Atom *atom;
	AtomInstance *atomInstance;
	AtomInstance *atomInstance2;
	EndInstance *endInstance;
	ChainAlignment *chainAlignment;
	struct hashtable *hash;
	struct List *list;
	struct List *list2;

	//Construct a list of chain instances.
	hash = create_hashtable(1, hashtable_key, hashtable_equalKey, NULL, NULL);
	list = constructEmptyList(0, (void (*)(void *))destructList);
	k = 0;
	for(i=0; i<atomsLength; i++) {
		atom = atoms[i];
		Atom_InstanceIterator *instanceIterator = atom_getInstanceIterator(atom);
		while((atomInstance = atom_getNext(instanceIterator)) != NULL) {
			k++;
			assert(atomInstance_getOrientation(atomInstance));
			if(hashtable_search(hash, atomInstance) == NULL) {
				list2 = constructEmptyList(atomsLength, NULL);
				for(j=0; j<atomsLength; j++) {
					list2->list[j] = NULL;
				}
				listAppend(list, list2);
				j = i;
				while(1) {
					assert(atomInstance_getOrientation(atomInstance));
					hashtable_insert(hash, atomInstance, atomInstance);
					list2->list[j++] = atomInstance;
					if(j == atomsLength) { //end of chain
						break;
					}
					endInstance = endInstance_getAdjacency(atomInstance_get3End(atomInstance));
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
		}
		atom_destructInstanceIterator(instanceIterator);
	}
	assert(k == (int32_t)hashtable_count(hash));

	//Now do actual construction;
	chainAlignment = malloc(sizeof(ChainAlignment));
	chainAlignment->matrix = malloc(sizeof(AtomInstance **) * atomsLength);
	for(i=0; i<atomsLength; i++) {
		chainAlignment->matrix[i] = malloc(sizeof(AtomInstance *) * list->length);
	}
	chainAlignment->columnNumber = atomsLength;
	chainAlignment->rowNumber = list->length;
	for(i=0; i<list->length; i++) {
		list2 = list->list[i];
		assert(list2->length == atomsLength);
		for(j=0; j<list2->length; j++) {
			chainAlignment->matrix[j][i] = list2->list[j];
		}
	}
	//Calculate the total length.
	chainAlignment->totalAlignmentLength = 0;
	for(i=0;i<atomsLength; i++) {
		chainAlignment->totalAlignmentLength += atom_getLength(atoms[i]);
	}
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
	int32_t i;
	for(i=0; i<chainAlignment->columnNumber; i++) {
		free(chainAlignment->matrix[i]);
	}
	free(chainAlignment->matrix);
	free(chainAlignment);
}

Event *copyConstructUnaryEvent(Event *event, EventTree *eventTree2) {
	assert(event_getChildNumber(event) == 1);

	double branchLength = event_getBranchLength(event);

	Event *parentEvent = event_getParent(event);
	while(eventTree_getEvent(eventTree2, event_getName(parentEvent)) == NULL) {
		branchLength += event_getBranchLength(parentEvent);
		parentEvent = event_getParent(parentEvent);
		assert(parentEvent != NULL);
	}
	parentEvent = eventTree_getEvent(eventTree2, event_getName(parentEvent));

	Event *childEvent = event_getChild(event, 0);
	while(eventTree_getEvent(eventTree2, event_getName(parentEvent)) == NULL) {
		assert(event_getChildNumber(childEvent) == 1);
		childEvent = event_getChild(parentEvent, 0);
		assert(childEvent != NULL);
	}
	childEvent = eventTree_getEvent(eventTree2, event_getName(childEvent));

	return event_construct2(event_getMetaEvent(event), branchLength, parentEvent, childEvent, eventTree2);
}

void usage() {
	fprintf(stderr, "cactus_tree, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
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
	int32_t i;
	Chain *chain;
	Atom *atom;
	struct avl_table *sortedChainAlignments;
	ChainAlignment *chainAlignment;
	struct List *list;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * netDiskName = NULL;
	char * netName = NULL;
	char * tempFileRootDirectory = NULL;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "netName", required_argument, 0, 'd' },
			{ "tempDirRoot", required_argument, 0, 'e' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:c:d:e:h", long_options, &option_index);

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
			case 'd':
				netName = stringCopy(optarg);
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
	assert(netName != NULL);
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
	logInfo("Net name : %s\n", netName);
	logInfo("Temp file root directory : %s\n", tempFileRootDirectory);

	//////////////////////////////////////////////
	//Set up the temp file root directory
	//////////////////////////////////////////////

	//initialiseTempFileTree(tempFileRootDirectory, 100, 4);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
	logInfo("Parsed the net to be refined\n");

	///////////////////////////////////////////////////////////////////////////
	//Construct the chain alignments for all the non-trivial chains.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	sortedChainAlignments = avl_create((int32_t (*)(const void *, const void *, void *))chainAlignment_cmpFn, NULL, NULL);
	Net_ChainIterator *chainIterator = net_getChainIterator(net);
	while((chain = net_getNextChain(chainIterator)) != NULL) {
		Atom **atomChain = chain_getAtomChain(chain, &i);
		avl_insert(sortedChainAlignments, chainAlignment_construct(atomChain, i));
		free(atomChain);
	}
	net_destructChainIterator(chainIterator);
	logInfo("Constructed the atom trees for the non-trivial chains in the net in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	//For each lone atom call the tree construction function.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	Net_AtomIterator *atomIterator = net_getAtomIterator(net);
	while((atom = net_getNextAtom(atomIterator)) != NULL) {
		if(atom_getChain(atom) == NULL) {
			avl_insert(sortedChainAlignments, chainAlignment_construct(&atom, 1));
		}
	}
	net_destructAtomIterator(atomIterator);

	logInfo("Constructed the atom trees for the trivial chains in the net in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	//For each chain call the tree construction function.
	///////////////////////////////////////////////////////////////////////////

	startTime = time(NULL);
	struct avl_traverser *chainAlignmentIterator = mallocLocal(sizeof(struct avl_traverser));
	avl_t_init(chainAlignmentIterator, sortedChainAlignments);
	list = constructEmptyList(0, NULL);
	i = INT32_MAX;
	while((chainAlignment = avl_t_next(chainAlignmentIterator)) != NULL) {
		assert(chainAlignment->totalAlignmentLength <= i);
		i = chainAlignment->totalAlignmentLength;
		listAppend(list, chainAlignment);
	}
	buildChainTrees((ChainAlignment **)list->list, list->length, net_getEventTree(net), tempFileRootDirectory);
	destructList(list);

	logInfo("Augmented the atom trees in: %i seconds\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	//Pass the end trees and augmented events to the child nets.
	///////////////////////////////////////////////////////////////////////////

	/*startTime = time(NULL);
	Net_AdjacencyComponentIterator *adjacencyComponentIterator = net_getAdjacencyComponentIterator(net);
	while((adjacencyComponent = net_getNextAdjacencyComponent(adjacencyComponentIterator)) != NULL) {
		net2 = adjacencyComponent_getNestedNet(adjacencyComponent);
		//add in the end trees and augment the event trees.
		AdjacencyComponent_EndIterator *endIterator = adjacencyComponent_getEndIterator(adjacencyComponent);
		while((end = adjacencyComponent_getNextEnd(endIterator)) != NULL) {
			end2 = net_getEnd(net2, end_getName(end));
			//copy the end instances.
			End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
			while((endInstance = end_getNext(instanceIterator)) != NULL) {
				if(end_getInstance(end2, endInstance_getName(endInstance)) == NULL) {
					//make sure the augmented event is in there.
					if(eventTree_getEvent(net_getEventTree(net2), event_getName(endInstance_getEvent(endInstance))) == NULL) {
						copyConstructUnaryEvent(endInstance_getEvent(endInstance), net_getEventTree(net2));
					}
					endInstance_construct(end, endInstance_getEvent(endInstance));
				}
			}
			//now copy the parent links.
			while((endInstance = end_getPrevious(instanceIterator)) != NULL) {
				if(endInstance_getParent(endInstance) != NULL) {
					endInstance_makeParentAndChild(
							end_getInstance(end2, endInstance_getName(endInstance_getParent(endInstance))),
							end_getInstance(end2, endInstance_getName(endInstance)));
				}
			}
			end_destructInstanceIterator(instanceIterator);
		}
		adjacencyComponent_destructEndIterator(endIterator);
	}
	net_destructAdjacencyComponentIterator(adjacencyComponentIterator);

	logInfo("Filled in end trees and augmented the event trees for the child nets in: %i seconds\n", time(NULL) - startTime);*/

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
