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

	//do allocation
	alignment = malloc(sizeof(void *) * chainAlignment->totalAlignmentLength);
	for(i=0; i<chainAlignment->totalAlignmentLength; i++) {
		alignment[i] = malloc(sizeof(char) *chainAlignment->rowNumber);
	}

	//fill in the array.
	for(j=0; j<chainAlignment->rowNumber; j++) {
		l = 0;
		for(i=0; i<chainAlignment->columnNumber; i++) {
			atomInstance = chainAlignment->matrix[i][j];
			if(atomInstance == NULL) {
				for(k=0; k<atom_getLength(chainAlignment->atoms[i]); k++) {
					alignment[l++][j]  = 'N';
				}
			}
			else {
				cA = atomInstance_getString(atomInstance);
				for(k=0; k<atomInstance_getLength(atomInstance); k++) {
					alignment[l++][j]  = cA[k];
				}
				free(cA);
			}
		}
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
							int32_t **atomBoundaries, char *eventTreeString, const char *tempDir) {
	/*
	 * Here's the function you need to fill in, I haven't defined the outputs yet - you get it working and then we can discuss.
	 *
	 * Arrays are all indexed first by the chain/concatenated atoms.
	 * Alignments are column/row indexed with rows as chain instances and columns as aligned bases.
	 * So for example: concatenatedAtoms[i][j][k] is the ith chain/concatenated atom, jth column (position in concatenated atom), kth row (chain instance).
	 *
	 * Names can be converted to strings with: netMisc_nameToString() and netMisc_nameToStringStatic() (See the API).
	 */


}

void buildChainTrees(ChainAlignment **chainAlignments, int32_t chainAlignmentNumber, EventTree *eventTree,
		const char *tempFilePath) {
	/*
	 * This is the function to do the tree construction in.
	 */
	int32_t i;

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
	buildChainTrees_Bernard(chainAlignmentNumber, concatenatedAtoms, _5Ends, _3Ends, leafEventLabels, atomBoundaries, eventTreeString, randomDir);

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

ChainAlignment *chainAlignment_construct(struct List *atoms) {
	int32_t i, j;
	Atom *atom;
	Atom *atom2;
	AtomInstance *atomInstance;
	AtomInstance *atomInstance2;
	ChainAlignment *chainAlignment;
	struct hashtable *hash;
	struct List *list;
	struct List *list2;

	//Construct a list of chain instances.
	hash = create_hashtable(1, hashtable_key, hashtable_equalKey, NULL, NULL);
	list = constructEmptyList(0, (void (*)(void *))destructList);
	for(i=0; i<atoms->length; i++) {
		atom = atoms->list[i];
		Atom_InstanceIterator *instanceIterator = atom_getInstanceIterator(atom);
		while((atomInstance = atom_getNext(instanceIterator)) != NULL) {
			if(hashtable_search(hash, atomInstance) == NULL) {
				list2 = constructEmptyList(atoms->length, NULL);
				for(j=0; j<atoms->length; j++) {
					list2->list[j] = NULL;
				}
				listAppend(list, list2);
				j = i;
				while(1) {
					hashtable_insert(hash, atomInstance, atomInstance);
					list2->list[j++] = atomInstance;
					atomInstance2 = endInstance_getAtomInstance(endInstance_getAdjacency(atomInstance_get3End(atomInstance)));
					if(atomInstance2 == NULL) {
						break;
					}
					atom2 = j < atoms->length ? atoms->list[j] : NULL;
					if(atomInstance_getAtom(atomInstance2) != atom2) {
						break;
					}
					atomInstance = atomInstance2;
				}
			}
		}
		atom_destructInstanceIterator(instanceIterator);
	}

	//Now do actual construction;
	chainAlignment = malloc(sizeof(ChainAlignment));
	chainAlignment->matrix = malloc(sizeof(AtomInstance **) * atoms->length);
	for(i=0; i<atoms->length; i++) {
		chainAlignment->matrix[i] = malloc(sizeof(AtomInstance *) * list->length);
	}
	chainAlignment->columnNumber = atoms->length;
	chainAlignment->rowNumber = list->length;
	for(i=0; i<list->length; i++) {
		list2 = list->list[i];
		for(j=0; j<list2->length; j++) {
			chainAlignment->matrix[j][i] = list2->list[j];
		}
	}
	//Calculate the total length.
	chainAlignment->totalAlignmentLength = 0;
	for(i=0;atoms->length; i++) {
		chainAlignment->totalAlignmentLength += atom_getLength(atoms->list[i]);
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
	//free(chainAlignment->leftEnds);
	//free(chainAlignment->rightEnds);
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
	 * For number 2, building the trees for the atoms it:
	 *
	 * (a) Gets the sequences for the chains.. and makes them into mafs.
	 *
	 * (b) Writes out the
	 *
	 * (c) Runs the program to build the trees from this input.
	 *
	 * (d) Reads the results.
	 */
	NetDisk *netDisk;
	Net *net;
	Net *net2;
	int32_t startTime;
	int32_t i;
	Chain *chain;
	Link *link;
	End *end;
	End *end2;
	EndInstance *endInstance;
	Atom *atom;
	AdjacencyComponent *adjacencyComponent;
	struct List *atoms;
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

	initialiseTempFileTree(tempFileRootDirectory, 100, 4);

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
		atoms = constructEmptyList(0, NULL);
		for(i=0; i<chain_getLength(chain); i++) {
			link = chain_getLink(chain, i);
			end = link_getLeft(link);
			atom = end_getAtom(end);
			if(atom != NULL) {
				listAppend(atoms, atom);
			}
		}
		link = chain_getLink(chain, chain_getLength(chain)-1);
		end = link_getRight(link);
		atom = end_getAtom(end);
		if(atom != NULL) {
			listAppend(atoms, atom);
		}
		avl_insert(sortedChainAlignments, chainAlignment_construct(atoms));
		destructList(atoms);
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
			atoms = constructEmptyList(0, NULL);
			listAppend(atoms, atom);
			avl_insert(sortedChainAlignments, chainAlignment_construct(atoms));
			destructList(atoms);
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
	while((chainAlignment = avl_t_next(chainAlignmentIterator)) != NULL) {
		listAppend(list, chainAlignment);
	}
	buildChainTrees((ChainAlignment **)list->list, list->length, net_getEventTree(net), tempFileRootDirectory);
	destructList(list);

	logInfo("Augmented the atom trees in: %i seconds\n", time(NULL) - startTime);

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

	logInfo("Filled in end trees and augmented the event trees for the child nets in: %i seconds\n", time(NULL) - startTime);

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
