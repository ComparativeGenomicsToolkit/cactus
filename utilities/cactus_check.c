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
 * The script checks the nets are structured as we expect. Returns non zero (and logs, if turned on) if it finds an error.
 *
 * The checks are not comprehensive, just what I've bothered to check.
 */


//Check for internal instances and trees, call this flag checkTrees
//Check internal instances adjacencies, call this flag checkInternalAdjacencies
//Stage 1, both checkTrees and checkInternalAdjacencies are false.
//Stage 2, checkTrees is true, checkInternalAdjacencies is false.
//Stage 3, checkTrees is true, checkInternalAdjacencies is true.

int32_t CHECK_TREES = 0;
int32_t CHECK_INTERNAL_ADJACENCIES = 0;

void checkNetsRecursively(Net *net, Net *parentNet);

void callCheckNetsRecursively(Net *net) {
	//Call problem recursively
	Net_AdjacencyComponentIterator *iterator = net_getAdjacencyComponentIterator(net);
	AdjacencyComponent *adjacencyComponent;
	while((adjacencyComponent = net_getNextAdjacencyComponent(iterator)) != NULL) {
		checkNetsRecursively(adjacencyComponent_getNestedNet(adjacencyComponent), net);
	}
	net_destructAdjacencyComponentIterator(iterator);
}

void checkAdjacencyComponents(Net *net, Net *parentNet) {
	/*
	 * Check child parent connection between nets in the reconstruction tree.
	 */
	AdjacencyComponent *adjacencyComponent = net_getParentAdjacencyComponent(net);
	//Check the adjacency component is not empty.
	assert(adjacencyComponent_getEndNumber(adjacencyComponent) != 0);

	//Check the connection first
	assert(adjacencyComponent_getNet(adjacencyComponent) == parentNet);
	assert(net == adjacencyComponent_getNestedNet(adjacencyComponent));
	assert(parentNet == adjacencyComponent_getNet(adjacencyComponent));
	assert(net_getName(net) == adjacencyComponent_getNestedNetName(adjacencyComponent));
	assert(net_getAdjacencyComponent(parentNet, adjacencyComponent_getNestedNetName(adjacencyComponent)) == adjacencyComponent);

	//Check the ends.
	AdjacencyComponent_EndIterator *endIterator = adjacencyComponent_getEndIterator(adjacencyComponent);
	End *end;
	while((end = adjacencyComponent_getNextEnd(endIterator)) != NULL) {
		assert(net_getEnd(parentNet, end_getName(end)) == end);
		assert(end_getNet(end) == parentNet);
		assert(net_getEnd(net, end_getName(end)) != NULL);
		assert(!end_isAtomEnd(net_getEnd(net, end_getName(end))));
		assert(end_getAdjacencyComponent(end) == adjacencyComponent);
	}
	adjacencyComponent_destructEndIterator(endIterator);

	//Check if part of link
	Link *link = adjacencyComponent_getLink(adjacencyComponent);
	if(link != NULL) {
		assert(link_getAdjacencyComponent(link) == adjacencyComponent);
		assert(adjacencyComponent_getEndNumber(adjacencyComponent) >= 2); //other checks are in check chains.
	}
}

void checkChains(Net *net) {
	/*
	 * Checks the chains are structured as we expect.
	 */
	Net_ChainIterator *chainIterator = net_getChainIterator(net);
	Chain *chain;
	while((chain = net_getNextChain(chainIterator)) != NULL) {
		int32_t i;
		for(i=0; i<chain_getLength(chain); i++) {
			Link *link = chain_getLink(chain, i);
			assert(link != NULL);
			assert(link_getChain(link) == chain);
			assert(end_getOrientation(link_getLeft(link)));
			assert(end_getOrientation(link_getRight(link)));
			AdjacencyComponent *adjacencyComponent = link_getAdjacencyComponent(link);
			assert(adjacencyComponent_getEndNumber(adjacencyComponent) >= 2);
			assert(adjacencyComponent_getEnd(adjacencyComponent, end_getName(link_getLeft(link))) == link_getLeft(link));
			assert(adjacencyComponent_getEnd(adjacencyComponent, end_getName(link_getRight(link))) == link_getRight(link));
			if(i > 1) {
				Link *link2 = chain_getLink(chain, i-1);
				assert(end_isAtomEnd(link_getRight(link2)));
				assert(end_isAtomEnd(link_getLeft(link)));
				assert(end_getAtom(link_getRight(link2)) == end_getAtom(link_getLeft(link)));
			}
		}
	}
	net_destructChainIterator(chainIterator);
}

void checkEnds(Net *net) {
	/*
	 * Checks the ends, including their adjacencies.
	 */
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	int32_t i, j;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		/*
		 * Check end is part of an adjacency component
		 */
		AdjacencyComponent *adjacencyComponent = end_getAdjacencyComponent(end);
		if(adjacencyComponent == NULL) {
			assert(net_getAdjacencyComponentNumber(end_getNet(end)) == 0);
		}
		else {
			assert(adjacencyComponent_getEnd(adjacencyComponent, end_getName(end)) == end);
		}

		/*
		 * Check stub/cap/atom-end status
		 */
		if(end_isAtomEnd(end)) {
			assert(end_getAtom(end) != NULL);
			assert(!end_isCap(end));
			assert(!end_isStub(end));
		}
		else {
			assert(end_getAtom(end) == NULL);
			assert((end_isCap(end) && !end_isStub(end)) || (!end_isCap(end) && end_isStub(end)));
		}

		End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
		EndInstance *endInstance, *endInstance2;
		while((endInstance = end_getNext(instanceIterator)) != NULL) {
			/*
			 * Check the end instance is connected.
			 */
			assert(endInstance_getEnd(endInstance) == end);

			/*
			 * Check end trees.
			 */
			if(CHECK_TREES) {
				//check is labelled with an event.
				assert(endInstance_getEvent(endInstance) != NULL);

				endInstance2 = endInstance_getParent(endInstance);
				if(endInstance2 == NULL) { //must be the root.
					assert(end_getRootInstance(end) != NULL);
					assert(endInstance == end_getRootInstance(end));
					//assert(endInstance_getEvent(endInstance) == eventTree_getRootEvent(net_getEventTree(net))); //checks root is aligned with the root.
				}
				else {
					//check parent and child are properly connected.
					j = 0;
					for(i=0; i<endInstance_getChildNumber(endInstance2); i++) {
						if(endInstance_getChild(endInstance2, i) == endInstance) {
							assert(!j);
							j = 1;
						}
					}
					assert(j);

					//check parent event is ancestral.
					assert(event_isAncestor(endInstance_getEvent(endInstance), endInstance_getEvent(endInstance2)));
				}
			}

			/*
			 * Skip checking internal instances if needed.
			 */
			if(endInstance_isInternal(endInstance) && !CHECK_INTERNAL_ADJACENCIES) {
				continue;
			}

			/*
			 * Check the required adjacency (1) between the end instances.
			 */
			endInstance2 = endInstance_getAdjacency(endInstance);
			assert(endInstance2 != NULL);
			assert(endInstance_getAdjacency(endInstance2) == endInstance);
			assert(endInstance_getEvent(endInstance) == endInstance_getEvent(endInstance2));
			assert(endInstance_getSequence(endInstance) == endInstance_getSequence(endInstance2));
			assert(endInstance_getStrand(endInstance) == endInstance_getStrand(endInstance2));
			if(endInstance_getCoordinate(endInstance) != INT32_MAX) { //if they have a coordinate
				assert(endInstance_getSide(endInstance) != endInstance_getSide(endInstance2));
				if(endInstance_getStrand(endInstance)) {
					if(!endInstance_getSide(endInstance)) {
						assert(endInstance_getCoordinate(endInstance) < endInstance_getCoordinate(endInstance2));
					}
					else {
						assert(endInstance_getCoordinate(endInstance) > endInstance_getCoordinate(endInstance2));
					}
				}
				else {
					if(endInstance_getSide(endInstance)) {
						assert(endInstance_getCoordinate(endInstance) < endInstance_getCoordinate(endInstance2));
					}
					else {
						assert(endInstance_getCoordinate(endInstance) > endInstance_getCoordinate(endInstance2));
					}
				}
			}
		}
	}
	net_destructEndIterator(endIterator);
}

void checkAtoms(Net *net) {
	Net_AtomIterator *atomIterator = net_getAtomIterator(net);
	Atom *atom;
	while((atom = net_getNextAtom(atomIterator)) != NULL) {
		assert(atom_getOrientation(atom));
		assert(atom_getLength(atom) > 0);

		/*
		 * Check we have two ends.
		 */
		End *leftEnd = atom_getLeftEnd(atom);
		End *rightEnd = atom_getRightEnd(atom);
		assert(leftEnd != NULL);
		assert(rightEnd != NULL);
		assert(end_getAtom(leftEnd) == atom);
		assert(end_getAtom(rightEnd) == atom);

		/*
		 * Check the instances.
		 */
		Atom_InstanceIterator *instanceIterator = atom_getInstanceIterator(atom);
		AtomInstance *atomInstance;
		while((atomInstance = atom_getNext(instanceIterator)) != NULL) {
			assert(atomInstance_getAtom(atomInstance) == atom);
			assert(atomInstance_getOrientation(atomInstance));
			/*
			 * Check atom instance has two instances.
			 */
			EndInstance *_5EndInstance = atomInstance_get5End(atomInstance);
			EndInstance *_3EndInstance = atomInstance_get3End(atomInstance);
			assert(_5EndInstance != NULL);
			assert(_3EndInstance != NULL);
			assert(endInstance_getAtomInstance(_5EndInstance) == atomInstance);
			assert(endInstance_getAtomInstance(_3EndInstance) == atomInstance);
			/*
			 * Check the coordinates.
			 */
			assert(atomInstance_getLength(atomInstance) == atom_getLength(atom));
			assert(atomInstance_getStrand(atomInstance) == endInstance_getStrand(_5EndInstance));
			assert(atomInstance_getStrand(atomInstance) == endInstance_getStrand(_3EndInstance));
			assert(endInstance_getSide(_5EndInstance));
			assert(!endInstance_getSide(_3EndInstance));
			if(atomInstance_getStart(atomInstance) != INT32_MAX) {
				assert(endInstance_getCoordinate(_5EndInstance) == atomInstance_getStart(atomInstance));
				assert(endInstance_getCoordinate(_3EndInstance) == atomInstance_getStart(atomInstance) + (atomInstance_getStrand(atomInstance) ? atomInstance_getLength(atomInstance) - 1 : -atomInstance_getLength(atomInstance) + 1));
			}
		}
	}
	net_destructAtomIterator(atomIterator);
}


void checkEvents(Net *net) {
	/*
	 * Checks all the events in the event tree, checks we have one root, and that
	 * all other events properly have a parent and are linked.
	 */
	EventTree *eventTree = net_getEventTree(net);
	Event *event;
	EventTree_Iterator *iterator = eventTree_getIterator(eventTree);
	while((event = eventTree_getNext(iterator)) != NULL) {
		Event *event2 = event_getParent(event);
		if(event2 == NULL) { //is the root event.
			assert(eventTree_getRootEvent(eventTree) == event);
		}
		else {
			assert(event_isAncestor(event, event2));
			int32_t i, j;
			j = 0;
			for(i=0; i<event_getChildNumber(event2); i++) {
				if(event_getChild(event2, i) == event) {
					assert(!j);
					j = 1;
				}
			}
			assert(j);
		}
	}
	eventTree_destructIterator(iterator);
}

void checkNetsRecursively(Net *net, Net *parentNet) {
	logInfo("Checking the net %s\n", netMisc_nameToStringStatic(net_getName(net)));
	if(parentNet != NULL) {
		checkAdjacencyComponents(net, parentNet);
	}
	checkChains(net);
	checkEnds(net);
	checkAtoms(net);
	checkEvents(net);
	callCheckNetsRecursively(net);
}

void usage() {
	fprintf(stderr, "cactus_tree, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --checkTrees : Check that each end and atom has a tree structure which is reconcilable with the event tree\n");
	fprintf(stderr, "-f --checkInternalAdjacencies : Checks that internal instances have adjacencies and that the operation structures which describe changes are correct\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
	NetDisk *netDisk;
	Net *net;

	/*
	 * Arguments/options
	 */
	char * logLevelString = NULL;
	char * netDiskName = NULL;
	char * netName = NULL;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "netName", required_argument, 0, 'd' },
			{ "checkTrees", no_argument, 0, 'e' },
			{ "checkInternalAdjacencies", no_argument, 0, 'f' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:c:d:efh", long_options, &option_index);

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
				CHECK_TREES = !CHECK_TREES;
				break;
			case 'f':
				CHECK_INTERNAL_ADJACENCIES = !CHECK_INTERNAL_ADJACENCIES;
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

	assert(netDiskName != NULL);
	assert(netName != NULL);

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

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
	logInfo("Parsed the top level net of the cactus tree to check\n");

	///////////////////////////////////////////////////////////////////////////
	// Recursive check the nets.
	///////////////////////////////////////////////////////////////////////////

	checkNetsRecursively(net, NULL);
	logInfo("Checked the nets/\n");

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	netDisk_destruct(netDisk);

	return 0;
}
