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
	Net_GroupIterator *iterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(iterator)) != NULL) {
		Net *nestedNet = group_getNestedNet(group);
		if(nestedNet != NULL) {
			checkNetsRecursively(nestedNet, net);
		}
	}
	net_destructGroupIterator(iterator);
}

void checkNonTerminalGroups(Net *net, Net *parentNet) {
	/*
	 * Check child parent connection between nets in the reconstruction tree.
	 */
	Group *group = net_getParentGroup(net);
	//Check the group is not empty.
	assert(group_getEndNumber(group) != 0);

	//Check the connection first
	assert(group_getNet(group) == parentNet);
	assert(net == group_getNestedNet(group));
	assert(parentNet == group_getNet(group));
	assert(net_getName(net) == group_getName(group));
	assert(net_getGroup(parentNet, group_getName(group)) == group);

	//Check the ends.
	Group_EndIterator *endIterator = group_getEndIterator(group);
	End *end;
	while((end = group_getNextEnd(endIterator)) != NULL) {
		assert(net_getEnd(parentNet, end_getName(end)) == end);
		assert(end_getNet(end) == parentNet);
		assert(net_getEnd(net, end_getName(end)) != NULL);
		assert(!end_isBlockEnd(net_getEnd(net, end_getName(end))));
		assert(end_getGroup(end) == group);
	}
	group_destructEndIterator(endIterator);

	//Check if part of link
	Link *link = group_getLink(group);
	if(link != NULL) {
		assert(link_getGroup(link) == group);
		assert(group_getEndNumber(group) >= 2); //other checks are in check chains.
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
			Group *group = link_getGroup(link);
			assert(group_getEndNumber(group) >= 2);
			assert(group_getEnd(group, end_getName(link_getLeft(link))) == link_getLeft(link));
			assert(group_getEnd(group, end_getName(link_getRight(link))) == link_getRight(link));
			if(i > 1) {
				Link *link2 = chain_getLink(chain, i-1);
				assert(end_isBlockEnd(link_getRight(link2)));
				assert(end_isBlockEnd(link_getLeft(link)));
				assert(end_getBlock(link_getRight(link2)) == end_getBlock(link_getLeft(link)));
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
		 * Check end is part of an group
		 */
		Group *group = end_getGroup(end);
		assert(group != NULL);
		assert(group_getEnd(group, end_getName(end)) == end);

		/*
		 * Check stub/cap/block-end status
		 */
		if(end_isBlockEnd(end)) {
			assert(end_getBlock(end) != NULL);
			assert(!end_isCap(end));
			assert(!end_isStub(end));
		}
		else {
			assert(end_getBlock(end) == NULL);
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
					assert(endInstance_getEvent(endInstance) == eventTree_getRootEvent(net_getEventTree(net))); //checks root is aligned with the root.
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
			 * Check the required adjacency between the end instances.
			 */
			endInstance2 = endInstance_getAdjacency(endInstance);
			assert(endInstance2 != NULL);
			assert(end_getGroup(endInstance_getEnd(endInstance2)) == group); //check they have the same group.
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

void checkBlocks(Net *net) {
	Net_BlockIterator *blockIterator = net_getBlockIterator(net);
	Block *block;
	while((block = net_getNextBlock(blockIterator)) != NULL) {
		assert(block_getOrientation(block));
		assert(block_getLength(block) > 0);

		/*
		 * Check we have two ends.
		 */
		End *leftEnd = block_getLeftEnd(block);
		End *rightEnd = block_getRightEnd(block);
		assert(leftEnd != NULL);
		assert(rightEnd != NULL);
		assert(end_getBlock(leftEnd) == block);
		assert(end_getBlock(rightEnd) == block);

		/*
		 * Check the instances.
		 */
		Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
		BlockInstance *blockInstance;
		while((blockInstance = block_getNext(instanceIterator)) != NULL) {
			assert(blockInstance_getBlock(blockInstance) == block);
			assert(blockInstance_getOrientation(blockInstance));
			/*
			 * Check block instance has two instances.
			 */
			EndInstance *_5EndInstance = blockInstance_get5End(blockInstance);
			EndInstance *_3EndInstance = blockInstance_get3End(blockInstance);
			assert(_5EndInstance != NULL);
			assert(_3EndInstance != NULL);
			assert(endInstance_getBlockInstance(_5EndInstance) == blockInstance);
			assert(endInstance_getBlockInstance(_3EndInstance) == blockInstance);
			/*
			 * Check the coordinates.
			 */
			assert(blockInstance_getLength(blockInstance) == block_getLength(block));
			assert(blockInstance_getStrand(blockInstance) == endInstance_getStrand(_5EndInstance));
			assert(blockInstance_getStrand(blockInstance) == endInstance_getStrand(_3EndInstance));
			assert(endInstance_getSide(_5EndInstance));
			assert(!endInstance_getSide(_3EndInstance));
			if(blockInstance_getStart(blockInstance) != INT32_MAX) {
				assert(endInstance_getCoordinate(_5EndInstance) == blockInstance_getStart(blockInstance));
				assert(endInstance_getCoordinate(_3EndInstance) == blockInstance_getStart(blockInstance) + (blockInstance_getStrand(blockInstance) ? blockInstance_getLength(blockInstance) - 1 : -blockInstance_getLength(blockInstance) + 1));
			}
		}
	}
	net_destructBlockIterator(blockIterator);
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

void checkBasesAccountedFor(Net *net) {
	int64_t totalBases = net_getTotalBaseLength(net);
	int64_t blockBases = 0.0;
	int64_t childBases = 0.0;
	Net_BlockIterator *blockIterator = net_getBlockIterator(net);
	Block *block;
	while((block = net_getNextBlock(blockIterator)) != NULL) {
		blockBases += block_getLength(block) * block_getInstanceNumber(block);
	}
	net_destructBlockIterator(blockIterator);
	Net_GroupIterator *iterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(iterator)) != NULL) {
		int64_t size = (int64_t)group_getTotalBaseLength(group);
		if(group_getNestedNet(group) != NULL) {
			assert(!group_isTerminal(group));
			assert(net_getTotalBaseLength(group_getNestedNet(group)) == size);
		}
		else {
			assert(group_isTerminal(group));
		}
		assert(size >= 0);
		childBases += size;
	}
	net_destructGroupIterator(iterator);
	if(blockBases + childBases != totalBases) {
		fprintf(stderr, "Got %i block bases, %i childBases and %i total bases\n", (int)blockBases, (int)childBases, (int)totalBases);
	}
	assert(blockBases + childBases == totalBases);
}

void checkNetsRecursively(Net *net, Net *parentNet) {
	logInfo("Checking the net %s\n", netMisc_nameToStringStatic(net_getName(net)));
	if(parentNet != NULL) {
		checkNonTerminalGroups(net, parentNet);
	}
	checkChains(net);
	checkEnds(net);
	checkBlocks(net);
	checkEvents(net);
	checkBasesAccountedFor(net);
	callCheckNetsRecursively(net);
}

void usage() {
	fprintf(stderr, "cactus_tree, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --checkTrees : Check that each end and block has a tree structure which is reconcilable with the event tree\n");
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
