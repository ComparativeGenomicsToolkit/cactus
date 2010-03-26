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


//Check for internal instances and trees, call the flag checkTrees
//Check internal instances adjacencies, call the flag checkInternalAdjacencies
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
		assert(end_isStubEnd(net_getEnd(net, end_getName(end))));
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
			//Check stub ends are not free stubs.
			if(end_isStubEnd(link_getLeft(link))) {
				assert(end_isAttached(link_getLeft(link)));
			}
			if(end_isStubEnd(link_getRight(link))) {
				assert(end_isAttached(link_getRight(link)));
			}
			//Check internal links link blocks.
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
	int32_t attachedEnds = 0;
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
			assert(!end_isStubEnd(end));
			//check not attached
			assert(end_isFree(end));
			assert(!end_isAttached(end));
		}
		else {
			//is stub
			assert(end_isStubEnd(end));
			assert(end_getBlock(end) == NULL);

			//check attachment
			if(end_isAttached(end)) {
				attachedEnds++;
				assert(!end_isFree(end));
			}
			else {
				assert(end_isFree(end));
			}
		}

		End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
		Cap *cap, *cap2;
		while((cap = end_getNext(instanceIterator)) != NULL) {
			/*
			 * Check the cap is connected.
			 */
			assert(cap_getEnd(cap) == end);

			/*
			 * Check end trees.
			 */
			if(CHECK_TREES) {
				//check is labelled with an event.
				assert(cap_getEvent(cap) != NULL);

				cap2 = cap_getParent(cap);
				if(cap2 == NULL) { //must be the root.
					assert(end_getRootInstance(end) != NULL);
					assert(cap == end_getRootInstance(end));
					assert(cap_getEvent(cap) == eventTree_getRootEvent(net_getEventTree(net))); //checks root is aligned with the root.
				}
				else {
					//check parent and child are properly connected.
					j = 0;
					for(i=0; i<cap_getChildNumber(cap2); i++) {
						if(cap_getChild(cap2, i) == cap) {
							assert(!j);
							j = 1;
						}
					}
					assert(j);

					//check parent event is ancestral.
					assert(event_isAncestor(cap_getEvent(cap), cap_getEvent(cap2)));
				}
			}

			/*
			 * Skip checking internal instances if needed.
			 */
			if(cap_isInternal(cap) && !CHECK_INTERNAL_ADJACENCIES) {
				continue;
			}

			/*
			 * Check the required adjacency between the caps.
			 */
			cap2 = cap_getAdjacency(cap);
			assert(cap2 != NULL);
			assert(end_getGroup(cap_getEnd(cap2)) == group); //check they have the same group.
			assert(cap_getAdjacency(cap2) == cap);
			assert(cap_getEvent(cap) == cap_getEvent(cap2));
			assert(cap_getSequence(cap) == cap_getSequence(cap2));
			assert(cap_getStrand(cap) == cap_getStrand(cap2));
			if(cap_getCoordinate(cap) != INT32_MAX) { //if they have a coordinate
				assert(cap_getSide(cap) != cap_getSide(cap2));
				if(cap_getStrand(cap)) {
					if(!cap_getSide(cap)) {
						assert(cap_getCoordinate(cap) < cap_getCoordinate(cap2));
					}
					else {
						assert(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
					}
				}
				else {
					if(cap_getSide(cap)) {
						assert(cap_getCoordinate(cap) < cap_getCoordinate(cap2));
					}
					else {
						assert(cap_getCoordinate(cap) > cap_getCoordinate(cap2));
					}
				}
			}
		}
	}
	net_destructEndIterator(endIterator);
	//Check we have an even number of attached ends.
	assert(attachedEnds % 2 == 0);
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
		Segment *segment;
		while((segment = block_getNext(instanceIterator)) != NULL) {
			assert(segment_getBlock(segment) == block);
			assert(segment_getOrientation(segment));
			/*
			 * Check segment has two instances.
			 */
			Cap *_5Cap = segment_get5Cap(segment);
			Cap *_3Cap = segment_get3Cap(segment);
			assert(_5Cap != NULL);
			assert(_3Cap != NULL);
			assert(cap_getSegment(_5Cap) == segment);
			assert(cap_getSegment(_3Cap) == segment);
			/*
			 * Check the coordinates.
			 */
			assert(segment_getLength(segment) == block_getLength(block));
			assert(segment_getStrand(segment) == cap_getStrand(_5Cap));
			assert(segment_getStrand(segment) == cap_getStrand(_3Cap));
			assert(cap_getSide(_5Cap));
			assert(!cap_getSide(_3Cap));
			if(segment_getStart(segment) != INT32_MAX) {
				assert(cap_getCoordinate(_5Cap) == segment_getStart(segment));
				assert(cap_getCoordinate(_3Cap) == segment_getStart(segment) + (segment_getStrand(segment) ? segment_getLength(segment) - 1 : -segment_getLength(segment) + 1));
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
	Block_InstanceIterator *segmentIterator;
	Segment *segment;
	while((block = net_getNextBlock(blockIterator)) != NULL) {
		segmentIterator = block_getInstanceIterator(block);
		while((segment = block_getNext(segmentIterator)) != NULL) {
			if(segment_getSequence(segment) != NULL) {
				blockBases += segment_getLength(segment);
			}
		}
		block_destructInstanceIterator(segmentIterator);
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
