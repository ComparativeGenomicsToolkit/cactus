#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "commonC.h"

/*
 * The script checks the nets are structured as we expect, essentially by
 * calling net_check for each net in the tree. We do a couple of other tests also..
 */

Cap *A(Cap *cap) {
	Cap *cap2 = cap_getParent(cap);
	while(1) {
		if(cap_getAdjacency(cap2) != NULL) {
			return cap2;
		}
		if(cap_getParent(cap2) == NULL) {
			assert(end_getRootInstance(cap_getEnd(cap2)) == cap2);
			return cap2;
		}
		cap2 = cap_getParent(cap2);
	}
}

Hash *computeHashA(Net *net) {
	/*
	 * For each attached, non-root cap c finds A(c) and puts it in a hash "hashA",
	 * such that hashA[c] = A(c).
	 * Where A(c) is the lowest common ancestor of c which is attached, else
	 * A(c) is the root ancestor of A(c).
	 */
	Hash *hashA = hash_construct();
	Event *rootEvent = eventTree_getRootEvent(net_getEventTree(net));
	Cap *cap;
	Net_CapIterator *capIterator = net_getCapIterator(net);

	while((cap = net_getNextCap(capIterator)) != NULL) {
		if(cap_getEvent(cap) != rootEvent && cap_getAdjacency(cap) != NULL) {
			Cap *cap2 = A(cap);
			assert(hash_search(hashA, cap) == NULL);
			hash_insert(hashA, cap, cap2);
		}
	}
	net_destructCapIterator(capIterator);

	return hashA;
}

Hash *computeHashInvA(Net *net) {
	/*
	 * For each attached or root cap c finds A'(c) and puts it in a hash "hashInvA",
	 *  such that hashInvA[c] = A'(c).
	 * Where A'(c) is the inverse of A(c) such that A'(c) = { d \in C | A(d) = c },
	 * and C is the set of all caps. The members of A'(c) are called bottom nodes.
	 *
	 */
	Hash *hashInvA = hash_construct();
	Event *rootEvent = eventTree_getRootEvent(net_getEventTree(net));
	Cap *cap;
	Net_CapIterator *capIterator = net_getCapIterator(net);

	while((cap = net_getNextCap(capIterator)) != NULL) {
		if(cap_getEvent(cap) != rootEvent && cap_getAdjacency(cap) != NULL) {
			Cap *cap2 = A(cap);
			struct List *list;
			if((list = hash_search(hashInvA, cap2)) == NULL) {
				list = constructEmptyList(0, NULL);
				hash_insert(hashInvA, cap2, list);
			}
			listAppend(list, cap);
		}
	}
	net_destructCapIterator(capIterator);

	return hashInvA;
}

static void checkFaces(Net *net) {
	/*
	 *
	 */
	//Compute hashes for A(c) and A'(c) (see legends for def).
	Hash *hashA = computeHashA(net);
	Hash *hashInvA = computeHashInvA(net);

	//Construct lifted edges using A(c) function hash.
	//Hash *liftedEdgesHash = computeLiftedEdges(net, hashA);

	//Constructs lifted edge/adjacency edge connected components, called modules.
	//Faces are simply the nodes in the modules (the top nodes) and the set of
	//bottom nodes.
	//struct List *modules = computeModules(net, liftedEdgesHash);

	//Check all trivial faces are not represented in the set of Faces stored in the net.


	//Check all non-trival faces are represented in the set of Faces stored in the net.
}

static void checkTreeIsTerminalNormalised(Net *net) {
	/*
	 * A cactus tree is terminally normalised if all leaf nets are terminal.
	 * We haven't made this part of the api checks yet, as I don't know if we'll keep it
	 * forever.
	 */
	if(net_isLeaf(net)) {
		assert(net_isTerminal(net));
		assert(net_getBlockNumber(net) == 0);
		//The following are defensive checks.
		Group *group;
		Net_GroupIterator *iterator = net_getGroupIterator(net);
		while((group = net_getNextGroup(iterator)) != NULL) {
			assert(group_isLeaf(group));
		}
		net_destructGroupIterator(iterator);
	}
}

static void checkBasesAccountedFor(Net *net) {
	/*
	 * Checks all the bases in a net end up in child net or a nested net.
	 */
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
			assert(!group_isLeaf(group));
			assert(net_getTotalBaseLength(group_getNestedNet(group)) == size);
		}
		else {
			assert(group_isLeaf(group));
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

static void checkNets(Net *net, int32_t recursive) {
	net_check(net);
	checkBasesAccountedFor(net);
	checkTreeIsTerminalNormalised(net);
	//Call problem recursively
	if(recursive) {
		Net_GroupIterator *iterator = net_getGroupIterator(net);
		Group *group;
		while((group = net_getNextGroup(iterator)) != NULL) {
			if(!group_isLeaf(group)) {
				checkNets(group_getNestedNet(group), 1);
			}
		}
		net_destructGroupIterator(iterator);
	}
}

void usage() {
	fprintf(stderr, "cactus_tree, version 0.2\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-e --recursive : Check all nets recursively\n");
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
	int32_t recursive = 0;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "recursive", no_argument, 0, 'e' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "a:c:eh", long_options, &option_index);

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
				recursive = 1;
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

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = netDisk_construct(netDiskName);
	logInfo("Set up the net disk\n");

	int32_t j;
	for (j = optind; j < argc; j++) {
		const char *netName = argv[j];
		logInfo("Processing the net named: %s", netName);

		///////////////////////////////////////////////////////////////////////////
		// Parse the basic reconstruction problem
		///////////////////////////////////////////////////////////////////////////

		net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
		logInfo("Parsed the top level net of the cactus tree to check\n");

		///////////////////////////////////////////////////////////////////////////
		// Recursive check the nets.
		///////////////////////////////////////////////////////////////////////////

		checkNets(net, recursive);
		logInfo("Checked the nets/\n");
	}

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	netDisk_destruct(netDisk);

	return 0;
}
