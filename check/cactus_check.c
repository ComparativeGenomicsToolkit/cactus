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

Cap *getTopNode(Cap *cap) {
	/*
	 * Gets the corresponding top node of a bottom node.
	 */
	assert(cap_getAdjacency(cap) != NULL);
	assert(end_getRootInstance(cap_getEnd(cap)) != cap);
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

Hash *hashTopNodes(Net *net) {
	/*
	 * For each bottom node finds the corresponding top node and returns a hash
	 * of bottom nodes to top nodes.
	 */
	Hash *topNodeHash = hash_construct();
	Event *rootEvent = eventTree_getRootEvent(net_getEventTree(net));
	Cap *cap;
	Net_CapIterator *capIterator = net_getCapIterator(net);

	while((cap = net_getNextCap(capIterator)) != NULL) {
		if(cap_getEvent(cap) != rootEvent && cap_getAdjacency(cap) != NULL) {
			Cap *cap2 = getTopNode(cap);
			assert(hash_search(topNodeHash, cap) == NULL);
			hash_insert(topNodeHash, cap, cap2);
		}
	}
	net_destructCapIterator(capIterator);

	return topNodeHash;
}

Hash *hashBottomNodes(Net *net) {
	/*
	 * For each top node finds the corresponding set of bottom nodes and returns a
	 * hash of top nodes to sets of bottom nodes.
	 */
	Hash *bottomNodesHash = hash_construct();
	Event *rootEvent = eventTree_getRootEvent(net_getEventTree(net));
	Cap *cap;
	Net_CapIterator *capIterator = net_getCapIterator(net);

	while((cap = net_getNextCap(capIterator)) != NULL) {
		if(cap_getEvent(cap) != rootEvent && cap_getAdjacency(cap) != NULL) {
			Cap *cap2 = getTopNode(cap);
			struct List *list;
			if((list = hash_search(bottomNodesHash, cap2)) == NULL) {
				list = constructEmptyList(0, NULL);
				hash_insert(bottomNodesHash, cap2, list);
			}
			listAppend(list, cap);
		}
	}
	net_destructCapIterator(capIterator);

	return bottomNodesHash;
}

Hash *computeLiftedEdges(Net *net, Hash *topNodeHash, Hash *bottomNodesHash) {
	/*
	 * For each top node finds the set of top nodes connected to it by a lifted
	 * edge. Returns a hash of top nodes to list of other top nodes
	 * connected by lifted edges.
	 */
	Hash *liftedEdgesHash = hash_construct();
	Hash_Iterator *iterator = hash_getIterator(bottomNodesHash);
	Cap *topCap;
	while((topCap = hash_getNext(iterator)) != NULL) {
		struct List *bottomNodes = hash_search(bottomNodesHash, topCap);
		struct List *liftedEdges = constructEmptyList(0, NULL);
		assert(hash_search(liftedEdgesHash, topCap) == NULL);
		hash_insert(liftedEdgesHash, topCap, liftedEdges);
		int32_t i;
		for(i=0; i<bottomNodes->length; i++) {
			Cap *bottomCap = bottomNodes->list[i];
			Cap *adjacentBottomCap = cap_getAdjacency(bottomCap);
			assert(adjacentBottomCap != NULL);
			Cap *adjacentTopCap = hash_search(topNodeHash, adjacentBottomCap);
			assert(adjacentTopCap != NULL);
			assert(hash_search(bottomNodesHash, adjacentTopCap) != NULL);
			listAppend(liftedEdges, adjacentTopCap);
		}
	}
	hash_destructIterator(iterator);
	return liftedEdgesHash;
}

void computeModulesP(Cap *topCap, Hash *liftedEdgesHash, struct List *module,
		Hash *modulesHash) {
	int32_t i;
	Cap *adjacentTopCap;
	if(hash_search(modulesHash, topCap) == NULL) {
		//Add to module
		hash_insert(modulesHash, topCap, module);
		listAppend(module, topCap);

		//Traverse the lifted edges
		struct List *liftedEdges = hash_search(liftedEdgesHash, topCap);
		for(i=0; i<liftedEdges->length; i++) {
			adjacentTopCap = liftedEdges->list[i];
			computeModulesP(adjacentTopCap, liftedEdgesHash, module, modulesHash);
		}

		//Traverse the direct adjaceny
		adjacentTopCap = cap_getAdjacency(topCap);
		if(adjacentTopCap != NULL) {
			computeModulesP(adjacentTopCap, liftedEdgesHash, module, modulesHash);
		}
	}
}

struct List *computeModules(Net *net, Hash *liftedEdges) {
	/*
	 * Finds the set of adjacency/lifted edge components, called modules,
	 * and returns them in a list.
	 */
	struct List *modules = constructEmptyList(0, (void (*)(void *))destructList);
	Hash *modulesHash = hash_construct();

	Hash_Iterator *iterator = hash_getIterator(liftedEdges);
	Cap *topCap;
	while((topCap = hash_getNext(iterator)) != NULL) {
		if(hash_search(modulesHash, topCap) == NULL) {
			struct List *module = constructEmptyList(0, NULL);
			computeModulesP(topCap, liftedEdges, module, modulesHash);
			listAppend(modules, module);
			assert(module->length >= 2);
		}
	}
	hash_destructIterator(iterator);
	hash_destruct(modulesHash);
	return modules;
}

int32_t isTrivialFaceP(Cap *cap, Hash *bottomNodesHash, Cap *otherCap) {
	struct List *bottomNodes = hash_search(bottomNodesHash, cap);
	assert(bottomNodes != NULL);
	int32_t i;
	//Check is connected to only the other cap by the lifted edge.
	for(i=0; i<bottomNodes->length; i++) {
		if(bottomNodes->list[i] != otherCap) {
			return 0;
		}
	}
	//check paths are disjoint.
	assert(bottomNodes->length > 0);
	Event *event1 = cap_getEvent(bottomNodes->list[0]);
	for(i=1; i<bottomNodes->length; i++) {
		Event *event2 = cap_getEvent(bottomNodes->list[i]);
		if(eventTree_getCommonAncestor(event1, event2) != cap_getEvent(cap)) {
			return 0;
		}
	}
	//if there is an adjacency, check it is to the other cap
	return cap_getAdjacency(cap) == NULL ? 1 : cap_getAdjacency(cap) == otherCap;
}

int32_t isTrivialFace(Net *net, struct List *module, Hash *bottomNodesHash) {
	/*
	 * Returns non-zero iff the module is trivial.
	 */
	if(module->length > 2) {
		return 0;
	}
	assert(module->length == 2);

	return isTrivialFaceP(module->list[0], bottomNodesHash, module->list[1]) &&
			isTrivialFaceP(module->list[1], bottomNodesHash, module->list[0]);
}

void checkTrivialFaceP(Net *net, Cap *topCap, Hash *bottomNodesHash) {
	assert(cap_getFace(topCap) == NULL); //top node of trivial face can not be in a cap.
	assert(cap_getFace(cap_getReverse(topCap)) == NULL); //defensive, again top node of trivial face can not be in a cap.
	//Currently we don't look at the bottom nodes
	struct List *bottomNodes = hash_search(bottomNodesHash, topCap);
	assert(bottomNodes != NULL);
	assert(net != NULL);
}

void checkTrivialFace(Net *net, struct List *module, Hash *bottomNodesHash) {
	/*
	 * Checks a trivial face.
	 */
	//Checks the top nodes do not have an associated fae.
	assert(module->length == 2);
	checkTrivialFaceP(net, module->list[0], bottomNodesHash);
	checkTrivialFaceP(net, module->list[1], bottomNodesHash);
}

void checkNonTrivialFace(Net *net, struct List *module, Hash *bottomNodesHash) {
	/*
	 * Checks a non-trivial face.
	 */
	//Checks the top nodes are all in one associated face.
	//Checks the set of bottom nodes for each face are in agreement.
	int32_t i, j, k, l;
	assert(module->length > 0);
	Cap *topCap = module->list[0];
	Face *face = cap_getFace(topCap);
	assert(face != NULL);
	for(i=1; i<module->length; i++) {
		topCap = module->list[i];
		assert(face == cap_getFace(topCap));
		j = 0; //face_getTopNodeIndex(face, topCap);
		assert(face_getTopNode(face, j) == topCap);
		struct List *bottomNodes = hash_search(bottomNodesHash, topCap);
		assert(bottomNodes != NULL);
		for(k=0; k<bottomNodes->length; k++) {
			Cap *bottomCap = bottomNodes->list[k];
		    l = 0; //face_getBottomNodeIndex(face, j, bottomCap);
		    assert(face_getBottomNode(face, j, l) == bottomCap);
		}
	}
}

void diffFaces(Net *net, struct List *modules, Hash *bottomNodesHash) {
	/*
	 * Compares the set of faces we compute with those stored.
	 */
	int32_t i;
	for(i=0; i<modules->length; i++) {
		struct List *module = modules->list[i];
		if(isTrivialFace(net, module, bottomNodesHash)) {
			checkTrivialFace(net, module, bottomNodesHash);
		}
		else {
			checkNonTrivialFace(net, module, bottomNodesHash);
		}
	}
}

static void checkFaces(Net *net) {
	/*
	 * Checks that the set of faces is as we expect - with a face created
	 * for each non-trivial face.
	 */
	return 0;
	if(net_builtFaces(net)) { //only check the faces if they have been built..
		Hash *topNodeHash = hashTopNodes(net);
		Hash *bottomNodesHash = hashBottomNodes(net);

		//Construct lifted edges
		Hash *liftedEdgesHash = computeLiftedEdges(net, topNodeHash, bottomNodesHash);

		//Constructs lifted edge/adjacency edge connected components, called modules.
		//Faces are simply the nodes in the modules (the top nodes) and the set of
		//bottom nodes.
		struct List *modules = computeModules(net, liftedEdgesHash);

		//Check all faces we have computed are the same as those computed by Daniel.
		diffFaces(net, modules, bottomNodesHash);
	}
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
	checkFaces(net);
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
