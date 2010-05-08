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
/*
Hash *computeLiftedEdges(Net *net, Hash *hashA, Hash *hashInvA) {
	Hash *liftedEdgesHash = hash_construct();
	Hash_Iterator *iterator = hash_getIterator(hashInvA);
	Cap *topCap;
	while((topCap = hash_getNext(iterator)) != NULL) {
		struct List *bottomCaps = hashInvA(topCap);
		struct List *liftedEdges = constructEmptyList(0);
		hash_insert(liftedEdgesHash, topCap, liftedEdges);
		int32_t i;
		for(i=0; i<bottomNodes->length; i++) {
			Cap *bottomCap = bottomCaps->list[i];
			Cap *adjacentBottomCap = cap_getAdjacency(bottomCap);
			assert(adjacentBottomCap != NULL);
			Cap *adjacentTopCap = hash_search(adjacentBottomCap, hashA);
			assert(adjacentTopCap != NULL);
			assert(hash_search(hashInvA, adjacentTopCap) != NULL);
			listAppend(liftedEdgesHash, adjacentTopCap);
		}
	}
	hash_destructIterator(iterator);
}

void computeModulesP(Cap *topCap, Hash *liftedEdges, struct List *module,
		Hash *modulesHash) {
	int32_t i;
	Cap *adjacentTopCap;
	if(hash_search(modulesHash, topCap) == NULL) {
		//Add to module
		hash_insert(modulesHash, topCap, module);
		listAppend(module, topCap);

		//Traverse the lifted edges
		struct List *liftedEdges = hash_search(liftedEdges, topCap);
		for(i=0; i<liftedEdges->length; i++) {
			adjacentTopCap = liftedEdges->list[i];
			computeModulesP(adjacentTopCap, liftedEdges, module, modulesHash);
		}

		//Traverse the direct adjaceny
		adjacentTopCap = cap_getAdjacency(topCap);
		if(adjacentTopCap != NULL) {
			computeModulesP(adjacentTopCap, liftedEdges, module, modulesHash);
		}
	}
}

struct List *computeModules(Net *net, Hash *liftedEdges) {
	struct List *modules = constructEmptyList(0, destructList);
	Hash *modulesHash = hash_construct();

	Hash_Iterator *iterator = hash_getIterator(hashInvA);
	Cap *topCap;
	while((topCap = hash_getNext(iterator)) != NULL) {
		if(hash_search(modulesHash, topCap) == NULL) {
			struct List *module = constructEmptyList();
			computeModulesP(topCap, liftedEdges, module, modulesHash);
			listAppend(modules, module);
			assert(module->length >= 2);
		}
	}
	hash_destructIterator(iterator);
	hash_destruct(modulesHash);
	return modules;
}

int32_t isDisjoint

int32_t isTrivialFace(Net *net, struct List *module, Hash *hashInvA) {
	if(module->length > 2) {
		return 0;
	}
	assert(module->length == 2);
	return isSimpleFace(net, module, hashInvA);
}

void checkTrivialFaceP(Net *net, Cap *topCap, Hash *hashInvA) {
	assert(cap_getFace(topCap) == NULL); //top node of trivial face can not be in a cap.
	assert(cap_getFace(cap_getReverse(topCap)) == NULL); //defensive, again top node of trivial face can not be in a cap.
	//Currently we don't look at the bottom nodes
	struct List *bottomNodes = hash_search(hashInvA, topCap);
	assert(bottomNodes != NULL);
	assert(net != NULL);
}

void checkTrivialFace(Net *net, struct List *module, Hash *hashInvA) {
	//Checks the top nodes do not have an associated fae.
	assert(module->length == 2);
	checkTrivialFaceP(net, module->list[0], hashInvA);
	checkTrivialFaceP(net, module->list[1], hashInvA);
}

void checkNonTrivialFace(Net *net, struct List *module, Hash *hashInvA) {
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
		j = face_getTopNodeIndex(face, topCap);
		assert(face_getTopNode(face, j) == topCap);
		struct List *bottomCaps = hash_search(hashInvA, topCap);
		assert(bottomCaps != NULL);
		for(k=0; k<bottomCaps->length; k++) {
			Cap *bottomCap = bottomCaps->list[k];
		    l = face_getBottomNodeIndex(face, j, bottomCap);
		    assert(face_getBottomNode(j, l) == bottomCap);
		}
	}
}

void diffFaces(Net *net, struct List *modules, Hash *hashInvA) {
	int32_t i;
	for(i=0; i<modules; i++) {
		struct List *module = modules->list[i];
		if(isTrivialFace(module, hashInvA)) {
			checkTrivialFace(net, module, hashInvA);
		}
		else {
			checkNonTrivialFace(net, module, hashInvA);
		}
	}
}

static void checkFaces(Net *net) {
	//Compute hashes for A(c) and A'(c) (see legends for def).
	Hash *hashA = computeHashA(net);
	Hash *hashInvA = computeHashInvA(net);

	//Construct lifted edges using A(c) function hash.
	Hash *liftedEdgesHash = computeLiftedEdges(net, hashA);

	//Constructs lifted edge/adjacency edge connected components, called modules.
	//Faces are simply the nodes in the modules (the top nodes) and the set of
	//bottom nodes.
	struct List *modules = computeModules(net, liftedEdgesHash);

	//Check all faces we have computed are the same as those computed by Daniel.
	diffFaces(net, modules);
}*/

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
