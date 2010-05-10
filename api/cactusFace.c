#include "cactusGlobalsPrivate.h"
#include "cactusFacePrivate.h"

/*
 * Private constructor
 */
Face * face_construct2(Name name, Net * net) {
	Face * face = calloc(1, sizeof(Face));
	face->net = net;
	face->name = name;
	net_addFace(net, face);
	return face;
}

/*
 * Basc constructor
 */
Face * face_construct(Net * net) {
	return face_construct2(netDisk_getUniqueID(net_getNetDisk(net)), net);
}

/*
 * Basic memory deallocator
 */
void face_destruct(Face * face)
{
	int32_t index;

	net_removeFace(face->net, face);

	if (face->cardinal > 0) {
		free(face->topNodes);
		for (index = 0; index < face->cardinal; index++)
			free(face->bottomNodes[index]);
		free(face->bottomNodes);
		free(face->bottomNodeNumbers);
		for (index = 0; index < face->cardinal; index++)
			free(face->derivedEdgeDestinations[index]);
		free(face->derivedEdgeDestinations);
		//Freeing face ends..
		for (index = 0; index < face_getCardinal(face); index++) {
			if(face->faceEnds[index] != NULL) {
				faceEnd_destruct(face->faceEnds[index]);
			}
		}
		free(face->faceEnds);
	}
	free(face);
}

Net *face_getNet(Face *face) {
	return face->net;
}

/*
 * Returns cardinal
 */
int32_t face_getCardinal(Face * face) {
	return face->cardinal;
}

/*
 * Get selected top node
 */
Cap * face_getTopNode(Face * face, int32_t index) {
	return face->topNodes[index];
}

Face_FaceEndIterator *face_getFaceEndIterator(Face *face) {
	Face_FaceEndIterator *iterator = mallocLocal(sizeof(Face_FaceEndIterator));
	iterator->face = face;
	iterator->index = 0;
	return iterator;
}

static FaceEnd *face_getFaceEnd(Face *face, int32_t index) {
	assert(index >= 0);
	assert(index < face_getCardinal(face));
	if(face->faceEnds[index] == NULL) {
		face->faceEnds[index] = faceEnd_construct(face, index);
	}
	return face->faceEnds[index];
}

FaceEnd *face_getNextFaceEnd(Face_FaceEndIterator *iterator) {
	if(iterator->index < face_getCardinal(iterator->face)) {
		return face_getFaceEnd(iterator->face, iterator->index++);
	}
	return NULL;
}

FaceEnd *face_getPrevious(Face_FaceEndIterator *iterator) {
	assert(iterator->index <= face_getCardinal(iterator->face));
	if(iterator->index > 0) {
		return face_getFaceEnd(iterator->face, --iterator->index);
	}
	return NULL;
}

Face_FaceEndIterator *face_copyFaceEndIterator(Face_FaceEndIterator *iterator) {
	Face_FaceEndIterator *iterator2 = mallocLocal(sizeof(Face_FaceEndIterator));
	iterator2->face = iterator->face;
	iterator2->index = iterator->index;
	return iterator2;
}

void face_destructFaceEndIterator(Face_FaceEndIterator *iterator) {
	free(iterator);
}

/*
 * Get selected derived destinations for selected top node
 */
Cap * face_getDerivedDestinationAtIndex(Face * face, int32_t topIndex, int32_t derivedIndex) {
#ifdef BEN_DEBUG
	assert(topIndex < face_getCardinal(face));
	assert(derivedIndex < face_getBottomNodeNumber(face, topIndex));
#endif
	return face->derivedEdgeDestinations[topIndex][derivedIndex];
}

/*
 * Get non-null derived destination of selected top node (useful for simple face)
 */
Cap * face_getDerivedDestination(Face * face, int32_t index) {
	int32_t bottomIndex;
#if BEN_DEBUG
	assert(index < face_getCardinal(face));
#endif

	for (bottomIndex = 0; bottomIndex < face_getBottomNodeNumber(face, index); bottomIndex++)
		if (face->derivedEdgeDestinations[index][bottomIndex])
			return face->derivedEdgeDestinations[index][bottomIndex];

	return NULL;
}

/*
 * Get the number of bottom nodes for the selected top node
 */
int32_t face_getBottomNodeNumber(Face * face, int32_t topIndex) {
	return face->bottomNodeNumbers[topIndex];
}

/*
 * Get selected bottom node from selected top node in face
 */
Cap * face_getBottomNode(Face * face, int32_t topNodeIndex, int32_t bottomNodeIndex) {
	return face->bottomNodes[topNodeIndex][bottomNodeIndex];
}

/*
 * Allocate arrays to allow for data
 */
void face_allocateSpace(Face * face, int32_t cardinal) {
	face->cardinal = cardinal;
	face->topNodes = calloc(cardinal, sizeof(Cap *));
	face->bottomNodes =
	    calloc(cardinal, sizeof(Cap **));
	face->bottomNodeNumbers =
	    calloc(cardinal, sizeof(int32_t));
	face->derivedEdgeDestinations =
	    calloc(cardinal, sizeof(Cap **));
	face->faceEnds =
		calloc(cardinal, sizeof(FaceEnd *));
}

/*
 * Sets the selected top node
 */
void face_setTopNode(Face * face, int32_t topIndex, Cap * topNode) {
	if (topNode)
		topNode = cap_getPositiveOrientation(topNode);
	face->topNodes[topIndex] = topNode;
	cap_setFace(topNode, face);
}

/*
 * Set bottom node count and allocate space to store pointers
 */
void face_setBottomNodeNumber(Face * face, int32_t topIndex, int32_t number) {
	face->bottomNodeNumbers[topIndex] = 0;
	if(number) {
		face->bottomNodes[topIndex] = calloc(number, sizeof(Cap *));
		face->derivedEdgeDestinations[topIndex] = calloc(number, sizeof(Cap *));
	} else {
		face->bottomNodes[topIndex] = NULL;
		face->derivedEdgeDestinations[topIndex] = NULL;
	}
}

/*
 * Sets the derived edge destination for a given top node in face
 */
void face_setDerivedDestination(Face * face, int32_t topIndex, int32_t bottomIndex, Cap * destination) {
	if (destination)
		destination = cap_getPositiveOrientation(destination);
	face->derivedEdgeDestinations[topIndex][bottomIndex] = destination;
}

/*
 * Adds bottom node to selected top node in face
 */
void face_addBottomNode(Face * face, int32_t topIndex, Cap * bottomNode) {
	if (bottomNode)
		bottomNode = cap_getPositiveOrientation(bottomNode);
	face->bottomNodes[topIndex][face->bottomNodeNumbers[topIndex]++] = bottomNode;
}

/*
 * Adds a top node with a simgle bottom node
 */
void face_engineerArtificialNodes(Face * face, Cap * topNode, Cap * bottomNode, int32_t nonDerived) {
	int32_t index = face->cardinal++;

	if (topNode)
		topNode = cap_getPositiveOrientation(topNode);

	face->topNodes =
	    realloc(face->topNodes,
		    face_getCardinal(face) * sizeof(Cap *));
	face_setTopNode(face, index, topNode);

	face->bottomNodeNumbers =
	    realloc(face->topNodes, face_getCardinal(face) * sizeof(int32_t));
	face->bottomNodeNumbers[index] = 1;

	face->bottomNodes =
	    realloc(face->topNodes,
		    face_getCardinal(face) * sizeof(Cap **));
	face->bottomNodes[index] = malloc(sizeof(Cap *));
	face_addBottomNode(face, index, bottomNode);

	face->derivedEdgeDestinations =
	    realloc(face->topNodes,
		    face_getCardinal(face) * sizeof(Cap *));
	face->derivedEdgeDestinations[index] = malloc(sizeof(Cap *));
	face_setDerivedDestination(face, index, 0,
	    cap_getParent(face_getTopNode(face, nonDerived)));
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private face functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void face_setNet(Face *face, Net *net) {
	net_removeFace(face_getNet(face), face);
	face->net = net;
	net_addFace(net, face);
}

/*
 * Serialization write function for a given top node in a face
 */
static void face_writeBinaryRepresentationAtIndex(Face * face, int32_t index,
				    void (*writeFn) (const void *ptr,
						     size_t size,
						     size_t count))
{
	int32_t index2;

	// Top node's names
	binaryRepresentation_writeName(cap_getName(face->topNodes[index]), writeFn);

	// Number of bottom nodes
	binaryRepresentation_writeInteger(face->bottomNodeNumbers[index], writeFn);

	for (index2 = 0; index2 < face->bottomNodeNumbers[index]; index2++) {
		// Names of bottom nodes
		binaryRepresentation_writeName(cap_getName(face->bottomNodes[index][index2]), writeFn);

		// destination of the derived edge
		if (face->derivedEdgeDestinations[index][index2])
			binaryRepresentation_writeName(cap_getName(face->derivedEdgeDestinations[index][index2]), writeFn);
		else
			binaryRepresentation_writeName(NULL_NAME, writeFn);
	}
}

/*
 * Serialisation write function
 */
void face_writeBinaryRepresentation(Face * face,
				    void (*writeFn) (const void *ptr,
						     size_t size,
						     size_t count))
{
	int32_t index;

	binaryRepresentation_writeElementType(CODE_FACE, writeFn);
	binaryRepresentation_writeName(face->name, writeFn);
	binaryRepresentation_writeInteger(face->cardinal, writeFn);

	for (index = 0; index < face->cardinal; index++)
		face_writeBinaryRepresentationAtIndex(face, index, writeFn);
}

/*
 * Serialisation read function for a given top node in a face
 */
static void face_loadFromBinaryRepresentationAtIndex(void **binaryString, Face * face, int32_t index, Net * net) {
	int32_t num, index2;
	Cap * cap;
	Name name;

	// Top nodes
	name = binaryRepresentation_getName(binaryString);
	cap = net_getCap(net, name);
	face->topNodes[index] = cap;
	cap_setFace(cap, face);

	// Number of bottom nodes
	num = binaryRepresentation_getInteger(binaryString);
	face->bottomNodeNumbers[index] = num;
	face->bottomNodes[index] = calloc(num, sizeof(Cap*));
	face->derivedEdgeDestinations[index] = calloc(num, sizeof(Cap*));

	// Names of bottom nodes
	for(index2 = 0; index2 < num; index2++) {
		name = binaryRepresentation_getName(binaryString);
		cap = net_getCap(net, name);
		face->bottomNodes[index][index2] = cap;

		// Derived edge destination
		name = binaryRepresentation_getName(binaryString); //name may be null_name
		cap = net_getCap(net, name);
	#ifdef BEN_DEBUG
		if(name == NULL_NAME)
			assert(cap == NULL);
	#endif
		face->derivedEdgeDestinations[index][index2] = cap;
	}
}

/*
 * Serialisation read function
 */
Face *face_loadFromBinaryRepresentation(void **binaryString, Net * net)
{
	Face *face = NULL;
	int32_t num, index;

	if (binaryRepresentation_peekNextElementType(*binaryString) ==
	    CODE_FACE) {
		binaryRepresentation_popNextElementType(binaryString);
		face = calloc(1, sizeof(Face));
		face->net = net;
		face->name = binaryRepresentation_getName(binaryString);
		num = binaryRepresentation_getInteger(binaryString);

		face->cardinal = num;
		face->topNodes = calloc(num, sizeof(Cap*));
		face->bottomNodes = calloc(num, sizeof(Cap**));
		face->bottomNodeNumbers = calloc(num, sizeof(int32_t));
		face->derivedEdgeDestinations = calloc(num, sizeof(Cap**));

		for (index = 0; index < num; index++)
			face_loadFromBinaryRepresentationAtIndex(binaryString, face, index, net);
	}

	return face;
}

FaceEnd *face_getFaceEndForTopNode(Face *face, Cap *cap) {
	int32_t i;
	for(i=0; i<face_getCardinal(face); i++) {
		if(face_getTopNode(face, i) == cap) {
			return face_getFaceEnd(face, i);
		}
	}
	assert(0);
	return NULL;
}

/*
 * Returns name of face (i.e. name of first cap)
 */
Name face_getName(Face * face) {
	return face->name;
}

/*
 * Function needed for unit testing
 */
Face *face_getStaticNameWrapper(Name name) {
	static Face face;
	face.name = name;
	return &face;
}

/*
 * Tests if all the top nodes are separate from the bottom nodes
 */
static int32_t face_isStronglyAcyclic(Face * face) {
	int32_t cardinal = face_getCardinal(face);
	int32_t topIndex1, topIndex, bottomIndex;
	Event * topEvent, * bottomEvent;

	for (topIndex1 = 0; topIndex1 < cardinal; topIndex1++) {
		topEvent = cap_getEvent(face_getTopNode(face, topIndex1));

		for (topIndex = 0; topIndex < cardinal; topIndex++) {
			for (bottomIndex = 0; bottomIndex < face_getBottomNodeNumber(face, topIndex); bottomIndex++) {
				bottomEvent = cap_getEvent(face_getBottomNode(face, topIndex, bottomIndex));
				if (bottomEvent == topEvent
				    || event_isAncestor(bottomEvent, topEvent))
					return false;
			}
		}
	}

	return true;
}

/*
 * Tests if all the descent paths coming out of one top node in face are edge disjoint
 */
static int32_t face_hasMergedDescentEdgesAtIndex(Face * face, int32_t topIndex) {
	int32_t bottomNodeNumber = face_getBottomNodeNumber(face, topIndex);
	int32_t index;
	Cap * topNode = face_getTopNode(face, topIndex);
	Cap * current;
	struct List * list = constructZeroLengthList(100, NULL);

	for (index = 0; index < bottomNodeNumber; index++) {
		current = face_getBottomNode(face, topIndex, index);
		while(current != topNode) {
			if (listContains(list, current)) {
				destructList(list);
				return true;
			} else
				listAppend(list, current);
			current = cap_getParent(current);
		}
	}

	destructList(list);
	return false;
}

/*
 * Tests is descent paths are edge-disjoint in face
 */
static int32_t face_hasSeparateDescentEdges(Face * face) {
	int32_t cardinal = face_getCardinal(face);
	int32_t index;

	for (index = 0; index < cardinal; index++)
		if (face_hasMergedDescentEdgesAtIndex(face, index))
			return false;

	return true;
}

/*
 * Returns the first attched ancestor of cap
 */
static Cap *face_getAttachedAncestor(Cap * cap)
{
	Cap *current = cap_getParent(cap);
	Cap *parent;

	while (current
	       && !cap_getAdjacency(cap)
	       && (parent = cap_getParent(current)))
		current = parent;

	return current;
}

/*
 * Tests if a given top node is part of a simple alternating cycle
 */
static int32_t face_breaksSimpleAlternatingPath(Face * face, int topIndex) {
	int32_t bottomNodeNumber = face_getBottomNodeNumber(face, topIndex);
	int32_t index;
	Cap * topNode = face_getTopNode(face, topIndex);
	Cap * partner = cap_getAdjacency(topNode);
	Cap * bottomNode;
	Cap * bottomNodePartner;
	Cap * bottomNodePartnerAncestor;
	Cap * derivedEdgeDestination = NULL;

	for (index = 0; index < bottomNodeNumber; index++) {
		bottomNode = face_getBottomNode(face, topIndex, index);
		bottomNodePartner = cap_getAdjacency(bottomNode);
#ifdef BEN_DEBUG
		if (!bottomNodePartner)
			abort();
#endif
		bottomNodePartnerAncestor = face_getAttachedAncestor(bottomNodePartner);
		if (bottomNodePartnerAncestor == partner)
			continue;
		else if (derivedEdgeDestination == NULL)
			derivedEdgeDestination = bottomNodePartnerAncestor;
		else
			return true;
	}

	return false;
}

/*
 * Tests if face is a simple alternating cycle
 */
static int32_t face_isSimpleAlternatingPath(Face * face) {
	int32_t cardinal = face_getCardinal(face);
	int32_t index;

	for(index = 0; index < cardinal; index++)
		if (face_breaksSimpleAlternatingPath(face, index))
			return false;

	return true;
}

/*
 * Tests if a face is simple
 */
int32_t face_isSimple(Face * face) {
	return face_isStronglyAcyclic(face)
	       && face_hasSeparateDescentEdges(face)
	       && face_isSimpleAlternatingPath(face);
}

/*
 * Tests if regular
 */
int32_t face_isRegular(Face * face)
{
	int32_t index;

	if (!face_isSimple(face))
		return false;

	for (index = 0; index < face_getCardinal(face); index++)
		if (face_getBottomNodeNumber(face, index) > 1)
			return false;

	return true;
}

/*
 * Tests if canonical
 */
int32_t face_isCanonical(Face * face)
{
	int32_t index;

	if (!face_isRegular(face))
		return false;

	for (index = 0; index < face_getCardinal(face); index++)
		if (face_getBottomNodeNumber(face, index) == 1)
			return false;

	return true;
}


/*
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * The following functions are all in aid of checking that the set of faces we have is well
 * formed.
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */

void face_check(Face *face) {
	assert(face_isSimple(face));
}

static Hash *hashbottomCaps(Net *net) {
	/*
	 * For each top node finds the corresponding set of bottom nodes and returns a
	 * hash of top nodes to sets of bottom nodes.
	 */
	Hash *bottomCapsHash = hash_construct();
	Event *rootEvent = eventTree_getRootEvent(net_getEventTree(net));
	Cap *cap;
	Net_CapIterator *capIterator = net_getCapIterator(net);

	while((cap = net_getNextCap(capIterator)) != NULL) {
		if(cap_getEvent(cap) != rootEvent && cap_getAdjacency(cap) != NULL) {
			Cap *cap2 = cap_getTopCap(cap);
			struct List *list;
			if((list = hash_search(bottomCapsHash, cap2)) == NULL) {
				list = constructEmptyList(0, NULL);
				hash_insert(bottomCapsHash, cap2, list);
			}
			listAppend(list, cap);
		}
	}
	net_destructCapIterator(capIterator);

	return bottomCapsHash;
}

static Hash *computeLiftedEdges(Net *net, Hash *bottomCapsHash) {
	/*
	 * For each top node finds the set of top nodes connected to it by a lifted
	 * edge. Returns a hash of top nodes to list of other top nodes
	 * connected by lifted edges.
	 */
	Hash *liftedEdgesHash = hash_construct();
	Hash_Iterator *iterator = hash_getIterator(bottomCapsHash);
	Cap *topCap;
	while((topCap = hash_getNext(iterator)) != NULL) {
		struct List *bottomCaps = hash_search(bottomCapsHash, topCap);
		struct List *liftedEdges = constructEmptyList(0, NULL);
		assert(hash_search(liftedEdgesHash, topCap) == NULL);
		hash_insert(liftedEdgesHash, topCap, liftedEdges);
		int32_t i;
		for(i=0; i<bottomCaps->length; i++) {
			Cap *bottomCap = bottomCaps->list[i];
			Cap *adjacentBottomCap = cap_getAdjacency(bottomCap);
			assert(adjacentBottomCap != NULL);
			Cap *adjacentTopCap = cap_getTopCap(adjacentBottomCap);
			assert(adjacentTopCap != NULL);
			assert(hash_search(bottomCapsHash, adjacentTopCap) != NULL);
			listAppend(liftedEdges, adjacentTopCap);
		}
	}
	hash_destructIterator(iterator);
	return liftedEdgesHash;
}

static void computeModulesP(Cap *topCap, Hash *liftedEdgesHash, struct List *module,
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

static struct List *computeModules(Net *net, Hash *liftedEdges) {
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

static int32_t isTrivialFaceP(Cap *cap, Hash *bottomCapsHash, Cap *otherCap) {
	struct List *bottomCaps = hash_search(bottomCapsHash, cap);
	assert(bottomCaps != NULL);
	int32_t i;
	//Check is connected to only the other cap by the lifted edge.
	for(i=0; i<bottomCaps->length; i++) {
		if(bottomCaps->list[i] != otherCap) {
			return 0;
		}
	}
	//check paths are disjoint.
	assert(bottomCaps->length > 0);
	Event *event1 = cap_getEvent(bottomCaps->list[0]);
	for(i=1; i<bottomCaps->length; i++) {
		Event *event2 = cap_getEvent(bottomCaps->list[i]);
		if(eventTree_getCommonAncestor(event1, event2) != cap_getEvent(cap)) {
			return 0;
		}
	}
	//if there is an adjacency, check it is to the other cap
	return cap_getAdjacency(cap) == NULL ? 1 : cap_getAdjacency(cap) == otherCap;
}

static int32_t isTrivialFace(Net *net, struct List *module, Hash *bottomCapsHash) {
	/*
	 * Returns non-zero iff the module is trivial.
	 */
	if(module->length > 2) {
		return 0;
	}
	assert(module->length == 2);

	return isTrivialFaceP(module->list[0], bottomCapsHash, module->list[1]) &&
			isTrivialFaceP(module->list[1], bottomCapsHash, module->list[0]);
}

static void checkTrivialFaceP(Net *net, Cap *topCap, Hash *bottomCapsHash) {
	assert(cap_getTopFace(topCap) == NULL); //top node of trivial face can not be in a cap.
	assert(cap_getTopFace(cap_getReverse(topCap)) == NULL); //defensive, again top node of trivial face can not be in a cap.
	struct List *bottomCaps = hash_search(bottomCapsHash, topCap);
	assert(bottomCaps != NULL);
	int32_t i;
	for(i=0; i<bottomCaps->length; i++) { //check nodes do not have a bottom cap..
		Cap *bottomCap = bottomCaps->list[i];
		assert(cap_getBottomFaceEnd(bottomCap) == NULL);
	}
	assert(net != NULL);
}

static void checkTrivialFace(Net *net, struct List *module, Hash *bottomCapsHash) {
	/*
	 * Checks a trivial face.
	 */
	//Checks the top nodes do not have an associated fae.
	assert(module->length == 2);
	checkTrivialFaceP(net, module->list[0], bottomCapsHash);
	checkTrivialFaceP(net, module->list[1], bottomCapsHash);
}

static void checkNonTrivialFace(Net *net, struct List *module, Hash *bottomCapsHash) {
	/*
	 * Checks a non-trivial face.
	 */
	//Checks the top nodes are all in one associated face.
	//Checks the set of bottom nodes for each face are in agreement.
	int32_t i, k;
	assert(module->length > 0);
	Cap *topCap = module->list[0];
	Face *face = cap_getTopFace(topCap);
	assert(face != NULL);
	assert(face_getCardinal(face) == module->length);

	for(i=0; i<module->length; i++) {
		topCap = module->list[i];
		FaceEnd *faceEnd = cap_getTopFaceEnd(topCap);
		assert(faceEnd != NULL);
		assert(face == faceEnd_getFace(faceEnd));
		assert(faceEnd_getTopNode(faceEnd) == topCap);
		struct List *bottomCaps = hash_search(bottomCapsHash, topCap);
		assert(bottomCaps != NULL);
		assert(faceEnd_getNumberOfBottomNodes(faceEnd) == bottomCaps->length);
		for(k=0; k<bottomCaps->length; k++) {
			Cap *bottomCap = bottomCaps->list[k];
			assert(cap_getBottomFaceEnd(bottomCap) == faceEnd);
		}
	}
}

static void diffFaces(Net *net, struct List *modules, Hash *bottomCapsHash) {
	/*
	 * Compares the set of faces we compute with those stored.
	 */
	int32_t i, j = 0;
	for(i=0; i<modules->length; i++) {
		struct List *module = modules->list[i];
		if(isTrivialFace(net, module, bottomCapsHash)) {
			checkTrivialFace(net, module, bottomCapsHash);
		}
		else {
			checkNonTrivialFace(net, module, bottomCapsHash);
			j++;
		}
	}
	assert(j == net_getFaceNumber(net)); //we should have checked exactly the number of
	//non-trivial faces.
}

void face_checkFaces(Net *net) {
	/*
	 * Checks that the set of faces is as we expect - with a face created
	 * for each non-trivial face.
	 */
	if(net_builtFaces(net)) { //only check the faces if they have been built..
		Hash *bottomCapsHash = hashbottomCaps(net);

		//Construct lifted edges
		Hash *liftedEdgesHash = computeLiftedEdges(net, bottomCapsHash);

		//Constructs lifted edge/adjacency edge connected components, called modules.
		//Faces are simply the nodes in the modules (the top nodes) and the set of
		//bottom nodes.
		struct List *modules = computeModules(net, liftedEdgesHash);

		//Check all faces we have computed are the same as those computed by Daniel.
		diffFaces(net, modules, bottomCapsHash);

		//Cleanup
		hash_destruct(bottomCapsHash);
		hash_destruct(liftedEdgesHash);
		destructList(modules);
	}
}

