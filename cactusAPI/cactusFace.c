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
		free(face->derivedEdgeDestinations);
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

/*
 * Get selected derived destination  of selected top node
 */
Cap * face_getDerivedDestination(Face * face, int32_t index) {
	return face->derivedEdgeDestinations[index];
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
	    calloc(cardinal, sizeof(Cap *));
}

/*
 * Sets the selected top node
 */
void face_setTopNode(Face * face, int32_t topIndex, Cap * topNode) {
	face->topNodes[topIndex] = topNode;
	cap_setFace(topNode, face);
}

/*
 * Set bottom node count and allocate space to store pointers
 */
void face_setBottomNodeNumber(Face * face, int32_t topIndex, int32_t number) {
	face->bottomNodeNumbers[topIndex] = 0;
	if(number)
		face->bottomNodes[topIndex] = calloc(number, sizeof(Cap *));
	else
		face->bottomNodes[topIndex] = NULL;
}

/*
 * Sets the derived edge destination for a given top node in face
 */
void face_setDerivedDestination(Face * face, int32_t topIndex, Cap * destination) {
	face->derivedEdgeDestinations[topIndex] = destination;
}

/*
 * Adds bottom node to selected top node in face
 */
void face_addBottomNode(Face * face, int32_t topIndex, Cap * bottomNode) {
	face->bottomNodes[topIndex][face->bottomNodeNumbers[topIndex]++] = bottomNode;
}

/*
 * Adds a top node with a simgle bottom node
 */
void face_engineerArtificialNodes(Face * face, Cap * topNode, Cap * bottomNode, int32_t nonDerived) {
	int32_t index = face->cardinal++;
	face->topNodes =
	    realloc(face->topNodes,
		    face_getCardinal(face) * sizeof(Cap *));
	face_setTopNode(face, index, topNode);
	face->bottomNodeNumbers =
	    realloc(face->topNodes, face_getCardinal(face) * sizeof(int32_t));
	face->bottomNodes =
	    realloc(face->topNodes,
		    face_getCardinal(face) * sizeof(Cap **));
	face->bottomNodes[index] = malloc(sizeof(Cap *));
	face_addBottomNode(face, index, bottomNode);
	face->derivedEdgeDestinations =
	    realloc(face->topNodes,
		    face_getCardinal(face) * sizeof(Cap *));
	face_setDerivedDestination(face, index,
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

	// Names of bottom nodes
	for (index2 = 0; index2 < face->bottomNodeNumbers[index]; index2++)
		binaryRepresentation_writeName(cap_getName(face->bottomNodes[index][index2]), writeFn);

	// destination of the derived edge
	if (face->derivedEdgeDestinations[index]) {
		binaryRepresentation_writeName(cap_getName(face->derivedEdgeDestinations[index]), writeFn);
	}
	else {
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
	face->bottomNodes[index] = calloc(num, sizeof(Cap()));

	// Names of bottom nodes
	for(index2 = 0; index2 < num; index2++) {
		name = binaryRepresentation_getName(binaryString);
		cap = net_getCap(net, name);
		face->bottomNodes[index][index2] = cap;
	}

	// Derived edge destination
	name = binaryRepresentation_getName(binaryString); //name may be null_name
	cap = net_getCap(net, name);
#ifdef BEN_DEBUG
	if(name == NULL_NAME) {
		assert(cap == NULL);
	}
#endif
	face->derivedEdgeDestinations[index] = cap;
}

/*
 * Serialisation read function
 */
Face *face_loadFromBinaryRepresentation(void **binaryString, Net * net)
{
	Face *face = calloc(1, sizeof(face));
	int32_t num, index;

	face = NULL;
	if (binaryRepresentation_peekNextElementType(*binaryString) ==
	    CODE_FACE) {
		binaryRepresentation_popNextElementType(binaryString);
		face->net = net;
		face->name = binaryRepresentation_getName(binaryString);
		num = binaryRepresentation_getInteger(binaryString);

		face->cardinal = num;
		face->topNodes = calloc(num, sizeof(Cap*));
		face->bottomNodes = calloc(num, sizeof(Cap**));
		face->bottomNodeNumbers = calloc(num, sizeof(int32_t));
		face->derivedEdgeDestinations = calloc(num, sizeof(Cap*));

		for (index = 0; index < num; index++)
			face_loadFromBinaryRepresentationAtIndex(binaryString, face, index, net);
	}

	return face;
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
 * Tests if a given node is a top node in face
 */
static int32_t face_isTopNode(Face * face, Cap * cap) {
	int32_t cardinal = face_getCardinal(face);
	int32_t index;

	for (index = 0; index < cardinal; index++)
		if (cap == face_getTopNode(face, index))
			return true;

	return false;
}

/*
 * Tests if all the top nodes are separate from the bottom nodes
 */
static int32_t face_isPseudoBipartite(Face * face) {
	int32_t cardinal = face_getCardinal(face);
	int32_t topIndex, bottomIndex;

	for (topIndex = 0; topIndex < cardinal; topIndex++)
		for (bottomIndex = 0; bottomIndex < face_getBottomNodeNumber(face, topIndex); bottomIndex++)
			if (face_isTopNode(face, face_getBottomNode(face, topIndex, bottomIndex)))
				return false;

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
	Cap * derivedDestination = NULL;

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
		else if (derivedDestination == NULL)
			derivedDestination = bottomNodePartnerAncestor;
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
	return face_isPseudoBipartite(face)
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

void face_check(Face *face) {
	assert(face_isSimple(face));
}

