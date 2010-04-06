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
	binaryRepresentation_writeName(cap_getName(face->derivedEdgeDestinations[index]), writeFn);
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
	name = binaryRepresentation_getName(binaryString);
	cap = net_getCap(net, name);
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
