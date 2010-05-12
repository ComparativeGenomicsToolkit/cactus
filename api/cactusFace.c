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
	assert(face->faceEnds != NULL);
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

FaceEnd *face_getPreviousFaceEnd(Face_FaceEndIterator *iterator) {
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

FaceEnd *face_getFaceEndForTopNode(Face *face, Cap *cap) {
	cap = cap_getPositiveOrientation(cap);
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
		bottomNodePartnerAncestor = cap_getTopCap(bottomNodePartner);
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


void face_check(Face *face) {
	//Checks the structure of the face
	Face_FaceEndIterator *faceEndIterator = face_getFaceEndIterator(face);
	FaceEnd *faceEnd;
	int32_t i=0;
	while((faceEnd = face_getNextFaceEnd(faceEndIterator)) != NULL) {
		faceEnd_check(faceEnd);
		assert(faceEnd_getFace(faceEnd) == face);
		i++;
	}
	assert(face_getCardinal(face) == i);
	face_destructFaceEndIterator(faceEndIterator);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private face functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

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
	assert(face->cardinal == 0);
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
