/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"
#include "cactusFacePrivate.h"

/*
 * Basic constructor
 */
Face * face_construct(Flower * flower) {
    Face * face = st_calloc(1, sizeof(Face));
    face->flower = flower;
    flower_addFace(flower, face);
    return face;
}

/*
 * Basic memory deallocator
 */
void face_destruct(Face * face) {
    int64_t index;

    flower_removeFace(face->flower, face);

    //Remove reference to faces in the caps..
    FaceEnd *faceEnd;
    Face_FaceEndIterator *iterator = face_getFaceEndIterator(face);
    while ((faceEnd = face_getNextFaceEnd(iterator)) != NULL) {
	if (cap_getTopFace(faceEnd_getTopNode(faceEnd)) == face)
	    cap_setTopFace(faceEnd_getTopNode(faceEnd), NULL);
    }
    face_destructFaceEndIterator(iterator);

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
            if (face->faceEnds[index] != NULL) {
                faceEnd_destruct(face->faceEnds[index]);
            }
        }
        free(face->faceEnds);
    }
    free(face);
}

Flower *face_getFlower(Face *face) {
    return face->flower;
}

/*
 * Returns cardinal
 */
int64_t face_getCardinal(Face * face) {
    return face->cardinal;
}

/*
 * Get selected top node
 */
Cap * face_getTopNode(Face * face, int64_t index) {
    return face->topNodes[index];
}

Face_FaceEndIterator *face_getFaceEndIterator(Face *face) {
    Face_FaceEndIterator *iterator = st_malloc(sizeof(Face_FaceEndIterator));
    iterator->face = face;
    iterator->index = 0;
    return iterator;
}

static FaceEnd *face_getFaceEnd(Face *face, int64_t index) {
    assert(index >= 0);
    assert(index < face_getCardinal(face));
    assert(face->faceEnds != NULL);
    if (face->faceEnds[index] == NULL) {
        face->faceEnds[index] = faceEnd_construct(face, index);
    }
    return face->faceEnds[index];
}

FaceEnd *face_getNextFaceEnd(Face_FaceEndIterator *iterator) {
    if (iterator->index < face_getCardinal(iterator->face)) {
        return face_getFaceEnd(iterator->face, iterator->index++);
    }
    return NULL;
}

FaceEnd *face_getPreviousFaceEnd(Face_FaceEndIterator *iterator) {
    assert(iterator->index <= face_getCardinal(iterator->face));
    if (iterator->index > 0) {
        return face_getFaceEnd(iterator->face, --iterator->index);
    }
    return NULL;
}

Face_FaceEndIterator *face_copyFaceEndIterator(Face_FaceEndIterator *iterator) {
    Face_FaceEndIterator *iterator2 = st_malloc(sizeof(Face_FaceEndIterator));
    iterator2->face = iterator->face;
    iterator2->index = iterator->index;
    return iterator2;
}

void face_destructFaceEndIterator(Face_FaceEndIterator *iterator) {
    free(iterator);
}

FaceEnd *face_getFaceEndForTopNode(Face *face, Cap *cap) {
    cap = cap_getPositiveOrientation(cap);
    int64_t i;
    for (i = 0; i < face_getCardinal(face); i++) {
        if (face_getTopNode(face, i) == cap) {
            return face_getFaceEnd(face, i);
        }
    }
    assert(0);
    return NULL;
}

/*
 * Tests if all the top nodes are separate from the bottom nodes
 */
static int64_t face_isStronglyAcyclic(Face * face) {
    int64_t cardinal = face_getCardinal(face);
    int64_t topIndex1, topIndex, bottomIndex;
    Event * topEvent, *bottomEvent;

    for (topIndex1 = 0; topIndex1 < cardinal; topIndex1++) {
        topEvent = cap_getEvent(face_getTopNode(face, topIndex1));

        for (topIndex = 0; topIndex < cardinal; topIndex++) {
            for (bottomIndex = 0; bottomIndex < face_getBottomNodeNumber(face,
                    topIndex); bottomIndex++) {
                bottomEvent = cap_getEvent(face_getBottomNode(face, topIndex,
                        bottomIndex));
                if (bottomEvent == topEvent || event_isAncestor(bottomEvent,
                        topEvent))
                    return false;
            }
        }
    }

    return true;
}

/*
 * Tests if all the descent paths coming out of one top node in face are edge disjoint
 */
static int64_t face_hasMergedDescentEdgesAtIndex(Face * face, int64_t topIndex) {
    int64_t bottomNodeNumber = face_getBottomNodeNumber(face, topIndex);
    int64_t index;
    Cap * topNode = face_getTopNode(face, topIndex);
    Cap * current;
    struct List * list = constructZeroLengthList(100, NULL);

    for (index = 0; index < bottomNodeNumber; index++) {
        current = face_getBottomNode(face, topIndex, index);
        while (current != topNode) {
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
static int64_t face_hasSeparateDescentEdges(Face * face) {
    int64_t cardinal = face_getCardinal(face);
    int64_t index;

    for (index = 0; index < cardinal; index++)
        if (face_hasMergedDescentEdgesAtIndex(face, index))
            return false;

    return true;
}

/*
 * Tests if a given top node is part of a simple alternating cycle
 */
static int64_t face_breaksSimpleAlternatingPath(Face * face, int topIndex) {
    int64_t bottomNodeNumber = face_getBottomNodeNumber(face, topIndex);
    int64_t index;
    Cap * topNode = face_getTopNode(face, topIndex);
    Cap * partner = cap_getAdjacency(topNode);
    Cap * bottomNode;
    Cap * bottomNodePartner;
    Cap * bottomNodePartnerAncestor;
    Cap * derivedEdgeDestination = NULL;

    for (index = 0; index < bottomNodeNumber; index++) {
        bottomNode = face_getBottomNode(face, topIndex, index);
        bottomNodePartner = cap_getAdjacency(bottomNode);
#ifndef NDEBUG
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
static int64_t face_isSimpleAlternatingPath(Face * face) {
    int64_t cardinal = face_getCardinal(face);
    int64_t index;

    for (index = 0; index < cardinal; index++)
        if (face_breaksSimpleAlternatingPath(face, index))
            return false;

    return true;
}

/*
 * Tests if a face is simple
 */
int64_t face_isSimple(Face * face) {
    return face_isStronglyAcyclic(face) && face_hasSeparateDescentEdges(face)
            && face_isSimpleAlternatingPath(face);
}

/*
 * Tests if regular
 */
int64_t face_isRegular(Face * face) {
    int64_t index;

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
int64_t face_isCanonical(Face * face) {
    int64_t index;

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
    int64_t i = 0;
    while ((faceEnd = face_getNextFaceEnd(faceEndIterator)) != NULL) {
        faceEnd_check(faceEnd);
        cactusCheck(faceEnd_getFace(faceEnd) == face);
        i++;
    }
    cactusCheck(face_getCardinal(face) == i);
    face_destructFaceEndIterator(faceEndIterator);
}

/*
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * The following functions are all in aid of checking that the set of faces we have is well
 * formed.
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */

static stHash *hashbottomCaps(Flower *flower) {
    /*
     * For each top node finds the corresponding set of bottom nodes and returns a
     * hash of top nodes to sets of bottom nodes.
     */
    stHash *bottomCapsHash = stHash_construct2(
            NULL, (void(*)(void *)) stList_destruct);
    Event *rootEvent = eventTree_getRootEvent(flower_getEventTree(flower));
    Cap *cap;
    Flower_CapIterator *capIterator = flower_getCapIterator(flower);

    while ((cap = flower_getNextCap(capIterator)) != NULL) {
        assert(cap_getOrientation(cap));
        if (cap_getEvent(cap) != rootEvent && cap_getAdjacency(cap) != NULL) {
            Cap *cap2 = cap_getTopCap(cap);
            assert(cap2 != NULL);
            assert(cap_getOrientation(cap2));
            stList *list;
            if ((list = stHash_search(bottomCapsHash, cap2)) == NULL) {
               list = stList_construct3(0, NULL);
                stHash_insert(bottomCapsHash, cap2,list);
            }
            stList_append(list, cap);
        }
    }
    flower_destructCapIterator(capIterator);

    return bottomCapsHash;
}

static stHash *computeLiftedEdges(stHash *bottomCapsHash) {
    /*
     * For each top node finds the set of top nodes connected to it by a lifted
     * edge. Returns a hash of top nodes tolist of other top nodes
     * connected by lifted edges.
     */
    stHash *liftedEdgesHash = stHash_construct2(
            NULL, (void(*)(void *)) stList_destruct);
    stHashIterator *iterator = stHash_getIterator(bottomCapsHash);
    Cap *topCap;
    while ((topCap = stHash_getNext(iterator)) != NULL) {
        assert(cap_getOrientation(topCap));
        stList *bottomCaps = stHash_search(bottomCapsHash, topCap);
        assert(bottomCapsHash != NULL);
        assert(stList_length(bottomCaps)> 0);
        stList *liftedEdges = stList_construct3(0, NULL);
        assert(stHash_search(liftedEdgesHash, topCap) == NULL);
        stHash_insert(liftedEdgesHash, topCap, liftedEdges);
        int64_t i;
        for (i = 0; i < stList_length(bottomCaps); i++) {
            Cap *bottomCap = stList_get(bottomCaps, i);
            Cap *adjacentBottomCap = cap_getAdjacency(bottomCap);
            assert(adjacentBottomCap != NULL);
            adjacentBottomCap = cap_getPositiveOrientation(adjacentBottomCap);
            Cap *adjacentTopCap = cap_getTopCap(adjacentBottomCap);
            assert(adjacentTopCap != NULL);
            assert(cap_getOrientation(adjacentTopCap));
            assert(stHash_search(bottomCapsHash, adjacentTopCap) != NULL);
            stList_append(liftedEdges, adjacentTopCap);
        }
        assert(stList_length(liftedEdges) == stList_length(bottomCaps));
    }
    stHash_destructIterator(iterator);
    return liftedEdgesHash;
}

static void computeModulesP(Cap *topCap, stHash *liftedEdgesHash,
        stList *module, stHash *modulesHash) {
    int64_t i;
    Cap *adjacentTopCap;
    assert(cap_getOrientation(topCap));
    if (stHash_search(modulesHash, topCap) == NULL) {
        //Add to module
        stHash_insert(modulesHash, topCap, module);
        stList_append(module, topCap);

        //Traverse the lifted edges
        stList *liftedEdges = stHash_search(liftedEdgesHash, topCap);
        if (liftedEdges != NULL) {
            assert(stList_length(liftedEdges)> 0);
            for (i = 0; i < stList_length(liftedEdges); i++) {
                adjacentTopCap = stList_get(liftedEdges, i);
                assert(cap_getOrientation(adjacentTopCap));
                computeModulesP(adjacentTopCap, liftedEdgesHash, module,
                        modulesHash);
            }
        }

        //Traverse the direct adjaceny
        adjacentTopCap = cap_getAdjacency(topCap);
        if (adjacentTopCap != NULL) {
            adjacentTopCap = cap_getPositiveOrientation(adjacentTopCap);
            assert(adjacentTopCap != topCap);
            computeModulesP(adjacentTopCap, liftedEdgesHash, module,
                    modulesHash);
        }
    } else {
        assert(stHash_search(modulesHash, topCap) == module);
    }
}

static stList *computeModules(stHash *liftedEdges) {
    /*
     * Finds the set of adjacency/lifted edge components, called modules,
     * and returns them in alist.
     */
    stList *modules =
            stList_construct3(0, (void(*)(void *)) stList_destruct);
    stHash *modulesHash = stHash_construct();

    stHashIterator *iterator = stHash_getIterator(liftedEdges);
    Cap *topCap;
    while ((topCap = stHash_getNext(iterator)) != NULL) {
        assert(cap_getOrientation(topCap));
        if (stHash_search(modulesHash, topCap) == NULL) {
            stList *module = stList_construct3(0, NULL);
            computeModulesP(topCap, liftedEdges, module, modulesHash);
            stList_append(modules, module);
            assert(stList_length(module) >= 1); //with a self loop can have a module of length 1.
        }
    }
    stHash_destructIterator(iterator);
    stHash_destruct(modulesHash);
    return modules;
}

static void checkFace(stList *module, stHash *bottomCapsHash) {
#ifndef NDEBUG
    /*
     * Checks a face.
     */
    //Checks the top nodes are all in one associated face.
    //Checks the set of bottom nodes for each face are in agreement.
    int64_t i, k;
    assert(stList_length(module)> 0);
    Cap *topCap = stList_get(module, 0);
    Face *face = cap_getTopFace(topCap);
    assert(face != NULL);
    assert(face_getCardinal(face) == stList_length(module));

    for (i = 0; i < stList_length(module); i++) {
        topCap = stList_get(module, i);
        FaceEnd *faceEnd = cap_getTopFaceEnd(topCap);
        assert(faceEnd != NULL);
        assert(face == faceEnd_getFace(faceEnd));
        assert(faceEnd_getTopNode(faceEnd) == topCap);
        stList *bottomCaps = stHash_search(bottomCapsHash, topCap);
        if (bottomCaps != NULL) { //could be null if top node has no lifted edges.
            for (k = 0; k < stList_length(bottomCaps); k++) {
                Cap *bottomCap = stList_get(bottomCaps, k);
                assert(cap_getBottomFaceEnd(bottomCap) == faceEnd);
            }
            //Temp debug output
            {
                st_uglyf(
                        "Number of caps in D's face end %" PRIi64 ", in B's face end %" PRIi64 ", topCap %" PRIi64 " \n",
                        faceEnd_getNumberOfBottomNodes(faceEnd),
                        stList_length(bottomCaps), topCap);
                FaceEnd_BottomNodeIterator *bottomNodeIterator =
                        faceEnd_getBottomNodeIterator(faceEnd);
                Cap *bottomCap;
                while ((bottomCap = faceEnd_getNextBottomNode(
                        bottomNodeIterator)) != NULL) {
                    st_uglyf("Bottom cap in Daniel's face: %" PRIi64 " %" PRIi64 "\n", bottomCap,
                            cap_getTopCap(bottomCap));
                }
                for (k = 0; k < stList_length(bottomCaps); k++) {
                    bottomCap = stList_get(bottomCaps, k);
                    st_uglyf("Bottom cap in Benedict's face: %" PRIi64 " %" PRIi64 "\n",
                            bottomCap, cap_getTopCap(bottomCap));
                }
            }
            //assert(faceEnd_getNumberOfBottomNodes(faceEnd) == stList_length(bottomCaps));
        } else {
            assert(faceEnd_getNumberOfBottomNodes(faceEnd) == 0);
        }
    }
#endif
}


void face_checkFaces(Flower *flower) {
    /*
     * Checks that the set of faces is as we expect - with a face created
     * for each non-trivial face.
     */
    if (flower_builtFaces(flower)) { //only check the faces if they have been built..
        stHash *bottomCapsHash = hashbottomCaps(flower);

        //Construct lifted edges
        stHash *liftedEdgesHash = computeLiftedEdges(bottomCapsHash);

        //Constructs lifted edge/adjacency edge connected components, called modules.
        //Faces are simply the nodes in the modules (the top nodes) and the set of
        //bottom nodes.
        stList *modules = computeModules(liftedEdgesHash);

        //Check all faces we have computed are the same as those computed by Daniel.
        int64_t i;
        for (i = 0; i < stList_length(modules); i++) {
            checkFace(stList_get(modules, i), bottomCapsHash);
        }
        assert(stList_length(modules) == flower_getFaceNumber(flower)); //we should have checked exactly the number of faces.

        //Cleanup
        stHash_destruct(bottomCapsHash);
        stHash_destruct(liftedEdgesHash);
        stList_destruct(modules);
    } else {
        //We do not like intermediate states.
        assert(flower_getFaceNumber(flower) == 0);
    }
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
Cap * face_getDerivedDestinationAtIndex(Face * face, int64_t topIndex,
        int64_t derivedIndex) {
#ifndef NDEBUG
    assert(topIndex < face_getCardinal(face));
    assert(derivedIndex < face_getBottomNodeNumber(face, topIndex));
#endif
    return face->derivedEdgeDestinations[topIndex][derivedIndex];
}

/*
 * Get non-null derived destination of selected top node (useful for simple face)
 */
Cap * face_getDerivedDestination(Face * face, int64_t index) {
    int64_t bottomIndex;
    assert(index < face_getCardinal(face));
    for (bottomIndex = 0; bottomIndex < face_getBottomNodeNumber(face, index); bottomIndex++)
        if (face->derivedEdgeDestinations[index][bottomIndex])
            return face->derivedEdgeDestinations[index][bottomIndex];

    return NULL;
}

/*
 * Get the number of bottom nodes for the selected top node
 */
int64_t face_getBottomNodeNumber(Face * face, int64_t topIndex) {
    return face->bottomNodeNumbers[topIndex];
}

/*
 * Get selected bottom node from selected top node in face
 */
Cap * face_getBottomNode(Face * face, int64_t topNodeIndex,
        int64_t bottomNodeIndex) {
    return face->bottomNodes[topNodeIndex][bottomNodeIndex];
}

/*
 * Allocate arrays to allow for data
 */
void face_allocateSpace(Face * face, int64_t cardinal) {
    assert(face->cardinal == 0);
    face->cardinal = cardinal;
    face->topNodes = st_calloc(cardinal, sizeof(Cap *));
    face->bottomNodes = st_calloc(cardinal, sizeof(Cap **));
    face->bottomNodeNumbers = st_calloc(cardinal, sizeof(int64_t));
    face->derivedEdgeDestinations = st_calloc(cardinal, sizeof(Cap **));
    face->faceEnds = st_calloc(cardinal, sizeof(FaceEnd *));
}

/*
 * Sets the selected top node
 */
void face_setTopNode(Face * face, int64_t topIndex, Cap * topNode) {
    if (topNode)
        topNode = cap_getPositiveOrientation(topNode);
    face->topNodes[topIndex] = topNode;
    cap_setTopFace(topNode, face);
}

/*
 * Set bottom node count and allocate space to store pointers
 */
void face_setBottomNodeNumber(Face * face, int64_t topIndex, int64_t number) {
    face->bottomNodeNumbers[topIndex] = 0;
    if (number) {
        face->bottomNodes[topIndex] = st_calloc(number, sizeof(Cap *));
        face->derivedEdgeDestinations[topIndex] = st_calloc(number,
                sizeof(Cap *));
    } else {
        face->bottomNodes[topIndex] = NULL;
        face->derivedEdgeDestinations[topIndex] = NULL;
    }
}

/*
 * Sets the derived edge destination for a given top node in face
 */
void face_setDerivedDestination(Face * face, int64_t topIndex,
        int64_t bottomIndex, Cap * destination) {
    if (destination)
        destination = cap_getPositiveOrientation(destination);
    face->derivedEdgeDestinations[topIndex][bottomIndex] = destination;
}

/*
 * Adds bottom node to selected top node in face
 */
void face_addBottomNode(Face * face, int64_t topIndex, Cap * bottomNode) {
    if (bottomNode)
        bottomNode = cap_getPositiveOrientation(bottomNode);
    face->bottomNodes[topIndex][face->bottomNodeNumbers[topIndex]++]
            = bottomNode;
}
