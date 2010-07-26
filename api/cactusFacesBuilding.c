#include "cactusGlobalsPrivate.h"

typedef struct _liftedEdge LiftedEdge;

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Ancestors function
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Utility function for the lifted edge hashtable
 */
static uint32_t buildFaces_hashfunction(const void *ptr) {
    Cap *key = (Cap *) ptr;
    return (uint32_t) cap_getName(key);
}

/*
 * Utility function for the lifted edge hashtable
 */
static int buildFaces_key_eq_fn(const void *ptrA, const void *ptrB) {
    return ptrA == ptrB;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Lifted edges functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct _liftedEdge {
    Cap *destination;
    Cap *bottomNode;
};

/*
 * Lifted edge destructor
 */
static void buildFaces_destructLiftedEdge(LiftedEdge * liftedEdge) {
    free(liftedEdge);
}
/*
 * Utility function for the lifted edge hashtable
 */
static void buildFaces_destructValue(void *ptr) {
    destructList((struct List *) ptr);
}

/*
 * Utility function for List struct
 */
static void buildFaces_destructListElem(void *ptr) {
    buildFaces_destructLiftedEdge((LiftedEdge *) ptr);
}

/*
 * Fill in a hashtable which to every node associates
 * a list of lifted edges
 */
static stHash *buildFaces_computeLiftedEdges(Net * net) {
    stHash *liftedEdgesTable = stHash_construct3(buildFaces_hashfunction,
            buildFaces_key_eq_fn, NULL, buildFaces_destructValue);
    Net_CapIterator *iter = net_getCapIterator(net);
    Cap *cap, *attachedAncestor;
    Cap *adjacency, *adjacencyAncestor;
    struct List *liftedEdges;
    LiftedEdge *liftedEdge;

    // Iterate through potential bottom nodes
    while ((cap = net_getNextCap(iter))) {
        // ... check if connected
        if ((adjacency = cap_getAdjacency(cap))) {
            // ... lift
            attachedAncestor = cap_getTopCap(cap);
            adjacencyAncestor = cap_getTopCap(cap_getPositiveOrientation(
                    adjacency));

#ifdef BEN_DEBUG
            assert((attachedAncestor && adjacencyAncestor) || (!attachedAncestor && !adjacencyAncestor));
#endif

            // If root node
            if (attachedAncestor == NULL)
                continue;

            // ... create lifted edge
            liftedEdge = st_malloc(sizeof(LiftedEdge));
            liftedEdge->destination = adjacencyAncestor;
            liftedEdge->bottomNode = cap;

#ifdef BEN_DEBUG
            // Self loop
            if (adjacencyAncestor == attachedAncestor)
                abort();
#endif

            // ... add it to the hashtable
            if ((liftedEdges
                    = stHash_search(liftedEdgesTable, attachedAncestor))) {
                listAppend(liftedEdges, liftedEdge);
            } else {
                liftedEdges = constructZeroLengthList(2,
                        buildFaces_destructListElem);
                listAppend(liftedEdges, liftedEdge);
                stHash_insert(liftedEdgesTable, attachedAncestor, liftedEdges);
            }
        }
    }

    net_destructCapIterator(iter);
    return liftedEdgesTable;
}

/*
 * Recursive function which fills a given list with the
 * connected nodes within a module
 */
static void buildFaces_fillTopNodeList(Cap * cap, struct List *list,
        stHash *liftedEdgesTable) {
    struct List *liftedEdges;
    int32_t index;

    // Limit of recursion
    if (listContains(list, cap))
        return;

    // Actual filling
    st_logInfo("Adding cap %p to face\n", cap);
    listAppend(list, cap);

    // Recursion through lifted edges
    if ((liftedEdges = stHash_search(liftedEdgesTable, cap)))
        for (index = 0; index < liftedEdges->length; index++)
            buildFaces_fillTopNodeList(
                    ((LiftedEdge *) liftedEdges-> list[index])->destination,
                    list, liftedEdgesTable);

    // Recursion through adjacency
    if (cap_getAdjacency(cap))
        buildFaces_fillTopNodeList(cap_getAdjacency(cap), list,
                liftedEdgesTable);

    // Remove from lifted edges table to prevent double usage
    // of end instances
    stHash_remove(liftedEdgesTable, cap);
    if (liftedEdges)
        destructList(liftedEdges);
}

/*
 * Constructs a face from a given Cap
 */
static void buildFaces_constructFromCap(Cap * startingCap,
        stHash *liftedEdgesTable, Net * net) {
    Face *face = face_construct(net);
    struct List *topNodes = constructZeroLengthList(16, NULL);
    struct List *liftedEdges;
    Cap *cap, *bottomNode, *ancestor;
    int32_t index, index2;

    printf("Constructing new face");

    // Establish list of top nodes
    buildFaces_fillTopNodeList(startingCap, topNodes, liftedEdgesTable);

#ifdef BEN_DEBUG
    // What, no top nodes!?
    if (topNodes->length == 0)
        abort();
#endif

    // Initialize data structure
    face_allocateSpace(face, topNodes->length);

    // For every top node
    for (index = 0; index < topNodes->length; index++) {
        cap = topNodes->list[index];
        face_setTopNode(face, index, cap);
        liftedEdges = stHash_search(liftedEdgesTable, cap);

        if (!liftedEdges) {
            face_setBottomNodeNumber(face, index, 0);
            continue;
        }

        face_setBottomNodeNumber(face, index, liftedEdges->length);
        // For every bottom node of that top node
        for (index2 = 0; index2 < liftedEdges->length; index2++) {
            bottomNode
                    = ((LiftedEdge *) liftedEdges-> list[index2])->bottomNode;
            face_addBottomNode(face, index, bottomNode);

#if BEN_DEBUG
            assert(cap_getAdjacency(bottomNode));
#endif
            ancestor = cap_getTopCap(cap_getPositiveOrientation(
                    cap_getAdjacency(bottomNode)));
            if (cap_getAdjacency(cap) != ancestor)
                face_setDerivedDestination(face, index, index2, ancestor);
            else
                face_setDerivedDestination(face, index, index2, NULL);

#ifdef BEN_DEBUG
            // If bottom nodes part of top nodes
            if (listContains(topNodes, cap_getPositiveOrientation(
                    ((LiftedEdge*) liftedEdges->list[index2])->bottomNode)))
                abort();
#endif
        }
    }

#ifdef BEN_DEBUG_ULTRA
    if (!buildFaces_isSimple(face))
    abort();
#endif

    // Clean up
    destructList(topNodes);
}

void net_reconstructFaces(Net * net) {
    net_destructFaces(net);
    stHash *liftedEdgesTable = buildFaces_computeLiftedEdges(net);
    Net_CapIterator *iter = net_getCapIterator(net);
    struct List *liftedEdges;
    Cap *current;

    while ((current = net_getNextCap(iter))) {
        if ((liftedEdges = stHash_search(liftedEdgesTable, current))
                && (liftedEdges->length >= 1)) {
            buildFaces_constructFromCap(current, liftedEdgesTable, net);
        }
    }
    stHash_destruct(liftedEdgesTable);
    net_destructCapIterator(iter);
}

void net_destructFaces(Net *net) {
    Face *face;
    while ((face = net_getFirstFace(net)) != NULL) {
        face_destruct(face);
    }
}

/*
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 * The following functions are all in aid of checking that the set of faces we have is well
 * formed.
 * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 */

static stHash *hashbottomCaps(Net *net) {
    /*
     * For each top node finds the corresponding set of bottom nodes and returns a
     * hash of top nodes to sets of bottom nodes.
     */
    stHash *bottomCapsHash = stHash_construct2(
            NULL, (void(*)(void *)) destructList);
    Event *rootEvent = eventTree_getRootEvent(net_getEventTree(net));
    Cap *cap;
    Net_CapIterator *capIterator = net_getCapIterator(net);

    while ((cap = net_getNextCap(capIterator)) != NULL) {
        assert(cap_getOrientation(cap));
        if (cap_getEvent(cap) != rootEvent && cap_getAdjacency(cap) != NULL) {
            Cap *cap2 = cap_getTopCap(cap);
            assert(cap2 != NULL);
            assert(cap_getOrientation(cap2));
            struct List *list;
            if ((list = stHash_search(bottomCapsHash, cap2)) == NULL) {
                list = constructEmptyList(0, NULL);
                stHash_insert(bottomCapsHash, cap2, list);
            }
            listAppend(list, cap);
        }
    }
    net_destructCapIterator(capIterator);

    return bottomCapsHash;
}

static stHash *computeLiftedEdges(stHash *bottomCapsHash) {
    /*
     * For each top node finds the set of top nodes connected to it by a lifted
     * edge. Returns a hash of top nodes to list of other top nodes
     * connected by lifted edges.
     */
    stHash *liftedEdgesHash = stHash_construct2(
            NULL, (void(*)(void *)) destructList);
    stHashIterator *iterator = stHash_getIterator(bottomCapsHash);
    Cap *topCap;
    while ((topCap = stHash_getNext(iterator)) != NULL) {
        assert(cap_getOrientation(topCap));
        struct List *bottomCaps = stHash_search(bottomCapsHash, topCap);
        assert(bottomCapsHash != NULL);
        assert(bottomCaps->length> 0);
        struct List *liftedEdges = constructEmptyList(0, NULL);
        assert(stHash_search(liftedEdgesHash, topCap) == NULL);
        stHash_insert(liftedEdgesHash, topCap, liftedEdges);
        int32_t i;
        for (i = 0; i < bottomCaps->length; i++) {
            Cap *bottomCap = bottomCaps->list[i];
            Cap *adjacentBottomCap = cap_getAdjacency(bottomCap);
            assert(adjacentBottomCap != NULL);
            adjacentBottomCap = cap_getPositiveOrientation(adjacentBottomCap);
            Cap *adjacentTopCap = cap_getTopCap(adjacentBottomCap);
            assert(adjacentTopCap != NULL);
            assert(cap_getOrientation(adjacentTopCap));
            assert(stHash_search(bottomCapsHash, adjacentTopCap) != NULL);
            listAppend(liftedEdges, adjacentTopCap);
        }
        assert(liftedEdges->length == bottomCaps->length);
    }
    stHash_destructIterator(iterator);
    return liftedEdgesHash;
}

static void computeModulesP(Cap *topCap, stHash *liftedEdgesHash,
        struct List *module, stHash *modulesHash) {
    int32_t i;
    Cap *adjacentTopCap;
    assert(cap_getOrientation(topCap));
    if (stHash_search(modulesHash, topCap) == NULL) {
        //Add to module
        stHash_insert(modulesHash, topCap, module);
        listAppend(module, topCap);

        //Traverse the lifted edges
        struct List *liftedEdges = stHash_search(liftedEdgesHash, topCap);
        if (liftedEdges != NULL) {
            assert(liftedEdges->length> 0);
            for (i = 0; i < liftedEdges->length; i++) {
                adjacentTopCap = liftedEdges->list[i];
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

static struct List *computeModules(stHash *liftedEdges) {
    /*
     * Finds the set of adjacency/lifted edge components, called modules,
     * and returns them in a list.
     */
    struct List *modules =
            constructEmptyList(0, (void(*)(void *)) destructList);
    stHash *modulesHash = stHash_construct();

    stHashIterator *iterator = stHash_getIterator(liftedEdges);
    Cap *topCap;
    while ((topCap = stHash_getNext(iterator)) != NULL) {
        assert(cap_getOrientation(topCap));
        if (stHash_search(modulesHash, topCap) == NULL) {
            struct List *module = constructEmptyList(0, NULL);
            computeModulesP(topCap, liftedEdges, module, modulesHash);
            listAppend(modules, module);
            assert(module->length >= 1); //with a self loop can have a module of length 1.
        }
    }
    stHash_destructIterator(iterator);
    stHash_destruct(modulesHash);
    return modules;
}

static void checkFace(struct List *module, stHash *bottomCapsHash) {
    /*
     * Checks a face.
     */
    //Checks the top nodes are all in one associated face.
    //Checks the set of bottom nodes for each face are in agreement.
    int32_t i, k;
    assert(module->length> 0);
    Cap *topCap = module->list[0];
    Face *face = cap_getTopFace(topCap);
    assert(face != NULL);
    assert(face_getCardinal(face) == module->length);

    for (i = 0; i < module->length; i++) {
        topCap = module->list[i];
        FaceEnd *faceEnd = cap_getTopFaceEnd(topCap);
        assert(faceEnd != NULL);
        assert(face == faceEnd_getFace(faceEnd));
        assert(faceEnd_getTopNode(faceEnd) == topCap);
        struct List *bottomCaps = stHash_search(bottomCapsHash, topCap);
        if (bottomCaps != NULL) { //could be null if top node has no lifted edges.
            for (k = 0; k < bottomCaps->length; k++) {
                Cap *bottomCap = bottomCaps->list[k];
                assert(cap_getBottomFaceEnd(bottomCap) == faceEnd);
            }
            //Temp debug output
            {
                st_uglyf(
                        "Number of caps in D's face end %i, in B's face end %i, topCap %i \n",
                        faceEnd_getNumberOfBottomNodes(faceEnd),
                        bottomCaps->length, topCap);
                FaceEnd_BottomNodeIterator *bottomNodeIterator =
                        faceEnd_getBottomNodeIterator(faceEnd);
                Cap *bottomCap;
                while ((bottomCap = faceEnd_getNextBottomNode(
                        bottomNodeIterator)) != NULL) {
                    st_uglyf("Bottom cap in Daniel's face: %i %i\n", bottomCap,
                            cap_getTopCap(bottomCap));
                }
                for (k = 0; k < bottomCaps->length; k++) {
                    bottomCap = bottomCaps->list[k];
                    st_uglyf("Bottom cap in Benedict's face: %i %i\n",
                            bottomCap, cap_getTopCap(bottomCap));
                }
            }
            //assert(faceEnd_getNumberOfBottomNodes(faceEnd) == bottomCaps->length);
        } else {
            assert(faceEnd_getNumberOfBottomNodes(faceEnd) == 0);
        }
    }
}

void face_checkFaces(Net *net) {
    /*
     * Checks that the set of faces is as we expect - with a face created
     * for each non-trivial face.
     */
    if (net_builtFaces(net)) { //only check the faces if they have been built..
        stHash *bottomCapsHash = hashbottomCaps(net);

        //Construct lifted edges
        stHash *liftedEdgesHash = computeLiftedEdges(bottomCapsHash);

        //Constructs lifted edge/adjacency edge connected components, called modules.
        //Faces are simply the nodes in the modules (the top nodes) and the set of
        //bottom nodes.
        struct List *modules = computeModules(liftedEdgesHash);

        //Check all faces we have computed are the same as those computed by Daniel.
        int32_t i;
        for (i = 0; i < modules->length; i++) {
            checkFace(modules->list[i], bottomCapsHash);
        }
        assert(modules->length == net_getFaceNumber(net)); //we should have checked exactly the number of faces.

        //Cleanup
        stHash_destruct(bottomCapsHash);
        stHash_destruct(liftedEdgesHash);
        destructList(modules);
    } else {
        //We do not like intermediate states.
        assert(net_getFaceNumber(net) == 0);
    }
}

