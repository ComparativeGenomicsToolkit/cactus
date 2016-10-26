#include <stdio.h>
#include <ctype.h>
#include "cactus.h"
#include "sonLib.h"

/*
 * Code to calculate a maximum likelihood (ML) string for a block using Felsenstein's pruning algorithm.
 */

/////
// Code to for creating a phylogenetic model of a given event tree with associated substitution matrices.
////
stMatrix *getSubMatrix(stTree *tree) {
    /*
     * Gets back the substitution matrix for the parent branch of a given node.
     */
    return ((void **) stTree_getClientData(tree))[0];
}

Event *getEvent(stTree *tree) {
    /*
     * Gets the event of the event tree that this node maps to.
     */
    return ((void **) stTree_getClientData(tree))[1];
}

static stTree *getPhylogeneticTree(Event *event, Event *eventToTreatAsParent,
        stMatrix *(*generateSubstitutionMatrix)(double)) {
    stTree *tree = stTree_construct();
    stMatrix *matrix = generateSubstitutionMatrix(
            event_getBranchLength(eventToTreatAsParent == NULL ? event : eventToTreatAsParent));
    void **attributes = st_malloc(sizeof(void *) * 2);
    attributes[0] = matrix;
    attributes[1] = event;
    stTree_setClientData(tree, attributes);
    for (int64_t i = 0; i < event_getChildNumber(event); i++) {
        if (eventToTreatAsParent != event_getChild(event, i)) {
            stTree_setParent(getPhylogeneticTree(event_getChild(event, i), NULL, generateSubstitutionMatrix), tree);
        }
    }
    return tree;
}

stTree *getPhylogeneticTreeRootedAtGivenEvent(Event *event, stMatrix *(*generateSubstitutionMatrix)(double)) {
    /*
     * Creates a stTree isomorphic to the eventTree that 'event' is part of, but rooted at 'event'.
     * Each node is the returned tree has two attributes, arranged in an array (see getSubMatrix and getEvent above).
     * The first is a substitution matrix giving substitution probabilities for bases along the incident parent branch of
     * the re-rooted tree.
     * The second is the event that it maps to in the original event tree.
     */
    stTree *tree = getPhylogeneticTree(event, NULL, generateSubstitutionMatrix); //This builds the subtree rooted at the given event
    stMatrix_destruct(getSubMatrix(tree)); //This cleans up the substitution matrix for the root of the remodeled tree.
    ((void **) stTree_getClientData(tree))[0] = generateSubstitutionMatrix(0.0); //And this parameterizes the substitution matrix of
    //the parent branch of the root to have zero length.

    //The following builds out the subtree of the eventTree not represented by tree
    Event *pEvent = NULL;
    stTree *tree2 = tree;
    while ((pEvent = event_getParent(event)) != NULL) {
        stTree *tree3 = getPhylogeneticTree(pEvent, event, generateSubstitutionMatrix);
        stTree_setParent(tree3, tree2);
        tree2 = tree3;
        event = pEvent;
    }
    return tree;
}

void cleanupPhylogeneticTreeP(stTree *tree) {
    for(int64_t i=0; i<stTree_getChildNumber(tree); i++) {
        cleanupPhylogeneticTreeP(stTree_getChild(tree, i));
    }
    stMatrix_destruct(getSubMatrix(tree));
    free(stTree_getClientData(tree));
}

void cleanupPhylogeneticTree(stTree *tree) {
    /*
     * Frees the phylogenetic tree created by getPhylogeneticTreeRootedAtGivenEvent, including the associated substitution matrices.
     */
    cleanupPhylogeneticTreeP(tree);
    stTree_destruct(tree);
}

stMatrix *generateJukesCantorMatrix(double distance) {
    /*
     * Generator function to make substitution matrices for the Jukes-Cantor model.
     */
    return stMatrix_jukesCantor(distance, 4);
}

/////
// The following calls the ML string for a block from a set of base probabilities.
/////

char indexToChar(int64_t i) {
    switch (i) {
    case 0:
        return 'A';
    case 1:
        return 'C';
    case 2:
        return 'G';
    case 3:
        return 'T';
    default:
        assert(0); //This should not happen
        return 'N';
    }
}

static char *getMaxLikelihoodString(double *baseProbs, int64_t length) {
    /*
     * For the "baseProbs" 2d array of base probabilities generates a ML string of bases.
     * The baseProbs array is organised as
     * [ Prob of A at position 0, Prob of C at position 0, Prob of G at position 0, Prob of T at position 0,
     *   Prob of A at position 1, Prob of C at position 1, Prob of G at position 1, Prob of T at position 1,
     *   ...
     *  etc.
     *  The returned string is a an upper case string of A, C, G and T.
     *  Length is the length of the string.
     *  In case of bases at a position with equal probability a (somewhat) random base is chosen.
     */
    char *mlString = st_malloc(sizeof(char) * (length+1));
    for (int64_t i = 0; i < length; i++) {
        int64_t k = 0;
        double m = baseProbs[i * 4];
        for (int64_t j = 1; j < 4; j++) {
            double n = baseProbs[i * 4 + j];
            if (n > m || (n == m && st_random() > 0.5)) {
                k = j;
                m = n;
            }
        }
        mlString[i] = indexToChar(k); //Convert the index of the ML base to a A,C,G,T character.
    }
    mlString[length] = '\0';
    return mlString;
}

///
// The following functions are the meat of the Felsenstein's algorithm implementation.
///

static double *transformBaseProbsBySubstitutionMatrix(double *baseProbs, int64_t length, stMatrix *substitutionMatrix) {
    /*
     * Updates the array of base probs, as described in getMaxLikelihoodString by multiplying the vector of base
     * probabilities at each position by the given substitution matrix.
     * Returns the input array.
     */
    for (int64_t i = 0; i < length; i++) {
        assert(stMatrix_n(substitutionMatrix) == 4);
        double *v = stMatrix_multiplySquareMatrixAndColumnVector(substitutionMatrix, &(baseProbs[i * 4]));
        memcpy(&(baseProbs[i * 4]), v, sizeof(double) * 4);
        free(v);
    }
    return baseProbs;
}

double *getEmptyBaseProbsString(int64_t length) {
    /*
     * Gets an array of base probs, as described in getMaxLikelihoodString,
     * for a block of 'length' positions, in which each position is initialised to 1.0.
     */
    double *baseProbs = st_calloc(length * 4, sizeof(double));
    for (int64_t i = 0; i < length * 4; i++) {
        baseProbs[i] = 1.0;
    }
    return baseProbs;
}

double *getBaseProbsString(char *string, int64_t length) {
    /*
     * Gets an array of base probs, as described in getMaxLikelihoodString, representing
     * the input string.
     */
    double *baseProbs = st_calloc(length * 4, sizeof(double)); //Gets the initial array initialised to 0.0 values
    for (int64_t i = 0; i < length; i++) {
        switch (toupper(string[i])) {
        case 'A':
            assert(baseProbs[i * 4] == 0.0);
            baseProbs[i * 4] = 1.0;
            break;
        case 'C':
            baseProbs[i * 4 + 1] = 1.0;
            break;
        case 'G':
            baseProbs[i * 4 + 2] = 1.0;
            break;
        case 'T':
            baseProbs[i * 4 + 3] = 1.0;
            break;
        default: //If N we treat marginalise over all possibilities.
            baseProbs[i * 4] = 1.0;
            baseProbs[i * 4 + 1] = 1.0;
            baseProbs[i * 4 + 2] = 1.0;
            baseProbs[i * 4 + 3] = 1.0;
            break;
        }
    }
    return baseProbs;
}

static void multiply(double *baseProbs1, double *baseProbs2, int64_t blockLength) {
    /*
     * Convenience function.
     * Updates baseProbs1, so that at each position i, baseProbs1[i] = baseProbs1[i] * baseProbs2[i], each
     * being the probability of a given base at a given position whose probability if the product of the initial probabilities.
     * Frees base probs2.
     */
    for (int64_t j = 0; j < blockLength * 4; j++) {
        baseProbs1[j] *= baseProbs2[j];
    }
    free(baseProbs2);
}

static double *computeBaseProbs(stTree *tree, stHash *eventsToStrings, int64_t blockLength) {
    /*
     * This is the Felsenstein's function to compute the probabilities of each base at each position of the block for the given root node of tree
     * (which is a phylogenetic tree and attached substitution matrices created by getSubstitutionTreeRootedAtGivenEvent).
     */
    //The code is recursive.
    if (stTree_getChildNumber(tree) > 0) { //Case root is an internal node.
        double *baseProbs = computeBaseProbs(stTree_getChild(tree, 0), eventsToStrings, blockLength);
        for (int64_t i = 1; i < stTree_getChildNumber(tree); i++) {
            multiply(baseProbs, computeBaseProbs(stTree_getChild(tree, i), eventsToStrings, blockLength), blockLength);
        }
        return transformBaseProbsBySubstitutionMatrix(baseProbs, blockLength, getSubMatrix(tree));
    } else { //Case root is a leaf
        double *baseProbs = getEmptyBaseProbsString(blockLength);
        stList *strings = stHash_search(eventsToStrings, getEvent(tree));
        if (strings != NULL) { //If there are strings associated with this event.
            for (int64_t i = 0; i < stList_length(strings); i++) {
                multiply(baseProbs,
                        transformBaseProbsBySubstitutionMatrix(getBaseProbsString(stList_get(strings, i), blockLength),
                                blockLength, getSubMatrix(tree)), blockLength);
            }
        }
        return baseProbs;
    }
}

////
// The following is used to soft-mask (make lower case) bases deemed to be repetitive in the source genomes.
////

void maskAncestralRepeatBases(Block *block, char *mlString) {
    /*
     * Soft masks the positions in the mlString that are deemed to be repetitive. A position is repetitive
     * if greater than 50% of the bases from which it is derived are not upper case.
     */

    int64_t *upperCounts = st_calloc(block_getLength(block), sizeof(int64_t)); //Counts of upper case bases at each position of the block.
    int64_t *nCounts = st_calloc(block_getLength(block), sizeof(int64_t)); //Counts of Ns at each position of the block.

    //Iterate through the sequences of the segments of a block and collate the number of upper case bases.
    Block_InstanceIterator *segmentIt = block_getInstanceIterator(block);
    Segment *segment;
    size_t numSegmentsWithSequence = 0;
    while ((segment = block_getNext(segmentIt)) != NULL) {
        if (segment_getSequence(segment) != NULL) {
            numSegmentsWithSequence++;
            char *string = segment_getString(segment);
            for (int64_t i = 0; i < block_getLength(block); i++) {
                char uC = toupper(string[i]);
                upperCounts[i] += uC == string[i] ? 1 : 0;
                nCounts[i] += (uC != 'A' && uC != 'C' && uC != 'G' && uC != 'T' ? 1 : 0);
            }
            free(string);
        }
    }
    block_destructInstanceIterator(segmentIt);

    //Convert any upper case character to lower case if the majority of bases
    //from which it is derived are not upper case.
    for (int64_t i = 0; i < block_getLength(block); i++) {
        if (nCounts[i] == numSegmentsWithSequence) {
            mlString[i] = 'N';
        }
        if (upperCounts[i] <= numSegmentsWithSequence / 2) {
            mlString[i] = tolower(mlString[i]);
        }
    }
    //Cleanup
    free(upperCounts);
    free(nCounts);
}

static stHash *hashEventsToSegmentStrings(Block *block) {
    /*
     * Returns a hash of events to the strings of segments with a given event.
     * The strings are stored in a list.
     */
    stHash *eventsToStrings = stHash_construct2(NULL, (void (*)(void *)) stList_destruct);
    Block_InstanceIterator *segmentIt = block_getInstanceIterator(block);
    Segment *segment;
    while ((segment = block_getNext(segmentIt)) != NULL) {
        if (segment_getSequence(segment) != NULL) {
            stList *strings = stHash_search(eventsToStrings, segment_getEvent(segment));
            if (strings == NULL) {
                strings = stList_construct3(0, free);
                stHash_insert(eventsToStrings, segment_getEvent(segment), strings);
            }
            stList_append(strings, segment_getString(segment));
        }
    }
    block_destructInstanceIterator(segmentIt);
    return eventsToStrings;
}

char *getMaximumLikelihoodString(stTree *tree, Block *block) {
    /*
     * Computes a maximum likelihood (ML) string for a given block.
     */
    char *mlString;
    if (block_getInstanceNumber(block) == 1
        && segment_getEvent(block_getFirst(block)) == getEvent(tree)) {
        // This block contains only one segment: the reference
        // segment. This is intended to be a "scaffold gap" of sorts
        // indicating that there is no direct support for the chosen
        // adjacency.
        mlString = malloc((block_getLength(block) + 1) * sizeof(char));
        memset(mlString, 'N', block_getLength(block));
        mlString[block_getLength(block)] = '\0';
    } else {
        stHash *eventsToStrings = hashEventsToSegmentStrings(block);
        double *baseProbs = computeBaseProbs(tree, eventsToStrings, block_getLength(block));
        mlString = getMaxLikelihoodString(baseProbs, block_getLength(block));
        //Cleanup
        free(baseProbs);
        stHash_destruct(eventsToStrings);
    }
    maskAncestralRepeatBases(block, mlString);
    return mlString;
}

