#include <stdio.h>
#include <ctype.h>
#include "cactus.h"
#include "sonLib.h"
#include "bioioC.h"

// OpenMP
#if defined(_OPENMP)
#include <omp.h>
#endif

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

// The whole alignment's tree, with the branch lengths reconstruction should use (see setReconstructionTree),
// and its nodes by name.  Set once, before any reconstruction, and only read after that.
// Likelihood vectors need these lengths too: they already carry the ancestor's uncertainty.
static stTree *reconstructionTree = NULL;
static stHash *reconstructionNodes = NULL;

static void addReconstructionNodes(stTree *tree) {
    if (stTree_getLabel(tree) != NULL) {
        stHash_insert(reconstructionNodes, (void *)stTree_getLabel(tree), tree);
    }
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        addReconstructionNodes(stTree_getChild(tree, i));
    }
}

void setReconstructionTree(const char *newick) {
    reconstructionTree = stTree_parseNewickString(newick);
    reconstructionNodes = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    addReconstructionNodes(reconstructionTree);
}

double getReconstructionBranchLength(Event *event) {
    /*
     * The event tree is the part of the whole tree that spans this alignment's genomes, with each
     * ancestor that is left with a single child removed and its branch added to the child's.  So the
     * branch above an event is the path from its node up to its parent's in the whole tree.
     */
    Event *parent = event_getParent(event);
    stTree *node = reconstructionNodes == NULL ? NULL : stHash_search(reconstructionNodes, (void *)event_getHeader(event));
    stTree *parentNode = node == NULL || parent == NULL ? NULL
                         : stHash_search(reconstructionNodes, (void *)event_getHeader(parent));
    if (parentNode != NULL) {
        double length = 0.0;
        for (; node != NULL && node != parentNode; node = stTree_getParent(node)) {
            length += stTree_getBranchLength(node); // INFINITY where the tree gives no length
        }
        if (node == parentNode && length != INFINITY) {
            return length;
        }
    }
    return event_getBranchLength(event);
}

static stTree *getPhylogeneticTree(Event *event, Event *eventToTreatAsParent,
        stMatrix *(*generateSubstitutionMatrix)(double)) {
    stTree *tree = stTree_construct();
    stMatrix *matrix = generateSubstitutionMatrix(
            getReconstructionBranchLength(eventToTreatAsParent == NULL ? event : eventToTreatAsParent));
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
    assert(stMatrix_n(substitutionMatrix) == 4);
    double v[4];
    for (int64_t i = 0; i < length; i++) {
        stMatrix_multiplySquareMatrixAndColumnVector2(substitutionMatrix, &(baseProbs[i * 4]), v);
        memcpy(&(baseProbs[i * 4]), v, sizeof(double) * 4);
    }
    return baseProbs;
}

double *getEmptyBaseProbsString(int64_t length) {
    /*
     * Gets an array of base probs, as described in getMaxLikelihoodString,
     * for a block of 'length' positions, in which each position is initialised to 1.0.
     */
    double *baseProbs = st_malloc(length * 4 * sizeof(double));
    for (int64_t i = 0; i < length * 4; i++) {
        baseProbs[i] = 1.0;
    }
    return baseProbs;
}

/////
// Likelihood vectors for ancestral sequences.
//
// A called ancestral base throws away how sure the call was, and the progressive alignment
// then treats that ancestor as a leaf when it builds the next one up.  Instead an ancestor can
// carry, at each position, P(bases of the genomes below it | its base) -- Felsenstein's upward
// message restricted to its children -- and its parent can use that in place of the one-hot
// vector of the called base.  Chained up the tree, this is exact pruning over all the leaves
// below each ancestor.  The vector must exclude the ancestor's outgroups: its parent sees
// those lineages again, and a posterior would count their evidence twice.
//
// Vectors are scaled so the largest of the four at each position is 1 (only their ratios
// matter), and are 1,1,1,1 where there is no information (e.g. scaffold gaps).
//
// File format, per sequence: a text line ">header\tlength\n", then length * 4 native float32s
// in A,C,G,T order.
/////

typedef struct {
    int64_t length;
    float *probs;
} LikelihoodVector;

static void likelihoodVector_destruct(LikelihoodVector *vector) {
    free(vector->probs);
    free(vector);
}

// Sequence -> LikelihoodVector for the input sequences whose bases are replaced by their
// likelihoods.  Set before the reference phase and only read during it.
static stHash *inputLikelihoods = NULL;

// The name the sequence is written under in the reference fasta (see getReferenceSequences.c),
// which is the name the parent will know it by.
static char *likelihoodSequenceName(Sequence *sequence) {
    const char *header = sequence_getHeader(sequence);
    if (strlen(header) > 0) {
        char *name = st_malloc(strlen(header) + 1);
        sscanf(header, "%s", name);
        return name;
    }
    return cactusMisc_nameToString(sequence_getName(sequence));
}

static void readLikelihoodFile(const char *path, stHash *nameToVector) {
    FILE *fh = st_fopen(path, "rb");
    int64_t bufSize = 1024;
    char *buf = st_malloc(bufSize);
    while (benLine(&buf, &bufSize, fh) != -1) {
        if (strlen(buf) == 0) {
            continue;
        }
        char *tab = strchr(buf, '\t');
        if (buf[0] != '>' || tab == NULL) {
            st_errAbort("Malformed header line in likelihood file %s: %s", path, buf);
        }
        *tab = '\0';
        LikelihoodVector *vector = st_malloc(sizeof(LikelihoodVector));
        if (sscanf(tab + 1, "%" SCNi64, &vector->length) != 1 || vector->length < 0) {
            st_errAbort("Malformed length in likelihood file %s for %s", path, buf + 1);
        }
        vector->probs = st_malloc(sizeof(float) * 4 * (vector->length > 0 ? vector->length : 1));
        if (fread(vector->probs, sizeof(float), 4 * vector->length, fh) != (size_t)(4 * vector->length)) {
            st_errAbort("Likelihood file %s is truncated in %s", path, buf + 1);
        }
        if (stHash_search(nameToVector, buf + 1) != NULL) {
            st_errAbort("Sequence %s appears twice in the likelihood files (second time in %s)", buf + 1, path);
        }
        stHash_insert(nameToVector, stString_copy(buf + 1), vector);
    }
    free(buf);
    st_fclose(fh, (char *)path);
}

void setInputAncestralLikelihoods(Flower *flower, stList *likelihoodFiles) {
    // name -> vector.  The vectors are handed over to inputLikelihoods as they are matched.
    stHash *nameToVector = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    for (int64_t i = 0; i < stList_length(likelihoodFiles); i++) {
        readLikelihoodFile(stList_get(likelihoodFiles, i), nameToVector);
    }

    // Resolve names to sequences once, so the base calling looks up by pointer
    inputLikelihoods = stHash_construct2(NULL, (void (*)(void *))likelihoodVector_destruct);
    Flower_SequenceIterator *seqIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(seqIt)) != NULL) {
        char *name = likelihoodSequenceName(sequence);
        LikelihoodVector *vector = stHash_search(nameToVector, name);
        if (vector != NULL) {
            if (vector->length != sequence_getLength(sequence)) {
                st_errAbort("Likelihood vector for %s has length %" PRIi64 " but the sequence has length %" PRIi64,
                            name, vector->length, sequence_getLength(sequence));
            }
            stHash_removeAndFreeKey(nameToVector, name); // frees the stored name, returns the vector
            stHash_insert(inputLikelihoods, sequence, vector);
        }
        free(name);
    }
    flower_destructSequenceIterator(seqIt);

    if (stHash_size(nameToVector) > 0) {
        stList *missing = stHash_getKeys(nameToVector);
        st_errAbort("%" PRIi64 " sequences in the likelihood files are not in the input, e.g. %s",
                    stList_length(missing), (char *)stList_get(missing, 0));
    }
    stHash_destruct(nameToVector);
    st_logInfo("Using likelihood vectors in place of the bases of %" PRIi64 " sequences\n", stHash_size(inputLikelihoods));
}

static double *getInputLikelihoods(Segment *segment) {
    /*
     * The loaded likelihoods of the segment, in the segment's orientation, or NULL if it has none.
     */
    if (inputLikelihoods == NULL) {
        return NULL;
    }
    LikelihoodVector *vector = stHash_search(inputLikelihoods, segment_getSequence(segment));
    if (vector == NULL) {
        return NULL;
    }
    bool strand = segment_getStrand(segment);
    int64_t start = segment_getStart(strand ? segment : segment_getReverse(segment)) -
                    sequence_getStart(segment_getSequence(segment));
    int64_t length = segment_getLength(segment);
    assert(start >= 0 && start + length <= vector->length);
    double *baseProbs = st_malloc(length * 4 * sizeof(double));
    for (int64_t i = 0; i < length; i++) {
        for (int64_t j = 0; j < 4; j++) {
            // on the reverse strand, position i is the complement of base (length-1-i), and in
            // A,C,G,T order the complement of base j is base 3-j
            baseProbs[i * 4 + j] = strand ? vector->probs[(start + i) * 4 + j]
                                          : vector->probs[(start + length - 1 - i) * 4 + 3 - j];
        }
    }
    return baseProbs;
}

double *getBaseProbsString(Segment *segment) {
    /*
     * Gets an array of base probs, as described in getMaxLikelihoodString, representing
     * the input string, or its likelihood vector if one was loaded.
     */
    double *likelihoods = getInputLikelihoods(segment);
    if (likelihoods != NULL) {
        return likelihoods;
    }
    char *string = segment_getString(segment);
    int64_t length = segment_getLength(segment);
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
    free(string);
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

static int getFirstSegmentMatchingEvent(const void *a, const void *b) {
    Event *e1 = (Event *)a, *e2 = segment_getEvent((Segment *)b);
    assert(e1 != NULL && e2 != NULL);
    // Must order events the same way sortByEvent does, as this searches the list
    // it produced
    return cactusMisc_nameCompare(event_getName(e1), event_getName(e2));
}

static double *computeBaseProbs(stTree *tree, stList *eventSortedSegments, int64_t blockLength) {
    /*
     * This is the Felsenstein's function to compute the probabilities of each base at each position of the block for the given root node of tree
     * (which is a phylogenetic tree and attached substitution matrices created by getSubstitutionTreeRootedAtGivenEvent).
     */
    //The code is recursive.
    if (stTree_getChildNumber(tree) > 0) { //Case root is an internal node.
        double *baseProbs = computeBaseProbs(stTree_getChild(tree, 0), eventSortedSegments, blockLength);
        int64_t i=1;
        // While there are no base probs, cos the subtree is empty replace base probs with those from another branch
        while(baseProbs == NULL && i < stTree_getChildNumber(tree)) {
            baseProbs = computeBaseProbs(stTree_getChild(tree, i++), eventSortedSegments, blockLength);
        }
        // Now that we have base probs combine the remaining branches
        while(i < stTree_getChildNumber(tree)) {
            double *baseProbs2 = computeBaseProbs(stTree_getChild(tree, i++), eventSortedSegments, blockLength);
            if(baseProbs2 != NULL) {
                multiply(baseProbs, baseProbs2, blockLength);
            }
        }
        return baseProbs == NULL ? NULL : transformBaseProbsBySubstitutionMatrix(baseProbs, blockLength, getSubMatrix(tree));
    } else { //Case root is a leaf
        Event *event = getEvent(tree);
        int64_t i = stList_binarySearchFirstIndex(eventSortedSegments, event, getFirstSegmentMatchingEvent);
        if(i == -1) {
            return NULL;
        }
        double *baseProbs = transformBaseProbsBySubstitutionMatrix(getBaseProbsString(stList_get(eventSortedSegments, i)),
                                                                   blockLength, getSubMatrix(tree));
        while(++i < stList_length(eventSortedSegments)) {
            Segment *segment = stList_get(eventSortedSegments, i);
            if(segment_getEvent(segment) != event) {
                break;
            }
            multiply(baseProbs, transformBaseProbsBySubstitutionMatrix(getBaseProbsString(segment), blockLength, getSubMatrix(tree)), blockLength);
        }
        return baseProbs;
    }
}

////
// The following is used to soft-mask (make lower case) bases deemed to be repetitive in the source genomes.
////

void maskAncestralRepeatBases(Block *block, stList *segments, char *mlString) {
    /*
     * Soft masks the positions in the mlString that are deemed to be repetitive. A position is repetitive
     * if greater than 50% of the bases from which it is derived are not upper case.
     */
    //assert(block_getInstanceNumber(block) == stList_length(segments));

    int64_t l = block_getLength(block), j = stList_length(segments);
    int64_t *upperCounts = st_calloc(l, sizeof(int64_t)); //Counts of upper case bases at each position of the block.
    int64_t *nCounts = st_calloc(l, sizeof(int64_t)); //Counts of Ns at each position of the block.

    //Iterate through the sequences of the segments of a block and collate the number of upper case bases.
    for(int64_t i=0; i<j; i++) {
        Segment *segment = stList_get(segments, i);
        assert(segment_getSequence(segment) != NULL);
        char *string = segment_getString(segment);
        for (int64_t k = 0; k < l; k++) {
            char uC = toupper(string[k]);
            upperCounts[k] += uC == string[k] ? 1 : 0;
            nCounts[k] += (uC != 'A' && uC != 'C' && uC != 'G' && uC != 'T' ? 1 : 0);
        }
        free(string);
    }

    //Convert any upper case character to lower case if the majority of bases
    //from which it is derived are not upper case.
    for (int64_t i = 0; i < l; i++) {
        if (nCounts[i] == j) {
            mlString[i] = 'N';
        }
        if (upperCounts[i] <= j / 2) {
            mlString[i] = tolower(mlString[i]);
        }
    }

    //Cleanup
    free(upperCounts);
    free(nCounts);
}

static int sortByEvent(const void *a, const void *b) {
    Event *e1 = segment_getEvent((Segment *)a), *e2 = segment_getEvent((Segment *)b);
    assert(e1 != NULL && e2 != NULL);
    // By name rather than by address: the addresses order the events differently
    // from one run to the next, and the order the segments are summed in decides
    // the ancestral bases.  Segment name breaks ties, so the order is total.
    int i = cactusMisc_nameCompare(event_getName(e1), event_getName(e2));
    if (i != 0) {
        return i;
    }
    return cactusMisc_nameCompare(segment_getName((Segment *)a), segment_getName((Segment *)b));
}

static stList *segmentsSortedByEvent(Block *block) {
    /*
     * Returns a list of segments in the block sorted by event
     */

    // Get the list of segments
    stList *segments = stList_construct();
    Block_InstanceIterator *segmentIt = block_getInstanceIterator(block);
    Segment *segment;
    while ((segment = block_getNext(segmentIt)) != NULL) {
        if (segment_getSequence(segment) != NULL) {
            stList_append(segments, segment);
        }
    }
    block_destructInstanceIterator(segmentIt);

    // Sort the segments by event
    stList_sort(segments, sortByEvent);

    return segments;
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
        stList *eventSortedSegments = segmentsSortedByEvent(block);
        double *baseProbs = computeBaseProbs(tree, eventSortedSegments, block_getLength(block));
        if(baseProbs == NULL) {
            baseProbs = getEmptyBaseProbsString(block_getLength(block));
        }
        mlString = getMaxLikelihoodString(baseProbs, block_getLength(block));
        maskAncestralRepeatBases(block, eventSortedSegments, mlString);
        //Cleanup
        free(baseProbs);
        stList_destruct(eventSortedSegments);
    }
    return mlString;
}

////
// The following computes and writes the likelihood vectors of the reference (ancestral) event.
////

static double *computeChildLikelihoods(stTree *tree, stList *eventSortedSegments, int64_t blockLength) {
    /*
     * The upward message at the root of the tree made by getPhylogeneticTreeRootedAtGivenEvent:
     * computeBaseProbs over its children, leaving out the child that carries the rest of the event
     * tree (the outgroups).  NULL if none of the children has a segment in the block.
     */
    Event *parentEvent = event_getParent(getEvent(tree));
    double *baseProbs = NULL;
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        stTree *child = stTree_getChild(tree, i);
        if (parentEvent != NULL && getEvent(child) == parentEvent) {
            continue;
        }
        double *childProbs = computeBaseProbs(child, eventSortedSegments, blockLength);
        if (childProbs == NULL) {
            continue;
        }
        if (baseProbs == NULL) {
            baseProbs = childProbs;
        } else {
            multiply(baseProbs, childProbs, blockLength);
        }
    }
    return baseProbs;
}

static void setReferenceLikelihoods(Block *block, Event *referenceEvent, stTree *tree, stHash *sequenceToVector) {
    stList *eventSortedSegments = segmentsSortedByEvent(block);
    int64_t length = block_getLength(block);
    double *baseProbs = NULL;
    bool computed = 0;
    for (int64_t k = 0; k < stList_length(eventSortedSegments); k++) {
        Segment *segment = stList_get(eventSortedSegments, k);
        if (segment_getEvent(segment) != referenceEvent) {
            continue;
        }
        LikelihoodVector *vector = stHash_search(sequenceToVector, segment_getSequence(segment));
        assert(vector != NULL);
        if (!computed) {
            computed = 1;
            baseProbs = computeChildLikelihoods(tree, eventSortedSegments, length);
            if (baseProbs == NULL) {
                break; // nothing below the ancestor here, so the vector stays at no information
            }
            for (int64_t i = 0; i < length; i++) {
                double m = 0.0;
                for (int64_t j = 0; j < 4; j++) {
                    m = baseProbs[i * 4 + j] > m ? baseProbs[i * 4 + j] : m;
                }
                for (int64_t j = 0; j < 4; j++) {
                    baseProbs[i * 4 + j] = m > 0.0 ? baseProbs[i * 4 + j] / m : 1.0;
                }
            }
        }
        bool strand = segment_getStrand(segment);
        int64_t start = segment_getStart(strand ? segment : segment_getReverse(segment)) -
                        sequence_getStart(segment_getSequence(segment));
        assert(start >= 0 && start + length <= vector->length);
        for (int64_t i = 0; i < length; i++) {
            for (int64_t j = 0; j < 4; j++) {
                // the inverse of the reverse strand mapping in getInputLikelihoods
                if (strand) {
                    vector->probs[(start + i) * 4 + j] = baseProbs[i * 4 + j];
                } else {
                    vector->probs[(start + length - 1 - i) * 4 + 3 - j] = baseProbs[i * 4 + j];
                }
            }
        }
    }
    free(baseProbs);
    stList_destruct(eventSortedSegments);
}

void writeAncestralLikelihoods(stList *flowerLayers, Flower *flower, Event *referenceEvent,
                               stMatrix *(*generateSubstitutionMatrix)(double), FILE *fileHandle) {
    stTree *tree = getPhylogeneticTreeRootedAtGivenEvent(referenceEvent, generateSubstitutionMatrix);

    // A vector for each sequence of the reference event, starting at no information
    stHash *sequenceToVector = stHash_construct2(NULL, (void (*)(void *))likelihoodVector_destruct);
    stList *sequences = stList_construct();
    Flower_SequenceIterator *seqIt = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while ((sequence = flower_getNextSequence(seqIt)) != NULL) {
        if (sequence_getEvent(sequence) == referenceEvent) {
            LikelihoodVector *vector = st_malloc(sizeof(LikelihoodVector));
            vector->length = sequence_getLength(sequence);
            vector->probs = st_malloc(sizeof(float) * 4 * (vector->length > 0 ? vector->length : 1));
            for (int64_t i = 0; i < 4 * vector->length; i++) {
                vector->probs[i] = 1.0;
            }
            stHash_insert(sequenceToVector, sequence, vector);
            stList_append(sequences, sequence);
        }
    }
    flower_destructSequenceIterator(seqIt);

    // Every ancestral base is in the reference segment of exactly one block, so the flowers can
    // fill in their parts of the vectors in parallel
    stList *flowers = stList_construct();
    for (int64_t i = 0; i < stList_length(flowerLayers); i++) {
        stList_appendAll(flowers, stList_get(flowerLayers, i));
    }
#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic)
#endif
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        // a flower has no block iterator, so visit each block from its 5' end
        Flower_EndIterator *endIt = flower_getEndIterator(stList_get(flowers, i));
        End *end;
        while ((end = flower_getNextEnd(endIt)) != NULL) {
            if (end_isBlockEnd(end)) {
                Block *block = block_getPositiveOrientation(end_getBlock(end));
                if (end_getPositiveOrientation(end) == block_get5End(block)) {
                    setReferenceLikelihoods(block, referenceEvent, tree, sequenceToVector);
                }
            }
        }
        flower_destructEndIterator(endIt);
    }
    stList_destruct(flowers);

    // The same sequences, under the same names, as getReferenceSequences writes
    for (int64_t i = 0; i < stList_length(sequences); i++) {
        sequence = stList_get(sequences, i);
        if (sequence_getLength(sequence) > 0 && !sequence_isTrivialSequence(sequence)) {
            LikelihoodVector *vector = stHash_search(sequenceToVector, sequence);
            char *name = likelihoodSequenceName(sequence);
            fprintf(fileHandle, ">%s\t%" PRIi64 "\n", name, vector->length);
            if (fwrite(vector->probs, sizeof(float), 4 * vector->length, fileHandle) != (size_t)(4 * vector->length)) {
                st_errAbort("Failed to write the likelihood vector of %s", name);
            }
            free(name);
        }
    }

    stList_destruct(sequences);
    stHash_destruct(sequenceToVector);
    cleanupPhylogeneticTree(tree);
}
