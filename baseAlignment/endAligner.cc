
void pairwiseAlignet(int32_t sequence1, const char *string1, int32_t sequence2,
        const char *string2, stList *alignedPairs) {
    vector<int> constraints;
    AlignerDPTable *fTable;
    AlignerDPTable *bTable;
    bfloat fProb;
    vector<double> gapPosteriorsX;
    vector<double> gapPosteriorsY;
    double divergence;
    struct PairwiseAlignmentInputParameters *pAIP =
            constructPairwiseAlignmentInputParameters();
    pAIP->matchThreshold = 0.01;
    XMLNode xMainNode = getDefaultModelParams();
    struct ModelParameters *modelParameters = readModelParameters(xMainNode);

    //Convert the sequences.
    char *convertedSeqX = convertSequence(seqX->string, convertFromDNA);
    char *convertedSeqY = convertSequence(seqY->string, convertFromDNA);

    //Do the alignment.
    vector<int> diagonals = computePairwiseAlignment(convertedSeqX,
            convertedSeqY, &modelParameters, pAIP, constraints, &fTable,
            &bTable, &fProb, gapPosteriorsX, gapPosteriorsY, &divergence);

    //Clean up.
    free(convertedSeqX);
    free(convertedSeqY);
    delete fTable;
    delete bTable;
    destructModelParameters(modelParameters);
    destructPairwiseAlignmentInputParameters(pAIP);
}

/*
 * Constructs a random spanning tree linking all the nodes in items into one component.
 */
void constructSpanningTree(int32_t numberOfSequences,
        stSortedSet *pairwiseAlignments) {
    stList *list = stList_construct();
    for (int32_t i = 0; i < numberOfSequences; i++) {
        stList_append(list, stIntTuple_construct(1, i));
    }
    while (stList_length(list) > 0) {
        stIntTuple *i = stList_removeItem(list, st_randomChoice(list));
        stIntTuple *j = st_randomChoice(list);
        stIntTuple *k = stIntTuple_construct(2, stIntTuple_getPosition(i, 0),
                stIntTuple_getPosition(j, 0));
        if (!stSortedSet_search(pairwiseAlignments, k)) {
            stSortedSet_insert(pairwiseAlignments, k);
        } else {
            stIntTuple_destruct(k);
        }
    }
    stList_destruct(list);
}

stSortedSet *makeEndAlignemnt(End *end, int32_t maxSeqLength,
        int32_t spanningTrees, struct PairwiseAlignmentInputParameters *pAIP) {
    //Get the sequences
    End_InstanceIterator *capIterator = end_getInstanceIterator(end);
    Cap *cap;
    stList *sequences = stList_construct3(0,
            (void(*)(void *)) destructAdjacencySequence);
    while ((cap = end_getNext(capIterator)) != NULL) {
        stList_append(sequences, getAdjacencySequenceAndCoordinates(cap,
                maxSeqLength));
    }
    end_destructInstanceIterator(capIterator);

    //Get the set of pairwise alignments (by constructing spanning trees)
    stSortedSet *pairwiseAlignments = stSortedSet_construct3((int(*)(
            const void *, const void *)) stIntTuple_cmpFn,
            (void(*)(void *)) stIntTuple_destruct);
    for (int i = 0; i < spanningTrees; i++) {
        constructSpanningTree(end_getInstanceNumber(end), pairwiseAlignments);
    }

    //Construct the alignments (using the Lunter code, output needs to be a set of weighted pairs above threshold)
    //and sort them by weight
    stSortedSetIterator *pairwiseAlignmentsIterator = stSortedSet_getIterator(
            pairwiseAlignments);
    stIntTuple *pairwiseAlignment;
    stList *alignedPairs = stList_construct();
    while ((pairwiseAlignment = stSortedSet_getNext(pairwiseAlignmentsIterator))
            != NULL) {
        int32_t sequence1 = stIntTuple_getPosition(pairwiseAlignment, 0);
        int32_t sequence2 = stIntTuple_getPosition(pairwiseAlignment, 1);
        AdjacencySequence *adjacencySequence1 =
                stList_get(sequences, sequence1);
        AdjacencySequence *adjacencySequence2 =
                stList_get(sequences, sequence2);
        pairwiseAlignet(sequence1, adjacencySequence1->string, sequence2,
                adjacencySequence2->string, alignedPairs);
    }
    stSortedSet_destructIterator(pairwiseAlignmentsIterator);

    //Sort the pairs (by weight)
    stList_sort(alignedPairs, stIntTuple_cmpFn);

    //Greedily construct poset and filter pairs..
    stPosetAlignment *posetAlignment = stPosetAlignment_construct(
            end_getInstanceNumber(end));
    stList *acceptedAlignedPairs = stList_construct3(0,
            (void(*)(void *)) stIntTuple_destruct);
    while (stList_length(alignedPairs) > 0) {
        stIntTuple *alignedPair = stList_pop(alignedPairs);
        int32_t sequence1 = stIntTuple_getPosition(alignedPair, 1);
        int32_t position1 = stIntTuple_getPosition(alignedPair, 2);
        int32_t sequence2 = stIntTuple_getPosition(alignedPair, 3);
        int32_t position2 = stIntTuple_getPosition(alignedPair, 4);
        if (stPosetAlignment_isPossible(posetAlignment, sequence1, position1, sequence2, position2)) {
            stPosetAlignment_add(posetAlignment, sequence1, position1, sequence2, position2);
            //Add a converted version to the accepted aligned pairs.

        }
        stIntTuple_destruct(alignedPair);
    }
    stList_destruct(alignedPairs);
    stPosetAlignment_destruct(posetAlignment);

    //Return the accepted pairs
    return acceptedAlignedPairs;
}

stSortedSet *makeEndAlignmentConsistent(stList *endAlignments) {
    //For each end alignment

    //Sort the pairs into a bucket for each sequence.

    //For each adjacency, get two alignment buckets.

    //Sum buckets cummulatively

    //Choose cutoff.. put all rejected pairs in a reject bucket..

    //filter all buckets by reject filter,

    //return set of filtered pairs (single list).
}

void cactus(stSortedSet *pairs, Net *net) {
    //Build pinch and adjacency graphs, fill out net.

}
