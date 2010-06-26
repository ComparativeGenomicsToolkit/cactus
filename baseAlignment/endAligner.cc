
struct _AlignedPair {
    Name sequence1;
    int32_t position1;
    bool strand1;
    Name sequence2;
    int32_t position2;
    bool strand2;
    int32_t score;
} AlignedPair;

AlignedPair alignedPair_construct(Name sequence1, int32_t position1, bool strand1,
                                  Name sequence2, int32_t position2, bool strand2, int32_t score) {
    AlignedPair *alignedPair = st_malloc(sizeof(AlignedPair));
    alignedPair->sequence1 = sequence1;
    alignedPair->position1 = position1;
    alignedPair->strand1 = strand1;
    alignedPair->sequence2 = sequence2;
    alignedPair->position2 = position2;
    alignedPair->strand2 = strand2;
    alignedPair->score = score;
    return alignedPair;
}

void alignedPair_destruct(AlignedPair *alignedPair) {
    free(alignedPair);
}

void pairwiseAlignment(int32_t sequence1, const char *string1, int32_t sequence2,
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

stList *makeEndAlignemnt(End *end, int32_t maxSeqLength,
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
        int32_t score = stIntTuple_getPosition(alignedPair, 0);
        int32_t sequence1 = stIntTuple_getPosition(alignedPair, 1);
        int32_t position1 = stIntTuple_getPosition(alignedPair, 2);
        int32_t sequence2 = stIntTuple_getPosition(alignedPair, 3);
        int32_t position2 = stIntTuple_getPosition(alignedPair, 4);
        AdjacencySequence *adjacencySequence1 = stList_get(sequences, sequence1);
        AdjacencySequence *adjacencySequence2 = stList_get(sequences, sequence2);
        if (stPosetAlignment_isPossible(posetAlignment, sequence1, position1, sequence2, position2)) {
            stPosetAlignment_add(posetAlignment, sequence1, position1, sequence2, position2);
            //Add a converted version to the accepted aligned pairs.
            stList_append(acceptedAlignedPairs, alignedPair_construct(adjacencySequence1->sequenceName, adjacencySequence1->start + (adjacencySequence1->strand ? position1 : -position1), adjacencySequence1->strand,
                    adjacencySequence2->sequenceName, adjacencySequence2->start + (adjacencySequence2->strand ? position2 : -position2), adjacencySequence2->strand), score);
        }
        stIntTuple_destruct(alignedPair);
    }
    stList_destruct(alignedPairs);
    stPosetAlignment_destruct(posetAlignment);

    //Return the accepted pairs
    return acceptedAlignedPairs;
}

/*
 * This builds an adjacency list structure for the the sequences. Every sequence-position
 * has a column in the hash with which it can be aligned with.
 */
static stHash *buildAdjacencyList(stList *alignedPairs) {
    stHash *hash = stHash_construct3((uint32_t (*)(const void *))stIntTuple_hashKey,
            (int (*)(const void *, const void *))stIntTuple_equalsFn,
            (void (*)(void *))stIntTuple_destruct, NULL);
    stListIterator *it = stList_getIterator(pairs);
    stIntTuple *pair;
    while((pair = stList_getNext(it)) != NULL) {
       stIntTuple *seqPos1 = stIntTuple_construct(2, stIntTuple_getPosition(pair, 0), stIntTuple_getPosition(pair, 1));
       stIntTuple *seqPos2 = stIntTuple_construct(2, stIntTuple_getPosition(pair, 2), stIntTuple_getPosition(pair, 3));
       stSortedSet *column1 = stHash_search(hash, seqPos1);
       assert(column1 != NULL);
       stSortedSet *column2 = stHash_search(hash, seqPos2);
       assert(column2 != NULL);
       if(column1 != column2) { //Merge the columns
           stSortedSetIterator *it2 = stSortedSet_getIterator(column2);
           stIntTuple *seqPos3;
           while((seqPos3 = stSortedSet_getNext(it2)) != NULL) {
               assert(stSortedSet_search(column1, seqPos3) == NULL);
               stSortedSet_insert(column1, seqPos3);
               assert(stHash_search(hash, seqPos3) == column2);
               stHash_insert(hash, seqPos3, column1);
               assert(stHash_search(hash, seqPos3) == column1);
           }
           stSortedSet_destructIterator(it2);
           stSortedSet_destruct(column2);
       }
       //Cleanup loop.
       stIntTuple_destruct(seqPos1);
       stIntTuple_destruct(seqPos2);
    }
    stList_destructIterator(it);
    return hash;
}

int32_t getCutOffPoint(stList *inducedAlignment1, stList *inducedAlignment2,
        AdjacencySequence *adjacencySequence) {
    //Cummulate the score of the columns with respect to the sequences.
    int64_t *cummulativeAlignmentScore1 = cummulateScore(inducedAlignment1, adjacencySequence);
    int64_t *cummulativeAlignmentScore2 = cummulateScoreReverse(inducedAlignment2, adjacencySequence);

    //Choose the cutoff point
    int64_t maxScore = cummulativeAlignmentScore2[0];
    int32_t j = 0;
    for(int32_t i=1; i<adjacencySequence->length; i++) {
        int64_t k = cummulativeAlignmentScore1[i-1] + cummulativeAlignmentScore2[i];
        if(k >= maxScore) {
            maxScore = k;
            j = i;
        }
    }
    if(cummulativeAlignmentScore1[adjacencySequence->length-1] > maxScore) {
        j = adjacencySequence->length-1;
    }
    //Cleanup
    free(cummulativeAlignmentScore1);
    free(cummulativeAlignmentScore2);

    return j;
}

/*
 * Removes
 */
void pruneAlignments(stHash *endAlignment1, stHash *endAlignment2,
        Cap *cap) {
    if(endAlignment1 == endAlignment2) {
        return; //We ignore self loops, which have been dealt with already by the
        //posed aligner.
    }
    //Build a list of the columns in end alignment1 and end alignment 2 that contains the
    //sequence we're interested in.
    AdjacencySequence *adjacencySequence = cap_getAdjacency(cap);
    if(adjacencySequence->length == 0) {
        adjacencySequence_destruct(adjacencySequence);
        return;
    }
    stList *inducedAlignment1 = getInducedAlignment(endAlignment1, adjacencySequence);
    stList *inducedAlignment2 = getInducedAlignment(endAlignment2, adjacencySequence);

    int32_t cutOffPoint = getCutOffPoint(inducedAlignment1, inducedAlignment1, adjacencySequence);

    //Remove positions in the overlap according to the cutoff
    stListIterator *it = stList_getIterator(inducedAlignment1);
    stSortedSet *column;
    while((column = stList_getNext(column)) != NULL) {
        st
    }

}

stList *makeEndAlignments(Net *net) {
    //Make the end alignments, representing each as an adjacency alignment.
    End *end;
    Net_EndIterator *endIterator = net_getEndIterator(net);
    stHash *endAlignments = stHash_construct2(NULL, (void (*)(void *))stList_destruct);
    while((end = net_getNextEnd(net)) != NULL) {
        stList *endAlignment = makeEndAlignments(end);
        stHash_insert(end, buildAdjacencyList(endAlignment)); //adjacency list construction destroys the pairs.
        assert(stList_length(endAlignment) == 0);
        stList_destruct(endAlignment);
    }
    net_destructEndIterator(endIterator);

    //Prune end alignments
    endIterator = net_getEndIterator(net);
    while((end = net_getNextEnd(net)) != NULL) {
        stHash *endAlignment1 = stHash_search(end, endAlignments);
        assert(endAlignment1 != NULL);
        Cap *cap;
        End_InstanceIterator *capIterator = end_getInstanceIterator(end);
        while((cap = end_getNext(capIterator)) != NULL) {
            Cap *adjacentCap = cap_getAdjacency(cap);
            assert(adjacentCap != NULL);
            stHash *endAlignment2 = stHash_search(cap_getEnd(cap_getReverse(cap)));
            assert(endAlignment2 != NULL);
            //Now traverse the alignments, pruning stuff
            pruneAlignments(endAlignment1, endAlignment2, cap);
        }
        end_destructInstanceIterator(capIterator);
    }
    net_destructEndIterator(endIterator);

    //Now convert to set of final aligned pairs to return.
    stList *alignedPairs = stList_construct3(0, (void (*)(void *))alignedPair_destruct);
    stList *endAlignmentsList = stHash_getValues(endAlignments);
    while(stList_length(endAlignmentsList) > 0) {
        stHash *endAlignment = stList_pop(endAlignmentsList);
        stList *values = stHash_getValues(endAlignment);
        stSortedSet *seen = stSortedSet_construct();
        while(stList_length(values) > 0) {
            stList *column = stList_pop(values);
            if(stSortedSet_search(seen, column) == NULL) {
                stList_appendAll(alignedPairs, column);
                stList_destruct(column);
            }
        }
        stList_destruct(values);
        stSortedSet_destruct(seen);
    }

    //Cleanup
    stHash_destruct(endAlignments);
    stList_destruct(endAlignmentsList);

    return alignedPairs;
}

void cactus(stSortedSet *pairs, Net *net) {
    //Build pinch and adjacency graphs, fill out net.

}
