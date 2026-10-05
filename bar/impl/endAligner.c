/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "endAligner.h"
#include <math.h>
#include "multipleAligner.h"
#include "adjacencySequences.h"
#include "pairwiseAligner.h"

AlignedPair *alignedPair_construct(int64_t subsequenceIdentifier1, int64_t position1, bool strand1,
        int64_t subsequenceIdentifier2, int64_t position2, bool strand2, int64_t score, int64_t rScore) {
    AlignedPair *alignedPair = st_malloc(sizeof(AlignedPair));
    alignedPair->subsequenceIdentifier = subsequenceIdentifier1;
    alignedPair->position = position1;
    alignedPair->strand = strand1;

    alignedPair->reverse = st_malloc(sizeof(AlignedPair));
    alignedPair->reverse->reverse = alignedPair;

    alignedPair->reverse->subsequenceIdentifier = subsequenceIdentifier2;
    alignedPair->reverse->position = position2;
    alignedPair->reverse->strand = strand2;

    alignedPair->score = score;
    alignedPair->reverse->score = rScore;

    return alignedPair;
}

void alignedPair_destruct(AlignedPair *alignedPair) {
    free(alignedPair); //We assume the reverse will be free independently.
}

static int alignedPair_cmpFnP(const AlignedPair *alignedPair1, const AlignedPair *alignedPair2) {
    int i = cactusMisc_nameCompare(alignedPair1->subsequenceIdentifier, alignedPair2->subsequenceIdentifier);
    if(i == 0) {
        i = alignedPair1->position > alignedPair2->position ? 1 : (alignedPair1->position < alignedPair2->position ? -1 : 0);
        if(i == 0) {
            i = alignedPair1->strand == alignedPair2->strand ? 0 : (alignedPair1->strand ? 1 : -1);
        }
    }
    return i;
}

int alignedPair_cmpFn(const AlignedPair *alignedPair1, const AlignedPair *alignedPair2) {
    int i = alignedPair_cmpFnP(alignedPair1, alignedPair2);
    if(i == 0) {
        i = alignedPair_cmpFnP(alignedPair1->reverse, alignedPair2->reverse);
    }
    return i;
}

struct _barDistances {
    stTree *tree;   // the reconstruction tree, or NULL
    stHash *nodes;  // its nodes by label
};

static void barDistances_addNodes(stHash *nodes, stTree *tree) {
    if (stTree_getLabel(tree) != NULL) {
        stHash_insert(nodes, (void *)stTree_getLabel(tree), tree);
    }
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        barDistances_addNodes(nodes, stTree_getChild(tree, i));
    }
}

BarDistances *barDistances_constructFromCactusParams(CactusParams *params) {
    BarDistances *distances = st_calloc(1, sizeof(BarDistances));
    if (cactusParams_has(params, 2, "reference", "reconstructionTree")) {
        char *newick = cactusParams_get_string(params, 2, "reference", "reconstructionTree");
        distances->tree = stTree_parseNewickString(newick);
        distances->nodes = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
        barDistances_addNodes(distances->nodes, distances->tree);
        free(newick);
    }
    return distances;
}

bool barDistances_hasReconstructionTree(BarDistances *distances) {
    return distances->tree != NULL;
}

double barDistances_get(BarDistances *distances, Event *event1, Event *event2) {
    stTree *node1 = distances->nodes == NULL ? NULL : stHash_search(distances->nodes, (void *)event_getHeader(event1));
    stTree *node2 = distances->nodes == NULL ? NULL : stHash_search(distances->nodes, (void *)event_getHeader(event2));
    if (node1 != NULL && node2 != NULL) {
        double d1 = 0.0;
        for (stTree *a = node1; a != NULL; a = stTree_getParent(a)) {
            double d2 = 0.0;
            for (stTree *b = node2; b != NULL; b = stTree_getParent(b)) {
                if (a == b) {
                    return d1 + d2;
                }
                d2 += stTree_getBranchLength(b);
            }
            d1 += stTree_getBranchLength(a);
        }
    }
    Event *ancestor = eventTree_getCommonAncestor(event1, event2);
    double d = 0.0;
    for (Event *e = event1; e != ancestor; e = event_getParent(e)) {
        d += event_getBranchLength(e);
    }
    for (Event *e = event2; e != ancestor; e = event_getParent(e)) {
        d += event_getBranchLength(e);
    }
    return d;
}

void barDistances_destruct(BarDistances *distances) {
    if (distances->tree != NULL) {
        stHash_destruct(distances->nodes);
        stTree_destruct(distances->tree);
    }
    free(distances);
}

void barModel_readGenomes(CactusParams *params, const char *section, const char *name, char **genomes) {
    genomes[0] = genomes[1] = NULL;
    if (!cactusParams_has(params, 4, "bar", section, name, "genomes")) {
        return;
    }
    char *value = cactusParams_get_string(params, 4, "bar", section, name, "genomes");
    stList *tokens = stString_split(value);
    if (stList_length(tokens) != 2) {
        st_errAbort("bar model %s/%s has genomes=\"%s\", which is not two genome names", section, name, value);
    }
    genomes[0] = stString_copy(stList_get(tokens, 0));
    genomes[1] = stString_copy(stList_get(tokens, 1));
    stList_destruct(tokens);
    free(value);
}

bool barModel_isForGenomes(char **genomes, const char *a, const char *b) {
    return genomes[0] != NULL && ((strcmp(genomes[0], a) == 0 && strcmp(genomes[1], b) == 0) ||
                                  (strcmp(genomes[0], b) == 0 && strcmp(genomes[1], a) == 0));
}

struct _pecanPairModels {
    int64_t modelNo;
    double *maxDistance; // NAN for a model with none
    char **genomes; // model m's two genomes at 2m and 2m+1, NULL for a model with none
    StateMachine **stateMachines;
    BarDistances *distances;
};

PecanPairModels *pecanPairModels_constructFromCactusParams(CactusParams *params) {
    int64_t modelNo = cactusParams_has(params, 3, "bar", "pecan", "pecanPairModels") ?
        cactusParams_get_int(params, 3, "bar", "pecan", "pecanPairModels") : 0;
    if (modelNo <= 0) {
        return NULL;
    }
    PecanPairModels *pairModels = st_calloc(1, sizeof(PecanPairModels));
    pairModels->modelNo = modelNo;
    pairModels->maxDistance = st_malloc(sizeof(double) * modelNo);
    pairModels->genomes = st_malloc(sizeof(char *) * 2 * modelNo);
    pairModels->stateMachines = st_malloc(sizeof(StateMachine *) * modelNo);
    int64_t keyed = 0;
    double lastMaxDistance = -INFINITY;
    for (int64_t m = 0; m < modelNo; m++) {
        char *name = stString_print("pairModel%" PRIi64, m);
        barModel_readGenomes(params, "pecan", name, pairModels->genomes + 2 * m);
        keyed += pairModels->genomes[2 * m] != NULL;
        pairModels->maxDistance[m] = NAN;
        if (cactusParams_has(params, 4, "bar", "pecan", name, "maxDistance")) {
            pairModels->maxDistance[m] = cactusParams_get_float(params, 4, "bar", "pecan", name, "maxDistance");
            if (pairModels->maxDistance[m] < lastMaxDistance) {
                st_errAbort("pecan pair models must be in ascending maxDistance, but %s has %f after %f", name,
                            pairModels->maxDistance[m], lastMaxDistance);
            }
            lastMaxDistance = pairModels->maxDistance[m];
        } else if (pairModels->genomes[2 * m] == NULL) {
            st_errAbort("pecan pair model %s has neither genomes nor a maxDistance", name);
        }
        char *json = cactusParams_get_string(params, 4, "bar", "pecan", name, "hmm");
        Hmm *hmm = hmm_jsonParse(json, strlen(json));
        if (hmm->type != fiveState) {
            st_errAbort("pecan pair model %s is a type %i HMM; bar's pecan uses the five-state model (type %i)",
                        name, (int)hmm->type, (int)fiveState);
        }
        pairModels->stateMachines[m] = hmm_getStateMachine(hmm);
        hmm_destruct(hmm);
        free(json);
        free(name);
    }
    pairModels->distances = barDistances_constructFromCactusParams(params);
    st_logInfo("bar: %" PRIi64 " pecan pair models (%" PRIi64 " for pairs of genomes), distances on the %s tree\n",
               modelNo, keyed, barDistances_hasReconstructionTree(pairModels->distances) ? "reconstruction" : "event");
    return pairModels;
}

void pecanPairModels_destruct(PecanPairModels *pairModels) {
    if (pairModels == NULL) {
        return;
    }
    for (int64_t m = 0; m < pairModels->modelNo; m++) {
        stateMachine_destruct(pairModels->stateMachines[m]);
    }
    for (int64_t i = 0; i < 2 * pairModels->modelNo; i++) {
        free(pairModels->genomes[i]);
    }
    free(pairModels->genomes);
    free(pairModels->stateMachines);
    free(pairModels->maxDistance);
    barDistances_destruct(pairModels->distances);
    free(pairModels);
}

/*
 * cPecan's per-pair callback: the model for the two sequences' genomes, else for the distance between them.
 */
typedef struct {
    PecanPairModels *pairModels;
    stList *events; // the genome of each sequence, by its index in the alignment
} PairModelChoice;

static StateMachine *choosePairModel(int64_t seqX, int64_t seqY, void *extraArgs) {
    PairModelChoice *choice = extraArgs;
    Event *eventX = stList_get(choice->events, seqX), *eventY = stList_get(choice->events, seqY);
    if (eventX == eventY) {
        return NULL;
    }
    PecanPairModels *pairModels = choice->pairModels;
    for (int64_t m = 0; m < pairModels->modelNo; m++) {
        if (barModel_isForGenomes(pairModels->genomes + 2 * m, event_getHeader(eventX), event_getHeader(eventY))) {
            return pairModels->stateMachines[m];
        }
    }
    double distance = barDistances_get(pairModels->distances, eventX, eventY);
    for (int64_t m = 0; m < pairModels->modelNo; m++) {
        if (distance <= pairModels->maxDistance[m]) { // false for a NAN maxDistance
            return pairModels->stateMachines[m];
        }
    }
    return NULL;
}

stSortedSet *makeEndAlignment(StateMachine *sM, End *end, int64_t spanningTrees, int64_t maxSequenceLength,
        bool useProgressiveMerging, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    return makeEndAlignmentWithPairModels(sM, NULL, end, spanningTrees, maxSequenceLength, useProgressiveMerging, gapGamma,
                                          pairwiseAlignmentBandingParameters);
}

stSortedSet *makeEndAlignmentWithPairModels(StateMachine *sM, PecanPairModels *pairModels, End *end, int64_t spanningTrees,
        int64_t maxSequenceLength, bool useProgressiveMerging, float gapGamma,
        PairwiseAlignmentParameters *pairwiseAlignmentBandingParameters) {
    //Make an alignment of the sequences in the ends

    //Get the adjacency sequences to be aligned.
    Cap *cap;
    End_InstanceIterator *it = end_getInstanceIterator(end);
    stList *sequences = stList_construct3(0, (void (*)(void *))adjacencySequence_destruct);
    stList *seqFrags = stList_construct3(0, (void (*)(void *))seqFrag_destruct);
    stList *events = stList_construct();
    stHash *endInstanceNumbers = stHash_construct2(NULL, free);
    while((cap = end_getNext(it)) != NULL) {
        if(cap_getSide(cap)) {
            cap = cap_getReverse(cap);
        }
        AdjacencySequence *adjacencySequence = adjacencySequence_construct(cap, maxSequenceLength);
        stList_append(sequences, adjacencySequence);
        assert(cap_getAdjacency(cap) != NULL);
        End *otherEnd = end_getPositiveOrientation(cap_getEnd(cap_getAdjacency(cap)));
        stList_append(seqFrags, seqFrag_construct(adjacencySequence->string, 0, end_getName(otherEnd)));
        stList_append(events, cap_getEvent(cap));
        //Increase count of seqfrags with a given end.
        int64_t *c = stHash_search(endInstanceNumbers, otherEnd);
        if(c == NULL) {
            c = st_calloc(1, sizeof(int64_t));
            assert(*c == 0);
            stHash_insert(endInstanceNumbers, otherEnd, c);
        }
        (*c)++;
    }
    end_destructInstanceIterator(it);

    //Get the alignment.
    PairModelChoice choice = { pairModels, events };
    MultipleAlignment *mA  = makeAlignmentWithPairStateMachines(sM, pairModels == NULL ? NULL : choosePairModel, &choice,
            seqFrags, spanningTrees, 100000000, useProgressiveMerging, gapGamma, pairwiseAlignmentBandingParameters);
    stList_destruct(events);

    //Build an array of weights to reweight pairs in the alignment.
    int64_t *pairwiseAlignmentsPerSequenceNonCommonEnds = st_calloc(stList_length(seqFrags), sizeof(int64_t));
    int64_t *pairwiseAlignmentsPerSequenceCommonEnds = st_calloc(stList_length(seqFrags), sizeof(int64_t));
    //First build array on number of pairwise alignments to each sequence, distinguishing alignments between sequences sharing
    //common ends.
    for(int64_t i=0; i<stList_length(mA->chosenPairwiseAlignments); i++) {
        stIntTuple *pairwiseAlignment = stList_get(mA->chosenPairwiseAlignments, i);
        int64_t seq1 = stIntTuple_get(pairwiseAlignment, 1);
        int64_t seq2 = stIntTuple_get(pairwiseAlignment, 2);
        assert(seq1 != seq2);
        SeqFrag *seqFrag1 = stList_get(seqFrags, seq1);
        SeqFrag *seqFrag2 = stList_get(seqFrags, seq2);
        int64_t *pairwiseAlignmentsPerSequence = seqFrag1->rightEndId == seqFrag2->rightEndId
                ? pairwiseAlignmentsPerSequenceCommonEnds : pairwiseAlignmentsPerSequenceNonCommonEnds;
        pairwiseAlignmentsPerSequence[seq1]++;
        pairwiseAlignmentsPerSequence[seq2]++;
    }
    //Now calculate score adjustments.
    double *scoreAdjustmentsNonCommonEnds = st_malloc(stList_length(seqFrags) * sizeof(double));
    double *scoreAdjustmentsCommonEnds = st_malloc(stList_length(seqFrags) * sizeof(double));
    for(int64_t i=0; i<stList_length(seqFrags); i++) {
        SeqFrag *seqFrag = stList_get(seqFrags, i);
        End *otherEnd = flower_getEnd(end_getFlower(end), seqFrag->rightEndId);
        assert(otherEnd != NULL);
        assert(stHash_search(endInstanceNumbers, otherEnd) != NULL);
        int64_t commonInstanceNumber = *(int64_t *)stHash_search(endInstanceNumbers, otherEnd);
        int64_t nonCommonInstanceNumber = stList_length(seqFrags) - commonInstanceNumber;

        assert(commonInstanceNumber > 0 && nonCommonInstanceNumber >= 0);
        assert(pairwiseAlignmentsPerSequenceNonCommonEnds[i] <= nonCommonInstanceNumber);
        assert(pairwiseAlignmentsPerSequenceNonCommonEnds[i] >= 0);
        assert(pairwiseAlignmentsPerSequenceCommonEnds[i] < commonInstanceNumber);
        assert(pairwiseAlignmentsPerSequenceCommonEnds[i] >= 0);

        //scoreAdjustmentsNonCommonEnds[i] = ((double)nonCommonInstanceNumber + commonInstanceNumber - 1)/(pairwiseAlignmentsPerSequenceNonCommonEnds[i] + pairwiseAlignmentsPerSequenceCommonEnds[i]);
        //scoreAdjustmentsCommonEnds[i] = scoreAdjustmentsNonCommonEnds[i];
        if(pairwiseAlignmentsPerSequenceNonCommonEnds[i] > 0) {
          scoreAdjustmentsNonCommonEnds[i] = ((double)nonCommonInstanceNumber)/pairwiseAlignmentsPerSequenceNonCommonEnds[i];
          assert(scoreAdjustmentsNonCommonEnds[i] >= 1.0);
          assert(scoreAdjustmentsNonCommonEnds[i] <= nonCommonInstanceNumber);
        }
        else {
          scoreAdjustmentsNonCommonEnds[i] = INT64_MIN;
        }
        if(pairwiseAlignmentsPerSequenceCommonEnds[i] > 0) {
          scoreAdjustmentsCommonEnds[i] = ((double)commonInstanceNumber-1)/pairwiseAlignmentsPerSequenceCommonEnds[i];
          assert(scoreAdjustmentsCommonEnds[i] >= 1.0);
          assert(scoreAdjustmentsCommonEnds[i] <= commonInstanceNumber-1);
        }
        else {
          scoreAdjustmentsCommonEnds[i] = INT64_MIN;
        }
    }

    //Convert the alignment pairs to an alignment of the caps..
    stSortedSet *sortedAlignment =
                stSortedSet_construct3((int (*)(const void *, const void *))alignedPair_cmpFn,
                (void (*)(void *))alignedPair_destruct);
    while(stList_length(mA->alignedPairs) > 0) {
        stIntTuple *alignedPair = stList_pop(mA->alignedPairs);
        assert(stIntTuple_length(alignedPair) == 5);
        int64_t seqIndex1 = stIntTuple_get(alignedPair, 1);
        int64_t seqIndex2 = stIntTuple_get(alignedPair, 3);
        AdjacencySequence *i = stList_get(sequences, seqIndex1);
        AdjacencySequence *j = stList_get(sequences, seqIndex2);
        assert(i != j);
        int64_t offset1 = stIntTuple_get(alignedPair, 2);
        int64_t offset2 = stIntTuple_get(alignedPair, 4);
        int64_t score = stIntTuple_get(alignedPair, 0);
        if(score <= 0) { //Happens when indel probs are included
            score = 1; //This is the minimum
        }
        assert(score > 0 && score <= PAIR_ALIGNMENT_PROB_1);
        SeqFrag *seqFrag1 = stList_get(seqFrags, seqIndex1);
        SeqFrag *seqFrag2 = stList_get(seqFrags, seqIndex2);
        assert(seqFrag1 != seqFrag2);
        double *scoreAdjustments = seqFrag1->rightEndId == seqFrag2->rightEndId ? scoreAdjustmentsCommonEnds : scoreAdjustmentsNonCommonEnds;
        assert(scoreAdjustments[seqIndex1] != INT64_MIN);
        assert(scoreAdjustments[seqIndex2] != INT64_MIN);
        AlignedPair *alignedPair2 = alignedPair_construct(
                i->subsequenceIdentifier, i->start + (i->strand ? offset1 : -offset1), i->strand,
                j->subsequenceIdentifier, j->start + (j->strand ? offset2 : -offset2), j->strand,
                score*scoreAdjustments[seqIndex1], score*scoreAdjustments[seqIndex2]); //Do the reweighting here.
        assert(stSortedSet_search(sortedAlignment, alignedPair2) == NULL);
        assert(stSortedSet_search(sortedAlignment, alignedPair2->reverse) == NULL);
        stSortedSet_insert(sortedAlignment, alignedPair2);
        stSortedSet_insert(sortedAlignment, alignedPair2->reverse);
        stIntTuple_destruct(alignedPair);
    }
    //Cleanup
    stList_destruct(seqFrags);
    stList_destruct(sequences);
    free(pairwiseAlignmentsPerSequenceNonCommonEnds);
    free(pairwiseAlignmentsPerSequenceCommonEnds);
    free(scoreAdjustmentsNonCommonEnds);
    free(scoreAdjustmentsCommonEnds);
    multipleAlignment_destruct(mA);
    stHash_destruct(endInstanceNumbers);

    return sortedAlignment;
}

void writeEndAlignmentToDisk(End *end, stSortedSet *endAlignment, FILE *fileHandle) {
    fprintf(fileHandle, "%" PRIi64 " %" PRIi64 "\n", end_getName(end), stSortedSet_size(endAlignment));
    stSortedSetIterator *it = stSortedSet_getIterator(endAlignment);
    AlignedPair *aP;
    while((aP = stSortedSet_getNext(it)) != NULL) {
        fprintf(fileHandle, "%" PRIi64 " %" PRIi64 " %i %" PRIi64 " ", aP->subsequenceIdentifier, aP->position, aP->strand, aP->score);
        aP = aP->reverse;
        fprintf(fileHandle, "%" PRIi64 " %" PRIi64 " %i %" PRIi64 "\n", aP->subsequenceIdentifier, aP->position, aP->strand, aP->score);
    }
    stSortedSet_destructIterator(it);
}

stSortedSet *loadEndAlignmentFromDisk(Flower *flower, FILE *fileHandle, End **end) {
    stSortedSet *endAlignment =
                stSortedSet_construct3((int (*)(const void *, const void *))alignedPair_cmpFn,
                (void (*)(void *))alignedPair_destruct);
    char *line = stFile_getLineFromFile(fileHandle);
    if(line == NULL) {
        *end = NULL;
        return NULL;
    }
    Name flowerName;
    int64_t lineNumber;
    int64_t i = sscanf(line, "%" PRIi64 " %" PRIi64 "", &flowerName, &lineNumber);
    if(i != 2 || lineNumber < 0) {
        st_errAbort("We encountered a mis-specified name in loading the first line of an end alignment from the disk: '%s'\n", line);
    }
    *end = flower_getEnd(flower, flowerName);
    if(*end == NULL) {
        st_errAbort("We encountered an end name that is not in the database: '%s'\n", line);
    }
    free(line);    
    for(int64_t i=0; i<lineNumber; i++) {
        line = stFile_getLineFromFile(fileHandle);
        if(line == NULL) {
            st_errAbort("Got a null line when parsing an end alignment\n");
        }
        int64_t sI1, sI2;
        int64_t p1, st1, p2, st2, score1, score2;
        int64_t i = sscanf(line, "%" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 " %" PRIi64 "", &sI1, &p1, &st1, &score1, &sI2, &p2, &st2, &score2);
        (void)i;
        if(i != 8) {
            st_errAbort("We encountered a mis-specified name in loading an end alignment from the disk: '%s'\n", line);
        }
        stSortedSet_insert(endAlignment, alignedPair_construct(sI1, p1, st1, sI2, p2, st2, score1, score2));
        free(line);
    }
    return endAlignment;
}

