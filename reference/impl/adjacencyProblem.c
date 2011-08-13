#include "cactus.h"
#include "sonLib.h"
#include "cactusMatchingAlgorithms.h"
#include "adjacencyProblem.h"

typedef struct _ReferenceInterval {
    int32_t _5Node;
    int32_t _3Node;
    struct _ReferenceInterval *nReferenceInterval;
} ReferenceInterval;

static ReferenceInterval *referenceInterval_construct(int32_t _5Node,
        int32_t _3Node, ReferenceInterval *nReferenceInterval) {
    ReferenceInterval *referenceInterval = st_malloc(sizeof(ReferenceInterval));
    assert(_5Node != _3Node);
    referenceInterval->_5Node = _5Node;
    referenceInterval->_3Node = _3Node;
    referenceInterval->nReferenceInterval = nReferenceInterval;
    return referenceInterval;
}

typedef struct _ReferenceIntervalInsertion {
    stIntTuple *chain;
    bool orientation;
    double score;
    ReferenceInterval *referenceInterval;
} ReferenceIntervalInsertion;

static ReferenceIntervalInsertion *referenceIntervalInsertion_construct(
        stIntTuple *chain, bool orientation, double score,
        ReferenceInterval *referenceInterval) {
    ReferenceIntervalInsertion *referenceIntervalInsertion = st_malloc(
            sizeof(ReferenceIntervalInsertion));
    referenceIntervalInsertion->chain = chain;
    referenceIntervalInsertion->orientation = orientation;
    referenceIntervalInsertion->score = score;
    referenceIntervalInsertion->referenceInterval = referenceInterval;
    return referenceIntervalInsertion;
}

static stList *getReferenceIntervalInsertions(stList *reference,
        stIntTuple *chain, double *z, int32_t nodeNumber) {
    /*
     * Calculates the scores of the inserting the chain at each possible position in the reference.
     */
    stList *referenceIntervalInsertions = stList_construct3(0, free);

    int32_t _5Node = stIntTuple_getPosition(chain, 0);
    int32_t _3Node = stIntTuple_getPosition(chain, 1);
    for (int32_t i = 0; i < stList_length(reference); i++) {
        ReferenceInterval *referenceInterval = stList_get(reference, i);
        double positiveScore = z[referenceInterval->_3Node * nodeNumber
                + _5Node];
        double negativeScore = z[referenceInterval->_3Node * nodeNumber
                + _3Node];
        while (referenceInterval != NULL) {
            positiveScore += z[referenceInterval->_5Node * nodeNumber + _3Node];
            negativeScore += z[referenceInterval->_5Node * nodeNumber + _5Node];
            referenceInterval = referenceInterval->nReferenceInterval;
        }
        referenceInterval = stList_get(reference, i);
        stList_append(
                referenceIntervalInsertions,
                referenceIntervalInsertion_construct(chain, 1, positiveScore,
                        referenceInterval));
        stList_append(
                referenceIntervalInsertions,
                referenceIntervalInsertion_construct(chain, 0, negativeScore,
                        referenceInterval));
        while (referenceInterval->nReferenceInterval != NULL) {
            referenceInterval = referenceInterval->nReferenceInterval;
            positiveScore = positiveScore + z[referenceInterval->_3Node
                    * nodeNumber + _5Node] - z[referenceInterval->_5Node
                    * nodeNumber + _3Node];
            negativeScore = negativeScore + z[referenceInterval->_3Node
                    * nodeNumber + _3Node] - z[referenceInterval->_5Node
                    * nodeNumber + _5Node];
            stList_append(
                    referenceIntervalInsertions,
                    referenceIntervalInsertion_construct(chain, 1,
                            positiveScore, referenceInterval));
            stList_append(
                    referenceIntervalInsertions,
                    referenceIntervalInsertion_construct(chain, 0,
                            negativeScore, referenceInterval));
        }
    }

    return referenceIntervalInsertions;
}

static void insert(ReferenceIntervalInsertion *referenceIntervalInsertion) {
    /*
     * Inserts a chain into a reference interval.
     */
    ReferenceInterval *newInterval = referenceInterval_construct(
            stIntTuple_getPosition(referenceIntervalInsertion->chain,
                    referenceIntervalInsertion->orientation ? 0 : 1),
            stIntTuple_getPosition(referenceIntervalInsertion->chain,
                    referenceIntervalInsertion->orientation ? 1 : 0),
            referenceIntervalInsertion->referenceInterval->nReferenceInterval);
    referenceIntervalInsertion->referenceInterval->nReferenceInterval
            = newInterval;
}

stList *makeReferenceGreedily(stList *stubs, stList *chainsList, double *z,
        int32_t nodeNumber, double *totalScore) {
    /*
     * Constructs a reference greedily, by picking a best possible insertion
     * at each of the |R| steps.
     */
    stList *reference = stList_construct3(0, free);

    /*
     * Make the stubs.
     */
    stListIterator *it = stList_getIterator(stubs);
    stIntTuple *stub;
    *totalScore = 0.0;
    while ((stub = stList_getNext(it)) != NULL) {
        int32_t _5Node = stIntTuple_getPosition(stub, 0);
        int32_t _3Node = stIntTuple_getPosition(stub, 1);
        *totalScore += z[_5Node * nodeNumber + _3Node];
        stList_append(reference,
                referenceInterval_construct(_5Node, _3Node, NULL));
    }
    stList_destructIterator(it);

    /*
     * Now greedily insert the chains.
     */
    stSortedSet *chains = stList_getSortedSet(chainsList, NULL);
    while (stSortedSet_size(chains) > 0) {
        ReferenceIntervalInsertion maxReferenceIntervalInsertion;
        maxReferenceIntervalInsertion.score = -1;
        stSortedSetIterator *it2 = stSortedSet_getIterator(chains);
        stIntTuple *chain;
        /*
         * Fine the best possible insertion.
         */
        while ((chain = stSortedSet_getNext(it2)) != NULL) {
            stList *referenceIntervalInsertions =
                    getReferenceIntervalInsertions(reference, chain, z,
                            nodeNumber);
            for (int32_t i = 0; i < stList_length(referenceIntervalInsertions); i++) {
                ReferenceIntervalInsertion *referenceIntervalInsertion =
                        stList_get(referenceIntervalInsertions, i);
                if (maxReferenceIntervalInsertion.score
                        < referenceIntervalInsertion->score) {
                    maxReferenceIntervalInsertion.score
                            = referenceIntervalInsertion->score;
                    maxReferenceIntervalInsertion.chain
                            = referenceIntervalInsertion->chain;
                    maxReferenceIntervalInsertion.orientation
                            = referenceIntervalInsertion->orientation;
                    maxReferenceIntervalInsertion.referenceInterval
                            = referenceIntervalInsertion->referenceInterval;
                }
            }
            stList_destruct(referenceIntervalInsertions);
        }
        stSortedSet_destructIterator(it2);

        /*
         * Update the reference
         */
        assert(maxReferenceIntervalInsertion.score != -1);
        stSortedSet_remove(chains, maxReferenceIntervalInsertion.chain);
        insert(&maxReferenceIntervalInsertion);
        *totalScore += maxReferenceIntervalInsertion.score;
    }
    stSortedSet_destruct(chains);

    return reference;
}

static void removeChainFromReference(stIntTuple *chain, stList *reference) {
    /*
     * Remove chain from the reference.
     */
    int32_t _5Node = stIntTuple_getPosition(chain, 0);
    int32_t _3Node = stIntTuple_getPosition(chain, 1);
    assert(_5Node != _3Node);
    for (int32_t i = 0; i < stList_length(reference); i++) {
        ReferenceInterval *referenceInterval = stList_get(reference, i);
        ReferenceInterval *pReferenceInterval = NULL;
        while (referenceInterval != NULL) {
            assert(referenceInterval->_3Node != referenceInterval->_5Node);
            if (referenceInterval->_5Node == _5Node
                    || referenceInterval->_5Node == _3Node) {
                assert(
                        referenceInterval->_3Node == _5Node
                                || referenceInterval->_3Node == _3Node);
                if (pReferenceInterval == NULL) {
                    stList_set(reference, i,
                            referenceInterval->nReferenceInterval);
                } else {
                    pReferenceInterval->nReferenceInterval
                            = referenceInterval->nReferenceInterval;
                }
                free(referenceInterval);
                return;
            } else {
                assert(
                        referenceInterval->_3Node != _5Node
                                && referenceInterval->_3Node != _3Node);
            }
            pReferenceInterval = referenceInterval;
            referenceInterval = referenceInterval->nReferenceInterval;
        }
    }
    assert(0);
}

void gibbsSamplingWithSimulatedAnnealing(stList *reference, stList *chains,
        double *z, int32_t permutations, double(*temperature)(double),
        bool pureGreedy) {
    /*
     * Update a reference by sampling.
     */
    int32_t nodeNumber = stList_length(reference) + stList_length(chains);
    for (int32_t k = 0; k < permutations; k++) {
        chains = stList_copy(chains, NULL);
        stList_shuffle(chains);
        for (int32_t i = 0; i < stList_length(chains); i++) {
            stIntTuple *chain = stList_get(chains, i);
            removeChainFromReference(chain, reference);

            /*
             * Now find a new place to insert it.
             */
            stList *referenceIntervalInsertions =
                    getReferenceIntervalInsertions(reference, chain, z,
                            nodeNumber);

            const int32_t l = stList_length(referenceIntervalInsertions);
            double conditionalSum = 0.0;
            double expConditionalSum = 0.0;
            double *normalScores = st_malloc(sizeof(double) * l);
            for (int32_t j = 0; j < l; j++) {
                ReferenceIntervalInsertion *referenceIntervalInsertion =
                        stList_get(referenceIntervalInsertions, j);
                normalScores[j] = referenceIntervalInsertion->score;
                conditionalSum += normalScores[j];
            }

            if (pureGreedy) {
                double maxScore = -1;
                int32_t insertPoint = -1;
                for (int32_t j = 0; j < l; j++) {
                    if (normalScores[j] > maxScore) {
                        maxScore = normalScores[j];
                        insertPoint = j;
                    }
                }
                assert(maxScore != -1);
                insert(stList_get(referenceIntervalInsertions, insertPoint));
            } else {
                for (int32_t j = 0; j < l; j++) {
                    normalScores[j] /= conditionalSum;
                    if (temperature != NULL) {
                        normalScores[j] = exp(
                                -1 / (temperature((double) k / permutations)
                                        * normalScores[j]));
                    }
                    expConditionalSum += normalScores[j];
                }
                double chosenScore = st_random() * expConditionalSum;
                bool b = 0;
                for (int32_t j = 0; j < l; j++) {
                    chosenScore -= normalScores[j];
                    if (chosenScore < 0) {
                        insert(stList_get(referenceIntervalInsertions, j));
                        b = 1;
                        break;
                    }
                }
                assert(b);
            }
            free(normalScores);
            stList_destruct(referenceIntervalInsertions);
        }
        stList_destruct(chains);
    }
}

stList *convertReferenceToAdjacencyEdges(stList *reference) {
    stList *edges = stList_construct3(0, (void(*)(void *)) stIntTuple_destruct);
    for (int32_t i = 0; i < stList_length(reference); i++) {
        ReferenceInterval *startReferenceInterval = stList_get(reference, i);
        ReferenceInterval *referenceInterval = startReferenceInterval;
        while (referenceInterval->nReferenceInterval != NULL) {
            stList_append(
                    edges,
                    constructEdge(referenceInterval->_3Node,
                            referenceInterval->nReferenceInterval->_5Node));
            referenceInterval = referenceInterval->nReferenceInterval;
        }
        assert(referenceInterval != NULL);
        assert(startReferenceInterval != NULL);
        stList_append(
                edges,
                constructEdge(startReferenceInterval->_5Node,
                        referenceInterval->_3Node));
    }
    return edges;
}

double calculateMaxZ(int32_t nodeNumber, double *z) {
    double maxZ = 0.0;
    for (int32_t i = 0; i < nodeNumber; i++) {
        for (int32_t j = i + 1; j < nodeNumber; j++) {
            assert(z[i * nodeNumber + j] <= z[j * nodeNumber + i] + 0.000001);
            assert(z[j * nodeNumber + i] <= z[i * nodeNumber + j] + 0.000001);
            maxZ += z[i * nodeNumber + j];
        }
    }
    return maxZ;
}

double calculateZScoreOfReference(stList *reference, int32_t nodeNumber, double *zMatrix) {
    double totalScore = 0.0;
    for (int32_t i = 0; i < stList_length(reference); i++) {
        ReferenceInterval *referenceInterval = stList_get(reference, i);
        while(referenceInterval->nReferenceInterval != NULL) {
            totalScore += zMatrix[referenceInterval->_3Node * nodeNumber + referenceInterval->nReferenceInterval->_5Node];
            referenceInterval = referenceInterval->nReferenceInterval;
        }
        ReferenceInterval *referenceInterval2 = stList_get(reference, i);
        totalScore += zMatrix[referenceInterval->_3Node * nodeNumber + referenceInterval2->_5Node];
    }
    return totalScore;
}

void logReference(stList *reference, int32_t nodeNumber, double *zMatrix,
        double totalScore, const char *message) {
    st_logDebug(
            "Reporting reference with %i intervals, %i nodes and %lf total score out of max score of %lf for %s\n",
            stList_length(reference), nodeNumber, totalScore,
            calculateMaxZ(nodeNumber, zMatrix), message);
    for (int32_t i = 0; i < stList_length(reference); i++) {
        ReferenceInterval *referenceInterval = stList_get(reference, i);
        st_logDebug("\tInterval : ");
        while (referenceInterval->nReferenceInterval != NULL) {
            st_logDebug(
                    "(%i %i %lf) ",
                    referenceInterval->_3Node,
                    referenceInterval->nReferenceInterval->_5Node,
                    zMatrix[referenceInterval->_3Node * nodeNumber
                            + referenceInterval->nReferenceInterval->_5Node]);
            referenceInterval = referenceInterval->nReferenceInterval;
        }
        ReferenceInterval *referenceInterval2 = stList_get(reference, i);
        assert(referenceInterval != NULL);
        assert(referenceInterval2 != NULL);
        st_logDebug(
                "(%i %i %lf) \n",
                referenceInterval2->_5Node,
                referenceInterval->_5Node,
                zMatrix[referenceInterval2->_5Node * nodeNumber
                        + referenceInterval->_3Node]);
    }
}

double calculateZScore(int32_t n, int32_t m, int32_t k, double theta) {
    assert(theta <= 1.0);
    assert(theta >= 0.0);
    if (theta == 0.0) {
        return n * m;
    }
    double beta = 1.0 - theta;
    return ((1 - pow(beta, n)) * pow(beta, k) * (1 - pow(beta, m))) / (theta
            * theta);
}

double exponentiallyDecreasingTemperatureFn(double d) {
    return 1000 * pow(100000, -d);
}

double constantTemperatureFn(double d) {
    return 1;
}

