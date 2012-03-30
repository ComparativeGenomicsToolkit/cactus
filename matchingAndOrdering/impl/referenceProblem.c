#include <stdlib.h>
#include <math.h>
#include "sonLib.h"

#include "stMatchingAlgorithms.h"
#include "stReferenceProblem.h"

typedef struct _ReferenceInterval {
    int32_t _5Node;
    int32_t _3Node;
    struct _ReferenceInterval *nReferenceInterval;
    struct _ReferenceInterval *pReferenceInterval;
} ReferenceInterval;

static ReferenceInterval *referenceInterval_construct(int32_t _5Node,
        int32_t _3Node, ReferenceInterval *nReferenceInterval,
        ReferenceInterval *pReferenceInterval) {
    ReferenceInterval *referenceInterval = st_malloc(sizeof(ReferenceInterval));
    assert(_5Node != _3Node);
    referenceInterval->_5Node = _5Node;
    referenceInterval->_3Node = _3Node;
    referenceInterval->nReferenceInterval = nReferenceInterval;
    referenceInterval->pReferenceInterval = pReferenceInterval;
    return referenceInterval;
}

static void referenceInterval_destruct(ReferenceInterval *referenceInterval) {
    if (referenceInterval->nReferenceInterval != NULL) {
        referenceInterval_destruct(referenceInterval->nReferenceInterval);
    }
    free(referenceInterval);
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

static double referenceIntervalInsertion_isConnectedLeft(
        ReferenceIntervalInsertion *referenceIntervalInsertion, double *z,
        int32_t nodeNumber) {
    int32_t _5Node = stIntTuple_getPosition(referenceIntervalInsertion->chain,
            referenceIntervalInsertion->orientation ? 0 : 1);
    assert(referenceIntervalInsertion->referenceInterval != NULL);
    assert(z[referenceIntervalInsertion->referenceInterval->_3Node
                    * nodeNumber + _5Node] >= 0.0);
    return z[referenceIntervalInsertion->referenceInterval->_3Node * nodeNumber
            + _5Node] > 0.0 ? 1 : 0;
}

static int32_t referenceInterval_get5Node(ReferenceInterval *referenceInterval) {
    if (referenceInterval->nReferenceInterval != NULL) {
        return referenceInterval->nReferenceInterval->_5Node;
    }
    while (referenceInterval->pReferenceInterval != NULL) {
        referenceInterval = referenceInterval->pReferenceInterval;
    }
    assert(referenceInterval != NULL);
    assert(referenceInterval->pReferenceInterval == NULL);
    return referenceInterval->_5Node;
}

static double referenceIntervalInsertion_isConnectedRight(
        ReferenceIntervalInsertion *referenceIntervalInsertion, double *z,
        int32_t nodeNumber) {
    int32_t _3Node = stIntTuple_getPosition(referenceIntervalInsertion->chain,
            referenceIntervalInsertion->orientation ? 1 : 0);
    int32_t _5Node = referenceInterval_get5Node(
            referenceIntervalInsertion->referenceInterval);
    assert(z[_5Node * nodeNumber + _3Node] >= 0.0);
    return z[_5Node * nodeNumber + _3Node] > 0.0 ? 1 : 0;
}

static double referenceIntervalInsertion_isCurrentlyConnected(
        ReferenceIntervalInsertion *referenceIntervalInsertion, double *z,
        int32_t nodeNumber) {
    assert(referenceIntervalInsertion->referenceInterval != NULL);
    int32_t _3Node = referenceIntervalInsertion->referenceInterval->_3Node;
    int32_t _5Node = referenceInterval_get5Node(
            referenceIntervalInsertion->referenceInterval);
    assert(z[_5Node * nodeNumber + _3Node] >= 0.0);
    return z[_5Node * nodeNumber + _3Node] > 0.0 ? 1 : 0;
}

static double referenceIntervalInsertion_connectionScore(
        ReferenceIntervalInsertion *referenceIntervalInsertion, double *z,
        int32_t nodeNumber) {
    return referenceIntervalInsertion_isConnectedLeft(
            referenceIntervalInsertion, z, nodeNumber)
            + referenceIntervalInsertion_isConnectedRight(
                    referenceIntervalInsertion, z, nodeNumber)
            - referenceIntervalInsertion_isCurrentlyConnected(
                    referenceIntervalInsertion, z, nodeNumber);
}

static stList *getReferenceIntervalInsertions(stList *reference,
        stIntTuple *chain, double *z, int32_t nodeNumber) {
    /*
     * Calculates the scores of the inserting the chain at each possible position in the reference.
     *
     * Now modified to avoid numerical bug derivied from doing running sum which involved subtraction.
     */
    stList *referenceIntervalInsertions = stList_construct3(0, free);

    int32_t _5Node = stIntTuple_getPosition(chain, 0);
    int32_t _3Node = stIntTuple_getPosition(chain, 1);
    for (int32_t i = 0; i < stList_length(reference); i++) {
        ReferenceInterval *referenceInterval = stList_get(reference, i);
        stList *intervals = stList_construct();
        while (referenceInterval != NULL) {
            stList_append(intervals, referenceInterval);
            referenceInterval = referenceInterval->nReferenceInterval;
        }
        double *positiveScores = st_malloc(
                stList_length(intervals) * sizeof(double));
        double *negativeScores = st_malloc(
                stList_length(intervals) * sizeof(double));
        referenceInterval = stList_get(reference, i);
        int32_t j = stList_length(intervals) - 1;
        assert(j >= 0);
        positiveScores[j] = z[referenceInterval->_5Node * nodeNumber + _3Node]; //Add the score of the 3 prime end of the interval to the 3 prime most insertion point
        negativeScores[j] = z[referenceInterval->_5Node * nodeNumber + _5Node];
        for (; j > 0; j--) {
            referenceInterval = stList_get(intervals, j);

            assert(z[referenceInterval->_5Node * nodeNumber + _3Node] >= 0);
            positiveScores[j - 1] = z[referenceInterval->_5Node * nodeNumber
                    + _3Node] + positiveScores[j];

            assert(z[referenceInterval->_5Node * nodeNumber + _5Node] >= 0);
            negativeScores[j - 1] = z[referenceInterval->_5Node * nodeNumber
                    + _5Node] + negativeScores[j];
        }
        double positiveScore = 0.0;
        double negativeScore = 0.0;
        for (j = 0; j < stList_length(intervals); j++) {
            referenceInterval = stList_get(intervals, j);

            assert(z[referenceInterval->_3Node * nodeNumber + _5Node] >= 0.0);
            assert(positiveScore + 0.0 == positiveScore);
            positiveScore += z[referenceInterval->_3Node * nodeNumber + _5Node];
            assert(positiveScore >= 0.0);

            assert(z[referenceInterval->_3Node * nodeNumber + _3Node] >= 0.0);
            assert(negativeScore + 0.0 == negativeScore);
            negativeScore += z[referenceInterval->_3Node * nodeNumber + _3Node];
            assert(negativeScore >= 0.0);

            positiveScores[j] += positiveScore;
            negativeScores[j] += negativeScore;

            stList_append(
                    referenceIntervalInsertions,
                    referenceIntervalInsertion_construct(chain, 1,
                            positiveScores[j], referenceInterval));
            stList_append(
                    referenceIntervalInsertions,
                    referenceIntervalInsertion_construct(chain, 0,
                            negativeScores[j], referenceInterval));
        }
        free(positiveScores);
        free(negativeScores);
        stList_destruct(intervals);
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
            referenceIntervalInsertion->referenceInterval->nReferenceInterval,
            referenceIntervalInsertion->referenceInterval);
    assert(referenceIntervalInsertion->referenceInterval != NULL);
    if (referenceIntervalInsertion->referenceInterval->nReferenceInterval
            != NULL) {
        referenceIntervalInsertion->referenceInterval->nReferenceInterval->pReferenceInterval
                = newInterval;
    }
    referenceIntervalInsertion->referenceInterval->nReferenceInterval
            = newInterval;
}

stList *makeReferenceGreedily(stList *stubs, stList *chainsList, double *z,
        int32_t nodeNumber, double *totalScore, bool fast) {
    /*
     * Constructs a reference greedily, by picking a best possible insertion
     * at each of the |R| steps.
     *
     * If fast=1 then it picks the best possible insertion for a single randomly selected
     * chain at each step.
     *
     */
    stList *reference = stList_construct3(0,
            (void(*)(void *)) referenceInterval_destruct);

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
                referenceInterval_construct(_5Node, _3Node, NULL, NULL));
    }
    stList_destructIterator(it);
    //assert(stList_length(stubs) > 0);
    //assert(stList_length(reference) > 0);

    /*
     * Now greedily insert the chains.
     */
    stSortedSet *chains = stList_getSortedSet(chainsList, NULL);
    while (stSortedSet_size(chains) > 0) {
        ReferenceIntervalInsertion maxReferenceIntervalInsertion;
        maxReferenceIntervalInsertion.score = -1;
        maxReferenceIntervalInsertion.chain = NULL;
        maxReferenceIntervalInsertion.orientation = 0;
        maxReferenceIntervalInsertion.referenceInterval = NULL;
        stSortedSetIterator *it2 = stSortedSet_getIterator(chains);
        assert(stList_length(chainsList) > 0);
        stIntTuple *chain;
        /*
         * Find the best possible insertion.
         */
        while ((chain = stSortedSet_getNext(it2)) != NULL) {
            stList *referenceIntervalInsertions =
                    getReferenceIntervalInsertions(reference, chain, z,
                            nodeNumber);
            assert(stList_length(referenceIntervalInsertions) > 0);
            for (int32_t i = 0; i < stList_length(referenceIntervalInsertions); i++) {
                ReferenceIntervalInsertion *referenceIntervalInsertion =
                        stList_get(referenceIntervalInsertions, i);
                if (maxReferenceIntervalInsertion.chain == NULL
                        || maxReferenceIntervalInsertion.score
                                < referenceIntervalInsertion->score
                        || (maxReferenceIntervalInsertion.score
                                == referenceIntervalInsertion->score
                                && referenceIntervalInsertion_connectionScore(
                                        referenceIntervalInsertion, z,
                                        nodeNumber)
                                        > referenceIntervalInsertion_connectionScore(
                                                &maxReferenceIntervalInsertion,
                                                z, nodeNumber))) {
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
            if(fast) {
                break;
            }
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
                assert(pReferenceInterval != NULL);
                pReferenceInterval->nReferenceInterval
                        = referenceInterval->nReferenceInterval;
                if (referenceInterval->nReferenceInterval) {
                    referenceInterval->nReferenceInterval->pReferenceInterval
                            = pReferenceInterval;
                }
                free(referenceInterval);
                return;
            } else {
                assert(
                        referenceInterval->_3Node != _3Node
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
    int32_t nodeNumber = (stList_length(reference) + stList_length(chains)) * 2;
    chains = stList_copy(chains, NULL);
    for (int32_t k = 0; k < permutations; k++) {
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
            double *connectionScores = st_malloc(sizeof(double) * l);
            for (int32_t j = 0; j < l; j++) {
                ReferenceIntervalInsertion *referenceIntervalInsertion =
                        stList_get(referenceIntervalInsertions, j);
                normalScores[j] = referenceIntervalInsertion->score;
                connectionScores[j]
                        = referenceIntervalInsertion_connectionScore(
                                referenceIntervalInsertion, z, nodeNumber);
                assert(normalScores[j] >= 0.0);
                conditionalSum += normalScores[j];
            }
            if (pureGreedy) {
                double maxScore = INT32_MIN;
                double maxConnectionScore = INT32_MIN;
                int32_t insertPoint = -1;
                for (int32_t j = 0; j < l; j++) {
                    if (normalScores[j] > maxScore || (normalScores[j]
                            == maxScore && connectionScores[j]
                            > maxConnectionScore)) {
                        maxScore = normalScores[j];
                        maxConnectionScore = connectionScores[j];
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
                for (int32_t j = 0; j < l; j++) {
                    chosenScore -= normalScores[j];
                    if (chosenScore < 0 || j == l - 1) {
                        insert(stList_get(referenceIntervalInsertions, j));
                        break;
                    }
                }
            }
            free(normalScores);
            free(connectionScores);
            stList_destruct(referenceIntervalInsertions);
        }
    }
    stList_destruct(chains);
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

double calculateZScoreOfReference(stList *reference, int32_t nodeNumber,
        double *zMatrix) {
    double totalScore = 0.0;
    for (int32_t i = 0; i < stList_length(reference); i++) {
        ReferenceInterval *referenceInterval = stList_get(reference, i);
        int32_t j = referenceInterval->_5Node;
        while (referenceInterval != NULL) {
            totalScore += zMatrix[referenceInterval->_3Node * nodeNumber + j];
            ReferenceInterval *referenceInterval2 =
                    referenceInterval->nReferenceInterval;
            while (referenceInterval2 != NULL) {
                totalScore += zMatrix[referenceInterval->_3Node * nodeNumber
                        + referenceInterval2->_5Node];
                referenceInterval2 = referenceInterval2->nReferenceInterval;
            }
            referenceInterval = referenceInterval->nReferenceInterval;
        }
    }
    return totalScore;
}

void logZScoreOfReference(stList *reference, int32_t nodeNumber,
        double *zMatrix) {
    for (int32_t i = 0; i < stList_length(reference); i++) {
        ReferenceInterval *referenceInterval = stList_get(reference, i);
        int32_t j = referenceInterval->_5Node;
        while (referenceInterval != NULL) {
            st_logDebug("The score of the adjacency %i %i %lf\n",
                    referenceInterval->_3Node, j,
                    zMatrix[referenceInterval->_3Node * nodeNumber + j]);
            ReferenceInterval *referenceInterval2 =
                    referenceInterval->nReferenceInterval;
            while (referenceInterval2 != NULL) {
                st_logDebug(
                        "The score of the adjacency %i %i %lf\n",
                        referenceInterval->_3Node,
                        referenceInterval2->_5Node,
                        zMatrix[referenceInterval->_3Node * nodeNumber
                                + referenceInterval2->_5Node]);
                referenceInterval2 = referenceInterval2->nReferenceInterval;
            }
            referenceInterval = referenceInterval->nReferenceInterval;
        }
    }
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
                referenceInterval->_3Node,
                referenceInterval2->_5Node,
                zMatrix[referenceInterval2->_5Node * nodeNumber
                        + referenceInterval->_3Node]);
    }
}

double calculateZScore(int32_t n, int32_t m, int32_t k, double theta) {
    assert(theta <= 1.0);
    assert(theta >= 0.0);
    if (theta == 0.0) {
        return ((double) n) * m;
    }
    double beta = 1.0 - theta;
    return ((1.0 - pow(beta, n)) / theta) * pow(beta, k)
            * ((1.0 - pow(beta, m)) / theta);
}

double exponentiallyDecreasingTemperatureFn(double d) {
    return 1000 * pow(100000, -d);
}

double constantTemperatureFn(double d) {
    return 1;
}

