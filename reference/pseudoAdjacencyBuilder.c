/*
 * pseudoAdjacencyBuilder.c
 *
 *  Created on: 1 Apr 2010
 *      Author: benedictpaten
 *
 * This source file contains the functions which choose the pseudo adjacencies
 * in a pseudo chromosome. This is the vital bit!
 */

#include "cactus.h"
#include "sonLib.h"
#include "adjacencyPairs.h"
#include "adjacencyPairsHash.h"
#include "pickAdjacencyPairs.h"
#include "mergeCycles.h"
#include "correctStubPairing.h"

static void fillInPseudoAdjacenciesP(End *end1,
        PseudoChromosome *pseudoChromosome, stHash *adjacencies) {
    /*
     * Traverses the connected component of the contig constructing the pseudo-adjacencies in the pseudo-chromosome.
     */
    AdjacencyPair *adjacencyPair = stHash_search(adjacencies, end1);
    assert(adjacencyPair != NULL);
    //construct the new pseudo-adjacency!
    End *end2 = adjacencyPair_getOtherEnd(adjacencyPair, end1);
    assert(end_getGroup(end1) == end_getGroup(end2));
    pseudoAdjacency_construct(end1, end2, pseudoChromosome);
    if (pseudoChromosome_get3End(pseudoChromosome) != end2) {
        assert(end_isBlockEnd(end2));
        fillInPseudoAdjacenciesP(end_getOtherBlockEnd(end2), pseudoChromosome,
                adjacencies);
    }
}

static void fillInPseudoAdjacencies(Flower *flower, Reference *reference,
        stHash *adjacencies) {
    /*
     * Walks through the list of pseudo-chromosomes, then
     * iterates through the adjacencies populating each psuedo chromosome with adjacencies.
     */
    PseudoChromosome *pseudoChromosome;
    Reference_PseudoChromosomeIterator *pseudoChromosomeIterator =
            reference_getPseudoChromosomeIterator(reference);
    while ((pseudoChromosome = reference_getNextPseudoChromosome(
            pseudoChromosomeIterator)) != NULL) {
        //Now we fill in the pseudo-adjacencies
        End *end = pseudoChromosome_get5End(pseudoChromosome);
        fillInPseudoAdjacenciesP(end, pseudoChromosome, adjacencies);
    }
    reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
}

void makePseudoAdjacencies(Flower *flower, Reference *reference) {
    /*
     * Makes the pseudo adjacencies for each pseudo chromosome in reference.
     */
	stHash *adjacenciesHash = pickAdjacencyPairs(flower, reference);
    mergeCycles(adjacenciesHash, flower);
    correctAttachedStubEndPairing(adjacenciesHash, flower, reference);
    fillInPseudoAdjacencies(flower, reference, adjacenciesHash);
    adjacencyHash_cleanUp(adjacenciesHash, flower);
}
