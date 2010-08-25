/*
 * correctSubPairing.c
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#include "sonLib.h"
#include "cactus.h"
#include "adjacencyPairs.h"
#include "adjacencyPairsHash.h"
#include "adjacencyPairComponents.h"
#include "adjacencySwitch.h"

void stubsToComponents_addComponent(stHash *stubsToComponents,
        stList *component) {
    End *end1, *end2;
    getAttachedStubEndsInComponent(component, &end1, &end2);
    assert(end1 != NULL);
    assert(end2 != NULL);
    stHash_insert(stubsToComponents, end1, component);
    stHash_insert(stubsToComponents, end2, component);
}

void stubsToComponents_removeComponent(stHash *stubsToComponents,
        stList *component) {
    End *end1, *end2;
    getAttachedStubEndsInComponent(component, &end1, &end2);
    assert(end1 != NULL);
    assert(end2 != NULL);
    stHash_remove(stubsToComponents, end1);
    stHash_remove(stubsToComponents, end2);
}

void correctAttachedStubEndPairing(stHash *adjacencies, Flower *flower,
        Reference *reference) {
    stList *components = adjacencyHash_getConnectedComponents(adjacencies,
            flower);
    //Put into a hash of stub ends to components..
    stHash *stubsToComponents = stHash_construct();
    for (int32_t i = 0; i < stList_length(components); i++) { //Get mispaired contigs, breaking them up.
        stubsToComponents_addComponent(stubsToComponents, stList_get(
                components, i));
    }
    PseudoChromosome *pseudoChromosome;
    Reference_PseudoChromosomeIterator *pseudoChromosomeIterator =
            reference_getPseudoChromosomeIterator(reference);
    while ((pseudoChromosome = reference_getNextPseudoChromosome(
            pseudoChromosomeIterator)) != NULL) {
        stList *component1 = stHash_search(stubsToComponents,
                pseudoChromosome_get5End(pseudoChromosome));
        stList *component2 = stHash_search(stubsToComponents,
                pseudoChromosome_get3End(pseudoChromosome));
        if (component1 != component2) {
            AdjacencySwitch *adjacencySwitch =
                    adjacencySwitch_getStrongestAdjacencySwitch(component1,
                            component2);
            assert(adjacencySwitch != NULL);
            //Make the swtich
            adjacencySwitch_switch(adjacencySwitch, adjacencies);

            //Get the the new components..
            stList *newComponent1 =
                    adjacencyHash_getConnectedComponent(adjacencies, flower,
                            adjacencyPair_getEnd1(
                                    adjacencySwitch_getAdjacencyPair1(
                                            adjacencySwitch)));
            stList *newComponent2 =
                    adjacencyHash_getConnectedComponent(adjacencies, flower,
                            adjacencyPair_getEnd1(
                                    adjacencySwitch_getAdjacencyPair2(
                                            adjacencySwitch)));
            assert(stList_length(component1) + stList_length(component2) ==
                    stList_length(newComponent1) + stList_length(newComponent2));

            //Update the list of components
            stList_removeItem(components, component1);
            stList_removeItem(components, component2);
            stList_append(components, newComponent1);
            stList_append(components, newComponent2);

            //Update the stubs to components hash
            stubsToComponents_removeComponent(stubsToComponents, component1);
            stubsToComponents_removeComponent(stubsToComponents, component2);
            stubsToComponents_addComponent(stubsToComponents, newComponent1);
            stubsToComponents_addComponent(stubsToComponents, newComponent2);

            //Destroy the old components.
            stList_destruct(component1);
            stList_destruct(component2);

            //Destroy the switch and adjacency pairs
            adjacencyPair_destruct(adjacencySwitch_getAdjacencyPair1(
                    adjacencySwitch));
            adjacencyPair_destruct(adjacencySwitch_getAdjacencyPair2(
                    adjacencySwitch));
            adjacencySwitch_destruct(adjacencySwitch);
        }
    }
    reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
    stHash_destruct(stubsToComponents);
    stList_destruct(components);
}
