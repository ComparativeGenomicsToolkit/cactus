/*
 * adjacencySwitch.h
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#ifndef ADJACENCYSWITCH_H_
#define ADJACENCYSWITCH_H_

#include "cactus.h"
#include "adjacencyPairs.h"
#include "sonLib.h"

typedef struct _adjacencySwitch AdjacencySwitch;

/*
 * Constructs an adjacencySwitch from two adjacency pairs.
 * Let an adjacency switch be a pair of adjacencies created
 * so that the four ends of the two initial adjacency pairs are paired into two, new and
 * distinct adjacency pairs. There are two possible adjacency switches for any AdjacencySwitch.
 */
AdjacencySwitch *adjacencySwitch_construct(AdjacencyPair *adjacencyPair1, AdjacencyPair *adjacencyPair2);

/*
 * Destroys the AdjacencySwitch, but not the two contained adjacency pairs.
 */
void adjacencySwitch_destruct(AdjacencySwitch *adjacencySwitch);

/*
 * Gets the first existing adjacency pair in the AdjacencySwitch.
 */
AdjacencyPair *adjacencySwitch_getAdjacencyPair1(AdjacencySwitch *adjacencySwitch);

/*
 * Gets the second existing adjacency pair in the AdjacencySwitch.
 */
AdjacencyPair *adjacencySwitch_getAdjacencyPair2(AdjacencySwitch *adjacencySwitch);

/*
 * Let the strength of an adjacency pair be given by the function adjacencyPair_strength().
 * Let the initial strength be the sum of the strengths of the two adjacency pairs in the
 * AdjacencySwitch.
 * Let the residual strength of an adjacency switch be the combined strength of the new
 * adjacency pairs. Let the maximum resdidual strength of an AdjacencySwitch be the maximum
 * over the two possible
 * adjacency switches of the AdjacencySwitch. The strength of an AdjacencySwitch
 * is the maximum residual strength minus its initial strength.
 */
double adjacencySwitch_getStrength(AdjacencySwitch *adjacencySwitch);

/*
 * Comparison function that compares the AdjacencySwitches by the the value of the
 * adjacencySwitch_getStrength
 * function.
 */
int32_t adjacencySwitch_cmpByStrength(AdjacencySwitch *adjacencySwitch, AdjacencySwitch *adjacencySwitch2);

/*
 * Creates two adjacency pairs for the maximum strength adjacency switch, adds them
 * to the adjacency hash and removes the original adjacency pairs in the adjacency switch
 * from the adjacency hash.
 */
void adjacencySwitch_switch(AdjacencySwitch *adjacencySwitch, stHash *adjacencyHash);

/*
 * Computes the strongest AdjacencySwitch between two components of adjacency pairs.
 * If no valid adjacency switch is found it returns NULL.
 */
AdjacencySwitch *adjacencySwitch_getStrongestAdjacencySwitch(stList *component, stList *component2);

#endif /* ADJACENCYSWITCH_H_ */
