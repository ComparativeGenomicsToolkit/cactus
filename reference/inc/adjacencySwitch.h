/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

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
 * Constructs an AdjacencySwitch from two adjacency pairs.
 * Let an adjacency switch be a pair of adjacencies created
 * so that the four ends of the two initial adjacency pairs are paired into two, new and
 * distinct adjacency pairs. There are two possible adjacency switches for any AdjacencySwitch, this is configured
 * by the switch state.
 */
AdjacencySwitch *adjacencySwitch_construct(AdjacencyPair *adjacencyPair1, AdjacencyPair *adjacencyPair2, bool switchState);

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
 * adjacency pairs. The strength of an AdjacencySwitch
 * is the residual strength.
 */
uint32_t adjacencySwitch_getStrength(AdjacencySwitch *adjacencySwitch);

/*
 * The adjacency switch may create zero, one or two pseudo-adjacencies. Returns this number.
 */
int32_t adjacencySwitch_getNumberOfPseudoAdjacencies(AdjacencySwitch *adjacencySwitch);

/*
 * Compares two adjacency switches. Adjacency switch 1 is greater than adjacency switch 2 iff:
 * if it has fewer pseudo-adjacencies or the same number of pseudo-adjacencies and greather strength.
 * It is equal in strength in the number of pseudo-adjacencies and strength are the same, otherwise
 * it is less than.
 */
int adjacencySwitch_compareStrengthAndPseudoAdjacencies(AdjacencySwitch *adjacencySwitch1, AdjacencySwitch *adjacencySwitch2);

/*
 * Creates two adjacency pairs for the maximum strength adjacency switch, adds them
 * to the adjacency hash and removes the original adjacency pairs in the adjacency switch
 * from the adjacency hash.
 */
void adjacencySwitch_switch(AdjacencySwitch *adjacencySwitch, stHash *adjacencyHash);

/*
 * Returns the strongest adjacency switch between the two components of adjacency pairs.
 */
AdjacencySwitch *adjacencySwitch_getStrongestAdjacencySwitch(stList *component, stList *component2, stHash *adjacencies);

#endif /* ADJACENCYSWITCH_H_ */
