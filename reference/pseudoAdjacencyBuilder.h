/*
 * pseudoAdjacencyBuilder.h
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#ifndef PSEUDOADJACENCYBUILDER_H_
#define PSEUDOADJACENCYBUILDER_H_

#include "cactus.h"

/*
 * Makes the pseudo adjacencies for each pseudo chromosome in reference.
 * In a top-down greedy definition of a reference this is the only function that can
 * be varied.
 */
void makePseudoAdjacencies(Flower *flower, Reference *reference);


#endif /* PSEUDOADJACENCYBUILDER_H_ */
