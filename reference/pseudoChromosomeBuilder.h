/*
 * pseudoChromosomeBuilder.h
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#ifndef PSEUDOCHROMOSOMEBUILDER_H_
#define PSEUDOCHROMOSOMEBUILDER_H_

#include "cactus.h"

/*
 * Makes a pseudo chromosome for each pair of attached ends in the reference. Assumes
 * that each attached end in the top level problem has only one cap.. as we currently do
 * with setup, if we merge these ends this will have to change.
 */
void makeTopLevelPseudoChromosomes(Flower *flower, Reference *reference);

/*
 * Gets the parent group of the flower and it's reference (which must be defined),
 * then for each set of attached ends in the flower it gathers the corresponding pairing
 * in the higher level reference and constructs a corresponding set of pseudo chromosomes,
 * one for each pair of attached ends.
 */
void makeIntermediateLevelPseudoChromosomes(Flower *flower, Reference *reference);

#endif /* PSEUDOCHROMOSOMEBUILDER_H_ */
