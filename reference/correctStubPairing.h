/*
 * correctStubPairing.h
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#ifndef CORRECTSTUBPAIRING_H_
#define CORRECTSTUBPAIRING_H_

#include "sonLib.h"
#include "cactus.h"

/*
 * If a-b and c-d are two components respectively containing stub ends a and b and c and d and the correct pairing
 * is a-c b-d then two adjacencies are broken and two are created to make this happen.
 * Destructs the list of components as we go.
 */
void correctAttachedStubEndPairing(stHash *adjacencies, Flower *flower,
		Reference *reference);

#endif /* CORRECTSTUBPAIRING_H_ */
