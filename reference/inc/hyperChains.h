/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * hyperChains.h
 *
 *  Created on: 25 Aug 2010
 *      Author: benedictpaten
 */

#ifndef HYPERCHAINS_H_
#define HYPERCHAINS_H_

#include "sonLib.h"
#include "cactus.h"

/*
 * See getListOfPressedFlowers for an explanation of a 'pressed flower'.
 * For flower A this method constructs a hash
 * of 'hyperlinks' between the ends in its set of pressed flowers.
 * A hyperlink between two ends a and b exists if there exists a chain in the set of
 * pressed flowers that has a one end and b on the other end.
 */
stHash *constructHyperChains(Flower *flower);

#endif /* HYPERCHAINS_H_ */
