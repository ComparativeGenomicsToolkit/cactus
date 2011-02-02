/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * terminalFlowers.h
 *
 *  Created on: 25 Aug 2010
 *      Author: benedictpaten
 */

#ifndef TERMINALFLOWERS_H_
#define TERMINALFLOWERS_H_

#include "cactus.h"
#include "sonLib.h"

/*
 * A flower A is a child of a flower B if it is connected to B by a group.
 * A flower A which is a child of flower B is either a 'tangle child' of B or
 * a 'link child' of B, dependent on if it is respectively connected by a link or tangle group.
 * A flower A is a descendant of a flower B if there is a path of successive
 * ancestors from B to A, such that each successive ancestor is the child of the previous ancestor.
 * A flower A which descends from some flower B is
 * a 'tangle descendant' of B if and only if each successive child on the path from B to
 * A is a tangle child of its parent.
 * For a flower A, the set of 'pressed flowers' for A is the set of all terminal tangle descendants
 * of A, or if A is terminal, A itself.
 * The set of pressed flowers can be viewed as a collapse of the multi-scale dimension of the cactus graph.
 */
stList *getListOfPressedFlowers(Flower *flower);

#endif /* TERMINALFLOWERS_H_ */
