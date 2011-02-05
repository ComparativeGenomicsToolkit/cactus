/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * referenceViewer.h
 *
 *  Created on: 27-Apr-2010
 *      Author: benedictpaten
 */

#ifndef REFERENCE_VIEWER_H_
#define REFERENCE_VIEWER_H_

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "commonC.h"
#include "referenceViewer.h"

///////////////////////
// Plotter interface.
//////////////////////

/*
 * Adds an end to the plot
 */
void addEnd(End *end, void *extraArgument);

/*
 * Adds a block edge to the plot.
 */
void addBlockEdge(End *_5End, End *_3End, void *extraArgument);

/*
 * Adds a pseudo-adjacency edge to the plot.
 */
void addPseudoAdjacencyEdge(End *_5End, End *_3End, void *extraArgument);

/*
 * Adds an adjacency edge to the plot.
 */
void addAdjacencyEdge(End *_5End, End *_3End, void *extraArgument);

/*
 * Adds a telomere edge.
 */
void addTelomereEdge(End *_5End, End *_3End, void *extraArgument);

///////////////////////
// End plotter interface.
//////////////////////

/*
 * Function which does the magic and calls the above functions,
 * in the order given (first ends, second block edges... etc. , last telomere edges).
 */
void makeReferenceGraph(Reference *reference, void *extraArgument);

#endif
