/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "referenceViewer.h"

/*
 * Functions to make plot in graph-viz friendly .dot format
 */

static void addEdge(End *end1, End *end2, const char *colour, const char *label, double length, double weight, const char *direction, void *extraArgument) {
    char *nameString1 = cactusMisc_nameToString(end_getName(end1));
    char *nameString2 = cactusMisc_nameToString(end_getName(end2));
    graphViz_addEdgeToGraph(nameString1, nameString2, (FILE *)extraArgument, 1 ? label : "", colour, length, weight, direction);
    free(nameString1);
    free(nameString2);
}

void addEnd(End *end, void *extraArgument) {
	const char *nameString = cactusMisc_nameToStringStatic(end_getName(end));
	const char *colour = "red";
	const char *shape = "circle";
	if(end_isBlockEnd(end)) {
		colour = "black";
		shape = "circle";
	}
	else {
		assert(end_isStubEnd(end));
		assert(end_isAttached(end));
	}
	graphViz_addNodeToGraph(nameString, (FILE *)extraArgument, nameString, 0.5, 0.5, shape, colour, 14);
}

void addBlockEdge(End *_5End, End *_3End, void *extraArgument) {
	addEdge(_5End, _3End, "black", NULL, 0.1, 100, "forward", extraArgument);
}

void addPseudoAdjacencyEdge(End *_5End, End *_3End, void *extraArgument) {
	addEdge(_5End, _3End, "green", NULL, 0.1, 100, "forward", extraArgument);
}

void addAdjacencyEdge(End *_5End, End *_3End, void *extraArgument) {
	const char *colour = "red";
	addEdge(_5End, _3End, colour, NULL, 1, 0.0, "forward", extraArgument);
}

void addTelomereEdge(End *_5End, End *_3End, void *extraArgument) {
	addEdge(_5End, _3End, "magenta", NULL, 0.1, 100, "forward", extraArgument);
}
