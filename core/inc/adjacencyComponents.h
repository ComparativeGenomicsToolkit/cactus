/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * adjacencyComponents.h
 *
 *  Created on: 10 Sep 2010
 *      Author: benedictpaten
 */

#ifndef ADJACENCYCOMPONENTS_H_
#define ADJACENCYCOMPONENTS_H_

#include "sonLib.h"
#include "pinchGraph.h"

/*
 * Returns non-zero if edge is a degree one edge that is not a stub.
 */
bool passThroughDegree1EdgesFn(struct PinchEdge *edge);

/*
 * Returns zero for all pinch edges.
 */
bool doNotPassThroughDegree1EdgesFn(struct PinchEdge *edge);

/*
 * Gets a list of adjacency components from the pinch graph. Adjacency components are components of adjacency 1-degree
 * edges. Each adjacency component is represented as a sorted set of pinch-vertices.
 */
stList *getAdjacencyComponents(struct PinchGraph *pinchGraph);

/*
 * Like getAdjacencyComponents. A pinch edge is counted as an adjacency edge if passThroughEdge returns non-zero.
 */
stList *getAdjacencyComponents2(struct PinchGraph *pinchGraph, bool (*passThroughEdge)(struct PinchEdge *));

/*
 * Creates a hash of nodes to adjacency components, each node being a member of one adjacency component.
 */
stHash *getVertexToAdjacencyComponentHash(struct PinchGraph *pinchGraph, stList *adjacencyComponents);

/*
 * Constructs a graph in which the adjacency components are nodes and the
 * edges are blocks. The graph is represented as a set of adjacency lists, one for adjacency component.
 */
stList *getAdjacencyComponentGraph(struct PinchGraph *pinchGraph, stList *adjacencyComponents, stHash *vertexToAdjacencyComponentsHash);


#endif /* ADJACENCYCOMPONENTS_H_ */
