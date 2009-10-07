#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "commonC.h"
#include "fastCMaths.h"
#include "bioioC.h"
#include "hashTableC.h"
#include "pinchGraph.h"
#include "cactus.h"
#include "pairwiseAlignment.h"

/*
 * Stuff to check the pinch graph.
 */


////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Test methods
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void checkPinchGraphDegree(struct PinchGraph *graph, int32_t maxDegree) {
	/*
	 * Checks no black edges have degree higher than given max degree (except those of stubs and caps which
	 * can not be undone).
	 */
#ifdef BEN_DEBUG
	int32_t i;
	struct PinchVertex *vertex;

	for(i=0; i<graph->vertices->length; i++) {
		vertex = graph->vertices->list[i];
		if(lengthBlackEdges(vertex) > 0) {
			if(isAStubOrCap(getFirstBlackEdge(vertex)) == FALSE) { //we don't check stubs/caps/
				assert(lengthBlackEdges(vertex) <= maxDegree);
			}
		}
	}
#endif
}


void checkPinchGraph(struct PinchGraph *graph) {
	/*
	 * Checks that the pinch graph is correctly formed. If not throws assertion
	 * errors.
	 */
#ifdef BEN_DEBUG

	int32_t i;
	struct PinchVertex *vertex;
	struct PinchVertex *vertex2;
	struct PinchEdge *edge;
	struct PinchEdge *edge2;
	void *greyEdgeIterator;
	void *blackEdgeIterator;

	for(i=0; i<graph->vertices->length; i++) {
		vertex = graph->vertices->list[i];

		//check the black edges first.
		if(vertex->vertexID == 0) {
			//special condition for sink vertex.
			assert(vertex->vertexID == i);
			assert(lengthBlackEdges(vertex) == 0);
		}
		else {
			assert(lengthBlackEdges(vertex) > 0);

			//check the edges are all okay
			blackEdgeIterator = getBlackEdgeIterator(vertex);
			while((edge = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
				//check not self edge
				assert(edge->to != edge->from);

				//check the vertices they point at are okay.
				assert(edge->from->vertexID == i);
				assert(edge->to->vertexID >= 0);
				assert(edge->to->vertexID < graph->vertices->length);

				//check reverse edge points at the same things.
				assert(edge->rEdge->to == edge->from);
				assert(edge->rEdge->from == edge->to);
				//check edge okay coordiates.
				assert(edge->segment->start <= edge->segment->end);
				//check reverse segment
				assert(edge->segment->rSegment->start == -edge->segment->end);
				assert(edge->segment->rSegment->end == -edge->segment->start);
				//check reverse edge segment is reverse segment
				assert(edge->rEdge->segment == edge->segment->rSegment);

				//check the edge is in the edges tree.
				assert(edge == getContainingBlackEdge(graph, edge->segment->contig, edge->segment->start));
				assert(edge == getContainingBlackEdge(graph, edge->segment->contig, edge->segment->end));

				//check its stub/cap-eyness.
				if(isAStubOrCap(edge)) { //is a cap and must have one end pointing at the source or free
					assert(lengthGreyEdges(edge->from) == 0 ||
						(lengthGreyEdges(edge->from) == 1 &&
						 getFirstGreyEdge(edge->from)->vertexID == 0) ||
						lengthGreyEdges(edge->to) == 0 ||
						(lengthGreyEdges(edge->to) == 1 &&
						 getFirstGreyEdge(edge->to)->vertexID == 0));
				}
				else {
					assert(lengthGreyEdges(edge->from) > 0); //can not be a cap instance
					assert(lengthGreyEdges(edge->to) > 0); //can not be a cap instance
					greyEdgeIterator = getGreyEdgeIterator(edge->from);
					while((vertex2 = getNextGreyEdge(edge->from, greyEdgeIterator)) != NULL) {
						assert(vertex2->vertexID != 0); //can not be connected to the source vertex
					}
					destructGreyEdgeIterator(greyEdgeIterator);
					greyEdgeIterator = getGreyEdgeIterator(edge->to);
					while((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
						assert(vertex2->vertexID != 0); //can not be connected to the source vertex
					}
					destructGreyEdgeIterator(greyEdgeIterator);
				}
			}
			destructBlackEdgeIterator(blackEdgeIterator);

			//check they all point consistently.
			void *blackEdgeIterator = getBlackEdgeIterator(vertex);
			edge = getNextBlackEdge(vertex, blackEdgeIterator);
			while((edge2 = getNextBlackEdge(vertex, blackEdgeIterator)) != NULL) {
				//check edges have the same length.
				assert(edge->segment->end - edge->segment->start == edge2->segment->end - edge2->segment->start);
				//check they point at the same vertex.
				assert(edge->to == edge2->to);
			}
			destructBlackEdgeIterator(blackEdgeIterator);
		}

		//check the grey edges.
		greyEdgeIterator = getGreyEdgeIterator(vertex);
		while((vertex2 = getNextGreyEdge(vertex, greyEdgeIterator)) != NULL) {
			//check the other vertex has a valid id.
			assert(vertex2->vertexID >= 0);
			assert(vertex2->vertexID < graph->vertices->length);

			//check it has a grey edge pointing back.
			assert(containsGreyEdge(vertex2, vertex));
		}
		destructGreyEdgeIterator(greyEdgeIterator);
	}
	logDebug("Checked the graph okay\n");

#endif
}
