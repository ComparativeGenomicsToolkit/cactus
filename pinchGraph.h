#ifndef PINCH_GRAPH_H_
#define PINCH_GRAPH_H_

#include "fastCMaths.h"
#include "commonC.h"
#include "hashTableC.h"
#include "pairwiseAlignment.h"
#include "avl.h"
//#include "net.h"
#include "cactus.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Segments
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct Segment {
	/*
	 * Type for representing a sub-sequence.
	 */
    //The sequence from which it comes.
	Name contig;
	//inclusive start coordinate
	int32_t start;
	//inclusive end coordinate.
	int32_t end;
	//the mirror segment in the opposite direction (this is always present
	//make switching between the forward and reverse complement really simple
	//and symmetric.
	struct Segment *rSegment;
};

void segment_recycle(struct Segment *segment, Name contig, int32_t start, int32_t end);

struct Segment *constructSegment(Name contig, int32_t start, int32_t end);

void destructSegment(struct Segment *segment);

void logSegment(struct Segment *segment);

int segmentComparator(struct Segment *segment1, struct Segment *segment2);

int segmentComparatorPointers(struct Segment **segment1, struct Segment **segment2);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic elements of the pinch graph.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchVertex {
	/*
	 * The basic vertex type of the pinch graph.
	 */
	//these are the pinch edges.
	struct avl_table *blackEdges;
	//these are just pointers to other vertices.
	struct avl_table *greyEdges;
	//a unique ID for each vertex.
	int32_t vertexID;
	bool isEnd;
	bool isDeadEnd;
};


struct PinchEdge {
	/*
	 * The basic 'segment' containing edge type of the pinch graph. (as opposed
	 * to grey edges which are not explicitly defined).
	 */
	struct PinchVertex *from;
	struct PinchVertex *to;
	struct Segment *segment;
	//we store everything forwards and backwards to make updating symmetrical.
	struct PinchEdge *rEdge;
};


struct PinchGraph {
	/*
	 * The type for holding the graph.
	 */
	struct avl_table *edges;
	struct List *vertices;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Vertex methods
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchVertex *constructPinchVertex(struct PinchGraph *graph, int32_t vertexID, bool isEnd, bool isDeadEnd);

void destructPinchVertex(struct PinchVertex *);


struct PinchVertex *mergeVertices(struct PinchGraph *graph, struct PinchVertex *vertex1, struct PinchVertex *vertex2);

void removeVertexFromGraphAndDestruct(struct PinchGraph *graph, struct PinchVertex *vertex);

/*
 * Returns true if the vertex is an end.
 */
bool vertex_isEnd(struct PinchVertex *vertex);

/*
 * Returns true if the vertex is a dead end.
 */
bool vertex_isDeadEnd(struct PinchVertex *vertex);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for accessing adjacency (grey) edges
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void *constructGreyEdges();

void destructGreyEdges(void *);

int32_t lengthGreyEdges(struct PinchVertex *vertex);

struct PinchVertex *getFirstGreyEdge(struct PinchVertex *vertex);

int32_t containsGreyEdge(struct PinchVertex *vertex, struct PinchVertex *vertex2);

void insertGreyEdge(struct PinchVertex *vertex, struct PinchVertex *vertex2);

void removeGreyEdge(struct PinchVertex *vertex, struct PinchVertex *vertex2);

struct PinchVertex *popGreyEdge(struct PinchVertex *vertex);

void *getGreyEdgeIterator(struct PinchVertex *vertex);

struct PinchVertex *getNextGreyEdge(struct PinchVertex *vertex, void *iterator);

void destructGreyEdgeIterator(void *iterator);

void connectVertices(struct PinchVertex *vertex, struct PinchVertex *vertex2);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Edge methods.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchEdge *constructPinchEdge(struct Segment *);

void destructPinchEdge(struct PinchEdge *);

void removePinchEdgeFromGraphAndDestruct(struct PinchGraph *graph, struct PinchEdge *edge);

void addPinchEdgeToGraph(struct PinchGraph *graph, struct PinchEdge *edge);

void connectPinchEdge(struct PinchEdge *edge, struct PinchVertex *from, struct PinchVertex *to);

int32_t edgeComparator(struct PinchEdge *edge1, struct PinchEdge *edge2, void *o);

int32_t isAStubOrCap(struct PinchEdge *edge);

/*
 * Types for choosing the sides of positions.
 */
#define LEFT 0
#define RIGHT 1

struct PinchEdge *getContainingBlackEdge(struct PinchGraph *graph, Name contig, int32_t position);

struct PinchEdge *getNextEdge(struct PinchGraph *graph, struct PinchEdge *edge, Net *net);

struct PinchVertex *splitEdge(struct PinchGraph *graph, Name contig, int32_t position, int32_t leftOrRight);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for accessing pinch (black) edges
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void *constructBlackEdges();

void destructBlackEdges(void *);

int32_t lengthBlackEdges(struct PinchVertex *vertex);

struct PinchEdge *getFirstBlackEdge(struct PinchVertex *vertex);

int32_t containsBlackEdge(struct PinchVertex *vertex, struct PinchEdge *edge);

void insertBlackEdge(struct PinchVertex *vertex, struct PinchEdge *edge);

void removeBlackEdge(struct PinchVertex *vertex, struct PinchEdge *edge);

struct PinchEdge *popBlackEdge(struct PinchVertex *vertex);

void *getBlackEdgeIterator(struct PinchVertex *vertex);

struct PinchEdge *getNextBlackEdge(struct PinchVertex *vertex, void *iterator);

void destructBlackEdgeIterator(void *iterator);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//The actual graph construction methods
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs an empty pinch graph.
 */
struct PinchGraph *pinchGraph_construct();

void destructPinchGraph(struct PinchGraph *);

void checkPinchGraph(struct PinchGraph *graph);

void checkPinchGraphDegree(struct PinchGraph *graph, int32_t maxDegree);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for 'pinching' the graph.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void pinchMergeSegment(struct PinchGraph *graph,
					   struct Segment *segment1,
					   struct Segment *segment2,
					   struct hashtable *vertexAdjacencyComponents);

void pinchMerge(struct PinchGraph *graph, struct PairwiseAlignment *pairwiseAlignment,
		void (*addFunction)(struct PinchGraph *pinchGraph, struct Segment *, struct Segment *, struct hashtable *, void *),
		void *extraParameter,
		struct hashtable *vertexAdjacencyComponents);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods to get the virtual segments that lie between two
//vertices connected by a grey edge.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct List *getGreyEdgeSegments(struct PinchVertex *vertex, struct PinchVertex *vertex2);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for writing out the basic pinch graph
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void writeOutPinchGraphWithChains(struct PinchGraph *pinchGraph,
								  struct List *chainsList,
								  struct List *groups,
								  FILE *fileHandle);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods from pinchGraphManipulation.c
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void removeOverAlignedEdges(struct PinchGraph *pinchGraph, float minimumTreeCoverage, int32_t maxDegree, int32_t extensionSteps, Net *net);

struct List *getRecursiveComponents(struct PinchGraph *pinchGraph, int32_t (*excludedEdgesFn)(void *));

struct List *getRecursiveComponents2(struct PinchGraph *pinchGraph, struct List *edgesToExclude);

void linkStubComponentsToTheSinkComponent(struct PinchGraph *pinchGraph);

void removeTrivialGreyEdgeComponents(struct PinchGraph *graph, struct List *listOfVertices, Net *net);

float treeCoverage(struct PinchVertex *vertex, Net *net,
		struct PinchGraph *pinchGraph);


#endif
