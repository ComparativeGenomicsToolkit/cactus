#ifndef PINCH_GRAPH_H_
#define PINCH_GRAPH_H_

#include "fastCMaths.h"
#include "commonC.h"
#include "hashTableC.h"
#include "pairwiseAlignment.h"
#include "avl.h"
//#include "flower.h"
#include "cactus.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Pieces
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct Piece {
	/*
	 * Type for representing a sub-sequence.
	 */
    //The sequence from which it comes.
	Name contig;
	//inclusive start coordinate
	int32_t start;
	//inclusive end coordinate.
	int32_t end;
	//the mirror piece in the opposite direction (this is always present
	//make switching between the forward and reverse complement really simple
	//and symmetric.
	struct Piece *rPiece;
};

void piece_recycle(struct Piece *piece, Name contig, int32_t start, int32_t end);

struct Piece *constructPiece(Name contig, int32_t start, int32_t end);

void destructPiece(struct Piece *piece);

void logPiece(struct Piece *piece);

int pieceComparator(struct Piece *piece1, struct Piece *piece2);

int pieceComparatorPointers(struct Piece **piece1, struct Piece **piece2);

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
	 * The basic 'piece' containing edge type of the pinch graph. (as opposed
	 * to grey edges which are not explicitly defined).
	 */
	struct PinchVertex *from;
	struct PinchVertex *to;
	struct Piece *piece;
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

void disconnectVertices(struct PinchVertex *vertex, struct PinchVertex *vertex2);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Edge methods.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct PinchEdge *constructPinchEdge(struct Piece *);

void destructPinchEdge(struct PinchEdge *);

void removePinchEdgeFromGraphAndDestruct(struct PinchGraph *graph, struct PinchEdge *edge);

void addPinchEdgeToGraph(struct PinchGraph *graph, struct PinchEdge *edge);

void connectPinchEdge(struct PinchEdge *edge, struct PinchVertex *from, struct PinchVertex *to);

int32_t edgeComparator(struct PinchEdge *edge1, struct PinchEdge *edge2, void *o);

int32_t isAStub(struct PinchEdge *edge);

/*
 * Types for choosing the sides of positions.
 */
#define LEFT 0
#define RIGHT 1

struct PinchEdge *getContainingBlackEdge(struct PinchGraph *graph, Name contig, int32_t position);

struct PinchEdge *getNextEdge(struct PinchGraph *graph, struct PinchEdge *edge, Flower *flower);

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
//Methods to get the virtual pieces that lie between two
//vertices connected by a grey edge.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct List *getGreyEdgePieces(struct PinchVertex *vertex, struct PinchVertex *vertex2);

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Methods for writing out the basic pinch graph
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void writeOutPinchGraphWithChains(struct PinchGraph *pinchGraph,
								  struct hashtable *edgeColours,
								  struct List *groups,
								  FILE *fileHandle);


#endif
