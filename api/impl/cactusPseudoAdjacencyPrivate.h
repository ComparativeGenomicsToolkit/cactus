#ifndef CACTUS_PSEUDO_ADJACENCY_PRIVATE_H_
#define CACTUS_PSEUDO_ADJACENCY_PRIVATE_H_

#include "cactusGlobals.h"

struct _pseudoAdjacency {
	Name name;
	End *_5End, *_3End;
	PseudoChromosome *pseudoChromosome;
};

/*
 * Constructs a new pseudo-adjacency.
 */
PseudoAdjacency *pseudoAdjacency_construct2(Name name,
		End *_5End, End *_3End,
		PseudoChromosome *pseudoChromosome);

/*
 * Destruct the pseudo-adjacency.
 */
void pseudoAdjacency_destruct(PseudoAdjacency *pseudoAdjacency);

/*
 * Gets name of pseudo adjacency.
 */
PseudoAdjacency *pseudoAdjacency_getStaticNameWrapper(Name name);

/*
 * Write a binary representation of the pseudo-adjacency to the write function.
 */
void pseudoAdjacency_writeBinaryRepresentation(PseudoAdjacency *pseudoAdjacency, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a pseudo-chromosome into memory from a binary representation.
 */
PseudoAdjacency *pseudoAdjacency_loadFromBinaryRepresentation(void **binaryString, PseudoChromosome *pseudoChromosome);

#endif
