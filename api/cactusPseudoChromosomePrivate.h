#ifndef CACTUS_PSEUDO_CHROMOSOME_PRIVATE_H_
#define CACTUS_PSEUDO_CHROMOSOME_PRIVATE_H_

#include "cactusGlobals.h"

struct _pseudoChromosome {
	Name name;
	stSortedSet *pseudoAdjacencies;
	End *_5End;
	End *_3End;
	Reference *reference;
};

/*
 * Constructor with name prespecified.
 */
PseudoChromosome *pseudoChromosome_construct2(Name name, Reference *reference,
		End *_5End, End *_3End);

/*
 * Destructor for pseudo-chromosome.
 */
void pseudoChromosome_destruct(PseudoChromosome *pseudoChromosome);

/*
 * Adds a pseudo-adjacency to the pseudo-chromosome structure.
 */
void pseudoChromosome_addPseudoAdjacency(PseudoChromosome *pseudoChromosome, PseudoAdjacency *pseudoAdjacency);

/*
 * Removes a pseudo-adjacency from the pseudo-chromosome structure.
 */
void pseudoChromosome_removePseudoAdjacency(PseudoChromosome *pseudoChromosome, PseudoAdjacency *pseudoAdjacency);

/*
 * Static name wrapper.
 */
PseudoChromosome *pseudoChromosome_getStaticNameWrapper(Name name);

/*
 * Write a binary representation of the pseudo-chromosome to the write function.
 */
void pseudoChromosome_writeBinaryRepresentation(PseudoChromosome *pseudoChromosome, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a pseudo-chromosome into memory from a binary representation.
 */
PseudoChromosome *pseudoChromosome_loadFromBinaryRepresentation(void **binaryString, Reference *reference);

#endif
