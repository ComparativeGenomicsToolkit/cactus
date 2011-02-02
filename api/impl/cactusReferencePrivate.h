#ifndef CACTUS_REFERENCE_PRIVATE_H_
#define CACTUS_REFERENCE_PRIVATE_H_

#include "cactusGlobals.h"

struct _reference {
	stSortedSet *pseudoChromosomes;
	Flower *flower;
};

/*
 * Adds a pseudo-chromosome to the reference structure.
 */
void reference_addPseudoChromosome(Reference *reference, PseudoChromosome *pseudoChromosome);

/*
 * Remove the pseudo-chromosome from the reference.
 */
void reference_removePseudoChromosome(Reference *reference, PseudoChromosome *pseudoChromosome);

/*
 * Write a binary representation of the reference to the write function.
 */
void reference_writeBinaryRepresentation(Reference *reference, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a reference into memory from a binary representation of the flower.
 */
Reference *reference_loadFromBinaryRepresentation(void **binaryString, Flower *flower);

#endif
