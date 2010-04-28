#ifndef CACTUS_REFERENCE_PRIVATE_H_
#define CACTUS_REFERENCE_PRIVATE_H_

#include "cactusGlobals.h"

struct _reference {
	struct avl_table *pseudoChromosomes;
	Net *net;
	Name name;
};

/*
 * Destructor for reference.
 */
void reference_destruct(Reference *reference);

/*
 * Adds a pseudo-chromosome to the reference structure.
 */
void reference_addPseudoChromosome(Reference *reference, PseudoChromosome *pseudoChromosome);

/*
 * Remove the pseudo-chromosome from the reference.
 */
void reference_removePseudoChromosome(Reference *reference, PseudoChromosome *pseudoChromosome);

/*
 * Gets another reference wrapped name.
 */
Reference *reference_getStaticNameWrapper(Name name);

/*
 * Write a binary representation of the reference to the write function.
 */
void reference_writeBinaryRepresentation(Reference *reference, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Loads a reference into memory from a binary representation of the net.
 */
Reference *reference_loadFromBinaryRepresentation(void **binaryString, Net *net);

/*
 * Sets the net associated with the reference.
 */
void reference_setNet(Reference *reference, Net *net);

#endif
