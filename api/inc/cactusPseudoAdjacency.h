#ifndef CACTUS_PSEUDO_ADJACENCY_H_
#define CACTUS_PSEUDO_ADJACENCY_H_

#include "cactusGlobals.h"

/*
 * Constructs a new pseudo-adjacency.
 */
PseudoAdjacency *pseudoAdjacency_construct(End *_5End, End *_3End,
		PseudoChromosome *pseudoChromosome);

/*
 * Gets the name of the pseudo adjacency.
 */
Name pseudoAdjacency_getName(PseudoAdjacency *pseudoAdjacency);

/*
 * Gets the 5End of the pseudo adjacency.
 */
End *pseudoAdjacency_get5End(PseudoAdjacency *pseudoAdjacency);

/*
 * Gets the 3End of the pseudo adjacency.
 */
End *pseudoAdjacency_get3End(PseudoAdjacency *pseudoAdjacency);

/*
 * Gets the pseudo-chromosome that contains the end, or NULL if not set.
 */
PseudoChromosome *pseudoAdjacency_getPseudoChromosome(PseudoAdjacency *pseudoAdjacency);

/*
 * Get the index (position) of the pseudo adjacency in its parent pseudo chromosome.
 */
int32_t pseudoAdjacency_getIndex(PseudoAdjacency *pseudoAdjacency);

#endif
