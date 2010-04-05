#ifndef CACTUS_REFERENCE_H_
#define CACTUS_REFERENCE_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Cactus reference ordering functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * A reference ordering structure for the non-stub ends in the net.
 *
 * There are three types of end in any net:
 * pseudo-telomere ends, stub ends and block ends.
 * Pseudo-telomere ends and stub ends are both inherited ends, in that they are
 * contained in any parent net. Pseudo-telomeres are always paired, such that
 * there is an even number of them.
 * Block ends are the ends of the blocks in the considered net
 * and become inherited in any child nets.
 *
 * A reference ordering orders all the pseudo-telomeres and block ends in a net
 * into a set of linear "pseudo chromosomes".
 * Let (a, b) represent an "pseudo-adjacency" between two ends, so called because
 * there may not necessarily exist an actual adjacency between caps representing
 * the two ends.
 *
 * A reference ordering can be represented as a 2-d list of lists of adjacencies.
 *
 * ( ( (L_{0, 0}, R_{0, 0}), (L_{0, 1}, R_{0, 1}) ... (L_{0, N}, R_{0, N}) ),
 *   ( (L_{1, 0}, R_{1, 0}), (L_{1, 1}, R_{1, 1}) ... (L_{1, O}, R_{1, O}) ),
 *   ...
 *   ( (L_{M, 0}, R_{M, 0}), (L_{M, 1}, R_{M, 1}) ... (L_{M, P}, R_{M, P}) ) )
 *
 * Where pair (L_{i, j}, R_{i, j}) represents the jth pseudo-adjacency,
 * where 0 <= j <= N, in the i_th  pseudo-chromosome, where 0 <= i <= M, which
 * contains N psuedo-adjacencies.
 *
 * A reference ordering R is valid for net N only if :
 * (1) the first, L_{i, 0}, and last, R_{i, N}, ends
 * of the first and last pseudo-adjacencies of each pseudo-chromosome
 * represent a pair of paired pseudo-telomeres.
 * Call this property of R "properly terminated"
 * (2) for each link (a, b) in the set of child chains in N there exists
 * a valid pseudo-adjacency, call this property of R "link respecting".
 * (3) for each pair of contiguous adjacencies
 *  ((L_{i, j}, R_{i, j}), (L_{i, j+1}, R_{i, j+1}))
 *  in a pseudo-chromosome, R_{i, j} L_{i, j+1} are opposite ends of the same block.
 *  Call this property of R "block respecting".
 *
 * A reference ordering R is valid for net N if and only if it contains only
 * the complete set of pseudo-telomeres and block ends for N and is properly terminated,
 * link respecting and block respecting.
 */
Reference *reference_construct(Net *net);

/*
 * Gets the name.
 */
Name reference_getName(Reference *reference);

/*
 * Gets net associated with reference.
 */
Net *reference_getNet(Reference *reference);

/*
 * Returns the number of pseudo-chromosomes in the reference.
 */
int32_t reference_getPseudoChromosomeNumber(Reference *reference);

/*
 * Get pseudo-chromosome with the given name.
 */
PseudoChromosome *reference_getPseudoChromosome(Reference *reference, Name name);

/*
 * Gets the first pseudo chromosome in the list of pseudo-chromosomes.
 */
PseudoChromosome *reference_getFirst(Reference *reference);

/*
 * Gets an iterator to iterate through the pseudo-chromosomes in the reference.
 */
Reference_PseudoChromosomeIterator *reference_getPseudoChromosomeIterator(Reference *reference);

/*
 * Gets the next pseudo-chromosome from the iterator.
 */
PseudoChromosome *reference_getNextPseudoChromosome(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator);

/*
 * Gets the previous pseudo-chromosome from the iterator.
 */
PseudoChromosome *reference_getPreviousPseudoChromosome(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator);

/*
 * Duplicates the iterator.
 */
Reference_PseudoChromosomeIterator *reference_copyPseudoChromosomeIterator(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator);

/*
 * Destructs the iterator.
 */
void reference_destructPseudoChromosomeIterator(Reference_PseudoChromosomeIterator *pseudoChromosomeIterator);

/*
 * Gets hash in which keys are ends and the corresponding psuedo-adjacencies are the
 * values. The hash is yours to cleanup and will be unique for each invocation of the function.
 */
Hash *reference_getEndToPseudoAdjacencyHash(Reference *reference);

#endif
