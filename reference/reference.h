/*
 * reference.h
 *
 *  Created on: 1 Apr 2010
 *      Author: benedictpaten
 */

#ifndef REFERENCE_H_
#define REFERENCE_H_

/////////////////////
//Core functions
/////////////////////

/*
 * Makes the pseudo adjacencies for each pseudo chromosome in reference.
 * In a top-down greedy definition of a reference this is the only function that can
 * be varied.
 */
void makePseudoAdjacencies(Net *net, Reference *reference);

/*
 * Makes a pseudo chromosome for each pair of attached ends in the reference. Assumes
 * that each attached end in the top level problem has only one cap.. as we currently do
 * with setup, if we merge these ends this will have to change.
 */
void makeTopLevelPseudoChromosomes(Net *net, Reference *reference);

/*
 * Gets the parent group of the net and it's reference (which must be defined),
 * then for each set of attached ends in the net it gathers the corresponding pairing
 * in the higher level reference and constructs a corresponding set of pseudo chromosomes,
 * one for each pair of attached ends.
 */
void makeIntermediateLevelPseudoChromosomes(Net *net, Reference *reference);

/*
 * After constructing a reference pseudo-adjacencies may link groups together
 * that were previously separate. To correct this, this function merges distinct groups
 * joined by novel pseudo-adjacencies.
 */
void mergeGroupsLinkedByPseudoAdjacencies(Net *net, Reference *reference);

/*
 * This is the core function called to create a reference with the given name for a genome.
 * If the reference already exists then it does nothing.
 */
void addReferenceToNet(Net *net);

////////////////////////
//Misc functions
////////////////////////

/*
 * Gets a group in which to place an end, using a bunch of heuristics.
 */
Group *getSpareGroup(Net *net);

/*
 * Ensures the given end is in all the child nets.
 */
void pushEndIntoChildNets(End *end);

/*
 * Returns the number of free stub ends in the group.
 */
int32_t getFreeStubEndNumber(Group *group);


#endif /* REFERENCE_H_ */
