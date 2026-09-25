/*
 * blockMLString.h
 *
 *  Created on: Nov 26, 2014
 *      Author: benedictpaten
 */

#ifndef BLOCKMLSTRING_H_
#define BLOCKMLSTRING_H_

char *getMaximumLikelihoodString(stTree *tree, Block *block);

stMatrix *generateJukesCantorMatrix(double distance);

stTree *getPhylogeneticTreeRootedAtGivenEvent(Event *event, stMatrix *(*generateSubstitutionMatrix)(double));

Event *getEvent(stTree *tree);

stMatrix *getSubMatrix(stTree *tree);

void cleanupPhylogeneticTree(stTree *tree);

void maskAncestralRepeatBases(Block *block, char *mlString);

/*
 * Take the branch lengths used to reconstruct the reference -- its bases (here) and its adjacencies
 * (getEventWeighting) -- from this tree, matched to events by name, rather than from the event tree.
 * The workflow may have lengthened the event tree's branches above ancestors (upweightAncestorDistances)
 * to make the alignment more sensitive around them, which says nothing about how far apart the genomes are.
 */
void setReconstructionTree(const char *newick);

/*
 * The length of the branch above the event for reconstructing the reference: the length of the path
 * from the event up to its parent in the reconstruction tree, if it was set and has both, else the
 * event tree's.  The two can differ in shape: the event tree drops ancestors left with a single child.
 */
double getReconstructionBranchLength(Event *event);

#endif /* BLOCKMLSTRING_H_ */
