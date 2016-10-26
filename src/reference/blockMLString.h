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

#endif /* BLOCKMLSTRING_H_ */
