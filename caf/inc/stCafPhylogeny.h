/*
 * stCafPhylogeny.h
 *
 *  Created on: Jun 2, 2014
 *      Author: benedictpaten
 */

#ifndef STCAFPHYLOGENY_H_
#define STCAFPHYLOGENY_H_

/*
 * Build tree for each block and then use it to partition homologies in the block into
 * those which occur before and after the speciation.
 */
void stCaf_buildTreesToRemoveAncientHomologies(stPinchThreadSet *threadSet, stHash *threadStrings, stSet *outgroupThreads);

#endif /* STCAFPHYLOGENY_H_ */
