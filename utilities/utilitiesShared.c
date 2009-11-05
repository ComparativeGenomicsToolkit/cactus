#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"

double calculateTotalContainedSequence(Net *net) {
	Net_EndIterator *endIterator = net_getEndIterator(net);
	End *end;
	double totalLength = 0.0;
	while((end = net_getNextEnd(endIterator)) != NULL) {
		if(!end_isAtomEnd(end)) {
			End_InstanceIterator *instanceIterator = end_getInstanceIterator(end);
			EndInstance *endInstance;
			while((endInstance = end_getNext(instanceIterator)) != NULL) {
				endInstance = endInstance_getStrand(endInstance) ? endInstance : endInstance_getReverse(endInstance);
				if(!endInstance_getSide(endInstance)) {
					EndInstance *endInstance2 = endInstance_getAdjacency(endInstance);
					while(end_isAtomEnd(endInstance_getEnd(endInstance2))) {
						AtomInstance *atomInstance = endInstance_getAtomInstance(endInstance2);
						assert(atomInstance != NULL);
						assert(atomInstance_get5End(atomInstance) == endInstance2);
						endInstance2 = endInstance_getAdjacency(atomInstance_get3End(atomInstance));
						assert(endInstance_getStrand(endInstance2));
						assert(endInstance_getSide(endInstance2));
					}
					assert(endInstance_getStrand(endInstance2));
					assert(endInstance_getSide(endInstance2));
					int32_t length = endInstance_getCoordinate(endInstance2) - endInstance_getCoordinate(endInstance) - 1;
					assert(length >= 0);
					totalLength += length;
				}
			}
			end_destructInstanceIterator(instanceIterator);
		}
	}
	net_destructEndIterator(endIterator);
	return totalLength;
}

double calculateChainSize(Chain *chain) {
	Atom **atoms;
	int32_t i, j;
	double k = 0.0;
	atoms = chain_getAtomChain(chain, &i);
	for(j=0; j<i; j++) {
		k += atom_getLength(atoms[j]);
	}
	free(atoms);
	for(i=0; i<chain_getLength(chain); i++) {
		k += calculateTotalContainedSequence(adjacencyComponent_getNestedNet(link_getAdjacencyComponent(chain_getLink(chain, i))));
	}
	return k;
}
