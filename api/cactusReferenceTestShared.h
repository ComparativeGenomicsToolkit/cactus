#include "cactusGlobalsPrivate.h"

static CactusDisk *cactusDisk = NULL;
static Net *net;
static bool nestedTest = 0;

static End *end1, *end2, *end3, *end4, *end5, *end6, *end7, *end8;
static Block *block1, *block2;
static Reference *reference;
static PseudoChromosome *pseudoChromosome1;
static PseudoChromosome *pseudoChromosome2;
static PseudoAdjacency *pseudoAdjacency1;
static PseudoAdjacency *pseudoAdjacency2;
static PseudoAdjacency *pseudoAdjacency3;
static PseudoAdjacency *pseudoAdjacency4;

static void cactusReferenceTestSharedTeardown() {
	if(cactusDisk != NULL) {
		cactusDisk_destruct(cactusDisk);
		testCommon_deleteTemporaryNetDisk();
		cactusDisk = NULL;
	}
}

static void cactusReferenceTestSharedSetup() {
	cactusReferenceTestSharedTeardown();
	cactusDisk = cactusDisk_construct(testCommon_getTemporaryNetDisk());
	net = net_construct(cactusDisk);

	end1 = end_construct(0, net);
	block1 = block_construct(1, net);
	end2 = block_get5End(block1);
	end3 = block_get3End(block1);
	block2 = block_construct(1, net);
	end4 = block_get5End(block2);
	end5 = block_get3End(block2);
	end6 = end_construct(0, net);
	end7 = end_construct(0, net);
	end8 = end_construct(0, net);

	reference = reference_construct(net);

	pseudoChromosome1 = pseudoChromosome_construct(reference, end1, end6);
	pseudoChromosome2 = pseudoChromosome_construct(reference, end7, end8);

	pseudoAdjacency1 = pseudoAdjacency_construct(end1, end2, pseudoChromosome1);
	pseudoAdjacency2 = pseudoAdjacency_construct(end3, end4, pseudoChromosome1);
	pseudoAdjacency3 = pseudoAdjacency_construct(end5, end6, pseudoChromosome1);

	pseudoAdjacency4 = pseudoAdjacency_construct(end7, end8, pseudoChromosome2);
}
