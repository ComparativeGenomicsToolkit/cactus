#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Public pseudo-adjacency functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

PseudoAdjacency *pseudoAdjacency_construct(End *_5End, End *_3End,
		PseudoChromosome *pseudoChromosome) {
	return pseudoAdjacency_construct2(cactusDisk_getUniqueID(flower_getNetDisk(reference_getNet(pseudoChromosome_getReference(pseudoChromosome)))),
			_5End, _3End, pseudoChromosome);
}

Name pseudoAdjacency_getName(PseudoAdjacency *pseudoAdjacency) {
	return pseudoAdjacency->name;
}

End *pseudoAdjacency_get5End(PseudoAdjacency *pseudoAdjacency) {
	return pseudoAdjacency->_5End;
}

End *pseudoAdjacency_get3End(PseudoAdjacency *pseudoAdjacency) {
	return pseudoAdjacency->_3End;
}

PseudoChromosome *pseudoAdjacency_getPseudoChromosome(PseudoAdjacency *pseudoAdjacency) {
	return pseudoAdjacency->pseudoChromosome;
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private pseudo-adjacency functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

PseudoAdjacency *pseudoAdjacency_construct2(Name name,
		End *_5End, End *_3End,
		PseudoChromosome *pseudoChromosome) {
	PseudoAdjacency *pseudoAdjacency = st_malloc(sizeof(PseudoAdjacency));
	assert(_5End != NULL);
	assert(_3End != NULL);
	assert(name != NULL_NAME);
	assert(pseudoChromosome != NULL);
	assert(end_getGroup(_5End) == end_getGroup(_3End));
	pseudoAdjacency->name = name;
	pseudoAdjacency->_5End = end_getPositiveOrientation(_5End);
	pseudoAdjacency->_3End = end_getPositiveOrientation(_3End);
	pseudoAdjacency->pseudoChromosome = pseudoChromosome;
	pseudoChromosome_addPseudoAdjacency(pseudoChromosome, pseudoAdjacency);
	return pseudoAdjacency;
}

void pseudoAdjacency_destruct(PseudoAdjacency *pseudoAdjacency) {
	pseudoChromosome_removePseudoAdjacency(pseudoAdjacency_getPseudoChromosome(pseudoAdjacency),
			pseudoAdjacency);
	free(pseudoAdjacency);
}

PseudoAdjacency pseudoAdjacency_getStaticNameWrapperP;
PseudoAdjacency *pseudoAdjacency_getStaticNameWrapper(Name name) {
	pseudoAdjacency_getStaticNameWrapperP.name = name;
	return &pseudoAdjacency_getStaticNameWrapperP;
}

void pseudoAdjacency_writeBinaryRepresentation(PseudoAdjacency *pseudoAdjacency, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeElementType(CODE_PSEUDO_ADJACENCY, writeFn);
	binaryRepresentation_writeName(pseudoAdjacency_getName(pseudoAdjacency), writeFn);
	binaryRepresentation_writeName(end_getName(pseudoAdjacency_get5End(pseudoAdjacency)), writeFn);
	binaryRepresentation_writeName(end_getName(pseudoAdjacency_get3End(pseudoAdjacency)), writeFn);
}

PseudoAdjacency *pseudoAdjacency_loadFromBinaryRepresentation(void **binaryString, PseudoChromosome *pseudoChromosome) {
	PseudoAdjacency *pseudoAdjacency = NULL;
	End *_5End, *_3End;
	Name name;
	if(binaryRepresentation_peekNextElementType(*binaryString) == CODE_PSEUDO_ADJACENCY) {
		binaryRepresentation_popNextElementType(binaryString);
		name = binaryRepresentation_getName(binaryString);
		_5End = flower_getEnd(reference_getNet(pseudoChromosome_getReference(pseudoChromosome)), binaryRepresentation_getName(binaryString));
		_3End = flower_getEnd(reference_getNet(pseudoChromosome_getReference(pseudoChromosome)), binaryRepresentation_getName(binaryString));
		pseudoAdjacency = pseudoAdjacency_construct2(name, _5End, _3End, pseudoChromosome);
	}
	return pseudoAdjacency;
}
