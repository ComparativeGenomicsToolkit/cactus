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
	return pseudoAdjacency_construct2(cactusDisk_getUniqueID(flower_getCactusDisk(reference_getFlower(pseudoChromosome_getReference(pseudoChromosome)))),
			_5End, _3End, pseudoChromosome, pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome));
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

int32_t pseudoAdjacency_getIndex(PseudoAdjacency *pseudoAdjacency) {
    return pseudoAdjacency->index;
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
		PseudoChromosome *pseudoChromosome, int32_t index) {
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
	pseudoAdjacency->index = index;
	pseudoChromosome_addPseudoAdjacency(pseudoChromosome, pseudoAdjacency);
	end_setPseudoAdjacency(_5End, pseudoAdjacency);
	end_setPseudoAdjacency(_3End, pseudoAdjacency);
	return pseudoAdjacency;
}

void pseudoAdjacency_destruct(PseudoAdjacency *pseudoAdjacency) {
	end_setPseudoAdjacency(pseudoAdjacency_get3End(pseudoAdjacency), NULL);
	end_setPseudoAdjacency(pseudoAdjacency_get5End(pseudoAdjacency), NULL);
	free(pseudoAdjacency);
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
		_5End = flower_getEnd(reference_getFlower(pseudoChromosome_getReference(pseudoChromosome)), binaryRepresentation_getName(binaryString));
		_3End = flower_getEnd(reference_getFlower(pseudoChromosome_getReference(pseudoChromosome)), binaryRepresentation_getName(binaryString));
		pseudoAdjacency = pseudoAdjacency_construct2(name, _5End, _3End, pseudoChromosome, pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome));
	}
	return pseudoAdjacency;
}
