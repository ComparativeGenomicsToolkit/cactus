#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Basic reference functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static int reference_constructP(const void *o1, const void *o2) {
    return cactusMisc_nameCompare(pseudoChromosome_getName(
            (PseudoChromosome *) o1), pseudoChromosome_getName(
            (PseudoChromosome *) o2));
}

Reference *reference_construct(Flower *flower) {
    Reference *reference = st_malloc(sizeof(Reference));
    //Setup the basic structure - a sorted set of pseudo-chromosomes.
    reference->pseudoChromosomes = stSortedSet_construct3(reference_constructP,
            NULL);
    //Link the reference and flower.
    reference->flower = flower;
    flower_setReference(flower, reference);
    return reference;
}

Flower *reference_getFlower(Reference *reference) {
    return reference->flower;
}

int32_t reference_getPseudoChromosomeNumber(Reference *reference) {
    return stSortedSet_size(reference->pseudoChromosomes);
}

PseudoChromosome *reference_getPseudoChromosome(Reference *reference, Name name) {
    PseudoChromosome *pseudoChromosome;
    pseudoChromosome = pseudoChromosome_getStaticNameWrapper(name);
    return stSortedSet_search(reference->pseudoChromosomes, pseudoChromosome);
}

PseudoChromosome *reference_getFirst(Reference *reference) {
    return stSortedSet_getFirst(reference->pseudoChromosomes);
}

Reference_PseudoChromosomeIterator *reference_getPseudoChromosomeIterator(
        Reference *reference) {
    return stSortedSet_getIterator(reference->pseudoChromosomes);
}

PseudoChromosome *reference_getNextPseudoChromosome(
        Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
    return stSortedSet_getNext(pseudoChromosomeIterator);
}

PseudoChromosome *reference_getPreviousPseudoChromosome(
        Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
    return stSortedSet_getPrevious(pseudoChromosomeIterator);
}

Reference_PseudoChromosomeIterator *reference_copyPseudoChromosomeIterator(
        Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
    return stSortedSet_copyIterator(pseudoChromosomeIterator);
}

void reference_destructPseudoChromosomeIterator(
        Reference_PseudoChromosomeIterator *pseudoChromosomeIterator) {
    stSortedSet_destructIterator(pseudoChromosomeIterator);
}

void reference_check(Reference *reference) {
    Flower *flower = reference_getFlower(reference);
    //Going ends --> pseudo adjacencies.
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        if (end_isAttached(end) || end_isBlockEnd(end)) {
            PseudoAdjacency *pseudoAdjacency = end_getPseudoAdjacency(end);
            assert(pseudoAdjacency != NULL);
            assert(pseudoAdjacency_get5End(pseudoAdjacency) == end || pseudoAdjacency_get3End(pseudoAdjacency) == end);
        } else {
            assert(end_getPseudoAdjacency(end) == NULL); //check free stub end is not in the pseudo chromosomes..
        }
    }
    flower_destructEndIterator(endIterator);
    assert(reference_getPseudoChromosomeNumber(reference)*2 == flower_getAttachedStubEndNumber(flower));

    //Going pseudo-adjacencies --> ends.
    Reference_PseudoChromosomeIterator *pseudoChromosomeIterator =
            reference_getPseudoChromosomeIterator(reference);
    PseudoChromosome *pseudoChromosome;
    int32_t i = 0;
    while ((pseudoChromosome = reference_getNextPseudoChromosome(
            pseudoChromosomeIterator)) != NULL) {
        //Here we check the structure of the pseudo chromosome also..
        assert(pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome) > 0); //must be at least one adjacency
        assert(pseudoChromosome_get5End(pseudoChromosome) == pseudoAdjacency_get5End(pseudoChromosome_getPseudoAdjacencyByIndex(pseudoChromosome, 0))); //check the 5 end matches the 5 end of the first pseudo adjacency.
        assert(pseudoChromosome_get3End(pseudoChromosome) == pseudoAdjacency_get3End(pseudoChromosome_getPseudoAdjacencyByIndex(pseudoChromosome, pseudoChromosome_getPseudoAdjacencyNumber(pseudoChromosome)-1))); //check the 5 end matches the 5 end of the first pseudo adjacency.

        //We check we are consistent with our parent reference
        if(flower_getParentGroup(flower) != NULL) {
            Group *parentGroup = flower_getParentGroup(flower);
            if(flower_getReference(group_getFlower(parentGroup)) != NULL) {
                End *parent5End = group_getEnd(parentGroup, end_getName(pseudoChromosome_get5End(pseudoChromosome)));
                End *parent3End = group_getEnd(parentGroup, end_getName(pseudoChromosome_get3End(pseudoChromosome)));
                assert(parent5End != NULL);
                assert(parent3End != NULL);
                PseudoAdjacency *parentPseudoAdjacency = end_getPseudoAdjacency(parent5End);
                assert(parentPseudoAdjacency != NULL);
                assert(pseudoAdjacency_get5End(parentPseudoAdjacency) == parent3End || pseudoAdjacency_get5End(parentPseudoAdjacency) == parent5End);
            }
        }

        PseudoChromsome_PseudoAdjacencyIterator *pseudoAdjacencyIterator =
                pseudoChromosome_getPseudoAdjacencyIterator(pseudoChromosome);
        PseudoAdjacency *pseudoAdjacency, *previousPseudoAdjacency = NULL;
        while ((pseudoAdjacency = pseudoChromosome_getNextPseudoAdjacency(
                pseudoAdjacencyIterator)) != NULL) {
            i++;
            End *_5End = pseudoAdjacency_get5End(pseudoAdjacency);
            End *_3End = pseudoAdjacency_get3End(pseudoAdjacency);
            assert(_5End == end_getPositiveOrientation(_5End)); //check they are positive orientation
            assert(_3End == end_getPositiveOrientation(_3End)); //check they are positive orientation
            assert(end_getPseudoAdjacency(_5End) == pseudoAdjacency); //check these are represented in the hash.. so that the mapping is unique.
            assert(end_getPseudoAdjacency(_3End) == pseudoAdjacency);
            assert(end_getGroup(_5End) != NULL); //check the groups are the same for both sides of the adjacency.
            assert(end_getGroup(_5End) == end_getGroup(_3End));
            assert(group_getLink(end_getGroup(_5End)) == group_getLink(end_getGroup(_3End))); //check if there is a link, they are in the same link.
            if (previousPseudoAdjacency != NULL) { //check the adjacency spans the block...
                assert(pseudoAdjacency_get3End(previousPseudoAdjacency) == end_getOtherBlockEnd(_5End));
                assert(_5End == end_getOtherBlockEnd(pseudoAdjacency_get3End(previousPseudoAdjacency))); //do it the other way, just for fun.
            }
            previousPseudoAdjacency = pseudoAdjacency;
        }
        pseudoChromosome_destructPseudoAdjacencyIterator(
                pseudoAdjacencyIterator);
    }
    reference_destructPseudoChromosomeIterator(pseudoChromosomeIterator);
    assert(i*2 == flower_getAttachedStubEndNumber(flower) + flower_getBlockEndNumber(flower)); //stHash_size(endsToPseudoAdjacencies));
}

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Private functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void reference_destruct(Reference *reference) {
    flower_removeReference(reference_getFlower(reference), reference);
    PseudoChromosome *pseudoChromosome;
    while ((pseudoChromosome = reference_getFirst(reference)) != NULL) {
        pseudoChromosome_destruct(pseudoChromosome);
    }
    stSortedSet_destruct(reference->pseudoChromosomes);
    free(reference);
}

void reference_addPseudoChromosome(Reference *reference,
        PseudoChromosome *pseudoChromosome) {
    assert(stSortedSet_search(reference->pseudoChromosomes, pseudoChromosome) == NULL);
    stSortedSet_insert(reference->pseudoChromosomes, pseudoChromosome);
}

void reference_removePseudoChromosome(Reference *reference,
        PseudoChromosome *pseudoChromosome) {
    assert(stSortedSet_search(reference->pseudoChromosomes, pseudoChromosome) != NULL);
    stSortedSet_remove(reference->pseudoChromosomes, pseudoChromosome);
}

void reference_writeBinaryRepresentation(Reference *reference, void(*writeFn)(
        const void * ptr, size_t size, size_t count)) {
    Reference_PseudoChromosomeIterator *iterator;
    PseudoChromosome *pseudoChromosome;
    binaryRepresentation_writeElementType(CODE_REFERENCE, writeFn);
    binaryRepresentation_writeInteger(reference_getPseudoChromosomeNumber(
            reference), writeFn);

    iterator = reference_getPseudoChromosomeIterator(reference);
    while ((pseudoChromosome = reference_getNextPseudoChromosome(iterator))
            != NULL) {
        pseudoChromosome_writeBinaryRepresentation(pseudoChromosome, writeFn);
    }
    reference_destructPseudoChromosomeIterator(iterator);
    binaryRepresentation_writeElementType(CODE_REFERENCE, writeFn);
}

Reference *reference_loadFromBinaryRepresentation(void **binaryString,
        Flower *flower) {
    Reference *reference = NULL;
    int32_t pseudoChromosomeNumber;

    if (binaryRepresentation_peekNextElementType(*binaryString)
            == CODE_REFERENCE) {
        binaryRepresentation_popNextElementType(binaryString);
        reference = reference_construct(flower);
        pseudoChromosomeNumber = binaryRepresentation_getInteger(binaryString);
        while (pseudoChromosomeNumber-- > 0) {
            pseudoChromosome_loadFromBinaryRepresentation(binaryString,
                    reference);
        }
        assert(binaryRepresentation_peekNextElementType(*binaryString)
            == CODE_REFERENCE);
        binaryRepresentation_popNextElementType(binaryString);
    }
    return reference;
}
