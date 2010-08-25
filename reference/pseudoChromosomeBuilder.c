
#include "sonLib.h"
#include "cactus.h"
#include "adjacencyPairs.h"

/*
 * Code to make pseudo chromosomes for a flower.
 */

stList *getAttachedStubEnds(Flower *flower) {
    /*
     * Get the top level attached ends.
     */
    stList *list = stList_construct();
    Flower_EndIterator *endIterator = flower_getEndIterator(flower);
    End *end;
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        if (end_isAttached(end)) {
            assert(end_isStubEnd(end));
            assert(!end_isFree(end));
            stList_append(list, end);
        }
    }
    flower_destructEndIterator(endIterator);
    assert(stList_length(list) % 2 == 0);
    return list;
}

static void makePseudoChromosomesFromPairs(stList *ends,
        Reference *reference) {
    /*
     * Make pseudo chromosomes from list of paired ends
     */
    assert(stList_length(ends) >= 2); //the reconstruction must contain 2
    //or more attached top level ends to build a reference genome currently.
    assert(stList_length(ends) % 2 == 0); //we must have an even number of ends.
    int32_t i;
    for (i = 0; i < stList_length(ends); i += 2) {
        End *end1 = stList_get(ends, i);
        End *end2 = stList_get(ends, i+1);
        assert(end1 != NULL);
        assert(end2 != NULL);
        assert(end_isStubEnd(end1));
        assert(end_isAttached(end1));
        assert(end_isStubEnd(end2));
        assert(end_isAttached(end2));
        pseudoChromosome_construct(reference, end1, end2);
    }
}

static void makePseudoChromosomes(Flower *flower, Reference *reference,
        int(*cmpFn)(End *, End *)) {
    /*
     * Uses the above functions to construct a set of pairs of ends and then construct the pseudo-chromsomes,
     * but without any pseudo-adjacencies.
     */
    stList *ends = getAttachedStubEnds(flower);
    /* Sort the ends and so construct pairing of attached ends to construct pseudo chromosomes*/
    stList_sort(ends, (int (*)(const void *, const void *))cmpFn);
    makePseudoChromosomesFromPairs(ends, reference);
    stList_destruct(ends);
}

/*
 * Code to make the pseudo chromosomes for a root problem.
 */

static Cap *makeTopLevelPseudoChromosomesP(End *end) {
    /*
     * Gets the leaf cap from a top level attached stub. There can only
     * be one such cap per top level attached stub, currently.
     */
    Cap *cap = NULL, *cap2;
    End_InstanceIterator *iterator = end_getInstanceIterator(end);
    while ((cap2 = end_getNext(iterator)) != NULL) {
        if (cap_getChildNumber(cap2) == 0) {
            assert(cap == NULL);
            cap = cap2;
        }
    }
    end_destructInstanceIterator(iterator);
    assert(cap != NULL);
    return cap;
}

static int makeTopLevelPseudoChromosomes_cmpEnds(End *end1, End *end2) {
    /*
     * Sorts the attached ends according to the sequences they are connected to.
     */
    Cap *cap1 = makeTopLevelPseudoChromosomesP(end1);
    Cap *cap2 = makeTopLevelPseudoChromosomesP(end2);
    Sequence *sequence1 = cap_getSequence(cap1);
    Sequence *sequence2 = cap_getSequence(cap2);
    assert(sequence1 != NULL);
    assert(sequence2 != NULL);
    int32_t i = cactusMisc_nameCompare(sequence_getName(sequence1),
            sequence_getName(sequence2));
    if (i == 0) {
        assert(cap_getSide(cap1) != cap_getSide(cap2));
        return cap_getSide(cap1) ? -1 : 1; //sort 5' to 3' (i.e. with the 5 prime with a lower index to the 3 prime side)
    } else {
        return i;
    }
}

void makeTopLevelPseudoChromosomes(Flower *flower, Reference *reference) {
    makePseudoChromosomes(flower, reference,
            makeTopLevelPseudoChromosomes_cmpEnds);
}

/*
 * Code to make the pseudo chromosomes for a non-root problem.
 */


static stHash *makeIntermediateLevelPseudoChromosomes_cmpEndsP = NULL;
static Flower *makeIntermediateLevelPseudoChromosomes_parentFlower = NULL;

static int makeIntermediateLevelPseudoChromosomes_cmpEnds(End *end1,
        End *end2) {
    /*
     * Sorts the attached ends according to the pairing in the higher level reference.
     */
    End *end3 = flower_getEnd(
            makeIntermediateLevelPseudoChromosomes_parentFlower, end_getName(
                    end1));
    End *end4 = flower_getEnd(
            makeIntermediateLevelPseudoChromosomes_parentFlower, end_getName(
                    end2));
#ifdef BEN_DEBUG
    assert(end3 != NULL);
    assert(end4 != NULL);
    assert(end_getOrientation(end3));
    assert(end_getOrientation(end4));
    assert(end_isAttached(end3) || end_isBlockEnd(end3));
    assert(end_isAttached(end4) || end_isBlockEnd(end4));
#endif

    PseudoAdjacency *pseudoAdjacency1 = stHash_search(
            makeIntermediateLevelPseudoChromosomes_cmpEndsP, end3);
    PseudoAdjacency *pseudoAdjacency2 = stHash_search(
            makeIntermediateLevelPseudoChromosomes_cmpEndsP, end4);
#ifdef BEN_DEBUG
    assert(pseudoAdjacency1 != NULL);
    assert(pseudoAdjacency2 != NULL);
#endif
    PseudoChromosome *pseudoChromosome1 = pseudoAdjacency_getPseudoChromosome(
            pseudoAdjacency1);
    PseudoChromosome *pseudoChromosome2 = pseudoAdjacency_getPseudoChromosome(
            pseudoAdjacency2);

    //Sort first by pseudo chromosome
    int32_t i = cactusMisc_nameCompare(pseudoChromosome_getName(
            pseudoChromosome1), pseudoChromosome_getName(pseudoChromosome2));
    if (i != 0) { //we order, so the chromosomes reflect the ordering of the parent adjacencies on the circle.
        return i;
    }
    assert(pseudoChromosome1 == pseudoChromosome2);
    //Then by pseudo-adjacency
    i = cactusMisc_nameCompare(pseudoAdjacency_getName(pseudoAdjacency1),
            pseudoAdjacency_getName(pseudoAdjacency2));
    if (i != 0) {
        return i;
    }
#ifdef BEN_DEBUG
    assert(pseudoAdjacency1 == pseudoAdjacency2);
    //if in the same pseudo adjacency then
    //return the proper pair according to ordering of ends in pair.
    assert(end3 == pseudoAdjacency_get5End(pseudoAdjacency1) || end3 == pseudoAdjacency_get3End(pseudoAdjacency1));
    assert(end4 == pseudoAdjacency_get5End(pseudoAdjacency1) || end4 == pseudoAdjacency_get3End(pseudoAdjacency1));
#endif
    return pseudoAdjacency_get5End(pseudoAdjacency1) == end3 ? -1 : 1;
}

void makeIntermediateLevelPseudoChromosomes(Flower *flower,
        Reference *reference) {
    Group *parentGroup = flower_getParentGroup(flower);
    assert(parentGroup != NULL);
    Flower *parentFlower = group_getFlower(parentGroup);
    Reference *parentReference = flower_getReference(parentFlower);
    assert(parentReference != NULL);

    makeIntermediateLevelPseudoChromosomes_parentFlower = parentFlower;
    makeIntermediateLevelPseudoChromosomes_cmpEndsP
            = reference_getEndToPseudoAdjacencyHash(parentReference);

    makePseudoChromosomes(flower, reference,
            makeIntermediateLevelPseudoChromosomes_cmpEnds);

    stHash_destruct(makeIntermediateLevelPseudoChromosomes_cmpEndsP);
}
