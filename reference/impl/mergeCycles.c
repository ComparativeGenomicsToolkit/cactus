#include "sonLib.h"
#include "cactus.h"
#include "adjacencyPairs.h"
#include "adjacencyPairsHash.h"
#include "adjacencyPairComponents.h"
#include "adjacencySwitch.h"

static void splitIntoContigsAndCycles(stList *components, stHash *hyperChains, stList **contigs,
        stList **cycles) {
    /*
     * Divides the list of components into rings and components terminated by a pair of stubs.
     */
    End *end1, *end2;
    *cycles = stList_construct();
    *contigs = stList_construct();
    stList *component;
    for (int32_t i = 0; i < stList_length(components); i++) {
        component = stList_get(components, i);
        getTopStubEndsInComponent(component, hyperChains, &end1, &end2);
        if (end1 != NULL) {
            assert(end2 != NULL);
            stList_append(*contigs, component);
        } else {
            assert(end2 == NULL);
            stList_append(*cycles, component);
        }
    }
}

void mergeCycles(stHash *adjacencies, stHash *hyperChains) {
    stList *contigs, *cycles;

    //Get the comopnents to merge.
    stList *components = adjacencyHash_getConnectedComponents(adjacencies, hyperChains);

    //Get the cycles and contigs
    splitIntoContigsAndCycles(components, hyperChains, &contigs, &cycles);

    //Iterate on the cycles merging them until we have none.
    while (stList_length(cycles) > 0) { //Now merge the cycles
        AdjacencySwitch *adjacencySwitch = NULL;
        stList *mergingComponent = NULL;
        stList *cycle = stList_pop(cycles);
        for (int32_t j = 0; j < stList_length(components); j++) {
            stList *component = stList_get(components, j);
            if (cycle != component) {
                AdjacencySwitch *adjacencySwitch2 =
                        adjacencySwitch_getStrongestAdjacencySwitch(cycle, component, adjacencies);
                if (adjacencySwitch2 != NULL) {
                    if(adjacencySwitch == NULL) {
                        adjacencySwitch = adjacencySwitch2;
                        mergingComponent = component;
                    } else {
                        if (adjacencySwitch_compareStrengthAndPseudoAdjacencies(
                                adjacencySwitch2, adjacencySwitch) > 0) {
                            adjacencySwitch_destruct(adjacencySwitch);
                            adjacencySwitch = adjacencySwitch2;
                            mergingComponent = component; 
                        } else {
                            adjacencySwitch_destruct(adjacencySwitch2);
                        }
                    }
                }
            }
        }
        assert(adjacencySwitch != NULL);

        //Do the switch
        adjacencySwitch_switch(adjacencySwitch, adjacencies);

        //Find the new component..
        stList *newComponent = adjacencyHash_getConnectedComponent(adjacencies,
                hyperChains, adjacencyPair_getEnd1(
                        adjacencySwitch_getAdjacencyPair1(adjacencySwitch)));
        assert(stList_length(cycle) + stList_length(mergingComponent) == stList_length(newComponent));

        //Update the list of components
        stList_removeItem(components, cycle);
        stList_removeItem(components, mergingComponent);
        stList_append(components, newComponent);
        //Update the list of cycles.
        if (stList_contains(cycles, mergingComponent)) {
            stList_removeItem(cycles, mergingComponent);
            stList_append(cycles, newComponent);
        }
        //Destroy the old components.
        stList_destruct(cycle);
        stList_destruct(mergingComponent);

        //Destroy the switch and adjacency pairs
        adjacencyPair_destruct(adjacencySwitch_getAdjacencyPair1(
                adjacencySwitch));
        adjacencyPair_destruct(adjacencySwitch_getAdjacencyPair2(
                adjacencySwitch));
        adjacencySwitch_destruct(adjacencySwitch);

    }
    stList_destruct(contigs);
    stList_destruct(cycles);
    stList_destruct(components);

#ifdef BEN_DEBUG
    components = adjacencyHash_getConnectedComponents(adjacencies, hyperChains);
    splitIntoContigsAndCycles(components, hyperChains, &contigs, &cycles);
    assert(stList_length(cycles) == 0);
    assert(stList_length(contigs) == stList_length(components));
    //assert(stList_length(contigs) > 0);
    stList_destruct(contigs);
    stList_destruct(cycles);
    stList_destruct(components);
#endif
}
