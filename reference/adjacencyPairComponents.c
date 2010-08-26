#include "sonLib.h"
#include "cactus.h"
#include "adjacencyPairs.h"

static void extendComponent(End *end, stList *component,
        stHash *componentsHash, stHash *adjacencies, stHash *hyperChains) {
    /*
     * Sub-function of get connected components, extends the component.
     */
    assert(end_getOrientation(end));
    if (stHash_search(componentsHash, end) == NULL) {
        stList_append(component, end);
        stHash_insert(componentsHash, end, component);
        AdjacencyPair *adjacencyPair = stHash_search(adjacencies, end);
        if (adjacencyPair != NULL) {
            extendComponent(adjacencyPair_getOtherEnd(adjacencyPair, end),
                    component, componentsHash, adjacencies, hyperChains);
        }
        if (stHash_search(hyperChains, end) != NULL) {
            extendComponent(stHash_search(hyperChains, end), component,
                    componentsHash, adjacencies, hyperChains);
        }
    }
}

stList *adjacencyHash_getConnectedComponent(stHash *adjacencies, stHash *hyperChains, End *end) {
    stList *component = stList_construct();
    stHash *componentsHash = stHash_construct();
    extendComponent(end, component, componentsHash, adjacencies, hyperChains);
    stHash_destruct(componentsHash);
    return component;
}

stList *adjacencyHash_getConnectedComponents(stHash *adjacencies, stHash *hyperChains) {
    stList *components =
            stList_construct3(0, (void(*)(void *)) stList_destruct);
    stHash *componentsHash = stHash_construct();
    stHashIterator *endIterator = stHash_getIterator(adjacencies);
    End *end;
    while ((end = stHash_getNext(endIterator)) != NULL) { //iterates over the positive oriented ends.
        assert(end_getOrientation(end));
        assert(end_isAttached(end)); //ignore free stubs
        assert(!end_isBlockEnd(end)); //We are only in terminal problems
        stList *component = stHash_search(componentsHash, end);
        if (component == NULL) {
            component = stList_construct();
            stList_append(components, component);
            extendComponent(end, component, componentsHash, adjacencies, hyperChains);
        }
    }
    stHash_destructIterator(endIterator);
#ifdef BEN_DEBUG
    assert(stHash_size(adjacencies) == stHash_size(componentsHash));
#endif
    stHash_destruct(componentsHash);
    return components;
}

void getTopStubEndsInComponent(stList *component, stHash *hyperChains, End **end1, End **end2) {
    int32_t i;
    *end1 = NULL;
    *end2 = NULL;
    for (i = 0; i < stList_length(component); i++) {
        End *end3 = stList_get(component, i);
        if (stHash_search(hyperChains, end3) == NULL) {
            assert(*end2 == NULL);
            if (*end1 == NULL) {
                *end1 = end3;
            } else {
                *end2 = end3;
            }
        }
    }
}
