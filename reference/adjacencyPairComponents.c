#include "sonLib.h"
#include "cactus.h"
#include "adjacencyPairs.h"

static void extendComponent(End *end, stList *component,
		stHash *componentsHash, stHash *adjacencies) {
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
					component, componentsHash, adjacencies);
		}
		if (end_isBlockEnd(end)) {
			extendComponent(end_getOtherBlockEnd(end), component,
					componentsHash, adjacencies);
		}
	}
}

stList *adjacencyHash_getConnectedComponent(stHash *adjacencies,
		Flower *flower, End *end) {
	stList *component = stList_construct();
	stHash *componentsHash = stHash_construct();
	extendComponent(end, component, componentsHash, adjacencies);
	stHash_destruct(componentsHash);
	return component;
}

stList *adjacencyHash_getConnectedComponents(stHash *adjacencies,
		Flower *flower) {
	stList *components =
			stList_construct3(0, (void(*)(void *)) stList_destruct);
	stHash *componentsHash = stHash_construct();
	End *end;
	Flower_EndIterator *endIterator = flower_getEndIterator(flower);
	while ((end = flower_getNextEnd(endIterator)) != NULL) { //iterates over the positive oriented ends.
		if (end_isBlockEnd(end) || end_isAttached(end)) { //ignore free stubs
			stList *component = stHash_search(componentsHash, end);
			if (component == NULL) {
				component = stList_construct();
				stList_append(components, component);
				extendComponent(end, component, componentsHash, adjacencies);
			}
		}
	}
	flower_destructEndIterator(endIterator);
#ifdef BEN_DEBUG
	assert(flower_getEndNumber(flower) - flower_getFreeStubEndNumber(flower) == stHash_size(componentsHash));
#endif
	stHash_destruct(componentsHash);
	return components;
}

void getAttachedStubEndsInComponent(stList *component, End **end1, End **end2) {
	int32_t i;
	*end1 = NULL;
	*end2 = NULL;
	for (i = 0; i < stList_length(component); i++) {
		End *end3 = stList_get(component, i);
		if (end_isStubEnd(end3)) {
			assert(end_isAttached(end3));
			assert(*end2 == NULL);
			if (*end1 == NULL) {
				*end1 = end3;
			} else {
				*end2 = end3;
			}
		}
	}
}
