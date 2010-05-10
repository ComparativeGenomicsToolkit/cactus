#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions on a sorted set and its iterator
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

SortedSet *sortedSet_construct(int32_t (*compareFn)(const void *, const void *, void *)) {
	return avl_create(compareFn, NULL, NULL);
}

void sortedSet_destruct(SortedSet *sortedSet, void (*destructElementFn)(void *, void *)) {
	avl_destroy(sortedSet, destructElementFn);
}

void sortedSet_insert(SortedSet *sortedSet, void *object) {
	avl_insert(sortedSet, object);
}

void *sortedSet_find(SortedSet *sortedSet, void *object) {
	return avl_find(sortedSet, object);
}

void sortedSet_delete(SortedSet *sortedSet, void *object) {
	avl_delete(sortedSet, object);
}

int32_t sortedSet_getLength(SortedSet *sortedSet) {
	return avl_count(sortedSet);
}

void *sortedSet_getFirst(SortedSet *items) {
	static SortedSet_Iterator iterator;
	avl_t_init(&iterator, items);
	return avl_t_first(&iterator, items);
}

void *sortedSet_getLast(SortedSet *items) {
	static SortedSet_Iterator iterator;
	avl_t_init(&iterator, items);
	return avl_t_last(&iterator, items);
}

SortedSet_Iterator *iterator_construct(SortedSet *items) {
	SortedSet_Iterator *iterator;
	iterator = mallocLocal(sizeof(SortedSet_Iterator));
	avl_t_init(iterator, items);
	return iterator;
}

void iterator_destruct(SortedSet_Iterator *iterator) {
	free(iterator);
}

void *iterator_getNext(SortedSet_Iterator *iterator) {
	return avl_t_next(iterator);
}

SortedSet_Iterator *iterator_copy(SortedSet_Iterator *iterator) {
	SortedSet_Iterator *copyIterator;
	copyIterator = mallocLocal(sizeof(SortedSet_Iterator));
	avl_t_copy(copyIterator, iterator);
	return copyIterator;
}

void *iterator_getPrevious(SortedSet_Iterator *iterator) {
	return avl_t_prev(iterator);
}
