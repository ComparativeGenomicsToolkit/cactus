#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions on a sorted set and its iterator
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

struct avl_table *sortedSet_construct(int32_t (*compareFn)(const void *, const void *, void *)) {
	return avl_create(compareFn, NULL, NULL);
}

void sortedSet_destruct(struct avl_table *sortedSet, void (*destructElementFn)(void *, void *)) {
	avl_destroy(sortedSet, destructElementFn);
}

void sortedSet_insert(struct avl_table *sortedSet, void *object) {
	avl_insert(sortedSet, object);
}

void *sortedSet_find(struct avl_table *sortedSet, void *object) {
	return avl_find(sortedSet, object);
}

void sortedSet_delete(struct avl_table *sortedSet, void *object) {
	avl_delete(sortedSet, object);
}

int32_t sortedSet_getLength(struct avl_table *sortedSet) {
	return avl_count(sortedSet);
}

void *sortedSet_getFirst(struct avl_table *items) {
	static struct avl_traverser iterator;
	avl_t_init(&iterator, items);
	return avl_t_first(&iterator, items);
}

void *sortedSet_getLast(struct avl_table *items) {
	static struct avl_traverser iterator;
	avl_t_init(&iterator, items);
	return avl_t_last(&iterator, items);
}

struct avl_traverser *iterator_construct(struct avl_table *items) {
	struct avl_traverser *iterator;
	iterator = mallocLocal(sizeof(struct avl_traverser));
	avl_t_init(iterator, items);
	return iterator;
}

void iterator_destruct(struct avl_traverser *iterator) {
	free(iterator);
}

void *iterator_getNext(struct avl_traverser *iterator) {
	return avl_t_next(iterator);
}

struct avl_traverser *iterator_copy(struct avl_traverser *iterator) {
	struct avl_traverser *copyIterator;
	copyIterator = mallocLocal(sizeof(struct avl_traverser));
	avl_t_copy(copyIterator, iterator);
	return copyIterator;
}

void *iterator_getPrevious(struct avl_traverser *iterator) {
	return avl_t_prev(iterator);
}
