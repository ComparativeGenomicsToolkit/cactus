#ifndef CACTUS_SORTED_SET_H_
#define CACTUS_SORTED_SET_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Sorted set functions
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Constructs a sorted set, using the given comparison function.
 */
struct avl_table *sortedSet_construct(int32_t (*compareFn)(const void *, const void *, void *));

/*
 * Destructs the sorted set, applying the destruct function to each element.
 */
void sortedSet_destruct(struct avl_table *sortedSet, void (*destructElementFn)(void *, void *));

/*
 * Inserts the object into the sorted set.
 */
void sortedSet_insert(struct avl_table *sortedSet, void *object);

/*
 * Finds the objects in the sorted set, or returns null.
 */
void *sortedSet_find(struct avl_table *sortedSet, void *object);

/*
 * Deletes the object in the sorted set.
 */
void sortedSet_delete(struct avl_table *sortedSet, void *object);

/*
 * Gets the number of elements in the sorted set.
 */
int32_t sortedSet_getLength(struct avl_table *sortedSet);

/*
 * Gets the first element (with lowest value), in the sorted set.
 */
void *sortedSet_getFirst(struct avl_table *items);

/*
 * Constructs an iterator for the sorted set.
 */
struct avl_traverser *iterator_construct(struct avl_table *items);

/*
 * Destructs an iterator for the sorted set.
 */
void iterator_destruct(struct avl_traverser *iterator);

/*
 * Gets next element in the sorted set.
 */
void *iterator_getNext(struct avl_traverser *iterator);

/*
 * Gets the previous element in the sorted set.
 */
void *iterator_getPrevious(struct avl_traverser *iterator);

/*
 * Copies the iterator.
 */
struct avl_traverser *iterator_copy(struct avl_traverser *iterator);

#endif
