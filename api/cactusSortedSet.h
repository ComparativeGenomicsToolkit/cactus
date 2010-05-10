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
SortedSet *sortedSet_construct(int32_t (*compareFn)(const void *, const void *, void *));

/*
 * Destructs the sorted set, applying the destruct function to each element.
 */
void sortedSet_destruct(SortedSet *sortedSet, void (*destructElementFn)(void *, void *));

/*
 * Inserts the object into the sorted set.
 */
void sortedSet_insert(SortedSet *sortedSet, void *object);

/*
 * Finds the objects in the sorted set, or returns null.
 */
void *sortedSet_find(SortedSet *sortedSet, void *object);

/*
 * Deletes the object in the sorted set.
 */
void sortedSet_delete(SortedSet *sortedSet, void *object);

/*
 * Gets the number of elements in the sorted set.
 */
int32_t sortedSet_getLength(SortedSet *sortedSet);

/*
 * Gets the first element (with lowest value), in the sorted set.
 */
void *sortedSet_getFirst(SortedSet *items);

/*
 * Gets the last element in the sorted set.
 */
void *sortedSet_getLast(SortedSet *items);

/*
 * Constructs an iterator for the sorted set.
 */
SortedSet_Iterator *iterator_construct(SortedSet *items);

/*
 * Destructs an iterator for the sorted set.
 */
void iterator_destruct(SortedSet_Iterator *iterator);

/*
 * Gets next element in the sorted set.
 */
void *iterator_getNext(SortedSet_Iterator *iterator);

/*
 * Gets the previous element in the sorted set.
 */
void *iterator_getPrevious(SortedSet_Iterator *iterator);

/*
 * Copies the iterator.
 */
SortedSet_Iterator *iterator_copy(SortedSet_Iterator *iterator);

#endif
