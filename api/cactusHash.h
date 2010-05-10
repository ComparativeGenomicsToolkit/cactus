#ifndef CACTUS_HASH_H_
#define CACTUS_HASH_H_

/*
 * cactusHash.h
 *
 *  Created on: 4 Apr 2010
 *      Author: benedictpaten
 */

#include "cactusGlobals.h"

/*
 * Constructs hash, with no destructors for keys or values.
 */
Hash *hash_construct();

/*
 * Constructs a hash with given destructors, if null then destructors are ignored
 */
Hash *hash_construct2(void (*destructKeys)(void *), void (*destructValues)(void *));

/*
 * Constructs a hash using the given comparison functions.
 */
Hash *hash_construct3(uint32_t (*hashKey)(void *), int32_t (*hashEqualsKey)(void *, void *),
		void (*destructKeys)(void *), void (*destructValues)(void *));

/*
 * Destructs a hash.
 */
void hash_destruct(Hash *hash);

/*
 * Insert element, overiding if already present.
 */
void hash_insert(Hash *hash, void *key, void *value);

/*
 * Search for value, returns null if not present.
 */
void *hash_search(Hash *hash, void *key);

/*
 * Removes element, returning removed element.
 */
void *hash_remove(Hash *hash, void *key);

/*
 * Returns the number of key/value pairs in the hash.
 */
int32_t hash_size(Hash *hash);

/*
 * Returns an iterator of the keys in the hash.
 */
Hash_Iterator *hash_getIterator(Hash *hash);

/*
 * Gets the next key from the iterator.
 */
void *hash_getNext(Hash_Iterator *iterator);

/*
 * Duplicates the iterator.
 */
Hash_Iterator *hash_copyIterator(Hash_Iterator *iterator);

/*
 * Destructs the iterator.
 */
void hash_destructIterator(Hash_Iterator *iterator);

#endif
