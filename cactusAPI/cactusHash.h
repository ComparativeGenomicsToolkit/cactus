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

#endif
