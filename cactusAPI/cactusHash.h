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
 * Constructs a hash.
 */
Hash *hash_construct();

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
