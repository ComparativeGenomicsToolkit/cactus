/*
 * cactusHash.c
 *
 *  Created on: 4 Apr 2010
 *      Author: benedictpaten
 */
#include "cactusGlobalsPrivate.h"

Hash *hash_construct() {
	return create_hashtable(1, hashtable_key, hashtable_equalKey, NULL, NULL);
}

void hash_destruct(Hash *hash) {
	hashtable_destroy(hash, 0, 0);
}

void hash_insert(Hash *hash, void *key, void *value) {
	hashtable_insert(hash, key, value);
}

void *hash_search(Hash *hash, void *key) {
	return hashtable_search(hash, key);
}

void *hash_remove(Hash *hash, void *key) {
	return hashtable_remove(hash, key, 0);
}
