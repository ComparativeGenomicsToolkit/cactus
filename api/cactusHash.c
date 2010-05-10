/*
 * cactusHash.c
 *
 *  Created on: 4 Apr 2010
 *      Author: benedictpaten
 */
#include "cactusGlobalsPrivate.h"

Hash *hash_construct() {
	return hash_construct3(hashtable_key, hashtable_equalKey, NULL, NULL);
}

Hash *hash_construct2(void (*destructKeys)(void *), void (*destructValues)(void *)) {
	return hash_construct3(hashtable_key, hashtable_equalKey, destructKeys, destructValues);
}

Hash *hash_construct3(uint32_t (*hashKey)(void *), int32_t (*hashEqualsKey)(void *, void *),
		void (*destructKeys)(void *), void (*destructValues)(void *)) {
	Hash *hash = malloc(sizeof(Hash));
	hash->hash = create_hashtable(0, hashKey, hashEqualsKey, destructKeys, destructValues);
	hash->destructKeys = destructKeys != NULL;
	hash->destructValues = destructValues != NULL;
	return hash;
}

void hash_destruct(Hash *hash) {
	hashtable_destroy(hash->hash, hash->destructValues, hash->destructKeys);
	free(hash);
}

void hash_insert(Hash *hash, void *key, void *value) {
	hashtable_insert(hash->hash, key, value);
}

void *hash_search(Hash *hash, void *key) {
	return hashtable_search(hash->hash, key);
}

void *hash_remove(Hash *hash, void *key) {
	return hashtable_remove(hash->hash, key, 0);
}

int32_t hash_size(Hash *hash) {
	return hashtable_count(hash->hash);
}

Hash_Iterator *hash_getIterator(Hash *hash) {
	return hashtable_iterator(hash->hash);
}

void *hash_getNext(Hash_Iterator *iterator) {
	if(iterator->e != NULL) {
		void *o = hashtable_iterator_key(iterator);
		hashtable_iterator_advance(iterator);
		return o;
	}
	return NULL;
}

Hash_Iterator *hash_copyIterator(Hash_Iterator *iterator) {
	Hash_Iterator *iterator2 = malloc(sizeof(Hash_Iterator));
	iterator2->h = iterator->h;
	iterator2->e = iterator->e;
	iterator2->parent = iterator->parent;
	iterator2->index = iterator->index;
	return iterator2;
}

void hash_destructIterator(Hash_Iterator *iterator) {
	free(iterator);
}
