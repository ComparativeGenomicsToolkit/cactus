/*
 * cactusHashPrivate.h
 *
 *  Created on: 05-Apr-2010
 *      Author: benedictpaten
 */

#ifndef CACTUS_HASH_PRIVATE_H_
#define CACTUS_HASH_PRIVATE_H_

#include "cactusGlobals.h"

struct _cactusHash {
	struct hashtable *hash;
	bool destructKeys, destructValues;
};

#endif /* CACTUSHASHPRIVATE_H_ */
