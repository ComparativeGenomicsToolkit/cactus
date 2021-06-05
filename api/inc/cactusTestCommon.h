/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_TEST_COMMON_H_
#define CACTUS_TEST_COMMON_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions shared by the test code.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Get a temporary directory for a test.
 */
char *testCommon_getTmpTestDir(const char *testName);

/*
 * Adds a thread with random nucleotides to the flower, and return its corresponding name in the pinch graph.
 */
 Name testCommon_addThreadToFlower(Flower *flower, char *header, int64_t length);

#endif
