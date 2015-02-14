/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_MISC_H_
#define CACTUS_MISC_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Useful utility functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

extern const char *CACTUS_CHECK_EXCEPTION_ID;

/*
 * Compares to names, giving an ordering to names (arbitrary but consistent).
 */
int64_t cactusMisc_nameCompare(Name name1, Name name2);

/*
 * Converts the string which holds the name (and nothing else), into a name.
 */
Name cactusMisc_stringToName(const char *stringName);

/*
 * Creates a new string (which must be freed) representing the name as a string.
 */
char *cactusMisc_nameToString(Name name);

/*
 * Creates a static string (which needn't be freed) representing the name as a string.
 */
const char *cactusMisc_nameToStringStatic(Name name);

/*
 * Creates a new string with orientation sign (which must be freed) representing the name as a string.
 */
char *cactusMisc_nameToStringWithOrientation(Name name, int64_t orientation);

/*
 * Creates a static string with orientation sign (which needn't be freed) representing the name as a string.
 */
const char *cactusMisc_nameToStringStaticWithOrientiation(Name name, int64_t orientation);

/*
 * Gets the default name of the reference event string.
 */
const char *cactusMisc_getDefaultReferenceEventHeader();

/*
 * Use a bulk get to efficiently precache the nested flowers of a set of parent flowers.
 */
void preCacheNestedFlowers(CactusDisk *cactusDisk, stList *flowers);

/*
 * Check a condition is true, if not throw an exception - a short hand to defining your own exception.
 */
void cactusCheck(bool condition);

/*
 * Check a condition is true, if not throw an exception with the given string.
 */
void cactusCheck2(bool condition, char *string, ...);

#endif
