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

/*
 * Computes the reverse complement character of a ACTGactg, returning other characters unmodified.
 */
char cactusMisc_reverseComplementChar(char c);

/*
 * Computes the reverse complement of the string, returning the r-c string in newly allocated memory that must be freed.
 */
char *cactusMisc_reverseComplementString(const char *string);

/*
 * Compares to names, giving an ordering to names (arbitrary but consistent).
 */
int32_t cactusMisc_nameCompare(Name name1, Name name2);

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
char *cactusMisc_nameToStringWithOrientation(Name name, int32_t orientation);

/*
 * Creates a static string with orientation sign (which needn't be freed) representing the name as a string.
 */
const char *cactusMisc_nameToStringStaticWithOrientiation(Name name, int32_t orientation);

/*
 * Gets the default name of the reference event string.
 */
const char *cactusMisc_getDefaultReferenceEventHeader();

/*
 * Encodes a list of flower names and prints the encoding to the file.
 */
void cactusMisc_encodeFlowersString(stList *flowerNames, FILE *fileHandle);

/*
 * Decodes a list of flower names and returns them from the filehandle.
 */
stList *cactusMisc_decodeFlowersString(CactusDisk *cactusDisk, FILE *fileHandle);

/*
 *  As parseFlowers, but reads the list from stdin.
 */
stList *cactusMisc_parseFlowersFromStdin(CactusDisk *cactusDisk);

/*
 * Use a bulk get to efficiently precache the nested flowers of a set of parent flowers.
 */
void preCacheNestedFlowers(CactusDisk *cactusDisk, stList *flowers);

#endif
