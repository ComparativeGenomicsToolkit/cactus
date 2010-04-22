#ifndef CACTUS_NET_MISC_H_
#define CACTUS_NET_MISC_H_

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
char netMisc_reverseComplementChar(char c);

/*
 * Computes the reverse complement of the string, returning the r-c string in newly allocated memory that must be freed.
 */
char *netMisc_reverseComplementString(const char *string);

/*
 * Compares to names, giving an ordering to names (arbitrary but consistent).
 */
int32_t netMisc_nameCompare(Name name1, Name name2);

/*
 * Converts the string which holds the name (and nothing else), into a name.
 */
Name netMisc_stringToName(const char *stringName);

/*
 * Creates a new string (which must be freed) representing the name as a string.
 */
char *netMisc_nameToString(Name name);

/*
 * Creates a static string (which needn't be freed) representing the name as a string.
 */
const char *netMisc_nameToStringStatic(Name name);

/*
 * Creates a new string with orientation sign (which must be freed) representing the name as a string.
 */
char *netMisc_nameToStringWithOrientation(Name name, int32_t orientation);

/*
 * Creates a static string with orientation sign (which needn't be freed) representing the name as a string.
 */
const char *netMisc_nameToStringStaticWithOrientiation(Name name, int32_t orientation);

#endif
