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
 * Creates a new string with orientation sign (which must be freed) representing the name as a string.
 */
char *cactusMisc_nameToStringWithOrientation(Name name, int64_t orientation);

/*
 * Gets the default name of the reference event string.
 */
const char *cactusMisc_getDefaultReferenceEventHeader();

/*
 * Check a condition is true, if not throw an exception - a short hand to defining your own exception.
 */
void cactusCheck(bool condition);

/*
 * Check a condition is true, if not throw an exception with the given string.
 */
void cactusCheck2(bool condition, char *string, ...);

/*
 * With jemalloc page retention on, hand every retained page back and run in stock decay for
 * the rest of the process once rss reaches CACTUS_RETENTION_OFF_MB.  Cheap enough for inner
 * loops: one clock read, /proc at most once a second.  A no-op without jemalloc, without the
 * variable, or once it has fired.  Called from caf's melting rounds, once per flower in bar
 * (every base aligner) and once per window in the poa engines.
 */
void cactus_retentionGuard(void);

/*
 * True when this machine uses strict overcommit accounting (vm.overcommit_memory = 2).  jemalloc
 * reads the same setting when it starts, and under it gives freed memory back by remapping it
 * PROT_NONE, which splits the process's mappings faster than the kernel merges them.  At
 * vm.max_map_count (65530 by default) every allocation that needs a new mapping fails, with RAM to
 * spare (https://github.com/jemalloc/jemalloc/issues/1328).  Measured on a 5-genome snake ancestor
 * (7.7M blocks after melting): mappings climbed ~1500/s through stCaf_finish to the limit, and
 * the job failed with most of the node's memory free.  Retaining pages never gives memory back, so
 * it never remaps: the same job then ran to the end with at most 1413 mappings.  Read once and
 * cached; false off Linux.
 */
bool cactus_strictOvercommit(void);

#endif
