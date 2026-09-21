/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <inttypes.h>

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Useful utility functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

int64_t cactusMisc_nameCompare(Name name1, Name name2) {
    return name1 > name2 ? 1 : (name1 < name2 ? -1 : 0);
}

Name cactusMisc_stringToName(const char *stringName) {
    assert(stringName != NULL);
    Name name;
    int64_t i = sscanf(stringName, NAME_STRING, &name);
    if (i != 1) {
        fprintf(stderr, "Can not get a valid name from the given string: %s\n", stringName);
        return NULL_NAME;
    }
    return name;
}

char *cactusMisc_nameToString(Name name) {
    char *cA;
    cA = st_malloc(sizeof(char) * 21);
    sprintf(cA, NAME_STRING, name);
    return cA;
}

const char *cactusMisc_getDefaultReferenceEventHeader() {
    return stString_print("reference");
}

const char *CACTUS_CHECK_EXCEPTION_ID = "CACTUS_CHECK_EXCEPTION_ID";

void cactusCheck(bool condition) {
    if (!condition) {
        //assert(0);
        stThrowNew(CACTUS_CHECK_EXCEPTION_ID, "Cactus check condition failed");
    }
}

void cactusCheck2(bool condition, char *string, ...) {
    if(!condition) {
        static char cA[100000];
        va_list ap;
        va_start(ap, string);
        vsprintf(cA, string, ap);
        va_end(ap);
        //assert(0);
        assert(strlen(cA) < 100000);
        stThrowNew(CACTUS_CHECK_EXCEPTION_ID, "Cactus check condition failed: %s", cA);
    }
}

/*
 * Retention guard.  Retention (dirty and muzzy decay both -1) never returns a page, so rss
 * follows cumulative churn rather than the working set and does not saturate: on salamander
 * Anc3 it reached 540 GiB before bar had allocated anything, all of it caf's, against an
 * 809 GiB request.  This is the way out.  The limit is a fraction of the job's request, set
 * by the workflow.  Once tripped it stays tripped -- there is nothing to switch back on.
 */
extern int mallctl(const char *, void *, size_t *, void *, size_t) __attribute__((weak));

static volatile int retentionGuardTripped = 0;
static volatile time_t retentionGuardLastCheck = 0;
static int64_t retentionGuardLimitMb = -1;   // -1 unread, 0 disabled

static int64_t retentionGuardRssMb(void) {
    FILE *f = fopen("/proc/self/statm", "r");
    if (f == NULL) return -1;
    long size = 0, resident = 0;
    if (fscanf(f, "%ld %ld", &size, &resident) != 2) resident = 0;
    fclose(f);
    return (int64_t)resident * (sysconf(_SC_PAGESIZE) / 1024) / 1024;
}

void cactus_retentionGuard(void) {
    if (retentionGuardTripped || retentionGuardLimitMb == 0 || mallctl == NULL) {
        return;
    }
    if (retentionGuardLimitMb < 0) {
        const char *env = getenv("CACTUS_RETENTION_OFF_MB");
        retentionGuardLimitMb = env != NULL ? atoll(env) : 0;
        if (retentionGuardLimitMb <= 0) {
            retentionGuardLimitMb = 0;
            return;
        }
    }
    time_t now = time(NULL);
    if (now == retentionGuardLastCheck) {   // at most once a second, racy on purpose
        return;
    }
    retentionGuardLastCheck = now;

    int64_t rss = retentionGuardRssMb();
    if (rss < retentionGuardLimitMb) {
        return;
    }
    retentionGuardTripped = 1;

    ssize_t stock[2] = { 10000, 0 };
    const char *names[2] = { "dirty_decay_ms", "muzzy_decay_ms" };
    unsigned narenas = 0;
    size_t nsz = sizeof(narenas);
    mallctl("arenas.narenas", &narenas, &nsz, NULL, 0);
    for (int w = 0; w < 2; w++) {
        char key[64];
        snprintf(key, sizeof(key), "arenas.%s", names[w]);
        mallctl(key, NULL, NULL, &stock[w], sizeof(stock[w]));
        for (unsigned i = 0; i < narenas; i++) {
            snprintf(key, sizeof(key), "arena.%u.%s", i, names[w]);
            mallctl(key, NULL, NULL, &stock[w], sizeof(stock[w]));
        }
    }
    for (unsigned i = 0; i < narenas; i++) {
        char key[64];
        snprintf(key, sizeof(key), "arena.%u.purge", i);
        mallctl(key, NULL, NULL, NULL, 0);
    }
    // GiB because every other memory figure this run prints is GiB, and the reader of this line is
    // being asked to change something: it is the one message the guard ever emits.
    st_logCritical("memory guard: rss reached %.1f GiB against a %.1f GiB limit, so jemalloc page "
                   "retention is now off for the rest of this run (%.1f GiB after; the difference "
                   "was pages jemalloc had not handed back, not live data).  Expect the rest of "
                   "the run to be 2-3x slower.  To avoid it next time, give the job more memory, "
                   "or set --consRetainPages 0 to run without retention from the start.\n",
                   rss / 1024.0, retentionGuardLimitMb / 1024.0, retentionGuardRssMb() / 1024.0);
}
