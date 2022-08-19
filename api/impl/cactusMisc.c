/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"
#include <ctype.h>
#include <stdio.h>

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

