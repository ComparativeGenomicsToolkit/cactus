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

const char *cactusMisc_nameToStringStatic(Name name) {
    static char cA[100];
    sprintf(cA, NAME_STRING, name);
    return cA;
}

const char *cactusMisc_getDefaultReferenceEventHeader() {
    static char cA[10];
    sprintf(cA, "reference");
    return cA;
}

void preCacheNestedFlowers(CactusDisk *cactusDisk, stList *flowers) {
    stList *nestedFlowerNames = stList_construct3(0, free);
    for (int64_t i = 0; i < stList_length(flowers); i++) {
        Flower *flower = stList_get(flowers, i);
        Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
        Group *group;
        while ((group = flower_getNextGroup(groupIt)) != NULL) {
            if (!group_isLeaf(group)) {
                int64_t *iA = st_malloc(sizeof(int64_t));
                iA[0] = group_getName(group);
                stList_append(nestedFlowerNames, iA);
            }
        }
        flower_destructGroupIterator(groupIt);
    }
    stList_destruct(cactusDisk_getFlowers(cactusDisk, nestedFlowerNames));
    stList_destruct(nestedFlowerNames);
}

const char *CACTUS_CHECK_EXCEPTION_ID = "CACTUS_CHECK_EXCEPTION_ID";

void cactusCheck(bool condition) {
    if (!condition) {
        assert(0);
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
        assert(0);
        assert(strlen(cA) < 100000);
        stThrowNew(CACTUS_CHECK_EXCEPTION_ID, "Cactus check condition failed: %s", cA);
    }
}

