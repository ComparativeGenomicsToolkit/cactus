/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"
#include <ctype.h>

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Useful utility functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

static char correctCase(char newChar, char oldChar) {
    return toupper(oldChar) == oldChar ? toupper(newChar) : tolower(newChar);
}

char cactusMisc_reverseComplementChar(char c) {
    switch (c) {
        case 'a':
        case 'A':
            return correctCase('t', c);
        case 'c':
        case 'C':
            return correctCase('g', c);
        case 'g':
        case 'G':
            return correctCase('c', c);
        case 't':
        case 'T':
            return correctCase('a', c);

            //Two redundant characters
        case 'w':
        case 'W':
        case 's':
        case 'S':
            return c; //Complement is the same for S and W

        case 'M':
        case 'm':
            return correctCase('k', c); //Complement is k
        case 'k':
        case 'K':
            return correctCase('m', c); //Complement is m

        case 'r':
        case 'R':
            return correctCase('y', c); //Complement is y
        case 'y':
        case 'Y':
            return correctCase('r', c); //Complement is r

            //Three redundant characters
        case 'b':
        case 'B':
            return correctCase('v', c); //Complement is v
        case 'v':
        case 'V':
            return correctCase('b', c); //Complement is b

        case 'd':
        case 'D':
            return correctCase('h', c); //Complement is h
        case 'h':
        case 'H':
            return correctCase('d', c); //Complement is d

        default: //Includes N, but also any other character.
            return c;
    }
}

char *cactusMisc_reverseComplementString(const char *string) {
    int32_t i, j;

    j = strlen(string);
    char *cA;

    cA = st_malloc(sizeof(char) * (j + 1));
    for (i = 0; i < j; i++) {
        cA[i] = cactusMisc_reverseComplementChar(string[j - 1 - i]);
    }
    cA[j] = '\0';
    return cA;
}

int32_t cactusMisc_nameCompare(Name name1, Name name2) {
    return name1 > name2 ? 1 : (name1 < name2 ? -1 : 0);
}

Name cactusMisc_stringToName(const char *stringName) {
    assert(stringName != NULL);
    Name name;
    int32_t i = sscanf(stringName, NAME_STRING, &name);
    if (i != 1) {
        fprintf(stderr, "Can not get a valid name from the given string: %s\n",
                stringName);
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

stList *parseFlowers(char **flowerNames, int32_t flowerNamesLength, CactusDisk *cactusDisk) {
    stList *flowerNamesList = stList_construct3(0, free);
    assert(flowerNamesLength % 2 == 0);
    for(int32_t i=0; i<flowerNamesLength; i+=2) {
        Name firstFlowerName = cactusMisc_stringToName(flowerNames[i]);
        int32_t flowerNumber;
        int32_t j = sscanf(flowerNames[i+1], "%i", &flowerNumber);
        assert(j == 1);
        assert(flowerNumber >= 1);
        for(j=0; j<flowerNumber; j++) {
            int64_t *iA = st_malloc(sizeof(int64_t));
            iA[0] = firstFlowerName + j;
            //cactusDisk_getFlower(cactusDisk, iA[0]);
            stList_append(flowerNamesList, iA);
        }
    }
    stList *flowers = cactusDisk_getFlowers(cactusDisk, flowerNamesList);
    stList_destruct(flowerNamesList);
    return flowers;
    /*stList *flowers = stList_construct();
    assert(flowerNamesLength % 2 == 0);
    for(int32_t i=0; i<flowerNamesLength; i+=2) {
        Name firstFlowerName = cactusMisc_stringToName(flowerNames[i]);
        int32_t flowerNumber;
        int32_t j = sscanf(flowerNames[i+1], "%i", &flowerNumber);
        assert(j == 1);
        for(j=0; j<flowerNumber; j++) {
            Flower *flower = cactusDisk_getFlower(cactusDisk, firstFlowerName + j);
            assert(flower != NULL);
            stList_append(flowers, flower);
        }
    }
    return flowers;*/
}

