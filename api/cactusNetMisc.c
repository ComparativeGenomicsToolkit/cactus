#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Useful utility functions.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

char netMisc_reverseComplementChar(char c) {
	switch(c) {
		case 'a':
			return 't';
		case 'c':
			return 'g';
		case 'g':
			return 'c';
		case 't':
			return 'a';
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		default:
			return c;
	}
}

char *netMisc_reverseComplementString(const char *string) {
	int32_t i, j;

	j = strlen(string);
	char *cA;

	cA = mallocLocal(sizeof(char) *(j+1));
	for(i=0; i<j; i++) {
		cA[i] = netMisc_reverseComplementChar(string[j-1-i]);
	}
	cA[j] = '\0';
	return cA;
}

int32_t netMisc_nameCompare(Name name1, Name name2) {
	return name1 - name2 > 0 ? 1 : (name1 < name2 ? -1 : 0); //cautious, to avoid 32bit overflows.
}

Name netMisc_stringToName(const char *stringName) {
	assert(stringName != NULL);
	Name name;
	int32_t i = sscanf(stringName, NAME_STRING, &name);
	if(i != 1) {
		fprintf(stderr, "Can not get a valid name from the given string: %s\n", stringName);
		return NULL_NAME;
	}
	return name;
}

char *netMisc_nameToString(Name name) {
	char *cA;
	cA = mallocLocal(sizeof(char)*21);
	sprintf(cA, NAME_STRING, name);
	return cA;
}

const char *netMisc_nameToStringStatic(Name name) {
	static char cA[100];
	sprintf(cA, NAME_STRING, name);
	return cA;
}

char *netMisc_nameToStringWithOrientation(Name name, int32_t orientation) {
	char *cA;
	cA = mallocLocal(sizeof(char)*22);
	sprintf(cA, orientation ? "%s" : "-%s", netMisc_nameToStringStatic(name));
	return cA;
}

const char *netMisc_nameToStringStaticWithOrientiation(Name name, int32_t orientation) {
	static char cA[100];
	sprintf(cA, orientation ? "%s" : "-%s", netMisc_nameToStringStatic(name));
	return cA;
}




