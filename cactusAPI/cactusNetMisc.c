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

	cA = malloc(sizeof(char) *(j+1));
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
	cA = malloc(sizeof(char)*21);
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
	cA = malloc(sizeof(char)*22);
	sprintf(cA, orientation ? "%s" : "-%s", netMisc_nameToStringStatic(name));
	return cA;
}

const char *netMisc_nameToStringStaticWithOrientiation(Name name, int32_t orientation) {
	static char cA[100];
	sprintf(cA, orientation ? "%s" : "-%s", netMisc_nameToStringStatic(name));
	return cA;
}

static int addAdjacenciesToLeafCapsP(Cap **cap1, Cap **cap2) {
	assert(cap_getStrand(*cap1) && cap_getStrand(*cap2));
	Sequence *sequence1 = cap_getSequence(*cap1);
	Sequence *sequence2 = cap_getSequence(*cap2);
	int32_t i = netMisc_nameCompare(sequence_getName(sequence1), sequence_getName(sequence2));
	if(i == 0) {
		int32_t j = cap_getCoordinate(*cap1);
		int32_t k = cap_getCoordinate(*cap2);
		i = j - k;
		if(i == 0) {
			assert(cap_getSegment(*cap1) == cap_getSegment(*cap2));
			j = cap_getSide(*cap1);
			k = cap_getSide(*cap2);
			assert((j && !k) || (!j && k));
			i = j ? -1 : 1;
		}
	}
	return i;
}

void netMisc_addAdjacenciesToLeafCaps(Net *net) {
	End *end;
	Cap *cap;
	Cap *cap2;
	Net_EndIterator *endIterator;
	End_InstanceIterator *instanceIterator;
	struct List *list;
	int32_t i;

	/*
	 * Build a list of caps, then sort them.
	 */
	list = constructEmptyList(0, NULL);
	endIterator = net_getEndIterator(net);
	while ((end = net_getNextEnd(endIterator)) != NULL) {
		instanceIterator = end_getInstanceIterator(end);
		while ((cap = end_getNext(instanceIterator)) != NULL) {
			if(!cap_isInternal(cap)) {
				assert(!cap_isAugmented(cap));
				if(!cap_getStrand(cap)) {
					cap = cap_getReverse(cap);
				}
				listAppend(list, cap);
			}
		}
		end_destructInstanceIterator(instanceIterator);
	}
	net_destructEndIterator(endIterator);
	assert((list->length % 2) == 0);

	/*
	 * Sort the caps.
	 */
	qsort(list->list, list->length, sizeof(void *),
			(int (*)(const void *v, const void *))addAdjacenciesToLeafCapsP);

	/*
	 * Now make the adjacencies.
	 */
	for(i=1; i<list->length; i+=2) {
		cap = list->list[i-1];
		cap2 = list->list[i];
		cap_makeAdjacent1(cap, cap2);
	}

	/*
	 * Clean up.
	 */
	destructList(list);
}



