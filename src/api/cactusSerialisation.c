/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "cactusGlobalsPrivate.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions for serialising the objects.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void binaryRepresentation_writeElementType(char elementCode, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	writeFn(&elementCode, sizeof(char), 1);
}

void binaryRepresentation_writeString(const char *name, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	int64_t i = strlen(name);
	writeFn(&i, sizeof(int64_t), 1);
	writeFn(name, sizeof(char), i);
}

void binaryRepresentation_writeInteger(int64_t i, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	writeFn(&i, sizeof(int64_t), 1);
}

void binaryRepresentation_writeName(Name name, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	binaryRepresentation_writeInteger(name, writeFn);
}

void binaryRepresentation_writeFloat(float f, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	writeFn(&f, sizeof(float), 1);
}

void binaryRepresentation_writeBool(bool i, void (*writeFn)(const void * ptr, size_t size, size_t count)) {
	writeFn(&i, sizeof(bool), 1);
}

char binaryRepresentation_peekNextElementType(void *binaryString) {
	return *((char *)binaryString);
}

char binaryRepresentation_popNextElementType(void **binaryString) {
	char *c;
	c = *binaryString;
	*binaryString = c + 1;
	return *c;
}

char *binaryRepresentation_getString(void **binaryString) {
	int64_t i;
	char *cA;
	i = binaryRepresentation_getInteger(binaryString);
	cA = st_malloc(sizeof(char)*(i+1));
	memcpy(cA, *binaryString, sizeof(char)*i);
	cA[i] = '\0';
	*binaryString = *((char **)binaryString) + i;
	return cA;
}

char *binaryRepresentation_getStringStatic_cA = NULL;
const char *binaryRepresentation_getStringStatic(void **binaryString) {
	if(binaryRepresentation_getStringStatic_cA != NULL) {
		free(binaryRepresentation_getStringStatic_cA);
	}
	binaryRepresentation_getStringStatic_cA = binaryRepresentation_getString(binaryString);
	return binaryRepresentation_getStringStatic_cA;
}

int64_t binaryRepresentation_getInteger(void **binaryString) {
	int64_t *i;
	i = *binaryString;
	*binaryString = i + 1;
	return *i;
}

Name binaryRepresentation_getName(void **binaryString) {
	return binaryRepresentation_getInteger(binaryString);
}

float binaryRepresentation_getFloat(void **binaryString) {
	float *i;
	i = *binaryString;
	*binaryString = i + 1;
	return *i;
}

bool binaryRepresentation_getBool(void **binaryString) {
	bool *i;
	i = *binaryString;
	*binaryString = i + 1;
	return *i;
}

int64_t binaryRepresentation_makeBinaryRepresentationP_i = 0;
void binaryRepresentation_makeBinaryRepresentationP(const void * ptr, size_t size, size_t count) {
	/*
	 * Records the cummulative size of the substrings written out in creating the flower.
	 */
	assert(ptr != NULL);
	binaryRepresentation_makeBinaryRepresentationP_i += size * count;
}

char *binaryRepresentation_makeBinaryRepresentationP2_vA = NULL;
void binaryRepresentation_makeBinaryRepresentationP2(const void * ptr, size_t size, size_t count) {
	/*
	 * Cummulates all the binary data into one array
	 */
	memcpy(binaryRepresentation_makeBinaryRepresentationP2_vA, ptr, size*count);
	binaryRepresentation_makeBinaryRepresentationP2_vA += size * count;
}

void *binaryRepresentation_makeBinaryRepresentation(void *object, void (*writeBinaryRepresentation)(void *, void (*writeFn)(const void * ptr, size_t size, size_t count)), int64_t *recordSize) {
	void *vA;
	binaryRepresentation_makeBinaryRepresentationP_i = 0;
	writeBinaryRepresentation(object, binaryRepresentation_makeBinaryRepresentationP);
	assert(binaryRepresentation_makeBinaryRepresentationP_i < INT64_MAX);
	vA = st_malloc(binaryRepresentation_makeBinaryRepresentationP_i);
	binaryRepresentation_makeBinaryRepresentationP2_vA = vA;
	writeBinaryRepresentation(object, binaryRepresentation_makeBinaryRepresentationP2);
	*recordSize = binaryRepresentation_makeBinaryRepresentationP_i;
	return vA;
}

void *binaryRepresentation_resizeObjectAsPowerOf2(void *vA, int64_t *recordSize) {
    if(*recordSize == 0) {
        *recordSize = 1;
    }
    int64_t finalSize = pow(2, log(*recordSize * 2)/log(2.0));
    assert(finalSize >= *recordSize);
    vA = realloc(vA, finalSize);
    if(vA == NULL) {
        st_errAbort("Could not realloc memory\n");
    }
    *recordSize = finalSize;
    return vA;
}
