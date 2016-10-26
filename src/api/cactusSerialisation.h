/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef CACTUS_SERIALISATION_H_
#define CACTUS_SERIALISATION_H_

#include "cactusGlobals.h"

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////
//Functions for serialising the objects.
////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

/*
 * Codes used to define which objects are being encoded/decoded from a binary database stream.
 */

#define CODE_META_EVENT 1
#define CODE_EVENT 2
#define CODE_POINTER 3
#define CODE_EVENT_TREE 4
#define CODE_META_SEQUENCE 5
#define CODE_SEQUENCE 6
#define CODE_ADJACENCY 7
#define CODE_PARENT 8
#define CODE_CAP 9
#define CODE_CAP_WITH_COORDINATES 10
#define CODE_CAP_WITH_COORDINATES_BUT_NO_SEQUENCE 11
#define CODE_END_WITHOUT_PHYLOGENY 12
#define CODE_END_WITH_PHYLOGENY 13
#define CODE_SEGMENT 14
#define CODE_BLOCK 15
#define CODE_GROUP 16
#define CODE_GROUP_END 17
#define CODE_LINK 18
#define CODE_CHAIN 19
#define CODE_FACE 20
#define CODE_FLOWER 21
#define CODE_REFERENCE 22
#define CODE_PSEUDO_CHROMOSOME 23
#define CODE_PSEUDO_ADJACENCY 24
#define CODE_CACTUS_DISK 25

/*
 * Writes a code for the element type.
 */
void binaryRepresentation_writeElementType(char elementCode, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Writes a string to the binary stream.
 */
void binaryRepresentation_writeString(const char *string, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Writes an integer to the binary stream
 */
void binaryRepresentation_writeInteger(int64_t i, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Writes an name to the binary stream
 */
void binaryRepresentation_writeName(Name name, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Writes a float to the binary stream.
 */
void binaryRepresentation_writeFloat(float f, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Writes an bool to the binary stream
 */
void binaryRepresentation_writeBool(bool i, void (*writeFn)(const void * ptr, size_t size, size_t count));

/*
 * Returns indicating which element is next, but does not increment the string pointer.
 */
char binaryRepresentation_peekNextElementType(void *binaryString);

/*
 * Returns indicating which element is next, while incrementing the string pointer.
 */
char binaryRepresentation_popNextElementType(void **binaryString);

/*
 * Parses out a string, returning it in a newly allocated string which must be freed.
 */
char *binaryRepresentation_getString(void **binaryString);

/*
 * Parses out a string, placing the memory in a buffer owned by the function. Thid buffer
 * will be overidden by the next call to the function.
 */
const char *binaryRepresentation_getStringStatic(void **binaryString);

/*
 * Parses an integer from binary string.
 */
int64_t binaryRepresentation_getInteger(void **binaryString);

/*
 * Parses a name from a binary string.
 */
Name binaryRepresentation_getName(void **binaryString);

/*
 * Parses a float from the binary string.
 */
float binaryRepresentation_getFloat(void **binaryString);

/*
 * Parses a bool from a binary string.
 */
bool binaryRepresentation_getBool(void **binaryString);

/*
 * Makes a binary representation of an object, using a passed function which writes
 * out the representation of the considered object.
 */
void *binaryRepresentation_makeBinaryRepresentation(void *object, void (*writeBinaryRepresentation)(void *, void (*writeFn)(const void * ptr, size_t size, size_t count)), int64_t *recordSize);

/*
 * Resizes a record as a power of 2.
 */
void *binaryRepresentation_resizeObjectAsPowerOf2(void *vA, int64_t *recordSize);


#endif
