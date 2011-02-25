/*
 * mafs.h
 *
 *  Created on: 24 Feb 2011
 *      Author: benedictpaten
 */

#ifndef MAFS_H_
#define MAFS_H_

void getMAFBlock(Block *block, FILE *fileHandle);

void getMAFsReferenceOrdered(Flower *flower, FILE *fileHandle, void (*getMafBlock)(Block *, FILE *));

void getMAFs(Flower *flower, FILE *fileHandle, void (*getMafBlock)(Block *, FILE *));

void makeMAFHeader(Flower *flower, FILE *fileHandle);

#endif /* MAFS_H_ */
