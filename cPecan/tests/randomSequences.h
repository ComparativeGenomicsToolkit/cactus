/*
 * randomSequences.h
 *
 *  Created on: 3 Mar 2012
 *      Author: benedictpaten
 */

#ifndef RANDOMSEQUENCES_H_
#define RANDOMSEQUENCES_H_

char getRandomChar();

char *getRandomSequence(int64_t length);

char *evolveSequence(const char *startSequence);

#endif /* RANDOMSEQUENCES_H_ */
