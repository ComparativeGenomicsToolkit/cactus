/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * treeStats.h
 *
 *  Created on: 15-Apr-2010
 *      Author: benedictpaten
 */

#ifndef TREESTATS_H_
#define TREESTATS_H_

/*
 * Writes a lot of stats.
 */
void reportCactusDiskStats(char *cactusDiskName, Flower *flower, FILE *fileHandle, bool perColumnStats, stSortedSet *includeSpecies, stSortedSet *excludeSpecies);

/*
 * Writes stats about blocks, including only those blocks for which include block is true.
 */
void reportBlockStatsP(Flower *flower, FILE *fileHandle, bool(*includeBlock)(
        Block *), const char *attribString, bool perColumnStats);

#endif /* TREESTATS_H_ */
