/*
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ST_CACTUS_SETUP_H_
#define ST_CACTUS_SETUP_H_

#include "cactus.h"

/*
 * Build the first flower, adding an event tree and sequence files.
 */
Flower *cactus_setup_first_flower(CactusDisk *cactusDisk, CactusParams *params,
                                  char *speciesTree, char *outgroupEvents, char *sequenceFilesAndEvents);

#endif /* ST_CACTUS_SETUP_H_ */
