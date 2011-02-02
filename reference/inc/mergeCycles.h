/*
 * Copyright (C) 2006-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

/*
 * mergeCycles.h
 *
 *  Created on: 24 Aug 2010
 *      Author: benedictpaten
 */

#ifndef MERGECYCLES_H_
#define MERGECYCLES_H_

#include "sonLib.h"
#include "cactus.h"

/*
 * Merges cycles not containing an attached stub end into components containing attached stub ends,
 * Updates adjacencies as we go and destroys the list of components.
 */
void mergeCycles(stHash *adjacencies, stHash *hyperChains);

#endif /* MERGECYCLES_H_ */
