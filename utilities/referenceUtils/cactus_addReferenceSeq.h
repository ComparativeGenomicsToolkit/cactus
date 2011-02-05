/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"

/*
 * Sep 27 2010: nknguyen@soe.ucsc.edu (based on the reference-parts in Benedict's cactus_MAFGenerator.c)
 * Adding the reference sequence into cactus structure
 */

typedef struct _referenceSequence {
    int32_t length;
    char *header;
    int32_t index;
    char *string;
} ReferenceSequence;

/*
 * Gets the sum of all the block lengths in the flower.
 */
/*static int32_t getTotalBlockLength(Flower *flower);

static ReferenceSequence *referenceSequence_construct(Flower *flower);

static void referenceSequence_destruct(ReferenceSequence *referenceSequence);

char *formatSequenceHeader(Sequence *sequence);

int32_t getNumberOnPositiveStrand(Block *block);


Sequence *getSequenceByHeader(Flower *flower, char *header);

void addReferenceSegmentToBlock(Block *block, ReferenceSequence *referenceSequence);

void block_metaSequence(Block *block, ReferenceSequence *referenceSequence);

void reference_walkDown(End *end, ReferenceSequence *referenceSequence);

void reference_walkUp(End *end, ReferenceSequence *referenceSequence);

void traverseReferenceBlocks(Flower *flower, ReferenceSequence *refseq);

MetaSequence *constructReferenceMetaSequence(Flower *flower, CactusDisk *cactusDisk, ReferenceSequence *refseq);

void constructReferenceSequence(MetaSequence *metaSequence, Flower *flower);
*/

char *getConsensusString(Block *block);
Flower *flower_addReferenceSequence(Flower *flower, CactusDisk *cactusDisk, char *header);


