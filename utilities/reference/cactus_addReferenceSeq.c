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

static int32_t getTotalBlockLength(Flower *flower) {
    /*
     * Gets the sum of all the block lengths in the flower.
     */
    int32_t length = 0;
    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    Block *block;
    while((block = flower_getNextBlock(blockIt)) != NULL) {
        if(block_getInstanceNumber(block) > 0){
            length += block_getLength(block);
        }
    }
    flower_destructBlockIterator(blockIt);
    Flower_GroupIterator *groupIt = flower_getGroupIterator(flower);
    Group *group;
    while((group = flower_getNextGroup(groupIt)) != NULL) {
        if(!group_isLeaf(group)) {
            length += getTotalBlockLength(group_getNestedFlower(group));
        }
    }
    flower_destructGroupIterator(groupIt);
    return length;
}

static ReferenceSequence *referenceSequence_construct(Flower *flower) {
    ReferenceSequence *referenceSequence = st_malloc(sizeof(ReferenceSequence));
    referenceSequence->index = 0;
    referenceSequence->header = stString_copy("reference");
    referenceSequence->length = getTotalBlockLength(flower);
    referenceSequence->string = st_malloc(sizeof(char)*(referenceSequence->length +1));
    //referenceSequence->string = "";
    //sscanf(sequenceHeader, "%s", referenceSequence->string);
    return referenceSequence;
}

static void referenceSequence_destruct(ReferenceSequence *referenceSequence) {
    free(referenceSequence->header);
    free(referenceSequence);
}

/*
 */

char *formatSequenceHeader(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) * (1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

int32_t getNumberOnPositiveStrand(Block *block) {
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment;
    int32_t i = 0;
    while ((segment = block_getNext(it)) != NULL) {
        if (segment_getChildNumber(segment) == 0) {
            if (segment_getStrand(segment)) {
                i++;
            }
        }
    }
    block_destructInstanceIterator(it);
    return i;
}

char *getConsensusString(Block *block) {
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment;
    while((segment = block_getNext(it)) != NULL) {
        if(segment_getSequence(segment) != NULL) {
            block_destructInstanceIterator(it);
            return segment_getString(segment);
        }
    }
    assert(0);
}

Sequence *getSequenceByHeader(Flower *flower, char *header){
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while((sequence = flower_getNextSequence(it)) != NULL){
        char *sequenceHeader = formatSequenceHeader(sequence);
        if(strcmp(sequenceHeader, header) == 0){
            flower_destructSequenceIterator(it);
            return sequence;
	}
    }
    flower_destructSequenceIterator(it);
    return NULL;    
}

void addReferenceSegmentToBlock(Block *block, ReferenceSequence *referenceSequence) {
    /*
     */
    Sequence *sequence = getSequenceByHeader(block_getFlower(block), "reference");
    assert(sequence != NULL);
    //Correct the orientation..
    if (getNumberOnPositiveStrand(block) == 0) {
        block = block_getReverse(block);
    }

    segment_construct2(block, referenceSequence->index, "+", sequence);
    referenceSequence->index += block_getLength(block);
}

void block_metaSequence(Block *block, ReferenceSequence *referenceSequence) {
    /*
     */
    //Correct the orientation..
    if (getNumberOnPositiveStrand(block) == 0) {
        block = block_getReverse(block);
    }

    if (block_getInstanceNumber(block) > 0) {
        char *instanceString = getConsensusString(block);
        //int length = strlen(instanceString) + strlen(referenceSequence->string);
        //referenceSequence->string = (char *)realloc(referenceSequence->string, (length+1)*sizeof(char));
        //st_logInfo("refString *%s*, instanceStr *%s*\n", referenceSequence->string, instanceString);
	referenceSequence->string = strcat(referenceSequence->string, instanceString);
        free(instanceString);
    }
}

void reference_walkDown(End *end, ReferenceSequence *referenceSequence);

void reference_walkUp(End *end, ReferenceSequence *referenceSequence) {
    assert(end != NULL);
    if (end_isBlockEnd(end)) {
        Block *block = end_getBlock(end);
        //st_logInfo("refseqStringLen %d, refseqLength %d\n", strlen(referenceSequence->string), referenceSequence->length);
        if(strlen(referenceSequence->string) == referenceSequence->length){
            addReferenceSegmentToBlock(block, referenceSequence);
	}else{
	    block_metaSequence(block,referenceSequence); 
        }
        //getMAFBlock(end_getBlock(end), referenceSequence);
        reference_walkDown(end_getOtherBlockEnd(end), referenceSequence);
    } else {
        assert(end_isAttached(end));
        Group *parentGroup = flower_getParentGroup(end_getFlower(end));
        if (parentGroup != NULL) {
            reference_walkUp(group_getEnd(parentGroup, end_getName(end)), referenceSequence);
        } else { //We reached the end of a pseudo-chromosome!
            assert(pseudoChromosome_get3End(pseudoAdjacency_getPseudoChromosome(end_getPseudoAdjacency(end))) == end);
        }
    }
}

void reference_walkDown(End *end, ReferenceSequence *referenceSequence) {
    assert(end != NULL);
    //assert(end_isAttached(end));
    Group *group = end_getGroup(end);
    if (group_isLeaf(group)) { //Walk across
        PseudoAdjacency *pseudoAdjacency = end_getPseudoAdjacency(end);
        assert(pseudoAdjacency != NULL);
        assert(pseudoAdjacency_get3End(pseudoAdjacency) == end || pseudoAdjacency_get5End(pseudoAdjacency) == end);
        if (pseudoAdjacency_get3End(pseudoAdjacency) == end) {
            end = pseudoAdjacency_get5End(pseudoAdjacency);
        } else {
            end = pseudoAdjacency_get3End(pseudoAdjacency);
        }
        //Now walk up
        reference_walkUp(end, referenceSequence);
    } else { //Walk down
        reference_walkDown(flower_getEnd(group_getNestedFlower(group), end_getName(end)), referenceSequence);
    }
}

MetaSequence *constructReferenceMetaSequence(Flower *flower, CactusDisk *cactusDisk, ReferenceSequence *refseq) {
    /*
     */
    st_logInfo("Getting reference metasequence...\n");
    Event *event = eventTree_getRootEvent(flower_getEventTree(flower));
    Name eventName = event_getName(event);
    MetaSequence *metaSequence;
    int32_t start = 0;

    Reference *reference = flower_getReference(flower);
    assert(reference != NULL);
    Reference_PseudoChromosomeIterator *it =
            reference_getPseudoChromosomeIterator(reference);
    PseudoChromosome *pseudoChromosome;
    while ((pseudoChromosome = reference_getNextPseudoChromosome(it)) != NULL) {
        End *end = pseudoChromosome_get5End(pseudoChromosome);
        assert(!end_isBlockEnd(end));
        reference_walkDown(end, refseq);
    }
    metaSequence = metaSequence_construct(start, refseq->length, refseq->string, refseq->header, eventName, cactusDisk);
    st_logInfo("Got metasequence: name *%s*, start %d, length %d, header *%s*\n", 
                cactusMisc_nameToString(metaSequence_getName(metaSequence)), metaSequence_getStart(metaSequence),
                metaSequence_getLength(metaSequence), metaSequence_getHeader(metaSequence));
    reference_destructPseudoChromosomeIterator(it);
    return metaSequence;
}

void constructReferenceSequence(MetaSequence *metaSequence, Flower *flower){
    sequence_construct(metaSequence, flower);
    Flower_GroupIterator *it = flower_getGroupIterator(flower);
    Group *group;
    while((group = flower_getNextGroup(it))!= NULL){
        Flower *nestedFlower = group_getNestedFlower(group);
	if(nestedFlower != NULL){
	    constructReferenceSequence(metaSequence, flower); 
	}
    }
    flower_destructGroupIterator(it);
}

void addReferenceSegments(Flower *flower, ReferenceSequence *refseq){
    Reference *reference = flower_getReference(flower);
    assert(reference != NULL);
    Reference_PseudoChromosomeIterator *it =
            reference_getPseudoChromosomeIterator(reference);
    PseudoChromosome *pseudoChromosome;
    while ((pseudoChromosome = reference_getNextPseudoChromosome(it)) != NULL) {
        End *end = pseudoChromosome_get5End(pseudoChromosome);
        assert(!end_isBlockEnd(end));
        reference_walkDown(end, refseq);
    }
    reference_destructPseudoChromosomeIterator(it);
}

void usage() {
    fprintf(stderr, "cactus_addReferenceSeq, version 0.0\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { 
                { "logLevel", required_argument, 0, 'a' }, 
                { "cactusDisk", required_argument, 0, 'c' }, 
		{ "flowerName", required_argument, 0, 'd' },
                { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:h", long_options,
                &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            case 'd':
                flowerName = stString_copy(optarg);
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(flowerName != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if (logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower name : %s\n", flowerName);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(
            cactusDiskDatabaseString);
    CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    Flower *flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(
            flowerName));
    st_logInfo("Parsed the top level flower of the cactus tree to check\n");

    ///////////////////////////////////////////////////////////////////////////
    // Recursive check the flowers.
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    ReferenceSequence *referenceSequence = NULL;
    referenceSequence = referenceSequence_construct(flower);
    st_logInfo("initialize referenceSequence: index %d, length %d, header *%s*, string *%s*\n", 
               referenceSequence->index, referenceSequence->length, referenceSequence->header
               ,referenceSequence->string);
    MetaSequence *metaSequence = constructReferenceMetaSequence(flower, cactusDisk, referenceSequence);
    constructReferenceSequence(metaSequence, flower);
    addReferenceSegments(flower, referenceSequence);
    referenceSequence_destruct(referenceSequence);
    st_logInfo("Added the reference sequence in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
