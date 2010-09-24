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
 * Object for representing the reference sequence.
 */
typedef struct _referenceSequence {
    int32_t length;
    char *header;
    int32_t index;
} ReferenceSequence;

static int32_t getTotalBlockLength(Flower *flower) {
    /*
     * Gets the sum of all the block lengths in the flower.
     */
    int32_t length = 0;
    Flower_BlockIterator *blockIt = flower_getBlockIterator(flower);
    Block *block;
    while((block = flower_getNextBlock(blockIt)) != NULL) {
        length += block_getLength(block);
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
    return referenceSequence;
}

static void referenceSequence_destruct(ReferenceSequence *referenceSequence) {
    free(referenceSequence->header);
    free(referenceSequence);
}

/*
 * The script outputs a maf file containing all the block in a flower and its descendants.
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

static void getMAFBlockP2(Segment *segment, FILE *fileHandle) {
    assert(segment != NULL);
    Sequence *sequence = segment_getSequence(segment);
    if (sequence != NULL) {
        char *sequenceHeader = formatSequenceHeader(sequence);
        int32_t start;
        if (segment_getStrand(segment)) {
            start = segment_getStart(segment) - sequence_getStart(sequence);
        } else { //start with respect to the start of the reverse complement sequence
            start = (sequence_getStart(sequence) + sequence_getLength(sequence)
                    - 1) - segment_getStart(segment);
        }
        int32_t length = segment_getLength(segment);
        char *strand = segment_getStrand(segment) ? "+" : "-";
        int32_t sequenceLength = sequence_getLength(sequence);
        char *instanceString = segment_getString(segment);
        fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n", sequenceHeader,
                start, length, strand, sequenceLength, instanceString);
        free(instanceString);
        free(sequenceHeader);
    }
}

static void getMAFBlockP(Segment *segment, FILE *fileHandle) {
    int32_t i;
    for (i = 0; i < segment_getChildNumber(segment); i++) {
        getMAFBlockP(segment_getChild(segment, i), fileHandle);
    }
    getMAFBlockP2(segment, fileHandle);
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

void getMAFBlock(Block *block, FILE *fileHandle, ReferenceSequence *referenceSequence) {
    /*
     * Outputs a MAF representation of the block to the given file handle.
     */
    //Correct the orientation..
    if (getNumberOnPositiveStrand(block) == 0) {
        block = block_getReverse(block);
    }
    if (block_getInstanceNumber(block) > 0) {
        //Add in the header
        if (block_getRootInstance(block) != NULL) {
            /* Get newick tree string with internal labels and no unary events */
            char *newickTreeString = block_makeNewickString(block, 1, 0);
            assert(newickTreeString != NULL);
            fprintf(fileHandle, "a score=%i tree='%s'\n",
                    block_getLength(block) * block_getInstanceNumber(block),
                    newickTreeString);
            free(newickTreeString);
        } else {
            fprintf(fileHandle, "a score=%i\n", block_getLength(block)
                    * block_getInstanceNumber(block));
        }
        //Now for the reference segment
        if (referenceSequence != NULL) {
            char *instanceString = getConsensusString(block);
            fprintf(fileHandle, "s\t%s\t%i\t%i\t%s\t%i\t%s\n",
                    referenceSequence->header, referenceSequence->index,
                    block_getLength(block), "+", referenceSequence->length,
                    instanceString);
            free(instanceString);
            referenceSequence->index += block_getLength(block);
        }
        //Now add the blocks in
        if (block_getRootInstance(block) != NULL) {
            assert(block_getRootInstance(block) != NULL);
            getMAFBlockP(block_getRootInstance(block), fileHandle);
            fprintf(fileHandle, "\n");
        } else {
            Block_InstanceIterator *iterator = block_getInstanceIterator(block);
            Segment *segment;
            while ((segment = block_getNext(iterator)) != NULL) {
                getMAFBlockP2(segment, fileHandle);
            }
            block_destructInstanceIterator(iterator);
            fprintf(fileHandle, "\n");
        }
    }
}

void getMAFSReferenceOrdered_walkDown(End *end, FILE *fileHandle, ReferenceSequence *referenceSequence);

void getMAFSReferenceOrdered_walkUp(End *end, FILE *fileHandle, ReferenceSequence *referenceSequence) {
    assert(end != NULL);
    if (end_isBlockEnd(end)) {
        getMAFBlock(end_getBlock(end), fileHandle, referenceSequence);
        getMAFSReferenceOrdered_walkDown(end_getOtherBlockEnd(end), fileHandle, referenceSequence);
    } else {
        assert(end_isAttached(end));
        Group *parentGroup = flower_getParentGroup(end_getFlower(end));
        if (parentGroup != NULL) {
            getMAFSReferenceOrdered_walkUp(group_getEnd(parentGroup,
                    end_getName(end)), fileHandle, referenceSequence);
        } else { //We reached the end of a pseudo-chromosome!
            assert(pseudoChromosome_get3End(pseudoAdjacency_getPseudoChromosome(end_getPseudoAdjacency(end))) == end);
        }
    }
}

void getMAFSReferenceOrdered_walkDown(End *end, FILE *fileHandle, ReferenceSequence *referenceSequence) {
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
        getMAFSReferenceOrdered_walkUp(end, fileHandle, referenceSequence);
    } else { //Walk down
        getMAFSReferenceOrdered_walkDown(flower_getEnd(group_getNestedFlower(
                group), end_getName(end)), fileHandle, referenceSequence);
    }
}

void getMAFsReferenceOrdered(Flower *flower, FILE *fileHandle, ReferenceSequence *referenceSequence) {
    /*
     * Outputs MAF representations of all the block in the flower and its descendants, ordered
     * according to the reference ordering.
     */
    Reference *reference = flower_getReference(flower);
    assert(reference != NULL);
    Reference_PseudoChromosomeIterator *it =
            reference_getPseudoChromosomeIterator(reference);
    PseudoChromosome *pseudoChromosome;
    while ((pseudoChromosome = reference_getNextPseudoChromosome(it)) != NULL) {
        End *end = pseudoChromosome_get5End(pseudoChromosome);
        assert(!end_isBlockEnd(end));
        getMAFSReferenceOrdered_walkDown(end, fileHandle, referenceSequence);
    }
    reference_destructPseudoChromosomeIterator(it);
}

void getMAFs(Flower *flower, FILE *fileHandle) {
    /*
     * Outputs MAF representations of all the block sin the flower and its descendants.
     */

    //Make MAF blocks for each block
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    while ((block = flower_getNextBlock(blockIterator)) != NULL) {
        getMAFBlock(block, fileHandle, NULL);
    }
    flower_destructBlockIterator(blockIterator);

    //Call child flowers recursively.
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    while ((group = flower_getNextGroup(groupIterator)) != NULL) {
        if (!group_isLeaf(group)) {
            getMAFs(group_getNestedFlower(group), fileHandle); //recursive call.
        }
    }
    flower_destructGroupIterator(groupIterator);
}

void makeMAFHeader(Flower *flower, FILE *fileHandle) {
    fprintf(fileHandle, "##maf version=1 scoring=N/A\n");
    char *cA = eventTree_makeNewickString(flower_getEventTree(flower));
    fprintf(fileHandle, "# cactus %s\n\n", cA);
    free(cA);
}

void usage() {
    fprintf(stderr, "cactus_mafGenerator, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr, "-e --outputFile : The file to write the MAFs in.\n");
    fprintf(stderr,
            "-f --orderByReference : Order the blocks by the reference ordering.\n");
    fprintf(
            stderr,
            "-g --includeReferenceSequence : Include a reference sequence in the blocks, must be used in conjunction with --orderByReference.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;
    char * outputFile = NULL;
    bool orderByReference = 0;
    bool includeReferenceSequence = 0;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "cactusDisk", required_argument,
                0, 'c' }, { "flowerName", required_argument, 0, 'd' }, {
                "outputFile", required_argument, 0, 'e' }, {
                "orderByReference", no_argument, 0, 'f' }, {
                "includeReferenceSequence", no_argument, 0, 'g' }, { "help",
                no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:c:d:e:fgh", long_options,
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
            case 'e':
                outputFile = stString_copy(optarg);
                break;
            case 'f':
                orderByReference = !orderByReference;
                break;
            case 'g':
                includeReferenceSequence = !includeReferenceSequence;
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
    assert(outputFile != NULL);
    if (includeReferenceSequence) {
        if (!orderByReference) {
            stExcept_new(
                    "MAF_GENERATOR_EXCEPTION",
                    "You have specified to include the reference sequence by not to order by reference, this is not possible currently");
        }
    }

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
    st_logInfo("Output MAF file : %s\n", outputFile);

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
    FILE *fileHandle = fopen(outputFile, "w");
    makeMAFHeader(flower, fileHandle);
    if (orderByReference) {
        st_logInfo("Ordering by reference\n");
        ReferenceSequence *referenceSequence = NULL;
        if(includeReferenceSequence) {
            referenceSequence = referenceSequence_construct(flower);
        }
        getMAFsReferenceOrdered(flower, fileHandle, referenceSequence);
        if(referenceSequence != NULL) {
            referenceSequence_destruct(referenceSequence);
        }
    } else {
        getMAFs(flower, fileHandle);
    }
    fclose(fileHandle);
    st_logInfo("Got the mafs in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
