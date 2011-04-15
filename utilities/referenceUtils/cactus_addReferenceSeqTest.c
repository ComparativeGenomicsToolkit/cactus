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
#include "cactus_addReferenceSeq.h"

/*
 */

void usage() {
    fprintf(stderr, "cactus_addReferenceSeq, version 0.0\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-b --name : name of the reference sequence\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

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

char *getChromName1(char *name, int num){
    char chrom[5];//assume there are less than 999 chroms! 
    sprintf(chrom, "%d", num);
    char *chromName = st_malloc(sizeof(char)*(strlen(name) + strlen(chrom) + 1));
    strcpy(chromName, name);
    strcat(chromName, "_");
    strcat(chromName, chrom);
    return chromName;
}

Sequence *getSequenceMatchesHeader(Flower *flower, char *header){
    //Returns the first Sequence whose name matches 'header'
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while((sequence = flower_getNextSequence(it)) != NULL){
        char *sequenceHeader = formatSequenceHeader(sequence);
        //if(strcmp(sequenceHeader, header) == 0){
        if(strstr(sequenceHeader, header) != NULL){
            flower_destructSequenceIterator(it);
            free(sequenceHeader);
            return sequence;
        }
    }
    flower_destructSequenceIterator(it);
    return NULL;
}

int32_t test(Flower *flower, char *name){
    //Test see if every block has 'reference' segment
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    int32_t totalSegmentLength = 0;
    //fprintf(stdout, "Number of blocks: %d\n", flower_getBlockNumber(flower));
    while((block = flower_getNextBlock(blockIterator)) != NULL) {
        Block_InstanceIterator *it = block_getInstanceIterator(block);
        Segment * segment;
        //fprintf(stdout, "\n%s\n", cactusMisc_nameToString(block_getName(block)));
        while((segment = block_getNext(it)) != NULL){
            Sequence *sequence = segment_getSequence(segment);
            if (sequence != NULL){
                char *header = formatSequenceHeader(sequence);
                //fprintf(stdout, "\t%s\n", cactusMisc_nameToString(segment_getName(segment)));
                //fprintf(stdout, "\t%s\n", header);
                if(strcmp(header, name) == 0){
                //if(strstr(header, name) != NULL){
                    totalSegmentLength += segment_getLength(segment);
		    break;
                }
                //fprintf(stdout, "\t%s\n", header);
            }
        }
        block_destructInstanceIterator(it);
    }
    flower_destructBlockIterator(blockIterator);
    
    Flower_GroupIterator *git = flower_getGroupIterator(flower);
    Group *group;
    while((group = flower_getNextGroup(git)) != NULL){
        Flower *nestedFlower = group_getNestedFlower(group);
        if(nestedFlower != NULL){
            totalSegmentLength += test(nestedFlower, name);
        }
    }
    flower_destructGroupIterator(git);

    //fprintf(stdout, "TotalSegmentLength %d\n\n", totalSegmentLength);
    return totalSegmentLength;
}

int main(int argc, char *argv[]) {
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    //char * flowerName = NULL;
    char * flowerName = "0";
    char *name = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { 
                { "logLevel", required_argument, 0, 'a' }, 
                { "name", required_argument, 0, 'b' }, 
                { "cactusDisk", required_argument, 0, 'c' }, 
		//{ "flowerName", required_argument, 0, 'd' },
                { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:h", long_options,
                &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                name = stString_copy(optarg);
                break;
            case 'c':
                cactusDiskDatabaseString = stString_copy(optarg);
                break;
            /*case 'd':
                flowerName = stString_copy(optarg);
                break;*/
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
    assert(name != NULL);
    assert(cactusDiskDatabaseString != NULL);

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
    st_logInfo("Sequence name : %s\n", name);

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
    //flower = flower_addReferenceSequence(flower, cactusDisk, name);
    //fprintf(stdout, "\n\nDONE adding reference sequence to flower *%s*\n\n", cactusMisc_nameToString(flower_getName(flower)));
    
    //Make sure that reference sequence has already been added to cactusDisk:
    if(getSequenceMatchesHeader(flower, name) == NULL){
        fprintf(stderr, "No reference sequence found in cactusDisk\n");
        exit(EXIT_FAILURE);
    }

    char *name1 = getChromName1(name, 1);
    int ref1 = test(flower, name1);
    fprintf(stdout, "*%s*: %d\n",name1, ref1);

    char *name2 = getChromName1(name, 2);
    int ref2 = test(flower, name2);
    fprintf(stdout, "%s: %d\n", name2, ref2);
    //test(flower, "reference_3");
    //test(flower, "reference_4");

    st_logInfo("Added the reference sequence in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
