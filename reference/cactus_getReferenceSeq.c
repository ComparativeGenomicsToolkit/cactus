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

/*
 */

void usage() {
    fprintf(stderr, "cactus_getReferenceSeq, version 0.0\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr,
            "-b --name : name of the reference sequence's event\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr,
            "-e --outputFile : Name of output fasta file\n");
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

void getReferenceSequences(FILE *fileHandle, Flower *flower, char *name){
   //get names of all the sequences in 'flower' that have their names start with 'name'
   Sequence *sequence;
   Flower_SequenceIterator * seqIterator = flower_getSequenceIterator(flower);
   while((sequence = flower_getNextSequence(seqIterator)) != NULL){
      char *sequenceHeader = formatSequenceHeader(sequence);
      st_logInfo("Sequence %s\n", sequenceHeader);
      if(strstr(sequenceHeader, name) != NULL){
          if(sequence_getLength(sequence) == 0){//skip empty sequences
              free(sequenceHeader);
              continue;
          }
          fprintf(fileHandle, ">%s\n", sequenceHeader);
          const int32_t lineLength = 10000000;
          const int32_t end = sequence_getLength(sequence) + sequence_getStart(sequence);
          for(int32_t i=sequence_getStart(sequence); i<end; i += lineLength) {
              char *cA = sequence_getString(sequence, i, end - i < lineLength ?  end - i : lineLength, 1);
              fprintf(fileHandle, "%s\n", cA);
              free(cA);
          }
      }
      free(sequenceHeader);
   }
   flower_destructSequenceIterator(seqIterator);
   return;
}


int main(int argc, char *argv[]) {
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;
    char * outputFile = NULL;
    char *referenceEventString = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { 
                { "logLevel", required_argument, 0, 'a' }, 
                { "referenceEventString", required_argument, 0, 'b' },
                { "cactusDisk", required_argument, 0, 'c' }, 
		{ "flowerName", required_argument, 0, 'd' },
		{ "outputFile", required_argument, 0, 'e' },
                { "help", no_argument, 0, 'h' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e:h", long_options,
                &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                referenceEventString = stString_copy(optarg);
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
    assert(referenceEventString != NULL);
    assert(cactusDiskDatabaseString != NULL);
    assert(outputFile != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower name : %s\n", flowerName);
    st_logInfo("Sequence name : %s\n", referenceEventString);
    st_logInfo("Output file : %s\n", outputFile);

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

    //int64_t startTime = time(NULL);
    //flower = flower_addReferenceSequence(flower, cactusDisk, name);
    //st_logInfo("Added the reference sequence in %i seconds/\n", time(NULL) - startTime);

    //Make sure that referenceSequence has already been added:
    if(getSequenceMatchesHeader(flower, referenceEventString) == NULL){
        fprintf(stderr, "No reference sequence found in cactusDisk\n");
        exit(EXIT_FAILURE); 
    }

    FILE *fileHandle = fopen(outputFile, "w");
    getReferenceSequences(fileHandle, flower, referenceEventString);
    fclose(fileHandle);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
