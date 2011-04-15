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
#include "cactusUtils.h"
//#include "cactus_addReferenceSeq.h"

/*
 *Feb 7 2011: modify to include the referene sequence in the flower
 *Jan 19 2011: Correct to only put chain-segments in the same bed-record if there are links between them
 *(any two adjacent segments in the record must belong to two different blocks). 
 *Start a new bed record if hit an end with a self-edge
 *
 *Dec 09 2010: Modify code to adapt the Benedict's List data-structure
 *Sep 08 2010: nknguyen@soe.ucsc.edu
 *Print to outputFile bed record of each chain of the interested/inputed species.
 */

bool isLinked(End *end1, End *end2){
    //Return true if there is a link between end1 and end2, otherwise return false
    Link *link = group_getLink(end_getGroup(end1));
    if (link != NULL){
        End *link3end = link_get3End(link);
        End *link5end = link_get5End(link);
        end1 = end_getPositiveOrientation(end1);
        end2 = end_getPositiveOrientation(end2);
        if ( (link3end == end1 && link5end == end2) || (link3end == end2 && link5end == end1) ){
            return true;
        }
    }
    return false;
}

void block_getBED(Block *block, FILE *fileHandle, char *species, char *chr, int start, int level) {
    /*
     */
    Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
    Segment *segment;
    while((segment = block_getNext(instanceIterator)) != NULL) {
        Sequence *sequence = segment_getSequence(segment);
        if(sequence != NULL) {
            char *sequenceHeader = formatSequenceHeader(sequence);
            if(strcmp(sequenceHeader, species) == 0){
                if(!segment_getStrand(segment)){//if segment on the neg strand, reverse it
                    segment = segment_getReverse(segment);
                }
                Cap *cap5 = segment_get5Cap(segment);
                Cap *cap3 = segment_get3Cap(segment);
                //int s = cap_getCoordinate(cap5) + start - 2;
                //int e = cap_getCoordinate(cap3) + start + 1 - 2;
                int s = cap_getCoordinate(cap5) + start - 1;
                int e = cap_getCoordinate(cap3) + start + 1 - 1;
                fprintf(fileHandle, "%s %d %d %s.%d %d %s %d %d %d %d %d %d\n", chr, s, e, "NA", level, 0, ".", s, e, 0, 1, e-s,0);
                //fprintf(fileHandle, "%s %d %d %d %d %s %d %d %d %d %d %d\n", chr, s, e, level, 0, "+", s, e, 0, 1, e-s,0);
            }
            free(sequenceHeader);
        }
    }
    block_destructInstanceIterator(instanceIterator);
}

void chain_getBEDs(Chain *chain, Cap *cap, FILE *fileHandle, char *species, char *chr, int start, int level) {
    /*
     */
    st_logInfo("chain_getBEDs\n");
    char *chainName = cactusMisc_nameToString(chain_getName(chain));
    int chromStart;
    int chromEnd;
    int blockCount;
    struct IntList *blockSizes = constructEmptyIntList(0);
    struct IntList *blockStarts = constructEmptyIntList(0);
    Block *block;
    End *currEnd;
    End *prevOtherEnd = NULL;

    if(isStubCap(cap)){ //startCap of the current flower's thread
        cap = cap_getAdjacency(cap);
    }//else: current cap is somewhere in the middle of the thread
    
    if(cap == NULL){ return; }
    while(!isStubCap(cap)){
        currEnd = cap_getEnd(cap);
        block = end_getBlock(currEnd);
        Chain *currchain = block_getChain(block);
        if (chain == currchain){
            if( prevOtherEnd == NULL || isLinked(currEnd, prevOtherEnd) ){
                intListAppend(blockStarts, cap_getCoordinate(cap));
                intListAppend(blockSizes, cap_getCoordinate(cap_getOtherSegmentCap(cap)) - cap_getCoordinate(cap) +1);
            }else{
                chain_getBEDs(chain, cap, fileHandle, species, chr, start, level);
                break;
            }
        }
        //st_logInfo("cap %s, %d\t", cactusMisc_nameToString(cap_getName(cap)), cap_getCoordinate(cap));
        prevOtherEnd = end_getOtherBlockEnd(currEnd);
        moveCapToNextBlock(&cap);
        //st_logInfo("movedTo %s, %d\n", cactusMisc_nameToString(cap_getName(cap)), cap_getCoordinate(cap));
    }
    
    blockCount = blockStarts->length;
    if(blockCount == 0){
       //fprintf(fileHandle, "No block found for chain %s\n", chainName);
        return;
    }
    chromStart = blockStarts->list[0];
    chromEnd = blockStarts->list[blockCount -1] + blockSizes->list[blockCount -1];
    int i;
    for(i=0; i< blockCount; i++){
        blockStarts->list[i] -= chromStart;
    }
    //chromStart = chromStart + start -2;
    //chromEnd = chromEnd + start -2;
    chromStart = chromStart + start -1;
    chromEnd = chromEnd + start -1;
    fprintf(fileHandle, "%s %d %d %s.%d %d %s %d %d %s %d ", chr, chromStart, chromEnd, chainName, level, 0, ".", chromStart, chromEnd, "0", blockCount);
    //fprintf(fileHandle, "%s %d %d %d %d %s %d %d %s %d ", chr, chromStart, chromEnd, level, 0, "+", chromStart, chromEnd, "0", blockCount);
    //Print blockSizes
    for(i=0; i< blockCount; i++){
        fprintf(fileHandle, "%d,", blockSizes->list[i] );
    }
    fprintf(fileHandle, " ");
    for(i=0; i< blockCount; i++){
        fprintf(fileHandle, "%d,", blockStarts->list[i] );
    }
    fprintf(fileHandle, "\n");
}

void getBEDs(Flower *flower, FILE *fileHandle, char *species, int level){
    char *chr;
    char *tok;
    int start = 0;
    //int start = 1;
    char sep[] = ".";

    assert(species != NULL);
    strtok(stString_copy(species), sep); //species e.g "hg18"
    chr = strtok(NULL, sep);
    if(chr == NULL){
        chr = "";
    }else{
        tok = strtok(NULL, sep);//chromsize
        if(tok != NULL){
            sscanf(strtok(NULL, sep), "%d", &start);
        }
    }
    //sscanf(strtok(NULL, sep), "%d", &len);
    //sscanf(strtok(NULL, sep), "%d", &strand);
    //fprintf(stderr, "chrom *%s*, chromsize %d, start %d\n", chr, chrsize, start);

    //Get beds for all chain of current level:
    struct List *startCaps;
    startCaps = flower_getThreadStarts(flower, species);
    st_logInfo("Number of start Caps at flower %s is %d\n",
                cactusMisc_nameToString(flower_getName(flower)), startCaps->length);
    for(int i=0; i< startCaps->length; i++){
        Cap *startcap = startCaps->list[i];
        Flower_ChainIterator *chainIterator = flower_getChainIterator(flower);
        Chain *chain;
        while((chain = flower_getNextChain(chainIterator)) != NULL){
            chain_getBEDs(chain, startcap, fileHandle, species, chr, start, level);
        }
        flower_destructChainIterator(chainIterator);
    }

    //Get beds for non-trivial chains
    Flower_BlockIterator *blockIterator = flower_getBlockIterator(flower);
    Block *block;
    while((block = flower_getNextBlock(blockIterator)) != NULL){
        if(block_getChain(block) == NULL){//non-trivial chain
            block_getBED(block, fileHandle, species, chr, start, level);
        }
    }
    flower_destructBlockIterator(blockIterator);

    //Call child flowers recursively.
    Flower_GroupIterator *groupIterator = flower_getGroupIterator(flower);
    Group *group;
    level ++;
    //fprintf(fileHandle, "Go to level %d\n", level);
    while((group = flower_getNextGroup(groupIterator)) != NULL) {
        Flower *nestedFlower = group_getNestedFlower(group);
        if(nestedFlower != NULL) {
            //fprintf(fileHandle, "level %d\n", level);
            getBEDs(group_getNestedFlower(group), fileHandle, species, level); //recursive call.
        }
    }
    flower_destructGroupIterator(groupIterator);
}

/*struct List *getSequences(Flower *flower, char *name){
   //get names of all the sequences in 'flower' that have their names start with 'name'
   Sequence *sequence;
   struct List *seqs = constructEmptyList(0, free);
   Flower_SequenceIterator * seqIterator = flower_getSequenceIterator(flower);
   while((sequence = flower_getNextSequence(seqIterator)) != NULL){
      char *sequenceHeader = formatSequenceHeader(sequence);
      //if 'sequenceHeader' starts with 'name'
      //if(strstr(sequenceHeader, name) == sequenceHeader){
      if(strstr(sequenceHeader, name) != NULL){
         listAppend(seqs, sequenceHeader);
      }
      //free(sequenceHeader);
   }
   flower_destructSequenceIterator(seqIterator);
   return seqs;
}*/

void usage() {
    fprintf(stderr, "cactus_bedGenerator, version 0.2\n");
    fprintf(stderr, "Prints to output file all segments of the target sequence that are in blocks that contain both query & target\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --species: species (sequence) of interest.If species = 'reference' then will build the reference sequence named 'reference'\n");
    fprintf(stderr, "-c --cactusDisk : The cactus database conf string\n");
    fprintf(stderr, "-d --flowerName : The name of the flower (the key in the database)\n");
    fprintf(stderr, "-e --outputFile : The file to write the BEDs in.\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
    CactusDisk *cactusDisk;
    Flower *flower;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = NULL;
    char * outputFile = NULL;
    char * species = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while(1) {
        static struct option long_options[] = {
            { "logLevel", required_argument, 0, 'a' },
            { "species", required_argument, 0, 'b' },
            { "cactusDisk", required_argument, 0, 'c' },
            { "flowerName", required_argument, 0, 'd' },
            { "outputFile", required_argument, 0, 'e' },
            { "help", no_argument, 0, 'h' },
            { 0, 0, 0, 0 }
        };

        int option_index = 0;

        int key = getopt_long(argc, argv, "a:b:c:d:e:h", long_options, &option_index);

        if(key == -1) {
            break;
        }

        switch(key) {
            case 'a':
                logLevelString = stString_copy(optarg);
                break;
            case 'b':
                species = stString_copy(optarg);
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

    assert(cactusDiskDatabaseString != NULL);
    assert(flowerName != NULL);
    assert(outputFile != NULL);
    assert(species != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower name : %s\n", flowerName);
    st_logInfo("Output BED file : %s\n", outputFile);

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    st_logInfo("Set up the flower disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Parse the basic reconstruction problem
    ///////////////////////////////////////////////////////////////////////////

    flower = cactusDisk_getFlower(cactusDisk, cactusMisc_stringToName(flowerName));
    st_logInfo("Parsed the top level flower of the cactus tree to check\n");

    ///////////////////////////////////////////////////////////////////////////
    // Recursive check the flowers.
    ///////////////////////////////////////////////////////////////////////////

    int64_t startTime = time(NULL);
    FILE *fileHandle = fopen(outputFile, "w");
    //fprintf(fileHandle, "track name=%s\n", species);
    if(strstr(species, "reference") != NULL && getSequenceMatchesHeader(flower, "reference") == NULL){
        //flower_addReferenceSequence(flower, cactusDisk, species);
        fprintf(stderr, "No reference sequence found in cactusDisk\n");
        exit(EXIT_FAILURE);
    }
    //addReferenceSequenceTest(flower);
    char *seq;
    struct List *seqs = getSequences(flower, species);    
    for(int i = 0; i < seqs->length; i++){
        seq = seqs->list[i];
        st_logInfo("Getting beds for sequence \"%s\"\n", seq);
        getBEDs(flower, fileHandle, seq, 0);
    }
    fclose(fileHandle);
    st_logInfo("Got the beds in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
