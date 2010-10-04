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
 *Sep 08 2010: nknguyen@soe.ucsc.edu
 *Print to outputFile bed record of each chain of the interested/inputed species.
 */

char *formatSequenceHeader(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if(strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) *(1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    }
    else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

bool isStubCap(Cap *cap){
    /*
     *Return true if cap is of a stubEnd, otherwise return false
     */
    assert(cap != NULL);
    End *end = cap_getEnd(cap);
    assert(end != NULL);
    return (end_isStubEnd(end)) ? true : false;
}

int flower_getThreadStarts(Flower *flower, char *name, Cap ***startCaps){
    /*
     *Get 3' end Stubs of the sequence by its name
     */
    int num = 0;
    Cap *cap;
    Flower_CapIterator *capIterator = flower_getCapIterator(flower);
    while((cap= flower_getNextCap(capIterator)) != NULL){
        if(isStubCap(cap)){//dead end or inherited end
            Sequence *sequence = cap_getSequence(cap);
            if(sequence == NULL){continue;}
            char *sequenceHeader = formatSequenceHeader(sequence);
            if(strcmp(sequenceHeader, name) == 0){//cap matched with name
                if(!cap_getStrand(cap)){//if cap is on negative strand - reverse it
                    cap = cap_getReverse(cap);
                }
                if(cap_getSide(cap)){ continue; }//if cap is the 5' end, ignore
                num++;
                if(num == 1){
                    (*startCaps) = st_malloc( sizeof(Cap *) );
                }else{
                    (*startCaps) = (Cap **)realloc( (*startCaps), num*sizeof(Cap *) );
                }
                *((*startCaps) + num -1) = cap;
            }
            free(sequenceHeader);
        }
    }
    flower_destructCapIterator(capIterator);
    return num;
}

void moveCapToNextBlock(Cap **cap){
    /*Move cap to the next block that its segment is adjacency to*/
    assert((*cap) != NULL);
    st_logInfo("cap %s, %d\t", cactusMisc_nameToString(cap_getName(*cap)), cap_getCoordinate(*cap));
    if(isStubCap(*cap)){
        st_logInfo("STUB\n");
        return;
    }
    Cap *otherCap = cap_getOtherSegmentCap(*cap);
    st_logInfo("other cap %s, %d\t", cactusMisc_nameToString(cap_getName(otherCap)), cap_getCoordinate(otherCap));
    *cap = cap_getAdjacency(otherCap);
    assert(*cap != NULL);
    st_logInfo("moved-cap %s, %d\n", cactusMisc_nameToString(cap_getName(*cap)), cap_getCoordinate(*cap));
    
    /*Cap *adjCap = cap_getAdjacency(*cap);
    assert(adjCap != NULL);
    if(cap_getEnd(*cap) == cap_getEnd(adjCap)){//DOUBLE CHECK ... self connected end
        *cap = adjCap;
    }else{
        *cap = cap_getAdjacency(cap_getOtherSegmentCap(*cap));
    //}*/
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
                int s = cap_getCoordinate(cap5) + start - 2;
                int e = cap_getCoordinate(cap3) + start + 1 - 2;
                fprintf(fileHandle, "%s %d %d %s.%d %d %s %d %d %d %d %d %d\n", chr, s, e, "NA", level, 0, "+", s, e, 0, 1, e-s,0);
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
    int blockCount = 0;
    int chromStart;
    int chromEnd;
    int *blockSizes = NULL;
    int *blockStarts = NULL;
    Block *block;

    cap = cap_getAdjacency(cap);
    
    if(cap == NULL){ return; }
    /*block = end_getBlock(cap_getEnd(cap));
    if(block == NULL){
       moveCapToNextBlock(&cap);
    }*/
    while(!isStubCap(cap)){
        block = end_getBlock(cap_getEnd(cap));
        Chain *currchain = block_getChain(block);
        if (chain == currchain){
            blockCount++;
            if(blockCount ==1){
                blockSizes = st_malloc( sizeof(int) );
                blockStarts = st_malloc( sizeof(int) );
            }else{
                blockSizes = (int *)realloc( blockSizes, blockCount*sizeof(int) );
                blockStarts = (int *)realloc( blockStarts, blockCount*sizeof(int) );
            }
            *(blockStarts + blockCount -1) = cap_getCoordinate(cap);
            *(blockSizes + blockCount -1) = cap_getCoordinate(cap_getOtherSegmentCap(cap)) - cap_getCoordinate(cap) + 1;
        }
        
        //st_logInfo("cap %s, %d\t", cactusMisc_nameToString(cap_getName(cap)), cap_getCoordinate(cap));
        moveCapToNextBlock(&cap);
        //st_logInfo("movedTo %s, %d\n", cactusMisc_nameToString(cap_getName(cap)), cap_getCoordinate(cap));
    }
    
    if(blockCount == 0){
       //fprintf(fileHandle, "No block found for chain %s\n", chainName);
        return;
    }
    chromStart = *(blockStarts);
    chromEnd = *(blockStarts + blockCount -1) + *(blockSizes + blockCount -1);
    int i;
    for(i=0; i< blockCount; i++){
        *(blockStarts + i) -= chromStart;
    }
    chromStart = chromStart + start -2;
    chromEnd = chromEnd + start -2;
    fprintf(fileHandle, "%s %d %d %s.%d %d %s %d %d %s %d ", chr, chromStart, chromEnd, chainName, level, 0, "+", chromStart, chromEnd, "0", blockCount);
    //fprintf(fileHandle, "%s %d %d %d %d %s %d %d %s %d ", chr, chromStart, chromEnd, level, 0, "+", chromStart, chromEnd, "0", blockCount);
    //Print blockSizes
    for(i=0; i< blockCount; i++){
        fprintf(fileHandle, "%d,", *(blockSizes + i));
    }
    fprintf(fileHandle, " ");
    for(i=0; i< blockCount; i++){
        fprintf(fileHandle, "%d,", *(blockStarts + i));
    }
    fprintf(fileHandle, "\n");
}

int32_t getSeqLength(Flower *flower, char *header){
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while((sequence = flower_getNextSequence(it)) != NULL){
        char *sequenceHeader = formatSequenceHeader(sequence);
        if(strcmp(sequenceHeader, header) == 0){
            flower_destructSequenceIterator(it);
            return sequence_getLength(sequence);
        }
    }
    flower_destructSequenceIterator(it);
    return 0;
}



void getBEDs(Flower *flower, FILE *fileHandle, char *species, int level){
    char *chr;
    char *tok;
    //int start = 0;
    int start = 1;
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
    Cap **startCaps;
    int capNum;
    capNum = flower_getThreadStarts(flower, species, &startCaps);
    st_logInfo("Number of start Caps at flower %s is %d\n", cactusMisc_nameToString(flower_getName(flower)), capNum);
    for(int i=0; i< capNum; i++){
        Cap *startcap = *(startCaps + i);
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

int getSequences(Flower *flower, char ***seqs, char *name){
   //get names of all the sequences in 'flower' that have their names start with 'name'
   int num = 0;
   Sequence *sequence;
   Flower_SequenceIterator * seqIterator = flower_getSequenceIterator(flower);
   while((sequence = flower_getNextSequence(seqIterator)) != NULL){
      char *sequenceHeader = formatSequenceHeader(sequence);
      //if 'sequenceHeader' starts with 'name'
      //if(strstr(sequenceHeader, name) == sequenceHeader){
      if(strstr(sequenceHeader, name) != NULL){
         if(num == 0){
            (*seqs) = st_malloc(sizeof(char *) + 1);
         }else{
            (*seqs) =  (char **)realloc( (*seqs), num*(sizeof(char *) +1));
         }
         *((*seqs) + num) = sequenceHeader;
         num++;
      }
      //free(sequenceHeader);
   }
   flower_destructSequenceIterator(seqIterator);
   return num;
}

void usage() {
    fprintf(stderr, "cactus_bedGenerator, version 0.2\n");
    fprintf(stderr, "Prints to output file all segments of the target sequence that are in blocks that contain both query & target\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --species: species (sequence) of interest\n");
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
    fprintf(fileHandle, "track name=%s\n", species);
    if(strstr(species, "reference") != NULL){
        flower_addReferenceSequence(flower, cactusDisk, species);
    }
    //addReferenceSequenceTest(flower);
    char **seqs;
    char *seq;
    int numseqs = getSequences(flower, &seqs, species);    
    for(int i = 0; i < numseqs; i++){
        seq = *(seqs + i);
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
