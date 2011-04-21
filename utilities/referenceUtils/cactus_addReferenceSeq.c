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

typedef struct _refseq {
    int32_t length;
    char *header;
    int32_t index;
    char *string;
} ReferenceSequence;

End *getPseudoAdjacentEnd(End *end){
    assert(end != NULL);
    PseudoAdjacency *pseudoAdjacency = end_getPseudoAdjacency(end);
    assert(pseudoAdjacency != NULL);
    assert(pseudoAdjacency_get3End(pseudoAdjacency) == end || pseudoAdjacency_get5End(pseudoAdjacency) == end);
    //assert( pseudoAdjacency_get3End(pseudoAdjacency) != pseudoAdjacency_get5End(pseudoAdjacency) );
    if (pseudoAdjacency_get3End(pseudoAdjacency) == end) {
        end = pseudoAdjacency_get5End(pseudoAdjacency);
    } else {
        end = pseudoAdjacency_get3End(pseudoAdjacency);
    }

    return end;
}

static int32_t pseudoChromosome_getLength(End *end){
    /*
     *Return the total number of bases of the pseudochromosome (and of all 
     *pseudochromosomes at lower flowers within this pseudochrom) that 'end' belongs to, 
     */
    int32_t len = 0;
    assert(end_isStubEnd(end));
    Group *group;
    End *inheritedEnd;

    group = end_getGroup(end);
    if(!group_isLeaf(group)){//has lower level
	inheritedEnd = flower_getEnd(group_getNestedFlower(group), end_getName(end));
        len += pseudoChromosome_getLength(inheritedEnd);
    }

    end = getPseudoAdjacentEnd(end);
    while(end_isBlockEnd(end)){
        Block *block = end_getBlock(end);
        //if (block_getInstanceNumber(block) > 0) {
	len += block_getLength(block);
        //}

	end = end_getOtherBlockEnd(end); 
        group = end_getGroup(end);
        if(!group_isLeaf(group)){//has lower level
	    inheritedEnd = flower_getEnd(group_getNestedFlower(group), end_getName(end));
            len += pseudoChromosome_getLength(inheritedEnd);
        }
        end = getPseudoAdjacentEnd(end);
    }
    return len;
}

static ReferenceSequence *referenceSequence_construct(Flower *flower, char *header, int length) {
    ReferenceSequence *refseq = st_malloc(sizeof(ReferenceSequence));
    refseq->index = 0;
    refseq->header = stString_copy(header);
    refseq->length = length;
    refseq->string = st_malloc(sizeof(char)*(length +1));
    strcpy(refseq->string, "");
    return refseq;
}

static void referenceSequence_destruct(ReferenceSequence *refseq) {
    free(refseq->string);
    free(refseq->header);
    free(refseq);
}

char *formatSequenceHeader1(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) * (1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}

/*int32_t getNumberOnPositiveStrand(Block *block) {
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

int32_t ref_getNumberOnPositiveStrand(Block *block) {
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
}*/

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
    return NULL;
}

Sequence *getSequenceByHeader(Flower *flower, char *header){
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while((sequence = flower_getNextSequence(it)) != NULL){
        char *sequenceHeader = formatSequenceHeader1(sequence);
        if(strcmp(sequenceHeader, header) == 0){
            flower_destructSequenceIterator(it);
            free(sequenceHeader);
            return sequence;
	}
    }
    flower_destructSequenceIterator(it);
    return NULL;    
}

/*
 *Return: the first cap in input 'end' that belongs to sequence 'header'
 *        NULL if not found
 */
Cap *end_getCapByHeader(End *end, char *header){
    End_InstanceIterator *it = end_getInstanceIterator(end);
    Cap *cap;
    while((cap = end_getNext(it)) != NULL){
        Sequence *sequence = cap_getSequence(cap);
	if(sequence !=NULL){
	    char *sequenceHeader = formatSequenceHeader1(sequence);
	    if(strcmp(sequenceHeader, header) == 0){
                end_destructInstanceIterator(it);
                free(sequenceHeader);
		return cap;
	    }
        }
    }
    end_destructInstanceIterator(it);
    return NULL;
}

Segment *addReferenceSegmentToBlock(End *end, ReferenceSequence *refseq) {
    /*
     */
    if(!end_getSide(end)){
        end = end_getReverse(end);
    }
    Block *block = end_getBlock(end);
    Sequence *sequence = getSequenceByHeader(block_getFlower(block), refseq->header);
    assert(sequence != NULL);
    //Correct the orientation..
    /*if (getNumberOnPositiveStrand(block) == 0) {
        block = block_getReverse(block);
    }*/

    //Adding segment to block
    Segment *segment = segment_construct2(block, refseq->index, "+", sequence);
    refseq->index += block_getLength(block);
    return segment;
}

void block_metaSequence(End *end, ReferenceSequence *refseq) {
    /*
     *Adding sequence of end's block to MetaSequence
     */
    if(!end_getSide(end)){
        end = end_getReverse(end);
    }
    Block *block = end_getBlock(end);
    //Correct the orientation..
    /*if (ref_getNumberOnPositiveStrand(block) == 0) {
        block = block_getReverse(block);
    }*/
    if (block_getInstanceNumber(block) > 0) {
        char *instanceString = getConsensusString(block);
	refseq->string = strcat(refseq->string, instanceString);
        free(instanceString);
    }else{//if block doesn't have any instance, add 'N'
        assert(block_getLength(block) == 1);
	refseq->string = strcat(refseq->string, "N");
    }
}

Cap *copyRefCapToLowerFlowers(Cap *cap){
    assert(cap != NULL);
    End *end = cap_getEnd(cap);
    Group *group = end_getGroup(end);
    Flower *nestedflower = group_getNestedFlower(group);
    if (nestedflower != NULL) { //has lower level
        End *inheritedEnd = flower_getEnd(nestedflower, end_getName(end));
	if(cap_getSide(cap) != end_getSide(inheritedEnd)){//make sure end & inheritedEnd are the same
	    inheritedEnd = end_getReverse(inheritedEnd);
	}
	cap = cap_copyConstruct(inheritedEnd, cap);
        copyRefCapToLowerFlowers(cap);
    }
    return cap;
}

void reference_walkDown(End *end, ReferenceSequence *refseq);

void reference_walkUp(End *end, ReferenceSequence *refseq) {
    assert(end != NULL);
    if (end_isBlockEnd(end)) {
        if(strlen(refseq->string) == refseq->length){
            //if(block_getInstanceNumber(block) > 0){
            Segment *segment = addReferenceSegmentToBlock(end, refseq);
            copyRefCapToLowerFlowers(segment_get5Cap(segment));
            copyRefCapToLowerFlowers(segment_get3Cap(segment));
            //}
	}else{
	    block_metaSequence(end, refseq); 
        }
        reference_walkDown(end_getOtherBlockEnd(end), refseq);
    } else {
        assert(end_isAttached(end));
        Group *parentGroup = flower_getParentGroup(end_getFlower(end));
        if (parentGroup != NULL) {
            reference_walkUp(group_getEnd(parentGroup, end_getName(end)), refseq);
        } else { //We reached the end of a pseudo-chromosome!
            assert(pseudoChromosome_get3End(pseudoAdjacency_getPseudoChromosome(end_getPseudoAdjacency(end))) == end);
            //adding last Stub (5')
            if(refseq->index > 0){
                Sequence *sequence = getSequenceByHeader(end_getFlower(end), refseq->header);
	        Cap *endCap = cap_construct2(end, refseq->index, true, sequence);
                copyRefCapToLowerFlowers(endCap);
            }
        }
    }
}

void reference_walkDown(End *end, ReferenceSequence *refseq) {
    assert(end != NULL);
    //assert(end_isAttached(end));
    Group *group = end_getGroup(end);
    if (group_isLeaf(group)) { //Walk across
        end = getPseudoAdjacentEnd(end);
        //Now walk up
        reference_walkUp(end, refseq);
    } else { //Walk down
        reference_walkDown(flower_getEnd(group_getNestedFlower(group), end_getName(end)), refseq);
    }
}

void addAdj(End *end, End *adjEnd, char *header){
    /*
     *Add adjacency between caps of 'header' sequence in 'end' and 'adjEnd'
     */
    Cap *cap = end_getCapByHeader(end, header);
    Cap *cap2 = end_getCapByHeader(adjEnd, header);
    assert(cap != NULL);
    assert(cap2 != NULL);
    if(cap_getSide(cap) == cap_getSide(cap2)){
        cap2 = cap_getReverse(cap2);
    }
    cap_makeAdjacent(cap, cap2); 
}

/*
 *1/ Walk along (across) current flower and add adjacencies between caps of the 
 *reference segments. Once reach the end of the pseudochromosome, 
 *return the last block-end cap
 *2/ Recursively add adjs of lower-level flower
 */
void addReferenceAdjacencies(End *end, char *header){
    End *adjEnd = getPseudoAdjacentEnd(end);//start of block
    //assert(end_isBlockEnd(adjEnd));    
    /*if(end_isStubEnd(adjEnd)){
        return;
    }*/    

    //while( end_isBlockEnd(adjEnd) ){
    while(true){
        Group *group = end_getGroup(end);
        if(!group_isLeaf(group)){//has lower level
	    End *inheritedEnd = flower_getEnd(group_getNestedFlower(group), end_getName(end));
            addReferenceAdjacencies(inheritedEnd, header);
        }
 
        addAdj(end, adjEnd, header);
        if(end_isStubEnd(adjEnd)){
            break;
        }else{
            end = end_getOtherBlockEnd(adjEnd);
            adjEnd = getPseudoAdjacentEnd(end);
        }
    }
    //add adjacency to the last stub:
    //addAdj(end, adjEnd, header);

    return;
}

MetaSequence *constructReferenceMetaSequence(End *end, CactusDisk *cactusDisk, ReferenceSequence *refseq) {
    /*
     *Traverse pseudochromosome (of 'end') and its lower levels to get the reference MetaSequence
     */
    st_logInfo("Getting reference MetaSequence...\n");
    Event *event = eventTree_getRootEvent(flower_getEventTree(end_getFlower(end)));
    Name eventName = event_getName(event);
    MetaSequence *metaSequence;
    int32_t start = 1; 

    reference_walkDown(end, refseq);
    assert(strlen(refseq->string) == refseq->length);
    metaSequence = metaSequence_construct(start, strlen(refseq->string), refseq->string, refseq->header, eventName, cactusDisk);
    return metaSequence;
}

void constructReferenceSequence(MetaSequence *metaSequence, Flower *flower){
    /*
     *Attach MetaSequence to cactus flowers
     */
    assert(flower != NULL);
    //add reference sequence to current flower
    sequence_construct(metaSequence, flower);

    //recursively add reference sequence to lower-level flowers
    Flower_GroupIterator *it = flower_getGroupIterator(flower);
    Group *group;
    while((group = flower_getNextGroup(it))!= NULL){
        Flower *nestedFlower = group_getNestedFlower(group);
	if(nestedFlower != NULL){
	    constructReferenceSequence(metaSequence, nestedFlower); 
	}
    }
    flower_destructGroupIterator(it);
}

char *getChromName(char *name, int num){
    /*
     *Add the order of the pseudochromosome to the reference sequence name. E.g 'reference.chr1'
     */
    char chrom[5];//assume there are less than 999 chroms! 
    sprintf(chrom, "%d", num);
    char *chromName = st_malloc(sizeof(char)*(strlen(name) + strlen(chrom) + strlen(".chr") + 1));
    strcpy(chromName, name);
    strcat(chromName, ".chr");
    strcat(chromName, chrom);
    return chromName;
}

/*
 *Adding the reference sequence of each pseudoChrom to cactus structure.
 *'header' is going to be set as the reference's name
 */
Flower *flower_addReferenceSequence(Flower *flower, CactusDisk *cactusDisk, char *header){
    Reference *reference = flower_getReference(flower);
    assert(reference != NULL);
    Reference_PseudoChromosomeIterator *it =
            reference_getPseudoChromosomeIterator(reference);
    PseudoChromosome *pseudoChromosome;
    int chromNum = 0;

    while ((pseudoChromosome = reference_getNextPseudoChromosome(it)) != NULL) {//Each pseudoChrom
        chromNum ++;
        char *chromHeader = getChromName(header, chromNum);
        End *end = pseudoChromosome_get5End(pseudoChromosome);
        assert(!end_isBlockEnd(end));

        int len = pseudoChromosome_getLength(end);
        //if(len == 0){continue;}
        assert(len != 0);
        ReferenceSequence *refseq = referenceSequence_construct(flower, chromHeader, len);
    
        st_logInfo("\nInitialize refseq: index %d, length %d, header *%s*, string *%s*\n", 
                   refseq->index, refseq->length, refseq->header ,refseq->string);
        
        //Construct the MetaSequence 
        MetaSequence *metaSequence = constructReferenceMetaSequence(end, cactusDisk, refseq);
        st_logInfo("Got metasequence: name *%s*, start %d, length %d, header *%s*\n", 
                    cactusMisc_nameToString(metaSequence_getName(metaSequence)), 
                    metaSequence_getStart(metaSequence), metaSequence_getLength(metaSequence), 
                    metaSequence_getHeader(metaSequence));

        st_logInfo("\nConstructing reference sequence...\n");
        constructReferenceSequence(metaSequence, flower);
        st_logInfo("Constructed reference sequence successfully.\n");

        //Add startStub (3' end)
        Sequence *sequence = getSequenceByHeader(flower, refseq->header);
	Cap *startcap = cap_construct2(end, refseq->index, true, sequence);
        cap_check(startcap);
	refseq->index ++;
        copyRefCapToLowerFlowers(startcap);

        //adding reference Segments to the blocks and creating inherited caps
        st_logInfo("Adding reference segments...\n");
        reference_walkDown(end, refseq); 
        st_logInfo("Added reference segments successfully.\n");
        
        //adding adjacencies to the reference caps 
        st_logInfo("Adding adjacencies to the reference caps...\n");
        addReferenceAdjacencies(end, chromHeader);
        st_logInfo("Added adjacencies to the reference caps successfully.\n");

        //write to Disk:
        cactusDisk_write(cactusDisk);
        
        //free memory:
        referenceSequence_destruct(refseq);
        free(chromHeader);
    }
    reference_destructPseudoChromosomeIterator(it);
    
    return flower;
}

void usage() {
    fprintf(stderr, "cactus_addReferenceSeq, version 0.0\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --name : Name of the reference sequence\n");
    fprintf(stderr,
            "-c --cactusDisk : The location of the flower disk directory (the databaseString)\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    //fprintf(stderr,
    //        "-d --flowerName : The name of the flower (the key in the database)\n");
}

int main(int argc, char *argv[]) {
    char * logLevelString = NULL;
    char * cactusDiskDatabaseString = NULL;
    char * flowerName = "0";
    char * name = NULL;

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
            //case 'd':
            //    flowerName = stString_copy(optarg);
            //    break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    //assert(flowerName != NULL);

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
    flower = flower_addReferenceSequence(flower, cactusDisk, name);
    //Flower *flower_addReferenceSequence(Flower *flower, CactusDisk *cactusDisk, char *header){

    //test(flower);
    st_logInfo("Added the reference sequence in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
