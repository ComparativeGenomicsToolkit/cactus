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

#define LINEMAXSIZE 3000

/*
 * nknguyen@soe.ucsc.edu
 * Jan 25 2011  
 ************************
 * This script finds homologs of each inputing species and return an xml tree
 * with found information. The steps involved are:
 *
 * 1/ Reads in a list of genes (bed format) and sorts them by chromosome,
 * then start coordinate, all relatively to the positive (forward) strand.
 * (This sorting helps to reduce the time traversing cactus graph to get to
 * the region that ovelaps with each gene).
 *
 * For each gene:
 * 2/ Maps exons of that gene to the sequence thread of the ref-species in cactus graph.
 * (Traversing the refspc's thread and recording the ends of blocks that ovelap with each Exon.
 * After mapping, each Gene contains a list of Exons; each Exon contains a list of Ends;
 * each End is the start of the block that overlaps with that Exon. 
 * Call these Genes ref-Genes, Exons ref-Exons, Ends ref-Ends
 *
 * 3/ Then, for each species (called currSpc): construct a list of MappedExons, in which
 * each MappedExon is the local maximal thread of currSpc that maps to ref-species with no gaps
 * or short (< 10% of the ref-exon's length) indels of multiple of three. Whenever there is a 
 * non-triplet gap, long gap, or other breaks (repeated segment (duplication), new exon, rearrangement...), 
 * start a new MappedExon. The final MappedExons are in the coorinate order of currSpc.
 *
 * The procedure to construct the MappedExon list is as followed:
 * Intialized mappedExonList = empty
 * 3.1/ For each ref-Exon:
 *        For each ref-End, find all the currSpc's threads in the corresponding Block.
 *          For each currSpcThread, go through the current mappedExonList and try to merge it with
 *          some existing mappedExon in the list. The requirement for merging is that there is NO
 *          insertion / deletion between the two segments (one of the mappedExon, and one is
 *          the currSpcThread). If cannot merge to any element of the list, then just append it to end of the 
 *          list as a new element.
 * 3.2/ After step 3.1, we now have a list of scrambled (unsorted) local maximal threads of currSpc that map to
 *      current gene's exons. We then sort these mappedExons in the order of the currSpc coordinates.
 * 3.3/ Now go through the list again and merge any two neighboring exons that has short indels of multiple of three.
 * 3.4/ Now look at all the mappedExons of each Exon, and merge them if the total indel size is of multiple of three.
 *
 * 4/CheckCompleteHomologyForward
 *   CheckCompleteHomologyBackward
 *   CheckPartialHomologyForward
 *   CheckPartialHomologyBackward
 *
 * Notes: 
 * Sequences of the reference species must have name that matches to:
 * "species.chr#". For example: "hg19.chr8"
 */

//============================ BED =====================================
struct bed{
    struct bed *next;
    char *chrom;
    int32_t chromStart;
    int32_t chromEnd;
    char *name;

    int score;
    //char strand[2];
    char *strand;
    int32_t thickStart;
    int32_t thickEnd;
    int32_t itemRgb;
    int32_t blockCount;
    struct IntList *blockSizes;
    struct IntList *chromStarts;
};

struct bed *constructbed(){ 
    struct bed *bed = st_malloc(sizeof(struct bed));
    bed->next = NULL;
    return bed;
}

struct bed *bedLoadAll(char *fileName){
    struct bed *list = NULL;
    struct bed *currbed = NULL;
    struct bed *prevbed = NULL;

    FILE *fp;
    char *line = st_malloc(sizeof(char)*LINEMAXSIZE);
    char *lend;
    struct List *lineList;

    fp = fopen(fileName, "r");
    assert(fp != NULL);
    
    while( fgets(line, LINEMAXSIZE, fp) != NULL ){//each line
        //catch if not enough buffer
        if ( (lend = strstr(line, "\n")) == NULL ){
            fprintf(stderr, "Input line is too long for buffer: \n*%s*\n", line);
	    exit(0);
        }else{//trim \n
            *lend = '\0';
        }
 
	//skip header line
	if(strstr(line, "track") == line){ continue; }

        lineList = splitString(line, "\t");

	//empty lines, skip
	if(lineList->length == 0){ continue; }

        if(lineList->length < 12){
	    fprintf(stderr, "Wrong input format. Need at least 12 fields for each bed\n");
	    exit(0);
        }else{
	    currbed = constructbed(); 
            currbed->chrom = stString_copy(lineList->list[0]);
            assert( sscanf(lineList->list[1], "%d", &(currbed->chromStart)) == 1);
            assert( sscanf(lineList->list[2], "%d", &(currbed->chromEnd)) == 1);
            currbed->name = stString_copy(lineList->list[3]);
            assert( sscanf(lineList->list[4], "%d", &(currbed->score)) == 1);
            currbed->strand = stString_copy(lineList->list[5]);
            assert( sscanf(lineList->list[6], "%d", &(currbed->thickStart)) == 1);
            assert( sscanf(lineList->list[7], "%d", &(currbed->thickEnd)) == 1);
            assert( sscanf(lineList->list[8], "%d", &(currbed->itemRgb)) == 1);
            assert( sscanf(lineList->list[9], "%d", &(currbed->blockCount)) == 1);
            currbed->blockSizes = splitIntString(stString_copy(lineList->list[10]), ",");
	    assert(currbed->blockSizes->length == currbed->blockCount);
            currbed->chromStarts = splitIntString(stString_copy(lineList->list[11]), ",");
	    assert(currbed->chromStarts->length == currbed->blockCount);
            //destructList(lineList);
	    if(list == NULL){ 
	        list = currbed; 
		prevbed = currbed;
	    }else{
                prevbed->next = currbed;
		prevbed = currbed;
	    }
	}
        st_logInfo("Gene: %s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\n", 
                    currbed->chrom, currbed->chromStart, currbed->chromEnd, currbed->name, 
                    currbed->score, currbed->strand, currbed->thickStart, currbed->thickEnd, 
                    currbed->itemRgb, currbed->blockCount);
    }

    fclose(fp);
    return list;
}
//============================ END OF BED ==========================

//=========================== GLOBAL STRUCTURES ===============
struct MappedExon{
    int exon;
    int blockIndex;
    int status; //status = 1 if complete, 0 if incomplete, 3 if incomplete because of missing data
    char *chr;
    int32_t start;
    int32_t length;
    int32_t exonStart;//block start relatively to the exon: exonStart = blockStart (5'end) - exonStart
    int32_t exonEnd;//block end relatively to the exon: exonEnd = exonEnd - (blockStart + blockLength)
    int32_t exonLength;
    int32_t exonTotalLength;
    int32_t insBases;
    int32_t delBases;
    char strand; 
    char *seq;
    bool isLeftStub;
    bool isRightStub;
    struct MappedExon *prev;
    struct MappedExon *next;
};

struct Gene{
    char *name;
    char *species;
    char strand;
    int exonCount;
    struct List *exons;//list of Threads
    struct IntList *introns;//list of boolean indicating whether each intron has duplication or not
};

struct Exon{
    int id;
    int32_t start;
    int32_t end;
    int hasDup;
    struct List *ends;
};

struct MappedEnd{
    Cap *cap;
    End *end;
    int32_t exonStart; //block start relatively to the exon: = blockStart (5'end) - exonStart
    int32_t exonEnd; //block start relatively to the exon: = exonEnd - (blockStart + blockLength)
};

struct CapList{
    Cap *cap;
    int32_t endSeqCoor;
    struct CapList *next;
    struct CapList *prev;
};

//============= CONSTRUCTOR ==============
struct MappedExon *constructMappedExon(int exon, int32_t exonStart, char *chr, 
    int32_t start, int32_t length, char strand, int blockIndex, char *seq, 
    int32_t exonEnd, int32_t exonLength, int32_t exonTotalLength){
    
    struct MappedExon *mexon = st_malloc(sizeof(struct MappedExon));
    mexon->exon = exon;
    mexon->exonStart = exonStart;
    mexon->exonEnd = exonEnd;
    mexon->chr = chr;
    mexon->start = start;
    mexon->length = length;
    mexon->exonLength = exonLength;
    mexon->exonTotalLength = exonTotalLength;
    mexon->strand = strand;
    mexon->blockIndex = blockIndex;
    mexon->seq = seq;
    mexon->insBases = 0;
    mexon->delBases = 0;
    mexon->isLeftStub = false;
    mexon->isRightStub = false;
    mexon->next = NULL;
    mexon->prev = NULL;
    return mexon;
}

struct Exon *constructExon(int id){ 
    struct Exon *exon = st_malloc(sizeof(struct Exon));
    exon->ends = constructEmptyList(0, free);
    exon->hasDup = 0;
    exon->id = id;
    return exon;
}

struct Gene *constructGene(char *name, char strand, int exonCount){
    struct Gene *gene = st_malloc(sizeof(struct Gene));
    gene->name = name;
    gene->exonCount = exonCount;
    gene->exons = constructEmptyList(0, free);
    gene->introns = constructEmptyIntList(0);
    gene->strand = strand;
    return gene;
}

struct MappedEnd *constructMappedEnd(Cap *cap, End *end, int32_t exonStart, int32_t exonEnd){
    struct MappedEnd *me = st_malloc(sizeof(struct MappedEnd));
    me->end = end;
    me->cap = cap;
    me->exonStart = exonStart;
    me->exonEnd = exonEnd;
    return me;
}

struct CapList *constructCapList(){
    struct CapList *caplist = st_malloc(sizeof(struct CapList));
    caplist->cap = NULL;
    caplist->next = NULL;
    caplist->prev = NULL;
    return caplist;
}

//=========== DESTRUCT FUNCTIONS ============
void destructExon(struct Exon *thread){
    destructList(thread->ends);
    free(thread);
}

void destructGene(struct Gene *gene){
    for(int i=0; i< gene->exons->length; i++){
        destructExon(gene->exons->list[i]);
    }
    free(gene->exons);
    destructIntList(gene->introns);
    free(gene);
}


//============================= UTILS FUNCTIONS ========================
bool visitedString(struct List *list, char *name){
    for(int i=0; i< list->length; i++){
        if(strcmp(list->list[i], name) == 0){
            return true;
        }
    }
    return false;
}

char *getChrom(char *name){
    assert(name != NULL);
    struct List *list = splitString(name, ".");
    if(list->length >=2){
        return stString_copy(list->list[1]);
    }else{
        return stString_copy(name);
    }
}

int32_t getSeqStartFromName(char *name){
    int32_t start = 0;
    struct List *list = splitString(name, ".");
    if(list->length >= 4){
        sscanf(list->list[3], "%d", &start);
    }
    return start;
}

int32_t getSeqEndFromName(char *name){
    int32_t start = 0;
    int32_t length = 0;
    struct List *list = splitString(name, ".");
    assert(list->length >= 5);
    sscanf(list->list[3], "%d", &start);
    sscanf(list->list[4], "%d", &length);
    return start + length;
}

int32_t mapCapCoor2(Cap *cap, struct List *list){
    int32_t start = 0;
    if(list->length >=4){
        sscanf(list->list[3], "%d", &start);
    }
    return cap_getCoordinate(cap) - 2 + start;
}

int32_t mapCapCoor(Cap *cap){
    //st_logInfo("mapCapCoor\n");
    char *seq = cap_getSequenceName(cap);
    struct List *list = splitString(seq, ".");
    int32_t coor = mapCapCoor2(cap, list);
    //destructList(list);
    return coor;
}

int capCmp(Cap *cap1, Cap *cap2){
    int cmp;
    char *seq1 = cap_getSequenceName(cap1);
    struct List *list1 = splitString(seq1, "."); 
    char *seq2 = cap_getSequenceName(cap2);
    struct List *list2 = splitString(seq2, ".");
    st_logInfo("capCmp: %s, %s\n", seq1, seq2);

    assert(list1->length >= 2 && list2->length >= 2);//name must have at least species.chr info
    assert(strcmp(list1->list[0], list2->list[0]) == 0);//from the same species
    
    char *chr1 = list1->list[1];
    char *chr2 = list2->list[1];
    cmp = strcmp(chr1, chr2);//compare two chromosomes
    if(cmp == 0){//if same chromosome, compare the coordinates
        int32_t coor1 = mapCapCoor2(cap1, list1);
        int32_t coor2 = mapCapCoor2(cap2, list2);
        if(coor1 < coor2){
            cmp = -1;
        }else if(coor1 == coor2){
            cmp = 0;
        }else{
            cmp = 1;
        }
    }
    //destructList(list1);
    //destructList(list2);
    return cmp;
}

struct CapList *insertCapList(struct CapList *caplist, Cap *cap){
    assert(cap != NULL);
    st_logInfo("insertCapList, cap: %s\n", cactusMisc_nameToString(cap_getName(cap)));
    //insert cap to a list of cap so that the caps are in order from small to large coordinates
    struct CapList *currcapL = caplist;
    struct CapList *new = constructCapList();
    new->cap = cap;

    int32_t seqLen = getSeqLength(end_getFlower(cap_getEnd(cap)), cap_getSequenceName(cap));
    new->endSeqCoor = mapCapCoor(cap) + seqLen + 1;    
    //st_logInfo("constructed 'NEW', capCoor = %d, mapCoor: %d, seqLen = %d, endSeqCoor = %d\n", cap_getCoordinate(cap), mapCapCoor(cap), seqLen, new->endSeqCoor);

    Cap *currcap;

    if(caplist->cap == NULL){//first cap
        caplist = new;
    }else{
        while(currcapL != NULL){
            currcap = currcapL->cap;
            assert(currcap != NULL);
            int cmp = capCmp(currcap, cap);
            if(cmp < 0){
                if(currcapL->next == NULL){//reach the end of list
                    currcapL->next = new;
                    new->prev = currcapL;
                    break;
                }else{
                    currcapL = currcapL->next;
                }
            }else if(cmp > 0){
                if(currcapL->prev == NULL){//currcapL is the first node in the linked list
                    new->next = currcapL;
                    currcapL->prev = new;
                    caplist = new;
                }else{
                    currcapL->prev->next = new;
                    new->prev = currcapL->prev;
                    new->next = currcapL;
                    currcapL->prev = new;
                }
                break;
            }else{
                break;
            }
        }
    }
    return caplist;
}

//=================== GET START OF THREAD (3' STUB) ====================================
struct CapList *flower_getThreadStart(Flower *flower, char *name){// from Flower = 0
    /*
     *Get 3' end Stub (the start) of the sequence by its name
     */
    struct CapList *caplist = constructCapList();
    Cap *cap;
    Flower_CapIterator *capIterator = flower_getCapIterator(flower);
    while((cap= flower_getNextCap(capIterator)) != NULL){
        if(isStubCap(cap)){//dead end or inherited end
            if( !cap_getStrand(cap) ){//convert cap to + strand
                cap = cap_getReverse(cap);
            }
            if( !cap_getSide(cap) ){//3' 
                Sequence *sequence = cap_getSequence(cap);
                if(sequence == NULL){continue;}
                char *sequenceHeader = formatSequenceHeader(sequence);
                if(strstr(sequenceHeader, name) != NULL){
                    caplist = insertCapList(caplist, cap);
                    /*char *strand = cap_getStrand(cap) ? "+" : "-";
                    char *side = cap_getSide(cap) ? "5" : "3";
                    st_logInfo("Found a thread start: %d, %s, %s\n", mapCapCoor(cap), strand, side);*/
                }
                free(sequenceHeader);
            }
        }
    }
    flower_destructCapIterator(capIterator);
    return caplist;
}

int block_checkDuplication(Block *block){
    //return 1 if has duplication, defined as at least two segments come 
    //from one species (regardless which chromosome of that species)
    //fprintf(stderr, "block_checkDup\t");
    int hasDup = 0;
    Segment *segment;
    struct List *spcList = constructEmptyList(0, free); //list of visited species
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    while((segment = block_getNext(it)) != NULL){
        Sequence *sequence = segment_getSequence(segment);
        //if(sequence == NULL){continue;}
        char *seqHeader = formatSequenceHeader(sequence);
        struct List *headerList = splitString(seqHeader, ".");
        assert(headerList->length > 0);
        if(!visitedString(spcList, headerList->list[0])){//has not visited this species yet
            listAppend(spcList, headerList->list[0]);
            //fprintf(stderr, "-added\t");
        }else{//has duplications
            //fprintf(stderr, "-dup!");
            hasDup = 1;
        }
        //destructList(headerList);
        free(seqHeader);
        if(hasDup == 1){break;}
    }
    //destructList(spcList);
    block_destructInstanceIterator(it);
    //fprintf(stderr, "\n");
    return hasDup;
}

//============== 
//Traverse the thread from inputed 'currCap' to find the Cap that is closest
//and upstream (relative to the positive strand) to inputed 'start'

Cap *walkDown(Cap *cap, int32_t start);

Cap *walkUp(Cap *cap, int32_t start){
    assert(cap != NULL);
    Segment *segment = cap_getSegment(cap);
    if (segment != NULL) {
        Cap *otherCap = cap_getOtherSegmentCap(cap);
        //int32_t otherCoor = cap_getCoordinate(otherCap); 
        int32_t otherCoor = mapCapCoor(otherCap);
        if(otherCoor < start){//still upstream of 'start', walkDown (or cross)
            cap = walkDown(otherCap, start);
        }
    } else {
        //assert(end_isAttached(cap_getEnd(cap)));
        Group *parentGroup = flower_getParentGroup(end_getFlower(cap_getEnd(cap)));
        if (parentGroup != NULL) {
            Cap *upperCap = flower_getCap(group_getFlower(parentGroup), cap_getName(cap));
            assert(upperCap != NULL);

            //Make sure always following the same strand
            if(cap_getStrand(cap) != cap_getStrand(upperCap)){
                upperCap = cap_getReverse(upperCap);
            }
            cap = walkUp(upperCap, start);
        }
    }
    return cap;
}

Cap *walkDown(Cap *cap, int32_t start){//cap is either stub or the 'otherSegmentCap'
    assert(cap != NULL);
    //int32_t currcoor = cap_getCoordinate(cap);
    int32_t currcoor = mapCapCoor(cap);
    if(start <= currcoor){//current cap is already downstream of (or is at) 'start'
        return cap;
    }

    Cap *adjCap = cap_getAdjacency(cap);
    //int32_t adjCoor = cap_getCoordinate(adjCap);
    int32_t adjCoor = mapCapCoor(adjCap);

    if(adjCoor <= start){//still upstream of 'coor', walk across
        cap = walkUp(adjCap, start);
    }else{//currcoor < start < adjCoor
        Group *group = end_getGroup(cap_getEnd(cap));
        assert(!group_isLeaf(group));//there must be lower level for currcoor < start < adjCoor
        Cap *lowerCap = flower_getCap(group_getNestedFlower(group), cap_getName(cap));
        if(cap_getStrand(cap) != cap_getStrand(lowerCap)){
            lowerCap = cap_getReverse(lowerCap);
        }
        cap = walkDown(lowerCap, start);
    }
    return cap;
}
//======================



//================ MAP GENE TO OTHER SPECIES (FINDING HOMOLOGS) ====================
bool intListContains(struct IntList *list, int32_t num){
    for(int i = 0; i < list->length; i++){
        if(list->list[i] == num){
            return true;
        }
    }
    return false;
}

void extendMappedExon(struct List *mappedExons, Cap *cap, char *chr, int currBlockIndex, 
                      int exon, int32_t exonStart, char *seqname, int exonEnd, int32_t exonTotalLength){
    char strand = cap_getStrand(cap) ? '+' : '-';
    struct MappedExon *mexon = NULL;
    int32_t start = mapCapCoor(cap);
    int32_t segmentLen = segment_getLength(cap_getSegment(cap));
    int i = -1;
    for(i=0; i< mappedExons->length; i++){
        mexon = mappedExons->list[i];
        //must have the same strand, same chr to be able to extend
        if(strcmp(chr, mexon->chr) == 0 && strand == mexon->strand && currBlockIndex == mexon->blockIndex + 1){
            if( (strand == '+' && mexon->start + mexon->exonLength == start) || 
                (strand == '-' && start + mexon->length == mexon->start)){
                mexon->length += segmentLen;
                mexon->exonLength += segmentLen;
                mexon->exonEnd = exonEnd;
                mexon->blockIndex = currBlockIndex;
                break;
            }
        }
    }
    //Could not extend to any of the previous mappedExon, create a new one
    if(i == mappedExons->length){
        mexon = constructMappedExon(exon, exonStart, chr, start, segmentLen, strand, currBlockIndex, 
                                    seqname, exonEnd, segmentLen, exonTotalLength);
        //mexon = constructMappedExon(exon, exonStart, chr, start, segmentLen, strand, currBlockIndex, seqname, lastExonBlockIndex);
        listAppend(mappedExons, mexon);
    }

    //add Stub information:
    if ( mexon != NULL && isStubCap(cap_getAdjacency(cap)) ){
        mexon->isLeftStub = true;
    }
    if ( mexon != NULL && isStubCap(cap_getAdjacency(cap_getOtherSegmentCap(cap))) ){
        mexon->isRightStub = true;
    }

    return; 
}

void getCoverageStatus(struct List *mappedExons){
    //check to see if each mapped-exon fully covers the exon (update status to 1) or not (status =0)
    //lastIndex is the index of the last block that mapped to the exon.
    struct MappedExon *mexon;
    for(int i=0; i< mappedExons->length; i++){
        mexon = mappedExons->list[i];
        //if(mexon->exonStart <= 0 && mexon->blockIndex == mexon->lastExonBlockIndex){//complete coverage
        if(mexon->exonStart <= 0 && mexon->exonEnd <= 0){//complete coverage
            mexon->status = 1;
        }else{
            mexon->status = 0;
        }
    }
    return;
}

struct List *exon_getCoverage(struct Exon *exon, char *species){
    /*
     * Finding maximal non-gap concatenated segments in 'chrom' that mapped to 'exon'
     */
    if (exon->ends == NULL || exon->ends->length == 0){return NULL;}
    struct List *mappedExons = constructEmptyList(0, free);
    int32_t exonTotalLength = exon->end - exon->start;
    for(int i=0; i< exon->ends->length; i++){//each cactusBlock that overlap with the exon, from left to right
        struct MappedEnd *mend = exon->ends->list[i];
        End_InstanceIterator *it = end_getInstanceIterator(mend->end);
        Cap *cap;
        while( (cap = end_getNext(it)) != NULL ){
            char *seqname = cap_getSequenceName(cap);
            if(seqname != NULL && strstr(seqname, species) != NULL){//found segment of interested species
                extendMappedExon(mappedExons, cap, getChrom(seqname), i, exon->id, 
                                 mend->exonStart, seqname, mend->exonEnd, exonTotalLength);
            }
        }
        end_destructInstanceIterator(it);
    }
    getCoverageStatus(mappedExons);
    return mappedExons;
}

int mappedExon_cmp(struct MappedExon *me1, struct MappedExon *me2){
    assert(me1 != NULL && me2 != NULL);
    int chrCmp = strcmp(me1->chr, me2->chr);
    if(chrCmp != 0){
        return chrCmp;
    }else if(me1->start < me2->start){
        return -1;
    }else if (me1->start == me2->start){
        if(me1->strand == '+'){
            if(me1->exon < me2->exon){
                return -1;
            }else if(me1->exon == me2->exon){
                return 0;
                /*if(me1->exonStart < me2->exonStart){
                    return -1;
                }else if(me1->exonStart == me2->exonStart){
                    return 0;
                }else{
                    return 1;
                }*/
            }else{
                return 1;
            }
        }else{
            if(me1->exon > me2->exon){
                return -1;
            }else if(me1->exon == me2->exon){
                return 0;
            }else{
                return 1;
            }
        }
    }else{//same block
        return 1;
    }
}

struct MappedExon *insertMappedExon(struct MappedExon *root,  struct MappedExon *me){
    assert(root != NULL && me != NULL);
    st_logInfo("insertMappedExon, chr: %s, start: %d, length: %d\n", me->chr, me->start, me->length);
    //insert cap to a list of cap so that the caps are in order from small to large coordinates
    struct MappedExon *currMexon = root;

    while(currMexon != NULL){
        int cmp = mappedExon_cmp(currMexon, me);
        if(cmp < 0){
            //st_logInfo("currME < me\n");
            if(currMexon->next == NULL){//reach the end of list
                currMexon->next = me;
                me->prev = currMexon;
                assert(me->next == NULL);
                break;
            }else{
                currMexon = currMexon->next;
            }
        }else if(cmp > 0){
            //st_logInfo("currME > me\n");
            if(currMexon->prev == NULL){//currMexon is the first node in the linked list
                me->next = currMexon;
                currMexon->prev = me;
                root = me;
                st_logInfo("root: %d, next: %d, nextPrev: %d\n", root->start, root->next->start, root->next->prev->start);
            }else{
                currMexon->prev->next = me;
                me->prev = currMexon->prev;
                me->next = currMexon;
                currMexon->prev = me;
            }
            break;
        }else{
            //st_logInfo("TWO MAPPED EXONS ARE EQUAL!!! POTENTIAL ERROR!!\n");
            break;
        }
    }
    return root;
}

//struct MappedExon *merge2MappedExons(struct MappedExon *curr, struct MappedExon *next, char strand, bool missingData){
struct MappedExon *merge2MappedExons(struct MappedExon *curr, struct MappedExon *next, char strand, int32_t insBases, int32_t delBases){
    assert(curr != NULL && next != NULL);
    //assert(curr != NULL && next != NULL && curr->next == next);
    st_logInfo("Merge %d With %d; insBases = %d, delBases = %d\n", curr->start, next->start, insBases, delBases);
    float maxIndel = 0;
    if(strand == '+'){//merge next to curr
        curr->length = next->start + next->length - curr->start;
        curr->exonLength = next->exonStart + next->exonLength - curr->exonStart;
        curr->exonEnd = next->exonEnd;
        curr->blockIndex = next->blockIndex;
        curr->isRightStub = next->isRightStub;
        
        curr->delBases += delBases + next->delBases;
        curr->insBases += insBases + next->insBases;
        maxIndel = curr->exonTotalLength * 0.1;
        st_logInfo("maxIndel: %f\n", maxIndel);
        if(curr->exonEnd <= 0 && curr->exonStart <=0 && curr->delBases <= maxIndel && curr->insBases <= maxIndel){
            curr->status = 1;
        }
        curr->next = next->next;
        if(curr->next != NULL){
            st_logInfo("add the reverse edge\n");
            (curr->next)->prev = curr;
        }else{
            st_logInfo("Next is NULL\n");
        }
        return curr;
    }else{//merge curr to next
        next->length = next->start - (curr->start - curr->length);
        next->exonLength = curr->exonStart + curr->exonLength - next->exonStart;
        next->blockIndex = curr->blockIndex;
        next->exonEnd = curr->exonEnd;
        next->isRightStub = curr->isRightStub;
        
        next->delBases += delBases + curr->delBases;
        next->insBases += insBases + curr->insBases;
        
        maxIndel = next->exonTotalLength * 0.1;
        st_logInfo("maxIndel: %f\n", maxIndel);
        if(next->exonEnd <= 0 && next->exonStart <=0 && next->delBases <= maxIndel && next->insBases <= maxIndel){
            next->status = 1;
        }
        next->prev = curr->prev;
        if(next->prev != NULL){
            st_logInfo("add the reverse edge\n");
            (next->prev)->next = next;
        }else{
            st_logInfo("Prev is NULL");
        }
        return next;
    }
}

struct MappedExon *mergeIncompleteMappedExons(struct MappedExon *mexons){
    assert(mexons != NULL);
    struct MappedExon *curr = mexons;
    struct MappedExon *next;
    struct MappedExon *root = mexons;
    int32_t delBases, insBases;
    
    //going from left to right, merge two neighbor incomplete mappedExons if
    //the incompleteness is caused by deletion/insertion of multiple of 3s
    while(curr->next != NULL){
        //st_logInfo("CURR: exonStart %d, CURRnext %d\n", curr->exonStart, curr->next->exonStart);
        next = curr->next;
        delBases = -1;
        insBases = -1;
        //current and next mappedExon are inComplete and on the same chromosome
        if(curr->status != 1 && next->status !=1 && strcmp(curr->chr, next->chr) == 0 &&
           curr->strand == next->strand && curr->exon == next->exon){
            /*st_logInfo("\tcurrStand: %c == nextStrand: %c, same chrom\n", curr->strand, next->strand);
            st_logInfo("\tif insertion case: curr->exonStart (%d) + curr->exonLength (%d) == next->exonStart (%d)\n", 
                         curr->exonStart, curr->exonLength, next->exonStart);
            st_logInfo("\tAlso if insertion: nextStart (%d)- (currStart(%d) + currLength (%d)) = %d is mul3\n", 
                         next->start, curr->start, curr->length, (next->start - curr->start - curr->length));
            */
            if(curr->strand == '+'){//positive strand
                insBases = next->start - (curr->start + curr->length);
                delBases = next->exonStart - (curr->exonStart + curr->exonLength);

                //if (next->start - curr->start ==  curr->length){ //deletion in current species
                //}else if(curr->exonStart + curr->exonLength == next->exonStart){//insertion
                if(delBases %3 == 0 && delBases <= curr->exonTotalLength*0.1 &&
                   insBases %3 == 0 && insBases <= curr->exonTotalLength*0.1){
                    curr = merge2MappedExons(curr, next, curr->strand, insBases, delBases);
                    if(curr->prev == NULL){ root = curr; }
                }else{
                    curr = curr->next;
                }
            }else{ //Negative strand
                delBases = curr->exonStart - (next->exonStart + next->exonLength);
                insBases = (next->start - next->length) - curr->start; 
                //if (next->start - curr->start == next->length){//deletion in current species
                //}else if (next->exonStart + next->exonLength == curr->exonStart){ //insertion in current species
                //if(indelBases %3 == 0 && abs(indelBases) <= curr->exonTotalLength*0.1){
                if(delBases %3 == 0 && delBases <= curr->exonTotalLength*0.1 &&
                   insBases %3 == 0 && insBases <= curr->exonTotalLength*0.1){
                    curr = merge2MappedExons(curr, next, curr->strand, insBases, delBases);
                    if(curr->prev == NULL){ root = curr; }
                }else{
                    curr = curr->next;
                }
            }
        }else{
            curr = curr->next;
        }
    }
    return root;
}

bool checkStrandConsistency(struct MappedExon *start, struct MappedExon *end){
    assert(start != NULL && end != NULL);
    struct MappedExon *curr = start;
    char strand = curr->strand;
    while(curr != end->next && curr != NULL){
        if(strand != curr->strand){
            return false; 
        }
        curr = curr->next;
    }
    return true;
}

struct MappedExon *mergeIncompleteMappedExons2(struct MappedExon *mexons){
    /*
     * This step is done after mergeIncompleMappedExon (which merged any two neighboring
     * incomplete mappedExon if the indel between them has a length of multiple of 3).
     * This second merge checks if the total indel length of all the continous mappedExons 
     * that mapped to an exon is a mult of 3. If so, merge all of them together.
     */
    assert(mexons != NULL);
    struct MappedExon *curr = mexons;
    struct MappedExon *merge = NULL;
    struct MappedExon *start = NULL;
    struct MappedExon *end = NULL;
    int currExon = -1;
    int indels = -1;
    int currIndel;
    int insBases = 0;
    int delBases = 0;
    bool repeatlast = false;
    bool hasLargeIndels = false;
    while(curr != NULL){
        st_logInfo("CURR: exonStart %d\n", curr->exonStart);
        if(curr->exon != currExon){//end prevExon, starts a new exon
            //prevExon contains incomplete mappedExons, but the total 
            //indels length is mult3, so we merge all the pieces together
            st_logInfo("exon: %d != prevExon: %d\n", curr->exon, currExon);
            if(indels%3 == 0 && !hasLargeIndels){
                st_logInfo("\tindels %d is mult3\n", indels);
                end = repeatlast ? curr : curr->prev;
                assert(end != NULL);
                //only merge if from start to end is the same contig. 
                //Otherwise there exists missing data and can't tell in that case.
                if(start != end && strcmp(start->seq, end->seq)==0 && checkStrandConsistency(start, end)){
                    st_logInfo("\t\tstart != end, same seq, same strand, merge\n");
                    merge = merge2MappedExons(start, end, start->strand, insBases, delBases);
                    st_logInfo("MERGE2, start = %d, end = %d\n", start->start, end->start);
                    if(merge->prev == NULL){ mexons = merge;}
                }else{
                    st_logInfo("\t\tstart == end OR diff seq: %s and %s OR diff strand: %c, %c, NOT merge\n", start->seq, end->seq, start->strand, end->strand);
                }
            }
            //reset the variables
            if(curr->status == 1){//already complete
                indels = -1;
                currExon = -1;
            }else{//the first mappedExon of the currExon
                start = curr;
                indels = 0;
                hasLargeIndels = false;
                currExon = curr->exon;
            }
            insBases = 0;
            delBases = 0;
        }else{//same exon
            st_logInfo("exon: %d == prevExon: %d\n", curr->exon, currExon);
            if( (curr->strand == '+' && curr->exonStart <= curr->prev->exonStart) ||
                (curr->strand == '-' && curr->exonStart >= curr->prev->exonStart) ||
                curr->strand != curr->prev->strand ){//reset
                st_logInfo("\treset\n");
                currExon = -1;
                indels = -1;
                insBases = 0;
                delBases = 0;
                continue;
            }else{//update the indels
                st_logInfo("\tupdate indels\n");
                if(curr->strand == '+'){//positive strand
                    if( curr->start - curr->prev->start ==  curr->prev->length ){ //deletion in current species
                        currIndel = (-1)*(curr->exonStart - (curr->prev->exonStart + curr->prev->exonLength));//deletion is (-)
                        delBases += abs(currIndel);
                    }else{
                        //assert(curr->prev->exonStart + curr->prev->length == curr->exonStart || 
                        //       strcmp(curr->prev->seq, curr->seq) != 0);//insertion in current species
                        currIndel = curr->start - (curr->prev->start + curr->prev->length);
                        insBases += currIndel;
                    }
                }else{//negative strand
                    if(curr->start - curr->prev->start == curr->length ){ //deletion in current species
                        currIndel = (-1)*(curr->prev->exonStart - (curr->exonStart + curr->exonLength));
                        delBases = abs(currIndel);
                    }else{
                        //assert(curr->exonStart + curr->length == curr->prev->exonStart ||
                        //       strcmp(curr->prev->seq, curr->seq) != 0);//insertion in current species
                        currIndel = (curr->start - curr->length) - curr->prev->start;
                        insBases += currIndel;
                    }
                }
                //reset if current indel (gap) size is larger than 10% of the exon Length
                if( abs(currIndel) > curr->exonTotalLength*0.1 ){
                    hasLargeIndels = true;
                }
                indels += currIndel;
                
                if(curr->next == NULL){//curr is the last MappedExon, but still need to be merged, so iterate one more time
                    st_logInfo("\tcurr->next is NULL, revisit curr\n");
                    currExon = -1; 
                    repeatlast = true;
                    continue;
                }
            }
        }
        curr = curr->next;
    }

    return mexons;
}

bool checkMissingExons(struct IntList *missingExons, int start, int end){
    //return true if all exon in [start, end] (start<= end, inclusive) in missingExons
    //otherwise return false
    assert(start <= end);
    for(int i= start; i<= end; i++ ){
        if(!intListContains(missingExons, i)){
            return false;
        }
    }
    return true;
}

bool gene_checkMissingData(struct MappedExon *mexons){
    //return true if there are more than one contigs come from the same chromosome
    //or if an exon got disrupted by missing data
    //LIMITATIONS: If a gene only contains one contig, and that contig doesn't interupt
    //any exons, but it starts in a middle of the gene and missed some exons at the begining
    //or end of the gene, this function will still return false. 
    //But in those cases, if a gene has partially complete homolog, then given out input data, 
    //there's a high chance that it's caused by missing data. So this missing data info can be
    //retrieved. If a gene is incomplete, then it doesn't matter if it has missingdata or not...
    struct MappedExon *curr = mexons;
    char *chr = curr->chr;
    char *seq = curr->seq;

    while(curr != NULL){
        if((curr->isLeftStub && curr->exonStart > 0) || (curr->isRightStub && curr->exonEnd > 0)){
            return true;
        }else if(strcmp(curr->chr, chr) == 0){//same chromosome
            if(strcmp(curr->seq, seq) != 0){//different contigs
                return true;
            }
        }else{
            chr = curr->chr;
        }
        seq = curr->seq;
        curr = curr->next;
    }
    return false;
}

struct IntList *gene_findMissingDataExons(struct MappedExon *mexons){
    struct IntList *mdexons = constructEmptyIntList(0);
    struct MappedExon *curr = mexons;
    while(curr != NULL){
        if( (curr->isRightStub && curr->exonEnd > 0) || (curr->isLeftStub && curr->exonStart > 0)){
            if(!intListContains(mdexons, curr->exon)){
                intListAppend(mdexons, curr->exon);
            }
        }
        curr = curr->next;
    }
    return mdexons;
}

void correctForAltSplice(struct MappedExon *mexons){
    struct MappedExon *curr = mexons;
    int32_t minExonLen;
    int32_t alnBases;
    while(curr != NULL){
        minExonLen = 0.9*curr->exonTotalLength;
        alnBases = curr->exonLength - curr->delBases;
        if(curr->status == 0 && alnBases >= minExonLen){
            if(curr->strand == '+'){
                if( (curr->prev == NULL || curr->exon != curr->prev->exon || 
                     (curr->exonStart < curr->prev->exonStart + curr->prev->exonLength)) &&
                    (curr->next == NULL || curr->exon != curr->next->exon ||
                     (curr->exonStart + curr->exonLength > curr->next->exonStart)) ){
                    curr->status = 1;
                }
            }else{
                if( (curr->prev == NULL || curr->exon != curr->prev->exon || 
                     (curr->prev->exonStart < curr->exonStart + curr->exonLength)) &&
                    (curr->next == NULL || curr->exon != curr->next->exon ||
                     (curr->next->exonStart + curr->next->exonLength > curr->exonStart)) ){
                    curr->status = 1;
                }
            }
        }
        curr = curr->next;
    }
    return;
}

bool checkCompleteHomologyForward(struct MappedExon *mexons, int exonCount){
    /*
     * A species is defined to have a complete homolog with the reference gene when:
     *  a/It contains a traversal of the reference gene's exons, in which the order and orientation of the reference exons are conserved.
     *    In another words, all the exon-adjacencies of the reference gene must be conserved in that species.
     *  b/ All the exons in the traversal are conserved.
     */
    st_logInfo("\nChecking complete homology forward:\n");
    int exon = -1;
    struct MappedExon *curr = mexons;

    char *chr = curr->chr;
    char strand = curr->strand;
    bool hasHomolog = NULL;

    //move from left to right, try to find a string of continous 0, 1, 2, ..., exonCount -1 with status == 1
    while(curr != NULL){
        st_logInfo("currExon: %d, start: %d, end: %d, exonStart: %d, exonEnd: %d, chr: %s\n", 
                    curr->exon, curr->start, curr->start + curr->length, curr->exonStart, curr->exonEnd, curr->chr);
        if (exon == exonCount -1 && curr->prev->status == 1 && strand == '+'){ //prevMapped exon is last exon, fully covered
            return true;
        }

        if( strcmp(curr->chr, chr) == 0 ){//same chromosome
            st_logInfo("\tcurrChr: %s == chr : %s\n", curr->chr, chr);
            if(curr->strand != strand){//switch strand, reset
                exon = -1;
                strand = curr->strand;
                continue;
            }else if(curr->exon < exon){//startOver at curr
                st_logInfo("\t\tcurrExon (%d) < exon (%d) \n", curr->exon, exon);
                exon = -1;
                continue;
            }else if (curr->exon == exon){//repeating exon
                st_logInfo("\t\tcurrExon (%d) = exon: (%d)\n", curr->exon, exon);
                if(curr->exon == 0 && curr->status == 1){//previous block is also exon 0, but since it's the first exon, we can start over the gene here
                    st_logInfo("\t\t\t exon == 0 and status = 1, set missingData to false\n");
                }else{//prevExon == currExon, and not the first exon <--- violate the rule. Stop. Reset.
                    st_logInfo("\t\t\treset\n");
                    exon = -1;
                }
            }else if (curr->exon == exon + 1){
                st_logInfo("\t\tcurrExon %d = exon%d + 1\n", curr->exon, exon);
                exon = (curr->status == 1) ? (exon + 1) : -1;
            }else{//curr->exon > exon + 1
                st_logInfo("\t\tcurrExon (%d ) > exon: (%d) + 1. Reset\n", curr->exon, exon);
                exon = -1;
            }
            curr = curr->next;
        }else{//diff chromosome
            exon = -1;
            chr = curr->chr;
        }
    }
   
    hasHomolog = (exon == exonCount - 1) ? true : false;
    return (hasHomolog && strand == '+');
}

bool checkCompleteHomologyBackward(struct MappedExon *mexons, int exonCount){
    /*
     * A species is defined to have a complete homolog with the reference gene when:
     *  a/It contains a traversal of the reference gene's exons, in which the order and orientation of the reference exons are conserved.
     *    In another words, all the exon-adjacencies of the reference gene must be conserved in that species.
     *  b/ All the exons in the traversal are conserved.
     */
    st_logInfo("\nChecking complete homology backward:\n");
    int exon = exonCount;
    struct MappedExon *curr = mexons;

    char *chr = curr->chr;
    char strand = curr->strand;
    bool hasHomolog = NULL;

    //move from left to right, try to find a string of continous exonCount -1, exonCount -2, ..., 2, 1, 0  with status == 1
    while(curr != NULL){
        if (exon == 0 && curr->prev->status == 1 && strand == '-'){ //prevMapped exon is exon 0, fully covered
            return true;
        }
        if( strcmp(curr->chr, chr) == 0 ){//same chromosome
            if(curr->strand != strand){//switch strand, reset
                exon = exonCount;
                strand = curr->strand;
                continue;
            }else if(curr->exon > exon){//startOver at curr
                exon = exonCount;
                continue;
            }else if (curr->exon == exon){//repeating exon
                //can start over if block is the last exon and status 1, otherwise reset
                if(curr->exon != exonCount -1 || curr->status != 1){
                    exon = exonCount;
                }
            }else if (curr->exon == exon - 1){//only move forward if current status is 1
                exon = (curr->status == 1) ? (exon - 1) : exonCount;
            }else{//curr->exon > exon + 1
                exon = exonCount;
            }
            curr = curr->next;
        }else{//diff chromosome
            exon = exonCount;
            chr = curr->chr;
        }
    }
   
    hasHomolog = (exon == 0) ? true : false;
    return (hasHomolog && strand == '-');
}

//bool checkPartialHomologyForward(struct MappedExon *mexons, struct IntList *missingExons, struct IntList *mdExons, int exonCount){
bool checkPartialHomologyForward(struct MappedExon *mexons, struct IntList *missingExons, int exonCount){
    /*
     * A 'partially-complete' homolog is similar to a 'complete' homolog, but with a few exceptions:
     *  a/ It is a traversal of the exons in the same order and orientation with the reference exons. 
     *     However, the exons in the traversal are allowed to have missing data.
     *  b/ Absent exons due to missing data are allowed
     *
     */
    st_logInfo("\nchecking partial homology forward:\n");
    int exon = -1;
    struct MappedExon *curr = mexons;

    char *chr = curr->chr;
    char strand = curr->strand;
    char *seq = curr->seq;
    bool hasHomolog = NULL;
    bool exonMissingData = false;
    
    //going from left to right search for sequence of 0, 1, 2 ... exonCount -1. Only move forward if the previous exon status is 1, 
    //or if there is missingData in previous exon. Whenever detect missing data within an exon, add it to the mdExons list (missingDataExons)

    if(curr->exon > 0 && checkMissingExons(missingExons, 0, curr->exon -1)){
        exon = curr->exon;
        curr = curr->next;
    }

    while(curr != NULL){
        st_logInfo("currExon: %d, start: %d, end: %d, exonStart: %d, exonEnd: %d, chr: %s, seq: %s\n", 
                    curr->exon, curr->start, curr->start + curr->length, curr->exonStart, curr->exonEnd, curr->chr, curr->seq);
        if( exon == exonCount -1 || (curr->exon < exon && strcmp(curr->seq, seq) != 0) ){
            if( strand == '+' && (curr->prev->status == 1 || exonMissingData) ){
                return true;
            }else if(curr->prev->isRightStub && curr->exonEnd > 0){
                /*if(!intListContains(mdExons, exon)){
                    intListAppend(mdExons, exon);
                }*/
                return true;
            }
        }

        if( strcmp(curr->chr, chr) == 0 ){//same chromosome
            st_logInfo("\tcurrChr: %s == chr : %s\n", curr->chr, chr);
            if(curr->strand != strand){//switch strand, reset
                exon = -1;
                strand = curr->strand;
                exonMissingData = false;
                continue;
            }else if(curr->exon < exon){//startOver at curr
                st_logInfo("\t\tcurrExon (%d) < exon (%d) \n", curr->exon, exon);
                exon = -1;
                exonMissingData = false;
                continue;
            }else if (curr->exon == exon){
                st_logInfo("\t\tcurrExon (%d) = exon: (%d)\n", curr->exon, exon);
                if(curr->exon == 0 && curr->status == 1){//previous block is also exon 0, but since it's the first exon, we can start over the gene here
                    st_logInfo("\t\t\t exon == 0 and status = 1, set missingData to false\n");
                    exonMissingData = false;
                }else{
                    if(curr->prev->exonStart + curr->prev->exonLength <= curr->exonStart){//keep moving on if curr frag on the right of prevFrag
                        st_logInfo("\t\t\t prevEstart (%d) + prevElen (%d) <= currEstart (%d) \n", curr->prev->exonStart, curr->prev->exonLength, curr->exonStart);

                        if(curr->prev->isRightStub && curr->isLeftStub && 
                           curr->strand == '+' && strcmp(curr->seq, seq) != 0){//detect missingData, set it to true
                            st_logInfo("\t\t\t\t prevIsRightStub, currIsLeftStub, and currStrand +, currSeq != prevSeq. Set missingData to True\n");
                            exonMissingData = true;
                            /*if(!intListContains(mdExons, exon)){
                                intListAppend(mdExons, exon);
                            }*/
                        }else{
                            st_logInfo("\t\t\t\t hasn't detected missingData yet\n");
                        }
                    }else{//reset if curr frag on the left of or overlap with prevFrag
                        exonMissingData = false;
                        exon = -1;
                        if(curr->exon == 0){continue;}
                    }
                }
            }else if (curr->exon == exon + 1){
                st_logInfo("\t\tcurrExon %d = exon%d + 1\n", curr->exon, exon);
                //take care of previous exon
                if(curr->exon != 0 && curr->prev->isRightStub && curr->prev->exonEnd > 0){
                    /*if(!intListContains(mdExons, exon)){
                        intListAppend(mdExons, exon);
                    }*/
                    exonMissingData = true;
                }
                
                if(curr->exon != 0 && curr->prev->status == 0 && !exonMissingData){//previous exon not perfect, and does not contain missing data, reset
                    exon = -1;
                }else{//currExon is 0, or prevExon was perfect, or prevExon not perfect but contain missingDtaa, keep moving
                    exon++;
                    if(curr->isLeftStub && curr->exonStart > 0){//current exon contain missing data
                        /*if(!intListContains(mdExons, exon)){
                            intListAppend(mdExons, exon);
                        }*/
                        exonMissingData = true;
                    }else{
                        exonMissingData = false;
                    }
                }
            }else{//curr->exon > exon + 1
                st_logInfo("\t\tcurrExon (%d ) > exon: (%d) + 1\n", curr->exon, exon);
                if( strcmp(curr->seq, seq) != 0 && checkMissingExons(missingExons, exon+1, curr->exon -1) ){
                    st_logInfo("\t\t\tmissing data detected..., set exon to %d\n", curr->exon);
                    exon = curr->exon;
                    if(curr->isLeftStub && curr->exonStart > 0){//current exon contain missing data
                        /*if(!intListContains(mdExons, exon)){
                            intListAppend(mdExons, exon);
                        }*/
                        exonMissingData = true;
                    }else{
                        exonMissingData = false;
                    }
                }else{//missed exons
                    st_logInfo("\t\t\tRESET\n");
                    exon = -1;
                    exonMissingData = false;
                }
            }

            if(curr->next == NULL){//the last mapped exon in the list
                if(exon != -1 && curr->isRightStub && curr->exonEnd > 0){
                    /*if(!intListContains(mdExons, exon)){
                        intListAppend(mdExons, exon);
                    }*/
                    exonMissingData = true;
                }
                if(curr->status == 1 || exonMissingData){
                    if(exon == exonCount -1){
                        hasHomolog = true;
                    }else if(exon < exonCount -1){
                        hasHomolog = checkMissingExons(missingExons, exon +1, exonCount -1);
                    }
                }else{
                    hasHomolog = false;
                }
            }

            seq = curr->seq;
            curr = curr->next;
        }else{//diff chromosome
            exon = -1;
            chr = curr->chr;
            exonMissingData = false;
        }
    }
    return(hasHomolog && strand == '+');
}

//bool checkPartialHomologyBackward(struct MappedExon *mexons, struct IntList *missingExons, struct IntList *mdExons, int exonCount){
bool checkPartialHomologyBackward(struct MappedExon *mexons, struct IntList *missingExons, int exonCount){
    /*
     * A 'partially-complete' homolog is similar to a 'complete' homolog, but with a few exceptions:
     *  a/ It is a traversal of the exons in the same order and orientation with the reference exons. 
     *     However, the exons in the traversal are allowed to have missing data.
     *  b/ Absent exons due to missing data are allowed
     *
     */
    st_logInfo("\nchecking partial homology backward:\n");
    int exon = exonCount;
    struct MappedExon *curr = mexons;

    char *chr = curr->chr;
    char strand = curr->strand;
    char *seq = curr->seq;
    bool hasHomolog = NULL;
    bool exonMissingData = false;
    
    if(curr->exon < exonCount -1 && checkMissingExons(missingExons, curr->exon + 1, exonCount -1)){
        exon = curr->exon;
        curr = curr->next;
    }
    while(curr != NULL){
        if( exon == 0 || (curr->exon > exon && strcmp(curr->seq, seq) != 0) ){
            if( strand == '-' && (curr->prev->status == 1 || exonMissingData) ){
                return true;
            }else if(curr->prev->isLeftStub && curr->exonStart > 0){
                /*if(!intListContains(mdExons, exon)){
                    intListAppend(mdExons, exon);
                }*/
                return true;
            }
        }

        if( strcmp(curr->chr, chr) == 0 ){//same chromosome
            if(curr->strand != strand){//switch strand, reset
                exon = exonCount;
                strand = curr->strand;
                exonMissingData = false;
                continue;
            }else if(curr->exon > exon){//startOver at curr
                exon = exonCount;
                exonMissingData = false;
                continue;
            }else if (curr->exon == exon){
                if(curr->exon == exonCount -1 && curr->status == 1){//previous block is also exon 0, but since it's the first exon, we can start over the gene here
                    exonMissingData = false;
                }else{
                    if(curr->exonStart + curr->exonLength <= curr->prev->exonStart){//keep moving on if curr frag on the left of prevFrag
                        if(curr->prev->isLeftStub && curr->isRightStub && curr->strand == '-' && strcmp(curr->seq, seq) != 0){//detect missingData, set it to true
                            exonMissingData = true;
                            /*if(!intListContains(mdExons, exon)){
                                intListAppend(mdExons, exon);
                            }*/
                        }
                    }else{//reset if curr frag on the right of or overlap with prevFrag
                        exonMissingData = false;
                        exon = exonCount;
                        if(curr->exon == exonCount -1){continue;}
                    }
                }
            }else if (curr->exon == exon - 1){
                //take care of previous exon
                if(curr->exon != exonCount -1 && curr->prev->isLeftStub && curr->prev->exonStart > 0){
                    /*if(!intListContains(mdExons, exon)){
                        intListAppend(mdExons, exon);
                    }*/
                    exonMissingData = true;
                }
                
                if(curr->exon != exonCount -1 && curr->prev->status == 0 && !exonMissingData){//previous exon not perfect, and does not contain missing data, reset
                    exon = exonCount;
                }else{//currExon is exontCount -1, or prevExon was perfect, or prevExon not perfect but contain missingDtaa, keep moving
                    exon --;
                    if(curr->isRightStub && curr->exonEnd > 0){//current exon contain missing data
                        /*if(!intListContains(mdExons, exon)){
                            intListAppend(mdExons, exon);
                        }*/
                        exonMissingData = true;
                    }else{
                        exonMissingData = false;
                    }
                }
            }else{//curr->exon < exon - 1
                if( strcmp(curr->seq, seq) != 0 && checkMissingExons(missingExons, curr->exon + 1, exon -1) ){
                    exon = curr->exon;
                    if(curr->isRightStub && curr->exonEnd > 0){//current exon contain missing data
                        /*if(!intListContains(mdExons, exon)){
                            intListAppend(mdExons, exon);
                        }*/
                        exonMissingData = true;
                    }else{
                        exonMissingData = false;
                    }
                }else{//missed exons
                    st_logInfo("\t\t\tRESET\n");
                    exon = exonCount;
                    exonMissingData = false;
                }
            }

            if(curr->next == NULL){//the last mapped exon in the list
                if(exon != exonCount && curr->isLeftStub && curr->exonStart > 0){
                    /*if(!intListContains(mdExons, exon)){
                        intListAppend(mdExons, exon);
                    }*/
                    exonMissingData = true;
                }
                if(curr->status == 1 || exonMissingData){
                    if(exon == 0){
                        hasHomolog = true;
                    }else if(exon > 0){
                        hasHomolog = checkMissingExons(missingExons, 0, exon - 1);
                    }
                }else{
                    hasHomolog = false;
                }
            }
            seq = curr->seq;
            curr = curr->next;
        }else{//diff chromosome
            exon = exonCount;
            chr = curr->chr;
            exonMissingData = false;
        }
    }
    return(hasHomolog && strand == '-');
}

void printIntList(FILE *fileHandle, struct IntList *list){
    if(list == NULL){return;}
    for(int i=0; i< list->length; i++){
        fprintf(fileHandle, "%d,", list->list[i]);
    }
    //fprintf(fileHandle, "\n");
    return;
}

bool hasIntronDups(struct IntList *list){
    assert(list != NULL);
    for(int i=0; i< list->length; i++){
        if(list->list[i] == 1){
            return true;
        }
    }
    return false;
}

bool hasExonDups(struct List *list){
    struct Exon *exon;
    for(int i=0; i< list->length; i++){
        exon = list->list[i];
        if(exon->hasDup == 1){
            return true;
        }
    }
    return false;
}

void printExonDups(FILE *fileHandle, struct List *list){
    struct Exon *exon;
    for(int i=0; i< list->length; i++){
        exon = list->list[i];
        fprintf(fileHandle, "%d,", exon->hasDup);
    }
    return;
}

void convertExonIdOrder(struct IntList *mdexons, int max){
    for(int i=0; i< mdexons->length; i++){
        mdexons->list[i] = max - mdexons->list[i];
    }
    return;
}

bool mapGene(FILE *fileHandle, struct Gene *refgene, struct List *spcList){
    //fprintf(stderr, "GENE: %s\n", refgene->name);
    fprintf(stdout, "\t<gene name=\"%s\" exonCount=\"%d\" strand=\"%c\">\n", refgene->name, refgene->exons->length, refgene->strand);
    
    //Intronic duplication
    fprintf(stdout, "\t\t<intronDups>");
    printIntList(stdout, refgene->introns);
    fprintf(stdout, "</intronDups>\n");
    //Exonic duplication
    fprintf(stdout, "\t\t<exonDups>");
    printExonDups(stdout, refgene->exons);
    fprintf(stdout, "</exonDups>\n");


    bool hasHomolog = true;
    bool hasPartialHomolog = false;
    int visitedSpecies = 0;
    char *geneMissingData = "no";

    for (int i=0; i< spcList->length; i++){//each species
        char *species = spcList->list[i];
        if(strstr(species, refgene->species) != NULL){visitedSpecies ++; continue;}//ignore the refspecies
        struct MappedExon *sortedMappedExons = NULL;
        struct IntList *missingExons = constructEmptyIntList(0); //list of exon that has 0 segments of current species
        struct IntList *mdExons; //list of exon that contains missing data 
        st_logInfo("\nSPECIES: %s\n", species);
        st_logInfo("Getting exon coverage...\n");
        for(int j=0; j< refgene->exons->length; j++){//each exon
            struct Exon *exon = refgene->exons->list[j];
            struct List *meList = exon_getCoverage(exon, species);
            if(meList->length == 0){
                intListAppend(missingExons, j);
            }

            //=========== DEBUG ===============
            st_logInfo("\n\tExon %d:\n", j);
            for (int k=0; k< meList->length; k++){
                struct MappedExon *me = meList->list[k];
                char *leftstub = me->isLeftStub ? "leftStub" : "";
                char *rightstub = me->isRightStub ? "rightStub" : "";
                st_logInfo("\t\tStart: %d, end: %d, exonStart: %d, exonEnd: %d, blockIndex: %d, status: %d, %s, %s\n", me->start, 
                            me->start + me->length, me->exonStart,me->exonEnd,  me->blockIndex, me->status, leftstub, rightstub);
            }
            //===========END DEBUG=================
            
            if(sortedMappedExons == NULL && meList->length > 0){
                sortedMappedExons = meList->list[0];
            }
            for(int k=0; k< meList->length; k++){
                sortedMappedExons = insertMappedExon(sortedMappedExons, meList->list[k]);
            }
        }

        //=========== DEBUG ===============
        st_logInfo("\nMissing exons: \t");
        for(int j=0; j< missingExons->length; j++){
            st_logInfo("%d, ", j);
        }
        st_logInfo("\n");
        st_logInfo("\nSorted MappedExons:\n");
        struct MappedExon *e = sortedMappedExons;
        while(e != NULL){
            st_logInfo("\tExon: %d, Start: %d, end: %d, exonStart: %d, exonEnd: %d, blockIndex: %d, chr: %s, seqs: %s, status: %d\n", 
                        e->exon, e->start, e->start + e->length, e->exonStart, e->exonEnd, e->blockIndex, e->chr, e->seq, e->status);
            e = e->next;
        }
        //===========END DEBUG=================
        if(sortedMappedExons == NULL){//doesn't found the gene in current species
            fprintf(stdout, "\t\t<species name=\"%s\" status=\"NA\" missingdataexons=\"\">\n", species);
            fprintf(stdout, "\t\t</species>\n");
            continue;
        }
        visitedSpecies ++;

        sortedMappedExons = mergeIncompleteMappedExons(sortedMappedExons);
        sortedMappedExons = mergeIncompleteMappedExons2(sortedMappedExons);
        
        if (gene_checkMissingData(sortedMappedExons)){
            geneMissingData = "yes";
        }

        mdExons = gene_findMissingDataExons(sortedMappedExons);
        
        correctForAltSplice(sortedMappedExons);

        if (checkCompleteHomologyForward(sortedMappedExons, refgene->exons->length) ||
            checkCompleteHomologyBackward(sortedMappedExons, refgene->exons->length)){
            fprintf(stdout, "\t\t<species name=\"%s\" status=\"C\" ", species);
        }else if (checkPartialHomologyForward(sortedMappedExons, missingExons, refgene->exons->length) ||
                  checkPartialHomologyBackward(sortedMappedExons, missingExons, refgene->exons->length)){
            fprintf(stdout, "\t\t<species name=\"%s\" status=\"MD\" ", species);
            hasPartialHomolog = true;
            geneMissingData = "yes";//HACK
        }else{
            fprintf(stdout, "\t\t<species name=\"%s\" status=\"N\" ", species);
            hasHomolog = false;
        }
        fprintf(stdout, "missingdataexons=\"");
        if(refgene->strand == '-'){
            convertExonIdOrder(mdExons, refgene->exonCount - 1);
        }
        printIntList(stdout, mdExons);
        fprintf(stdout, "\">\n");
        
        //=============DEBUG============
        st_logInfo("\nMerge Incomplete MappedExons:\n");
        e = sortedMappedExons;
        fprintf(stdout, "\t\t\t<exonString>\n");
        while(e != NULL){
            st_logInfo("\tExon: %d, Start: %d, end: %d, exonStart: %d, exonEnd: %d, blockIndex: %d, chr: %s, seq: %s, status: %d\n", 
                        e->exon, e->start, e->start + e->length, e->exonStart, e->exonEnd, e->blockIndex, e->chr, e->seq, e->status);
            int exonid = refgene->strand == '+' ? e->exon : refgene->exonCount - 1 - e->exon;
            fprintf(stdout, "\t\t\t\t<exon id=\"%d\" status=\"%d\" strand=\"%c\" insBases=\"%d\" delBases=\"%d\"></exon>\n", exonid, e->status, e->strand, e->insBases, e->delBases);
            //fprintf(stdout, "\t\t\t\t<exon id=\"%d\" status=\"%d\" strand=\"%c\"></exon>\n", e->exon, e->status, e->strand);
            e = e->next;
        }
        fprintf(stdout, "\t\t\t</exonString>\n");
        fprintf(stdout, "\t\t</species>\n");
        //fprintf(fileHandle, "\n");
        //============END DEBUG ============
    }
   
    char *intronDups = hasIntronDups(refgene->introns) ? "yes" : "no";//if has dup within intronic region, intronDups = 1, else = 0
    char *exonDups = hasExonDups(refgene->exons) ? "yes" : "no";
    char *checkHomolog = hasHomolog ? "C" : "N";
    if(visitedSpecies == 1){//only the reference species
        checkHomolog = "NA";
        fprintf(fileHandle, "%s\t%s\n", refgene->name, checkHomolog);
    }else{
        //if(hasHomolog && hasPartialHomolog){//some species has partial homologs
        if( hasHomolog && (hasPartialHomolog || visitedSpecies < spcList->length) ){//some species has partial homologs
            checkHomolog = "MD";
        }
        fprintf(fileHandle, "%s\t%s\t%s\t%s\n", refgene->name, checkHomolog, exonDups, intronDups);
    }

    fprintf(stdout, "\t\t<missingData>%s</missingData>\n", geneMissingData);
    //fprintf(stdout, "\t\t<intronDups>%s</intronDups>\n", intronDups);
    //fprintf(stdout, "\t\t<exonDups>%s</exonDups>\n", exonDups);
    fprintf(stdout, "\t\t<status>%s</status>\n", checkHomolog);
    fprintf(stdout, "\t</gene>\n");
    
    return hasHomolog;
}

//================ MAP GENE TO REFERENCE SEQUENCE ================
//overlap_walkUp and overlap_walkDown traverse cactus thread and find blocks that are overlapped
//with the 'start' and 'end' coordinates. If inputed 'exon' is not NULL, record the list of 5'caps
//of overlapped blocks. 
//If hasDup is specified (not NULL), overlap_walkUp/Down only checks for and return
//as soon as find any duplication within [start, end)
Cap *overlap_walkDown(Cap *cap, int32_t start, int32_t end, struct Exon *exon, int *hasDup);

Cap *overlap_walkUp(Cap *cap, int32_t start, int32_t end, struct Exon *exon, int *hasDup){
    assert((exon != NULL && hasDup == NULL) || (exon == NULL && hasDup != NULL));
    int32_t coor = mapCapCoor(cap);
    //st_logInfo("overlap_walkUp: cap %d, start %d, end %d\n", coor, start, end);
    if(end <= coor){//cases of missing data
        (*hasDup) = 0;
        return cap;
    }
    //coor < end:
    Segment *segment = cap_getSegment(cap);
    if(segment != NULL){
        Cap *otherCap = cap_getOtherSegmentCap(cap);
        int checkDup = block_checkDuplication(segment_getBlock(segment));
        if(hasDup != NULL){
            *hasDup = checkDup;
            if((*hasDup) == 1){ 
                return cap; 
            }//the whole region has some duplication, return
        }else{
            if(checkDup == 1){exon->hasDup = 1;}
            int32_t exonStart = mapCapCoor(cap) - start; 
            int32_t exonEnd = end - (mapCapCoor(cap) + segment_getLength(segment)); 
            struct MappedEnd *mappedEnd = constructMappedEnd(cap, cap_getEnd(cap), exonStart, exonEnd);
            listAppend(exon->ends, mappedEnd);
        }
        if(mapCapCoor(otherCap) < end -1){//NOT YET end of exon(/intron/region), keep walking
            cap = overlap_walkDown(otherCap, start, end, exon, hasDup);
        }//else: end of exon
    }else{
        Group *parentGroup = flower_getParentGroup(end_getFlower(cap_getEnd(cap)));
        if (parentGroup != NULL) {
            Cap *upperCap = flower_getCap(group_getFlower(parentGroup), cap_getName(cap));
            assert(upperCap != NULL);
            //Make sure always following the same strand
            if(cap_getStrand(cap) != cap_getStrand(upperCap)){
                upperCap = cap_getReverse(upperCap);
            }
            cap = overlap_walkUp(upperCap, start, end, exon, hasDup);
        }
    }
    return cap;
}

Cap *overlap_walkDown(Cap *cap, int32_t start, int32_t end, struct Exon *exon, int *hasDup){
    //int32_t coor = cap_getCoordinate(cap);
    //int32_t coor = mapCapCoor(cap);
    //st_logInfo("overlap_walkDown: cap %d, start %d, end %d\n", coor, start, end);
    //assert(start <= coor && coor < end -1);

    Cap *adjCap = cap_getAdjacency(cap);
    Group *group = end_getGroup(cap_getEnd(cap));
    if(group_isLeaf(group)){
        cap = overlap_walkUp(adjCap, start, end, exon, hasDup);
    }else{
        Cap *lowerCap = flower_getCap(group_getNestedFlower(group), cap_getName(cap));
        if(cap_getStrand(cap) != cap_getStrand(lowerCap)){
            lowerCap = cap_getReverse(lowerCap);
        }
        cap = overlap_walkDown(lowerCap, start, end, exon, hasDup);
    }
    return cap;
}

void mapGeneToRef(struct Gene *mg, Cap *cap, struct bed *gene){
    st_logInfo("mapGeneToRef: gene %s, startCap %d\n", gene->name, mapCapCoor(cap));
    int32_t eStart, eEnd;
    int32_t iStart, iEnd;
    struct Exon *exon;

    for(int i=0; i< gene->blockCount; i++){//each exon
        int exonid = i;
        eStart = gene->chromStart + gene->chromStarts->list[i];
        eEnd = eStart + gene->blockSizes->list[i];
        //fprintf(stderr, "exon %d, eStart: %d, eEnd: %d, cap: %d\n", i, eStart, eEnd, mapCapCoor(cap));
        assert(mg->exons->length >= i);
        if(mg->exons->length == i){
            exon = constructExon(exonid);
            exon->start = eStart;
            exon->end = eEnd;
        }else{//mg->exons->legth > i
            exon = mg->exons->list[i];
        }
        cap = overlap_walkUp(cap, eStart, eEnd, exon, NULL);
        listAppend(mg->exons, exon);

        if(i < gene->blockCount -1){//not the last exon
            //check if intron has duplications...
            iStart = eEnd;
            iEnd = gene->chromStart + gene->chromStarts->list[i+1];
            //fprintf(stderr, "Done exon, start intron %d, iStart %d, iEnd %d, cap: %d\n", i, iStart, iEnd, mapCapCoor(cap));
            int hasDup = 0;
            //fprintf(stderr, "\n\nINTRON %d\n", i);
            overlap_walkUp(cap, iStart, iEnd, NULL, &hasDup);
            //Cap *cap2 = overlap_walkUp(cap, iStart, iEnd, NULL, &hasDup);
            //fprintf(stderr, "Done intron %d, cap2: %d\n", i, mapCapCoor(cap2));
            intListAppend(mg->introns, hasDup);
            
            //move to next exon:
            cap = walkUp(cap, iEnd);
            //fprintf(stderr, "walk to nextExon, cap: %d\n", mapCapCoor(cap));
        }
    }
    return;
}

//==================================
struct List *flower_getSpecies(Flower *flower){
    /*
     * Return the list of all species (e.g: hg18, panTro2, ponAbe2...) 
     * that are in the inputed flower
     */
    struct List *list = constructEmptyList(0, free);
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *seq;
    while((seq = flower_getNextSequence(it)) != NULL){
        char *seqHeader = formatSequenceHeader(seq);
        struct List *headerList = splitString(seqHeader, ".");
        assert(headerList->length > 0);
        if (!visitedString(list, headerList->list[0])){
            listAppend(list, headerList->list[0]);
        }
    }
    flower_destructSequenceIterator(it);
    return list;
}

//============================ MAP ALL GENES =========================
void mapGenes(Flower *flower, FILE *fileHandle, struct bed *gene, char *refspecies){
    st_logInfo("Flower %s\n", cactusMisc_nameToString(flower_getName(flower)));
    struct List *spcList = flower_getSpecies(flower);
    st_logInfo("Found %d spcies\n", spcList->length);
    char *geneChrom = NULL;
    struct CapList *refStartCapList = NULL;
    struct CapList *capNode;

    fprintf(stdout, "<geneHomologs>\n");
    while(gene != NULL){//The inputed genes were sorted by their chroms, then startCoors
        //Get the start of the target sequence: 
        st_logInfo("\nGene %s:\n", gene->name);
        
        //The sequence name of the refspecies must matches the pattern "spc.chr#", for ex: hg19.chr19
        if(geneChrom == NULL || strcmp(geneChrom, gene->chrom) != 0){//new chromosome 
            geneChrom = gene->chrom;
            st_logInfo("New chromosome: %s\n", geneChrom);
            struct List *strlist = constructEmptyList(0, free);
            listAppend(strlist, refspecies);
            listAppend(strlist, geneChrom);
            char *refname = str_joinList(strlist, ".");
            //destructList(strlist);
            st_logInfo("Finding thread starts of %s\n", refname);
            refStartCapList = flower_getThreadStart(flower, refname);
        }
        
        //DEBUG
        /*st_logInfo("Printing current refStartCapList:\n");
        capNode = refStartCapList;
        while(capNode != NULL){
            char *currStrand = cap_getStrand(capNode->cap) ? "+": "-";
            char *currSide = cap_getSide(capNode->cap) ? "5" : "3";
            st_logInfo("cap: %d, %s, %s\n", mapCapCoor(capNode->cap), currStrand, currSide);
            capNode = capNode->next;
        }*/
        //END DEBUG

        char refstrand = strcmp(gene->strand, "+") == 0 ? '+' : '-';
        struct Gene *refGene = constructGene(gene->name, refstrand, gene->blockCount);//mapped gene
        refGene->species = refspecies;
        capNode = refStartCapList;
        Cap *startCap;

        while(capNode != NULL){
            startCap = capNode->cap;
            st_logInfo("Current startCap is at: %d\n", mapCapCoor(startCap));
            if( (gene->chromEnd <= mapCapCoor(startCap)) ||
                (capNode->endSeqCoor <= gene->chromStart) ){//if not overlap, move to next sequence
                st_logInfo("current sequence does not overlap with gene, move to next Start\n");
                capNode = capNode->next;
                continue;
            }
            st_logInfo("Move to the 5'-most cap of the gene:\n");
            if(cap_getSide(startCap)){ 
                startCap = walkUp(startCap, gene->chromStart);
            }else{
                startCap = walkDown(startCap, gene->chromStart);
            }
            capNode->cap = startCap;
            
            /*char *strand = cap_getStrand(startCap) ? "+" : "-";
            char *side = cap_getSide(startCap)? "5" : "3";
            st_logInfo("startCap now is at: %d, %s, %s\n", mapCapCoor(startCap), strand, side);*/
            
            //Traverse the reference seq 
            mapGeneToRef(refGene, startCap, gene);

            //DEBUG
            /*st_logInfo("\nmapGeneToRef: %s, %s\n", refGene->name, refGene->species);
            for(int e = 0; e< refGene->exons->length; e++){
                st_logInfo("Exon %d: ", e);
                struct Exon *exon = refGene->exons->list[e];
                for(int f=0; f < exon->ends->length; f++){
                    struct MappedEnd *myend = exon->ends->list[f];
                    st_logInfo("\tEnd %s, index%d\n", cactusMisc_nameToString(end_getName(myend->end)), f);
                }
            }*/
            //END DEBUG
            
            capNode = capNode->next;
        }
        //mapGene(fileHandle, refGene, spcList);
        if( mapGene(fileHandle, refGene, spcList) ){
            //fprintf(stdout, "%s\tYES\n", gene->name);
        }else{
            //fprintf(stdout, "%s\tNO\n", gene->name);
        }
        gene = gene->next;
    }
    fprintf(stdout, "</geneHomologs>\n");
    return;
}

void usage() {
   fprintf(stderr, "cactus_genemapHomolog, version 0.2\n");
   fprintf(stderr, "-a --st_logLevel : Set the st_log level\n");
   fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
   fprintf(stderr, "-o --outputFile : output file, in xml format.\n");
   fprintf(stderr, "-s --species : Name of reference species, e.g 'hg18.chr11.134452384.116118685.50084.1'\n");
   fprintf(stderr, "-g --geneFile : File that contains Beds of transcripts/genes/CDS.\n");
   fprintf(stderr, "-h --help : Print this help screen\n");
}

//============================== MAIN =========================================
int main(int argc, char *argv[]) {
   Flower *flower;

   /*
    * Arguments/options
    */
   char * st_logLevelString = NULL;
   char * cactusDiskDatabaseString = NULL;
   char * flowerName = "0";
   char * outputFile = NULL;
   char * species = NULL;
   char * geneFile = NULL;

   ///////////////////////////////////////////////////////////////////////////
   // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
   ///////////////////////////////////////////////////////////////////////////

   while(1) {
      static struct option long_options[] = {
         { "genePslFile", required_argument, 0, 'g' },
         { "species", required_argument, 0, 's' },
         { "st_logLevel", required_argument, 0, 'a' },
         { "cactusDisk", required_argument, 0, 'c' },
         { "outputFile", required_argument, 0, 'o' },
         { "help", no_argument, 0, 'h' },
         { 0, 0, 0, 0 }
      };

      int option_index = 0;

      int key = getopt_long(argc, argv, "s:g:o:a:c:h", long_options, &option_index);

      if(key == -1) {
         break;
      }

      switch(key) {
         case 'a':
            st_logLevelString = stString_copy(optarg);
            break;
         case 'c':
            cactusDiskDatabaseString = stString_copy(optarg);
            break;
         case 'o':
            outputFile = stString_copy(optarg);
            break;
         case 's':
            species = stString_copy(optarg);
            break;
         case 'g':
            geneFile = stString_copy(optarg);
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
   assert(outputFile != NULL);
   assert(species != NULL);
   assert(geneFile != NULL);

   //////////////////////////////////////////////
   //Set up st_logging
   //////////////////////////////////////////////

   if(st_logLevelString != NULL && strcmp(st_logLevelString, "INFO") == 0) {
      st_setLogLevel(ST_LOGGING_INFO);
   }
   if(st_logLevelString != NULL && strcmp(st_logLevelString, "DEBUG") == 0) {
      st_setLogLevel(ST_LOGGING_DEBUG);
   }

   //////////////////////////////////////////////
   //Log (some of) the inputs
   //////////////////////////////////////////////

   st_logInfo("Flower disk name : %s\n", cactusDiskDatabaseString);
   st_logInfo("Output file : %s\n", outputFile);
   st_logInfo("Species: %s\n", species);
   st_logInfo("GenePslFile: %s\n", geneFile);

   //////////////////////////////////////////////
   //Load the database
   //////////////////////////////////////////////

   stKVDatabaseConf *kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
   CactusDisk *cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
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
   struct bed *gene = bedLoadAll(geneFile);
   mapGenes(flower, fileHandle, gene, species);
   fclose(fileHandle);
   st_logInfo("Map genes in %i seconds/\n", time(NULL) - startTime);

   ///////////////////////////////////////////////////////////////////////////
   // Clean up.
   ///////////////////////////////////////////////////////////////////////////

   cactusDisk_destruct(cactusDisk);

   return 0;
}
