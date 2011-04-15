#include <assert.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>
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
 * nknguyen@soe.ucsc.edu
 * Dec 1 2010, edit output to be in the 'annotated maf' format (with i and e
 * rows, for details see:http://genome.ucsc.edu/FAQ/FAQformat.html#format5
 * Note: haven't dealt with missing data
 */

//======= Global structures ======
struct MafSegment {
    Segment *segment;
    char *name;
    int32_t srcSize;
    int32_t insertSize; //size of the insert between previous block and current block if there is any
    int32_t gapSize; //size of the deletion between previous block and current block if there is any
    int32_t gapStart;
    char strand;
    bool empty; //(white space)
    bool unAligned; //double line
    bool missingData; //missing data
    struct MafSegment *prev;//previous segment
    struct MafSegment *next;//next segment
};

//====== Initialization functions ========
struct MafSegment *constructMafSegment(Segment *segment){
    struct MafSegment *mafSegment;
    mafSegment = st_malloc(sizeof(struct MafSegment));
    mafSegment->segment = segment;
    mafSegment->name = NULL;
    mafSegment->srcSize = 0;
    mafSegment->insertSize = 0;
    mafSegment->gapSize = 0;
    mafSegment->empty = false;
    mafSegment->unAligned = false;
    mafSegment->missingData = false;
    mafSegment->prev = NULL;
    mafSegment->next = NULL;
    return mafSegment;
}

void destructMafSegment(struct MafSegment *ms){
    free(ms);
}

//=================
int32_t getSrcSize(Segment *segment){
    int32_t srcSize = 0;
    if(segment != NULL){
	Sequence *sequence = segment_getSequence(segment);
	assert(sequence != NULL);
	srcSize = sequence_getLength(sequence);
    }
    return srcSize;
}

char *getSegmentName(Segment *segment){
    assert(segment != NULL);
    Cap *cap = segment_get5Cap(segment);
    return cap_getSequenceName(cap);
}

//========================
int32_t getSegmentStart(Segment *segment){
    int32_t start = -1;
    Sequence *sequence = segment_getSequence(segment);
    if (sequence != NULL) {
        if (segment_getStrand(segment)) {
            start = segment_getStart(segment) - sequence_getStart(sequence);
        } else { //start with respect to the start of the reverse complement sequence
            start = (sequence_getStart(sequence) + sequence_getLength(sequence)
                    - 1) - segment_getStart(segment);
        }
    }
    return start;
}

bool checkContinuity(struct MafSegment *ms1, struct MafSegment *ms2, int32_t gapSize){
    /*
     * Return true if ms1 and ms2 have same strand and:
     *[ms1Start-ms1End](abuts)[ms2Start-ms2End]  or [ms2S-ms2E][ms1S-ms2E]
     */

    //if(ms1->segment == NULL || ms2->segment == NULL){return true;}
    assert(ms1->segment != NULL);
    assert(ms2->segment != NULL);
    
    if(strcmp(ms1->name, ms2->name) != 0){//030111
        return false;
    }
    char strand1 = segment_getStrand(ms1->segment) ? '+' : '-';
    char strand2 = segment_getStrand(ms2->segment) ? '+' : '-';
    if(strand1 == strand2){
        int32_t start1 = getSegmentStart(ms1->segment);
        int32_t end1 = start1 + segment_getLength(ms1->segment);
        int32_t start2 = getSegmentStart(ms2->segment);
        //int32_t end2 = start2 + segment_getLength(ms2->segment);
        //if(start2 == end1 || start1 == end2){
        //if((strand1 == '+' && end1 + gapSize == start2) || //  |--ms1-->|---ms2---> 
        //   (strand1 == '-' && end2 + gapSize == start1)){  //  <--ms1---|<--ms2----|
        if(end1 + gapSize == start2) {
            return true;
        }
    }
    return false;
}

struct IntList *getRefMatchedColumns(struct List *refrow, Block *block){
    //Return indices (base 1) of reference cells (columns) that belong to the same block with 'block'
    //If reference segment belong to opposite-strand block, then return
    //index*(-1)
    assert(refrow->list != NULL);
    assert(block != NULL);
    struct IntList *indexList = constructEmptyIntList(0);
    int32_t i;
    for (i = 0; i< refrow->length; i++){
        struct MafSegment *refms = refrow->list[i];
        Segment *segment = refms->segment;
	if(segment == NULL){continue;}
        if(block == segment_getBlock(segment)){
	    intListAppend(indexList, i + 1);
	}else if(block_getReverse(block) == segment_getBlock(segment)){
	    intListAppend(indexList, (i+1)*(-1));
	}
    }
    return indexList;
}

struct List *getInitializedRow(int32_t length){
    //create a list of length 'length' of MafSegments (a row)
    struct List *row = constructEmptyList(0, free);
    for(int32_t i =0; i< length; i++){
        struct MafSegment *mafSegment = constructMafSegment( NULL );
        if (i >= 1){
            struct MafSegment *prevms = row->list[i-1];
            prevms->next = mafSegment;
            mafSegment->prev = prevms;
        }
        listAppend(row, mafSegment);   
    }
    return row;
}

bool checkInsert(struct List *row, int32_t c, struct IntList *prevCols, int32_t insertSize){
    /*
     *c = index of the 'matched' cell of current segment(base 1)
     *prevCols = indices of 'matched' cells of previous segment (that
     *aligned)(base1). If there is a column pc in prevCols so that
     *pc is immediately to left (when thread goes from left to right) 
     *or immediately to the right of c (when thread goes from right to left), 
     *returns true. Otherwise return false.
     */
    struct MafSegment *leftms, *rightms;
    int32_t i, pc, left, right;
    for(i = 0; i < prevCols->length; i++){
        pc = prevCols->list[i]; //previous column
        //if (c*pc > 0 && c == pc+1 ){ //Insertion. Note: if there is an insertion right before inversion, ignore
        if( ((pc > 0 && c>0) || (pc <0 && c <0)) && pc + 1 == c){
            if(c > 0){
                left = pc -1;
                right = c-1;
            }else{
                left = c*(-1) -1;
                right = pc*(-1) -1;
            }
            leftms = row->list[left];
            rightms = row->list[right];
            if(leftms == NULL || rightms == NULL){continue;}
            if(leftms->segment == NULL || rightms->segment == NULL){continue;}
            if(strcmp(leftms->name, rightms->name) != 0){continue;} //030111
            if (checkContinuity(leftms, rightms, insertSize)){
                rightms->insertSize = insertSize;
                return true;
            }
        }
    }
    return false;
}

//void fillInDoubleLine(struct List *row, int32_t c, struct IntList *prevCols, 
//                      int32_t gapStart, int32_t gapSize){
void fillInDoubleLine(struct List *row, int32_t c, struct IntList *prevCols, 
                      int32_t gapSize){
    assert(row != NULL && row->length >= c);
    bool hasDoubleline;
    int32_t pc, i, j;
    int32_t left, right;
    struct MafSegment * ms, *leftms, *rightms;

    for(i = 0; i< prevCols->length; i++){
        pc = prevCols->list[i];
        //if( pc*c > 0 && pc + 1 < c){
        if( ((pc > 0 && c>0) || (pc <0 && c <0)) && pc + 1 < c){
            if (c > 0){
                left = pc -1;
                right = c -1;
            }else{
                left = c*(-1) -1;
                right = pc*(-1) -1;
            }
            leftms = row->list[left];
            rightms = row->list[right];
            if(leftms == NULL || rightms == NULL){continue;}
            if(leftms->segment == NULL || rightms->segment == NULL){continue;}
            if(strcmp(leftms->name, rightms->name) != 0){continue;} //030111
            if (!checkContinuity(leftms, rightms, gapSize)){continue;}

            //st_logInfo("checkingDoubleLine: pc: %d, c: %d, left: %d, right: %d\n", pc, c, left, right);
            hasDoubleline = true;
            for(j= left+1; j < right; j++){//all cells in between pc and c must be gaps
                ms = row->list[j];
                if(ms->segment != NULL || ms->gapSize > 0){
                    hasDoubleline = false;
                    break;
                }
            }
            if(hasDoubleline){//fill in gap-cells
                for(j= left+1; j < right; j++){//all cells in btw pc and c are gaps
                    ms = row->list[j];
                    ms->srcSize = leftms->srcSize;
                    ms->name = leftms->name;
                    ms->gapSize = gapSize;
                    ms->strand = segment_getStrand(leftms->segment) ? '+' : '-';
                    ms->gapStart = getSegmentStart(leftms->segment) +
                                   segment_getLength(leftms->segment);
                    //ms->gapStart = gapStart;
                    ms->unAligned = true;
                }
                break;
            }
        }
    }
    return;
}

bool hasInsert(struct IntList *insertSizes){
    bool check = false;
    for(int32_t i=0; i< insertSizes->length; i++){
        if(insertSizes->list[i] > 0){
	    check = true;
            break;
	}
    }
    return check;
}

void fillInDeletion(struct List *refrow, struct List *row, int32_t c, struct IntList *prevCols){
    /*
     * If there exists a column pc in prevCols so that pc + 1 < c and
     * [pc+1, c-1] are empty cells, then we mark those cells as a deletion
     * (btwn pc and c)
     */
    assert(row != NULL && row->length >= c);
    bool hasDeletion;
    int32_t pc, i, j;
    int32_t left, right;
    int32_t gapSize = 0;
    struct MafSegment *ms, *refms, *leftms, *rightms;

    for(i = 0; i< prevCols->length; i++){
        pc = prevCols->list[i];
        gapSize = 0;
        //if( pc*c > 0 && pc + 1 < c){
        if( ((pc > 0 && c>0) || (pc <0 && c <0)) && pc + 1 < c){
            if (c > 0){
                left = pc -1; //convert back to base 0
                right = c -1; //convert back to base 0
            }else{
                left = c*(-1) -1;
                right = pc*(-1) -1;
            }
            assert(left < right);
            leftms = row->list[left];
            rightms = row->list[right];
            assert(leftms != NULL && rightms != NULL);
            if(leftms->segment == NULL || rightms->segment == NULL){continue;}
            if(strcmp(leftms->name, rightms->name) != 0){continue;}//030111
            
            //st_logInfo("checkingDeletion: pc: %d, c: %d, pc*c: %d, left: %d, right: %d\n", pc, c, pc*c, left, right);
            hasDeletion = true;
            for(j= left+1; j < right; j++){//all cells in between pc and c must be gaps
                ms = row->list[j];
                assert(ms != NULL);
                if(ms->segment != NULL || ms->gapSize > 0){
                    hasDeletion = false;
                    break;
                }
                refms = refrow->list[j];
	        gapSize += segment_getLength(refms->segment);
            }
            if(hasDeletion){//fill in gap-cells
                for(j= left+1; j < right; j++){//all cells in btw pc and c are gaps
                    ms = row->list[j];
                    ms->srcSize = leftms->srcSize;
                    ms->name = leftms->name;
                    ms->gapSize = gapSize;
                    ms->strand = segment_getStrand(leftms->segment) ? '+' : '-';
                    //if(c > 0){//goes from left to right
                        ms->gapStart = getSegmentStart(leftms->segment) +
                                       segment_getLength(leftms->segment);
                    /*}else{//goes from right to left
                        ms->gapStart = getSegmentStart(rightms->segment) +
                                       segment_getLength(rightms->segment);
                    }*/
                }
                break;
            }
        }
    }

    return;
}

bool block_hasRef(Block *block, char *name){
    bool hasRef = false;
    Block_InstanceIterator *it = block_getInstanceIterator(block);
    Segment *segment;
    while((segment = block_getNext(it)) != NULL){
        Cap *cap = segment_get5Cap(segment);
        char *currname = cap_getSequenceName(cap);
        if(strstr(currname, name)){//input 'block' contains segment of the reference
            hasRef = true;
            break;
        }
    }
    block_destructInstanceIterator(it);
    return hasRef;
}

void addMafSegment(struct MafSegment *ms, Segment *segment){
    ms->segment = segment;
    ms->srcSize = getSrcSize(segment);
    ms->name = getSegmentName(segment);
    return;
}

int32_t putSegmentToCell(Cap *cap, struct List *rows, struct List *refrow, 
                         char *refname, int32_t prevUnaligned, Cap *prevCap){
    Segment *segment = cap_getSegment(cap);
    Block *block = segment_getBlock(segment);

    //check if currentSegment aligns to anywhere on the ref species at all
    bool isAligned = block_hasRef(block, refname);

    struct IntList *cols = getRefMatchedColumns(refrow, block);
    //st_logInfo("Number of matched columns: %d\n", cols->length);

    struct List *r;
    int32_t c, i, j;
    //int32_t insert;
    struct MafSegment *mafsegment;
    struct IntList *prevCols = NULL;
    Segment *prevSegment = NULL;
 
    //st_logInfo("putSegmentToCell: %d, segmentLength: %d, prevUnaligned: %d\n",
    //          cap_getCoordinate(cap), segment_getLength(segment), prevUnaligned);

    if (cols->length > 0){//has matched columns
        if(prevCap != NULL){
	    prevSegment = cap_getSegment(prevCap);
	    assert(prevSegment != NULL);
	    prevCols =  getRefMatchedColumns(refrow, segment_getBlock(prevSegment));
        }//else: prevCap == NULL: beginning of thread... ignored..

        for(i = 0; i < cols->length; i++){//each match
	    c = cols->list[i];

            Segment *segment2 = segment;
	    if(c <0){//inversion
	        segment2 = segment_getReverse(segment);
		c = c*(-1);
	    }
            c -= 1;//change back to base-0

	    bool needNewRow = true;
	    for(j=0; j< rows->length; j++){//check to see if can fill segment into existing rows
		r = rows->list[j];
		mafsegment = r->list[c];
		if(mafsegment->segment == NULL && mafsegment->gapSize ==0){
		    //if(!checkDeletion(r, c, refrow)){//not deletion
                    addMafSegment(mafsegment, segment2);
                    needNewRow = false;
		    //}
		    break;
		}
	    }
	    if(needNewRow){//haven't found a cell for segment yet
                r = getInitializedRow(refrow->length);
                st_logInfo("Adding row #%d, length %d\n", rows->length, r->length);

		mafsegment = r->list[c];
                addMafSegment(mafsegment, segment2);
		listAppend(rows, r);
	    }

            //Check for insertion:
            bool hasInsert = false;
            if (prevCols != NULL && prevUnaligned >0){
                hasInsert = checkInsert(r, cols->list[i], prevCols, prevUnaligned);
            }
            //insert = hasInsert ? prevUnaligned : 0;
            //mafsegment->insertSize = insert;

            if( prevCols != NULL ){
                //if prevUnaligned == 0 && prevCols->length == 0: previous
                //bases aligns to somewhere else on the ref spc, but not
                //current 'refrow'
                if( prevUnaligned == 0){//check for deletion
                    fillInDeletion(refrow, r, cols->list[i], prevCols);
                }else if(!hasInsert){//doubleLine
                    //int32_t gapStart = getSegmentStart(prevSegment) + segment_getLength(prevSegment);
                    //fillInDoubleLine(r, cols->list[i], prevCols, gapStart, prevUnaligned);
                    fillInDoubleLine(r, cols->list[i], prevCols, prevUnaligned);
                }
            }
        }

        /*if(prevCap != NULL){
            destructIntList(cols);
            destructIntList(prevCols);
        }*/
        prevUnaligned = 0;
    }else if (isAligned){//current segment aligns somewhere on ref spc, but not current refrow
        prevUnaligned = 0;
    }else{//current segment does not algin to anywhere on the refrence spc
        prevUnaligned += segment_getLength(segment);
    }
    //st_logInfo(", newUnaligned %d\n", prevUnaligned);
    return prevUnaligned;
}

void walkDown(Cap *cap, struct List *rows, struct List *refrow, char *refname, int32_t prevUnaligned, Cap *prevCap);

void walkUp(Cap *cap, struct List *rows, struct List *refrow, char *refname, int32_t prevUnaligned, Cap *prevCap) {
    assert(cap != NULL);
    //st_logInfo("walkUp: %d, %s\n", cap_getCoordinate(cap), cactusMisc_nameToString(cap_getName(cap)));

    Segment *segment = cap_getSegment(cap);
    if (segment != NULL) {
        prevUnaligned = putSegmentToCell(cap, rows, refrow, refname, prevUnaligned, prevCap);//add segment to cell
        //prevUnaligned resets when segment aligns to one or more reference cell(s)
        if(prevUnaligned == 0){//reset prevCap to current cap
            prevCap = cap;
        }
        walkDown(cap_getOtherSegmentCap(cap), rows, refrow, refname, prevUnaligned, prevCap);
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
            walkUp(upperCap, rows, refrow, refname, prevUnaligned, prevCap);
        }
    }
}

void walkDown(Cap *cap, struct List *rows, struct List *refrow, char *refname, int32_t prevUnaligned, Cap *prevCap) {
    assert(cap != NULL);
    //st_logInfo("walkDown: %d\n", cap_getCoordinate(cap));
    //assert(end_isAttached(end));
    Group *group = end_getGroup(cap_getEnd(cap));
    if (group_isLeaf(group)) { //Walk across
        cap = cap_getAdjacency(cap);
        //Now walk up
        walkUp(cap, rows, refrow, refname, prevUnaligned, prevCap);
    } else { //Walk down
        Cap *lowerCap = flower_getCap(group_getNestedFlower(group), cap_getName(cap));
        if(cap_getStrand(cap) != cap_getStrand(lowerCap)){
            lowerCap = cap_getReverse(lowerCap);
        }
        walkDown(lowerCap, rows, refrow, refname, prevUnaligned, prevCap);
    }
}

void fillInEmptyCells(struct List *threadRows, struct List *refRow){
    struct MafSegment *rms;
    struct MafSegment *ms;
    for(int32_t j = 0; j < threadRows->length; j++){
        struct List *row = threadRows->list[j];
        for(int32_t i = 0; i< refRow->length; i++){
	    assert(refRow->length == row->length);
	    ms = row->list[i];
	    if(ms->segment == NULL && ms->gapSize == 0){
                rms = refRow->list[i];
		ms->empty = true;
		ms->gapSize = segment_getLength(rms->segment);
	    }
        } 
    } 
    return;
}

struct List *getRows(Flower *flower, char *name, struct List *refRows, char *refname){
    /*
     *Get rows for species 'name'
     */
    assert(refRows->list != NULL);
    struct List *rowsList = constructEmptyList(0, free);//rowsList->(ref)rows->row1, row2, ...
    struct List *refRow;
    Cap *cap;
    int i, j;
    
    for(i = 0; i < refRows->length; i++){//for each row of the reference species
        refRow = refRows->list[i];
        st_logInfo("\tGetting rows for source species %s that map to refRow %d\n", name, i);   
        st_logInfo("refRow Length: %d\n", refRow->length);

        //Get the starts of all the threads of current species
        //struct List *startCaps = getLeftThreadStarts(name, refRow);
        struct List *startCaps = flower_getThreadStarts(flower, name);
        if(startCaps->length == 0){
            st_logInfo("Could not find any %s sequence that aligns with the reference\n",name);
        }
        struct List *rows = constructEmptyList(0, free);
        for(j = 0; j < startCaps->length; j++){//each thread in the current species
            cap = startCaps->list[j];
            st_logInfo("\nCap %d: %s, sequence %s, coor: %d\n", j, cactusMisc_nameToString(cap_getName(cap)), 
                                                            cap_getSequenceName(cap), cap_getCoordinate(cap));
            walkDown(cap, rows, refRow, refname, 0, NULL);
        }

        //free startCaps list, but not delete the Caps themselves
        /*if(startCaps->length > 0){
            free(startCaps->list);
        }
        free(startCaps);*/
        //destructList(startCaps);
        fillInEmptyCells(rows, refRow);
        listAppend(rowsList, rows);    
        st_logInfo("\tDone getting rows for %s, refRow %d. Number of rows: %d\n", name, i, rows->length);
    }
    return rowsList;
}

//====================
void refWalkDown(Cap *cap, struct List *row);

void refWalkUp(Cap *cap, struct List *row) {
    assert(cap != NULL);
    st_logInfo("refWalkUp, cap %d, seq: %s\n", cap_getCoordinate(cap), cap_getSequenceName(cap));
    Segment *segment = cap_getSegment(cap);
    if (segment != NULL) {
        struct MafSegment *mafSegment = constructMafSegment( segment );
        mafSegment->srcSize = getSrcSize(segment);
        mafSegment->name = getSegmentName(segment);
        if(row->length >= 1){
            struct MafSegment *prevMs = row->list[row->length -1];
            prevMs->next = mafSegment;
            mafSegment->prev = prevMs;
        }
        listAppend(row, mafSegment);
        st_logInfo("\totherSegmentCap, cap %d, seq: %s\n", cap_getCoordinate(cap_getOtherSegmentCap(cap)), cap_getSequenceName(cap_getOtherSegmentCap(cap)));
        refWalkDown(cap_getOtherSegmentCap(cap), row);
    } else {
        //assert(end_isAttached(end));
        Group *parentGroup = flower_getParentGroup(end_getFlower(cap_getEnd(cap)));
        if (parentGroup != NULL) {
            Cap *upperCap = flower_getCap(group_getFlower(parentGroup), cap_getName(cap));
            if(cap_getStrand(cap) != cap_getStrand(upperCap)){
                upperCap = cap_getReverse(upperCap);
            }
            st_logInfo("\tupperCap, cap %d, seq: %s\n", cap_getCoordinate(upperCap), cap_getSequenceName(upperCap));
            refWalkUp(upperCap, row);
        }
    }
}

void refWalkDown(Cap *cap, struct List *row) {
    assert(cap != NULL);
    st_logInfo("refWalkDown, cap %d, seq: %s\n", cap_getCoordinate(cap), cap_getSequenceName(cap));
    //assert(end_isAttached(end));
    Group *group = end_getGroup(cap_getEnd(cap));
    if (group_isLeaf(group)) { //Walk across
        cap = cap_getAdjacency(cap);
        st_logInfo("\tadjCap, cap %d, seq: %s\n", cap_getCoordinate(cap), cap_getSequenceName(cap));
        //Now walk up
        refWalkUp(cap, row);
    } else { //Walk down
        Cap *lowerCap = flower_getCap(group_getNestedFlower(group), cap_getName(cap));
        if(cap_getStrand(cap) != cap_getStrand(lowerCap)){
            lowerCap = cap_getReverse(lowerCap);
        }
        st_logInfo("\tlowerCap, cap %d, seq: %s\n", cap_getCoordinate(lowerCap), cap_getSequenceName(lowerCap));
        refWalkDown(lowerCap, row);
    }
}

struct List *getReferenceRows(Flower *flower, char *name){
    /*
     * Each refRow represents a thread of the reference species
     * in the inputed flower. A thread could be a chromosome or a contig...
     */
    struct List *refRows = constructEmptyList(0, free);//list of rows of species 'name'
    Cap *cap;
    struct List *startCaps = flower_getThreadStarts(flower, name);
    for(int i = 0; i < startCaps->length; i++){
        cap = startCaps->list[i];
        struct List *row = constructEmptyList(0, free);
        refWalkDown(cap, row);
        listAppend(refRows, row);    
    }
    //free the startCaps list, but not the (Caps) themselves
    /*free(startCaps->list);
    free(startCaps);*/
    return refRows;
}

//================= PRINT MAF FOR EACH BLOCK ==================
char getLeftInfo(struct MafSegment *ms, int *count){
    assert(ms != NULL);
    char status;
    
    if(ms->prev == NULL){//start new sequence (or blue bar)
        status = 'N';
    }else{//(ms->segment != NULL) 
        if( ms->insertSize > 0 ){//insertion
            status = 'I';
            *count = ms->insertSize;
        }else if(ms->prev->segment != NULL){//prev block is not a gap
            if( checkContinuity(ms->prev, ms, 0) ){//continuous
                status = 'C';
            }else{//blue bar
                status = 'N';
            }
        }else{//ms->prev is a gap
            if(ms->prev->unAligned){
                status = 'I';
                *count = ms->prev->gapSize;
            }else if(ms->prev->missingData){
                status = 'M';
                *count = ms->prev->gapSize;
            }else if(ms->prev->empty){
                status = 'N';
            }else{
                status = 'C';
            }
        }
    }
    return status;
}

char getRightInfo(struct MafSegment *rightms, int *count){
    char status;
    
    if(rightms == NULL){//start new sequence (or blue bar)
        status = 'N';
    }else{
        if(rightms->segment != NULL){
            status = getLeftInfo(rightms, count);
        }else{
            if(rightms->unAligned){
                status = 'I';
                *count = rightms->gapSize;
            }else if (rightms->missingData){//NOTE!!! NEED TO COME BACK AND DEAL WITH MISSING DATA PROPERLY
                status = 'M';
                *count = rightms->gapSize;
            }else if(rightms->empty){
                status = 'N';
            }else{//deletion
                status = 'C';
            }
        }
    }
    return status;
}

void printIrow(struct MafSegment *ms, char *name, FILE *fh){
    assert(ms->segment != NULL);
    char leftStatus;
    int32_t leftCount = 0;
    char rightStatus;
    int32_t rightCount = 0;
    leftStatus = getLeftInfo( ms, &leftCount);
    rightStatus = getRightInfo( ms->next, &rightCount);
    fprintf(fh, "i\t%s\t%c\t%d\t%c\t%d\n", name, leftStatus, leftCount,
                                           rightStatus, rightCount);
    return;
}

void printErow(struct MafSegment *ms, char *name, FILE *fh){
    int32_t count = 0;
    char status = getRightInfo(ms, &count);
    fprintf(fh, "e\t%s\t%d\t%d\t%c\t%d\t%c\n", name, ms->gapStart, ms->gapSize,
                                               ms->strand, ms->srcSize, status);
}

void printMafBlockRow(struct MafSegment *mafSegment, int rownum, FILE *fh){
//void printMafBlockRow(struct MafSegment *mafSegment, char *name, FILE *fh){
    assert(mafSegment != NULL);
    Segment *segment = mafSegment->segment;
    char *name;
    if(segment != NULL){
        if (rownum < 0 ){//reference sequence
            name = mafSegment->name;
        }else{
            name = appendIntToName(mafSegment->name, rownum);
        }

        int32_t totalLen = mafSegment->srcSize;
	int32_t start = getSegmentStart(segment);
	char strand = segment_getStrand(segment) ? '+' : '-';
	int32_t len = segment_getLength(segment);//number of bases in the row
	fprintf(fh, "s\t%s\t%d\t%d\t%c\t%d\t%s\n", name, start, len, strand, totalLen, segment_getString(segment));
        printIrow(mafSegment, name, fh);
    }else{//gap, write 'e' row
        if(! mafSegment->empty ){
            name = appendIntToName(mafSegment->name, rownum);
            printErow(mafSegment, name, fh);
        }
    }
    return;
}
/*
void printRefDup(struct MafSegment *mafSegment, FILE *fh){
    assert(mafSegment != NULL);
    Segment *segment = mafSegment->segment;
    int32_t len = segment_getLength(segment);//number of bases in the row
    if(segment != NULL){
        Cap *cap = segment_get5Cap(segment);
        char *seqName = cap_getSequenceName(cap);
	int rownum = 0;
	Block *block = segment_getBlock(segment);
	Segment *currsegment;
	Block_InstanceIterator *it = block_getInstanceIterator(block);
	while((currsegment = block_getNext(it)) != NULL){
	    Cap *currcap = segment_get5Cap(currsegment);
	    char *currname = cap_getSequenceName(currcap);
            if(strcmp(currname, seqName) == 0 && cap != currcap){
	        currname = appendIntToName(seqName, rownum);
		Sequence *sequence = segment_getSequence(currsegment);
		assert(sequence != NULL);
		int32_t totalLen = sequence_getLength(sequence);//length of the chromsome the row belongs to
		int32_t start = getSegmentStart(currsegment);
		char strand = segment_getStrand(currsegment) ? '+' : '-';
		fprintf(fh, "s\t%s\t%d\t%d\t%c\t%d\t%s\n", currname, start, len, strand, totalLen, segment_getString(currsegment));
	        rownum++;
	    }
	}
	block_destructInstanceIterator(it);
    }
    return;
}*/

void printMafBlocks(struct List *refrow, int32_t c, struct List *spcRows, FILE *fh){
    int32_t i, j, h;
    for(i=0; i< refrow->length; i++){//each block
        //st_logInfo("\tColumn %d:\t", i);
        fprintf(fh, "\na\n");
        struct MafSegment *refms = refrow->list[i];
        printMafBlockRow(refms, -1, fh);//print the reference row
        //printRefDup(refms, fh);
	for(j=0; j < spcRows->length; j++){//each species
            struct List *currSpcRows = spcRows->list[j];
	    struct List *rows = currSpcRows->list[c];//correspondant row(s) to refrow
	    for(h=0; h < rows->length; h++){//each thread of current species that aligns to refrow
                //if(j==0 && h==0){continue;}
	        struct List *row = rows->list[h];
	        struct MafSegment *ms = row->list[i];
		printMafBlockRow(ms, h, fh);
            }
	}
    }
    return;
}

//=============== END PRINTING MAF FOR EACH BLOCK =================

/*void destructSpcRows(struct List *spcRows){
    for(int i=0; i < spcRows->length; i++){//for each species
        struct List *rowsList = spcRows->list[i];
        for (int j=0; j < rowsList->length; j++){//each ref-row
            struct List *rows = rowsList->list[j];
            for (int k=0; k < rows->length; k++){//each row
                //destructMyList(rows->list->list[k]);
            }
            free(rows->list);
            free(rows);
        }
        free(rowsList->list);
        free(rowsList);
    }
    free(spcRows->list);
    free(spcRows);
}*/

void getAugmentedMafs(Flower *flower, FILE *fh, char *species){
    /*
     *Get agumented mafs for the inputed flower. Agumented maf means new rows
     *or duplications are included in the maf records. Print to output file the
     *annotated mafs
     */
    struct List *spcList = splitString(species, " ");
    assert(spcList->length > 0);
    char *refSpc = spcList->list[0];

    //Get the reference row
    st_logInfo("Getting the reference Rows (%s)\n", refSpc);
    struct List *refRows = getReferenceRows(flower, refSpc);
    st_logInfo("Done. There are %d reference Rows.\n", refRows->length);
    assert(refRows->length > 0);

    struct List *spcRows = constructEmptyList(0, free);//list of rows of other species
    struct List *currSpcRows;
    if(refRows->length > 0){
        //Get rows for other species
        //for(int i=1; i < spcList->length; i++){//each species
        for(int i=0; i < spcList->length; i++){//each species
            st_logInfo("Getting Rows for %s\n", spcList->list[i]);
	    currSpcRows = getRows(flower, spcList->list[i], refRows, refSpc);
	    assert(currSpcRows->length == refRows->length);
            listAppend(spcRows, currSpcRows);
        }
        //destructList(spcList);
 
	for(int i=0; i < refRows->length; i++){//each reference row
            //st_logInfo("\nRow: %d\n", i);
            printMafBlocks(refRows->list[i], i, spcRows, fh);
            //destructMyList(refRows->list->list[i]);

            /*fprintf(fh, "\na\n");
	    printMafRow(refRows->list->list[i], fh);
	    for(int s=0; s < spcRows->length; s++){//each species
                currSpcRows = spcRows->list[s];
		struct MyList *crows = currSpcRows->list->list[i];//correspondant row(s) to current refrow (ith)
		for(int k=0; k< crows->list->length; k++){
	            struct MyList *row = crows->list->list[k];
		    printMafRow(row, fh);
		}
            }*/
	}
    }else{
        fprintf(stderr, "Could not find the reference sequence (species): %s\n", refSpc);
    }
    /*free(refRows->list);
    free(refRows);*/
    //destructSpcRows(spcRows);
    return;
}

void makeMAFHeader(Flower *flower, FILE *fileHandle) {
    fprintf(fileHandle, "##maf version=1 scoring=N/A\n");
    char *cA = eventTree_makeNewickString(flower_getEventTree(flower));
    fprintf(fileHandle, "# cactus %s\n\n", cA);
    free(cA);
}

void usage(){
    fprintf(stderr, "cactus_augmentedMaf, version 0.0\n");
    fprintf(stderr, "-a --logLevel: Set the log level\n");
    //fprintf(stderr, "-b --referenceSpecies: name of the species whose sequences to be the reference\n");
    fprintf(stderr, "-b --species: list of the species to be included in the outputed mafs.Eg: \"hg18 panTro2 rheSus2 ponAbe2\". The first species is the one whose sequences to be the reference\n");
    fprintf(stderr, "-c --cactusDisk: location of the flower disk directory\n");
    fprintf(stderr, "-d --flowerName: name of the starting flower (key in the database)\n");
    fprintf(stderr, "-e --outputFile: name of the file to write the Mafs in\n");
    fprintf(stderr, "-h --help: print this help screen\n");
}

int main(int argc, char *argv[]){
    char *logLevelString = NULL;
    char *cactusDiskDatabaseString = NULL;
    char *flowerName = NULL;
    char *species = NULL;
    char *outputFile = NULL;

    while(1){
        static struct option long_options[] = { 
	    {"logLevel", required_argument, 0, 'a'},
	    {"species", required_argument, 0, 'b'},
	    {"cactusDisk", required_argument, 0, 'c'},
	    {"flowerName", required_argument, 0, 'd'},
	    {"outputFile", required_argument, 0, 'e'},
	    {"help", no_argument, 0, 'h'},
	    {0, 0, 0, 0}
	};
	int option_index = 0;
	int key = getopt_long(argc, argv, "a:b:c:d:e:h", long_options, &option_index);
	if (key == -1){ break; }
	switch(key){
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

    assert(flowerName != NULL);
    assert(outputFile != NULL);
    assert(species != NULL);

    if (logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

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

    //if(strstr(species, "reference") != NULL){
    //    flower = flower_addReferenceSequence(flower, cactusDisk, "reference");
    //}

    //Make sure that referenceSequence has already been added:
    if(strstr(species, "reference") != NULL && getSequenceMatchesHeader(flower, "reference") == NULL){
        fprintf(stderr, "No reference sequence found in cactusDisk\n");
        exit(EXIT_FAILURE);
    }

    FILE *fh = fopen(outputFile, "w");
    makeMAFHeader(flower, fh);
       
    getAugmentedMafs(flower, fh, species);
    fprintf(fh, "\n");

    fclose(fh);
    st_logInfo("Got the mafs in %i seconds/\n", time(NULL) - startTime);

    ///////////////////////////////////////////////////////////////////////////
    // Clean up.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}

