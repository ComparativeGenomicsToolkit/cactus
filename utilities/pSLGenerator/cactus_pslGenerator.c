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
#include "common.h"
#include "psl.h"
/*
 * nknguyen@soe.ucsc.edu
 * 04/28/2010
 * This script generates pair-wise alignment for 
 */

//====================== GLOBAL STRUCTURES ====================================
struct Thread{
   Cap **caps;//array of caps in order
   int size; //number of caps
};

struct Align{
   int *qIndices;//array of q caps indices (in qthread->caps) in the alignment
   int *tIndices;//array of t caps indices (in tthread->caps) in the alignment
   int size; //number of segments in the alignment
};

struct Stubs{
   Cap **qstubs;
   Cap **tstubs;
   int qnum;
   int tnum;
};

//========================== PROTOTYPES =======================================
bool isStubCap(Cap *cap);
int getPSL(struct Align *align, struct Thread *qThread, struct Thread *tThread, char *query, char *target, FILE *fileHandle);
void getPSLNet(FILE *fileHandle, char *query, char *target, int start, int end, Cap *qstartCap, bool exhaust);
Cap *net_getChildCap(Net *net, Cap *pcap);
//int getPSL(Net *net, struct Align *align, struct Thread *qThread, struct Thread *tThread, char *query, char *target, FILE *fileHandle, int startI);
//void getPSLNet(Net *net, FILE *fileHandle, char *query, char *target, int start, int end, Cap *qstartCap);

//========================= INITIALIZATION FUNCTIONS ==========================
struct Thread *setThread(){
   struct Thread *thread = AllocA(struct Thread);
   thread->size = 0;
   return thread;
}
struct Align *setAlign(int *qI, int *tI, int size){
   struct Align *align = AllocA(struct Align);
   align->qIndices = qI;
   align->tIndices = tI;
   align->size = size;
   return align;
}
struct Stubs *setStubs(){
   struct Stubs *stubs = AllocA(struct Stubs);
   stubs->qnum = 0; 
   stubs->tnum = 0; 
   return stubs;
}
void resetStubs(struct Stubs **stubs){
   (*stubs)->qnum = 0;
   (*stubs)->tnum = 0;
}
void pslInitialize(struct psl *psl){
        psl->match = 0;
        psl->misMatch = 0;
        psl->repMatch = 0;
        psl->nCount = 0;
        psl->qNumInsert = 0;
        psl->qBaseInsert = 0;
        psl->tNumInsert = 0;
        psl->tBaseInsert = 0;
        psl->blockCount = 0;
   return;
}
//============================ PRINT FUNCTIONS ================================
//For testing purpose - not necessary for main program
void printCap(FILE *fileHandle, Cap *cap){
   fprintf(fileHandle, "\t%s", netMisc_nameToString(cap_getName(cap)));
   if(isStubCap(cap)){ return; }
   char side = cap_getSide(cap) ? '5' : '3';
   char strand = cap_getStrand(cap) ? '+' : '-';
   int coor = cap_getCoordinate(cap);
   Segment *segment = cap_getSegment(cap);
   Cap *otherCap = cap_getSide(cap) ? segment_get3Cap(segment): segment_get5Cap(segment);
   int otherCoor = cap_getCoordinate(otherCap);
   fprintf(fileHandle, " %c%c%d => %d\t", side, strand, coor, otherCoor);
}
void printStubs(FILE *fileHandle, struct Stubs *stubs){
   int i;
   Cap *cap;
   fprintf(fileHandle, "Query:\n");
   fprintf(fileHandle, "3-end Stubs, %d total:\t", stubs->qnum);
   for(i=0; i< stubs->qnum; i++){
      cap = *(stubs->qstubs + i);
      fprintf(fileHandle, "\t");
      printCap(fileHandle, cap);
      fprintf(fileHandle, " end%s", netMisc_nameToString(end_getName(cap_getEnd(cap))));
   }
   fprintf(fileHandle, "\n");
   fprintf(fileHandle, "Target:\n");
   fprintf(fileHandle, "3-end Stubs, %d total:\t", stubs->tnum);
   for(i=0; i< stubs->tnum; i++){
      cap = *(stubs->tstubs + i);
      fprintf(fileHandle, "\t");
      printCap(fileHandle, cap);
      fprintf(fileHandle, " end%s", netMisc_nameToString(end_getName(cap_getEnd(cap))));
   }
   fprintf(fileHandle, "\n");
}
void printThread(FILE *fileHandle, struct Thread *thread){
   fprintf(fileHandle, "printing Thread\n");
   if(thread == NULL){
      fprintf(fileHandle, "NULL\n");
   }
   int i;
   fprintf(fileHandle, "Thread size %d\n", thread->size);
   for(i=0; i< thread->size; i++){
      Cap *cap = *(thread->caps +i);
      //Segment *segment = cap_getSegment(cap);
      //Block *block = segment_getBlock(segment);
      //fprintf(fileHandle, "\t%s", netMisc_nameToString(block_getName(block)));
      fprintf(fileHandle, "\t%s", netMisc_nameToString(end_getName(cap_getEnd(cap))));
   }
   fprintf(fileHandle, "\nCaps:\t");
   for(i=0; i< thread->size; i++){
      Cap *cap = *(thread->caps +i);
      fprintf(fileHandle, "\t");
      printCap(fileHandle, cap);
   }
   fprintf(fileHandle, "\n");
}
void printAlign(FILE *fileHandle, struct Align *align){
   fprintf(fileHandle, "printing Align\n");
   fprintf(fileHandle, "Align size %d\n", align->size);
   int i;
   fprintf(fileHandle, "qIndices:");
   for(i=0; i< align->size; i++){
      fprintf(fileHandle, "\t%d", *(align->qIndices +i));
   }
   fprintf(fileHandle, "\n");
   fprintf(fileHandle, "tIndices:");
   for(i=0; i< align->size; i++){
      fprintf(fileHandle, "\t%d", *(align->tIndices +i));
   }
   fprintf(fileHandle, "\n");
}

//============================= UTILS FUNCTIONS ========================
void makePSLHeader(Net *net, FILE *fileHandle) {
   pslWriteHead(fileHandle);
} 

char *formatSequenceHeader(Sequence *sequence) {
   const char *sequenceHeader = sequence_getHeader(sequence);
   if(strlen(sequenceHeader) > 0) {
      char *cA = malloc(sizeof(char) *(1 + strlen(sequenceHeader)));
      sscanf(sequenceHeader, "%s", cA);
      return cA;
   }
   else {
      return netMisc_nameToString(sequence_getName(sequence));
   }
}

Sequence *net_getSequenceByName(Net *net, char *name){
   Sequence *sequence;
   Net_SequenceIterator * seqIterator = net_getSequenceIterator(net);
   while((sequence = net_getNextSequence(seqIterator)) != NULL){
      char *sequenceHeader = formatSequenceHeader(sequence);
      //if 'sequenceHeader' starts with 'name'
      if(strstr(sequenceHeader, name) == sequenceHeader){
         free(sequenceHeader);
         return sequence;
      }
   }
   net_destructSequenceIterator(seqIterator);
   return NULL;
}

int getSequences(Net *net, char ***seqs, char *name){
   //get all the sequences in net that have name start with 'name'
   int num = 0;
   Sequence *sequence;
   Net_SequenceIterator * seqIterator = net_getSequenceIterator(net);
   while((sequence = net_getNextSequence(seqIterator)) != NULL){
      char *sequenceHeader = formatSequenceHeader(sequence);
      //if 'sequenceHeader' starts with 'name'
      if(strstr(sequenceHeader, name) == sequenceHeader){
         if(num == 0){
            (*seqs) = AllocA(char *);
         }else{
            (*seqs) = needMoreMem((*seqs), num*sizeof(char*), (num+1)*sizeof(char *));
         } 
         *((*seqs) + num) = sequenceHeader;
         num++;
      }
      //free(sequenceHeader);
   }
   net_destructSequenceIterator(seqIterator);
   return num;
}

bool isStubCap(Cap *cap){
   End *end = cap_getEnd(cap);
   return (end_isStubEnd(end)) ? true : false;
}

bool isPlus(char strand){
   return (strand == '+')? true : false;
}

bool isCommonBlock(Block *block, char *query, char *target){
   //Return true if block contains both query and target
   if(block == NULL){ return false; }
   bool hasT = false;
   bool hasQ = false;
   Segment *segment;
   Block_InstanceIterator *segmentIterator= block_getInstanceIterator(block);
   while((segment= block_getNext(segmentIterator)) != NULL){
      Sequence *sequence = segment_getSequence(segment);
      if(sequence == NULL){
         continue;
      }
      char *sequenceHeader = formatSequenceHeader(sequence);
      if(strcmp(sequenceHeader, query) == 0){
         hasQ = true;
      }else 
      if(strcmp(sequenceHeader, target) == 0){
         hasT = true;
      }
         //free(sequenceHeader);
   }
   block_destructInstanceIterator(segmentIterator);
   if (hasQ && hasT){
      return true;
   }else{
      return false;
   }
}

void moveCapToNextBlock(Cap **cap){
   /*Move cap to the next block that its segment is adjacency to*/
   Cap *adjCap = cap_getAdjacency(*cap);
   if(cap_getEnd(*cap) == cap_getEnd(adjCap)){//self connected end
      *cap = adjCap;
   }else{
      if(cap_getSide(*cap)){//5' end
         *cap = cap_getAdjacency(segment_get3Cap(cap_getSegment(*cap)));
      }else{
         *cap = cap_getAdjacency(segment_get5Cap(cap_getSegment(*cap)));
      }
   }
}

void moveCapToPrevBlock(Cap **cap){
   /*Move cap to the next block that its segment is adjacency to*/
   Cap *adjCap = cap_getAdjacency(*cap);
   if(cap_getEnd(*cap) == cap_getEnd(adjCap)){//self connected end
      *cap = adjCap;
   }else{
      if(cap_getSide(adjCap)){//5' end
         *cap = segment_get3Cap(cap_getSegment(adjCap));
      }else{
         *cap = segment_get5Cap(cap_getSegment(adjCap));
      }
   }
}

void moveCapToNextCommonBlock(Cap **cap, char *name){
   bool commonBlock = false;
   Cap *currCap;
   while(!commonBlock && !isStubCap(*cap)){
      moveCapToNextBlock(cap);
      End *end = cap_getEnd(*cap);
      End_InstanceIterator *capIterator= end_getInstanceIterator(end);
      while((currCap= end_getNext(capIterator)) != NULL){
         Sequence *sequence = cap_getSequence(currCap);
         if(sequence == NULL){ continue; }
         char *sequenceHeader = formatSequenceHeader(sequence);
	 if(strcmp(sequenceHeader, name) == 0){
	    free(sequenceHeader);
	    commonBlock = true;
	    break;
	 }
      }
      end_destructInstanceIterator(capIterator);
   }
   return;
}

bool isMatch(struct Thread *qThread, struct Thread *tThread, int qi, int ti){
   /*Return true if qsegment & tsegment come from same block, false otherwise*/
   Cap *qcap = *(qThread->caps + qi);
   Cap *tcap = *(tThread->caps + ti);
   End *qend = cap_getEnd(qcap);
   End *tend = cap_getEnd(tcap);
   if(qend == NULL || tend == NULL){
      return false;
   }
   return qend == tend ? true : false;
}

void fixNegStrandPsl(struct psl *psl){
   int32_t temp;
   //convert the negative strand psl into correct coordinates
   if(!isPlus(psl->strand[0])){//query is on negative strand
      temp = psl->qStart + 1;
      psl->qStart = psl->qEnd + 1;
      psl->qEnd = temp;
   }
   if(!isPlus(psl->strand[1])){//target is on negative strand
      temp = psl->tStart + 1;
      psl->tStart = psl->tEnd + 1;
      psl->tEnd = temp;
   }
   return; 
}

void addStub(Cap ***caps, Cap *cap, int *num){
   (*num)++;
   if(*num == 1){
      *caps = AllocA(Cap *);
      **caps = cap;
   }else{
      *caps = needMoreMem(*caps, (*num-1)*sizeof(Cap *), (*num)*sizeof(Cap *));
      *(*caps + *num -1) = cap;
   }
   return;
}

//========== CONSTRUCTING A PSL - ADDING THE PARTS ============================
int32_t getBlockStart(Cap *cap){
   /*Get the start coordinate of input segment
    *psl starts from 0, while catus start at an added stub, which is at position 1
    *Therefore need to substract 2 to convert to the right coordinate
    */
   int32_t start = cap_getCoordinate(cap);//coordinate of cap on positive strand
   if(cap_getStrand(cap)){
      return start -2;
   }else{
      Sequence *sequence = cap_getSequence(cap);
      int32_t seqLen = sequence_getLength(sequence);
      return seqLen - 1 - (start -2);
   }
}

void addPSLSizes(struct psl *psl, Cap *qcap, Cap *tcap){
   Sequence *qsequence = cap_getSequence(qcap);
   Sequence *tsequence = cap_getSequence(tcap);
   psl->qSize = sequence_getLength(qsequence);
   psl->tSize = sequence_getLength(tsequence);
}

void addPSLStarts(struct psl *psl, Cap *qcap, Cap *tcap){
   //Get threads' starts relatively to the forward strands
   int32_t qstart = cap_getCoordinate(qcap);
   int32_t tstart = cap_getCoordinate(tcap);
   psl->qStart = qstart-2;
   psl->tStart = tstart-2;
}

void addPSLStrand(struct psl *psl, Cap *qcap, Cap *tcap){
   psl->strand[0] = cap_getStrand(qcap) ? '+' : '-';
   psl->strand[1] = cap_getStrand(tcap) ? '+' : '-';
   //psl->strand[2] = NULL;
}

void addPSLEnds(struct psl *psl){
   int32_t i = psl->blockCount -1;
   int32_t qend = *(psl->qStarts +i) + *(psl->blockSizes +i);
   int32_t tend = *(psl->tStarts +i) + *(psl->blockSizes +i);
   psl->qEnd = isPlus(psl->strand[0]) ? qend : psl->qSize - 1 - qend; 
   psl->tEnd = isPlus(psl->strand[1]) ? tend : psl->tSize - 1 - tend; 
}

void addPSLInserts(struct psl *psl, int32_t currIndex){
   int32_t qdiff;
   int32_t tdiff;
   if(currIndex == 0){
      return;
   }
   qdiff = *(psl->qStarts + currIndex) - *(psl->qStarts + currIndex -1) - *(psl->blockSizes + currIndex -1);
   tdiff = *(psl->tStarts + currIndex) - *(psl->tStarts + currIndex -1) - *(psl->blockSizes + currIndex -1);
   if(qdiff > 0){
      psl->qNumInsert ++;
      psl->qBaseInsert += qdiff;
   }
   if(tdiff > 0){
      psl->tNumInsert ++;
      psl->tBaseInsert += tdiff;
   }
}

void addBlockToPSL(struct psl *psl, Cap *qcap, Cap *tcap){
   /*
    */
   psl->blockCount++;
   Segment *segment = cap_getSegment(qcap);
   int32_t length = segment_getLength(segment);
   int32_t currIndex = psl->blockCount -1;
   //blockSizes
   psl->blockSizes = (unsigned *)realloc(psl->blockSizes, psl->blockCount*sizeof(unsigned));
   *(psl->blockSizes + currIndex) = length;
   //qStarts
   psl->qStarts = (unsigned *)realloc(psl->qStarts, psl->blockCount*sizeof(unsigned));
   *(psl->qStarts + currIndex) = getBlockStart(qcap);
   //tStarts
   psl->tStarts = (unsigned *)realloc(psl->tStarts, psl->blockCount*sizeof(unsigned));
   *(psl->tStarts + currIndex) = getBlockStart(tcap);
   psl->match += length;
   //Inserts
   addPSLInserts(psl, currIndex);
}

void addPSLBlocks(struct psl *psl, struct Align *align, struct Thread *qThread, struct Thread *tThread, FILE *fileHandle){
   //return number of blocks get added
   int i, qi, ti;
   Cap *qcap;
   Cap *tcap;
   bool qStrand;
   bool tStrand;
   char pslStrands[2];
   pslStrands[0] = psl->strand[0];
   pslStrands[1] = psl->strand[1];
   
   i = 0;
   while(i <align->size){
      qi = *(align->qIndices + i);
      ti = *(align->tIndices + i);
      qcap = *(qThread->caps + qi);
      tcap = *(tThread->caps + ti);
      qStrand = cap_getStrand(qcap);
      tStrand = cap_getStrand(tcap);
      addBlockToPSL(psl, qcap, tcap);//add to qstarts, tstarts   
      i++;
   }
   return;
}

//int getPSL(Net *net, struct Align *align, struct Thread *qThread, struct Thread *tThread, char *query, char *target, FILE *fileHandle, int startI){
int getPSL(struct Align *align, struct Thread *qThread, struct Thread *tThread, char *query, char *target, FILE *fileHandle){
   /*
    *Get PSL
    */
   //char *pslDesc = NULL;
   //FILE *err = fopen("dump", "w");
   if(align == NULL || qThread == NULL || tThread == NULL){
      return 0;
   }
   int blockSpace = 32; //default number of blocks 
   char strand[3];
   //char *netName = netMisc_nameToString(net_getName(net));
   //struct psl * psl = pslNew(netName, 0, 0, 0, target, 0, 0, 0, strand, blockSpace, 0);
   struct psl * psl = pslNew(query, 0, 0, 0, target, 0, 0, 0, strand, blockSpace, 0);
   pslInitialize(psl);
   
   int numBlocks = 0;
   int qi = align->qIndices[0];
   int ti = align->tIndices[0];
   Cap *qcap = *(qThread->caps + qi);
   Cap *tcap = *(tThread->caps + ti);
   addPSLSizes(psl, qcap, tcap);
   addPSLStarts(psl, qcap, tcap);
   addPSLStrand(psl, qcap, tcap);
   addPSLBlocks(psl, align, qThread, tThread, fileHandle);
   addPSLEnds(psl);
   if(psl->blockCount >= 1){
      char sep = '\t';
      char lastStep = '\n';
      fixNegStrandPsl(psl);
      if(!isPlus(psl->strand[1])){
         pslRc(psl);
      }
      //if(pslCheck(pslDesc, err, psl) == 0){
         pslOutput(psl, fileHandle, sep, lastStep);
      //}
   }
   numBlocks = psl->blockCount;
   pslFree(&psl);
   return numBlocks;
}

//=============================== ALIGNING THREADS ============================
void initializeIndices(int **qIndices, int **tIndices, int size, int qi, int ti){
   (*qIndices) = AllocN(int, size);
   (*tIndices) = AllocN(int, size);
   *((*qIndices) + size -1) = qi;
   *((*tIndices) + size -1) = ti;
   return;
}
/*struct Align **alignThreads(struct Thread *qThread, struct Thread *tThread, int qi, int ti, int *size){
   struct Align **aligns = NULL;
   struct Align **prevAligns;
   int *qIndices;
   int *tIndices;
   int curr_ti, k, l;
   int count = 0;
   int prevSize = 0;
   bool checkMatch;
   if(qi == qThread->size -1 || ti == tThread->size -1){//at the end of one (or both) thread(s)
      if(isMatch(qThread, tThread, qi, ti)){
         aligns = AllocA(struct Align *);
         qIndices = AllocA(int);
         tIndices = AllocA(int);
	 *qIndices = qi;
	 *tIndices = ti;
         *aligns = setAlign(qIndices, tIndices, 1);
	 *size = 1;
	 return aligns;
      }else{
         return NULL;
      }
   }
   for(curr_ti = ti; curr_ti < tThread->size - 1; curr_ti++){
      checkMatch = isMatch(qThread, tThread, qi, curr_ti);
      if(checkMatch){
         prevSize = 0;
         prevAligns = alignThreads(qThread, tThread, qi + 1, curr_ti + 1, &prevSize);
	 if(prevSize == 0){
            if(count == 0){
               aligns = AllocA(struct Align *);
	    }else{
               aligns = needMoreMem(aligns, count*sizeof(struct Align *), (count +1)*sizeof(struct Align *));
	    }
            qIndices = AllocA(int);
            tIndices = AllocA(int);
	    *qIndices = qi;
	    *tIndices = ti;
	    *(aligns + count) = setAlign(qIndices, tIndices, 1);
            count += 1;
         }else{
            if(count == 0){
               aligns = AllocN(struct Align *, prevSize);
	    }else{
               aligns = needMoreMem(aligns, count*sizeof(struct Align *), (count +prevSize)*sizeof(struct Align *));
	    }
	    for(k=0; k < prevSize; k++){
	       struct Align *pA = *(prevAligns +k);
               qIndices = AllocN(int, 1+pA->size);
               tIndices = AllocN(int, 1+pA->size);
	       *qIndices = qi;
	       *tIndices = curr_ti;
               for(l=0; l< pA->size ; l++){
	          *(qIndices +l +1) = *(pA->qIndices + l);
	          *(tIndices +l +1) = *(pA->tIndices + l);
	       }
	       *(aligns + count +k) = setAlign(qIndices, tIndices, pA->size +1);
	    }
	    count += prevSize;
         }
         freeMem(prevAligns);
      }
   }
   *size = count;
   return aligns;
}*/
void aligns_getMem(struct Align ***aligns, int oldsize, int addsize){
   if(oldsize == 0){
      (*aligns) = AllocN(struct Align *, addsize);
   }else{
      (*aligns) = needMoreMem((*aligns), oldsize*sizeof(struct Align *), (oldsize + addsize)*sizeof(struct Align *));
   }
   return;
}
bool align_compare(struct Align *a1, struct Align *a2){
   //If either a1 or a2 is NULL, return false
   if(a1 == NULL || a2 == NULL){
      return false;
   }
   int i, q1, q2, t1, t2;
   bool isEqual = true;
   if(a1->size == a2->size){ 
      for(i=0; i< a1->size; i++){
         q1 = *(a1->qIndices + i);
         q2 = *(a2->qIndices + i);
         t1 = *(a1->tIndices + i);
         t2 = *(a2->tIndices + i);
         if(q1 != q2 || t1 != t2){
	    isEqual = false;
	    break;
         } 
      }
   }else{
      isEqual = false;
   }
   return isEqual;
}
struct Align *align_getBest(struct Align **aligns, int num, struct Thread *qthread){
   int i, j, b;
   int numbase = 0;
   int max = 0;
   struct Align *align;
   struct Align *best;
   Cap *cap;
   Cap **caps = qthread->caps;
   //fprintf(stderr, "Getting Best alignment\n");
   for(i=0; i< num; i++){//each align
      align = *(aligns + i);
      //if(align == NULL){ continue; } 
      numbase = 0;
      for(j=0; j< align->size; j++){//each block in the align
         cap = *(caps + *(align->qIndices + j));
         Segment *segment = cap_getSegment(cap);
         numbase += segment_getLength(segment);
      }
      if(max < numbase){
         max = numbase;
         best = align;
	 b = i;
      }
   }
   return best;
}

struct Align *alignThreads(struct Thread *qThread, struct Thread *tThread){
   /* Dynamic programming - linear space to get the best alignment of qThread and tThread*/
   struct Align ***pA = NULL;//alignment of last & current rows
   struct Align *align;
   struct Align **aligns;
   int *qIndices;
   int *tIndices;
   int i1, i2, i, j, qi, ti, prevSize, count, k;
   int qsize = qThread->size;
   int tsize = tThread->size;
   
   if(qsize < 1 || tsize < 1){ return NULL; }
   
   //Allocate memery for pA
   pA = AllocN(struct Align **, 2);
   *pA = AllocN(struct Align *, tsize + 1);//start at position -1 of target
   *(pA + 1) = AllocN(struct Align *, tsize + 1);
   //Initialize pA
   for(i = 0; i < 2; i++){
      for(j = 0; j < tsize + 1; j++){//start at position -1 (nothing) of query
         *(*(pA + i) + j) = NULL;
      }
   }
   //Aligning
   for(i = 1; i < qsize + 1; i++){
      i1 = (i-1)%2; //i - 1
      i2 = i%2; //i
      qi = i-1;
      for(j = 1; j < tsize + 1; j++){
	 ti = j -1;
	 count = 0;
	 aligns = NULL;
	 //q0..qi-1, t0...ti-1
         align = *(*(pA + i1) + j-1); 
         if(isMatch(qThread, tThread, qi, ti)){
            if(align == NULL){
               prevSize = 0;
            }else{
               prevSize = align->size;
            }
            aligns_getMem(&aligns, count, 1);
            qIndices = NULL;
            tIndices = NULL;
            initializeIndices(&qIndices, &tIndices, 1 + prevSize, qi, ti);
            for(k=0; k < prevSize ; k++){
	       *(qIndices + k) = *(align->qIndices + k);
	       *(tIndices + k) = *(align->tIndices + k);
            }
            *(aligns + count) = setAlign(qIndices, tIndices, prevSize +1);
            count += 1;
	 }else{
            if(align != NULL){
               aligns_getMem(&aligns, count, 1);
               *(aligns + count) = align;
               count += 1;
            }
	 }
	 //q0..qi-1, t0...ti
         align = *(*(pA + i1) + j);
	 if(align != NULL){
            aligns_getMem(&aligns, count, 1);
            *(aligns + count) = align;
            count += 1;
	 }
	 //q0..qi, t0...ti-1
         align = *(*(pA + i2) + j-1);
	 if(align != NULL){
            aligns_getMem(&aligns, count, 1);
            *(aligns + count) = align;
            count += 1;
	 }
         if(count > 0){
	    *(*(pA + i2) + j) = align_getBest(aligns, count, qThread);
         }else{
	    *(*(pA + i2) + j) = NULL;
         }
      }
   }
   return *(*(pA + qsize%2) + tsize);
}

struct Align **alignThreads_exhaust(struct Thread *qThread, struct Thread *tThread, int *size){
   fprintf(stderr, "Exhaustly aligning ...\n");
   struct Align ****pA = NULL;//alignment of last & current rows
   int **sizes;
   struct Align **aligns;
   struct Align **prevAligns;
   struct Align *align;
   int *qIndices;
   int *tIndices;
   int i1, i2, i, j, qi, ti, prevSize, count, k, l;
   int qsize = qThread->size;
   int tsize = tThread->size;
   
   if(qsize < 1 || tsize < 1){ return NULL; }
   
   //Allocate memory for pA
   pA = AllocN(struct Align ***, 2);
   *pA = AllocN(struct Align **, tsize + 1);//start at position -1 of target
   *(pA + 1) = AllocN(struct Align **, tsize + 1);
   //Initialize pA
   for(i = 0; i < 2; i++){
      for(j = 0; j < tsize + 1; j++){//start at position -1 (nothing) of query
         *(*(pA + i) + j) = NULL;
      }
   }
   //Allocate memory for sizes
   sizes = AllocN(int*, 2);
   *sizes = AllocN(int, tsize + 1);//start at position -1 of target
   *(sizes + 1) = AllocN(int, tsize + 1);
   
   //Aligning
   for(i = 1; i < qsize + 1; i++){
      i1 = (i-1)%2; //i - 1
      i2 = i%2; //i
      qi = i-1;
      for(j = 1; j < tsize + 1; j++){
	 ti = j -1;
	 aligns = NULL;
	 count = 0;
	 //q0..qi-1, t0...ti-1
         prevAligns = *(*(pA + i1) + j-1); 
	 prevSize = *(*(sizes + i1) + j-1);
         fprintf(stderr, "qi, ti %d %d, prevSize (qi-1, ti-1)= %d\n", qi, ti, prevSize);
         if(isMatch(qThread, tThread, qi, ti)){
            if(prevSize == 0){ 
               aligns_getMem(&aligns, count, 1);
               initializeIndices(&qIndices, &tIndices, 1, qi, ti);
                (*aligns) = setAlign(qIndices, tIndices, 1);
               count++;
            }else{
               aligns_getMem(&aligns, count, prevSize);
	       for(l = 0; l < prevSize; l++){
	          align = *(prevAligns + l);
                  initializeIndices(&qIndices, &tIndices, 1 + align->size, qi, ti);
                  for(k=0; k < align->size ; k++){
	             *(qIndices + k) = *(align->qIndices + k);
	             *(tIndices + k) = *(align->tIndices + k);
                  }
                  *(aligns + l) = setAlign(qIndices, tIndices, align->size +1);
               }
	    count += prevSize;
            }
	 }else{
            if(prevSize > 0){
               aligns_getMem(&aligns, count, prevSize);
	       aligns = prevAligns;
               count += prevSize;
            }
	 }
	 //q0..qi-1, t0...ti
         prevAligns = *(*(pA + i1) + j);
	 prevSize = *(*(sizes + i1) + j);
         fprintf(stderr, "qi -1, ti, size = %d\n", prevSize);
	 if(prevSize > 0){
	    for(l = 0; l < prevSize; l++){
	       align = *(prevAligns + l);
	       for(k = 0; k < count; k++){
	          if(align_compare(align, *(aligns + k))){//already added this alignment
		     break;
		  }
	       }
	       if(k == count){//haven't added this alignment
                  aligns_getMem(&aligns, count, 1);
                  *(aligns + count) = align;
		  count ++;
	       }
	    }
	 }
	 //q0..qi, t0...ti-1
         prevAligns = *(*(pA + i2) + j -1);
	 prevSize = *(*(sizes + i2) + j -1);
         fprintf(stderr, "qi, ti -1, size = %d\n", prevSize);
	 if(prevSize > 0){
	    for(l = 0; l < prevSize; l++){
	       align = *(prevAligns + l);
	       for(k = 0; k < count; k++){
	          if(align_compare(align, *(aligns + k))){//already added this alignment
		     break;
		  }
	       }
	       if(k == count){//haven't added this alignment
                  aligns_getMem(&aligns, count, 1);
                  *(aligns + count) = align;
		  count ++;
	       }
	    }
	 }
         if(count > 0){
	    *(*(pA + i2) + j) = aligns;
         }else{
	    *(*(pA + i2) + j) = NULL;
         }
         *(*(sizes + i2) + j) = count;
      }
   }
   *size = *(*(sizes + qsize%2) + tsize);
   return *(*(pA + qsize%2) + tsize);
}

//========================== GETTING THREADS ==================================
End *end_getOppEnd(End *end){
   /*Get the Other End of the block end belongs to*/
   End *oppEnd = NULL;
   Block *block = end_getBlock(end);
   if(block != NULL){
      if(end == block_get5End(block)){
         oppEnd = block_get3End(block);
      }else{
         oppEnd = block_get5End(block);
      }
   }
   return oppEnd;
}
Cap *cap_getOppCap(Cap *cap){
   /*Get the Other End of the block end belongs to*/
   Cap *oppCap = NULL;
   Segment *segment = cap_getSegment(cap);
   if(segment != NULL){
      if(cap_getSide(cap)){//5' end
         oppCap = segment_get3Cap(segment);
      }else{
         oppCap = segment_get5Cap(segment);
      }
   }
   return oppCap;
}

End *traverseQuery(Cap *cap, struct Thread **thread, char *query, char *target, int start, int end, FILE *fileHandle, bool exhaust){
   cap = cap_getAdjacency(cap);
   int coor;
   bool past = false;
   Cap *prevcap = cap;
   int size = (*thread)->size;
   End *startEnd = NULL; 
   if(!isCommonBlock(end_getBlock(cap_getEnd(cap)), query, target)){
      moveCapToNextCommonBlock(&cap, target);
   }
   //traverse the net, start from 'cap', to get the thread
   while(!isStubCap(cap)){
      End *cend = cap_getEnd(cap);
      Block *block = end_getBlock(cend);
      int blockLen = block_getLength(block);
      coor = cap_getCoordinate(cap);//Cap coordinate is always the coordinate on + strand
      if(coor >= end){ //past range
         if(past){//if already visited this cap before, break
            break;
         }else{//if past the range, but hasn't visited any lower nets
            cap = prevcap;
            past = true;
            continue;
         } 
      }
      if(coor + blockLen < start){ //curr cap is upstream of 'start'
         if(!past){//hasn't visited this cap yet
            prevcap = cap;
            moveCapToNextCommonBlock(&cap, target);
	    continue;
         }
      }
      past = true;
      if(block_getChain(block)){//if block belongs to a chain
	 if(startEnd == NULL){
	    startEnd = cend;
	 }
         if(size == 0){
            (*thread)->caps = AllocA(Cap *);
         }else{
            (*thread)->caps = needMoreMem((*thread)->caps, size*sizeof(Cap *), (size+1)*sizeof(Cap *));
	 }
	 *((*thread)->caps + size) = cap;
	 size++;
         //Traverse lower level nets if exists 
	 Group *group = end_getGroup(end_getOppEnd(cend));
	 Net *nestedNet = group_getNestedNet(group);
         net_check(nestedNet);
	 if(nestedNet != NULL){//recursive call
            Cap *childCap = net_getChildCap(nestedNet, cap_getOppCap(cap));
            if(childCap != NULL){
               getPSLNet(fileHandle, query, target, start, end, childCap, exhaust);
	    }
	 }
      }
      moveCapToNextCommonBlock(&cap, target);
   }
   (*thread)->size = size;
   return startEnd;
}
bool cmpWithEnd(Cap *cap, char *name, int end){
   bool isLess = false;
   Cap *currCap;
   End_InstanceIterator *capIterator = end_getInstanceIterator(cap_getEnd(cap));
   while((currCap = end_getNext(capIterator)) != NULL){
      Sequence *sequence = cap_getSequence(cap);
      char *sequenceHeader = formatSequenceHeader(sequence);
      if(strcmp(sequenceHeader, name) == 0){
         if(cap_getCoordinate(currCap) < end){
	    isLess = true;
	    break;
	 }
      }   
   }
   end_destructInstanceIterator(capIterator);
   return isLess;
}
void traverseTarget(Cap *cap, struct Thread **thread, char *query, char *target, int start, int end){
   int size = (*thread)->size;
   //traverse the net, start from 'cap', to get the thread
   while(!isStubCap(cap)){
      End *cend = cap_getEnd(cap);
      Block *block = end_getBlock(cend);
      if(!cmpWithEnd(cap, target, end)){ break; } //past range
      if(block_getChain(block)){//if block belongs to a chain
         if(size == 0){
            (*thread)->caps = AllocA(Cap *);
         }else{
            (*thread)->caps = needMoreMem((*thread)->caps, size*sizeof(Cap *), (size+1)*sizeof(Cap *));
	 }
	 *((*thread)->caps + size) = cap;
	 size++;
      }
      moveCapToNextCommonBlock(&cap, query);
   }
   (*thread)->size = size;
   return;
}

//=================== GET START OF THREAD (3' STUB) ====================================
Cap *net_getThreadStart(Net *net, char *name){// from Net = 0
   //Get 3' end Stubs of the sequence by its name
   Cap *cap;
   Net_CapIterator *capIterator = net_getCapIterator(net);
   while((cap= net_getNextCap(capIterator)) != NULL){
      if(isStubCap(cap) && !cap_getSide(cap)){//3' dead end or inherited end
         Sequence *sequence = cap_getSequence(cap);
         if(sequence == NULL){continue;}
         char *sequenceHeader = formatSequenceHeader(sequence);
         if(strcmp(sequenceHeader, name) == 0){
            break;
         }
         free(sequenceHeader);
      }
   }
   net_destructCapIterator(capIterator);
   return cap;
}

Cap *net_getChildCap(Net *net, Cap *pcap){
   //Find stub in net that is correspond to input pcap, which is from parent net
   Cap *cap;
   Net_CapIterator *capIterator = net_getCapIterator(net);
   while((cap= net_getNextCap(capIterator)) != NULL){
      if(isStubCap(cap)){//inherited end
         if(cap_getName(pcap) == cap_getName(cap)){
            net_destructCapIterator(capIterator);
	    //Always return the cap on the positive strand
	    return cap_getStrand(cap) ? cap : cap_getReverse(cap);
         }
      }
   }
   net_destructCapIterator(capIterator);
   return NULL;
}

//===================== GETTING PSLs FROM CHAINS FOR CURRENT NET ==============
//void getPSLNet(Net *net, FILE *fileHandle, char *query, char *target, int start, int end, Cap *qstartCap){
void getPSLNet(FILE *fileHandle, char *query, char *target, int start, int end, Cap *qstartCap, bool exhaust){
   /*
    *Get PSLs for net, current level
    */
   //Each thread is an array of ordered caps obtained by traversing the net
   int size;
   int i=0;
   Cap *cap;
   //char *netName = netMisc_nameToString(net_getName(net));
   //fprintf(stderr, "\nNET: %s\n", netName);
   struct Thread *qThread = setThread();
   struct Align *align;
   End *startEnd = traverseQuery(qstartCap, &qThread, query, target, start, end, fileHandle, exhaust);
   if(startEnd == NULL){
      return;
   }
   End_InstanceIterator *capIterator = end_getInstanceIterator(startEnd);
   while((cap = end_getNext(capIterator)) != NULL){
      Sequence *sequence = cap_getSequence(cap);
      if(sequence == NULL){ continue; }
      char *sequenceHeader = formatSequenceHeader(sequence);
      if(strcmp(sequenceHeader, target) == 0){
         struct Thread *tThread = setThread();
         traverseTarget(cap, &tThread, query, target, start, end);
         i++;
         //fprintf(stderr, "qThreadSize = %d, tThreadSize = %d\n", qThread->size, tThread->size);
         if(!exhaust){
            align = alignThreads(qThread, tThread);
            if(align != NULL){
               getPSL(align, qThread, tThread, query, target, fileHandle);
	       freeMem(align);
            }
	 }else{
            size = 0;
            struct Align **aligns;
            aligns = alignThreads_exhaust(qThread, tThread, &size);
            int a;
            for(a = 0; a < size; a++){
               //printAlign(stderr, *(aligns +a));
                  getPSL(*(aligns + a), qThread, tThread, query, target, fileHandle);
            }
	    if(size > 0){ freeMem(aligns); }
         }
         if(tThread != NULL){ freeMem(tThread); }
      }
   }
   if(qThread != NULL){ freeMem(qThread); }
   end_destructInstanceIterator(capIterator);
   return;
}

//========================= TANGLE PSLs =======================================
void block_getPSL(Net *net, Cap *qcap, Cap *tcap, char *query, char *target, FILE *fileHandle){
   /*
    *print PSL for 2 aligned segments (applied for tangle cases)
    */
   char *pslDesc;
   FILE *err = fopen("dump", "w");
   char sep = '\t';
   char lastStep = '\n';
   int blockSpace = 1; //default number of blocks 
   char strand[3];
   //char *netName = netMisc_nameToString(net_getName(net));
   //struct psl * psl = pslNew(netName, 0, 0, 0, target, 0, 0, 0, strand, blockSpace, 0);
   struct psl * psl = pslNew(query, 0, 0, 0, target, 0, 0, 0, strand, blockSpace, 0);
   pslInitialize(psl);
   
   addPSLSizes(psl, qcap, tcap);
   addPSLStarts(psl, qcap, tcap);
   addPSLStrand(psl, qcap, tcap);
   addBlockToPSL(psl, qcap, tcap);//add to qstarts, tstarts   
   addPSLEnds(psl);
   if(psl->blockCount >= 1){
      fixNegStrandPsl(psl);
      if(!isPlus(psl->strand[1])){
         pslRc(psl);
      }
      if(pslCheck(pslDesc, err, psl) == 0){
         pslOutput(psl, fileHandle, sep, lastStep);
      }
   }
   pslFree(&psl);
   return;
}

bool block_inRange(Block *block, char *name, int s, int e){
   bool inrange = false;
   Segment *segment;
   int bstart, bend, temp;
   Block_InstanceIterator *segmentIterator = block_getInstanceIterator(block);
   while((segment = block_getNext(segmentIterator)) != NULL){
      Sequence *sequence = segment_getSequence(segment);
      if(sequence == NULL){continue;}
      char *sequenceHeader = formatSequenceHeader(sequence);
      if(strcmp(sequenceHeader, name) == 0){
         bstart = cap_getCoordinate(segment_get5Cap(segment));
         bend = cap_getCoordinate(segment_get3Cap(segment));
         if(bstart > bend){
            temp = bend;
            bend = bstart;
            bstart = temp;
         }
         if(bstart < e && s < bend){
	    inrange = true;
	    break;
	 }
      }   
   }
   block_destructInstanceIterator(segmentIterator);
   return inrange;
}
void getPSLTangle(Net *net, FILE *fileHandle, char *query, char *target, int s, int e){
   End *end;
   Cap *cap;
   Cap *qcap;
   Cap *tcap;
   int i, j;
   struct Stubs *stubs = setStubs();
   Net_EndIterator *endIterator = net_getEndIterator(net);
   while((end= net_getNextEnd(endIterator)) != NULL){
      Block *block = end_getBlock(end);
     if(block == NULL){ continue; }
     if(block_inRange(block, query, s, e)){
      if(end_isBlockEnd(end) && !block_getChain(block)){//blocks that don't belong to any chain
         End_InstanceIterator *capIterator = end_getInstanceIterator(end);
         while((cap=end_getNext(capIterator)) != NULL){ 
            Sequence *sequence = cap_getSequence(cap);
            if(sequence == NULL){continue;}
            char *sequenceHeader = formatSequenceHeader(sequence);
            if(strcmp(sequenceHeader, query) == 0){
               if(cap_getSide(cap)){//5' cap
                  addStub(&(stubs->qstubs), cap, &(stubs->qnum));
               }
            }else if(strcmp(sequenceHeader, target) == 0){
               if(cap_getSide(cap)){//5' cap
                  addStub(&(stubs->tstubs), cap, &(stubs->tnum));
               }
            }
            //free(sequenceHeader);
         }
         if(stubs->qnum > 0 && stubs->tnum > 0){
            for(i=0; i< stubs->qnum; i++){
	       qcap = *(stubs->qstubs + i);
               for(j=0; j< stubs->tnum; j++){
	          tcap = *(stubs->tstubs + j);
                  block_getPSL(net, qcap, tcap,query, target, fileHandle);
	       }
	    }
	 }
	 if(stubs->qnum > 0 || stubs->tnum > 0){
	    resetStubs(&stubs);
         }
         end_destructInstanceIterator(capIterator);
      }
     }
   }
   net_destructEndIterator(endIterator);

   Net_GroupIterator *groupIterator = net_getGroupIterator(net);
   Group *group;
   while((group = net_getNextGroup(groupIterator)) != NULL) {
      Net *nestedNet = group_getNestedNet(group);
      if(nestedNet != NULL) {
         getPSLTangle(nestedNet, fileHandle, query, target, s, e);
      }
   }
   net_destructGroupIterator(groupIterator);
}

//============================ GETTING ALL THE PSLs =========================
void getPSLs(Net *net, FILE *fileHandle, char *query, char *target, int *starts, int *ends, int size, bool tangle, bool exhaust) {
   /*
    *Print to output file PSLs of net and all (nested) nets in the lower levels
    */
   char *netName = netMisc_nameToString(net_getName(net));
   fprintf(stderr, "\nNET: %s\n", netName);
   int i, start, end;
   Cap *qstartCap = net_getThreadStart(net, query);
   if(size == 0){//No reference - set limit to the whole sequence - Won't do exhaust!
      start = 2;
      Sequence *qseq = net_getSequenceByName(net, query);
      end = sequence_getLength(qseq) + 2;
      fprintf(stderr, "Getting psl for range: <%d -  %d>\n", start, end);
      getPSLNet(fileHandle, query, target, start, end, qstartCap, false);
      if(tangle){
         getPSLTangle(net, fileHandle, query, target, start, end);
      }
   }else{
      for(i=0; i< size; i++){
         fprintf(stderr, "Getting psl for range: <%d -  %d>\n", *(starts +i), *(ends +i));
         getPSLNet(fileHandle, query, target, *(starts + i), *(ends + i), qstartCap, exhaust);
         if(tangle){
            getPSLTangle(net, fileHandle, query, target, *(starts + i), *(ends + i));
         }
      }
   }
   return;
}
int getRanges(struct psl *refpsl, int offset, int **starts, int **ends){
   int i, start, end;
   int size = 0;
   while(refpsl != NULL){
      start = refpsl->tStart - offset +2;
      end = refpsl->tEnd - offset +2;
      if(size == 0){
         (*starts) = AllocN(int, 1);
         (*ends) = AllocN(int, 1);
         *(*starts) = start;
         *(*ends) = end;
         size++;
      }else{
         for(i=0; i < size; i++){
            //overlap with existing limits, merge (not exhaust merge)
	    if(start < *((*ends) + i) && *((*starts) +i) < end){
               if(start < *((*starts) + i)){
                  *((*starts) + i) = start;
               }
               if(end > *((*ends) + i)){
                  *((*ends) + i) = end;
               }
	       break;
	    }
         }
	 if(i == size){//didnt overlap with existing limit --> new limit - add to the limits array
	    size++;
            (*starts) = needMoreMem((*starts), (size -1)*sizeof(int), size*sizeof(int));
            (*ends) = needMoreMem((*ends), (size -1)*sizeof(int), size*sizeof(int));
            *((*starts) + size - 1) = start;
            *((*ends) + size - 1) = end;
	 }
      }
      refpsl = refpsl->next;
   }
   return size;
}

void getAllPSLs(Net *net, FILE *fileHandle, char *query, char *target, struct psl *refpsl, int offset, bool tangle, bool exhaust) {
   char **qseqs;
   char **tseqs;
   char *qName;
   char *tName;
   int q, t;
   net = group_getNestedNet(net_getFirstGroup(net));
   //Look for all query and target sequences with name start with 'query' and 'target'
   int qnum = getSequences(net, &qseqs, query);
   int tnum = getSequences(net, &tseqs, target);

   //Find the set of limits from refpslList
   int *starts = NULL; 
   int *ends = NULL;
   int size = getRanges(refpsl,offset, &starts, &ends);
   for(q = 0; q < qnum; q++){
      qName = *(qseqs + q);
      fprintf(stderr, "Current Query: %s\n", qName);
      for(t = 0; t < tnum; t++){
         tName = *(tseqs + t);
         fprintf(stderr, "\tCurrent Target: %s\n", tName);
         getPSLs(net, fileHandle, qName, tName, starts, ends, size, tangle, exhaust);
      }
   }
   if(size > 0){
      freeMem(starts);
      freeMem(ends);
   }
   if(qnum > 0){
      freeMem(qseqs);
   }
   if(tnum > 0){
      freeMem(tseqs);
   }
   return;
}



void usage() {
   fprintf(stderr, "cactus_pslGenerator, version 0.2\n");
   fprintf(stderr, "-a --logLevel : Set the log level\n");
   fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
   //fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
   fprintf(stderr, "-e --outputFile : The file to write the PSLs in.\n");
   fprintf(stderr, "-q --query : Name of the query sequence.\n");
   fprintf(stderr, "-t --target : Name of the target sequence.\n");
   fprintf(stderr, "-r --ref : The file that has refseq psls with 'query' as the target.\n");
   fprintf(stderr, "-o --offset : Where query start on the genome sequence.\n");
   fprintf(stderr, "-g --tangle : if specified, will include the tangle groups\n");
   fprintf(stderr, "-x --exhaust : if specified, will exhaustly return all possible pairwise alignments.\n");
   fprintf(stderr, "Very slow - do not use for large regions.");
   fprintf(stderr, "If not specified, return the best alignment. Note, if no ref is specified, then will not do exhaust\n");
   fprintf(stderr, "-h --help : Print this help screen\n");
}

//============================== MAIN =========================================
int main(int argc, char *argv[]) {
   NetDisk *netDisk;
   Net *net;

   /*
    * Arguments/options
    */
   char * logLevelString = NULL;
   char * netDiskName = NULL;
   char * netName = NULL;
   char * outputFile = NULL;
   char * query = NULL;
   char * target = NULL;
   char * ref = NULL;
   int offset = 0;
   bool tangle = false;
   bool exhaust = false;

   ///////////////////////////////////////////////////////////////////////////
   // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
   ///////////////////////////////////////////////////////////////////////////

   while(1) {
      static struct option long_options[] = {
         //{ "netName", required_argument, 0, 'd' },
         { "exhaust", no_argument, 0, 'x' },
         { "tangle", no_argument, 0, 'g' },
         { "offset", required_argument, 0, 'o' },
         { "ref", required_argument, 0, 'r' },
         { "query", required_argument, 0, 'q' },
         { "target", required_argument, 0, 't' },
         { "logLevel", required_argument, 0, 'a' },
         { "netDisk", required_argument, 0, 'c' },
         { "outputFile", required_argument, 0, 'e' },
         { "help", no_argument, 0, 'h' },
         { 0, 0, 0, 0 }
      };

      int option_index = 0;

      int key = getopt_long(argc, argv, "o:r:q:t:a:c:d:e:fxgh", long_options, &option_index);

      if(key == -1) {
         break;
      }

      switch(key) {
         case 'a':
            logLevelString = stringCopy(optarg);
            break;
         case 'c':
            netDiskName = stringCopy(optarg);
            break;
         /*case 'd':
            netName = stringCopy(optarg);
            break;*/
         case 'e':
            outputFile = stringCopy(optarg);
            break;
         case 'q':
            query = stringCopy(optarg);
            break;
         case 't':
            target = stringCopy(optarg);
            break;
         case 'r':
            ref = stringCopy(optarg);
            break;
         case 'o':
            sscanf(stringCopy(optarg), "%d", &offset);
            break;
         case 'g':
            tangle = true;
            break;
         case 'x':
            exhaust = true;
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

   assert(netDiskName != NULL);
   //assert(netName != NULL);
   assert(outputFile != NULL);
   assert(query != NULL);
   assert(target != NULL);

   //////////////////////////////////////////////
   //Set up logging
   //////////////////////////////////////////////

   if(logLevelString != NULL && strcmp(logLevelString, "INFO") == 0) {
      setLogLevel(LOGGING_INFO);
   }
   if(logLevelString != NULL && strcmp(logLevelString, "DEBUG") == 0) {
      setLogLevel(LOGGING_DEBUG);
   }

   //////////////////////////////////////////////
   //Log (some of) the inputs
   //////////////////////////////////////////////

   logInfo("Net disk name : %s\n", netDiskName);
   logInfo("Net name : %s\n", netName);
   logInfo("Output PSL file : %s\n", outputFile);
   logInfo("Query: %s\n", query);
   logInfo("Target: %s\n", target);

   //////////////////////////////////////////////
   //Load the database
   //////////////////////////////////////////////

   netDisk = netDisk_construct(netDiskName);
   logInfo("Set up the net disk\n");

   ///////////////////////////////////////////////////////////////////////////
   // Parse the basic reconstruction problem
   ///////////////////////////////////////////////////////////////////////////
   netName = stringCopy("0");
   net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
   logInfo("Parsed the top level net of the cactus tree to check\n");

   ///////////////////////////////////////////////////////////////////////////
   // Recursive check the nets.
   ///////////////////////////////////////////////////////////////////////////

   int64_t startTime = time(NULL);
   FILE *fileHandle = fopen(outputFile, "w");
   makePSLHeader(net, fileHandle);
   struct psl *refpsl = NULL;
   if(ref != NULL){
      refpsl = pslLoadAll(ref);
   }
   getAllPSLs(net, fileHandle, query, target, refpsl, offset, tangle, exhaust);
   fclose(fileHandle);
   logInfo("Got the psls in %i seconds/\n", time(NULL) - startTime);

   ///////////////////////////////////////////////////////////////////////////
   // Clean up.
   ///////////////////////////////////////////////////////////////////////////

   netDisk_destruct(netDisk);

   return 0;
}
