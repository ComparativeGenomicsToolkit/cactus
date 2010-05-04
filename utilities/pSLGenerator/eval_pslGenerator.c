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
 * The script outputs a psl file containing psl records for all nets.
 */
void getPSL(Cap **qcap, Cap **tcap, FILE *fileHandle, char *query, char *target, char *netName, char *prevStrands);

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

bool isStubCap(Cap *cap){
   return (end_isStubEnd(cap_getEnd(cap))) ? true : false;
}

bool isPlus(char strand){
   if(strand == '+'){
      return true;
   }else{
      return false;
   }
}

void moveCapToNextBlock(Cap **cap){
   /*Move cap to the next block that its segment is adjacency to*/
   //if(cap_getSegment(*cap) != NULL){
      if(cap_getSide(*cap)){//5' end
         *cap = cap_getAdjacency(segment_get3Cap(cap_getSegment(*cap)));
      }else{
         *cap = cap_getAdjacency(segment_get5Cap(cap_getSegment(*cap)));
      }
   //}
}

void moveCapToPrevBlock(Cap **cap){
   /*Move cap to the previous block that its segment is adjacency to*/
   *cap = cap_getAdjacency(*cap);
   //if(cap_getSegment(*cap) != NULL){
      if(cap_getSide(*cap)){//5' end
         *cap = segment_get3Cap(cap_getSegment(*cap));
      }else{
         *cap = segment_get5Cap(cap_getSegment(*cap));
      }
   //}
}

void moveCapToStubEnd(Cap **cap){
   /*move cap to the other end of the sequence (stub)*/
   if(!isStubCap(*cap)){
      while(!isStubCap(*cap)){
         moveCapToNextBlock(cap);
      }
   }
}

bool compareStrands(bool qstrand, bool tstrand, char *prevStrands){
   /*compare qstrand & tstrand with prevStrands
    *++ and -- are equivalent; +- and -+ are equivalent
    */
   char *strands; 
   char *prevstrands;
   if(prevStrands == NULL){ 
      return false; 
   }
   if((!qstrand && !tstrand) ||(qstrand && tstrand)){  
      strands = "++";
   }else{
      strands = "-+";
   }
   bool qPrevStrand = isPlus(prevStrands[0]);
   bool tPrevStrand = isPlus(prevStrands[1]);
   if((!qPrevStrand && !tPrevStrand) ||(qPrevStrand && tPrevStrand)){  
      prevstrands = "++";
   }else{
      prevstrands = "-+";
   }
   if(strcmp(strands, prevstrands) == 0){
      return true;
   }else{
      return false;
   }
}

void fixNegStrandPsl(struct psl *psl){
   //convert the negative strand psl into correct coordinates
   if(!isPlus(psl->strand[0])){//query is on negative strand
      int32_t temp = psl->qSize - psl->qStart;
      psl->qStart = psl->qSize - psl->qEnd;
      psl->qEnd = temp;
   }
   if(!isPlus(psl->strand[1])){//target is on negative strand
      int32_t temp = psl->tSize - psl->tStart;
      psl->tStart = psl->tSize - psl->tEnd;
      psl->tEnd = temp;
   }
   return; 
}

int32_t getSegmentStart(Segment *segment, bool isNotReverse){
   /*Get the start coordinate of input segment
    *cactus starts from 1 while psl starts from 0. Cactus' segment-start
    *also 1 off from the real start because of the cap. Therefore minus 2
    */
   if(isNotReverse){	
      return segment_getStrand(segment) ? segment_getStart(segment)-2 
                  : segment_getStart(segment) - segment_getLength(segment)-1;
   }else{
      Sequence *sequence = segment_getSequence(segment);
      //convert to cactus coordinate relative to the negative strand
      int32_t start = sequence_getLength(sequence) - segment_getStart(segment) +1;
      return segment_getStrand(segment) ? start - segment_getLength(segment) +1
                                        : start;
   }
}

int net_getStubsBySequenceName(Net *net, Cap ***caps, char *name, bool capSide){
   /*Get caps (end-instance) by the sequence it is from, and the specified capSide (true if 5' and false if 3')*/
   //return number of caps found;
   Cap *cap = NULL;
   int numCaps = 0;
   Net_CapIterator *capIterator = net_getCapIterator(net);
   while((cap= net_getNextCap(capIterator)) != NULL){
      if(isStubCap(cap)){//dead end or inherited end
         Sequence *sequence = cap_getSequence(cap);
         char *sequenceHeader = formatSequenceHeader(sequence);
         if(strcmp(sequenceHeader, name) == 0){
            free(sequenceHeader);
            if((cap_getSide(cap) && capSide) ||(!(cap_getSide(cap)) && !capSide)){
               //break;
               numCaps++;
	       if(numCaps ==1){
                  *caps = malloc(sizeof(Cap *));
		  **caps = cap;
               }else{
	          *caps = (Cap **)realloc(*caps, numCaps*sizeof(Cap *));
	          *(*caps + numCaps-1) = cap;
	       }
            }
         }
      }
   }
   net_destructCapIterator(capIterator);
   return numCaps;
}

int get3Caps(Net *net, Cap ***caps, char *name){
   //get cap from 3' deadEnd or inheritedEnd if exist
   //return number of caps found
   int numCaps = net_getStubsBySequenceName(net, caps, name, false);
   return numCaps;
}

int get5Caps(Net *net, Cap ***caps, char *name){
   //get cap from 5' deadEnd or inheritedEnd if exist
   //return number of caps found
   int numCaps = net_getStubsBySequenceName(net, caps, name, true);
   return numCaps;
}

void addPSLSizes(struct psl *psl, Cap *qcap, Cap *tcap){
   Sequence *qsequence = cap_getSequence(qcap);
   Sequence *tsequence = cap_getSequence(tcap);
   psl->qSize = sequence_getLength(qsequence);
   psl->tSize = sequence_getLength(tsequence);
}

void addPSLStarts(struct psl *psl, Cap *qcap, Cap *tcap){
   //Get threads' starts relatively to the forward strands
   int32_t qsize = sequence_getLength(cap_getSequence(qcap));
   int32_t tsize = sequence_getLength(cap_getSequence(tcap));
   int32_t qstart = cap_getCoordinate(qcap);
   int32_t tstart = cap_getCoordinate(tcap);
   psl->qStart = cap_getStrand(qcap)? qstart-1 : qsize - qstart +2;
   psl->tStart = cap_getStrand(tcap)? tstart-1 : tsize - tstart +2;
}

void addPSLStrand(struct psl *psl, Cap *qcap, Cap *tcap){
   psl->strand[0] = cap_getStrand(qcap) ? '+' : '-';
   psl->strand[1] = cap_getStrand(tcap) ? '+' : '-';
   //psl->strand[2] = NULL;
}

void addPSLEnds(struct psl *psl, Cap *qcap, Cap *tcap){
   int32_t qsize = sequence_getLength(cap_getSequence(qcap));
   int32_t tsize = sequence_getLength(cap_getSequence(tcap));
   int32_t qend = cap_getCoordinate(qcap);
   int32_t tend = cap_getCoordinate(tcap);
   psl->qEnd = cap_getStrand(qcap)? qend-2 : qsize - qend +1;
   psl->tEnd = cap_getStrand(tcap)? tend-2 : tsize - tend +1;
}

void addPSLInserts(struct psl *psl, int32_t currIndex){
   int32_t qdiff;
   int32_t tdiff;
   if(currIndex == 0){
      qdiff = *(psl->qStarts) - psl->qStart;
      tdiff = *(psl->tStarts) - psl->tStart;
   }else{
      qdiff = *(psl->qStarts + currIndex) - *(psl->qStarts + currIndex -1) - *(psl->blockSizes + currIndex -1);
      tdiff = *(psl->tStarts + currIndex) - *(psl->tStarts + currIndex -1) - *(psl->blockSizes + currIndex -1);
   }
   if(qdiff > 0){
      psl->qNumInsert ++;
      psl->qBaseInsert += qdiff;
   }
   if(tdiff > 0){
      psl->tNumInsert ++;
      psl->tBaseInsert += tdiff;
   }
}

void addPSLLastInserts(struct psl *psl){
   int32_t currIndex = psl->blockCount -1;
   int32_t qendGap = psl->qEnd - *(psl->qStarts +currIndex) - *(psl->blockSizes +currIndex);
   int32_t tendGap = psl->tEnd - *(psl->tStarts +currIndex) - *(psl->blockSizes +currIndex);
   if(qendGap >0){
      psl->qNumInsert++;
      psl->qBaseInsert += qendGap;
   }
   if(tendGap >0){
      psl->tNumInsert++;
      psl->tBaseInsert += tendGap;
   }
}

void addBlockToPSL(struct psl *psl, Segment *qsegment, Segment *tsegment){
   /*
    */
   psl->blockCount++;
   int32_t length = segment_getLength(qsegment);
   int32_t currIndex = psl->blockCount -1;
   //blockSizes
   psl->blockSizes = (unsigned *)realloc(psl->blockSizes, psl->blockCount*sizeof(unsigned));
   *(psl->blockSizes + currIndex) = length;
   //qStarts
   psl->qStarts = (unsigned *)realloc(psl->qStarts, psl->blockCount*sizeof(unsigned));
   *(psl->qStarts + currIndex) = isPlus(psl->strand[0]) ? getSegmentStart(qsegment, true) : getSegmentStart(qsegment, false);
   //tStarts
   psl->tStarts = (unsigned *)realloc(psl->tStarts, psl->blockCount*sizeof(unsigned));
   *(psl->tStarts + currIndex) = isPlus(psl->strand[1]) ? getSegmentStart(tsegment, true) : getSegmentStart(tsegment, false);
   psl->match += length;
   //Inserts
   addPSLInserts(psl, currIndex);
}

int32_t block_getSegmentNumber(Block *block, char *name){
   //get number of segments in input block that come from the sequence 'name'
   Block_InstanceIterator *segmentIterator = block_getInstanceIterator(block);
   Segment * segment;
   int32_t count = 0;
   while((segment = block_getNext(segmentIterator)) != NULL){
      Sequence *sequence = segment_getSequence(segment);
                char *sequenceHeader = formatSequenceHeader(sequence);
                if(strcmp(sequenceHeader, name) == 0){
                   count++;
      }
                free(sequenceHeader);
   }
   block_destructInstanceIterator(segmentIterator);
   return count;
}

Segment *goToNextCommonBlock (Cap **cap, char *name){
   Block *block;
   bool commonBlock = false;
   Segment *segment;
   while(!commonBlock && !isStubCap(*cap)){
      segment = cap_getSegment(*cap);
      block = segment_getBlock(segment);
      //if(block_getChain(block) != NULL){//only look at chain-blocks
      Block_InstanceIterator *segmentIterator = block_getInstanceIterator(block);
      Segment *segment2;
      while((segment2 = block_getNext(segmentIterator)) != NULL){
         Sequence *sequence = segment_getSequence(segment2);
         char *sequenceHeader = formatSequenceHeader(sequence);
         if(strcmp(sequenceHeader, name) == 0){
            free(sequenceHeader);
            commonBlock = true;
            break;
         }
      }
      block_destructInstanceIterator(segmentIterator);
      //}
      moveCapToNextBlock(cap);
   }
   //free(block);
   return segment;
}

void goToNextCommonBlock2(Segment **qsegment, Segment **tsegment, Cap **qcap, Cap **tcap, char *query, char *target){ 
   Cap * qtempcap = *qcap;
   Cap * ttempcap = *tcap;
   Segment *qtempseg = *qsegment;
   Segment *ttempseg = *tsegment;
   int32_t qcount = block_getSegmentNumber(segment_getBlock(*qsegment), query);
   int32_t tcount = block_getSegmentNumber(segment_getBlock(*tsegment), target);
   if(qcount > tcount){
      while((segment_getBlock(qtempseg) != segment_getBlock(*tsegment)) && qtempseg != NULL && !isStubCap(qtempcap)){
         qtempseg = goToNextCommonBlock(&qtempcap, target);
   }
      *qcap = qtempcap;
      *qsegment = qtempseg;
   }else if(qcount < tcount){
      while((segment_getBlock(ttempseg) != segment_getBlock(*qsegment)) && ttempseg != NULL && !isStubCap(ttempcap)){
         ttempseg = goToNextCommonBlock(&ttempcap, target);
      }
      *tcap = ttempcap;
      *tsegment = ttempseg;
   }else{
      while( ((segment_getBlock(qtempseg) != segment_getBlock(*tsegment)) && qtempseg != NULL && !isStubCap(qtempcap)) &&
             ((segment_getBlock(ttempseg) != segment_getBlock(*qsegment)) && ttempseg != NULL && !isStubCap(ttempcap)) ){
         qtempseg = goToNextCommonBlock(&qtempcap, target);
         ttempseg = goToNextCommonBlock(&ttempcap, query);
      }
      if(segment_getBlock(qtempseg) == segment_getBlock(*tsegment) || ttempseg == NULL){
               *qcap = qtempcap;
               *qsegment = qtempseg;
      }else{
               *tcap = ttempcap;
               *tsegment = ttempseg;
      }
   }
   return; 
}

void addPSLBlocks(struct psl *psl, Cap **qcap, Cap **tcap, char *query, char *target, char *prevStrands, FILE *fileHandle){
   Segment *qsegment;
   Segment *tsegment;
   if(prevStrands == NULL){
      *qcap = cap_getAdjacency(*qcap);
      *tcap = cap_getAdjacency(*tcap);
   }
   while( !isStubCap(*qcap) && !isStubCap(*tcap) && //EXIT 'while' only  when q or t is stub
          !compareStrands(cap_getStrand(*qcap), cap_getStrand(*tcap), prevStrands) ){//or when switches back to previous strands
      qsegment = goToNextCommonBlock(qcap, target);
      tsegment = goToNextCommonBlock(tcap, query);
      if(qsegment == NULL || tsegment == NULL){
         return;
      }
      if(segment_getBlock(qsegment) != segment_getBlock(tsegment)){
         if(isStubCap(*qcap) || isStubCap(*tcap)){
            //tangle net.. one sequence reach the end already and not come back...
            break;
         }
         goToNextCommonBlock2(&qsegment, &tsegment, qcap, tcap, query, target);
      }
      char pslStrands[2];
      pslStrands[0] = psl->strand[0];
      pslStrands[1] = psl->strand[1];
      if(!compareStrands(segment_getStrand(qsegment),segment_getStrand(tsegment),pslStrands)){//switch to new strands
         moveCapToPrevBlock(qcap);
         moveCapToPrevBlock(tcap);
         getPSL(qcap, tcap, fileHandle, query, target, psl->qName, pslStrands);
      }
      if(segment_getBlock(qsegment) == segment_getBlock(tsegment)){
         addBlockToPSL(psl, qsegment, tsegment);//add to qstarts, tstarts   
      }else{
         //printf("Potentially early-terminated net ... Tangled!!!\n");
      }
   }
   //if exited 'while' NOT because need to switch back to previous strands:
   if(isStubCap(*qcap) || isStubCap(*tcap)){
      moveCapToStubEnd(qcap);
      moveCapToStubEnd(tcap);
   }
   addPSLEnds(psl, *qcap, *tcap);
   addPSLLastInserts(psl);
   return;
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

//==============MAP LOCATION OF PSL BLOCKS TO THE REFERENCE SEQUENCES================
int32_t mapCoor(int32_t refsize, int32_t start, int32_t coor, bool strand){
   /*Map 'coor' to the coordinates of the reference sequence
    *refsize: size of reference sequence, start: location of current seq on the ref-seq
    *strand: + if current seq is on the forward strand of the ref-seq, - otherwise
    */
   return strand? coor + start : refsize - (coor + start);   
}
//Given that all of cactus' input sequences are on the forward strand
void mapPslCoor(struct psl *psl, int32_t qrefsize, int32_t qstart, int32_t qsize, int32_t trefsize, int32_t tstart, int32_t tsize){
   int i;
   psl->qSize = qrefsize;
   psl->tSize = trefsize;
   psl->qStart = mapCoor(qrefsize, qstart, psl->qStart, true);
   psl->qEnd = mapCoor(qrefsize, qstart, psl->qEnd, true);
   psl->tStart = mapCoor(trefsize, tstart, psl->tStart, true);
   psl->tEnd = mapCoor(trefsize, tstart, psl->tEnd, true);
   if(!isPlus(psl->strand[0])){
      qstart = qrefsize - (qstart + qsize);
   }
   if(!isPlus(psl->strand[1])){
      tstart = trefsize - (tstart + tsize);
   }
   for(i=0; i< psl->blockCount; i++){
      *(psl->qStarts + i) = mapCoor(qrefsize, qstart, *(psl->qStarts +i), true);
      *(psl->tStarts + i) = mapCoor(trefsize, tstart, *(psl->tStarts +i), true);
   }
   return;
}
//=============END MAPPING=========

void getPSL(Cap **qcap, Cap **tcap, FILE *fileHandle, char *query, char *target, char *netName, char *prevStrands){
   /*
    *Get PSL that starts with qcap & tcap, current level
    */
   //char *psldesc;
   //FILE *out = fopen("err", "a");
   int blockSpace = 32; //default number of blocks 
   char strand[3];
   //struct psl * psl = pslNew(query, 0, 0, 0, target, 0, 0, 0, strand, blockSpace, 0);
   struct psl * psl = pslNew("chr6", 0, 0, 0, "chr6", 0, 0, 0, strand, blockSpace, 0);
   //struct psl * psl = pslNew(netName, 0, 0, 0, target, 0, 0, 0, strand, blockSpace, 0);
   pslInitialize(psl);
   addPSLSizes(psl, *qcap, *tcap);
   addPSLStarts(psl, *qcap, *tcap);
   addPSLStrand(psl, *qcap, *tcap);
   addPSLBlocks(psl, qcap, tcap, query, target, prevStrands,fileHandle);
   if(psl->blockCount >= 1){
      char sep = '\t';
      char lastStep = '\n';
      fixNegStrandPsl(psl);
      //mapPslCoor(psl, 173908612, 30906992, 7482, 171115067, 30227334, 7439);
      //pslCheck(psldesc, out, psl);
      if(!isPlus(psl->strand[1])){
         pslRc(psl);
      }
      pslOutput(psl, fileHandle, sep, lastStep);
   }
   pslFree(&psl);
   return;
}
void getPSLNet(Net *net, FILE *fileHandle, char *query, char *target){
   /*
    *Get PSLs for net, current level
    */
   int i = 0;
   int j = 0;
   char *netName = netMisc_nameToString(net_getName(net));
   //Get the start of query & target segments at current level
   //Start from 3' stubs (which then go to the 5' starts of the interested sequence)
   Cap **qcaps = NULL;
   Cap **tcaps = NULL;
   Cap *qcap ;
   Cap *tcap;
   int qnumCaps = get3Caps(net, &qcaps, query);
   int tnumCaps = get3Caps(net, &tcaps, target);
   //while (i < tnumCaps && i < qnumCaps){
   while (i < qnumCaps){
      while (j < tnumCaps){
         qcap = *(qcaps + i);
         tcap = *(tcaps + j);
         getPSL(&qcap, &tcap, fileHandle, query, target, netName, NULL);
         j++;
      }
      i++;
   }
   free(qcaps);
   free(tcaps);
}
void getPSLs(Net *net, FILE *fileHandle, char *query, char *target) {
   /*
    *Print to output file PSLs of net and all (nested) nets in the lower levels
    */
   //Make PSL record for net's blocks
   getPSLNet(net, fileHandle, query, target);

   //Call child nets recursively.
   Net_GroupIterator *groupIterator = net_getGroupIterator(net);
   Group *group;
   while((group = net_getNextGroup(groupIterator)) != NULL) {
      Net *nestedNet = group_getNestedNet(group);
      if(nestedNet != NULL) {
         getPSLs(group_getNestedNet(group), fileHandle, query, target); //recursive call.
      }
   }
   net_destructGroupIterator(groupIterator);
}
void makePSLHeader(Net *net, FILE *fileHandle) {
   pslWriteHead(fileHandle);
   //fprintf(fileHandle, "track name=\"hcrm:Rat-hlaB\" description=\"hlaB - human,chimp,rat,mouse: RAT\" color=193,72,52 visibility=pack\n");
   //fprintf(fileHandle, "track name=\"hc:CHIMP\" description=\"HCG27-MICA - human,chimp: CHIMP\" color=193,72,52 visibility=pack\n");
}

void usage() {
   fprintf(stderr, "cactus_pslGenerator, version 0.2\n");
   fprintf(stderr, "-a --logLevel : Set the log level\n");
   fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
   fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
   fprintf(stderr, "-e --outputFile : The file to write the PSLs in.\n");
   fprintf(stderr, "-q --query : Name of the query sequence.\n");
   fprintf(stderr, "-t --target : Name of the target sequence.\n");
   fprintf(stderr, "-h --help : Print this help screen\n");
}

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

   ///////////////////////////////////////////////////////////////////////////
   // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
   ///////////////////////////////////////////////////////////////////////////

   while(1) {
      static struct option long_options[] = {
         { "query", required_argument, 0, 'q' },
         { "target", required_argument, 0, 't' },
         { "logLevel", required_argument, 0, 'a' },
         { "netDisk", required_argument, 0, 'c' },
         { "netName", required_argument, 0, 'd' },
         { "outputFile", required_argument, 0, 'e' },
         { "help", no_argument, 0, 'h' },
         { 0, 0, 0, 0 }
      };

      int option_index = 0;

      int key = getopt_long(argc, argv, "a:c:d:e:q:t:fh", long_options, &option_index);

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
         case 'd':
            netName = stringCopy(optarg);
            break;
         case 'e':
            outputFile = stringCopy(optarg);
            break;
         case 'q':
            query = stringCopy(optarg);
            break;
         case 't':
            target = stringCopy(optarg);
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
   assert(netName != NULL);
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

   net = netDisk_getNet(netDisk, netMisc_stringToName(netName));
   logInfo("Parsed the top level net of the cactus tree to check\n");

   ///////////////////////////////////////////////////////////////////////////
   // Recursive check the nets.
   ///////////////////////////////////////////////////////////////////////////

   int64_t startTime = time(NULL);
   FILE *fileHandle = fopen(outputFile, "w");
   makePSLHeader(net, fileHandle);
   getPSLs(net, fileHandle, query, target);
   fclose(fileHandle);
   logInfo("Got the psls in %i seconds/\n", time(NULL) - startTime);

   ///////////////////////////////////////////////////////////////////////////
   // Clean up.
   ///////////////////////////////////////////////////////////////////////////

   netDisk_destruct(netDisk);

   return 0;
}
