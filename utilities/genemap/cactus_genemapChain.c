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

#define LINEMAXSIZE 3000

/*
 * nknguyen@soe.ucsc.edu
 * 05/26/2010
 * Edit: Sep 14 2010 - transcripts or CDS inputed in BED format instead of PSL
 * Edit: Nov 09 2010 - clean up unnecessary functions
 ************************
 * This script sees how genes (or transcripts) map to cactus graph.
 * We want to compare the gene structure with the chain structure of cactus
 * Input a list of reference genes, and returns an XML tree how these genes
 * mapped with the chain/block structure
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


//============================= UTILS FUNCTIONS ========================

/*void moveCapToNextBlock(Cap **cap){
   //Move cap to the next block that its segment is adjacency to
   if(isStubCap(*cap)){
      return;
   }
   Cap *adjCap = cap_getAdjacency(*cap);
   if(cap_getEnd(*cap) == cap_getEnd(adjCap)){//DOUBLE CHECK ... self connected end
      *cap = adjCap;
   }else{
      *cap = cap_getAdjacency(cap_getOtherSegmentCap(*cap));
   }
}*/

void mapBlockToExon(Cap *cap, int level, FILE *fileHandle){
   fprintf(fileHandle, "\t\t\t<block>\n");
   Block *block = end_getBlock(cap_getEnd(cap));
   Chain *chain = block_getChain(block);
   int start = cap_getCoordinate(cap);
   int end = cap_getCoordinate(cap_getOtherSegmentCap(cap)) +1;
   fprintf(fileHandle, "\t\t\t\t<blockName>%s</blockName>\n", cactusMisc_nameToString(block_getName(block)));
   if(chain != NULL){
      fprintf(fileHandle, "\t\t\t\t<chainName>%s</chainName>\n", cactusMisc_nameToString(chain_getName(chain)));
   }else{
      fprintf(fileHandle, "\t\t\t\t<chainName>NA</chainName>\n");
   }
   fprintf(fileHandle, "\t\t\t\t<level>%d</level>\n", level);
   fprintf(fileHandle, "\t\t\t\t<start>%d</start>\n", start);
   fprintf(fileHandle, "\t\t\t\t<end>%d</end>\n", end);
   fprintf(fileHandle, "\t\t\t</block>\n");
   st_logInfo("mapBlockToExon: start: %d, end: %d\n", start, end);
}

int mapGene(Cap *cap, int level, int exon, struct bed *gene, FILE *fileHandle){
   /*
    *Following cactus adjacencies, starting from 'cap', find regions that overlap with 
    *exons of input gene. Report chain relations of these regions with the exons.
    *cap: current cap. Level = chain level. exon = exon number. gene = bed record of gene
    */
   int32_t exonStart, exonEnd;
   if(isStubCap(cap)){
      Group *group = end_getGroup(cap_getEnd(cap));
      Flower *nestedFlower = group_getNestedFlower(group);
      if(nestedFlower != NULL){//recursive call
         Cap *childCap = flower_getCap(nestedFlower, cap_getName(cap));
         assert(childCap != NULL);
         exon = mapGene(childCap, level + 1, exon, gene, fileHandle);
         exonStart = gene->chromStarts->list[exon] + gene->chromStart;
         exonEnd = exonStart + gene->blockSizes->list[exon];
      }
   }

   cap = cap_getAdjacency(cap);
   Cap *nextcap;
   int32_t capCoor;
   exonStart = gene->chromStarts->list[exon] + gene->chromStart;
   exonEnd = exonStart + gene->blockSizes->list[exon];
   Block *block = end_getBlock(cap_getEnd(cap));  
 
   if(block == NULL){
      moveCapToNextBlock(&cap);
   }
   while(!isStubCap(cap) && exon < gene->blockCount){
      End *cend = cap_getEnd(cap);
      capCoor = cap_getCoordinate(cap);//Cap coordinate is always the coordinate on + strand
      nextcap = cap_getAdjacency(cap_getOtherSegmentCap(cap));
      st_logInfo("capCoor: %d, nextCap: %d, eStart: %d, eEnd: %d. Exon: %d\n", 
                  capCoor, cap_getCoordinate(nextcap), exonStart, exonEnd, exon);

      //keep moving if nextBlock Start is still upstream of current exon
      if(cap_getCoordinate(nextcap) <= exonStart){
         moveCapToNextBlock(&cap);
         st_logInfo("Still upstream, nextcap <= exonStart. Move to next chainBlock\n");
      }else if(capCoor >= exonEnd){//Done with current exon, move to next
         st_logInfo("Done with current exon, move to next one\n\n");
         fprintf(fileHandle, "\t\t</exon>\n");//end previous exon
         exon++;
         if(exon < gene->blockCount){
            exonStart = gene->chromStarts->list[exon] + gene->chromStart;
            exonEnd = exonStart + gene->blockSizes->list[exon];
            fprintf(fileHandle, "\t\t<exon id=\"%d\" start=\"%d\" end=\"%d\">\n", exon, exonStart, exonEnd);
         }
      }else{//current exon overlaps with current block Or with lower level flower
         Cap *oppcap = cap_getOtherSegmentCap(cap);
         st_logInfo("Current exon overlaps with current block or with lower flower\n");
         if(cap_getCoordinate(oppcap) >= exonStart && exonEnd > capCoor){
            mapBlockToExon(cap, level, fileHandle);
            if(exonEnd <= cap_getCoordinate(oppcap) + 1){
               st_logInfo("Done with current exon, move to next one\n\n");
               fprintf(fileHandle, "\t\t</exon>\n");//end previous exon
               exon++;
	       if(exon < gene->blockCount){
		  exonStart = gene->chromStarts->list[exon] + gene->chromStart;
		  exonEnd = exonStart + gene->blockSizes->list[exon];
		  fprintf(fileHandle, "\t\t<exon id=\"%d\" start=\"%d\" end=\"%d\">\n", exon, exonStart, exonEnd);
	       }
               continue;
            }
         }
         //Traverse lower level flowers if exists
         Group *group = end_getGroup(end_getOtherBlockEnd(cend));
         Flower *nestedFlower = group_getNestedFlower(group);
         if(nestedFlower != NULL){//recursive call
            Cap *childCap = flower_getCap(nestedFlower, cap_getName(cap_getOtherSegmentCap(cap)));
            assert(childCap != NULL);
            exon = mapGene(childCap, level + 1, exon, gene, fileHandle);
            exonStart = gene->chromStarts->list[exon] + gene->chromStart;
            exonEnd = exonStart + gene->blockSizes->list[exon];
         }
         moveCapToNextBlock(&cap);
      }
   }
   return exon;
}

//=================== GET START OF THREAD (3' STUB) ====================================
/*struct List *flower_getThreadStart(Flower *flower, char *name){// from Flower = 0
    //Get 3' end Stub (the start) of the sequence by its name
   struct List *list = constructEmptyList(0, free);
   Cap *cap;
   Flower_CapIterator *capIterator = flower_getCapIterator(flower);
   while((cap= flower_getNextCap(capIterator)) != NULL){
      if(isStubCap(cap) && !cap_getSide(cap)){//3' dead end or inherited end
         Sequence *sequence = cap_getSequence(cap);
         if(sequence == NULL){continue;}
         char *sequenceHeader = formatSequenceHeader(sequence);
         if(strstr(sequenceHeader, name) != NULL){
            listAppend(list, cap);
         }
         free(sequenceHeader);
      }
   }
   flower_destructCapIterator(capIterator);
   return list;
}*/


//============================ MAP ALL GENES =========================
void printOpeningTag(const char *tag, FILE *fileHandle) {
        /*
         * Creates an opening XML tag.
         */
        fprintf(fileHandle, "<%s>", tag);
}

void printClosingTag(const char *tag, FILE *fileHandle) {
        /*
         * Creates a closing XML tag.
         */
        fprintf(fileHandle, "</%s>", tag);
}

void mapGenes(Flower *flower, FILE *fileHandle, struct bed *gene, char *species){
   st_logInfo("Flower %s\n", cactusMisc_nameToString(flower_getName(flower)));
   printOpeningTag("geneMap", fileHandle);
   fprintf(fileHandle, "\n");
   
   int level = 0;//Flower level
   while(gene != NULL){
      //Get the start of the target sequence: 
      st_logInfo("Gene %s:\n", gene->name);
      Cap *startCap;
      struct List *capList = flower_getThreadStarts(flower, species);
      for(int i=0; i < capList->length; i++){
          startCap = capList->list[i];
          st_logInfo("Cap %d, %s\n", i, cactusMisc_nameToString(cap_getName(startCap)));
	  //Traverse cactus and get regions that overlap with exons of the gene, report the involved chains relations
	  fprintf(fileHandle, "\t<gene name=\"%s\" target=\"%s\" start=\"%d\" end=\"%d\" exonCount=\"%d\" strand=\"%c\">\n", 
                                 gene->name, species, gene->chromStart, gene->chromEnd, gene->blockCount, gene->strand[0]);
	  fprintf(fileHandle, "\t\t<exon id=\"0\" start=\"%d\" end=\"%d\">\n", 
                                 gene->chromStart, gene->chromStart + gene->blockSizes->list[0]);
	  
          mapGene(startCap, level, 0, gene, fileHandle);
	  
          fprintf(fileHandle, "\t</gene>\n");
      }
      gene = gene->next;
   }
   printClosingTag("geneMap", fileHandle);
   return;
}

void usage() {
   fprintf(stderr, "cactus_geneMap, version 0.2\n");
   fprintf(stderr, "-a --st_logLevel : Set the st_log level\n");
   fprintf(stderr, "-c --cactusDisk : The location of the flower disk directory\n");
   fprintf(stderr, "-e --outputFile : output file, in xml format.\n");
   fprintf(stderr, "-s --species : Name of species of interest, e.g 'hg18' or 'hg18.chr11.134452384.116118685.50084.1'\n");
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
