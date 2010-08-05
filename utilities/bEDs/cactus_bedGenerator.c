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
 *Prints to output file all segments of the target sequence that are in blocks that contain both query & target
 *For example, if option 'query' is specified as 'hg19', and 'target' is as 'panTro2'. The program will look for
 *all blocks in the input net and all recursive nets that contain segments of both hg19 and panTro2, and print 
 *out all panTro2's (the target) BED records.
 *Each BED record is just the Start and End of one segment.
 *
 *In case option 'ref' is specified:
 *   'ref' is a file that contains all the PSL records of reference transcripts/genes, whose target is the 'query' above.
 *   In our example, ref's target is 'hg19'.
 *   The program will only print BEDs of target segments that are in the same blocks with those query segments that are
 *   overlap with the range found in ref PSLs.
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

void printCap(FILE *fileHandle, Cap *cap){
   char side = cap_getSide(cap) ? '5' : '3';
   char strand = cap_getStrand(cap) ? '+' : '-';
   int coor = cap_getCoordinate(cap);
   Segment *segment = cap_getSegment(cap);
   Cap *otherCap = cap_getSide(cap) ? segment_get3Cap(segment): segment_get5Cap(segment);
   int otherCoor = cap_getCoordinate(otherCap);
   fprintf(fileHandle, "%s", cactusMisc_nameToString(cap_getName(cap)));
   fprintf(fileHandle, " %c%c%d (%d)", side, strand, coor, otherCoor);
}

bool block_inRange(Block *block, char *name, int s, int e){
   bool inrange = false;
   Segment *segment;
   int bstart, bend, temp;
   Block_InstanceIterator *segmentIterator = block_getInstanceIterator(block);
   while((segment = block_getNext(segmentIterator)) != NULL){
      Sequence *sequence = segment_getSequence(segment);
      if(sequence != NULL){
         char *sequenceHeader = formatSequenceHeader(sequence);
         if(strstr(sequenceHeader, name) == sequenceHeader){
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
   }
   block_destructInstanceIterator(segmentIterator);
   return inrange;
}

bool block_hasQueryTarget(Block *block, char *query, char *target){
   //Return TRUE if contains segments of BOTH query and target
   Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
   Segment *segment;
   bool hasQuery = false;
   bool hasTarget = false;
   while((segment = block_getNext(instanceIterator)) != NULL){
      Sequence *sequence = segment_getSequence(segment);
      if(sequence != NULL){ 
         char *sequenceHeader = formatSequenceHeader(sequence);
         if(strstr(sequenceHeader, query) == sequenceHeader){
            hasQuery = true;
         }else if(strstr(sequenceHeader, target) == sequenceHeader){
            hasTarget = true;
         }
         free(sequenceHeader);
      }
   }
   block_destructInstanceIterator(instanceIterator);
   if(hasTarget && hasQuery){
      return true;
   }else{
      return false;
   }
}

void getBEDBlock(Block *block, FILE *fileHandle, char *target, int tstart) {
	/*
	 * Outputs target BED records of the block to the given file handle.
	 */
	Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
	Segment *segment;
	while((segment = block_getNext(instanceIterator)) != NULL) {
		Sequence *sequence = segment_getSequence(segment);
		if(sequence != NULL) {
			char *sequenceHeader = formatSequenceHeader(sequence);
         		if(strstr(sequenceHeader, target) == sequenceHeader){
                                Cap *cap5 = segment_get5Cap(segment);
                                Cap *cap3 = segment_get3Cap(segment);
				int32_t coor5 = cap_getCoordinate(cap5);
				int32_t coor3 = cap_getCoordinate(cap3);
				int32_t start, end;
				if(coor5 < coor3){
					start = coor5 + tstart -2;
					end = coor3 + tstart +1 -2;
				}else{
					start = coor3 + tstart -2;
					end = coor5 + tstart +1 -2;
				}

				fprintf(fileHandle, "chr\t%d\t%d\n", start, end);
				
			}
			free(sequenceHeader);
		}
	}
	block_destructInstanceIterator(instanceIterator);
}

void getBED(Flower *net, FILE *fileHandle, char *query, char *target, int tstart, int s, int e) {
	/*
	 * Outputs target BED records of all the blocks in the net and its descendants.
	 */

	Flower_BlockIterator *blockIterator = flower_getBlockIterator(net);
	Block *block;
	while((block = flower_getNextBlock(blockIterator)) != NULL) {
        	//if(block_getChain(block) == NULL){//ONLY print out NON CHAIN blocks
        	if(block_hasQueryTarget(block, query, target)){
                        if( (s < 0 || e < 0) || block_inRange(block, query, s, e) ){
				getBEDBlock(block, fileHandle, target, tstart);
			}
		}
        	//}
	}
	flower_destructBlockIterator(blockIterator);

	//Call child nets recursively.
	Flower_GroupIterator *groupIterator = flower_getGroupIterator(net);
	Group *group;
	while((group = flower_getNextGroup(groupIterator)) != NULL) {
		Flower *nestedNet = group_getNestedNet(group);
		if(nestedNet != NULL) {
			getBED(group_getNestedNet(group), fileHandle, query, target, tstart, s, e); //recursive call.
		}
	}
	flower_destructGroupIterator(groupIterator);
}

void getBEDs(Flower *net, FILE *fileHandle, char *query, char *target, int tstart, int offset, struct psl *refpsl){
   int start, end;
   while(refpsl != NULL){
      start = refpsl->tStart - offset +2;
      end = refpsl->tEnd - offset +2;
      getBED(net, fileHandle, query, target, tstart, start, end);
      refpsl = refpsl->next;
   }
}

void usage() {
	fprintf(stderr, "cactus_bedGenerator, version 0.2\n");
	fprintf(stderr, "Prints to output file all segments of the target sequence that are in blocks that contain both query & target\n");
	fprintf(stderr, "-a --logLevel : Set the log level\n");
	fprintf(stderr, "-c --netDisk : The location of the net disk directory\n");
	fprintf(stderr, "-d --netName : The name of the net (the key in the database)\n");
	fprintf(stderr, "-e --outputFile : The file to write the BEDs in.\n");
	fprintf(stderr, "-r --ref : file that contains psl records of reference genes, whose target is the same with query.\n");
	fprintf(stderr, "-s --tstart : target sequence position on the chromosome\n");
	fprintf(stderr, "-o --qstart : query sequence position on the chromosome\n");
	fprintf(stderr, "-q --query\n");
	fprintf(stderr, "-t --target\n");
	fprintf(stderr, "-h --help : Print this help screen\n");
}

int main(int argc, char *argv[]) {
	CactusDisk *netDisk;
	Flower *net;

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
	int tstart = 0;
	int qstart = 0;

	///////////////////////////////////////////////////////////////////////////
	// (0) Parse the inputs handed by genomeCactus.py / setup stuff.
	///////////////////////////////////////////////////////////////////////////

	while(1) {
		static struct option long_options[] = {
			{ "logLevel", required_argument, 0, 'a' },
			{ "netDisk", required_argument, 0, 'c' },
			{ "netName", required_argument, 0, 'd' },
			{ "outputFile", required_argument, 0, 'e' },
			{ "query", required_argument, 0, 'q' },
			{ "target", required_argument, 0, 't' },
			{ "tstart", required_argument, 0, 's' },
			{ "qstart", required_argument, 0, 'o' },
			{ "ref", required_argument, 0, 'r' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "r:o:a:c:d:e:q:t:s:h", long_options, &option_index);

		if(key == -1) {
			break;
		}

		switch(key) {
			case 'a':
				logLevelString = stString_copy(optarg);
				break;
			case 'c':
				netDiskName = stString_copy(optarg);
				break;
			case 'd':
				netName = stString_copy(optarg);
				break;
			case 'e':
				outputFile = stString_copy(optarg);
				break;
			case 'q':
				query = stString_copy(optarg);
				break;
			case 't':
				target = stString_copy(optarg);
				break;
			case 's':
				sscanf( stString_copy(optarg), "%d", &tstart);
				break;
			case 'o':
				sscanf( stString_copy(optarg), "%d", &qstart);
				break;
			case 'r':
				ref = stString_copy(optarg);
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

	st_logInfo("Net disk name : %s\n", netDiskName);
	st_logInfo("Net name : %s\n", netName);
	st_logInfo("Output BED file : %s\n", outputFile);

	//////////////////////////////////////////////
	//Load the database
	//////////////////////////////////////////////

	netDisk = cactusDisk_construct(netDiskName);
	st_logInfo("Set up the net disk\n");

	///////////////////////////////////////////////////////////////////////////
	// Parse the basic reconstruction problem
	///////////////////////////////////////////////////////////////////////////

	net = cactusDisk_getNet(netDisk, cactusMisc_stringToName(netName));
	st_logInfo("Parsed the top level net of the cactus tree to check\n");

	///////////////////////////////////////////////////////////////////////////
	// Recursive check the nets.
	///////////////////////////////////////////////////////////////////////////

	int64_t startTime = time(NULL);
	FILE *fileHandle = fopen(outputFile, "w");
        if(ref == NULL){//no limit, print all bed records
           getBED(net, fileHandle, query, target, tstart, -1, -1);
        }else{
           struct psl *refpsl = pslLoadAll(ref);
	   getBEDs(net, fileHandle, query, target, tstart, qstart, refpsl);
        }
	fclose(fileHandle);
	st_logInfo("Got the beds in %i seconds/\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	cactusDisk_destruct(netDisk);

	return 0;
}
