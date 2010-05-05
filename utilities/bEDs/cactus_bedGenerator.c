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
 * The script outputs a maf file containing all the block in a net and its descendants.
 */

/*
 * Global variables.
 */
bool includeTreesInMafBlocks = 0;

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
void printCap(FILE *fileHandle, Cap *cap){
   char side = cap_getSide(cap) ? '5' : '3';
   char strand = cap_getStrand(cap) ? '+' : '-';
   int coor = cap_getCoordinate(cap);
   Segment *segment = cap_getSegment(cap);
   Cap *otherCap = cap_getSide(cap) ? segment_get3Cap(segment): segment_get5Cap(segment);
   int otherCoor = cap_getCoordinate(otherCap);
   fprintf(fileHandle, "%s", netMisc_nameToString(cap_getName(cap)));
   fprintf(fileHandle, " %c%c%d (%d)", side, strand, coor, otherCoor);
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

bool block_hasQueryTarget(Block *block, char *query, char *target){
   //Return TRUE if contains segments of BOTH query and target
   Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
   Segment *segment;
   bool hasQuery = false;
   bool hasTarget = false;
   while((segment = block_getNext(instanceIterator)) != NULL){
      Sequence *sequence = segment_getSequence(segment);
      if(sequence == NULL){continue;}
      if(sequence != NULL){ 
         char *sequenceHeader = formatSequenceHeader(sequence);
         if(strcmp(sequenceHeader, query) == 0){
            hasQuery = true;
         }else if(strcmp(sequenceHeader, target) == 0){
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
	 * Outputs a MAF representation of the block to the given file handle.
	 */
	Block_InstanceIterator *instanceIterator = block_getInstanceIterator(block);
	Segment *segment;
	while((segment = block_getNext(instanceIterator)) != NULL) {
		Sequence *sequence = segment_getSequence(segment);
		if(sequence != NULL) {
			char *sequenceHeader = formatSequenceHeader(sequence);
         		if(strcmp(sequenceHeader, target) == 0){
                                Cap *cap5 = segment_get5Cap(segment);
                                Cap *cap3 = segment_get3Cap(segment);
				int32_t coor5 = cap_getCoordinate(cap5);
				int32_t coor3 = cap_getCoordinate(cap3);
				//int32_t seqlength = sequence_getLength(sequence);
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

void getBED(Net *net, FILE *fileHandle, char *query, char *target, int tstart, int s, int e) {
	/*
	 * Outputs MAF representations of all the block sin the net and its descendants.
	 */

	//Make MAF blocks for each block
	Net_BlockIterator *blockIterator = net_getBlockIterator(net);
	Block *block;
	while((block = net_getNextBlock(blockIterator)) != NULL) {
        	//if(block_getChain(block) == NULL){//ONLY print out NON CHAIN blocks
        	if(block_hasQueryTarget(block, query, target)){
                        if(s < 0 || e < 0){
				getBEDBlock(block, fileHandle, target, tstart);
                	}else if(block_inRange(block, query, s, e)){
				getBEDBlock(block, fileHandle, target, tstart);
			}
		}
        	//}
	}
	net_destructBlockIterator(blockIterator);

	//Call child nets recursively.
	Net_GroupIterator *groupIterator = net_getGroupIterator(net);
	Group *group;
	while((group = net_getNextGroup(groupIterator)) != NULL) {
		Net *nestedNet = group_getNestedNet(group);
		if(nestedNet != NULL) {
			getBED(group_getNestedNet(group), fileHandle, query, target, tstart, s, e); //recursive call.
		}
	}
	net_destructGroupIterator(groupIterator);
}

void getBEDs(Net *net, FILE *fileHandle, char *query, char *target, int tstart, int offset, struct psl *refpsl){
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
	fprintf(stderr, "-r --ref : psl of reference genes\n");
	fprintf(stderr, "-s --tstart : target sequence position on the chromosome\n");
	fprintf(stderr, "-o --qstart : query sequence position on the chromosome\n");
	fprintf(stderr, "-q --query\n");
	fprintf(stderr, "-t --target\n");
	fprintf(stderr, "-f --includeTrees : Include trees for each MAF block inside of a comment line.\n");
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
			{ "includeTrees", no_argument, 0, 'f' },
			{ "help", no_argument, 0, 'h' },
			{ 0, 0, 0, 0 }
		};

		int option_index = 0;

		int key = getopt_long(argc, argv, "r:o:a:c:d:e:q:t:s:fh", long_options, &option_index);

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
			case 's':
				sscanf( stringCopy(optarg), "%d", &tstart);
				break;
			case 'o':
				sscanf( stringCopy(optarg), "%d", &qstart);
				break;
			case 'r':
				ref = stringCopy(optarg);
				break;
			case 'f':
				includeTreesInMafBlocks = !includeTreesInMafBlocks;
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
	logInfo("Output MAF file : %s\n", outputFile);

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
        if(ref == NULL){//no limit, print all bed records
           getBED(net, fileHandle, query, target, tstart, -1, -1);
        }else{
           struct psl *refpsl = pslLoadAll(ref);
	   getBEDs(net, fileHandle, query, target, tstart, qstart, refpsl);
        }
	fclose(fileHandle);
	logInfo("Got the mafs in %i seconds/\n", time(NULL) - startTime);

	///////////////////////////////////////////////////////////////////////////
	// Clean up.
	///////////////////////////////////////////////////////////////////////////

	netDisk_destruct(netDisk);

	return 0;
}
