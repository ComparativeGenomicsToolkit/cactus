/*
 * Copyright (C) 2009-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <getopt.h>

#include "cactus.h"
#include "avl.h"
#include "commonC.h"
#include "hashTableC.h"

///////////////////////////////////////////////////////
///Common functions used by cactus utilities-scripts///
///////////////////////////////////////////////////////

char *formatSequenceHeader(Sequence *sequence) {
    const char *sequenceHeader = sequence_getHeader(sequence);
    if (strlen(sequenceHeader) > 0) {
        char *cA = st_malloc(sizeof(char) * (1 + strlen(sequenceHeader)));
        sscanf(sequenceHeader, "%s", cA);
        return cA;
    } else {
        return cactusMisc_nameToString(sequence_getName(sequence));
    }
}


Sequence *getSequenceByHeader(Flower *flower, char *header){
    /*
     *Iterates through the Sequences in 'flower' and return 
     *the first Sequence whose name is 'header'
     */
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence = NULL;
    while((sequence = flower_getNextSequence(it)) != NULL){
        char *sequenceHeader = formatSequenceHeader(sequence);
        if(strcmp(sequenceHeader, header) == 0){
            free(sequenceHeader);
            break;
        }
        free(sequenceHeader);
    }
    flower_destructSequenceIterator(it);
    return sequence;
}

Sequence *getSequenceMatchesHeader(Flower *flower, char *header){
    //Returns the first Sequence whose name matches 'header'
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence = NULL;
    while((sequence = flower_getNextSequence(it)) != NULL){
        char *sequenceHeader = formatSequenceHeader(sequence);
        if(strstr(sequenceHeader, header) != NULL){
            free(sequenceHeader);
            break;
        }
        free(sequenceHeader);
    }
    flower_destructSequenceIterator(it);
    return sequence;
}

struct List *getSequences(Flower *flower, char *name){
   //get names of all the sequences in 'flower' that have their names start with 'name'
   Sequence *sequence;
   struct List *seqs = constructEmptyList(0, free);
   Flower_SequenceIterator * seqIterator = flower_getSequenceIterator(flower);
   while((sequence = flower_getNextSequence(seqIterator)) != NULL){
      char *sequenceHeader = formatSequenceHeader(sequence);
      //if(strstr(sequenceHeader, name) == sequenceHeader){
      if(strstr(sequenceHeader, name) != NULL){
         listAppend(seqs, sequenceHeader);
      }
      free(sequenceHeader);
   }
   flower_destructSequenceIterator(seqIterator);
   return seqs;
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

char *cap_getSequenceName(Cap *cap){
    Sequence *sequence = cap_getSequence(cap);
    char *name = NULL;
    if(sequence != NULL){
        name = formatSequenceHeader(sequence);
    }
    return name;
}

char *appendIntToName(char *name, int num){
    char index[10];//assume max num is 999! 
    sprintf(index, "%d", num);
    char *newName = st_malloc(sizeof(char)*(strlen(name) + strlen(index) + 2));
    strcpy(newName, name);
    strcat(newName, "_");
    strcat(newName, index);
    return newName;
}

int32_t getSeqLength(Flower *flower, char *header){
    Flower_SequenceIterator *it = flower_getSequenceIterator(flower);
    Sequence *sequence;
    while((sequence = flower_getNextSequence(it)) != NULL){
        char *sequenceHeader = formatSequenceHeader(sequence);
        if(strcmp(sequenceHeader, header) == 0){
            flower_destructSequenceIterator(it);
            free(sequenceHeader);
            return sequence_getLength(sequence);
        }
        free(sequenceHeader);
    }
    flower_destructSequenceIterator(it);
    return 0;
}

void moveCapToNextBlock(Cap **cap){
    /*Move cap to the next block that its segment is adjacency to*/
    assert((*cap) != NULL);
    if(isStubCap(*cap)){
        st_logInfo("STUB\n");
        return;
    }
    Cap *otherCap = cap_getOtherSegmentCap(*cap);
    *cap = cap_getAdjacency(otherCap);
    assert(*cap != NULL);
    st_logInfo("moved-cap %s, %d\n", cactusMisc_nameToString(cap_getName(*cap)), cap_getCoordinate(*cap));
}

struct List *flower_getThreadStarts(Flower *flower, char *name){
    /*
     *Get 3' end Stubs of the forward strand of sequences by its name
     *I.e (Each of these stubs is the first cap at the 5' end of each thread with 
     *name including 'name')
     */
    Cap *cap;
    struct List *startCaps = constructEmptyList(0, free);
    Flower_CapIterator *capIterator = flower_getCapIterator(flower);
    while((cap= flower_getNextCap(capIterator)) != NULL){
        if(isStubCap(cap)){//dead end or inherited end
            char *sequenceHeader = cap_getSequenceName(cap);
            if(sequenceHeader == NULL){continue;}
            if(strstr(sequenceHeader, name) != NULL){//cap matched with name
                if(!cap_getStrand(cap)){//if cap is on negative strand - reverse it
                    cap = cap_getReverse(cap);
                }
                if(cap_getSide(cap)){ continue; }//if cap is the 5' end, ignore
                listAppend(startCaps, cap);
            }
            free(sequenceHeader);
        }
    }
    flower_destructCapIterator(capIterator);
    return startCaps;
}


struct List *splitString(char *str, char *delim){
    struct List *list = constructEmptyList(0, free);
    char *tok;
    tok = strtok(stString_copy(str), delim);
    while (tok != NULL){
        listAppend(list, tok);
        tok = strtok(NULL, delim);
    }
    return list;
}

struct IntList *splitIntString(char *str, char *delim){
    struct IntList *list = constructEmptyIntList(0);
    char *tok;
    int32_t i;
    tok = strtok(stString_copy(str), delim);
    while (tok != NULL){
        assert( sscanf(tok, "%d", &i) == 1);
        intListAppend(list, i);
        tok = strtok(NULL, delim);
    }
    return list;
}

char *str_joinList(struct List *strList, char *sep){
    /*
     *Return a concatenated string of a list of string, which separated by 'sep' 
     */
    if(strList->length == 0){return NULL;}
    int32_t strLength = 0;
    for(int i=0; i < strList->length; i++){
        strLength += strlen(strList->list[i]);
    }
    //adding the strList->length -1 separators and end of str
    strLength += (strList->length - 1)*strlen(sep) + 1;
    char *str = st_malloc(sizeof(char)*strLength);

    strcpy(str, strList->list[0]);
    for(int i=1; i < strList->length; i++){
        strcat(str, sep);
        strcat(str, strList->list[i]);
    }
    return str;
}




