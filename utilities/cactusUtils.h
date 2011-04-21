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

/////////////////////////////////////////////////////
//Common functions used by cactus utilities-scripts// 
/////////////////////////////////////////////////////

/*
 *Return the fasta header of the sequence if any. Otherwise return the sequence's internal name.
 */
char *formatSequenceHeader(Sequence *sequence);


/*
 *Return the first sequence in flower whose name is 'header'
 */
Sequence *getSequenceByHeader(Flower *flower, char *header);

/*
 *Return the first sequence in flower whose name matches 'header'
 */
Sequence *getSequenceMatchesHeader(Flower *flower, char *header);

/*
 *Return the list of Sequences in flower whose names mathces 'name'.
 */
struct List *getSequences(Flower *flower, char *name);

/*
 *Return true if cap is a stubEnd, otherwise return False
 */
bool isStubCap(Cap *cap);

/*
 *Return the sequence header (or it's internal name if header is not available) of the input cap
 */
char *cap_getSequenceName(Cap *cap);

/*
 *Append the integer 'num' to string 'name'. Return the concatenated string.
 */
char *appendIntToName(char *name, int num);

/*
 *Return length of the sequence in flower whose name is header
 */
int32_t getSeqLength(Flower *flower, char *header);

/*
 *Move cap to the next block that its segment is adjacency to
 */
void moveCapToNextBlock(Cap **cap);

/*
 *Get 3'-end stubs (caps) of the forward strand of sequences by its name
 *(i.e each of these stubs is the first cap at the 5' end of each thread whose name matches 'name'
 */
struct List *flower_getThreadStarts(Flower *flower, char *name);

/*
 *Split a string using 'delim'. Return the list of token strings
 */
struct List *splitString(char *str, char *delim);

struct IntList *splitIntString(char *str, char *delim);

char *str_joinList(struct List *strList, char *sep);


