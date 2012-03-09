#include <stdio.h>
#include <ctype.h>
#include "cactus.h"
#include "sonLib.h"

////////////////////////////////////
////////////////////////////////////
//Calculate the bases of a reference sequence block
////////////////////////////////////
////////////////////////////////////


static int32_t *collateCounts(stList *strings, int32_t blockLength, int32_t **upperCounts) {
    //Matrix to store the number of occurrences of each base type, for each column in the block
    int32_t *baseCounts = st_calloc(blockLength * 5, sizeof(int32_t));
    //Array storing number of bases that are upper case letters (non repetitive)..
    *upperCounts = st_calloc(blockLength, sizeof(int32_t));
    for (int32_t j = 0; j < stList_length(strings); j++) {
        char *string = stList_get(strings, j);
        for (int32_t i = 0; i < blockLength; i++) {
            (*upperCounts)[i] += toupper(string[i]) == string[i] ? 1 : 0;
            switch (toupper(string[i])) {
                case 'A':
                    baseCounts[i * 5]++;
                    break;
                case 'C':
                    baseCounts[i * 5 + 1]++;
                    break;
                case 'G':
                    baseCounts[i * 5 + 2]++;
                    break;
                case 'T':
                    baseCounts[i * 5 + 3]++;
                    break;
                default:
                    baseCounts[i * 5 + 4]++;
            }
        }
    }
    return baseCounts;
}

static char getMajorityBase(int32_t *baseCounts, int32_t *outgroupBaseCounts) {
    static const char bases[5] = { 'a', 'c', 'g', 't' };
    static char seq[5];

    int32_t maxBaseCount = 0;
    for (int32_t j = 0; j < 4; j++) {
        int32_t k = baseCounts[j] * 2 + outgroupBaseCounts[j];
        if (maxBaseCount < k) {
            maxBaseCount = k;
        }
    }
    //st_uglyf("max base count %i\n", maxBaseCount);
    if (maxBaseCount == 0) {
        return 'n';
    } else {
        int32_t k = 0;
        for (int32_t j = 0; j < 4; j++) {
            if (baseCounts[j] * 2 + outgroupBaseCounts[j] == maxBaseCount) {
                seq[k++] = bases[j];
            }
        }
        seq[k] = '\0';
        return seq[st_randomInt(0, strlen(seq))];
    }
}

char *getConsensusStringP(stList *strings, stList *outgroupStrings, int32_t blockLength) {
    int32_t *upperCounts = NULL;
    int32_t *baseCounts = collateCounts(strings, blockLength, &upperCounts);
    int32_t *upperCountsOutgroup = NULL;
    int32_t *baseCountsOutgroup = collateCounts(outgroupStrings, blockLength, &upperCountsOutgroup);

    char *string = st_malloc(sizeof(char) * (blockLength + 1));
    string[blockLength] = '\0';

    for (int32_t i = 0; i < blockLength; i++) {
        /*
         * Logic to choose base..
         */
        char base = getMajorityBase(&(baseCounts[i * 5]), &(baseCountsOutgroup[i * 5]));

        /*
         * Now choose if it should be upper case.
         */
        string[i] = (upperCounts[i] + upperCountsOutgroup[i]) > ((double) stList_length(strings) + stList_length(
                outgroupStrings)) / 2 ? toupper(base) : base;
    }
    free(baseCounts);
    free(upperCounts);

    return string;
}

char *getConsensusString(Block *block, Name outgroupEventName) {
    /*
     * Returns a consensus string for a block.
     */

    //Get the strings.
    Segment *segment;
    Block_InstanceIterator *instanceIt = block_getInstanceIterator(block);
    stList *strings = stList_construct3(0, free);
    stList *outgroupStrings = stList_construct3(0, free);
    while ((segment = block_getNext(instanceIt)) != NULL) {
        if (segment_getSequence(segment) != NULL) {
            if (event_getName(segment_getEvent(segment)) == outgroupEventName) {
                stList_append(outgroupStrings, segment_getString(segment));
            } else {
                stList_append(strings, segment_getString(segment));
            }
        }
    }
    block_destructInstanceIterator(instanceIt);

    char *string = getConsensusStringP(strings, outgroupStrings, block_getLength(block));
    stList_destruct(strings);
    stList_destruct(outgroupStrings);
    return string;
}

void reverseComplementString(char *string) {

}
