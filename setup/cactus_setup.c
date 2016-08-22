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
#include <inttypes.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <dirent.h>
#include <math.h>
#include <ctype.h>

const char *CACTUS_SETUP_EXCEPTION = "CACTUS_SETUP_EXCEPTION";

#include "bioioC.h"
#include "cactus.h"

void usage() {
    fprintf(stderr, "cactus_setup [fastaFile]xN, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --cactusDisk : The location of the flower disk directory\n");
    fprintf(stderr, "-f --speciesTree : The species tree, which will form the skeleton of the event tree\n");
    fprintf(stderr, "-g --outgroupEvents : Leaf events in the species tree identified as outgroups\n");
    fprintf(stderr, "-i --makeEventHeadersAlphaNumeric : Remove non alpha-numeric characters from event header names\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(stderr, "-d --debug : Run some extra debug checks at the end\n");
}

/*
 * Plenty of global variables!
 */
int64_t isComplete = 1;
char * cactusDiskDatabaseString = NULL;
char * cactusSequencesPath = NULL;
CactusDisk *cactusDisk;
Flower *flower;
EventTree *eventTree;
Event *event;
int64_t totalSequenceNumber = 0;

void checkBranchLengthsAreDefined(stTree *tree) {
    if (isinf(stTree_getBranchLength(tree))) {
        st_errAbort("Got a non defined branch length in the input tree: %s.\n", stTree_getNewickTreeString(tree));
    }
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        checkBranchLengthsAreDefined(stTree_getChild(tree, i));
    }
}

char *makeAlphaNumeric(const char *string) {
    char *cA = stString_copy(string);
    int64_t j = 0;
    for (int64_t i = 0; i < strlen(string); i++) {
        if (isalpha(string[i]) || isdigit(string[i])) {
            cA[j++] = string[i];
        }
    }
    cA[j] = '\0';
    return cA;
}

void makeEventHeadersAlphaNumericFn(stTree *tree) {
    char *cA = makeAlphaNumeric(stTree_getLabel(tree));
    stTree_setLabel(tree, cA);
    free(cA);
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        makeEventHeadersAlphaNumericFn(stTree_getChild(tree, i));
    }
}

void processSequence(const char *fastaHeader, const char *string, int64_t length) {
    /*
     * Processes a sequence by adding it to the flower disk.
     */
    End *end1;
    End *end2;
    Cap *cap1;
    Cap *cap2;
    MetaSequence *metaSequence;
    Sequence *sequence;

    //Now put the details in a flower.
    metaSequence = metaSequence_construct(2, length, string, fastaHeader, event_getName(event), cactusDisk);
    sequence = sequence_construct(metaSequence, flower);

    end1 = end_construct2(0, isComplete, flower);
    end2 = end_construct2(1, isComplete, flower);
    cap1 = cap_construct2(end1, 1, 1, sequence);
    cap2 = cap_construct2(end2, length + 2, 1, sequence);
    cap_makeAdjacent(cap1, cap2);
    totalSequenceNumber++;
}

void setCompleteStatus(const char *fileName) {
    isComplete = 0;
    int64_t i = strlen(fileName);
    if (i >= 9) {
        const char *cA = fileName + i - 9;
        if (strcmp(cA, ".complete") == 0) {
            isComplete = 1;
            st_logInfo("The file %s is specified complete, the sequences will be attached\n", fileName);
            return;
        }
    }
    if (i >= 12) {
        const char *cA = fileName + i - 12;
        if (strcmp(cA, ".complete.fa") == 0) {
            isComplete = 1;
            st_logInfo("The file %s is specified complete, the sequences will be attached\n", fileName);
            return;
        }
    }
    st_logInfo("The file %s is specified incomplete, the sequences will not be attached\n", fileName);
}

int main(int argc, char *argv[]) {
    /*
     * Open the database.
     * Construct a flower.
     * Construct an event tree representing the species tree.
     * For each sequence contruct two ends each containing an cap.
     * Make a file for the sequence.
     * Link the two caps.
     * Finish!
     */

    int64_t key, j;
    struct List *stack;
    FILE *fileHandle = NULL;
    int64_t totalEventNumber;
    Group *group;
    Flower_EndIterator *endIterator;
    End *end;
    bool makeEventHeadersAlphaNumeric = 0;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * speciesTree = NULL;
    char * outgroupEvents = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'a' }, { "cactusDisk", required_argument, 0, 'b' }, {"cactusSequencesPath", required_argument, 0, 'f'}, {
                "speciesTree", required_argument, 0, 'g' }, { "outgroupEvents", required_argument, 0, 'h' },
                { "help", no_argument, 0, 'i' }, { "makeEventHeadersAlphaNumeric", no_argument, 0, 'j' }, { 0, 0, 0, 0 } };

        int option_index = 0;

        key = getopt_long(argc, argv, "a:b:f:hg:i", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = optarg;
                break;
            case 'b':
                cactusDiskDatabaseString = optarg;
                break;
            case 'f':
                cactusSequencesPath = optarg;
                break;
            case 'g':
                speciesTree = optarg;
                break;
            case 'h':
                outgroupEvents = optarg;
                break;
            case 'i':
                usage();
                return 0;
            case 'j':
                makeEventHeadersAlphaNumeric = 1;
                break;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    //assert(logLevelString == NULL || strcmp(logLevelString, "CRITICAL") == 0 || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
    assert(cactusDiskDatabaseString != NULL);
    assert(speciesTree != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Flower disk name : %s\n", cactusDiskDatabaseString);

    for (j = optind; j < argc; j++) {
        st_logInfo("Sequence file/directory %s\n", argv[j]);
    }

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    stKVDatabaseConf *kvDatabaseConf = kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskDatabaseString);
    if (stKVDatabaseConf_getType(kvDatabaseConf) == stKVDatabaseTypeTokyoCabinet || stKVDatabaseConf_getType(kvDatabaseConf)
            == stKVDatabaseTypeKyotoTycoon) {
        assert(stKVDatabaseConf_getDir(kvDatabaseConf) != NULL);
        cactusDisk = cactusDisk_construct2(kvDatabaseConf, cactusSequencesPath);
    } else {
        cactusDisk = cactusDisk_construct(kvDatabaseConf, 1);
    }
    st_logInfo("Set up the flower disk\n");

    //////////////////////////////////////////////
    //Construct the flower
    //////////////////////////////////////////////

    if (cactusDisk_getFlower(cactusDisk, 0) != NULL) {
        cactusDisk_destruct(cactusDisk);
        st_logInfo("The first flower already exists\n");
        return 0;
    }
    flower = flower_construct2(0, cactusDisk);
    assert(flower_getName(flower) == 0);
    st_logInfo("Constructed the flower\n");

    //////////////////////////////////////////////
    //Construct the event tree
    //////////////////////////////////////////////

    st_logInfo("Going to build the event tree with newick string: %s\n", speciesTree);
    stTree *tree = stTree_parseNewickString(speciesTree);
    st_logInfo("Parsed the tree\n");
    if (makeEventHeadersAlphaNumeric) {
        makeEventHeadersAlphaNumericFn(tree);
    }
    stTree_setBranchLength(tree, INT64_MAX);
    checkBranchLengthsAreDefined(tree);
    eventTree = eventTree_construct2(flower); //creates the event tree and the root even
    totalEventNumber = 1;
    st_logInfo("Constructed the basic event tree\n");

    //now traverse the tree
    stack = constructEmptyList(0, NULL);
    listAppend(stack, eventTree_getRootEvent(eventTree));
    listAppend(stack, tree);
    j = optind;
    while (stack->length > 0) {
        tree = stack->list[--stack->length];
        event = stack->list[--stack->length];
        assert(tree != NULL);
        totalEventNumber++;
        if (stTree_getChildNumber(tree) > 0) {
            event = event_construct3(stTree_getLabel(tree), stTree_getBranchLength(tree), event, eventTree);
            for (int64_t i = stTree_getChildNumber(tree) - 1; i >= 0; i--) {
                listAppend(stack, event);
                listAppend(stack, stTree_getChild(tree, i));
            }
        } else {
            assert(j < argc);
            assert(stTree_getLabel(tree) != NULL);

            assert(stTree_getBranchLength(tree) != INFINITY);
            event = event_construct3(stTree_getLabel(tree), stTree_getBranchLength(tree), event, eventTree);

            char *fileName = argv[j];

            if (!stFile_exists(fileName)) {
                st_errAbort("File does not exist: %s\n", fileName);
            }

            if (stFile_isDir(fileName)) {
                st_logInfo("Processing directory: %s\n", fileName);
                stList *filesInDir = stFile_getFileNamesInDirectory(fileName);
                for (int64_t i = 0; i < stList_length(filesInDir); i++) {
                    char *absChildFileName = stFile_pathJoin(fileName, stList_get(filesInDir, i));
                    assert(stFile_exists(absChildFileName));
                    setCompleteStatus(absChildFileName); //decide if the sequences in the file should be free or attached.
                    fileHandle = fopen(absChildFileName, "r");
                    fastaReadToFunction(fileHandle, processSequence);
                    fclose(fileHandle);
                    free(absChildFileName);
                }
                stList_destruct(filesInDir);
            } else {
                st_logInfo("Processing file: %s\n", fileName);
                setCompleteStatus(fileName); //decide if the sequences in the file should be free or attached.
                fileHandle = fopen(fileName, "r");
                fastaReadToFunction(fileHandle, processSequence);
                fclose(fileHandle);
            }
            j++;
        }
    }
    char *eventTreeString = eventTree_makeNewickString(eventTree);
    st_logInfo(
            "Constructed the initial flower with %" PRIi64 " sequences and %" PRIi64 " events with string: %s\n",
            totalSequenceNumber, totalEventNumber, eventTreeString);
    assert(event_getSubTreeBranchLength(eventTree_getRootEvent(eventTree)) >= 0.0);
    free(eventTreeString);
    //assert(0);

    //////////////////////////////////////////////
    //Label any outgroup events.
    //////////////////////////////////////////////

    if (outgroupEvents != NULL) {
        stList *outgroupEventsList = stString_split(outgroupEvents);
        for (int64_t i = 0; i < stList_length(outgroupEventsList); i++) {
            char *outgroupEvent = makeEventHeadersAlphaNumeric ? makeAlphaNumeric(stList_get(outgroupEventsList, i)) : stString_copy(stList_get(outgroupEventsList, i));
            Event *event = eventTree_getEventByHeader(eventTree, outgroupEvent);
            if (event == NULL) {
                st_errAbort("Got an outgroup string that does not match an event, outgroup string %s", outgroupEvent);
            }
            if (event_getChildNumber(event) != 0) {
                st_errAbort("Attempting to label an internal node as an outgroup, outgroup string %s", outgroupEvent);
            }
            assert(!event_isOutgroup(event));
            event_setOutgroupStatus(event, 1);
            assert(event_isOutgroup(event));
            free(outgroupEvent);
        }
        stList_destruct(outgroupEventsList);
    }

    //////////////////////////////////////////////
    //Construct the terminal group.
    //////////////////////////////////////////////

    if (flower_getEndNumber(flower) > 0) {
        group = group_construct2(flower);
        endIterator = flower_getEndIterator(flower);
        while ((end = flower_getNextEnd(endIterator)) != NULL) {
            end_setGroup(end, group);
        }
        flower_destructEndIterator(endIterator);
        assert(group_isLeaf(group));

        // Create a one link chain if there is only one pair of attached ends..
        group_constructChainForLink(group);
        assert(!flower_builtBlocks(flower));
    } else {
        flower_setBuiltBlocks(flower, 1);
    }

    ///////////////////////////////////////////////////////////////////////////
    // Write the flower to disk.
    ///////////////////////////////////////////////////////////////////////////

    //flower_check(flower);
    cactusDisk_write(cactusDisk);
    st_logInfo("Updated the flower on disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup.
    ///////////////////////////////////////////////////////////////////////////

    return 0; //Exit without clean up is quicker, enable cleanup when doing memory leak detection.

    stTree_destruct(tree);
    cactusDisk_destruct(cactusDisk);
    stKVDatabaseConf_destruct(kvDatabaseConf);

    return 0;
}
