#include "cactus.h"
#include "sonLib.h"
#include "bioioC.h"
#include <stdio.h>
#include <ctype.h>

void checkBranchLengthsAreDefined(stTree *tree) {
    if (isinf(stTree_getBranchLength(tree))) {
        st_errAbort("Got a non defined branch length in the input tree: %s.\n", stTree_getNewickTreeString(tree));
    }
    if (stTree_getBranchLength(tree) == 0.0) {
        stTree_setBranchLength(tree, DBL_MIN);
        st_logCritical("0-length branch found for %s. Resetting to very small value\n", stTree_getLabel(tree));
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

bool getCompleteStatus(const char *fileName) {
    int64_t i = strlen(fileName);
    if (i >= 9) {
        const char *cA = fileName + i - 9;
        if (strcmp(cA, ".complete") == 0) {
            st_logInfo("The file %s is specified complete, the sequences will be attached\n", fileName);
            return 1;
        }
    }
    if (i >= 12) {
        const char *cA = fileName + i - 12;
        if (strcmp(cA, ".complete.fa") == 0) {
            st_logInfo("The file %s is specified complete, the sequences will be attached\n", fileName);
            return 1;
        }
    }
    st_logInfo("The file %s is specified incomplete, the sequences will not be attached\n", fileName);
    return 0;
}

typedef struct _processSequenceVars {
    bool isComplete;
    Event *event;
    int64_t totalSequenceNumber;
    Flower *flower;
    CactusDisk *cactusDisk;
} ProcessSequenceVars;

void processSequence(void* destination, const char *fastaHeader, const char *string, int64_t length) {
    /*
     * Processes a sequence by adding it to the flower disk.
     */
    //Now put the details in a flower.
    ProcessSequenceVars *p = destination;
    Sequence *sequence = sequence_construct(2, length, string, fastaHeader, p->event, p->cactusDisk);
    flower_addSequence(p->flower, sequence);

    End *end1 = end_construct2(0, p->isComplete, p->flower);
    End *end2 = end_construct2(1, p->isComplete, p->flower);
    Cap *cap1 = cap_construct2(end1, 1, 1, sequence);
    Cap *cap2 = cap_construct2(end2, length + 2, 1, sequence);
    cap_makeAdjacent(cap1, cap2);
    p->totalSequenceNumber++;
}

static int64_t assignSequences(CactusDisk *cactusDisk, Flower *flower, EventTree *eventTree, char *sequenceFilesAndEvents) {
    stList *sequenceFilesAndEventsList = stString_split(sequenceFilesAndEvents);
    if (stList_length(sequenceFilesAndEventsList) % 2 != 0) {
        stList_destruct(sequenceFilesAndEventsList);
        st_errAbort("Sequences weren't provided in a proper "
                    "'event seq' space-separated format");
    }
    ProcessSequenceVars p; // Struct to pass around storing variables for
    // making sequences
    p.totalSequenceNumber = 0;
    p.flower = flower;
    p.cactusDisk = cactusDisk;

    for (int64_t i = 0; i < stList_length(sequenceFilesAndEventsList); i += 2) {
        char *eventName = stList_get(sequenceFilesAndEventsList, i);
        char *fileName = stList_get(sequenceFilesAndEventsList, i+1);

        st_logInfo("Assigning sequence %s to %s\n", fileName, eventName);

        if (!stFile_exists(fileName)) {
            st_errAbort("File does not exist: %s\n", fileName);
        }

        // Set the global "event" variable, which is needed for the
        // function provided to fastaReadToFunction.
        p.event = eventTree_getEventByHeader(eventTree, eventName);
        if (p.event == NULL) {
            st_errAbort("No such event: %s", eventName);
        }
        if (stFile_isDir(fileName)) {
            st_logInfo("Processing directory: %s\n", fileName);
            stList *filesInDir = stFile_getFileNamesInDirectory(fileName);
            for (int64_t j = 0; j < stList_length(filesInDir); j++) {
                char *absChildFileName = stFile_pathJoin(fileName, stList_get(filesInDir, j));
                assert(stFile_exists(absChildFileName));
                p.isComplete = getCompleteStatus(absChildFileName); //decide if the sequences in the file should be free or attached.
                FILE *fileHandle = fopen(absChildFileName, "r");
                fastaReadToFunction(fileHandle, &p, processSequence);
                fclose(fileHandle);
                free(absChildFileName);
            }
            stList_destruct(filesInDir);
        } else {
            st_logInfo("Processing file: %s\n", fileName);
            p.isComplete = getCompleteStatus(fileName); //decide if the sequences in the file should be free or attached.
            FILE *fileHandle = fopen(fileName, "r");
            fastaReadToFunction(fileHandle, &p, processSequence);
            fclose(fileHandle);
        }
    }
    stList_destruct(sequenceFilesAndEventsList);

    return p.totalSequenceNumber;
}

static int64_t constructEvents(Event *parentEvent, stTree *tree, EventTree *eventTree) {
    Event *myEvent = NULL; // To distinguish from the global "event" variable.
    assert(tree != NULL);
    int64_t totalEventNumber = 1;
    myEvent = event_construct3(stTree_getLabel(tree),
                               stTree_getBranchLength(tree), parentEvent,
                               eventTree);
    for (int64_t i = 0; i < stTree_getChildNumber(tree); i++) {
        totalEventNumber += constructEvents(myEvent, stTree_getChild(tree, i), eventTree);
    }
    return totalEventNumber;
}

Flower *cactus_setup_first_flower(CactusDisk *cactusDisk, CactusParams *params,
                                  char *speciesTree, char *outgroupEvents, char *sequenceFilesAndEvents) {
    /*
     * Flow is:
     * Construct a flower.
     * Construct an event tree representing the species tree.
     * For each sequence contruct two ends each containing a cap.
     * Make a file for the sequence.
     * Link the two caps.
     * Finish!
     */

    //////////////////////////////////////////////
    //Construct the flower
    //////////////////////////////////////////////

    if (cactusDisk_getFlower(cactusDisk, 0) != NULL) {
        cactusDisk_destruct(cactusDisk);
        st_logInfo("The first flower already exists\n");
        return 0;
    }
    Flower *flower = flower_construct2(0, cactusDisk);
    assert(flower_getName(flower) == 0);
    st_logInfo("Constructed the first flower\n");

    //////////////////////////////////////////////
    //Construct the event tree
    //////////////////////////////////////////////

    st_logInfo("Going to build the event tree with newick string: %s\n", speciesTree);
    stTree *tree = stTree_parseNewickString(speciesTree);
    st_logInfo("Parsed the tree\n");
    bool makeEventHeadersAlphaNumeric = cactusParams_get_int(params, 2, "setup", "makeEventHeadersAlphaNumeric");
    if (makeEventHeadersAlphaNumeric) {
        makeEventHeadersAlphaNumericFn(tree);
    }
    stTree_setBranchLength(tree, INT64_MAX);
    checkBranchLengthsAreDefined(tree);
    EventTree *eventTree = eventTree_construct2(cactusDisk); //creates the event tree and the root event
    st_logInfo("Constructed the basic event tree\n");

    // Construct a set of outgroup names so that ancestral outgroups
    // get recognized.
    stSet *outgroupNameSet = stSet_construct3(stHash_stringKey,
                                              stHash_stringEqualKey,
                                              free);
    if(outgroupEvents != NULL) {
        stList *outgroupNames = stString_split(outgroupEvents);
        for(int64_t i = 0; i < stList_length(outgroupNames); i++) {
            char *outgroupName = stList_get(outgroupNames, i);
            stSet_insert(outgroupNameSet, stString_copy(outgroupName));
        }
        stList_destruct(outgroupNames);
    }

    //now traverse the tree
    int64_t totalEventNumber = constructEvents(eventTree_getRootEvent(eventTree), tree, eventTree) + 1; // Plus one
    // for the root event

    //////////////////////////////////////////////
    //Construct the sequences and associate them with events
    //////////////////////////////////////////////

    int64_t totalSequenceNumber = assignSequences(cactusDisk, flower, eventTree, sequenceFilesAndEvents);

    //////////////////////////////////////////////
    //Log the constructed event tree and sequences
    //////////////////////////////////////////////

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
        Group *group = group_construct2(flower);
        Flower_EndIterator *endIterator = flower_getEndIterator(flower);
        End *end;
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

    //////////////////////////////////////////////
    //Cleanup.
    //////////////////////////////////////////////

    stTree_destruct(tree);
    stSet_destruct(outgroupNameSet);

    return flower;
}
