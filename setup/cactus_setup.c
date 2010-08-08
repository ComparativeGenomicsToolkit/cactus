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

#include "bioioC.h"
#include "cactus.h"

void usage() {
    fprintf(stderr, "cactus_setup [fastaFile]xN, version 0.2\n");
    fprintf(stderr, "-a --logLevel : Set the log level\n");
    fprintf(stderr, "-b --netDisk : The location of the net disk directory\n");
    fprintf(
            stderr,
            "-f --speciesTree : The species tree, which will form the skeleton of the event tree\n");
    fprintf(stderr, "-h --help : Print this help screen\n");
    fprintf(stderr, "-d --debug : Run some extra debug checks at the end\n");
}

/*
 * Plenty of global variables!
 */
int32_t isComplete = 1;
char * netDiskName = NULL;
CactusDisk *netDisk;
Flower *net;
EventTree *eventTree;
Event *event;
int32_t totalSequenceNumber = 0;

void fn(const char *fastaHeader, const char *string, int32_t length) {
    /*
     * Processes a sequence by adding it to the net disk.
     */
    End *end1;
    End *end2;
    Cap *cap1;
    Cap *cap2;
    MetaSequence *metaSequence;
    Sequence *sequence;

    //Now put the details in a net.
    metaSequence = metaSequence_construct(2, length, string, fastaHeader,
            event_getName(event), netDisk);
    sequence = sequence_construct(metaSequence, net);
    //isComplete = 0;
    end1 = end_construct2(0, isComplete, net);
    end2 = end_construct2(1, isComplete, net);
    cap1 = cap_construct2(end1, 1, 1, sequence);
    cap2 = cap_construct2(end2, length + 2, 1, sequence);
    cap_makeAdjacent(cap1, cap2);
    totalSequenceNumber++;
}

void setCompleteStatus(const char *fileName) {
    isComplete = 1;
    int32_t i = strlen(fileName);
    if (i > 11) {
        const char *cA = fileName + i - 11;
        if (strcmp(cA, ".incomplete") == 0) {
            isComplete = 0;
            st_logInfo(
                    "The file %s is specified incomplete, the sequences will be free\n",
                    fileName);
            return;
        }
    }
    st_logInfo(
            "The file %s is specified complete, the sequences will be attached\n",
            fileName);
}

int main(int argc, char *argv[]) {
    /*
     * Open the database.
     * Construct a net.
     * Construct an event tree representing the species tree.
     * For each sequence contruct two ends each containing an cap.
     * Make a file for the sequence.
     * Link the two caps.
     * Finish!
     */

    int32_t key, j;
    struct List *stack;
    struct BinaryTree *binaryTree;
    FILE *fileHandle;
    bool debug = 0;
    int32_t totalEventNumber;
    Group *group;
    Flower_EndIterator *endIterator;
    End *end;

    /*
     * Arguments/options
     */
    char * logLevelString = NULL;
    char * speciesTree = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // (0) Parse the inputs handed by genomeCactus.py / setup stuff.
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel",
                required_argument, 0, 'a' }, { "netDisk", required_argument, 0,
                'b' }, { "speciesTree", required_argument, 0, 'f' }, { "help",
                no_argument, 0, 'h' }, { "debug", no_argument, 0, 'd' }, { 0,
                0, 0, 0 } };

        int option_index = 0;

        key = getopt_long(argc, argv, "a:b:f:h:d", long_options, &option_index);

        if (key == -1) {
            break;
        }

        switch (key) {
            case 'a':
                logLevelString = optarg;
                break;
            case 'b':
                netDiskName = optarg;
                break;
            case 'f':
                speciesTree = optarg;
                break;
            case 'h':
                usage();
                return 0;
            case 'd':
                debug = 1;
                break;
            default:
                usage();
                return 1;
        }
    }

    ///////////////////////////////////////////////////////////////////////////
    // (0) Check the inputs.
    ///////////////////////////////////////////////////////////////////////////

    assert(logLevelString == NULL || strcmp(logLevelString, "CRITICAL") == 0 || strcmp(logLevelString, "INFO") == 0 || strcmp(logLevelString, "DEBUG") == 0);
    assert(netDiskName != NULL);
    assert(speciesTree != NULL);

    //////////////////////////////////////////////
    //Set up logging
    //////////////////////////////////////////////

    if (strcmp(logLevelString, "INFO") == 0) {
        st_setLogLevel(ST_LOGGING_INFO);
    }
    if (strcmp(logLevelString, "DEBUG") == 0) {
        st_setLogLevel(ST_LOGGING_DEBUG);
    }

    //////////////////////////////////////////////
    //Log (some of) the inputs
    //////////////////////////////////////////////

    st_logInfo("Net disk name : %s\n", netDiskName);

    for (j = optind; j < argc; j++) {
        st_logInfo("Sequence file/directory %s\n", argv[j]);
    }

    //////////////////////////////////////////////
    //Load the database
    //////////////////////////////////////////////

    netDisk = cactusDisk_construct(netDiskName);
    st_logInfo("Set up the net disk\n");

    //////////////////////////////////////////////
    //Construct the net
    //////////////////////////////////////////////

    if (cactusDisk_getFlowerNumberOnDisk(netDisk) != 0) {
        cactusDisk_destruct(netDisk);
        st_logInfo("The first net already exists\n");
        return 0;
    }
    net = flower_construct(netDisk);
    st_logInfo("Constructed the net\n");

    //////////////////////////////////////////////
    //Construct the event tree
    //////////////////////////////////////////////

    st_logInfo("Going to build the event tree with newick string: %s\n",
            speciesTree);
    binaryTree = newickTreeParser(speciesTree, 0.0, 0);
    binaryTree->distance = INT32_MAX;
    eventTree = eventTree_construct2(net); //creates the event tree and the root even
    totalEventNumber = 1;
    st_logInfo("Constructed the basic event tree\n");

    //now traverse the tree
    stack = constructEmptyList(0, NULL);
    listAppend(stack, eventTree_getRootEvent(eventTree));
    listAppend(stack, binaryTree);
    j = optind;
    while (stack->length > 0) {
        binaryTree = stack->list[--stack->length];
        event = stack->list[--stack->length];
        assert(binaryTree != NULL);
        totalEventNumber++;
        if (binaryTree->internal) {
            event = event_construct(metaEvent_construct(binaryTree->label,
                    netDisk), binaryTree->distance, event, eventTree);
            listAppend(stack, event);
            listAppend(stack, binaryTree->right);
            listAppend(stack, event);
            listAppend(stack, binaryTree->left);
        } else {
            assert(j < argc);
            event = event_construct(metaEvent_construct(binaryTree->label,
                    netDisk), binaryTree->distance, event, eventTree);

            struct stat info;//info about the file.
            exitOnFailure(stat(argv[j], &info),
                    "Failed to get information about the file: %s\n", argv[j]);
            if (S_ISDIR(info.st_mode)) {
                st_logInfo("Processing directory: %s\n", argv[j]);
                struct dirent *file;//a 'directory entity' AKA file
                DIR *dh = opendir(argv[j]);
                while ((file = readdir(dh)) != NULL) {
                    if (file->d_name[0] != '.') {
                        struct stat info2;
                        char *cA = pathJoin(argv[j], file->d_name);
                        //ascertain if complete or not
                        exitOnFailure(
                                stat(cA, &info2),
                                "Failed to get information about the file: %s\n",
                                file->d_name);
                        setCompleteStatus(file->d_name); //decide if the sequences in the file should be free or attached.
                        if (!S_ISDIR(info2.st_mode)) {
                            st_logInfo("Processing file: %s\n", cA);
                            fileHandle = fopen(cA, "r");
                            fastaReadToFunction(fileHandle, fn);
                            fclose(fileHandle);
                        }
                        free(cA);
                    }
                }
                closedir(dh);
            } else {
                st_logInfo("Processing file: %s\n", argv[j]);
                fileHandle = fopen(argv[j], "r");
                setCompleteStatus(argv[j]); //decide if the sequences in the file should be free or attached.
                fastaReadToFunction(fileHandle, fn);
                fclose(fileHandle);
            }
            j++;
        }
    }
    char *eventTreeString = eventTree_makeNewickString(eventTree);
    st_logInfo(
            "Constructed the initial net with %i sequences and %i events with string: %s\n",
            totalSequenceNumber, totalEventNumber, eventTreeString);
    assert(event_getSubTreeBranchLength(eventTree_getRootEvent(eventTree)) >= 0.0);
    free(eventTreeString);
    //assert(0);

    //////////////////////////////////////////////
    //Construct the terminal group.
    //////////////////////////////////////////////

    group = group_construct2(net);
    endIterator = flower_getEndIterator(net);
    while ((end = flower_getNextEnd(endIterator)) != NULL) {
        end_setGroup(end, group);
    }
    flower_destructEndIterator(endIterator);
    assert(group_isLeaf(group));

    ///////////////////////////////////////////////////////////////////////////
    // Create a one link chain if there is only one pair of attached ends..
    ///////////////////////////////////////////////////////////////////////////

    group_constructChainForLink(group);

    ///////////////////////////////////////////////////////////////////////////
    // Set the blocksBuilt flag to true for this top level net.
    ///////////////////////////////////////////////////////////////////////////

    assert(!flower_builtBlocks(net));
    flower_setBuiltBlocks(net, 1);
    assert(flower_builtBlocks(net));

    ///////////////////////////////////////////////////////////////////////////
    // Write the net to disk.
    ///////////////////////////////////////////////////////////////////////////

#ifdef BEN_DEBUG //Check we've done okay.
    flower_check(net);
#endif
    cactusDisk_write(netDisk);
    st_logInfo("Updated the net on disk\n");

    ///////////////////////////////////////////////////////////////////////////
    // Cleanup.
    ///////////////////////////////////////////////////////////////////////////

    cactusDisk_destruct(netDisk);

    ///////////////////////////////////////////////////////////////////////////
    // Debug
    ///////////////////////////////////////////////////////////////////////////

    if (debug) {
        netDisk = cactusDisk_construct(netDiskName);
        net = cactusDisk_getFlower(netDisk, 0);
        assert(net != NULL);
        assert(flower_getSequenceNumber(net) == totalSequenceNumber);
        assert(flower_getEndNumber(net) == 2*totalSequenceNumber);
        assert(eventTree_getEventNumber(flower_getEventTree(net)) == totalEventNumber);
        cactusDisk_destruct(netDisk);
    }

    return 0;
}
