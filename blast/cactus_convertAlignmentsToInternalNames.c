#include <getopt.h>
#include <errno.h>
#include "cactus.h"
#include "sonLib.h"
#include "pairwiseAlignment.h"
#include "bioioC.h"

static void usage(void)
{
    fprintf(stderr, "cactus_convertAlignmentsToInternalNames --cactusDisk cactusDisk alignmentsFile outputFile\n");
}

static void convertHeadersToNames(struct PairwiseAlignment *pA, stHash *headerToName)
{
    Name *name = NULL;
    if((name = stHash_search(headerToName, pA->contig1)) == NULL) {
        fprintf(stderr, "Error: sequence %s is not loaded into the cactus "
                "database\n", pA->contig1);
        exit(1);
    }
    pA->contig1 = cactusMisc_nameToString(*name);
    // Coordinates have to be shifted by 2 to keep compatibility with
    // cactus coordinates.
    pA->start1 += 2;
    pA->end1 += 2;
    if((name = stHash_search(headerToName, pA->contig2)) == NULL) {
        fprintf(stderr, "Error: sequence %s is not loaded into the cactus "
                "database\n", pA->contig2);
        exit(1);
    }
    pA->contig2 = cactusMisc_nameToString(*name);
    pA->start2 += 2;
    pA->end2 += 2;
}

int main(int argc, char *argv[])
{
    char *cactusDiskString = NULL;
    CactusDisk *cactusDisk;
    stKVDatabaseConf *kvDatabaseConf;
    stHash *headerToName;
    stList *flowers;
    Flower_EndIterator *endIt;
    End_InstanceIterator *capIt;
    End *end;
    FILE *alignmentFile;
    FILE *outputFile;
    struct option longopts[] = { {"cactusDisk", required_argument, NULL, 'c' },
                                 {0, 0, 0, 0} };
    int flag;
    while((flag = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
        switch(flag) {
        case 'c':
            cactusDiskString = stString_copy(optarg);
            break;
        case '?':
        default:
            usage();
            return 1;
        }
    }
    assert(argc == optind + 2);

    headerToName = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                     free, free);

    // Load a header->cactus ID map from the cactus DB
    if (cactusDiskString == NULL) {
        st_errAbort("--cactusDisk option must be provided");
    }
    kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    assert(stList_length(flowers) == 1);
    endIt = flower_getEndIterator(stList_get(flowers, 0));
    while((end = flower_getNextEnd(endIt)) != NULL) {
        capIt = end_getInstanceIterator(end);
        Cap *cap;
        while((cap = end_getNext(capIt)) != NULL) {
            const char *header;
            Name name;
            Name *heapName;
            if(!cap_getStrand(cap)) {
                cap = cap_getReverse(cap);
            }
            if(cap_getSide(cap)) {
                continue;
            }
            name = cap_getName(cap);
            header = sequence_getHeader(cap_getSequence(cap));
            heapName = st_malloc(sizeof(Name));
            *heapName = name;
            // The casts to void * are just to drop the const
            // qualifier. The functions won't modify the header.
            if(stHash_search(headerToName, (void *) header)) {
                // There is already a header -> cap name map, check
                // that it has the same name.
                Name *otherName = stHash_search(headerToName, (void *) header);
                fprintf(stderr, "Collision with header %s: name %" PRIi64
                        " otherName: %" PRIi64 "\n", header, name, *otherName);
                assert(*otherName == name);
            }
            stHash_insert(headerToName, (void *) header, heapName);
        }
        end_destructInstanceIterator(capIt);
    }

    // Scan over the given alignment file and convert the headers to
    // cactus Names.
    alignmentFile = fopen(argv[optind], "r");
    if(alignmentFile == NULL) {
    	st_errnoAbort("error opening alignment file %s", argv[optind]);
    }
    outputFile = fopen(argv[optind + 1], "w");
    if(outputFile == NULL) {
    	st_errnoAbort("error opening output file %s", argv[optind + 1]);
    }
    for(;;) {
        struct PairwiseAlignment *pA = cigarRead(alignmentFile);
        if(pA == NULL) {
            // Signals end of cigar file.
            break;
        }
        convertHeadersToNames(pA, headerToName);
        checkPairwiseAlignment(pA);
        cigarWrite(outputFile, pA, TRUE);
    }

    // Cleanup.
    flower_destructEndIterator(endIt);
    cactusDisk_destruct(cactusDisk);
}
