#include <getopt.h>
#include <errno.h>
#include <error.h>
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
    Name *name;
    if((name = stHash_search(headerToName, pA->contig1)) == NULL) {
        fprintf(stderr, "Error: sequence %s is not loaded into the cactus "
                "database\n", pA->contig1);
        exit(1);
    }
    pA->contig1 = cactusMisc_nameToString(*name);
    if((name = stHash_search(headerToName, pA->contig2)) == NULL) {
        fprintf(stderr, "Error: sequence %s is not loaded into the cactus "
                "database\n", pA->contig2);
        exit(1);
    }
    pA->contig2 = cactusMisc_nameToString(*name);
}

int main(int argc, char *argv[])
{
    char *cactusDiskString;
    CactusDisk *cactusDisk;
    stKVDatabaseConf *kvDatabaseConf;
    stHash *headerToName;
    stList *flowers;
    Flower_SequenceIterator *flowerIterator;
    Sequence *sequence;
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
    kvDatabaseConf = stKVDatabaseConf_constructFromString(cactusDiskString);
    cactusDisk = cactusDisk_construct(kvDatabaseConf, 0);
    flowers = flowerWriter_parseFlowersFromStdin(cactusDisk);
    assert(stList_length(flowers) == 1);
    flowerIterator = flower_getSequenceIterator(stList_get(flowers, 0));
    while((sequence = flower_getNextSequence(flowerIterator)) != NULL) {
        const char *header;
        Name name;
        Name *heapName;
        header = sequence_getHeader(sequence);
        name = sequence_getName(sequence);
        heapName = st_malloc(sizeof(Name));
        *heapName = name;
        // The casts to void * are just to drop the const
        // qualifier. The functions won't modify the header.
        assert(stHash_search(headerToName, (void *) header) == NULL);
        stHash_insert(headerToName, (void *) header, heapName);
    }

    // Scan over the given alignment file and convert the headers to
    // cactus Names.
    alignmentFile = fopen(argv[optind], "r");
    if(alignmentFile == NULL) {
        error(1, errno, "error opening alignment file %s", argv[optind]);
    }
    outputFile = fopen(argv[optind + 1], "w");
    if(outputFile == NULL) {
        error(1, errno, "error opening output file %s", argv[optind + 1]);
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
    flower_destructSequenceIterator(flowerIterator);
    cactusDisk_destruct(cactusDisk);
}
