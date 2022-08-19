/*
 * paf_view: Remove duplicates from a paf file
 *
 *  Released under the MIT license, see LICENSE.txt
 *
 * Overview:
 * (1) Load paf records
 * (2) For each paf in order of file add to a dictionary of sequence coordinates.
 * (3) If a record does not have the same coordinates as a previous entry omit print it, else omit it
*/

#include "paf.h"
#include "inc/paf.h"
#include <getopt.h>
#include <time.h>
#include "bioioC.h"

void usage() {
    fprintf(stderr, "paf_dedupe [options], version 0.1\n");
    fprintf(stderr, "Pretty print PAF alignments\n");
    fprintf(stderr, "-i --inputFile : Input paf file to invert. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-a --checkInverse : Also deduplicate alignments that are the same, but with query and target reversed\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

uint64_t paf_hash_key(const void *k) {
    Paf *p = (Paf *)k;
    uint64_t key = p->query_start + p->query_end + p->target_start + p->target_end;
    // Use the hash from <https://stackoverflow.com/a/12996028>
    // We can't just -ull these until C++11, if we're in C++
    key = (key ^ (key >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    key = (key ^ (key >> 27)) * UINT64_C(0x94d049bb133111eb);
    return key ^ (key >> 31);
}

int paf_equal_key(const void *k, const void *k2) {
    Paf *p = (Paf *)k, *p2 = (Paf *)k2;
    return strcmp(p->query_name, p2->query_name) == 0
            && strcmp(p->target_name, p2->target_name) == 0
            && p->same_strand == p2->same_strand
            && p->target_start == p2->target_start
            && p->target_end == p2->target_end
            && p->query_start == p2->query_start
            && p->query_end == p2->query_end; // Pafs are equal if they have the same query and target coordinates.
}

int main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    bool check_inverse=0;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "checkInverse", required_argument, 0, 'a' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:ha", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                inputFile = optarg;
                break;
            case 'o':
                outputFile = optarg;
                break;
            case 'a':
                check_inverse = 1;
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);

    //////////////////////////////////////////////
    // Remove duplicate paf records
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");
    stHash *pafs = stHash_construct3(paf_hash_key, paf_equal_key, NULL, (void (*)(void *))paf_destruct);
    Paf *paf;
    while((paf = paf_read(input)) != NULL) {
        // Get the query sequence
        Paf *pPaf = stHash_search(pafs, paf);
        if(check_inverse && pPaf == NULL) { // In case we want to check if we already output the inverse
            paf_invert(paf); // Invert
            pPaf = stHash_search(pafs, paf); // See if inverse is in there
            paf_invert(paf); // Invert back
            paf_check(paf); // Check is okay
        }
        if(pPaf == NULL) {  // If duplicate is not already in there
            stHash_insert(pafs, paf, paf);  // Add the paf
            paf_write(paf, output); // Write the paf to the output
        }
        else {
            // If debug output report info on dupe
            // Print the alignment levels
            if(st_getLogLevel() >= debug) {
                st_logDebug("Got duplicate pafs:\n");
                char *paf_string = paf_print(paf);
                char *paf_string_2 = paf_print(pPaf);
                st_logDebug("\t\tdupe (1) - : %s\n", paf_string);
                st_logDebug("\t\tdupe (2) - : %s\n", paf_string_2);
                free(paf_string); free(paf_string_2);
            }
            paf_destruct(paf); // Just clean it up and don't output it
        }
    }

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(inputFile != NULL) {
        fclose(input);
    }
    if(outputFile != NULL) {
        fclose(output);
    }
    stHash_destruct(pafs);

    st_logInfo("Paf dedupe is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    //while(1);
    //assert(0);

    return 0;
}

