/*
 * paf_trim: Trim bases from the prefix and suffix matches of a paf alignment
 *   (1) Load local alignments file (PAF)
 *   (2) Trim percentage (2% by default) of prefix and suffix matches.
 *   (3) Output local alignments file (PAF)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "paf.h"
#include <getopt.h>
#include <time.h>

static float trim_end_fraction = 0.02;

 void usage() {
     fprintf(stderr, "paf_trim [options], version 0.1\n");
     fprintf(stderr, "Trims the ends of a PAF file\n");
     fprintf(stderr, "-i --inputFile : Input paf file to invert. If not specified reads from stdin\n");
     fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
     fprintf(stderr, "-t --trimFraction : Fraction (from 0 to 1) of aligned bases"
                     "to trim from each end of the alignment (default:%f)\n", trim_end_fraction);
     fprintf(stderr, "-l --logLevel : Set the log level\n");
     fprintf(stderr, "-h --help : Print this help message\n");
 }

 int main(int argc, char *argv[]) {
     time_t startTime = time(NULL);

     /*
      * Arguments/options
      */
     char *logLevelString = NULL;
     char *inputFile = NULL;
     char *outputFile = NULL;

     ///////////////////////////////////////////////////////////////////////////
     // Parse the inputs
     ///////////////////////////////////////////////////////////////////////////

     while (1) {
         static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                 { "inputFile", required_argument, 0, 'i' },
                                                 { "outputFile", required_argument, 0, 'o' },
                                                 { "trimFraction", required_argument, 0, 't' },
                                                 { "help", no_argument, 0, 'h' },
                                                 { 0, 0, 0, 0 } };

         int option_index = 0;
         int64_t key = getopt_long(argc, argv, "l:i:o:ht:", long_options, &option_index);
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
             case 't':
                 trim_end_fraction = atof(optarg);
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
     st_logInfo("Trim fraction using : %f\n", trim_end_fraction);

     //////////////////////////////////////////////
     // Invert the paf
     //////////////////////////////////////////////

     FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
     FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

     Paf *paf;
     while((paf = paf_read(input)) != NULL) {
         paf_trim_end_fraction(paf, trim_end_fraction); // the invert routine
         paf_check(paf);
         paf_write(paf, output);
         paf_destruct(paf);
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

     st_logInfo("Paf invert is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

     //while(1);
     //assert(0);

     return 0;
 }
