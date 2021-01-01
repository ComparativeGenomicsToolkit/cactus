/*
 * Released under the MIT license, see LICENSE.txt
 */

#include "CuTest.h"
#include "sonLib.h"
#include "cactus_params_parser.h"

/*
 * TODOs:
 *
 * make branch
 *
 * setup parsing inputs
 *
 * refactor setup, caf, bar, reference to use XML object and be callable
 * from this code
 */

int main(int argc, char *argv[]) {
    char *params_file = "/Users/benedictpaten/CLionProjects/cactus/src/cactus/cactus_progressive_config.xml";

    CactusParams *p = cactusParams_load(params_file);

    cactusParams_destruct(p); // Cleanup

    return 0;
}
