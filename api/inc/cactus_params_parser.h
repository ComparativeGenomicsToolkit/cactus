/*
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef ST_CACTUS_PARAMS_PARSER_H_
#define ST_CACTUS_PARAMS_PARSER_H_

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>

/*
 * Cactus parameters object.
 */
typedef struct _cactusParams {
    xmlDocPtr doc; // The underlying xml document representing the parameters
    xmlNodePtr root; // The root node
    xmlNodePtr cur; // The node of the xml tree we search from to retrieve parameters.
    // can be set by cactusParams_set_root(CactusParams *p, int, ...), by default is set
    // to the root of the tree.
} CactusParams;

/*
 * Cleanup the CactusParams.
 */
void cactusParams_destruct(CactusParams *p);

/*
 * Load the CactusParams.
 */
CactusParams *cactusParams_load(char *file_name);

/*
 * Set the root node of the params tree.
 * e.g. cactusParams_set_root(p, 2, "caf", "divergence") would set the root
 * to the cactusWorkflowConfig->caf->divergence node.
 */
void cactusParams_set_root(CactusParams *p, int num, ...);

/*
 * Get a string parameter.
 */
char *cactusParams_get_string(CactusParams *p, int, ...);

/*
 * Get an integer parameter.
 */
int64_t cactusParams_get_int(CactusParams *p, int, ...);

/*
 * Get a float parameter.
 */
double cactusParams_get_float(CactusParams *p, int, ...);

/*
 * Get a sequence of integers
 */
int64_t *cactusParams_get_ints(CactusParams *p, int64_t *length, int, ...);

#endif /* ST_CACTUS_PARAMS_PARSER_H_ */
