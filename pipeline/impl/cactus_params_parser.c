/*
 * Released under the MIT license, see LICENSE.txt
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <stdio.h>
#include <stdarg.h>

#include "sonLib.h"
#include "cactus.h"
#include "cactus_params_parser.h"

void cactusParams_destruct(CactusParams *p) {
    xmlFreeDoc(p->doc);
    free(p);
}

CactusParams *cactusParams_load(char *file_name) {
    CactusParams *p = st_calloc(1, sizeof(CactusParams));

    // Parse the XML file
    p->doc = xmlParseFile(file_name);

    if (p->doc == NULL ) {
        fprintf(stderr,"ERROR: Cactus XML params file: %s not parsed successfully. \n", file_name);
        free(p);
        return NULL;
    }

    p->root = xmlDocGetRootElement(p->doc);

    if (p->root == NULL) {
        fprintf(stderr,"ERROR: Empty Cactus params file: %s\n", file_name);
        cactusParams_destruct(p);
        return NULL;
    }

    if (xmlStrcmp(p->root->name, (const xmlChar *) "cactusWorkflowConfig")) {
        fprintf(stderr,"ERROR: Cactus XML params file: root node != cactusWorkflowConfig");
        cactusParams_destruct(p);
        return NULL;
    }

    // Set the current root node pointer to the actual root of the xml tree.
    p->cur = p->root;

    return p;
}

/*
 * Gets the child node with the given name of the current node.
 */
static xmlNodePtr get_child_node(xmlNodePtr cur, const char *child_node_name) {
    cur = cur->xmlChildrenNode;
    while (cur != NULL) {
        if ((!xmlStrcmp(cur->name, (const xmlChar *) child_node_name))) {
            return cur;
        }
        cur = cur->next;
    }
    return NULL;
}

/*
 * Get the descendant node with the given path.
 *
 * For example, get_descendant_node(cur, 2, "one", "two")
 * would get the great grandchild of cur with parent labelled "one" and
 * node label "two".
 */
static xmlNodePtr get_descendant_node(xmlNodePtr cur, int num, va_list args) {
    /* access all the arguments assigned to valist */
    for (int64_t i = 0; i < num; i++) {
        const char *node_name = va_arg(args, char *);
        cur = get_child_node(cur, node_name);
        if(cur == NULL) {
            fprintf(stderr,"ERROR: Cactus XML param node name %s not found\n", node_name);
            return NULL;
        }
    }

    return cur;
}

void cactusParams_set_root(CactusParams *p, int num, ...) {
    va_list args;
    va_start(args, num);
    p->cur = get_descendant_node(p->root, num, args);
    va_end(args);
}

static const char *cactusParams_get_string2(CactusParams *p, int num, va_list args) {
    xmlNodePtr c = get_descendant_node(p->cur, num-1, args);

    if(c == NULL) {
        st_errAbort("ERROR: Failed to get string from cactus XML");
    }
    const char *attribute_name = va_arg(args, char *);
    char *v =  (char *)xmlGetProp(c, (const xmlChar *)attribute_name);

    if(v == NULL) {
        st_errAbort("ERROR: Failed to get attribute: %s from cactus XML", attribute_name);
    }

    return v;
}

const char *cactusParams_get_string(CactusParams *p, int num, ...) {
    va_list args;
    va_start(args, num);
    const char *c = cactusParams_get_string2(p, num, args);
    va_end(args);

    return c;
}

int64_t cactusParams_get_int(CactusParams *p, int num, ...) {
    va_list args;
    va_start(args, num);

    const char *c = cactusParams_get_string2(p, num, args);
    int64_t j;
    int i = sscanf(c, "%" PRIi64 "", &j);
    assert(i == 1);

    va_end(args);
    return j;
}

double cactusParams_get_float(CactusParams *p, int num, ...) {
    va_list args;
    va_start(args, num);

    const char *c = cactusParams_get_string2(p, num, args);
    float j;
    int i = sscanf(c, "%f", &j);
    assert(i == 1);

    va_end(args);
    return j;
}
