#ifndef TREE_LIB_H
#define TREE_LIB_H

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align.h"
#include "buildtree.h"   
#include "cluster.h"   
#include "distancemat.h"
#include "tree.h"
#include "util.h"

#include "treelib_strings.h"

char * msa2tree(char **mfa, unsigned int num_species);

#endif
