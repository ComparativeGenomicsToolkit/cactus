/*  Last edited: Aug 24 14:28 1999 (klh) */
/**********************************************************************
 ** FILE: buildtree.h
 ** NOTES:
 **  Contains functions for building trees from distance matrices
 **  (and vice versa)
 **********************************************************************/

#ifndef _BUILDTREE
#define _BUILDTREE

#include <float.h>
#include "util.h"
#include "tree.h"
#include "cluster.h"

/******************* structure definitions ****************************/


/********************** function prototypes ***************************/


/**********************************************************************
 FUNCTION: export_distances_buildtree
 DESCRIPTION: 
   Returns the distance matrix induced from the given tree, i.e. by summing
   the branch paths between two nodes to obtain their distance
 ARGS: 
   A Tree
   A DistanceMatrix
 RETURNS: 

 NOTES: 
   This function does not create the memory for the distance matrix, it
   merely fills in the given matrix
 **********************************************************************/
void export_distances_buildtree( struct Tree *, struct DistanceMatrix *);


/**********************************************************************
 FUNCTION: find_path_buildtree
 DESCRIPTION: 
   Given a leaf node, recursively calculates the branch-length distance
   from the node to all other leaves in the tree, and places it in the
   appropriate part of the DistanceMatrix
 ARGS: 
   unsigned int (the node number from which we are finding all distances)
   Tnode (the current node under consideration)
   DistanceMatrix
   Boolean array (to store which nodes have already been considered)
   unsigned int (the size of this boolean array)
 RETURNS: 
 NOTES: 
 **********************************************************************/
void find_path_buildtree( unsigned int,
			  struct Tnode *,
			  struct DistanceMatrix *,
			  Distance,
			  unsigned int *);


/**********************************************************************
 FUNCTION: leaf_find_buildtree
 DESCRIPTION: 
   Finds the leaf nodes descended from the given interior node, and
   calculates the distance from these nodes to all other nodes, by
   way of another function call
 ARGS: 
   Tnode
   DistanceMatrix
   Boolean array (to store which nodes have already been considered)
   unsigned int (the size of this boolean array)
 RETURNS: 
 NOTES: 
 **********************************************************************/
void leaf_find_buildtree( struct Tnode *, 
			  struct DistanceMatrix *, 
			  unsigned int *,
			  unsigned int);



/**********************************************************************
 FUNCTION: neighbour_joining_buildtree
 DESCRIPTION: 
   Returns a phylogenetic tree of the sequences in the 
   given alignment, using Saitou and Nei's neighbour-joining 
   algorithm
 ARGS: 
    A ClusterGroup pointer (cluster.h)
    Boolean, for whether to calc information needed for later bootstrapping
 RETURNS:
    A Tree (trees.h)
 NOTES: The function allocates all the memory necessary for the tree.
   The caller should call free_tree (tree.h) to free this memory when
   the tree is no longer needed
 **********************************************************************/
struct Tree *neighbour_joining_buildtree( struct ClusterGroup *,
					  unsigned int);


/**********************************************************************
 FUNCTION: UPGMA_buildtree
 DESCRIPTION: 
   Returns a phylogenetic tree of the sequences in the 
   given alignment, using the Unweighted Pair-Group method based on
   Arithmentic Averages (UPGMA)
 ARGS: 
   A ClusterGroup pointer (cluster.h)
   Boolean, for whether to calc information needed for later bootstrapping
 RETURNS: 
   struct Tree (trees.h)
 NOTES: The function allocates all the memory necessary for the tree.
   The caller should call free_Tnode (tree.h) to free this memory when
   the tree is no longer needed

   This algorithm produces a rooted tree; hence the returned Tree
   will have the root in child[0]; child[1] and child[2] (used to
   represent the trichotomy of unrooted trees) will be NULL 
 **********************************************************************/
struct Tree *UPGMA_buildtree( struct ClusterGroup *, unsigned int);


#endif
