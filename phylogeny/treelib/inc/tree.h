/*  Last edited: Jun 21 12:03 2001 (klh) */
/**********************************************************************
 ** FILE: tree.h
 ** NOTES:
 **  Functions and types for the manipulation of trees
 **********************************************************************/

#ifndef _TREE
#define _TREE

#include <stdio.h>

#include "util.h"
#include "sequence.h"
#include "cluster.h"

#define LINE_BUFFER_SIZE 40


/******************* structure definitions ****************************/

/******************
   This structure represents a tree-node. Everything is self-explanatory,
   except that much of the tree-building code relies upon a correspondence
   between node-numbers of leaf nodes, and sequence/cluster numbers. This 
   is particulary important when re-constructing distance matrices from 
   trees, and for bootstrapping trees. This may be troublesome if a tree is
   read in from a file... 
******************/


struct Tnode {
  struct Tnode *left;
  struct Tnode *right;
  struct Tnode *parent;
  double distance;
  unsigned int nodenumber;
  struct Cluster *clust;
  unsigned int bootstrap;
  unsigned int *child_ids;
};


/***********************
 This structure represents an unrooted binary tree of three or more 
   nodes; Two node trees, although trivial, can be modelled by having
   one of the three children as null. A rooted tree can be represented
   by a Tnode only, or by having two of the children of the trichotomy
   as NULL.
*************************/

struct Tree {
  struct Tnode *child[3];
  unsigned int numnodes;
};


/********************** function prototypes ***************************/

/********************************************************************* 
 FUNCTION: assign_nodesnumbers_Tnode
 DESCRIPTION: 
   This function assigns node numbers to the internal nodes of the 
   given Tnode, starting from the given, and returns the next free 
   nodenumner 
 RETURNS: unsigned int
 ARGS: 
   struct Tnode * (tree)
   unsigned int (starting node number)
 NOTES: 
 *********************************************************************/
unsigned int assign_nodenumbers_Tnode( struct Tnode *, unsigned int);


/********************************************************************* 
 FUNCTION: clone_Tnode
 DESCRIPTION: 
   This function makes a complete copy of the tree rooted at the given
   node and returns it
 RETURNS: Tnode *
 ARGS: 
   struct Tnode *
 NOTES: 
 *********************************************************************/
struct Tnode *clone_Tnode( struct Tnode *);

/********************************************************************* 
 FUNCTION: clone_Tree
 DESCRIPTION: 
   This function makes a complete copy of the Tree and returns it
 RETURNS: struct Tree *
 ARGS: 
   struct Tree *
 NOTES: 
 *********************************************************************/
struct Tree *clone_Tree( struct Tree *);


/********************************************************************* 
 FUNCTION: compare_to_bootstrap_sample_Tnode
 DESCRIPTION: 
   Updates the bootstrap values of the given Tnode, according to
   the topology of the given sample Tnode
 RETURNS:
 ARGS:
   Destination Tnode
   Sample Tnode
   The number of leaf nodes in each tree
   Boolean for whether the tree is binary or not
 NOTES:
   This function assumes that the given trees have been created by 
   calling either neighbourjoin_buildtree or UPGMA_buildtree with 
   the bootstrap boolean arguement set to true. A fatal error results
   if this is not the case. This is due to the fact that the information
   needed for the tree comparisons is determined when the trees are 
   constructed
 *********************************************************************/
void compare_to_bootstrap_sample_Tnode( struct Tnode *, struct Tnode *,
					unsigned int, unsigned int);

/********************************************************************* 
 FUNCTION: empty_Tree
 DESCRIPTION: 
   Creates and returns a tree with null nodes
 RETURNS: struct Tree * (trees.h)
 ARGS: 
 NOTES:
 *********************************************************************/
struct Tree *empty_Tree( void );

/********************************************************************* 
 FUNCTION: free_Tnode
 DESCRIPTION: 
   This function releases the memory used by this Tnode and all of its
   children
 RETURNS: A null pointer
 ARGS: 
   struct Tnode *
 NOTES: 
 *********************************************************************/
void *free_Tnode( struct Tnode *);

/********************************************************************* 
 FUNCTION: free_Tree
 DESCRIPTION: 
   This function releases the memory used br the Tnode chain in the
   given Tree
 RETURNS: A null pointer
 ARGS: 
   struct Tree
 NOTES: 
 *********************************************************************/
void *free_Tree( struct Tree *);

/********************************************************************* 
 FUNCTION: get_root_Tnode
 DESCRIPTION: This function takes a three node tree (struct Tree), and
    returns a Tnode as the 'root' of the tree. The root is inserted
    somewhat arbitraily between the three top-level nodes
 RETURNS: struct Tnode *
 ARGS: struct Tree
 NOTES: 
   The nodes of the given tree are cloned for use. This means
   that the old tree is still available and safe on return, and must
   be freed when finished with. The rooted tree returned by this 
   function must be freed by a call to free_Tnode_tree
 *********************************************************************/
struct Tree *get_root_Tnode( struct Tree * );

/********************************************************************* 
 FUNCTION: new_interior_Tnode (unsigned int)
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new interior (non-leaf) tree node, filling it with what it 
   knows, and returning it.
 RETURNS: struct Tnode * (trees.h)
 ARGS: 
   A integer to hold the node number.
 NOTES:
   This function s for creating internal nodes, which have no string
   id but just a number identifier
 *********************************************************************/
struct Tnode *new_interior_Tnode( unsigned int );

/********************************************************************* 
 FUNCTION: new_leaf_Tnode (int, char *)
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new tree node, filling it with what it knows, and 
   returning it.
 RETURNS: struct Tnode * (trees.h)
 ARGS: 
   A integer to hold the sequence number associated with the node
   The name of the node
 NOTES: 
   This function is for creating leaf nodes, which have a name.
 *********************************************************************/
struct Tnode *new_leaf_Tnode( unsigned int, struct Cluster *);

/********************************************************************* 
 FUNCTION: read_newhampshire_Tnode
 DESCRIPTION: 
   Construct a Tnode from the given file handle
 RETURNS: The total number of nodes constructed as a result of the call
 ARGS: 
   File handle (assumed in New Hampshire format)
 NOTES:
   Numbering of internal nodes is abandoned in favour of giving leaf
   nodes numbers from 0 to the number of leaves in the tree. This makes
   them handy indices into a distance matrix for example
 *********************************************************************/
unsigned int read_newhampshire_Tnode( FILE *, struct Tnode **, 
				      struct Tnode *, unsigned int);

/********************************************************************* 
 FUNCTION: read_newhampshire_Tree
 DESCRIPTION: 
   Constructs a tree from the given file handle
 RETURNS:  Tree *
 ARGS: A handle to the file (assumed in New Hamshire format)
 NOTES:
   The convention that the numbering of the interior nodes starts
   after all leaf nodes have been numbered is is violated with this
   method; the numbering is performed in a bottom-up, left-to-right 
   fashion, so that the nodes in the left subtree all have number-ids
   strictly less than those in the right subtree

 *********************************************************************/
struct Tree *read_newhampshire_Tree( FILE *);

/********************************************************************* 
 FUNCTION: scale_bootstraps_Tnode
 DESCRIPTION: 
   This function traverses the given Tnode, dividing the bootstrap
   values by the given number
 RETURNS: 
 ARGS: 
   Tnode
   A the number of bootstrap iterations that were performed
 NOTES:

 *********************************************************************/
void scale_bootstraps_Tnode( struct Tnode *, unsigned int);


/********************************************************************* 
 FUNCTION: scale_bootstraps_Tree
 DESCRIPTION: 
   This function traverses the gieven tree, dividing the bootstrap
   values by the given number
 RETURNS: 
 ARGS: 
   Tree
   A the number of bootstrap iterations that were performed
 NOTES:

 *********************************************************************/
void scale_bootstraps_Tree( struct Tree *, unsigned int);


/********************************************************************* 
 FUNCTION: update_bootstraps_Tree
 DESCRIPTION: 
   Updates the bootstrap values of the destination Tree, according to
   the topology of the given sample tree
 RETURNS:
 ARGS:
   Destination tree
   Sample tree
   The number of leaf nodes in the trees
 NOTES:
 *********************************************************************/
void update_bootstraps_Tree( struct Tree *, struct Tree *, unsigned int);

/********************************************************************* 
 FUNCTION: update_bootstraps_Tnode
 DESCRIPTION: 
   Updates the bootstrap values of the destination Tnode, according to
   the topology of the given sample Tnode
 RETURNS:
 ARGS:
   Destination Tnode
   Sample Tnode
   The number of leaf nodes in the tree
   Boolean for whether the trees are binary or not
 NOTES:
 *********************************************************************/
void update_bootstraps_Tnode( struct Tnode *, struct Tnode *, 
			      unsigned int, unsigned int);

/********************************************************************* 
 FUNCTION: write_clustering_data_Tnode
 DESCRIPTION: 
   This routine prints a text description of the clustering details
   of the given Tnode. It was written for the old implementation,
   where leaves were named "leaf_1", "leaf_2" etc, and it would
   often be the case that each leaf would contain several sequences.
   With the current implementation, this function is not used; a 
   new-hampshire output of each cluster is printed in-situ, preserving
   sequence names from the original alignment
 RETURNS:
 ARGS: 
   File handle
   TNode *
 NOTES:
 *********************************************************************/
void write_clustering_data_Tnode( FILE *, struct Tnode *);


/********************************************************************* 
 FUNCTION: write_debug_Tnode
 DESCRIPTION: 
   Writes the given Tnode to the give file handle in 'debug' format
 RETURNS:
 ARGS: 
   File handle
   TNode *
   Integer offset
 NOTES:
 *********************************************************************/
void write_debug_Tnode( FILE *, struct Tnode *, unsigned int);

/********************************************************************* 
 FUNCTION: write_debug_Tree
 DESCRIPTION: 
   prints the given Tree in a format suitable for debugging
 RETURNS:
 ARGS: 
   File handle
   Tree *
 NOTES:
 *********************************************************************/
void write_debug_Tree( FILE *, struct Tree *);

/********************************************************************* 
 FUNCTION: write_MUL_flattened_Tnode
 DESCRIPTION: 
   Prints the given tree as a MUL format alignment, with the sequences
   in 'tree order'
 RETURNS:
 ARGS: 
   File handle
   TNode *
 NOTES:
 *********************************************************************/
void write_MUL_flattened_Tnode( FILE *, struct Tnode *);

/********************************************************************* 
 FUNCTION: write_MUL_flattened_Tree
 DESCRIPTION: 
   Prints the given tree as a MUL format sequence alignment, with
   the sequences in 'tree order'
 RETURNS:
 ARGS: 
   File handle
   Tree *
 NOTES:
 *********************************************************************/
void write_MUL_flattened_Tree( FILE *, struct Tree *);


/********************************************************************* 
 FUNCTION: write_newhampshire_Tnode
 DESCRIPTION: 
   prints the given Tree in 'New Hampshire' text format to the given 
   file handle
 RETURNS:
 ARGS: 
   File handle
   TNode *
 NOTES:
 *********************************************************************/
void write_newhampshire_Tnode( FILE *, struct Tnode *, unsigned int);

/********************************************************************* 
 FUNCTION: write_newhampshire_Tree
 DESCRIPTION: 
   prints the given Tnode in 'New Hampshire' text format to the given 
   file handle
 RETURNS:
 ARGS: 
   File handle
   TNode *
 NOTES:
 *********************************************************************/
void write_newhampshire_Tree( FILE *, struct Tree *, unsigned int);



#endif
