/*  Last edited: Mar 18 11:55 2002 (klh) */
/**********************************************************************
 ** FILE: tree.c
 ** NOTES:
 **  Functions for the manipulation of trees
 **********************************************************************/

#include "tree.h"



/********************************************************************* 
 FUNCTION: assign_nodenumbers_Tnode
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
unsigned int assign_nodenumbers_Tnode( struct Tnode *node, unsigned int start) {

  if (node != NULL) {
    start += assign_nodenumbers_Tnode( node->left, start );
    start += assign_nodenumbers_Tnode( node->right, start );

    node->nodenumber = start++;
  }

  return start; 
}



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
struct Tnode *clone_Tnode( struct Tnode *source) {
  struct Tnode *dest = NULL;

  if (source != NULL) {
    dest = (struct Tnode *) malloc_util( sizeof( struct Tnode ) );
    
    dest->distance = source->distance;
    dest->nodenumber = source->nodenumber;
    dest->bootstrap = source->bootstrap;
    dest->clust = clone_Cluster( source->clust );
    if ( (dest->left = clone_Tnode( source->left )) != NULL)
      dest->left->parent = dest;
    if ( (dest->right = clone_Tnode( source->right )) != NULL)
      dest->right->parent = dest;
  }

  return dest; 
}



/********************************************************************* 
 FUNCTION: clone_Tree
 DESCRIPTION: 
   This function makes a complete copy of the Tree and returns it
 RETURNS: struct Tree *
 ARGS: 
   struct Tree *
 NOTES: 
 *********************************************************************/
struct Tree *clone_Tree( struct Tree *source) {
  struct Tree *dest = NULL;

  if ( source != NULL ) {
    dest = (struct Tree *) malloc_util( sizeof( struct Tree ) );
    dest->numnodes = source->numnodes;
    dest->child[0] = clone_Tnode( source->child[0] );
    dest->child[1] = clone_Tnode( source->child[1] );
    dest->child[2] = clone_Tnode( source->child[2] );
  }

  return dest;

}




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
void compare_to_bootstrap_sample_Tnode( struct Tnode *dest, 
					struct Tnode *sample,
					unsigned int numleaves,
					unsigned int is_binary) {
  unsigned int matchcounter, i;

  if (dest != NULL) {
    if (dest->child_ids != NULL) {

      /* we have to explore every non-terminal node of sample to see
	 if the bit pattern is the same, or a mirror image. If yes, then
	 this sub-tree exists in sample. The mirror image is to take care
	 of the fact that trichotomius trees may be isomorphic but 
	 'centred' at a different node.
	 
	 For non-trichotomous (i.e. rooted binary) trees, it is erroneous
	 to allow these 'mirror image' cases. It is easy to found out what
	 sort of trees we have by examining the child fields in one of the
	 trees. If only one of them is non-null, then we have a rooted 
	 binary tree.
      */

      if (sample != NULL) {
	if (sample->child_ids != NULL) {
	  matchcounter = 0;
	  for(i=0; i < numleaves; i++) {
	    if (dest->child_ids[i] == sample->child_ids[i]) {
	      matchcounter++;
	    }
	  }
	  if ((matchcounter == numleaves) || (matchcounter == 0 && ! is_binary)) {
	    dest->bootstrap++;
	    /* printf ("Incrementing node %d...\n", dest->nodenumber ); */
	  }
	  
	  compare_to_bootstrap_sample_Tnode( dest, sample->left, numleaves, is_binary );
	  compare_to_bootstrap_sample_Tnode( dest, sample->right, numleaves, is_binary );
	}
      }
    }
  }
}




/********************************************************************* 
 FUNCTION: empty_Tree
 DESCRIPTION: 
   Creates and returns a tree with null nodes
 RETURNS: struct Tree * (trees.h)
 ARGS: 
 NOTES:
 *********************************************************************/

struct Tree *empty_Tree( void ) {
  struct Tree *ret;

  ret = (struct Tree *) malloc_util(sizeof( struct Tree ));
  ret->child[0] = NULL;
  ret->child[1] = NULL;
  ret->child[2] = NULL;
  ret->numnodes = 0;

  return ret;
}



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

void *free_Tnode( struct Tnode *tn ) {
  if ( tn != NULL ) {
    tn->clust = free_Cluster( tn->clust );
    tn->left = free_Tnode( tn->left );
    tn->right = free_Tnode( tn->right );
    if (tn->child_ids != NULL) {
      tn->child_ids = free_util( tn->child_ids );
    }
    tn = free_util( tn );
  }
  return tn; 

}



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

void *free_Tree( struct Tree *t) {
  if ( t != NULL ) {
    t->child[0] = free_Tnode( t->child[0] );
    t->child[1] = free_Tnode( t->child[1] );
    t->child[2] = free_Tnode( t->child[2] );
    t = free_util( t );
  }
  return t;
}



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
   function must be freed by a call to free_Tnode
 *********************************************************************/

struct Tree *get_root_Tnode( struct Tree *source ) {
  struct Tnode *focal, *root, *children[3];
  struct Tree *ret;
  unsigned int rootleft, focalleft, focalright;
  Distance maxdist;

  /***** Method **************
     0. Clone the given tree
     1. Create the imaginary node between the three nodes in the Tree
     2. Create a root node;
  **************************/

  children[0] = clone_Tnode( source->child[0] );
  children[1] = clone_Tnode( source->child[1] );
  children[2] = clone_Tnode( source->child[2] );

  focal = new_interior_Tnode( source->numnodes );
  root = new_interior_Tnode( source->numnodes + 1 );

  /* arbitrarity choose halfway along the longest branch between
     the three nodes as the position for the root */

  maxdist = children[0]->distance;
  rootleft = 0;
  focalleft = 1;
  focalright = 2;
  if (children[1]->distance > maxdist) {
    rootleft = 1;
    focalleft = 0;
    focalright = 2;
  }   
  if (children[2]->distance > maxdist) {
    rootleft = 2;
    focalleft = 0;
    focalright = 1;
  }

  /* sort out distances; root node has zero distances */

  children[rootleft]->distance = children[rootleft]->distance * 0.5;
  focal->distance = children[rootleft]->distance;

  /* now sort out the links */

  root->right = focal;
  root->left = children[rootleft];
  focal->parent = children[rootleft]->parent = root;

  focal->left = children[focalleft];
  focal->right = children[focalright];
  children[focalleft]->parent = children[focalright] = focal;

  ret = empty_Tree();
  ret->child[0] = root;
  ret->numnodes = source->numnodes + 2;

  return ret;
}



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

struct Tnode *new_interior_Tnode( unsigned int label ) {
  struct Tnode *newNode;

  newNode = (struct Tnode *) malloc_util(sizeof(struct Tnode));
  newNode->left = NULL;
  newNode->right = NULL;
  newNode->parent = NULL;
  newNode->distance = 0.0;
  newNode->nodenumber = label;
  newNode->clust = NULL;
  newNode->bootstrap = 0;
  newNode->child_ids = NULL;

  return newNode;
}





/********************************************************************* 
 FUNCTION: new_leaf_Tnode(int, char *)
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

struct Tnode *new_leaf_Tnode(unsigned int label, struct Cluster *given) {
  struct Tnode *newNode;

  newNode = new_interior_Tnode( label );
  newNode->clust = given;

  return newNode;
}



/********************************************************************* 
 FUNCTION: read_newhampshire_Tnode
 DESCRIPTION: 
   Construct a Tnode from the given file handle
 RETURNS: The total number of nodes constructed as a result of the call
 ARGS: 
   File handle (assumed in New Hampshire format)
 NOTES: Due to the fact that the neighbourjoining implementation presented
   in these modules has 'clusters' of sequences at the leaf nodes;
   This means that the written tree will no longer be binary, so this
   method wil no longer read in the trees produced by 
   write_newhampshire_Tree. I should look into ways around this
 *********************************************************************/

unsigned int read_newhampshire_Tnode( FILE *handle, 
				      struct Tnode **nodeptrptr,
				      struct Tnode *parent,
				      unsigned int nodecounter ) {
  
  char c;
  unsigned int index;
  struct Sequence *newseq;
  double distance;
  

  if ( !fscanf( handle, "%1s", &c) )
      fatal_util("fscanf failed");

  if ( c == '(' ) {
    /* we do not know the node number until we have parsed the children,
       so give it the value zero for now */

    *nodeptrptr = new_interior_Tnode( 0 );
  
    nodecounter += read_newhampshire_Tnode( handle, 
					    &((*nodeptrptr)->left),
					    *nodeptrptr,
					    nodecounter );

    if (!fscanf( handle, "%1s", &c)) /* should be , */
        fatal_util("fscanf failed");
    if ( c != ',')
      fatal_util( "Parse error: ',' expected");

    nodecounter += read_newhampshire_Tnode( handle, 
					    &((*nodeptrptr)->right),
					    *nodeptrptr,
					    nodecounter );
  
    /* (*nodeptrptr)->nodenumber = nodecounter++; */

    if(!fscanf( handle, "%1s", &c))  /* should be ) */
        fatal_util("fscanf failed");
    if ( c != ')')
      fatal_util( "Parse error: ')' expected");
    if (!fscanf( handle, "%1s", &c))  /* should be : */
        fatal_util("fscanf failed");
    if ( c != ':')
      fatal_util( "Parse error: ':' expected");
    if (!fscanf( handle, "%lf", &distance ))
      fatal_util( "Parse error: floating point number expexted");

  }
  else { /* Must be the first char of an identifier */
    ungetc( c, handle );

    newseq = empty_Sequence();
    newseq->name = (char *) malloc_util( MAX_NAME_LENGTH * sizeof( char ) );

    for( index=0; (newseq->name[index] = fgetc( handle )) != ':'; index++);
    newseq->name[index] = '\0';

    if (!fscanf( handle, "%lf", &distance ))
      fatal_util( "Parse error: floating point number expexted");
    
    *nodeptrptr = new_leaf_Tnode( nodecounter++, single_Sequence_Cluster( newseq ) );
  }

  (*nodeptrptr)->parent = parent;
  (*nodeptrptr)->distance = distance;

  return nodecounter;
}



/********************************************************************* 
 FUNCTION: read_newhampshire_Tree
 DESCRIPTION: 
   Constructs a tree from the given file handle
 RETURNS:  Tree *
 ARGS: A handle to the file (assumed in New Hamshire format)
 NOTES:
   Numbering of internal nodes is abandoned in favour of giving leaf
   nodes numbers from 0 to the number of leaves in the tree. This makes
   them handy indices into a distance matrix for example
 *********************************************************************/

struct Tree *read_newhampshire_Tree( FILE *handle ) {
  unsigned int numnodes = 0; 
  char c;

  struct Tree *thetree = empty_Tree();

  if(!fscanf( handle, "("))
      fatal_util("fscanf failed");

  numnodes += read_newhampshire_Tnode( handle, 
				       &thetree->child[0],
				       NULL,
				       numnodes);
  
  if(!fscanf( handle, "%1s", &c )) /* should be , */
      fatal_util("fscanf failed");
  if ( c != ',')
    fatal_util( "Parse error: ',' expected");

  numnodes += read_newhampshire_Tnode( handle, 
				       &thetree->child[1],
				       NULL,
				       numnodes);

  if(!fscanf( handle, "%1s", &c )) /* should be , */
      fatal_util("fscanf failed");
  if ( c != ',')
    fatal_util( "Parse error: ',' expected");

  numnodes += read_newhampshire_Tnode( handle,
				       &thetree->child[2],
				       NULL,
				       numnodes);

  if(!fscanf( handle, "%1s", &c))
      fatal_util("fscanf failed");
  if ( c != ')')
    fatal_util( "Parse error: ')' expected");
  if(!fscanf( handle, "%1s", &c))
      fatal_util("fscanf failed");
  if ( c != ';')
    fatal_util( "Parse error: ';' expected");

  numnodes += assign_nodenumbers_Tnode( thetree->child[0], numnodes );
  numnodes += assign_nodenumbers_Tnode( thetree->child[1], numnodes );
  numnodes += assign_nodenumbers_Tnode( thetree->child[2], numnodes );
  thetree->numnodes = numnodes;

  return thetree;

}



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
void scale_bootstraps_Tnode( struct Tnode *node, unsigned int iters) {

  if (node != NULL) {
    node->bootstrap = (int) (((double) node->bootstrap / (double) iters) * 100.0);
    scale_bootstraps_Tnode( node->left, iters);
    scale_bootstraps_Tnode( node->right, iters);
  }

}


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
void scale_bootstraps_Tree( struct Tree *thetree, unsigned int iters) {

  if (thetree != NULL) {
    /* The first node will always be defined...(he says) */
    scale_bootstraps_Tnode( thetree->child[0], iters );
    if (thetree->child[1] != NULL) {
      scale_bootstraps_Tnode( thetree->child[1], iters );
      if (thetree->child[2] != NULL) {
	scale_bootstraps_Tnode( thetree->child[2], iters );
      }
    }
  }
}



/********************************************************************* 
 FUNCTION: update_bootstraps_Tree
 DESCRIPTION: 
   Updates the bootstrap values of the destination tree, according to
   the topology of the given sample tree
 RETURNS:
 ARGS:
   Destination tree
   Sample tree
   the number of leaf nodes in the tree
 NOTES:
   This function assumes that the given trees have been created by 
   calling either neighbourjoin_buildtree or UPGMA_buildtree with 
   the bootstrap boolean arguement set to true. A fatal error results
   if this is not the case

   Another thing to note is that the method uses node numbers for
   comparisons. This means that both ethe sample and reference trees
   must have been built from the same initial list of nodes, which
   is fine for bootstrapping, because the leaf nodes are stored in
   the order in which appear in the alignment
 *********************************************************************/
void update_bootstraps_Tree( struct Tree *dest, struct Tree *sample, 
				  unsigned int numleaves) {
  unsigned int is_binary,i,j;

  is_binary = ( dest->child[1] == NULL && dest->child[2] == NULL )?1:0; 

  for (i=0; i < 3; i++) {
    for (j=0; j < 3; j++) {
      update_bootstraps_Tnode( dest->child[i], 
				    sample->child[j],
				    numleaves,
				    is_binary );
    }
  }
}



/********************************************************************* 
 FUNCTION: update_bootstraps_Tnode
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
   if this is not the case
 *********************************************************************/
void update_bootstraps_Tnode( struct Tnode *dest, 
			      struct Tnode *sample,
			      unsigned int numleaves,
			      unsigned int is_binary) {


  compare_to_bootstrap_sample_Tnode( dest, sample, numleaves, is_binary);
  if (dest != NULL) {
    update_bootstraps_Tnode( dest->left, sample, numleaves, is_binary);
    update_bootstraps_Tnode( dest->right, sample, numleaves, is_binary);
  }
}




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
void write_clustering_data_Tnode( FILE *handle, struct Tnode *node) {
  unsigned int i;

  if (node != NULL) {
    if (node->left == NULL && node->right == NULL && node->clust != NULL) {
      /* we have a leaf node */
      if (node->clust->clustersize == 0 || node->clust->members == NULL) 
	fatal_util("Fatal Error: encountered a leaf node with no cluster members");
      else if (node->clust->clustersize > 1) {
	fprintf( handle, "Cluster_%d:\n", node->nodenumber);
	for( i=0; i < node->clust->clustersize; i++) {
	  if ( i == 0 )
 	    fprintf( handle, "\t" );
	  else if ( i % 4 == 0 )
	    fprintf( handle, "\n\t" );
	  else 
	    fprintf( handle, ", ");
	  fprintf( handle, "%s", node->clust->members[i]->name);
	}
	fprintf( handle, "\n\n");
      }
    }
    else {
      write_clustering_data_Tnode( handle, node->left );
      write_clustering_data_Tnode( handle, node->right );
    }
  }
}



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

void write_debug_Tnode( FILE *handle, struct Tnode *node, unsigned int offset) {
  unsigned int i,j;

  if (node != NULL) {
    /* We need to determine whether the node is a leaf or internal;
       since in this implementation internal nodes do not have names,
       it is sufficient to check the nodes name for nullness; however,
       this precludes the possibilities of internal nodes being given
       names  in the future, hence the check for internalness is made
       on the basis of the nullness of the children.
    */
    
    if ( node->left == NULL && node->right == NULL ) {
      /* this is a leaf node, so its cluster must have members; pain if not */
      if (node->clust->clustersize == 0 || node->clust->members == NULL)
	fatal_util( "Fatal Error: encountered a leaf node with no cluster info"); 
      else {
	for (i=0; i < node->clust->clustersize; i++) {
	  for (j=0; j < offset; j++) fprintf( handle, " ");
	  /* all leaves in the cluster are printed at the same offset */
	  fprintf( handle, 
//		   "%d:%s:%.5f\n", 
		   "%d:%s:%g\n", 
		   node->nodenumber, 
		   node->clust->members[i]->name,
		   node->distance); 
	}
      }
    }
    else if ( node->left != NULL && node->right != NULL ) {
      for (j=0; j < offset; j++) fprintf( handle, " ");
//      fprintf(handle, "Node %d:%.5f\n", node->nodenumber, node->distance);
      fprintf(handle, "Node %d:%g\n", node->nodenumber, node->distance);
      write_debug_Tnode( handle, node->left, offset+2);
      write_debug_Tnode( handle, node->right, offset+2);
    }
    /* else do nothing */
    
    fflush( handle );
  }
}



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

void write_debug_Tree( FILE *handle, struct Tree *thetree) {
  if (thetree != NULL) {
      write_debug_Tnode( handle, thetree->child[0], 0 );
      write_debug_Tnode( handle, thetree->child[1], 0 );
      write_debug_Tnode( handle, thetree->child[2], 0 );

      fflush( handle );
  }

}



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

void write_MUL_flattened_Tnode( FILE *handle, struct Tnode *node) {
  unsigned int i,j;

  if (node != NULL) {
    write_MUL_flattened_Tnode( handle, node->left );
    if (node->clust != NULL) {
      for( i=0; i < node->clust->clustersize; i++ ) {
	fprintf( handle, "%-24s", node->clust->members[i]->name );
	for (j=0; j < node->clust->members[i]->length; j++)
	  fprintf( handle, "%c", node->clust->members[i]->seq[j] );
	fprintf( handle, "\n");
      }
    }
    write_MUL_flattened_Tnode( handle, node->right );
  }

}



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

void write_MUL_flattened_Tree( FILE *handle, struct Tree *tr) {
  if (tr != NULL) {
    write_MUL_flattened_Tnode( handle, tr->child[0] );
    write_MUL_flattened_Tnode( handle, tr->child[1] );
    write_MUL_flattened_Tnode( handle, tr->child[2] );
  }
  fflush( handle );

}



/********************************************************************* 
 FUNCTION: write_newhampshire_Tnode
 DESCRIPTION: 
   prints the given Tree in 'New Hampshire' text format to the given 
   file handle
 RETURNS:
 ARGS: 
   File handle
   TNode *
   Whether or not to show bootstrap values
 NOTES:
 *********************************************************************/

void write_newhampshire_Tnode( FILE *handle, struct Tnode *node, 
			       unsigned int show_bootstraps  ) {
  if (node != NULL) {
    /* We need to determine whether the node is a leaf or internal;
       since in this implementation internal nodes do not have names,
       it is sufficient to check the nodes name for nullness; however,
       this precludes the possibilities of internal nodes being given
       names  in the future, hence the check for internalness is made
       on the basis of the nullness of the children.
    */
    
    if ( node->left == NULL && node->right == NULL) {
      /* this is a leaf node, so its cluster must have members; pain if not */
      if (node->clust->clustersize == 0 || node->clust->members == NULL)
	fatal_util( "Fatal Error: encountered a leaf node with no cluster info"); 
      else if (node->clust->clustersize == 1) 
//	fprintf( handle, "%s:%.5f", node->clust->members[0]->name, node->distance );
	fprintf( handle, "%s:%g", node->clust->members[0]->name, node->distance );
      else {
	/* if there is more than one sequence belonging to the cluster, then this piece
	   of code will generate som internal nodes in the output tree with no bootstrap
	   values. Such is life... */
        unsigned int i;

	for (i=0; i < node->clust->clustersize - 1; i++) {
//	 fprintf( handle, "(\n%s:%.5f,\n", node->clust->members[i]->name, 0.0 ); 
	 fprintf( handle, "(%s:%g,", node->clust->members[i]->name, 0.0 ); 
	}
//	fprintf( handle, "%s:%.5f)\n", node->clust->members[i]->name, 0.0);
	fprintf( handle, "%s:%g)", node->clust->members[i]->name, 0.0);
	for (i=0; i < node->clust->clustersize - 2; i++) {
//	 fprintf( handle, ":%.5f)\n", 0.0 ); 
	 fprintf( handle, ":%g)", 0.0 ); 
	}

//	fprintf( handle, ":%.5f", node->distance);
	fprintf( handle, ":%g", node->distance);
	/* fprintf( handle, "Cluster_%d:%.5f", node->nodenumber, node->distance ); */
      }
    }
    else if ( node->left != NULL && node->right != NULL ) {
//      fprintf( handle, "(\n");
      fprintf( handle, "(");
      write_newhampshire_Tnode( handle, node->left, show_bootstraps );
//      fprintf( handle, ",\n" );
      fprintf( handle, "," );
      write_newhampshire_Tnode( handle, node->right, show_bootstraps );
      if (show_bootstraps) {
//	fprintf( handle, ")\n%d:%.5f", node->bootstrap, node->distance );
	fprintf( handle, ")%d:%g", node->bootstrap, node->distance );
      }
      else {
//	fprintf( handle, ")\n:%.5f", node->distance);
	fprintf( handle, "):%g", node->distance);
      }
    }
    /* else do nothing */
  }
  fflush( handle );
}



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

void write_newhampshire_Tree( FILE *handle, struct Tree *thetree,
			      unsigned int show_bootstraps) {

  /* write_newhampshire_Tnode always places parenttheses around the
     sub-tree if there is more than one sub-node, but this is not
     appropriate if there there are no other sub-trees to draw at
     this level. Hence there is a special case for a single
     sub-tree (as returned by the UPGMA method) */

  if (thetree != NULL) {
    if (thetree->child[0] != NULL) {
      if (thetree->child[1] == NULL) {
	/* draw rooted tree */

	if (thetree->child[0]->left != NULL && thetree->child[0]->right != NULL) {
//	  fprintf( handle, "(\n");
	  fprintf( handle, "(");
	  write_newhampshire_Tnode( handle, thetree->child[0]->left, show_bootstraps );
//	  fprintf( handle, ",\n");
	  fprintf( handle, ",");
	  write_newhampshire_Tnode( handle, thetree->child[0]->right, show_bootstraps );
//	  fprintf( handle, ");\n");
	  fprintf( handle, ");");
	}
	else {
	  /* this is a leaf node, and leaf nodes may contain a single sequence
	     or a cluser of sequences. If this leaf contains a single sequence,
	     then we have a tree of one sequence, in which case we print an
	     error, because trees of one sequence do not make sense */
	  
	  if (thetree->child[0]->clust->clustersize == 1) 
	    fatal_util( "Cannot build a tree with a single sequence %s",
		     thetree->child[0]->clust->members[0]->name);
	  else {
	    unsigned int i;
	    
	    for (i=0; i < thetree->child[0]->clust->clustersize - 1; i++)
//	      fprintf( handle, "(\n%s:%.5f,\n", thetree->child[0]->clust->members[i]->name, 0.0 ); 
	      fprintf( handle, "(%s:%g,", thetree->child[0]->clust->members[i]->name, 0.0 ); 

//	    fprintf( handle, "%s:%.5f", thetree->child[0]->clust->members[i]->name, 0.0);
	    fprintf( handle, "%s:%g", thetree->child[0]->clust->members[i]->name, 0.0);

	    for (i=0; i < thetree->child[0]->clust->clustersize - 2; i++) 
//	      fprintf( handle, ")\n:%.5f)\n", 0.0 ); 
	      fprintf( handle, "):%g)", 0.0 ); 

//	    fprintf( handle, ");\n"); 
	    fprintf( handle, ");"); 
	  }
	}
      }
      else {
//	fprintf( handle, "(\n");
	fprintf( handle, "(");
	write_newhampshire_Tnode( handle, thetree->child[0], show_bootstraps );
//	fprintf( handle, ",\n");
	fprintf( handle, ",");
	write_newhampshire_Tnode( handle, thetree->child[1], show_bootstraps );
	if (thetree->child[2] != NULL) {
//	  fprintf( handle, ",\n");
	  fprintf( handle, ",");
	  write_newhampshire_Tnode( handle, thetree->child[2], show_bootstraps );
	}
//	fprintf( handle, ");\n");
	fprintf( handle, ");");
      }
    }
  }
  fflush( handle );
}


