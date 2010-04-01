/*  Last edited: Mar 18 11:52 2002 (klh) */
/**********************************************************************
 ** FILE: buildtree.c
 ** NOTES:
 **  Contains functions for building trees from distance matrices
 **  (and vice versa)
 **********************************************************************/

#include "buildtree.h"

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
void export_distances_buildtree( struct Tree *thetree, 
				 struct DistanceMatrix *mat) {

  struct Tnode *root;
  unsigned int *considered;

  /* first, create the imaginary node that sits between the three nodes
     of the unrooted tree */

  root = new_interior_Tnode( thetree->numnodes );

  /* now set up the connections */

  root->parent = thetree->child[0];
  root->left = thetree->child[1];
  root->right = thetree->child[2];
  root->left->parent = root;
  root->right->parent = root;
  root->parent->parent = root; /* This one seems odd but we need it */
  root->distance = root->parent->distance;

  /* A boolean array to store which nodes have been considered */

  considered = (unsigned int *) 
    malloc_util( (thetree->numnodes + 1) * sizeof( unsigned int ) );

  leaf_find_buildtree( root->parent, mat, considered, thetree->numnodes + 1);
  leaf_find_buildtree( root->left, mat, considered, thetree->numnodes + 1);
  leaf_find_buildtree( root->right, mat, considered, thetree->numnodes + 1);


  considered = free_util( considered );

  /* When freeing the temporary root we have created, we must nullify
     the parent and children. We don't want to cause a cascading delete
     of the tree */

  root->parent = NULL;
  root->left = NULL;
  root->right = NULL;
  root = free_Tnode( root );

}




/**********************************************************************
 FUNCTION: find_path_buildtree
 DESCRIPTION: 
   Given a leaf node, recursively calculates the branch-length distance
   from the node to all other leaves in the tree, and places it in the
   appropriate part of the DistanceMatrix
 ARGS: 
   unsigned int (the number of the node from which we are finding all distances)
   Tnode (the current node under consideration)
   DistanceMatrix
   Boolean array (to store which nodes have already been considered)
   unsigned int (the size of this boolean array)
 RETURNS: 
 NOTES: 
 **********************************************************************/
void find_path_buildtree( unsigned int home,
			  struct Tnode *node,
			  struct DistanceMatrix *mat,
			  Distance dist,
			  unsigned int *considered) {

  if (node->left == NULL && node->right == NULL) {
    if (home < node->nodenumber) {
      mat->data[node->nodenumber][home] = dist;
    }
    else {
      mat->data[home][node->nodenumber] = dist;
    }
  }
  else {
    /* we have an internal node */
    considered[node->nodenumber] = 1;

    if (node->left != NULL) {
      if (! considered[node->left->nodenumber]) {
	find_path_buildtree( home, 
			     node->left, 
			     mat, 
			     dist + node->left->distance,
			     considered );
      }
    }
    if (node->right != NULL) {
      if (! considered[node->right->nodenumber]) {
	find_path_buildtree( home, 
			     node->right, 
			     mat, 
			     dist + node->right->distance,
			     considered );
      }
    }
    if (node->parent != NULL) {
      if (! considered[node->parent->nodenumber]) {
	find_path_buildtree( home, 
			     node->parent, 
			     mat, 
			     dist + node->distance,
			     considered );
      }
    }
  }

}


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
 RETURNS: 
 NOTES: 
 **********************************************************************/
void leaf_find_buildtree( struct Tnode *node, 
			  struct DistanceMatrix *mat, 
			  unsigned int *considered,
			  unsigned int size) {

  unsigned int i;

  if (node->left != NULL || node->right != NULL) {
    if (node->left != NULL) 
      leaf_find_buildtree( node->left, mat, considered, size );
    if (node->right != NULL) 
      leaf_find_buildtree( node->right, mat, considered, size );
  }
  else {

    mat->data[node->nodenumber][node->nodenumber] = 0.0;
    for(i=0; i < size; i++) considered[i] = 0;
    considered[node->nodenumber] = 1;
    find_path_buildtree( node->nodenumber, 
			 node->parent, 
			 mat, 
			 node->distance, 
			 considered);
  }

} 




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
struct Tree *neighbour_joining_buildtree( struct ClusterGroup *group,
					  unsigned int bootstrap) { 
  unsigned int numseqs, i, j;       /*** The current pair of nodes   ***/
  unsigned int k, m, nodecount;     /*** loop counters               ***/
  unsigned int row, column;         /*** matrix indices              ***/
  unsigned int mini = 0, minj = 0;  /*** neighbouring nodes          ***/
  unsigned int nextfreenode;        /*** incremental labels to nodes ***/ 
  unsigned int leftovers[3];        /*** three remaining nodes       ***/
  struct DistanceMatrix *mat;
  struct Tree *theTree;
  struct Tnode **nodes;             /*** starts off holding leaves    ***/
  struct Tnode *newnode;
  double fnumseqs;                  /*** divisor for sums            ***/
  double dmj, dmi, ri, minsofar, dist, dij, dist_i, dist_j, dist_k;
  Distance *r;                        /*** stores the r values         ***/


  /* METHOD ***********************************
     I have implemented the neightbour-joining algorithm as descibed
     by Durbin, Eddy, Krogh and Mitchison (1998 pp 171-172); see the
     text for reasoning behind variable names. I have inserted several
     changes from the basic algorithm for reasons of efficiency

     To keep track of dead nodes, we use the 'nodes' pointer array.
     Whenever we amalgamate nodes, we overwrite one of the entries
     with the new node and set the other to NULL. The NULL entries
     in the nodes array thus mark the dead nodes; this is important
     when updating the distance matrix because we must not add terms
     from nodes that are gone.

     Note the complex indexing of the distance matrix, since it is
     assumed that it is triangularised at the top-right. Two Ways
     around this: 
     1. Abstract the indexing into a distancemat method
     2. Assume a square symmetrical matrix
     Method one may have performance impact; method two is more viable
     but need to be much more cafeful with the updating of the matrix
     at the end of each iteration

     I have also used a modified version of Bill Bruno's idea to attempt
     to eliminate negative branch lengths in generated trees.

  ********************************************/



  /******* intialisation ********************/

  nodes = (struct Tnode **) malloc_util( group->numclusters * sizeof(struct Tnode *));
  for( i=0; i < group->numclusters; i++) { 
    nodes[i] = new_leaf_Tnode( i, clone_Cluster( group->clusters[i]) );
  }

  mat = group->matrix;
  numseqs = nextfreenode =  mat->size;
  
  theTree = empty_Tree();

  fnumseqs = (double) numseqs;
  

  if (numseqs > 2) {

    r = (Distance *) malloc_util( numseqs * sizeof( Distance ) );
    /* Calculate r[i] for all i */
    for( i=0; i < numseqs; i++ ) {
      ri = 0.0;
      for (k=0; k < numseqs; k++) {
	if (k > i) ri += mat->data[k][i];
	else ri += mat->data[i][k];
      }
      r[i] = ri / (fnumseqs - 2.0);
    }
    
    
    /******* main loop ************************/
    
    for (nodecount=0; nodecount < numseqs-3; nodecount++) {

      /* do the intialisation necessary for each iteration here */

      minsofar = FLT_MAX;  /* from float.h */

      /******* for each pair of matrix entries *********************/

      for( i=0; i < numseqs; i++ ) {
	if (nodes[i] == NULL) continue;
	for( j=0; j < i; j++ ) {
	  if (nodes[j] == NULL) continue;
	  
	  dist = mat->data[i][j] - (r[i] + r[j]);
	  if (dist < minsofar) {
	    minsofar = dist;
	    mini = i;
	    minj = j;
	  }
	}
      }
      
      /* printf("i = %d, j = %d\n", mini, minj); */
      
      /* we have the neighbouring i, j; lets calc distances and make the new node */

      dij = mat->data[mini][minj];
      dist_i = (dij + r[mini] - r[minj]) * 0.5;
      dist_j = dij - dist_i;

      /* Adjustment to allow for negative branch lengths */
      if (dist_i < 0.0) {
	dist_i = 0.0;
	dist_j = dij;
	if (dist_j < 0.0)
	  dist_j = 0.0;
      }
      else if (dist_j < 0.0) {
	dist_j = 0.0;
	dist_i = dij;
	if (dist_i < 0.0)
	  dist_i = 0.0;
      }

      nodes[mini]->distance = dist_i;
      /* printf("DistanceI = %f\n", nodes[mini]->distance); */
      
      nodes[minj]->distance = dist_j;
      /* printf("DistanceJ = %f\n", nodes[minj]->distance); */
      
      newnode = new_interior_Tnode( nextfreenode++ );
      newnode->left = nodes[mini];
      newnode->right = nodes[minj];
      nodes[mini]->parent = newnode;
      nodes[minj]->parent = newnode;
      nodes[mini] = newnode;
      nodes[minj] = NULL;
      if (bootstrap) { 
	/* we need to create and load the 'bit' field of child ids */
	
	newnode->child_ids = (unsigned int *) 
	  malloc_util( group->numclusters * sizeof( unsigned int ) );
	for (m=0; m < group->numclusters; m++) {
	  if ( (newnode->left->child_ids != NULL && newnode->left->child_ids[m]) ||
	       (newnode->right->child_ids != NULL && newnode->right->child_ids[m]) ||
	       (newnode->left->nodenumber == m || 
		newnode->right->nodenumber == m ||
		newnode->nodenumber == m)) {
	    newnode->child_ids[m] = 1;
	  }
	  else {
	    newnode->child_ids[m] = 0;
	  }
	}
      }
      
      /* now update the distance matrix; This needs hackery to make sure that the
	 indexing is correct */

      r[mini] = 0.0;  /* This is the only r[i] that requires wholesale changes */
      for( m=0; m < numseqs; m++ ) {
	if (nodes[m] == NULL) continue;
	
	if (m != mini) {
	  
	  if (m > minj) dmj = mat->data[m][minj];
	  else dmj = mat->data[minj][m];
	  
	  if (m > mini) {
	    row = m;
	    column = mini;
	  }
	  else {
	    row = mini;
	    column = m;
	  }
	  dmi = mat->data[row][column];
	  
	  /* we can actually adjust r[m] here, by using the form:
	     rm = ((rm * numseqs) - dmi - dmj + dmk) / (numseqs-1)
	  */
	  
	  /* Note: in Bill Bruno's method for negative branch elimination, then if either
	     dist_i is positive and dist_j is 0, or dist_i is zero and dist_j is positive
	     (after adjustment) then the matrix entry is formed from the distance to the
	     node in question (m) to the node with the zero branch length (whichever it was).
	     I think my code already has the same effect; this is certainly true if dij is
	     equal to dist_i + dist_j, which it should have been fixed to
	  */
	  
	  mat->data[row][column] = (dmi + dmj - dij) * 0.5;
	  r[m] = ((r[m] * (fnumseqs - 2.0)) - dmi - dmj + mat->data[row][column]) / (fnumseqs - 3.0); 
	  r[mini] += mat->data[row][column];
	}
      }
      
      fnumseqs -= 1.0;
      r[mini] /= fnumseqs - 2.0;
      
    }
    /******* end of main loop ******************/
    
    
    /* Now there are just 3 nodes left. Need to locate those three nodes */
    /* The following looks dangerous, because we only have room for three nodes in the tree;
       However, all nodes except three should be NULL, or else something is seriously wrong */
    
    for(k=0, m=0; k < numseqs; k++) {
      if (nodes[k] != NULL) { 
	theTree->child[m] = nodes[k];
	leftovers[m++] = k;
	nodes[k] = NULL; /* so that the nodes are not released when input is */
      }
    }
    

    /* Now to get rid of those negative branch lengths (see you in ~70 lines....) */
    
    dist_i = theTree->child[0]->distance =
      (mat->data[leftovers[1]][leftovers[0]] +
       mat->data[leftovers[2]][leftovers[0]] -
       mat->data[leftovers[2]][leftovers[1]]) * 0.5;
    dist_j = theTree->child[1]->distance = mat->data[leftovers[1]][leftovers[0]] - theTree->child[0]->distance;
    dist_k = theTree->child[2]->distance = mat->data[leftovers[2]][leftovers[0]] - theTree->child[0]->distance;


    if (dist_i < 0.0) {
      dist_i = 0.0;
      dist_j = mat->data[leftovers[1]][leftovers[0]];
      dist_k = mat->data[leftovers[2]][leftovers[0]];
      if (dist_j < 0.0) {
	dist_j = 0.0;
	dist_k = ( mat->data[leftovers[2]][leftovers[0]] +
		   mat->data[leftovers[2]][leftovers[1]] ) * 0.5;
	if (dist_k < 0.0) 
	  dist_k = 0.0;
      }
      else if (dist_k < 0.0) {
	dist_k = 0.0;
	dist_j = 
	  (mat->data[leftovers[1]][leftovers[0]] +
	   mat->data[leftovers[2]][leftovers[1]] ) * 0.5;
	if (dist_j < 0.0)
	  dist_j = 0.0;
      }
    }
    else if (dist_j < 0.0) {
      dist_j = 0.0;
      dist_i = mat->data[leftovers[1]][leftovers[0]];
      dist_k = mat->data[leftovers[2]][leftovers[1]];
      if (dist_i < 0.0) {
	dist_i = 0.0;
	dist_k = ( mat->data[leftovers[2]][leftovers[0]] +
		   mat->data[leftovers[2]][leftovers[1]] ) * 0.5;
	if (dist_k < 0.0) 
	  dist_k = 0.0;
      }
      else if (dist_k < 0.0) {
	dist_k = 0.0;
	dist_i = 
	  (mat->data[leftovers[1]][leftovers[0]] +
	   mat->data[leftovers[2]][leftovers[0]] ) * 0.5;
	if (dist_i < 0.0)
	  dist_i = 0.0;
      }
    }
    else if (dist_k < 0.0) {
      dist_k = 0.0;
      dist_i = mat->data[leftovers[2]][leftovers[0]];
      dist_j = mat->data[leftovers[2]][leftovers[1]];
      if (dist_i < 0.0) {
	dist_i = 0.0;
	dist_j = ( mat->data[leftovers[1]][leftovers[0]] +
		   mat->data[leftovers[2]][leftovers[1]] ) * 0.5;
	if (dist_j < 0.0) 
	  dist_j = 0.0;
      }
      else if (dist_j < 0.0) {
	dist_j = 0.0;
	dist_i = 
	  (mat->data[leftovers[1]][leftovers[0]] +
	   mat->data[leftovers[2]][leftovers[0]] ) * 0.5;
	if (dist_i < 0.0)
	  dist_i = 0.0;
      }
    }



    theTree->child[0]->distance = dist_i; 
    theTree->child[1]->distance = dist_j;
    theTree->child[2]->distance = dist_k;


    r = free_util( r );
  }
  else {
    /* deal with the trivial case of less than three leaves */

    for( i=0; i < numseqs; i++ ) {
      theTree->child[i] = nodes[i];
      nodes[i] = NULL;
    }

    if (numseqs == 2) {
      theTree->child[0]->distance = mat->data[1][0] * 0.5;
      theTree->child[1]->distance = mat->data[1][0] * 0.5;
    }
  }
 
  /* nextfreenode records how many calls to new_Tnode_tree we made; therefore, at this stage, it
     holds the total number of nodes (including leaves) in the tree */
  theTree->numnodes = nextfreenode;
  
  nodes = free_util( nodes );
  
  /* The caller of the function should free the tree when finished with it */
  
  return theTree;
}



/**********************************************************************
 FUNCTION: UPGMA_buildtree
 DESCRIPTION: 
   Returns a phylogenetic tree of the sequences in the 
   given alignment, using the Unweighted Pair-Group method based on
   Arithmentic Averages (UPGMA)
 ARGS: 
   A ClusterGroup pointer (cluster.h)
   A boolean,for whether to record information needed later for bootstrapping
 RETURNS: 
   struct Tree (trees.h)
 NOTES: The function allocates all the memory necessary for the tree.
   The caller should call free_Tnode (tree.h) to free this memory when
   the tree is no longer needed

   This algorithm produces a rooted tree; hence the returned Tree
   will have the root in child[0]; child[1] and child[2] (used to
   represent the trichotomy of unrooted trees) will be NULL 
 **********************************************************************/
struct Tree *UPGMA_buildtree(struct ClusterGroup *group,
			     unsigned int bootstrap) { 

  unsigned int numseqs, i, j;       /*** The current pair of nodes   ***/
  unsigned int m, nodecount;        /*** loop counters               ***/
  unsigned int row, column;         /*** matrix indices              ***/
  unsigned int mini = 0, minj = 0;  /*** neighbouring nodes          ***/
  unsigned int nextfreenode;        /*** incremental labels to nodes ***/ 
  unsigned int *subtreesizes;
  Distance *heights;
  double minsofar, newnodeheight, dmi, dmj;
  
  struct DistanceMatrix *mat = NULL;
  struct Tree *theTree = NULL;
  struct Tnode **nodes = NULL;           /*** starts off holding leaves    ***/
  struct Tnode *newnode = NULL;

  /* METHOD ***********************************
     I have implemented the UPGMA algorithm as presented 
     by Durbin, Eddy, Krogh and Mitchison (1998 pp 171-172).

     To keep track of merged nodes, we use the 'nodes' pointer array.
     Whenever we amalgamate nodes, we overwrite one of the entries
     with the new node and set the other to NULL. The NULL entries
     in the nodes array thus mark the dead nodes; this is important
     when updating the distance matrix because we must not add terms
     from nodes that are gone.

     Note the complex indexing of the distance matrix, since it is
     assumed that it is triangularised at the top-right. Two Ways
     around this: 
     1. Abstract the indexing into a distancemat method
     2. Assume a square symmetrical matrix
     Method one may have performance impact; method two is more viable
     but need to be much more cafeful with the updating of the matrix
     at the end of each iteration

  ********************************************/



  /******* intialisation ********************/

  nodes = (struct Tnode **) malloc_util( group->numclusters * sizeof(struct Tnode *));
  heights = (Distance *) malloc_util( group->numclusters * sizeof(Distance));
  subtreesizes = (unsigned int *) malloc_util( group->numclusters * sizeof(unsigned int));
  for( i=0; i < group->numclusters; i++) { 
    nodes[i] = new_leaf_Tnode( i, group->clusters[i] );
    /* clusters[i] in group has conceptually moved into the Tnode, so... */
    group->clusters[i] = NULL;
    
    heights[i] = 0.0;
    subtreesizes[i] = 1;
  }

  mat = group->matrix;
  numseqs = nextfreenode =  mat->size;

  theTree = empty_Tree();
   
  /******* main loop ************************/

  if (numseqs == 1)
    newnode = nodes[0];
  else {
    /*** new code ***/
    Distance *min_of_row = (Distance *) malloc_util( numseqs * sizeof( Distance ));
    int *min_of_row_idx = (int *) malloc_util( numseqs * sizeof( int ));

    /* calculate in advance the minimum column for each row */
    for(i=1; i < numseqs; i++) {
      min_of_row[i] = FLT_MAX;
      min_of_row_idx[i] = -1;

      for(j=0; j < i; j++) {

	if (mat->data[i][j] < min_of_row[i]) {
	  min_of_row[i] = mat->data[i][j];
	  min_of_row_idx[i] = j;
	}
      }
    }


    for (nodecount=0; nodecount < numseqs-1; nodecount++) {

      /* do the intialisation necessary for each iteration here */
      
      minsofar = FLT_MAX;  /* from float.h */
      
      for( i=1; i < numseqs; i++) {
	if (nodes[i] == NULL) continue;
	if ( min_of_row[i] < minsofar ) {
	  minsofar = min_of_row[i];
	  mini = i;
	  minj = min_of_row_idx[i];
	}
      }

      /* we have the neighbouring i, j; lets calc distances and make the new node */
      
      newnodeheight = mat->data[mini][minj] * 0.5;
      nodes[mini]->distance = newnodeheight - heights[mini];
      nodes[minj]->distance = newnodeheight - heights[minj];
      
      newnode = new_interior_Tnode( nextfreenode++ );
      newnode->left = nodes[mini];
      newnode->right = nodes[minj];
      nodes[mini]->parent = newnode;
      nodes[minj]->parent = newnode;
      if (bootstrap) { 
	/* we need to create and load the 'bit' field of child ids */
	
	newnode->child_ids = (unsigned int *) 
	  malloc_util( group->numclusters * sizeof( unsigned int ) );
	for (i=0; i < group->numclusters; i++) {
	  if ( (newnode->left->child_ids != NULL && 
		(newnode->left->child_ids[i] || newnode->left->child_ids[i])) ||	       
	       (newnode->right->child_ids != NULL &&
		(newnode->right->child_ids[i] || newnode->right->child_ids[i])) ||
		(newnode->left->nodenumber == i || 
		 newnode->right->nodenumber == i ||
		 newnode->nodenumber == i)) {
	    newnode->child_ids[i] = 1; 
	  }
	  else {
	    newnode->child_ids[i] = 0;
	  }
	}
      }

      /* now update the distance matrix; This needs hackery to make sure that the
	 indexing is correct */
      
      min_of_row[mini] = FLT_MAX;
      min_of_row_idx[mini] = -1;

      for( m=0; m < numseqs; m++ ) {
	if (nodes[m] == NULL) continue;
	
	/*                  dmini,m*|Ci| + dminj,m*|Cj|
	     dmini,m =      ---------------------------   
			            |Ci| + |Cj|
	*/
	
	if (m > minj) dmj = mat->data[m][minj];
	else dmj = mat->data[minj][m];
	
	if (m > mini) {
	  row = m;
	  column = mini;
	}
	else {
	  row = mini;
	  column = m;
	}

	dmi = mat->data[row][column];
	
	mat->data[row][column] = 
	  (( dmi * subtreesizes[mini])+( dmj * subtreesizes[minj])) / 
	  (subtreesizes[mini] + subtreesizes[minj]);

	if (mat->data[row][column] < min_of_row[row] 
	    && nodes[column] != NULL
	    && column != minj
	    && row != column) {
	  min_of_row[row] = mat->data[row][column];
	  min_of_row_idx[row] = column; 
	}

	/* node minj is being nullified, so any row that points to it
	   is min should be updated */
	if (m != mini && m != 0) {
	  if (min_of_row_idx[m] == minj || min_of_row_idx[m] == mini) {
	    min_of_row[m] = FLT_MAX;
	    min_of_row_idx[m] = -1;
	    for(j=0; j < m; j++) {
	      if (nodes[j] == NULL  || j == minj) continue;
	      if (mat->data[m][j] < min_of_row[m]) {
		min_of_row[m] = mat->data[m][j];
		min_of_row_idx[m] = j;
	      }
	    }
	  }
	}
      }

      heights[mini] = newnodeheight;
      subtreesizes[mini] = subtreesizes[mini] + subtreesizes[minj] + 1; 
      
      nodes[mini] = newnode;
      nodes[minj] = NULL;
      
    }

    free_util( min_of_row );
    free_util( min_of_row_idx );

  }
  /******* end of main loop ******************/
  
    
  theTree->child[0] = newnode;
  
  /* nextfreenode records how many calls to new_Tnode_tree we made; therefore, at this stage, it
     holds the total number of nodes (including leaves) in the tree */
  
  theTree->numnodes = nextfreenode;
  nodes = free_util( nodes );
  heights = free_util( heights );
  subtreesizes = free_util( subtreesizes );

  /* The caller of the function should free the tree when finished with it */
  
  return theTree;

}
