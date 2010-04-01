/*  Last edited: Feb  1 18:18 2002 (klh) */
/**********************************************************************
 ** FILE: cluster.c
 ** NOTES:
 **  A DistanceMatrix should always be part of a Cluster 
 **  It makes no sense to have a set of pairwise distances without the 
 **  associated sequences (even if we just store their names)
 **********************************************************************/

#include "cluster.h"



/*********************************************************************
  FUNCTION: alignment_to_ClusterGroup
  DESCRIPTION: 
    This function returns a ClusterGroup, given an Alignment. 
    if the secind arg is true, In doing indentical sequences in the 
    alignment are merged. If bootstrapping is required, the consensus
    alignment can be extracted from the ClusterGroup using 
    get_consensus_from_ClusterGroup
  RETURNS: struct ClusterGroup
  ARGS: 
    1. A source Alignment pointer
    2. A boolean specifying whether duplicate sequences should be
       merged.
  NOTES: 
*********************************************************************/
struct ClusterGroup *alignment_to_ClusterGroup( struct Alignment *aln,
						unsigned int remove_duplicates) {
  unsigned int i, j, numclusters;
  struct ClusterGroup *group;
  struct Cluster **newclusts;

  group = empty_ClusterGroup();

  /* Need to create a cluster for every sequence in the given cluster */

  newclusts = (struct Cluster **) malloc_util( aln->numseqs * sizeof( struct Cluster *) );


  for( i=0; i < aln->numseqs; i++) {
    newclusts[i] = single_Sequence_Cluster( clone_Sequence(aln->seqs[i]));
  }
  numclusters = aln->numseqs;  

  if (remove_duplicates) {
    for( i=0; i < aln->numseqs; i++) {
      if (newclusts[i] == NULL) continue;
      for( j=i+1; j < aln->numseqs; j++) {
	if (newclusts[j] == NULL) continue;
	
	if (strncmp( aln->seqs[i]->seq, aln->seqs[j]->seq, aln->length) == 0) {
	  /* these two clusters contain the same sequence, so we can merge them */
	  newclusts[j] = merge_Cluster( newclusts[i], newclusts[j] );
	  numclusters--;
	}
      }
    }
  }


  /* newclusts will now be a sparse array containing clusters of
     identical sequences */

  group->numclusters = numclusters;
  group->clusters = (struct Cluster **) malloc_util( numclusters * sizeof( struct Cluster *) );
  for(i=0, j=0; i < aln->numseqs; i++) {
    if (newclusts[i] != NULL) {
      group->clusters[j++] = newclusts[i];
    }
  }

  newclusts = free_util( newclusts );

  return group;
}


/********************************************************************* 
 FUNCTION: clone_Cluster
 DESCRIPTION: 
   This function makes a complete copy of the given Cluster
   and returns it
 RETURNS: struct Cluster *
 ARGS: 
   struct Cluster *
 NOTES: 
*********************************************************************/
struct Cluster *clone_Cluster( struct Cluster *source) {
  unsigned int i;
  struct Cluster *dest = NULL; 
  
  if (source != NULL) {
    dest = empty_Cluster();
    dest->clustersize = source->clustersize;
    dest->members = (struct Sequence **) malloc_util( dest->clustersize
						      * sizeof( struct Sequence * ));
    for( i=0; i < source->clustersize; i++) {
      dest->members[i] = clone_Sequence( source->members[i] );
    }
    dest->consensus = clone_Sequence( source->consensus );
  }
  
  return dest;
}



/*********************************************************************
  FUNCTION: consensus_aln_from_ClusterGroup
  DESCRIPTION: 
    This function creates an alignment by taking the consensus 
    sequences from each Cluster in the given ClusterGroup
  RETURNS: struct ClusterGroup
  ARGS: 
    1. A source Alignment pointer
    2. A boolean specifying whether duplicate sequences should be
       merged.
  NOTES: 
*********************************************************************/
struct Alignment *consensus_aln_from_ClusterGroup( struct ClusterGroup *grp) {
  unsigned int i;
  struct Alignment *cons;

  cons = (struct Alignment *) malloc_util( sizeof( struct Alignment ));
  cons->numseqs = grp->numclusters;
  cons->seqs = (struct Sequence **) malloc_util( grp->numclusters * sizeof( struct Sequence *));

  for (i=0; i < grp->numclusters; i++) {
    cons->seqs[i] = clone_Sequence( grp->clusters[i]->consensus );
  }
  /* All consensus seqs will be aligned, so we can pick any one to get 
     the alignment width */
  cons->length = cons->seqs[0]->length;

  return cons;
}





/********************************************************************* 
  FUNCTION: empty_Cluster
  DESCRIPTION: 
    This function handles the simple task of allocating the space
    for a new Cluster.
  RETURNS: struct Cluster *
  ARGS: 
  NOTES: 
*********************************************************************/
struct Cluster *empty_Cluster( void ) {
  struct Cluster *newclust;
  
  newclust = (struct Cluster *) malloc_util( sizeof( struct Cluster ));
  newclust->clustersize = 0;
  newclust->members = NULL;
  newclust->consensus = NULL;
  newclust->matrix = NULL; 
  
  return newclust;
  
}



/********************************************************************* 
  FUNCTION: empty_ClusterGroup
  DESCRIPTION: 
    This function handles the simple task of allocating the space
    for a new ClusterGroup
  RETURNS: struct Cluster *
  ARGS: 
  NOTES: 
*********************************************************************/
struct ClusterGroup *empty_ClusterGroup( void ) {
  struct ClusterGroup *group;
  
  group = (struct ClusterGroup *) malloc_util( sizeof(struct ClusterGroup));
  group->numclusters = 0;
  group->clusters = NULL;
  group->matrix = NULL;
  
  return group;
}




/********************************************************************* 
  FUNCTION: free_Cluster
  DESCRIPTION: 
    This function releases the memory used by this Cluster and all of its
    members
  RETURNS: A null pointer
  ARGS: 
    struct Cluster *
  NOTES: 
    In the majority of cases, all sequences in the Cluster come from
    an alignment, and if this alignment is subsequently needed (e.g.
    for bootstrapping) then then the seqs should not be freed. To 
    prevent this, the members field should be set to null by the caller
    to prevent freeing og the alignment
*********************************************************************/
void *free_Cluster( struct Cluster *given ) {
  unsigned int i;
  
  if (given != NULL) {
    if (given->members != NULL) {
      for( i=0; i < given->clustersize; i++) {
	given->members[i] =  free_Sequence( given->members[i] );
      }
      given->members = free_util( given->members );
    }
    given->matrix = free_DistanceMatrix( given->matrix );
    given->consensus = free_Sequence( given->consensus );
    given = free_util( given );
  }
  return given;
}




/********************************************************************* 
  FUNCTION: free_ClusterGroup
  DESCRIPTION: 
    This function releases the memory used by this Cluster and all of its
    members
  RETURNS: A null pointer
  ARGS: 
    struct Cluster *
  NOTES: 
*********************************************************************/
void *free_ClusterGroup( struct ClusterGroup *given ) {
  unsigned int i;
  
  if (given != NULL) {
    if (given->clusters != NULL) {
      for( i=0; i < given->numclusters; i++ ) {
	given->clusters[i] = free_Cluster( given->clusters[i] );
      }
      given->clusters = free_util( given->clusters );
    }
    given->matrix = free_DistanceMatrix( given->matrix );
    given = free_util( given );
  }
  
  return given;
}




/********************************************************************* 
  FUNCTION: merge_Cluster
  DESCRIPTION: 
    Adds the sequences in second arg to first arg, freeing the second
    arg, returning the result of this freeing (hopefully NULL);
  RETURNS: The result of freeing the second cluster (NULL if all is well)
  ARGS: 
    Destination Cluster *, 
    Source Cluster *
  NOTES:
*********************************************************************/
void *merge_Cluster( struct Cluster *dest, struct Cluster *source) {
  unsigned int i;
  
  /* take the sequences in source and add them onto the destination list */
  
  dest->members = (struct Sequence **) 
    realloc_util( dest->members, 
		  (dest->clustersize + source->clustersize) * sizeof(struct Sequence *));
  for ( i=0; i < source->clustersize; i++) {
    dest->members[ dest->clustersize++ ] = source->members[i];
    source->members[i] = NULL;
  }
  
  /* Need to update the consensus sequence with respect to the merge.
     At this stage, I am only merging identical clusters, so the
     consensus sequence is already correct in the destination 
  */

  source = free_Cluster( source );
  
  return source;
}



/********************************************************************* 
 FUNCTION: single_Sequence_Cluster
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new Cluster with the single given Sequence.
 RETURNS: struct Cluster *
 ARGS: 
   A pointer to a Sequence, or NULL for an empty Cluster
 NOTES: 
*********************************************************************/

struct Cluster *single_Sequence_Cluster( struct Sequence *seq) {
  struct Cluster *newclust;

  newclust = empty_Cluster();

  if (seq != NULL) {
    newclust->clustersize = 1;
    newclust->members = (struct Sequence **) malloc_util( sizeof( struct Sequence *) );
    newclust->members[0] = seq;
    newclust->consensus = clone_Sequence( seq );
    /* Note how the consensus sequence will end up with the same name as
       the first representative of trhe cluster. This is fine, because
       when the consensus is used, it is assumed that it is not a real
       sequence so its name is meaningless 
    */
  }

  return newclust;
}




/********************************************************************* 
 FUNCTION: single_Cluster_ClusterGroup
 DESCRIPTION: 
   This function takes the given cluster and very simlpy makes a 
   single-Cluster ClusterGroup from it
 RETURNS: struct ClusterGroup *
 ARGS: 
   A pointer to a Cluster
 NOTES: 
*********************************************************************/

struct ClusterGroup *single_Cluster_ClusterGroup( struct Cluster *clust ) {
  struct ClusterGroup *group;

  group = empty_ClusterGroup();

  if (clust != NULL) {
    group->numclusters = 1;
    group->clusters = (struct Cluster **) malloc_util( sizeof(struct Cluster *));
    group->clusters[0] = clust;
  }
  
  return group;
}








