/*  Last edited: Aug 25 15:19 1999 (klh) */
/**********************************************************************
 ** FILE: cluster.h
 ** NOTES:
 **  A DistanceMatrix should always be part of a Cluster 
 **  It makes no sense to have a set of pairwise distances without the 
 **  associated sequences (even if we just store their names)
 **********************************************************************/

#ifndef _CLUSTER
#define _CLUSTER

#include "sequence.h"
#include "distancemat.h"



/******************* structure definitions ****************************/

struct Cluster {
  unsigned int clustersize;
  struct Sequence **members;
  struct Sequence *consensus;
  struct DistanceMatrix *matrix;
};


/*
  Clusters contain groups of identical sequence. I intend to investigate
  methods where clusters contain groups of similar (not necessarily
  identical) sequences. The DistanceMatric field, although not currently
  used, will allow for the building of trees from these clusters
*/


struct ClusterGroup {
  unsigned int numclusters;
  struct Cluster **clusters;
  struct DistanceMatrix *matrix;
};


/********************** function prototypes ***************************/



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
struct ClusterGroup *alignment_to_ClusterGroup( struct Alignment *,unsigned int);

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
struct Cluster *clone_Cluster( struct Cluster *);

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
struct Alignment *consensus_aln_from_ClusterGroup( struct ClusterGroup *);

/********************************************************************* 
 FUNCTION: empty_Cluster
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new Cluster.
 RETURNS: struct Cluster *
 ARGS: 
 NOTES: 
 *********************************************************************/
struct Cluster *empty_Cluster( void );

/********************************************************************* 
 FUNCTION: empty_ClusterGroup
 DESCRIPTION: 
   This function handles the simple task of allocating the space
   for a new ClusterGroup
 RETURNS: struct Cluster *
 ARGS: 
 NOTES: 
 *********************************************************************/
struct ClusterGroup *empty_ClusterGroup( void );

/********************************************************************* 
 FUNCTION: free_Cluster
 DESCRIPTION: 
   This function releases the memory used by this Cluster and all of its
   members
 RETURNS: A null pointer
 ARGS: 
   struct Cluster *
 NOTES: 
 *********************************************************************/
void *free_Cluster( struct Cluster *);

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
void *free_ClusterGroup( struct ClusterGroup *);

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
void *merge_Cluster( struct Cluster *, struct Cluster *);

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
struct Cluster *single_Sequence_Cluster( struct Sequence *);

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
struct ClusterGroup *single_Cluster_ClusterGroup( struct Cluster *);


#endif
