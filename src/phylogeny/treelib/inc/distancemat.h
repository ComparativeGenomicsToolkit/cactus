/*  Last edited: Aug 25 15:25 1999 (klh) */
/**********************************************************************
 ** FILE: distancemat.h
 ** NOTES:
 **   Functions and types for the manipulation of Distance Matrices
 **********************************************************************/

#ifndef _DISTANCEMAT
#define _DISTANCEMAT

#include <math.h>

#include "util.h"
#include "align.h"
#include "time.h"


/******************* structure definitions ****************************/

typedef float Distance;

struct DistanceMatrix {
  Distance **data;
  int size;
};


/********************** function prototypes ***************************/


/*********************************************************************
 FUNCTION: calc_DistanceMatrix
 DESCRIPTION: 
   Produces a distance matrix from the given multiple alignment
 RETURNS: struct DistanceMatrix
 ARGS:
   A DistanceMatrix to fill in
   A multiple alignment
   A boolean indicating whether or not random columns should be used
     for purposes of bootstrapping
   A boolean indicating whether the Kimura distance adjustment is to 
       be used or not.
 NOTES:
   0. the given DistanceMatrix and Alignment should be of the same order

   1. The matrix produced is in bottom-left triangular format; don't you
   go trying to access that top-right section (I'm warning you...)

   2. At the moment, the function calculates distance based on sequence
   identity, using Kimura's function if that option is raised.

   3. If use_rand_cols is true, then the matrix is constructed using
   random sampling  of columns, for the purposes of bootstrapping. At 
   the moment, the native function 'rand' is used to do this, suitable 
   seeded by time (by the caller). This may prove unsatisfactory...

   4. Where no information is available to determine the distance 
   between two sequences, a value of twice the maximum observed 
   distance is assigned (inspiration from ISMB99 poster by Huson,
   Smith and Warnow).

 *********************************************************************/
void calc_DistanceMatrix(struct DistanceMatrix *, 
			 struct Alignment *,
			 unsigned int,
			 unsigned int );

/*********************************************************************
 FUNCTION: clone_DistanceMatrix
 DESCRIPTION: 
   Produces a brand new DistanceMatrix, identical to the source
 RETURNS: struct DistanceMatrix
 ARGS: 
   A source distane matrix
 NOTES: 
   1. The matrix produced is in bottom-left triangular format; don't you
   go trying to access that top-right section (I'm warning you...)
 *********************************************************************/
struct DistanceMatrix *clone_DistanceMatrix( struct DistanceMatrix *);

/*********************************************************************
 FUNCTION: empty_DistanceMatrix
 DESCRIPTION: 
   Produces an empty distance matrixof the given size, uninitialised
 RETURNS: struct DistanceMatrix
 ARGS: 
   The size of the matrix to be created
 NOTES: 
   1. The matrix produced is in bottom-left triangular format; don't you
   go trying to access that top-right section (I'm warning you...)
 *********************************************************************/
struct DistanceMatrix *empty_DistanceMatrix( unsigned int );

/*********************************************************************
 FUNCTION: free_DistanceMatrix
 DESCRIPTION: 
   Frees the memory for the given distance matrix
 RETURNS:
 ARGS: 
   struct DistanceMatrix *
 NOTES: 
 *********************************************************************/
void *free_DistanceMatrix( struct DistanceMatrix *);

/********************************************************************** 
 FUNCTION: index_DistanceMatrix
 DESCRIPTION: 
   indexes the given distance matrix with the given indices,
   returning the appropraite distance.
 RETURNS: distance (float)
 ARGS: 
   A distance matrix *
   row index
   column index
 NOTES: 
   This function is necessary to account for the fact that the distance 
   matrix may be implemented as a symmtrical or triangular matrix.
   It therefore abstracts the internals of the distance matrix, at the
   cost of a function call for each lookup (is this wise...?)
 **********************************************************************/
Distance index_DistanceMatrix( struct DistanceMatrix *, unsigned int, unsigned int );

/*********************************************************************
 FUNCTION: print_DistanceMatrix
 DESCRIPTION: 
   Prints the given distance matrix.
 RETURNS:
 ARGS: 
   struct DistanceMatrix *
 NOTES: 
   A DistanceMatrix does not exist in isolation in practice but as
   part of a Cluster (this is to maintain the tight coupling between 
   the matrix and the sequences for which it is expressing the distances). 
   Therefore, to read or write a useful distance
   matrix (for compatibility with the phylip package for example)
   use write_phylip_Cluster
 *********************************************************************/
void print_DistanceMatrix( FILE *handle,  struct DistanceMatrix * );


/********************************************************************* 
 FUNCTION: read_phylip_DistanceMatrix
 DESCRIPTION: 
   This function creates a DistanceMatrix from the given input file.
   It also crates a dummy alignment (sequences with just names) and
   puts it in the given Alignment pointer
 RETURNS: struct Cluster *
 ARGS: 
   A file handle
   A pointer to an Alignment pointer
 NOTES: 
   The file is assumed to be the distance matrix file format  used
   by the phlip package:

     4
  Name_1  0.0000   0.6776   0.6786  0.2342
  Name_2  0.6776   0.0000   0.1111  0.9999
  Name_3  0.6786   0.1111   0.0000  0.4444
  Name_4  0.2342   0.9999   0.4444  0.0000
 *********************************************************************/
struct DistanceMatrix *read_phylip_DistanceMatrix( FILE *, struct Alignment **);

/********************************************************************* 
 FUNCTION: write_phylip_DistanceMatrix
 DESCRIPTION: 
   This function takes the given DistanceMatrix and writes it to the
   given file handle in phylip format. The alignment is needed for the
   Sequence names
   format
 RETURNS: 
 ARGS: 
   A file handle
   A DistanceMatrix pointer (cluster.h)
   An Alignment pointer
 NOTES: 
   The file is written in the distance matrix file format used
   by the phlip package:

     4
  Name_1  0.0000   0.6776   0.6786  0.2342
  Name_1  0.6776   0.0000   0.1111  0.9999
  Name_1  0.6786   0.1111   0.0000  0.4444
  Name_1  0.2342   0.9999   0.4444  0.0000
*********************************************************************/

void write_phylip_DistanceMatrix( FILE *, struct DistanceMatrix *, struct Alignment *);


#endif
