/*  Last edited: Feb  1 18:15 2002 (klh) */
/**********************************************************************
 ** FILE: distancemat.c
 ** NOTES:
 **   Functions and types for the manipulation of Distance Matrices
 **********************************************************************/

#include "distancemat.h"


/*********************** static variables *****************************/

/* This is a table of estimated PAMs for a range of percentage-differences,
   ranging from 75% dissimilarity, going up in 0.1% steps, to 93%
   dissimilarity. For percentage dissimilarty outside this range, we either
   use Kimura's formula (for < 75) or give an arbitrarily high distance 
   (for > 93) */

static int dayhoff_pams[]={
  195,   /* 75.0% observed d; 195 PAMs estimated = 195% estimated d */
  196,   /* 75.1% observed d; 196 PAMs estimated */
                  197,    198,    199,    200,    200,    201,    202,  203,    
  204,    205,    206,    207,    208,    209,    209,    210,    211,  212,    
  213,    214,    215,    216,    217,    218,    219,    220,    221,  222,    
  223,    224,    226,    227,    228,    229,    230,    231,    232,  233,    
  234,    236,    237,    238,    239,    240,    241,    243,    244,  245,    
  246,    248,    249,    250,    /* 250 PAMs = 80.3% observed d */          
                                  252,    253,    254,    255,    257,  258,    
  260,    261,    262,    264,    265,    267,    268,    270,    271,  273,    
  274,    276,    277,    279,    281,    282,    284,    285,    287,  289,    
  291,    292,    294,    296,    298,    299,    301,    303,    305,  307,    
  309,    311,    313,    315,    317,    319,    321,    323,    325,  328,    
  330,    332,    335,    337,    339,    342,    344,    347,    349,  352,    
  354,    357,    360,    362,    365,    368,    371,    374,    377,  380,    
  383,    386,    389,    393,    396,    399,    403,    407,    410,  414,    
  418,    422,    426,    430,    434,    438,    442,    447,    451,  456,    
  461,    466,    471,    476,    482,    487,    493,    498,    504,  511,    
  517,    524,    531,    538,    545,    553,    560,    569,    577,  586,    
  595,    605,    615,    626,    637,    649,    661,    675,    688,  703,    
  719,    736,    754,    775,    796,    819,    845,    874,    907,  945,
         /* 92.9% observed; 945 PAMs */    
  988    /* 93.0% observed; 988 PAMs */
};




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

void calc_DistanceMatrix( struct DistanceMatrix *mat,
			  struct Alignment *aln,
			  unsigned int use_rand_cols,
			  unsigned int use_kimura ) {

  /* this function will take alignment and return a distance matrix */
  /* This gives a clear separation between tree making and distances */

  unsigned int i, j, k, table_index, num_undefined_distances, mem_increment;
  unsigned int *columnlist;
  Distance residuecount, distance, max_observed_distance;
  Distance **undefined_distances;

  columnlist = (unsigned int *) malloc_util( aln->length * sizeof(unsigned int));
  for (i=0; i < aln->length; i++) {
    if (use_rand_cols) {
      /* generate random column here */
      columnlist[i] = (unsigned int) rand() % aln->length;
    }
    else {
      columnlist[i] = i;
    }
  }

  max_observed_distance = 0.0;
  mem_increment = 10;
  undefined_distances = NULL;
  num_undefined_distances = 0;

  for( i=0; i < aln->numseqs; i++) {
    mat->data[i][i] = 0.0;
    for( j=0; j < i; j++ ) {
      residuecount = distance = 0.0;
      mat->data[i][j] = 0.0;

      for( k=0; k < aln->length; k++) {
	if ( aln->seqs[i]->seq[columnlist[k]] == '.' || 
	     aln->seqs[j]->seq[columnlist[k]] == '.' ||
	     aln->seqs[i]->seq[columnlist[k]] == '-' || 
	     aln->seqs[j]->seq[columnlist[k]] == '-' ||
	     aln->seqs[i]->seq[columnlist[k]] == ' ' ||
	     aln->seqs[j]->seq[columnlist[k]] == ' ')
	  continue;

	/* neither character is a gap, so proceed */
	residuecount += 1.0;
	if ( aln->seqs[i]->seq[columnlist[k]] != aln->seqs[j]->seq[columnlist[k]]) 
	  distance += 1.0;
      }

      /* if residue count was zero here, there must have been a gap in every position;
	 in this case, %identity is undefined so a decision must be made to consider the
	 sequences 100% identical or 100% different. I go with the standard approach of
	 assigning twoce the maximum observed distance */

      if (residuecount > 0) {
	distance = distance / residuecount;
      }
      else {
	distance = -1.0;
	if (num_undefined_distances % mem_increment == 0) {
	  if (num_undefined_distances == 0) {
	    undefined_distances = (Distance **) 
	      malloc_util( sizeof( Distance *) * mem_increment);
	  }
	  else {
	    undefined_distances = (Distance **) 
	      realloc_util( undefined_distances, 
			    (num_undefined_distances + mem_increment) * sizeof(Distance *));
	  }
	}
	undefined_distances[num_undefined_distances++] = &(mat->data[i][j]);
	/* distances should be zero, but may be slightly less due to floating point 
	   arithmetic */
      }

      /* Use Kimura's formula to convert percentage dissimilarity to a distance */
      
      if (use_kimura) {
	if ( distance < 0.75) {
	  if (distance > 0.0) 
	    distance = - log( 1.0 - distance - (distance * distance * 0.20) );
	}
	else {
	  if (distance > 0.930) {
    	    distance = 10.0;
	  }
	  else {
	    table_index = (int) ((distance*1000.0) - 750.0);
	    distance = (Distance) dayhoff_pams[ table_index ];
	    distance /= 100.0;
	  }
	}
      }
      if (distance > max_observed_distance) {
	max_observed_distance = distance;
      }
      mat->data[i][j] = distance;
    }
  }
  
  /* before we go, lets find those undefined distances and replace them with 
     twice the maximum observed distance */

  for (i=0; i < num_undefined_distances; i++ ) {
    *undefined_distances[i] = 2 * max_observed_distance;
  }

  columnlist = free_util( columnlist );
  if (undefined_distances != NULL) {
    undefined_distances = free_util( undefined_distances );
  }
}



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
struct DistanceMatrix *clone_DistanceMatrix( struct DistanceMatrix *source) {
  unsigned int i,j;
  struct DistanceMatrix *dest;

  if (source != NULL) {
    dest = empty_DistanceMatrix( source->size );
    
    for( i=0; i < dest->size; i++) {
      for( j=0; j <= i; j++ ) {
	dest->data[i][j] = source->data[i][j];	
      }
    }
  }
  else {
    dest = NULL;
  }

  return dest;
}



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
struct DistanceMatrix *empty_DistanceMatrix( unsigned int size) {
  unsigned int i;
  struct DistanceMatrix *mat;

  mat = (struct DistanceMatrix *) malloc_util(sizeof(struct DistanceMatrix));
  mat->size = size;
  mat->data = (Distance **) malloc_util( mat->size * sizeof(Distance *) );

  for( i=0; i < mat->size; i++)
    mat->data[i] = (Distance *) malloc_util( (i+1) * sizeof(Distance) );
 
  return mat;
}



/*********************************************************************
 FUNCTION: free_DistanceMatrix
 DESCRIPTION: 
   Frees the memory for the given distance matrix
 RETURNS:
 ARGS: 
   struct DistanceMatrix *
 NOTES: 
 *********************************************************************/

void *free_DistanceMatrix( struct DistanceMatrix *mat ) {
  int i;

  if ( mat != NULL ) {
    if (mat->data != NULL) {
      for( i=0; i < mat->size; i++ ) {
	if (mat->data[i] != NULL)
	  mat->data[i] = free_util( mat->data[i] );
      }
      mat->data = free_util( mat->data );
    }
    mat = free_util( mat );
  }

  return mat;
}




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

Distance index_DistanceMatrix( struct DistanceMatrix *mat, 
			     unsigned int i, 
			     unsigned int j) {
  if (i > j) 
    return mat->data[i][j];
  else 
    return mat->data[j][i];
}


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

void print_DistanceMatrix( FILE *handle, struct DistanceMatrix *mat ) {
  unsigned int row, column;

  fprintf( handle, "Size:%d\n", mat->size);
  
  for(row=0; row < mat->size; row++) {
    fprintf( handle, "%5d", row);
    for(column=0; column <= row; column++)
	fprintf( handle, "%10.5f", mat->data[row][column]);
    fprintf( handle, "\n");
  }
  fflush( handle );
}




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

struct DistanceMatrix *read_phylip_DistanceMatrix( FILE *handle, struct Alignment **aln_loc) {
  struct DistanceMatrix *mat;
  unsigned int size, i, j;
  char identifier[11];
  double dist;

  /* The size of the matrix will be on the first line on its own */
  if (! fscanf( handle, "%d", &size ))
    fatal_util( "Parse error: The first line should contain the size of matrix");

  *aln_loc = (struct Alignment *) malloc_util( sizeof(struct Alignment )); 

  (*aln_loc)->numseqs = size;
  (*aln_loc)->seqs = (struct Sequence **) 
    malloc_util( size * sizeof( struct Sequence *));
  (*aln_loc)->length = 0;
  mat = empty_DistanceMatrix( size );


  for (i=0; i < size; i++) {
    /* The name should be exactly 10 chars, and the scanf should place a \0
       at the end, making 11 */
    fscanf( handle, "%s", identifier );
    /* Right; the rest of the line will consist of exactly 'size' floating
       point numbers */
    (*aln_loc)->seqs[i] = empty_Sequence();
    (*aln_loc)->seqs[i]->name = (char *) malloc_util( 11 * sizeof(char));
    strcpy( (*aln_loc)->seqs[i]->name, identifier );
    for (j=0; j  < size; j++) {
      fscanf( handle, "%lf", &dist);
      if (j <= i) 
	mat->data[i][j] = (Distance) dist;
    }
  } 

  return mat;
}



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

void write_phylip_DistanceMatrix( FILE *handle, 
				  struct DistanceMatrix *mat,
				  struct Alignment *align) {

  unsigned int row, column;

  fprintf( handle, "\t%d\n", align->numseqs);
  
  for(row=0; row < align->numseqs; row++) {
    fprintf( handle, "%10.10s", align->seqs[row]->name);
    for(column=0; column < align->numseqs; column++) {
      if (row > column )
	fprintf( handle, "%10.5f", mat->data[row][column]);
      else 
	fprintf( handle, "%10.5f", mat->data[column][row]);
    }
    fprintf( handle, "\n");
  }
  fflush( handle );
}
