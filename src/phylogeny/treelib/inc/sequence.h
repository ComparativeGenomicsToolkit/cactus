/*  Last edited: Nov 15 12:24 1999 (klh) */
/**********************************************************************
 ** FILE: sequence.h
 ** NOTES:
 **   Functions and structures for the manipulation of protein 
 **   sequences. Only minimal functionality is needed
 **********************************************************************/

#ifndef _SEQUENCE
#define _SEQUENCE

#include <string.h>
#include "util.h"


#define MAX_NAME_LENGTH 25

/******************* structure definitions ****************************/

struct Sequence {
  unsigned int length;
  char *name;
  char *seq;
  char *sec_struct;
  char *surf_acc;
  char *trans_mem;
  char *post_prob;
  char *lig_bind;
};

/********************** function prototypes ***************************/


/**********************************************************************
 FUNCTION: clone_Sequence
 DESCRIPTION: 
   Performs a deep copy of the given Sequence and returns it
 ARGS:
   A Sequenced pointer (sequence.h)
 RETURNS:
   A pointer to a Sequence structure
 NOTES:
 **********************************************************************/
struct Sequence *clone_Sequence( struct Sequence *);

/**********************************************************************
 FUNCTION: empty_Sequence
 DESCRIPTION: 
   Creates and returns an new, empty Sequence object
 ARGS: 
 RETURNS:
   A pointer to a Sequence structure
 NOTES:
 **********************************************************************/
struct Sequence *empty_Sequence( void );

/**********************************************************************
 FUNCTION: free_Sequence
 DESCRIPTION: 
   frees the memory occupied by the given Sequence reference
 ARGS: 
   A pointer to a struct Sequence
 RETURNS: 
   The NULL pointer
 NOTES:
 **********************************************************************/
void *free_Sequence( struct Sequence *);


#endif
