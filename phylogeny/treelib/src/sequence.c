/*  Last edited: Feb  1 18:21 2002 (klh) */
/**********************************************************************
 ** FILE: sequence.c
 ** NOTES:
 **   Functions and structures for the manipulation of protein 
 **   sequences. Only minimal functionality is needed
 **********************************************************************/

#include "sequence.h"


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

struct Sequence *clone_Sequence( struct Sequence *source) {
  unsigned int i;
  struct Sequence *dest = NULL;

  if (source != NULL) {
    dest = empty_Sequence();

    dest->name = (char *) malloc_util( (strlen(source->name)+1) * sizeof(char));
    strcpy( dest->name, source->name );
    dest->length = source->length;
    if (source->length) 
      dest->seq = (char *) malloc_util( dest->length * sizeof( char ));
    for( i=0; i < dest->length; i++) {
      dest->seq[i] = source->seq[i];
    }
    if (source->sec_struct != NULL) {
      dest->sec_struct = (char *) malloc_util( dest->length * sizeof( char ));
      for (i=0; i < dest->length; i++) {
	dest->sec_struct[i] = source->sec_struct[i];
      }
    }
    if (source->surf_acc != NULL) {
      dest->surf_acc = (char *) malloc_util( dest->length * sizeof( char ));
      for (i=0; i < dest->length; i++) {
	dest->surf_acc[i] = source->surf_acc[i];
      }
    }
    if (source->post_prob != NULL) {
      dest->post_prob = (char *) malloc_util( dest->length * sizeof( char ));
      for (i=0; i < dest->length; i++) {
	dest->post_prob[i] = source->post_prob[i];
      }
    }
    if (source->lig_bind != NULL) {
      dest->lig_bind = (char *) malloc_util( dest->length * sizeof( char ));
      for (i=0; i < dest->length; i++) {
	dest->lig_bind[i] = source->lig_bind[i];
      }
    }
  }

  return dest;
}



/**********************************************************************
 FUNCTION: empty_Sequence
 DESCRIPTION: 
   Creates and returns an new, empty Sequence object
 ARGS: 
 RETURNS:
   A pointer to a Sequence structure
 NOTES:
 **********************************************************************/
struct Sequence *empty_Sequence( void ) {
  struct Sequence *theseq;

  theseq = (struct Sequence *) malloc_util( sizeof(struct Sequence));
  theseq->length = 0;
  theseq->name = NULL;
  theseq->seq = NULL;
  theseq->sec_struct = NULL;
  theseq->surf_acc = NULL;
  theseq->post_prob = NULL;
  theseq->lig_bind = NULL;

  return theseq;
}




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
void *free_Sequence( struct Sequence *theseq ) {

  if (theseq != NULL) {
    if (theseq->name != NULL) 
      theseq->name = free_util( theseq->name);
    if (theseq->seq != NULL) 
      theseq->seq = free_util( theseq->seq);
    if (theseq->sec_struct != NULL)
      theseq->sec_struct = free_util( theseq->sec_struct);
    if (theseq->surf_acc != NULL)
      theseq->surf_acc = free_util( theseq->surf_acc);
    if (theseq->post_prob != NULL)
      theseq->post_prob = free_util( theseq->post_prob);
    if (theseq->lig_bind != NULL)
      theseq->lig_bind = free_util( theseq->lig_bind);

    theseq = free_util( theseq );
  }
  return theseq;
}






