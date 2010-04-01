/*  Last edited: Feb  1 17:52 2002 (klh) */
/**********************************************************************
 ** FILE: util.c
 ** NOTES:
 **   This file contains general utiliy functions used throughout the 
 **   application, such as those for memory management and error
 **   messaging. I have used it as a place to put other general stuff
 **   until I have somewhere better to put it.
 **********************************************************************/

#include "util.h"



/********************************************************************* 
 FUNCTION: calloc_util
 DESCRIPTION: 
   A wrapper for the stdlib.h function calloc; exits if memory cannot
   be allocated
 RETURNS: A generic pointer to the reallocated memory
 ARGS: 
   The number of bytes to be allocated
 NOTES:
 *********************************************************************/

void *calloc_util( size_t numobjs, size_t size) {
  void *ret;
	
  if ((ret = calloc( numobjs, size )) == NULL)
    fatal_util("calloc_util: Out of memory");

  return ret;	
}



/********************************************************************* 
 FUNCTION: malloc_util
 DESCRIPTION: 
   A wrapper for the stdlib.h function malloc; exits if memory cannot
   be allocated
 RETURNS: A generic pointer to the reallocated memory
 ARGS: 
   The number of bytes to be allocated
 NOTES:
 *********************************************************************/

void *malloc_util(size_t numbytes) {
  void *ret;

  if ((ret = malloc( numbytes )) == NULL)
    fatal_util("malloc_util: out of memory when requesting %d bytes", numbytes);

  return ret;	
}



/********************************************************************* 
 FUNCTION: realloc_util
 DESCRIPTION: 
   A wrapper for the stdlib.h function realloc; exits if memory cannot
   be allocated
 RETURNS: A generic pointer to the reallocated memory
 ARGS:
   A pointer to the memory to be reallocated
   The number of bytes to be allocated
 NOTES:
 *********************************************************************/

void *realloc_util(void *ptr, size_t bytes) {
  void *ret = NULL;

  if (ptr == NULL)
    fatal_util("Call to realloc_util with a null pointer");
  else { 
    if ((ret = realloc(ptr, bytes)) == NULL)
      fatal_util("realloc_util: out of memory when requesting %d bytes", bytes);
  }

  return ret;
}  



/********************************************************************* 
 FUNCTION: free_util
 DESCRIPTION: 
   A wrapper for the stdlib.h function free; warns if the free could 
   not be performed
 RETURNS: A (hopefully null) pointer
 ARGS:
   A pointer to the memory to be deallocated
 NOTES:
 *********************************************************************/

void *free_util( void *ptr ) {
  if (ptr == NULL)
    warning_util("Call to free_util with null pointer");
  else {
    free(ptr);
    ptr = NULL;
  }
  return ptr;
}



/********************************************************************* 
 FUNCTION: fatal_util
 DESCRIPTION: 
   Prints the given formatted error message and exits
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/

void fatal_util( char *fmt, ... ) {
  va_list args;
  
  va_start( args, fmt );
  fprintf( stderr, "\nA Fatal Error occurred: ");
  vfprintf( stderr, fmt, args);
  fprintf( stderr,"\n");
  va_end( args );
  exit(1);
}



/********************************************************************* 
 FUNCTION: warning_util
 DESCRIPTION: 
   Prints the given formatted warning to stderr
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/

void warning_util( char *fmt, ...) {
  va_list args;
	
  va_start( args, fmt );
  fprintf( stderr, "\nWARNING: " );
  vfprintf( stderr, fmt, args );
  fprintf( stderr, "\n" );
  va_end( args );
}



