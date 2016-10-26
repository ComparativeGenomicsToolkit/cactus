/*  Last edited: Feb  1 15:32 2002 (klh) */
/**********************************************************************
 ** FILE: util.h
 ** NOTES:
 **   This file contains general utiliy functions used throughout the 
 **   application, such as those for memory management and error
 **   messaging. I have used it as a place to put other general stuff
 **   until I have somewhere better to put it.
 **********************************************************************/

#ifndef _UTIL
#define _UTIL

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>


#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif


/********************** function prototypes ***************************/

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
void *calloc_util( size_t, size_t );


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
void *malloc_util( size_t );


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
void *realloc_util( void *, size_t );


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
void *free_util( void * );


/********************************************************************* 
 FUNCTION: fatal_util
 DESCRIPTION: 
   Prints the given formatted error message and exits
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/
void fatal_util( char *, ... );


/********************************************************************* 
 FUNCTION: warning_util
 DESCRIPTION: 
   Prints the given formatted warning to stderr
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/
void warning_util( char *, ... );


#endif
