/***********************************************************************/
/*                                                                     */
/* File    : global.h                                                  */
/* Author  : Torsten Rupp                                              */
/* Contents: global definitions                                        */
/* Systems : Linux                                                     */
/*                                                                     */
/***********************************************************************/

/****************************** Includes *******************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#include "global.h"

/****************** Conditional compilation switches *******************/

/***************************** Constants *******************************/

/**************************** Datatypes ********************************/

/**************************** Variables ********************************/

/****************************** Macros *********************************/

/**************************** Functions ********************************/

#ifdef __cplusplus
extern "C" {
#endif

void __halt(const char   *filename,
            unsigned int lineNb,
            int          exitcode,
            const char   *format,
            ...
           )
{
  va_list arguments;

  assert(filename != NULL);
  assert(format != NULL);

  va_start(arguments,format);
  vfprintf(stderr,format,arguments);
  va_end(arguments);
  fprintf(stderr," - halt in file %s, line %d\n", filename, lineNb);
  exit(exitcode);
}

void __abort(const char   *filename,
             unsigned int lineNb,
             const char   *prefix,
             const char   *format,
             ...
            )
{
  va_list arguments;

  assert(filename != NULL);
  assert(format != NULL);

  if (prefix != NULL) fprintf(stderr,"%s", prefix);
  va_start(arguments,format);
  vfprintf(stderr,format,arguments);
  va_end(arguments);
  fprintf(stderr," - program aborted in file %s, line %d\n", filename, lineNb);
  abort();
}

#ifdef __cplusplus
}
#endif

/* end of file */
