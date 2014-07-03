/***********************************************************************\
*
* $Source: /home/torsten/cvs/bar/bar/strings.c,v $
* $Revision: 1.7 $
* $Author: torsten $
* Contents: dynamic string functions
* Systems: all
*
\***********************************************************************/

/****************************** Includes *******************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <sys/types.h>
#include <regex.h>
#include <assert.h>

#include "global.h"
#ifndef NDEBUG
  #include <pthread.h>
  #include "lists.h"
#endif /* not NDEBUG */

#include "treelib_strings.h"

/****************** Conditional compilation switches *******************/
#define HALT_ON_INSUFFICIENT_MEMORY
#define _FILL_MEMORY

/***************************** Constants *******************************/
#define STRING_START_LENGTH 64   // string start length
#define STRING_DELTA_LENGTH 32   // string delta increasing/decreasing

LOCAL const char *DEFAULT_TRUE_STRINGS[] =
{
  "1",
  "true",
  "yes",
  "on",
};
LOCAL const char *DEFAULT_FALSE_STRINGS[] =
{
  "0",
  "false",
  "no",
  "off",
};

#define DEBUG_FILL_BYTE 0xFE

/***************************** Datatypes *******************************/
typedef enum
{
  FORMAT_LENGTH_TYPE_INTEGER,
  FORMAT_LENGTH_TYPE_LONG,
  FORMAT_LENGTH_TYPE_LONGLONG,
  FORMAT_LENGTH_TYPE_DOUBLE,
  FORMAT_LENGTH_TYPE_QUAD,
  FORMAT_LENGTH_TYPE_POINTER,
} FormatLengthType;

typedef struct
{
  char             token[16];
  unsigned int     length;
  bool             alternateFlag;
  bool             zeroPaddingFlag;
  bool             leftAdjustedFlag;
  bool             blankFlag;
  bool             signFlag;
  unsigned int     width;
  unsigned int     precision;
  FormatLengthType lengthType;
  char             quoteChar;
  char             conversionChar;
} FormatToken;

struct __String
{
  ulong length;
  ulong maxLength;
  char  *data;
  #ifndef NDEBUG
    ulong checkSum;
  #endif /* not NDEBUG */
};

#ifndef NDEBUG
  typedef struct DebugStringNode
  {
    LIST_NODE_HEADER(struct DebugStringNode);

    const char      *fileName;
    ulong           lineNb;
    const char      *prevFileName;
    ulong           prevLineNb;
    struct __String *string;
  } DebugStringNode;

  typedef struct
  {
    LIST_HEADER(DebugStringNode);
  } DebugStringList;
#endif /* not NDEBUG */

/***************************** Variables *******************************/
#ifndef NDEBUG
  pthread_once_t  debugStringInitFlag = PTHREAD_ONCE_INIT;
  pthread_mutex_t debugStringLock;
  DebugStringList debugAllocStringList;
  DebugStringList debugFreeStringList;
#endif /* not NDEBUG */

/****************************** Macros *********************************/

#ifndef NDEBUG
  #define CHECK_VALID(string) \
    do { \
      if (string != NULL) \
      { \
        if (((ulong)(string)->length^(ulong)(string)->maxLength^(ulong)(string)->data) != (string)->checkSum) \
        { \
          HALT_INTERNAL_ERROR("Invalid checksum 0x%08x in string %p, length %ld (max. %ld) (expected 0x%08x)!",\
                              (string)->checkSum,\
                              string,\
                              (ulong)(string)->length^(ulong)(string)->maxLength^(ulong)(string)->data,\
                              (string)->length,\
                              (string)->maxLength\
                             ); \
        } \
      } \
    } while (0)
  #define UPDATE_VALID(string) \
    do { \
      if (string != NULL) \
      { \
        (string)->checkSum = (ulong)(string)->length^(ulong)(string)->maxLength^(ulong)(string)->data; \
        } \
    } while (0)
#else /* NDEBUG */
  #define CHECK_VALID(string) \
    do { \
    } while (0)
  #define UPDATE_VALID(string) \
    do { \
    } while (0)
#endif /* not NDEBUG */

/***************************** Forwards ********************************/

/***************************** Functions *******************************/

#ifdef __cplusplus
  extern "C" {
#endif

#ifndef NDEBUG
LOCAL void debugStringInit(void)
{
  pthread_mutex_init(&debugStringLock,NULL);
  List_init(&debugAllocStringList);
  List_init(&debugFreeStringList);
}
#endif /* not NDEBUG */

/***********************************************************************\
* Name   : allocString
* Purpose: allocate a new string
* Input  : -
* Output : -
* Return : allocated string or NULL on insufficient memory
* Notes  : -
\***********************************************************************/

LOCAL inline struct __String* allocString(void)
{
  struct __String *string;

  string = (struct __String*)malloc(sizeof(struct __String));
  if (string == NULL)
  {
    #ifdef HALT_ON_INSUFFICIENT_MEMORY
      HALT_INSUFFICIENT_MEMORY(NULL);
    #else /* not HALT_ON_INSUFFICIENT_MEMORY */
      return NULL;
    #endif /* HALT_ON_INSUFFICIENT_MEMORY */
  }
  string->data = (char*)malloc(STRING_START_LENGTH);
  if (string->data == NULL)
  {
    free(string);
    #ifdef HALT_ON_INSUFFICIENT_MEMORY
      HALT_INSUFFICIENT_MEMORY(NULL);
    #else /* not HALT_ON_INSUFFICIENT_MEMORY */
      return NULL;
    #endif /* HALT_ON_INSUFFICIENT_MEMORY */
  }

  string->length    = 0;
  string->maxLength = STRING_START_LENGTH;
  string->data[0]   = '\0';
  #ifndef NDEBUG
    #ifdef FILL_MEMORY
      memset(&string->data[1],DEBUG_FILL_BYTE,STRING_START_LENGTH-1);
    #endif /* FILL_MEMORY */
  #endif /* not NDEBUG */

  UPDATE_VALID(string);

  return string;
}

/***********************************************************************\
* Name   : deleteString
* Purpose: delete a string
* Input  : string - string to delete
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void deleteString(struct __String *string)
{
  assert(string != NULL);
  assert(string->data != NULL);

  free(string->data);
  free(string);
}

/***********************************************************************\
* Name   : allocTmpString
* Purpose: allocate a temporary string
* Input  : -
* Output : -
* Return : allocated string or NULL on insufficient memory
* Notes  : temporary strings are not included in list of allocated
*          strings
\***********************************************************************/

LOCAL inline struct __String* allocTmpString(void)
{
  return allocString();
}

/***********************************************************************\
* Name   : assignTmpString
* Purpose: assign data of string to string
* Input  : toString   - to string
*          fromString - from strong
* Output : -
* Return : -
* Notes  : fromString will become invalid after operation!
\***********************************************************************/

LOCAL inline void assignTmpString(struct __String *toString, struct __String *fromString)
{
  assert(toString != NULL);
  assert(toString->data != NULL);
  assert(fromString != NULL);
  assert(fromString->data != NULL);

  free(toString->data);

  toString->length    = fromString->length;
  toString->maxLength = fromString->maxLength;
  toString->data      = fromString->data;
  #ifndef NDEBUG
    fromString->length    = 0L;
    fromString->maxLength = 0L;
    fromString->data      = NULL;
  #endif /* not NDEBUG */

  UPDATE_VALID(toString);
}

/***********************************************************************\
* Name   : ensureStringLength
* Purpose: ensure min. length of string
* Input  : string    - string
*          newLength - new min. length of string
* Output : -
* Return : TRUE if string length is ok, FALSE on insufficient memory
* Notes  : -
\***********************************************************************/

LOCAL void ensureStringLength(struct __String *string, ulong newLength)
{
  char  *newData;
  ulong newMaxLength;

  if ((newLength + 1) > string->maxLength)
  {
    newMaxLength = ((newLength + 1) + STRING_DELTA_LENGTH - 1) & ~(STRING_DELTA_LENGTH - 1);
    assert(newMaxLength >= (newLength + 1));
    newData = realloc(string->data,newMaxLength);
//??? error message?
    if (newData == NULL)
    {
      perror("insufficient memory for allocating string - program halted");
      exit(128);
    }
    #ifndef NDEBUG
      #ifdef FILL_MEMORY
        memset(&newData[string->maxLength],DEBUG_FILL_BYTE,newMaxLength-string->maxLength);
      #endif /* FILL_MEMORY */
    #endif /* not NDEBUG */

    string->data      = newData;
    string->maxLength = newMaxLength;
  }
}

/***********************************************************************\
* Name   : parseNextFormatToken
* Purpose: parse next format token
* Input  : format - format string
* Output : formatToken - format token
* Return : next char after format specifier
* Notes  : -
\***********************************************************************/

LOCAL const char *parseNextFormatToken(const char *format, FormatToken *formatToken)
{
  #define ADD_CHAR(formatToken,ch) \
    do \
    { \
      assert(formatToken->length<sizeof(formatToken->token)); \
      formatToken->token[formatToken->length] = ch; formatToken->length++; \
    } while (0)

  assert(format != NULL);
  assert(formatToken != NULL);

  formatToken->length           = 0;
  formatToken->alternateFlag    = FALSE;
  formatToken->zeroPaddingFlag  = FALSE;
  formatToken->leftAdjustedFlag = FALSE;
  formatToken->blankFlag        = FALSE;
  formatToken->signFlag         = FALSE;
  formatToken->width            = 0;
  formatToken->precision        = 0;
  formatToken->lengthType       = FORMAT_LENGTH_TYPE_INTEGER;
  formatToken->quoteChar        = '\0';

  /* format start character */
  assert((*format) == '%');
  ADD_CHAR(formatToken,(*format));
  format++;

  /* flags */
  while (   ((*format) != '\0')
         && (   ((*format) == '#')
             || ((*format) == '0')
             || ((*format) == '-')
             || ((*format) == ' ')
             || ((*format) == '+')
            )
        )
  {
    ADD_CHAR(formatToken,(*format));
    switch (*format)
    {
      case '#': formatToken->alternateFlag    = TRUE; break;
      case '0': formatToken->zeroPaddingFlag  = TRUE; break;
      case '-': formatToken->leftAdjustedFlag = TRUE; break;
      case ' ': formatToken->blankFlag        = TRUE; break;
      case '+': formatToken->blankFlag        = TRUE; break;
      #ifndef NDEBUG
        default:
          HALT_INTERNAL_ERROR_UNHANDLED_SWITCH_CASE();
          break; /* not reached */
      #endif /* NDEBUG */
    }
    format++;
  }

  /* width, precision */
  while (   ((*format) != '\0')
         && isdigit((int)(*format))
        )
  {
    ADD_CHAR(formatToken,(*format));

    formatToken->width=formatToken->width*10+((*format)-'0');
    format++;
  }

  /* precision */
  if (((*format) != '\0') && ((*format) == '.'))
  {
    ADD_CHAR(formatToken,(*format));
    format++;
    while (isdigit((int)(*format)))
    {
      ADD_CHAR(formatToken,(*format));

      formatToken->precision=formatToken->precision*10+((*format)-'0');
      format++;
    }
  }

  /* quoting character */
  if (((*format) != '\0') && ((*(format+1) == 's') || (*((format+1)) == 'S')))
  {
    formatToken->quoteChar = (*format);
    format++;
  }

  /* length modifier */
  if ((*format) != '\0')
  {
    if      (((*format) == 'h') && (*((format+1)) == 'h'))
    {
      ADD_CHAR(formatToken,(*(format+0)));
      ADD_CHAR(formatToken,(*(format+1)));

      formatToken->lengthType = FORMAT_LENGTH_TYPE_INTEGER;
      format+=2;
    }
    else if ((*format) == 'h')
    {
      ADD_CHAR(formatToken,(*format));

      formatToken->lengthType = FORMAT_LENGTH_TYPE_INTEGER;
      format++;
    }
    else if (((*format) == 'l') && (*((format+1)) == 'l'))
    {
      ADD_CHAR(formatToken,(*(format+0)));
      ADD_CHAR(formatToken,(*(format+1)));

      formatToken->lengthType = FORMAT_LENGTH_TYPE_LONGLONG;
      format+=2;
    }
    else if ((*format) == 'l')
    {
      ADD_CHAR(formatToken,(*format));

      formatToken->lengthType = FORMAT_LENGTH_TYPE_LONG;
      format++;
    }
    else if ((*format) == 'q')
    {
      ADD_CHAR(formatToken,(*format));

      formatToken->lengthType = FORMAT_LENGTH_TYPE_QUAD;
      format++;
    }
    else if ((*format) == 'j')
    {
      ADD_CHAR(formatToken,(*format));

      formatToken->lengthType = FORMAT_LENGTH_TYPE_INTEGER;
      format++;
    }
    else if ((*format) == 'z')
    {
      ADD_CHAR(formatToken,(*format));

      formatToken->lengthType = FORMAT_LENGTH_TYPE_INTEGER;
      format++;
    }
    else if ((*format) == 't')
    {
      ADD_CHAR(formatToken,(*format));

      formatToken->lengthType = FORMAT_LENGTH_TYPE_INTEGER;
      format++;
    }
  }

  /* conversion character */
  if ((*format) != '\0')
  {
    switch (*format)
    {
      case 'S':
        ADD_CHAR(formatToken,'s');
        formatToken->conversionChar = 'S';
        break;
      default:
        ADD_CHAR(formatToken,(*format));
        formatToken->conversionChar = (*format);
        break;
    }
    format++;
  }

  ADD_CHAR(formatToken,'\0');
  
  return format;

  #undef ADD_CHAR
}

/***********************************************************************\
* Name   : formatString
* Purpose: format a string and append (like printf)
* Input  : String    - string
*          format    - format string
*          arguments - arguments
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

/***********************************************************************\
* Name   : parseString
* Purpose: parse a string (like scanf)
* Input  : String       - string
*          format       - format string
*          arguments    - arguments
*          stringQuotes - string chars or NULL
* Output : nextIndex - index of next character in string not parsed
*                      (can be NULL)
* Return : TRUE if parsing sucessful, FALSE otherwise
* Notes  : Additional conversion chars:
*            S - string
*            y - boolean
*          Not implemented conversion chars:
*            p
*            n
\***********************************************************************/

LOCAL bool parseString(const struct __String *string,
                       ulong                 index,
                       const char            *format,
                       va_list         arguments,
                       const char            *stringQuotes,
                       ulong                 *nextIndex
                      )
{
  FormatToken formatToken;
  union
  {
    int             *i;
    long int        *l;
    long long int   *ll;
    float           *f;
    double          *d;
    char            *c;
    char            *s;
    void            *p;
    bool            *b;
    struct __String *string;
  } value;
  char        buffer[64];
  ulong       z;
  const char  *stringQuote;
  bool        foundFlag;

  CHECK_VALID(string);

  while ((*format) != '\0')
  {
    /* skip white spaces in format */
    while (((*format) != '\0') && isspace(*format))
    {
      format++;
    }

    /* skip white-spaces in string */
    while ((index < string->length) && isspace(string->data[index]))
    {
      index++;
    }

    if ((*format) != '\0')
    {
      if ((*format) == '%')
      {
        /* get format token */
        format = parseNextFormatToken(format,&formatToken);

        /* parse string and store values */
        switch (formatToken.conversionChar)
        {
          case 'i':
          case 'd':
          case 'u':
            /* get data */
            z = 0;
            if ((index < string->length) && ((string->data[index] == '+') || (string->data[index] == '-')))
            {
              buffer[z] = string->data[index];
              z++;
              index++;
            }
            while (   (index < string->length)
                   && (z < sizeof(buffer)-1)
                   && isdigit(string->data[index])
                  )
            {
              buffer[z] = string->data[index];
              z++;
              index++;
            }
            buffer[z] = '\0';

            /* convert */
            if (z > 0)
            {
              switch (formatToken.lengthType)
              {
                case FORMAT_LENGTH_TYPE_INTEGER:
                  value.i = va_arg(arguments,int*);
                  if (value.i != NULL) (*value.i) = strtol(buffer,NULL,10);
                  break;
                case FORMAT_LENGTH_TYPE_LONG:
                  value.l = va_arg(arguments,long int*);
                  if (value.l != NULL) (*value.l) = strtol(buffer,NULL,10);
                  break;
                case FORMAT_LENGTH_TYPE_LONGLONG:
                  value.ll = va_arg(arguments,long long int*);
                  if (value.ll != NULL) (*value.ll) = strtoll(buffer,NULL,10);
                  break;
                default:
                  #ifndef NDEBUG
                    HALT_INTERNAL_ERROR_UNHANDLED_SWITCH_CASE();
                  #endif /* NDEBUG */
                  break; /* not reached */
              }
            }
            else
            {
              return FALSE;
            }
            break;
          case 'c':
            /* convert */
            if (index < string->length)
            {
              value.c = va_arg(arguments,char*);
              if (value.c != NULL) (*value.c) = string->data[index];
              index++;
            }
            else
            {
              return FALSE;
            }
            break;
          case 'o':
            /* get data */
            z = 0;
            while (   (index < string->length)
                   && (z < sizeof(buffer)-1)
                   && (string->data[index] >= '0')
                   && (string->data[index] <= '7')
                  )
            {
              buffer[z] = string->data[index];
              z++;
              index++;
            }
            buffer[z] = '\0';

            /* convert */
            if (z > 0)
            {
              switch (formatToken.lengthType)
              {
                case FORMAT_LENGTH_TYPE_INTEGER:
                  value.i = va_arg(arguments,int*);
                  if (value.i != NULL) (*value.i) = strtol(buffer,NULL,8);
                  break;
                case FORMAT_LENGTH_TYPE_LONG:
                  value.l = va_arg(arguments,long int*);
                  if (value.l != NULL) (*value.l) = strtol(buffer,NULL,10);
                  break;
                case FORMAT_LENGTH_TYPE_LONGLONG:
                  value.ll = va_arg(arguments,long long int*);
                  if (value.ll != NULL) (*value.ll) = strtoll(buffer,NULL,10);
                  break;
                default:
                  #ifndef NDEBUG
                    HALT_INTERNAL_ERROR_UNHANDLED_SWITCH_CASE();
                  #endif /* NDEBUG */
                  break; /* not reached */
              }
            }
            else
            {
              return FALSE;
            }
            break;
          case 'x':
          case 'X':
            /* get data */
            if (((index+1) < string->length) && (string->data[index+0] == '0') && (string->data[index+0] == 'x'))
            {
              index+=2;
            }
            z = 0;
            while (   (index < string->length)
                   && (z < sizeof(buffer)-1)
                   && isdigit(string->data[index])
                  )
            {
              buffer[z] = string->data[index];
              z++;
              index++;
            }
            buffer[z] = '\0';

            /* convert */
            if (z > 0)
            {
              switch (formatToken.lengthType)
              {
                case FORMAT_LENGTH_TYPE_INTEGER:
                  value.i = va_arg(arguments,int*);
                  if (value.i != NULL) (*value.i) = strtol(buffer,NULL,16);
                  break;
                case FORMAT_LENGTH_TYPE_LONG:
                  value.l = va_arg(arguments,long int*);
                  if (value.l != NULL) (*value.l) = strtol(buffer,NULL,16);
                  break;
                case FORMAT_LENGTH_TYPE_LONGLONG:
                  value.ll = va_arg(arguments,long long int*);
                  if (value.ll != NULL) (*value.ll) = strtol(buffer,NULL,16);
                  break;
                default:
                  #ifndef NDEBUG
                    HALT_INTERNAL_ERROR_UNHANDLED_SWITCH_CASE();
                  #endif /* NDEBUG */
                  break; /* not reached */
              }
            }
            else
            {
              return FALSE;
            }
            break;
          case 'e':
          case 'E':
          case 'f':
          case 'F':
          case 'g':
          case 'G':
          case 'a':
          case 'A':
            /* get data */
            z = 0;
            if ((index < string->length) && ((string->data[index] == '+') || (string->data[index] == '-')  || (string->data[index] == '.')))
            {
              buffer[z] = string->data[index];
              z++;
              index++;
            }
            while (   (index < string->length)
                   && (z < sizeof(buffer)-1)
                   && isdigit(string->data[index])
                  )
            {
              buffer[z] = string->data[index];
              z++;
              index++;
            }
            if ((index < string->length) && (string->data[index] == '.'))
            {
              buffer[z] = '.';
              z++;
              index++;
              while (   (index < string->length)
                     && (z < sizeof(buffer)-1)
                     && isdigit(string->data[index])
                    )
              {
                buffer[z] = string->data[index];
                z++;
                index++;
              }
            }
            buffer[z] = '\0';

            /* convert */
            if (z > 0)
            {
              switch (formatToken.lengthType)
              {
                case FORMAT_LENGTH_TYPE_INTEGER:
                  value.f = va_arg(arguments,float*);
                  if (value.f != NULL) (*value.f) = strtod(buffer,NULL);
                  break;
                case FORMAT_LENGTH_TYPE_LONG:
                  value.d = va_arg(arguments,double*);
                  if (value.d != NULL) (*value.d) = strtod(buffer,NULL);
                  break;
                case FORMAT_LENGTH_TYPE_LONGLONG:
                  HALT_INTERNAL_ERROR_STILL_NOT_IMPLEMENTED(NULL);
                  break;
                default:
                  #ifndef NDEBUG
                    HALT_INTERNAL_ERROR_UNHANDLED_SWITCH_CASE();
                  #endif /* NDEBUG */
                  break; /* not reached */
              }
            }
            else
            {
              return FALSE;
            }
            break;
          case 's':
            /* get and copy data */
            value.s = va_arg(arguments,char*);
            assert(formatToken.width > 0);

            if (index < string->length)
            {
              z = 0;
              while (   (index < string->length)
                     && (formatToken.blankFlag || !isspace(string->data[index]))
                     && (string->data[index] != (*format))
                    )
              {
                if (string->data[index] == '\\')
                {
                  index++;
                  if (index < string->length)
                  {
                    if ((formatToken.width == 0) || (z < formatToken.width-1))
                    {
                      if (value.s != NULL) value.s[z] = string->data[index];
                      z++;
                    }
                    index++;
                  }
                }
                else
                {
                  /* check for string quote */
                  stringQuote = NULL;
                  if ((formatToken.quoteChar != '\0') && (formatToken.quoteChar == string->data[index])) stringQuote = &formatToken.quoteChar;
                  if ((stringQuote == NULL) && (stringQuotes != NULL)) stringQuote = strchr(stringQuotes,string->data[index]);

                  if (stringQuote != NULL)
                  {
                    do
                    {
                      /* skip quote-char */
                      index++;

                      /* get string */
                      while ((index < string->length) && (string->data[index] != (*stringQuote)))
                      {
                        if (string->data[index] == '\\')
                        {
                          index++;
                          if (index < string->length)
                          {
                            if ((formatToken.width == 0) || (z < formatToken.width-1))
                            {
                              if (value.s != NULL) value.s[z] = string->data[index];
                              z++;
                            }
                            index++;
                          }
                        }
                        else
                        {
                          if (z < (formatToken.width-1))
                          {
                            if (value.s != NULL) value.s[z] = string->data[index];
                            z++;
                          }
                          index++;
                        }
                      }

                      /* skip quote-char */
                      if (index < string->length)
                      {
                        index++;
                      }

                      stringQuote = NULL;
                      if (index < string->length)
                      {
                        if ((formatToken.quoteChar != '\0') && (formatToken.quoteChar == string->data[index])) stringQuote = &formatToken.quoteChar;
                        if ((stringQuote == NULL) && (stringQuotes != NULL)) stringQuote = strchr(stringQuotes,string->data[index]);
                      }
                    }
                    while (stringQuote != NULL);
                  }
                  else
                  {
                    if (z < (formatToken.width-1))
                    {
                      if (value.s != NULL) value.s[z] = string->data[index];
                      z++;
                    }
                    index++;
                  }
                }
              }
              if (value.s != NULL) value.s[z] = '\0';
            }
            else
            {
              return FALSE;
            }
            break;
          case 'p':
          case 'n':
            break;
          case 'S':
            /* get and copy data */
            value.string = va_arg(arguments,String);
            CHECK_VALID(value.string);

            if (index < string->length)
            {
              String_clear(value.string);           
              z = 0;
              while (   (index < string->length)
                     && (formatToken.blankFlag || !isspace(string->data[index]))
    // NUL in string here a problem?
                     && (string->data[index] != (*format))
                    )
              {
                if (string->data[index] == '\\')
                {
                  index++;
                  if (index < string->length)
                  {
                    if ((formatToken.width == 0) || (z < formatToken.width-1))
                    {
                      String_appendChar(value.string,string->data[index]);
                      z++;
                    }
                    index++;
                  }
                }
                else
                {
                  /* check for string quote */
                  stringQuote = NULL;
                  if ((formatToken.quoteChar != '\0') && (formatToken.quoteChar == string->data[index])) stringQuote = &formatToken.quoteChar;
                  if ((stringQuote == NULL) && (stringQuotes != NULL)) stringQuote = strchr(stringQuotes,string->data[index]);

                  if (stringQuote != NULL)
                  {
                    do
                    {
                      /* skip quote-char */
                      index++;

                      /* get string */
                      while ((index < string->length) && (string->data[index] != (*stringQuote)))
                      {
                        if (string->data[index] == '\\')
                        {
                          index++;
                          if (index < string->length)
                          {
                            if ((formatToken.width == 0) || (z < formatToken.width-1))
                            {
                              String_appendChar(value.string,string->data[index]);
                              z++;
                            }
                            index++;
                          }
                        }
                        else
                        {
                          if ((formatToken.width == 0) || (z < formatToken.width-1))
                          {
                            String_appendChar(value.string,string->data[index]);
                            z++;
                          }
                          index++;
                        }
                      }

                      /* skip quote-char */
                      if (index < string->length)
                      {
                        index++;
                      }

                      /* check for string quote */
                      stringQuote = NULL;
                      if (index < string->length)
                      {
                        if ((formatToken.quoteChar != '\0') && (formatToken.quoteChar == string->data[index])) stringQuote = &formatToken.quoteChar;
                        if ((stringQuote == NULL) && (stringQuotes != NULL)) stringQuote = strchr(stringQuotes,string->data[index]);
                      }
                    }
                    while (stringQuote != NULL);
                  }
                  else
                  {
                    if ((formatToken.width == 0) || (z < formatToken.width-1))
                    {
                      String_appendChar(value.string,string->data[index]);
                      z++;
                    }
                    index++;
                  }
                }
              }
            }
            else
            {
              return FALSE;
            }
            break;
  #if 0
  still not implemented
          case 'b':
            /* binaray value */
            switch (formatToken.lengthType)
            {
              case FORMAT_LENGTH_TYPE_INTEGER:
                {
                  unsigned int bits;

                  bits = va_arg(arguments,unsigned int);
                }
                break;
              case FORMAT_LENGTH_TYPE_LONG:
                {
                  unsigned long bits;

                  bits = va_arg(arguments,unsigned long);
                }
                break;
              case FORMAT_LENGTH_TYPE_LONGLONG:
                {
                  #if defined(_LONG_LONG) || defined(HAVE_LONG_LONG)
                    unsigned long long bits;

                    bits = va_arg(arguments,unsigned long long);
                  }
                  #else /* not _LONG_LONG || HAVE_LONG_LONG */
                    HALT_INTERNAL_ERROR("long long not supported");
                  #endif /* _LONG_LONG || HAVE_LONG_LONG */
                break;
              #ifndef NDEBUG
                default:
                  HALT_INTERNAL_ERROR_UNHANDLED_SWITCH_CASE();
                  break; /* not reached */
              #endif /* NDEBUG */
            }
            bufferCount = 0;
            FLUSH_OUTPUT();
  HALT_NOT_YET_IMPLEMENTED();
            break;
  #endif /* 0 */
          case 'y':
            /* get data */
            z = 0;
            while (   (index < string->length)
                   && !isspace(string->data[index])
                  )
            {
              if (z < sizeof(buffer)-1)
              {
                buffer[z] = string->data[index];
                z++;
              }
              index++;
            }
            buffer[z] = '\0';

            /* convert */
            if (z > 0)
            {
              value.b = va_arg(arguments,bool*);
              foundFlag = FALSE;
              z = 0;
              while (!foundFlag && (z < SIZE_OF_ARRAY(DEFAULT_TRUE_STRINGS)))
              {
                if (strcmp(buffer,DEFAULT_TRUE_STRINGS[z]) == 0)
                {
                  if (value.b != NULL) (*value.b) = TRUE;
                  foundFlag = TRUE;
                }
                z++;
              }
              z = 0;
              while (!foundFlag && (z < SIZE_OF_ARRAY(DEFAULT_FALSE_STRINGS)))
              {
                if (strcmp(buffer,DEFAULT_FALSE_STRINGS[z]) == 0)
                {
                  if (value.b != NULL) (*value.b) = FALSE;
                  foundFlag = TRUE;
                }
                z++;
              }

              if (!foundFlag)
              {
                return FALSE;
              }
            }
            else
            {
              return FALSE;
            }
            break;
          case '*':
            /* skip value */
            z = 0;
            while (   (index < string->length)
                   && !isspace(string->data[index])
                   && (string->data[index] != (*format))
                  )
            {
              index++;
            }
            break;
          case '%':
            if ((index >= string->length) || (string->data[index] != '%'))
            {
              return FALSE;
            }
            index++;
            break;
          default:
            HALT_INTERNAL_ERROR_UNHANDLED_SWITCH_CASE(NULL);
            break;
        }
      }
      else
      {
        if ((index >= string->length) || (string->data[index] != (*format)))
        {
          return FALSE;
        }
        index++;
        format++;
      }
    }
  }
  if (nextIndex != NULL)
  {
    (*nextIndex) = (index >= string->length)?STRING_END:index;
  }
  else
  {
    if (index < string->length)
    {
      return FALSE;
    }
  }

  return TRUE;
}

/***********************************************************************\
* Name   : getUnitFactor
* Purpose: get unit factor
* Input  : stringUnits     - string units
*          stringUnitCount - number of string units
*          string          - string start
*          unitString      - unit string start
* Output : nextIndex - index of next character in string not parsed or
*                      STRING_END if string completely parsed (can be
*                      NULL)
* Return : -
* Notes  : -
\***********************************************************************/

LOCAL ulong getUnitFactor(const StringUnit stringUnits[],
                          uint             stringUnitCount,
                          const char       *string,
                          const char       *unitString,
                          long             *nextIndex
                         )
{
  uint  z;
  ulong factor;

  z = 0;
  while ((z < stringUnitCount) && (strcmp(unitString,stringUnits[z].name) != 0))
  {
    z++;
  }
  if (z < stringUnitCount)
  {
    factor = stringUnits[z].factor;
    if (nextIndex != NULL) (*nextIndex) = STRING_END;
  }
  else
  {
    factor = 1L;
    if (nextIndex != NULL) (*nextIndex) = (ulong)(unitString-string);
  }

  return factor;
}

/*---------------------------------------------------------------------*/

#ifdef NDEBUG
String String_new(void)
#else /* not DEBUG */
String __String_new(const char *fileName, ulong lineNb)
#endif /* NDEBUG */
{
  struct __String *string;
  #ifndef NDEBUG
    DebugStringNode *debugStringNode;;
  #endif /* not NDEBUG */

  string = allocString();
  if (string == NULL)
  {
    return NULL;
  }

  #ifndef NDEBUG
    pthread_once(&debugStringInitFlag,debugStringInit);

    pthread_mutex_lock(&debugStringLock);
    debugStringNode = debugFreeStringList.head;
    while ((debugStringNode != NULL) && (debugStringNode->string != string))
    {
      debugStringNode = debugStringNode->next;
    }
    if (debugStringNode != NULL)
    {
      List_remove(&debugFreeStringList,debugStringNode);
    }
    else
    {
      debugStringNode = LIST_NEW_NODE(DebugStringNode);
      if (debugStringNode == NULL)
      {
        HALT_INSUFFICIENT_MEMORY();
      }
    }
    debugStringNode->fileName     = fileName;
    debugStringNode->lineNb       = lineNb;
    debugStringNode->prevFileName = NULL;
    debugStringNode->prevLineNb   = 0;
    debugStringNode->string   = string;
    List_append(&debugAllocStringList,debugStringNode);
    pthread_mutex_unlock(&debugStringLock);
  #endif /* not NDEBUG */

  UPDATE_VALID(string);

  return string;
}

#ifdef NDEBUG
String String_newCString(const char *s)
#else /* not NDEBUG */
String __String_newCString(const char *fileName, ulong lineNb, const char *s)
#endif /* NDEBUG */
{
  String string;

  #ifdef NDEBUG
    string = String_new();
  #else /* not DEBUG */
    string = __String_new(fileName,lineNb);
  #endif /* NDEBUG */
  String_setCString(string,s);

  UPDATE_VALID(string);

  return string;
}

#ifdef NDEBUG
String String_newChar(char ch)
#else /* not NDEBUG */
String __String_newChar(const char *fileName, ulong lineNb, char ch)
#endif /* NDEBUG */
{
  String string;

  #ifdef NDEBUG
    string = String_new();
  #else /* not DEBUG */
    string = __String_new(fileName,lineNb);
  #endif /* NDEBUG */
  String_setChar(string,ch);

  UPDATE_VALID(string);

  return string;
}

#ifdef NDEBUG
String String_newBuffer(const void *buffer, ulong bufferLength)
#else /* not NDEBUG */
String __String_newBuffer(const char *fileName, ulong lineNb, const void *buffer, ulong bufferLength)
#endif /* NDEBUG */
{
  String string;

  #ifdef NDEBUG
    string = String_new();
  #else /* not DEBUG */
    string = __String_new(fileName,lineNb);
  #endif /* NDEBUG */
  String_setBuffer(string,buffer,bufferLength);

  UPDATE_VALID(string);

  return string;
}

#ifdef NDEBUG
String String_duplicate(const String fromString)
#else /* not NDEBUG */
String __String_duplicate(const char *fileName, ulong lineNb, const String fromString)
#endif /* NDEBUG */
{
  struct __String *string;

  CHECK_VALID(fromString);

  if (fromString != NULL)
  {
    assert(fromString->data != NULL);

    #ifdef NDEBUG
      string = String_new();
    #else /* not DEBUG */
      string = __String_new(fileName,lineNb);
    #endif /* NDEBUG */
    if (string == NULL)
    {
      return NULL;
    }

    ensureStringLength(string,fromString->length);
    memcpy(&string->data[0],&fromString->data[0],fromString->length);
    string->data[fromString->length] ='\0';
    string->length = fromString->length;

    UPDATE_VALID(string);
  }
  else
  {
    string = NULL;
  }

  return string;
}

#ifdef NDEBUG
String String_copy(String *string, const String fromString)
#else /* not NDEBUG */
String __String_copy(const char *fileName, ulong lineNb, String *string, const String fromString)
#endif /* NDEBUG */
{
  CHECK_VALID(fromString);

  if (fromString != NULL)
  {
    assert(fromString->data != NULL);

    if ((*string) == NULL)
    {
      #ifdef NDEBUG
        (*string) = String_new();
      #else /* not DEBUG */
        (*string) = __String_new(fileName,lineNb);
      #endif /* NDEBUG */
      if ((*string) == NULL)
      {
        return NULL;
      }
    }

    ensureStringLength((*string),fromString->length);
    memcpy(&(*string)->data[0],&fromString->data[0],fromString->length);
    (*string)->data[fromString->length] ='\0';
    (*string)->length = fromString->length;

    UPDATE_VALID(*string);
  }
  else
  {
    if ((*string) != NULL)
    {
      (*string)->data[0] ='\0';
      (*string)->length = 0;

      UPDATE_VALID(*string);
    }
  }

  return (*string);
}

#ifdef NDEBUG
void String_delete(String string)
#else /* not NDEBUG */
void __String_delete(const char *fileName, ulong lineNb, String string)
#endif /* NDEBUG */
{
  #ifndef NDEBUG
    DebugStringNode *debugStringNode;;
  #endif /* not NDEBUG */

  CHECK_VALID(string);

  #ifndef NDEBUG
    pthread_mutex_lock(&debugStringLock);
    debugStringNode = debugFreeStringList.head;
    while ((debugStringNode != NULL) && (debugStringNode->string != string))
    {
      debugStringNode = debugStringNode->next;
    }
    if (debugStringNode != NULL)
    {
      fprintf(stderr,"DEBUG WARNING: multiple free of string %p at %s, %ld and previoulsy at %s, %ld which was allocated at %s, %ld!\n",
              string,
              fileName,
              lineNb,
              debugStringNode->prevFileName,
              debugStringNode->prevLineNb,
              debugStringNode->fileName,
              debugStringNode->lineNb
             );
      HALT_INTERNAL_ERROR("");
    }
    pthread_mutex_unlock(&debugStringLock);
  #endif /* not NDEBUG */

  if (string != NULL)
  {
    assert(string->data != NULL);

    #ifndef NDEBUG
      pthread_once(&debugStringInitFlag,debugStringInit);

      pthread_mutex_lock(&debugStringLock);
      debugStringNode = debugAllocStringList.head;
      while ((debugStringNode != NULL) && (debugStringNode->string != string))
      {
        debugStringNode = debugStringNode->next;
      }
      if (debugStringNode != NULL)
      {
        List_remove(&debugAllocStringList,debugStringNode);
        debugStringNode->prevFileName = fileName;
        debugStringNode->prevLineNb   = lineNb;
        List_append(&debugFreeStringList,debugStringNode);
      }
      else
      {
        fprintf(stderr,"DEBUG WARNING: string '%s' not found in debug list at %s, line %ld\n",
                string->data,
                fileName,
                lineNb
               );
        HALT_INTERNAL_ERROR("");
      }
      pthread_mutex_unlock(&debugStringLock);
    #endif /* not NDEBUG */

    free(string->data);
    free(string);
  }
}

String String_clear(String string)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    string->data[0] = '\0';
    string->length = 0;

    UPDATE_VALID(string);
  }

  return string;
}

String String_erase(String string)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    memset(string->data,0,string->maxLength);
    string->length = 0;

    UPDATE_VALID(string);
  }

  return string;
}

String String_set(String string, const String sourceString)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (sourceString != NULL)
    {
      assert(sourceString->data != NULL);

      ensureStringLength(string,sourceString->length);
      memcpy(&string->data[0],&sourceString->data[0],sourceString->length);
      string->data[sourceString->length] = '\0';
      string->length = sourceString->length;
    }
    else
    {
      string->data[0] = '\0';
      string->length = 0;
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_setCString(String string, const char *s)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (s != NULL)
    {
      String_setBuffer(string,s,strlen(s));
    }
    else
    {
      string->data[0] = '\0';
      string->length = 0;
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_setChar(String string, char ch)
{
  CHECK_VALID(string);

  String_setBuffer(string,&ch,1); 

  UPDATE_VALID(string);

  return string;
}

String String_setBuffer(String string, const void *buffer, ulong bufferLength)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (buffer != NULL)
    {
      assert(buffer != NULL);

      ensureStringLength(string,bufferLength);
      memcpy(&string->data[0],buffer,bufferLength);
      string->data[bufferLength] = '\0';
      string->length = bufferLength;
    }
    else
    {
      string->data[0] = '\0';
      string->length = 0;
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_sub(String string, const String fromString, ulong fromIndex, long fromLength)
{
  ulong n;

  CHECK_VALID(string);
  CHECK_VALID(fromString);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (fromString != NULL)
    {
      assert(fromString->data != NULL);

      if (fromIndex < fromString->length)
      {
        if (fromIndex == STRING_END)
        {
          n = MIN(fromString->length,fromString->length-fromIndex);
        }
        else
        {
          n = MIN(fromLength,fromString->length-fromIndex);
        }
        ensureStringLength(string,n);
        memcpy(&string->data[0],&fromString->data[fromIndex],n);
        string->data[n] ='\0';
        string->length = n;
      }
    }
    else
    {
      string->data[0] = '\0';
      string->length = 0;
    }

    UPDATE_VALID(string);
  }

  return string;
}

char *String_subCString(char *s, const String fromString, ulong fromIndex, long fromLength)
{
  ulong n;

  assert(s != NULL);

  CHECK_VALID(fromString);

  if (fromLength > 0)
  {
    if (fromString != NULL)
    {
      assert(fromString->data != NULL);

      if (fromIndex < fromString->length)
      {
        n = MIN(fromLength,fromString->length-fromIndex);
        memcpy(s,&fromString->data[fromIndex],n);
        s[n] = '\0';
      }
      else
      {
        s[0] = '\0';
      }
    }
    else
    {
      s[0] = '\0';
    }
  }

  return s;
}

char *String_subBuffer(char *buffer, const String fromString, ulong fromIndex, long fromLength)
{
  ulong n;

  assert(buffer != NULL);

  CHECK_VALID(fromString);

  if (fromLength > 0)
  {
    if (fromString != NULL)
    {
      assert(fromString->data != NULL);

      if (fromIndex < fromString->length)
      {
        n = MIN(fromLength,fromString->length-fromIndex);
        memcpy(&buffer[0],&fromString->data[fromIndex],n);
        memset(&buffer[n],0,fromLength-n);
      }
      else
      {
        memset(buffer,0,fromLength);
      }
    }
    else
    {
      memset(buffer,0,fromLength);
    }
  }

  return buffer;
}

String String_append(String string, const String appendString)
{
  ulong n;

  CHECK_VALID(string);
  CHECK_VALID(appendString);

  if (string != NULL)
  {
    assert(string->data != NULL);
    if (appendString != NULL)
    {
      n = string->length+appendString->length;
      ensureStringLength(string,n);
      memcpy(&string->data[string->length],&appendString->data[0],appendString->length);
      string->data[n] = '\0';
      string->length = n;
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_appendSub(String string, const String fromString, ulong fromIndex, long fromLength)
{
  ulong n;

  CHECK_VALID(string);
  CHECK_VALID(fromString);

  if (string != NULL)
  {
    assert(string->data != NULL);
    if (fromString != NULL)
    {
      assert(fromString->data != NULL);

      if (fromIndex < fromString->length)
      {
        if (fromIndex == STRING_END)
        {
          n = MIN(fromString->length,fromString->length-fromIndex);
        }
        else
        {
          n = MIN(fromLength,fromString->length-fromIndex);
        }
        ensureStringLength(string,string->length+n);
        memcpy(&string->data[string->length],&fromString->data[fromIndex],n);
        string->data[string->length+n] ='\0';
        string->length += n;
      }
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_appendBuffer(String string, const char *buffer, ulong bufferLength)
{
  ulong n;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);
    if (buffer != NULL)
    {
      n = string->length+bufferLength;
      ensureStringLength(string,n);
      memcpy(&string->data[string->length],buffer,bufferLength);
      string->data[n] = '\0';
      string->length = n;
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_appendCString(String string, const char *s)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (s != NULL)
    {
      String_appendBuffer(string,s,strlen(s));
    }
  }

  return string;
}

String String_appendChar(String string, char ch)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    String_appendBuffer(string,&ch,1);
  }

  return string;
}

String String_insert(String string, ulong index, const String insertString)
{
  ulong n;

  CHECK_VALID(string);
  CHECK_VALID(insertString);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (insertString != NULL)
    {
      if      (index == STRING_END)
      {
        n = string->length+insertString->length;
        ensureStringLength(string,n);
        memcpy(&string->data[string->length],&insertString->data[0],insertString->length);
        string->data[n] = '\0';
        string->length = n;
      }
      else if (index <= string->length)
      {
        n = string->length+insertString->length;
        ensureStringLength(string,n);
        memmove(&string->data[index+insertString->length],&string->data[index],string->length-index);
        memcpy(&string->data[index],&insertString->data[0],insertString->length);
        string->data[n] = '\0';
        string->length = n;
      }
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_insertSub(String string, ulong index, const String fromString, ulong fromIndex, long fromLength)
{
  ulong n;

  CHECK_VALID(string);
  CHECK_VALID(fromString);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (fromString != NULL)
    {
      if (fromIndex < fromString->length)
      {
        if (fromIndex == STRING_END)
        {
          n = MIN(fromString->length,fromString->length-fromIndex);
        }
        else
        {
          n = MIN(fromLength,fromString->length-fromIndex);
        }

        if      (index == STRING_END)
        {
          ensureStringLength(string,string->length+n);
          memcpy(&string->data[string->length],&fromString->data[fromIndex],n);
          string->data[string->length+n] = '\0';
          string->length += n;
        }
        else if (index <= string->length)
        {
          ensureStringLength(string,string->length+n);
          memmove(&string->data[index+n],&string->data[index],string->length-index);
          memcpy(&string->data[index],&fromString->data[fromIndex],n);
          string->data[string->length+n] = '\0';
          string->length += n;
        }
      }
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_insertBuffer(String string, ulong index, const char *buffer, ulong bufferLength)
{
  ulong n;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (buffer != NULL)
    {
      if      (index == STRING_END)
      {
        n = string->length+bufferLength;
        ensureStringLength(string,n);
        memcpy(&string->data[string->length],buffer,bufferLength);
        string->data[n] = '\0';
        string->length = n;
      }
      else if (index <= string->length)
      {
        n = string->length+bufferLength;
        ensureStringLength(string,n);
        memmove(&string->data[index+bufferLength],&string->data[index],string->length-index);
        memcpy(&string->data[index],buffer,bufferLength);
        string->data[n] = '\0';
        string->length = n;
      }
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_insertCString(String string, ulong index, const char *s)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (s != NULL)
    {
      String_insertBuffer(string,index,s,strlen(s));
    }
  }

  return string;
}

String String_insertChar(String string, ulong index, char ch)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    String_insertBuffer(string,index,&ch,1);
  }

  return string;
}

String String_remove(String string, ulong index, ulong length)
{
  ulong n;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if      (index == STRING_END)
    {
      n = (string->length > length)?string->length-length:0;
      string->data[n] = '\0';
      string->length = n;
    }
    else if (index < string->length)
    {
      if ((index + length) < string->length)
      {
        memmove(&string->data[index],&string->data[index + length],string->length - (index + length));
        n = string->length - length;
      }
      else
      {
        n = index;
      }
      string->data[n] = '\0';
      string->length = n;
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_replace(String string, ulong index, ulong length, const String insertString)
{
  CHECK_VALID(string);
  CHECK_VALID(insertString);

  String_remove(string,index,length);
  String_insert(string,index,insertString);

  return string;
}

String String_replaceBuffer(String string, ulong index, ulong length, const char *buffer, ulong bufferLength)
{
  CHECK_VALID(string);

  String_remove(string,index,length);
  String_insertBuffer(string,index,buffer,bufferLength);

  return string;
}

String String_replaceCString(String string, ulong index, ulong length, const char *s)
{
  CHECK_VALID(string);

  String_remove(string,index,length);
  String_insertCString(string,index,s);

  return string;
}

String String_replaceChar(String string, ulong index, ulong length, char ch)
{
  CHECK_VALID(string);

  String_remove(string,index,length);
  String_insertChar(string,index,ch);

  return string;
}

ulong String_length(const String string)
{
  CHECK_VALID(string);

  return (string != NULL)?string->length:0;
}

bool String_empty(const String string)
{
  CHECK_VALID(string);

  return (string != NULL)?(string->length == 0):TRUE;
}

char String_index(const String string, ulong index)
{
  char ch;

  CHECK_VALID(string);

  if (string != NULL)
  {
    if      (index == STRING_END)
    {
      ch = (string->length > 0)?string->data[string->length-1]:'\0';
    }
    else if (index < string->length)
    {
      ch = string->data[index];
    }
    else
    {
      ch = '\0';
    }
  }
  else
  {
    ch = '\0';
  }

  return ch;
}

const char *String_cString(const String string)
{
  CHECK_VALID(string);

  return (string != NULL)?&string->data[0]:NULL;
}

int String_compare(const String          string1,
                   const String          string2,
                   StringCompareFunction stringCompareFunction,
                   void                  *stringCompareUserData
                  )
{
  ulong n;
  ulong z;
  int   result;

  assert(string1 != NULL);
  assert(string2 != NULL);

  CHECK_VALID(string1);
  CHECK_VALID(string2);

  result = 0;
  n = MIN(string1->length,string2->length);
  z = 0;
  if (stringCompareFunction != NULL)
  {
    while ((result == 0) && (z < n))
    {
      result = stringCompareFunction(stringCompareUserData,string1->data[z],string2->data[z]);
      z++;
    }
  }
  else
  {
    while ((result == 0) && (z < n))
    {
      if      (string1->data[z] < string2->data[z]) result = -1;
      else if (string1->data[z] > string2->data[z]) result =  1;
      z++;
    }
  }
  if (result == 0)
  {
    if      (string1->length < string2->length) result = -1;
    else if (string1->length > string2->length) result =  1;
  }

  return result;
}

bool String_equals(const String string1, const String string2)
{
  bool  equalFlag;
  ulong z;

  if ((string1 != NULL) && (string2 != NULL))
  {
    CHECK_VALID(string1);
    CHECK_VALID(string2);

    if (string1->length == string2->length)
    {
      equalFlag = TRUE;
      z         = 0;
      while (equalFlag && (z < string1->length))
      {
        equalFlag = (string1->data[z] == string2->data[z]);
        z++;
      }
    }
    else
    {
      equalFlag = FALSE;
    }
  }
  else
  {
    equalFlag = (string1 == NULL) && (string2 == NULL);
  }

  return equalFlag;
}

bool String_equalsBuffer(const String string, const char *buffer, ulong bufferLength)
{
  bool  equalFlag;
  ulong z;

  assert(string != NULL);
  assert(buffer != NULL);

  if ((string != NULL) && (buffer != NULL))
  {
    CHECK_VALID(string);

    if (string->length == bufferLength)
    {
      equalFlag = TRUE;
      z         = 0;
      while (equalFlag && (z < string->length))
      {
        equalFlag = (string->data[z] == buffer[z]);
        z++;
      }
    }
    else
    {
      equalFlag = FALSE;
    }
  }
  else
  {
    equalFlag = (string == NULL) && (bufferLength == 0);
  }

  return equalFlag;
}

bool String_equalsCString(const String string, const char *s)
{
  bool equalFlag;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (s != NULL)
    {
      equalFlag = String_equalsBuffer(string,s,strlen(s));
    }
    else
    {
      equalFlag = (string->length == 0);
    }
  }
  else
  {
    equalFlag = (s == NULL);
  }

  return equalFlag;
}

bool String_equalsChar(const String string, char ch)
{
  bool equalFlag;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    equalFlag = ((string->length == 1) && (string->data[0] == ch));
  }
  else
  {
    equalFlag = FALSE;
  }

  return equalFlag;
}

bool String_subEquals(const String string1, const String string2, long index, ulong length)
{
  long  i;
  bool  equalFlag;
  ulong z;

  assert(string1 != NULL);
  assert(string2 != NULL);

  if ((string1 != NULL) && (string2 != NULL))
  {
    CHECK_VALID(string1);
    CHECK_VALID(string2);

    i = (index != STRING_END)?index:(long)string1->length-(long)length;
    if (   (i >= 0)
        && ((i+length) <= string1->length)
        && (length <= string2->length)
       )
    {
      equalFlag = TRUE;
      z         = 0;
      while (equalFlag && (z < length))
      {
        equalFlag = (string1->data[i+z] == string2->data[z]);
        z++;
      }
    }
    else
    {
      equalFlag = FALSE;
    }
  }
  else
  {
    equalFlag = (string1 == NULL) && (string2 == NULL);
  }

  return equalFlag;
}

bool String_subEqualsBuffer(const String string, const char *buffer, ulong bufferLength, long index, ulong length)
{
  long  i;
  bool  equalFlag;
  ulong z;

  assert(string != NULL);
  assert(buffer != NULL);

  if ((string != NULL) && (buffer != NULL))
  {
    CHECK_VALID(string);

    i = (index != STRING_END)?index:(long)string->length-(long)length;
    if (   (i >= 0)
        && ((i+length) <= string->length)
        && (length <= bufferLength)
       )
    {
      equalFlag = TRUE;
      z         = 0;
      while (equalFlag && (z < length))
      {
        equalFlag = (string->data[i+z] == buffer[z]);
        z++;
      }
    }
    else
    {
      equalFlag = FALSE;
    }
  }
  else
  {
    equalFlag = (string == NULL) && (bufferLength == 0);
  }

  return equalFlag;
}

bool String_subEqualsCString(const String string, const char *s, long index, ulong length)
{
  bool equalFlag;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (s != NULL)
    {
      equalFlag = String_subEqualsBuffer(string,s,strlen(s),index,length);
    }
    else
    {
      equalFlag = (string->length == 0);
    }
  }
  else
  {
    equalFlag = (s == NULL);
  }

  return equalFlag;
}

bool String_subEqualsChar(const String string, char ch, long index)
{
  bool equalFlag;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    equalFlag = ((index < string->length) && (string->data[index] == ch));
  }
  else
  {
    equalFlag = FALSE;
  }

  return equalFlag;
}

long String_find(const String string, ulong index, const String findString)
{
  long findIndex;
  long  i;
  ulong z;

  assert(string != NULL);
  assert(findString != NULL);

  CHECK_VALID(string);
  CHECK_VALID(findString);

  findIndex = -1;

  i = (index != STRING_BEGIN)?index:0;
  while (((i+findString->length) <= string->length) && (findIndex < 0))
  {
    z = 0;
    while ((z < findString->length) && (string->data[i+z] == findString->data[z]))
    {
      z++;
    }
    if (z >=  findString->length) findIndex = i;

    i++;
  }

  return findIndex;
}

long String_findCString(const String string, ulong index, const char *s)
{
  long findIndex;
  long sLength;
  long  i;
  ulong z;

  assert(string != NULL);
  assert(s != NULL);

  CHECK_VALID(string);

  findIndex = -1;

  sLength = strlen(s);
  i = (index != STRING_BEGIN)?index:0;
  while (((i+sLength) <= string->length) && (findIndex < 0))
  {
    z = 0;
    while ((z < sLength) && (string->data[i+z] == s[z]))
    {
      z++;
    }
    if (z >= sLength) findIndex = i;

    i++;
  }

  return findIndex;
}

long String_findChar(const String string, ulong index, char ch)
{
  long i;

  assert(string != NULL);

  CHECK_VALID(string);

  i = (index != STRING_BEGIN)?index:0;
  while ((i < string->length) && (string->data[i] != ch))
  {
    i++;
  }

  return (i < string->length)?i:-1;
}

long String_findLast(const String string, long index, String findString)
{
  long  findIndex;
  long  i;
  ulong z;

  assert(string != NULL);
  assert(findString != NULL);

  CHECK_VALID(string);

  findIndex = -1;

  i = (index != STRING_END)?index:string->length-1;
  while ((i >= 0) && (findIndex < 0))
  {
    z = 0;
    while ((z < findString->length) && (string->data[i+z] == findString->data[z]))
    {
      z++;
    }
    if (z >= findString->length) findIndex = i;

    i--;
  }

  return findIndex;
}

long String_findLastCString(const String string, long index, const char *s)
{
  long findIndex;
  long sLength;
  long  i;
  ulong z;

  assert(string != NULL);
  assert(s != NULL);

  CHECK_VALID(string);

  findIndex = -1;

  sLength = strlen(s);
  i = (index != STRING_END)?index:string->length-1;
  while ((i >= 0) && (findIndex < 0))
  {
    z = 0;
    while ((z < sLength) && (string->data[i+z] == s[z]))
    {
      z++;
    }
    if (z >=  sLength) findIndex = i;

    i--;
  }

  return findIndex;
}

long String_findLastChar(const String string, long index, char ch)
{
  long i;

  assert(string != NULL);

  CHECK_VALID(string);

  i = (index != STRING_END)?index:string->length-1;
  while ((i >= 0) && (string->data[i] != ch))
  {
    i--;
  }

  return (i >= 0)?i:-1;
}

String String_iterate(                      String string,
                      StringIterateFunction stringIterateFunction,
                      void                  *stringIterateUserData
                     )
{
  ulong      z;
  const char *s;
  ulong      n;

  assert(stringIterateFunction != NULL);

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    z = 0;
    while (z < string->length)
    {
      s = stringIterateFunction(stringIterateUserData,string->data[z]);
      if (s != NULL)
      {
        n = strlen(s);
        ensureStringLength(string,string->length+n-1);
        memmove(&string->data[z+n],&string->data[z+1],string->length-(z+1));
        memcpy(&string->data[z],s,n);
        string->data[string->length+n-1] = '\0';
        string->length += n-1;

        z += n;
      }
      else
      {
        z += 1;
      }
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_toLower(String string)
{
  ulong z;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    for (z = 0; z < string->length; z++)
    {
      string->data[z] = tolower(string->data[z]);
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_toUpper(String string)
{
  ulong z;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    for (z = 0; z < string->length; z++)
    {
      string->data[z] = toupper(string->data[z]);
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_trim(String string, const char *chars)
{
  CHECK_VALID(string);

  String_trimRight(string,chars);
  String_trimLeft(string,chars);

  return string;
}

String String_trimRight(String string, const char *chars)
{
  ulong n;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    n = string->length;
    while ((n > 0) && (strchr(chars,string->data[n - 1]) != NULL))
    {
      n--;
    }
    string->data[n] = '\0';
    string->length = n;

    UPDATE_VALID(string);
  }

  return string;
}

String String_trimLeft(String string, const char *chars)
{
  ulong z,n;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    z = 0;
    while ((z < string->length) && (strchr(chars,string->data[z]) != NULL))
    {
      z++;
    }
    if (z > 0)
    {
      n = string->length - z;
      memmove(&string->data[0],&string->data[z],n);
      string->data[n] = '\0';
      string->length = n;
    }

    UPDATE_VALID(string);
  }

  return string;
}

String String_escape(String string, const char *chars, char escapeChar)
{
  String s;
  ulong  z;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    s = allocTmpString();
    for (z = 0; z < string->length; z++)
    {
      if ((string->data[z] == escapeChar) || ((chars != NULL) && (strchr(chars,string->data[z]) != NULL)))
      {
        String_appendChar(s,escapeChar);
      }
      String_appendChar(s,string->data[z]);
    }
    assignTmpString(string,s);
  }

  return string;
}

String String_unescape(String string, char escapeChar)
{
  String s;
  ulong  z;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    s = allocTmpString();
    z = 0;
    while (z < string->length)
    {
      if (string->data[z] == escapeChar)
      {
        z++;
        if (z < string->length)
        {
          String_appendChar(s,string->data[z]);
          z++;
        }
      }
      else
      {
        String_appendChar(s,string->data[z]);
      }
      z++;
    }
    assignTmpString(string,s);
  }

  return string;
}

String String_quote(String string, char quoteChar)
{
  String s;
  ulong  z;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    s = allocTmpString();
    String_appendChar(s,quoteChar);
    for (z = 0; z < string->length; z++)
    {
      if (string->data[z] == quoteChar)
      {
        String_appendChar(s,'\\');
      }
      String_appendChar(s,string->data[z]);
    }
    String_appendChar(s,quoteChar);
    assignTmpString(string,s);
  }

  return string;
}

String String_unquote(String string, const char *quoteChars)
{
  const char *t0,*t1;
  char       quoteChar;
  String     s;
  ulong      z;

  CHECK_VALID(string);

  if (string != NULL)
  {
    assert(string->data != NULL);

    if (string->length > 0)
    {
      t0 = strchr(quoteChars,string->data[0]);
      t1 = strchr(quoteChars,string->data[string->length-1]);
      if ((t0 != NULL) && (t1 != NULL) && ((*t0) == (*t1)))
      {
        quoteChar = (*t0);

        s = allocTmpString();
        z = 1;
        while (z < string->length-1)
        {
          if (string->data[z] == quoteChar)
          {
            z++;
            if (z < string->length-1)
            {
              String_appendChar(s,string->data[z]);
              z++;
            }
          }
          else
          {
            String_appendChar(s,string->data[z]);
          }
          z++;
        }
        assignTmpString(string,s);
      }
    }
  }

  return string;
}

String String_padRight(String string, ulong length, char ch)
{
  ulong n;

  assert(string != NULL);

  CHECK_VALID(string);

  if (string != NULL)
  {
    if (string->length < length)
    {
      n = length-string->length;
      ensureStringLength(string,length);
      memset(&string->data[string->length],ch,n);
      string->data[length] = '\0';
      string->length = length;

      UPDATE_VALID(string);
    }
  }

  return string;
}

String String_padLeft(String string, ulong length, char ch)
{
  ulong n;

  assert(string != NULL);

  CHECK_VALID(string);

  if (string != NULL)
  {
    if (string->length < length)
    {
      n = length-string->length;
      ensureStringLength(string,length);
      memmove(&string->data[n],&string->data[0],string->length);
      memset(&string->data[0],ch,n);
      string->data[length] = '\0';
      string->length = length;

      UPDATE_VALID(string);
    }
  }

  return string;
}

String String_fillChar(String string, ulong length, char ch)
{
  CHECK_VALID(string);

  if (string != NULL)
  {
    ensureStringLength(string,length);
    memset(&string->data[0],ch,length);
    string->data[length] = '\0';
    string->length = length;

    UPDATE_VALID(string);
  }

  return string;
}

void String_initTokenizer(StringTokenizer *stringTokenizer,
                          const String    string,
                          ulong           index,
                          const char      *separatorChars,
                          const char      *stringQuotes,
                          bool            skipEmptyTokens
                         )
{
  assert(stringTokenizer != NULL);
  assert(string != NULL);

  CHECK_VALID(string);

  stringTokenizer->data            = string->data;
  stringTokenizer->length          = string->length;
  stringTokenizer->index           = index;
  stringTokenizer->separatorChars  = separatorChars;
  stringTokenizer->stringQuotes    = stringQuotes;
  stringTokenizer->skipEmptyTokens = skipEmptyTokens;
  #ifdef NDEBUG
    stringTokenizer->token         = String_new();
  #else /* not DEBUG */
    stringTokenizer->token         = __String_new(__FILE__,__LINE__);
  #endif /* NDEBUG */
}

void String_initTokenizerCString(StringTokenizer *stringTokenizer,
                                 const char      *s,
                                 const char      *separatorChars,
                                 const char      *stringQuotes,
                                 bool            skipEmptyTokens
                                )
{
  assert(stringTokenizer != NULL);
  assert(s != NULL);

  stringTokenizer->data            = s;
  stringTokenizer->length          = strlen(s);
  stringTokenizer->index           = 0;
  stringTokenizer->separatorChars  = separatorChars;
  stringTokenizer->stringQuotes    = stringQuotes;
  stringTokenizer->skipEmptyTokens = skipEmptyTokens;
  #ifdef NDEBUG
    stringTokenizer->token         = String_new();
  #else /* not DEBUG */
    stringTokenizer->token         = __String_new(__FILE__,__LINE__);
  #endif /* NDEBUG */
}

void String_doneTokenizer(StringTokenizer *stringTokenizer)
{
  assert(stringTokenizer != NULL);
  assert(stringTokenizer->data != NULL);

  String_delete(stringTokenizer->token);
}

bool String_getNextToken(StringTokenizer *stringTokenizer, String *const token, long *tokenIndex)
{
  const char *s;

  assert(stringTokenizer != NULL);

  /* check index */
  if (stringTokenizer->index >= stringTokenizer->length)
  {
    return FALSE;
  }

  if (stringTokenizer->skipEmptyTokens)
  {
    /* skip separator chars */
    while (   (stringTokenizer->index < stringTokenizer->length)
           && (strchr(stringTokenizer->separatorChars,stringTokenizer->data[stringTokenizer->index]) != NULL)
          )
    {
      stringTokenizer->index++;
    }
    if (stringTokenizer->index >= stringTokenizer->length) return FALSE;
  }

  /* get token */
  if (tokenIndex != NULL) (*tokenIndex) = stringTokenizer->index;
  String_clear(stringTokenizer->token);
  if (stringTokenizer->stringQuotes != NULL)
  {
    while (   (stringTokenizer->index < stringTokenizer->length)
           && (strchr(stringTokenizer->separatorChars,stringTokenizer->data[stringTokenizer->index]) == NULL)
          )
    {
      s = strchr(stringTokenizer->stringQuotes,stringTokenizer->data[stringTokenizer->index]);
      if (s != NULL)
      {
        stringTokenizer->index++;
        while (   (stringTokenizer->index < stringTokenizer->length)
               && (stringTokenizer->data[stringTokenizer->index] != (*s))
              )
        {
          String_appendChar(stringTokenizer->token,stringTokenizer->data[stringTokenizer->index]);
          stringTokenizer->index++;
        }
        if (stringTokenizer->index < stringTokenizer->length) stringTokenizer->index++;
      }
      else
      {
        String_appendChar(stringTokenizer->token,stringTokenizer->data[stringTokenizer->index]);
        stringTokenizer->index++;
      }
    }
  }
  else
  {
    while (   (stringTokenizer->index < stringTokenizer->length)
           && (strchr(stringTokenizer->separatorChars,stringTokenizer->data[stringTokenizer->index]) == NULL)
          )
    {
      String_appendChar(stringTokenizer->token,stringTokenizer->data[stringTokenizer->index]);
      stringTokenizer->index++;
    }
  }
  if (token != NULL) (*token) = stringTokenizer->token;

  /* skip token separator */
  if (   (stringTokenizer->index < stringTokenizer->length)
      && (strchr(stringTokenizer->separatorChars,stringTokenizer->data[stringTokenizer->index]) != NULL)
     )
  {
    stringTokenizer->index++;
  }

  return TRUE;
}

bool String_scan(const String string, ulong index, const char *format, ...)
{
  va_list arguments;
  bool    result;

  assert(string != NULL);
  assert((index == STRING_BEGIN) || (index == STRING_END) || (index < string->length));
  assert(format != NULL);

  CHECK_VALID(string);

  va_start(arguments,format);
  result = parseString(string,index,format,arguments,NULL,NULL);
  va_end(arguments);

  return result;
}

bool String_parse(const String string, ulong index, const char *format, ulong *nextIndex, ...)
{
  va_list arguments;
  bool    result;

  assert(string != NULL);
  assert((index == STRING_BEGIN) || (index == STRING_END) || (index < string->length));
  assert(format != NULL);

  CHECK_VALID(string);

  va_start(arguments,nextIndex);
  result = parseString(string,index,format,arguments,STRING_QUOTES,nextIndex);
  va_end(arguments);

  return result;
}

int String_toInteger(const String convertString, ulong index, long *nextIndex, const StringUnit stringUnits[], uint stringUnitCount)
{
  int  n;
  char *nextData;

  assert(convertString != NULL);

  CHECK_VALID(convertString);

  if (index < convertString->length)
  {
    n = strtol(&convertString->data[index],&nextData,0);
    if ((ulong)(nextData-convertString->data) < convertString->length)
    {
      n = n*getUnitFactor(stringUnits,stringUnitCount,convertString->data,nextData,nextIndex);
    }
    else
    {
      if (nextIndex != NULL) (*nextIndex) = STRING_END;
    }
  }
  else
  {
    n = 0;
    if (nextIndex != NULL) (*nextIndex) = index;
  }

  return n;
}

int64 String_toInteger64(const String convertString, ulong index, long *nextIndex, const StringUnit stringUnits[], uint stringUnitCount)
{
  int64 n;
  char  *nextData;

  assert(convertString != NULL);

  if (index < convertString->length)
  {
    n = strtoll(&convertString->data[index],&nextData,0);
    if ((ulong)(nextData-convertString->data) < convertString->length)
    {
      n = n*(int64)getUnitFactor(stringUnits,stringUnitCount,convertString->data,nextData,nextIndex);
    }
    else
    {
      if (nextIndex != NULL) (*nextIndex) = STRING_END;
    }
  }
  else
  {
    n = 0LL;
    if (nextIndex != NULL) (*nextIndex) = index;
  }

  return n;
}

double String_toDouble(const String convertString, ulong index, long *nextIndex, const StringUnit stringUnits[], uint stringUnitCount)
{
  double n;
  char   *nextData;

  assert(convertString != NULL);

  CHECK_VALID(convertString);

  if (index < convertString->length)
  {
    n = strtod(&convertString->data[index],&nextData);
    if ((ulong)(nextData-convertString->data) < convertString->length)
    {
      n = n*(double)getUnitFactor(stringUnits,stringUnitCount,convertString->data,nextData,nextIndex);
    }
    else
    {
      if (nextIndex != NULL) (*nextIndex) = STRING_END;
    }
  }
  else
  {
    n = 0.0;
    if (nextIndex != NULL) (*nextIndex) = index;
  }

  return n;
}

bool String_toBoolean(const String convertString, ulong index, long *nextIndex, const char *trueStrings[], uint trueStringCount, const char *falseStrings[], uint falseStringCount)
{
  bool       n;
  bool       foundFlag;
  const char **strings;
  uint       stringCount;
  uint       z;

  assert(convertString != NULL);

  CHECK_VALID(convertString);

  n = FALSE;

  if (index < convertString->length)
  {
    foundFlag = FALSE;
    if (!foundFlag)
    {
      if (trueStrings != NULL)
      {
        strings     = trueStrings;
        stringCount = trueStringCount;
      }
      else
      {
        strings     = DEFAULT_TRUE_STRINGS;
        stringCount = SIZE_OF_ARRAY(DEFAULT_TRUE_STRINGS);
      }
      z = 0;
      while (!foundFlag && (z < stringCount))
      {
        if (strcmp(&convertString->data[index],strings[z]) == 0)
        {
          n = TRUE;
          foundFlag = TRUE;
        }
        z++;
      }
    }
    if (!foundFlag)
    {
      if (falseStrings != NULL)
      {
        strings     = falseStrings;
        stringCount = falseStringCount;
      }
      else
      {
        strings     = DEFAULT_FALSE_STRINGS;
        stringCount = SIZE_OF_ARRAY(DEFAULT_FALSE_STRINGS);
      }
      z = 0;
      while (!foundFlag && (z < stringCount))
      {
        if (strcmp(&convertString->data[index],strings[z]) == 0)
        {
          n = FALSE;
          foundFlag = TRUE;
        }
        z++;
      }
    }
    if (!foundFlag)
    {
      if (nextIndex != NULL) (*nextIndex) = index;
    }
  }
  else
  {
    n = 0;
    if (nextIndex != NULL) (*nextIndex) = index;
  }
  
  return n;
}

String String_toString(String string, const String convertString, ulong index, long *nextIndex, const char *stringQuotes)
{
  char *stringQuote;

  if (index < convertString->length)
  {
    while ((index < convertString->length) && !isspace(convertString->data[index]))
    {
      stringQuote = (stringQuotes != NULL)?strchr(stringQuotes,convertString->data[index]):NULL;
      if (stringQuote != NULL)
      {
        do
        {
          /* skip string-char */
          index++;
          /* get string */
          while ((index < convertString->length) && (convertString->data[index] != (*stringQuote)))
          {
            if (convertString->data[index] == '\\')
            {
              if (index < convertString->length)
              {
                String_appendChar(string,convertString->data[index]);
                index++;
              }
            }
            else
            {
              String_appendChar(string,convertString->data[index]);
              index++;
            }
          }
          /* skip string-char */
          index++;
          /* next string char */
          stringQuote = ((stringQuotes != NULL) && (index < convertString->length))?strchr(stringQuotes,convertString->data[index]):NULL;
        }
        while (stringQuote != NULL);
      }
      else
      {
        String_appendChar(string,convertString->data[index]);
        index++;
      }
    }
    if (nextIndex != NULL) (*nextIndex) = index;
  }
  else
  {
    if (nextIndex != NULL) (*nextIndex) = index;
  }

  return string;
}

char* String_toCString(const String string)
{
  char *cString;

  assert(string != NULL);

  CHECK_VALID(string);

  cString = (char*)malloc(string->length+1);
  if (cString == NULL)
  {
    #ifdef HALT_ON_INSUFFICIENT_MEMORY
      HALT_INSUFFICIENT_MEMORY(NULL);
    #else /* not HALT_ON_INSUFFICIENT_MEMORY */
      return NULL;
    #endif /* HALT_ON_INSUFFICIENT_MEMORY */
  }

  memcpy(cString,string->data,string->length);
  cString[string->length] = '\0';

  return cString;
}

#ifndef NDEBUG
void String_debug(void)
{
  DebugStringNode *debugStringNode;

  pthread_once(&debugStringInitFlag,debugStringInit);

  pthread_mutex_lock(&debugStringLock);
  for (debugStringNode = debugAllocStringList.head; debugStringNode != NULL; debugStringNode = debugStringNode->next)
  {
    fprintf(stderr,"DEBUG WARNING: string %p '%s' allocated at %s, line %ld\n",
            debugStringNode->string,
            debugStringNode->string->data,
            debugStringNode->fileName,
            debugStringNode->lineNb
           );
  }
  pthread_mutex_unlock(&debugStringLock);
}

void String_debugDone(void)
{
  pthread_once(&debugStringInitFlag,debugStringInit);

  pthread_mutex_lock(&debugStringLock);
  List_done(&debugFreeStringList,NULL,NULL);
  List_done(&debugFreeStringList,NULL,NULL);
  pthread_mutex_unlock(&debugStringLock);
}

#endif /* not NDEBUG */

#ifdef __cplusplus
  }
#endif

/* end of file */
