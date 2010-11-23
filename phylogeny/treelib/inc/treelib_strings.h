/***********************************************************************\
*
* $Source: /home/torsten/cvs/bar/bar/strings.h,v $
* $Revision: 1.4 $
* $Author: torsten $
* Contents: dynamic string functions
* Systems: all
*
\***********************************************************************/

#ifndef __STRINGS__
#define __STRINGS__

/****************************** Includes *******************************/
#include <stdlib.h>
#include <stdarg.h>

#include "global.h"

/****************** Conditional compilation switches *******************/

/***************************** Constants *******************************/

#define STRING_BEGIN 0
#define STRING_END   -1

#define STRING_WHITE_SPACES " \t\f\v\n\r"
#define STRING_QUOTE        "'"
#define STRING_QUOTES       "\"'"

/***************************** Datatypes *******************************/

/***************************** Variables *******************************/

/* string */
typedef struct __String* String;

/* internal tokenizer data */
typedef struct
{
  const char *data;                     // string data
  ulong      length;                    // string length
  long       index;                     // index in string
  const char *separatorChars;           // token separator characters
  const char *stringQuotes;             // string quote characters
  bool       skipEmptyTokens;           // TRUE for skipping empty tokens
  String     token;                     // next token
} StringTokenizer;

/* comparison, iteration functions */
typedef int(*StringCompareFunction)(void *userData, char ch1, char ch2);
typedef const char*(*StringIterateFunction)(void *userData, char ch);

typedef struct
{
  const char *name;
  uint64     factor;
} StringUnit;

/****************************** Macros *********************************/

#ifndef NDEBUG
  #define String_new()                          __String_new(__FILE__,__LINE__)
  #define String_newCString(s)                  __String_newCString(__FILE__,__LINE__,s)
  #define String_newChar(ch)                    __String_newChar(__FILE__,__LINE__,ch)
  #define String_newBuffer(buffer,bufferLength) __String_newBuffer(__FILE__,__LINE__,buffer,bufferLength)
  #define String_duplicate(fromString)          __String_duplicate(__FILE__,__LINE__,fromString)
  #define String_copy(string,fromString)        __String_copy(__FILE__,__LINE__,string,fromString)
  #define String_delete(string)                 __String_delete(__FILE__,__LINE__,string)
#endif /* not NDEBUG */

/***************************** Forwards ********************************/

/***************************** Functions *******************************/

#ifdef __cplusplus
  extern "C" {
#endif

/***********************************************************************\
* Name   : String_new, String_newCString, String_newChar,
*          String_newBuffer
* Purpose: create new string
* Input  : s            - C-string
*          ch           - character
*          buffer       - buffer
*          bufferLength - length of buffer
* Output : -
* Return : string or NULL
* Notes  : -
\***********************************************************************/

#ifdef NDEBUG
String String_new(void);
String String_newCString(const char *s);
String String_newChar(char ch);
String String_newBuffer(const void *buffer, ulong bufferLength);
#else /* not NDEBUG */
String __String_new(const char *fileName, ulong lineNb);
String __String_newCString(const char *fileName, ulong lineNb, const char *s);
String __String_newChar(const char *fileName, ulong lineNb, char ch);
String __String_newBuffer(const char *fileName, ulong lineNb, const void *buffer, ulong bufferLength);
#endif /* NDEBUG */

/***********************************************************************\
* Name   : String_duplicate
* Purpose: duplicate string
* Input  : fromString - string to duplicate
* Output : -
* Return : new string (copy)
* Notes  : -
\***********************************************************************/

#ifdef NDEBUG
String String_duplicate(const String fromString);
#else /* not NDEBUG */
String __String_duplicate(const char *fileName, ulong lineNb, const String fromString);
#endif /* NDEBUG */

/***********************************************************************\
* Name   : String_copy
* Purpose: copy string
* Input  : string     - reference to string variable (if referenced
*                       string is NULL, a new string is allocated)
*          fromString - string to copy
* Output : -
* Return : new string (copy)
* Notes  : -
\***********************************************************************/

#ifdef NDEBUG
String String_copy(String *string, const String fromString);
#else /* not NDEBUG */
String __String_copy(const char *fileName, ulong lineNb, String *string, const String fromString);
#endif /* NDEBUG */

/***********************************************************************\
* Name   : String_delete
* Purpose: delete string
* Input  : string - string to delete
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

#ifdef NDEBUG
void String_delete(String string);
#else /* not NDEBUG */
void __String_delete(const char *fileName, ulong lineNb, String string);
#endif /* NDEBUG */

/***********************************************************************\
* Name   : String_clear
* Purpose: clear string
* Input  : string - string to clear
* Output : -
* Return : cleared string (empty)
* Notes  : -
\***********************************************************************/

String String_clear(String string);

/***********************************************************************\
* Name   : String_erase
* Purpose: erase string content
* Input  : string - string to erase
* Output : -
* Return : erased string (empty)
* Notes  : -
\***********************************************************************/

String String_erase(String string);

/***********************************************************************\
* Name   : String_set, String_setCString, String_setChar,
*          Stirng_setBuffer
* Purpose: set string (copy string)
* Input  : string       - string to set
*          sourceString - source string
*          s            - C-string
*          ch           - character
*          buffer       - buffer
*          bufferLength - length of buffer
* Output : -
* Return : string
* Notes  : -
\***********************************************************************/

String String_set(String string, const String sourceString);
String String_setCString(String string, const char *s);
String String_setChar(String string, char ch);
String String_setBuffer(String string, const void *buffer, ulong bufferLength);

/***********************************************************************\
* Name   : String_sub
* Purpose: get sub-string from string
* Input  : string/buffer - string/buffer to set
*          fromString    - string to get sub-string from
*          index         - start index (0..n-1)
*          length        - length of sub-string (0..n) or STRING_END
*                          (String_sub only!)
* Output : -
* Return : new sub-string/buffer
* Notes  : -
\***********************************************************************/

String String_sub(String string, const String fromString, ulong index, long length);
char *String_subCString(char *s, const String fromString, ulong index, long length);
char *String_subBuffer(char *buffer, const String fromString, ulong index, long length);

/***********************************************************************\
* Name   : String_append, String_appendSub, String_appendBuffer,
*          String_appendCString, String_appendChar
* Purpose: append to string
* Input  : string                - string
*          appendString/s/buffer - string to append
*          ch                    - character to append
*          bufferLength          - buffer length
* Output : -
* Return : string
* Notes  : -
\***********************************************************************/

String String_append(String string, String appendString);
String String_appendSub(String string, const String fromString, ulong fromIndex, long fromLength);
String String_appendBuffer(String string, const char *buffer, ulong bufferLength);
String String_appendCString(String string, const char *s);
String String_appendChar(String string, char ch);

/***********************************************************************\
* Name   : String_insert, String_insertSub, String_insertCString,
*          String_insertChar
* Purpose: insert into string
* Input  : string                - string
*          index                 - index where to insert
*          insertString/s/buffer - string to insert
*          ch                    - character to insert
*          bufferLength          - buffer length
* Output : -
* Return : string
* Notes  : -
\***********************************************************************/

String String_insert(String string, ulong index, const String insertString);
String String_insertSub(String string, ulong index, const String fromString, ulong fromIndex, long fromLength);
String String_insertBuffer(String string, ulong index, const char *buffer, ulong bufferLength);
String String_insertCString(String string, ulong index, const char *s);
String String_insertChar(String string, ulong index, char ch);

/***********************************************************************\
* Name   : String_remove
* Purpose: remove part of string
* Input  : string - string
*          index  - index of first character to remove or STRING_END to
*                   remove characters from end
*          length - length to remove
* Output : -
* Return : string
* Notes  : -
\***********************************************************************/

String String_remove(String string, ulong index, ulong length);

/***********************************************************************\
* Name   : String_replace, String_replaceCString, String_replaceChar
* Purpose: replace part of string with other string
* Input  : string                - string                
*          index                 - index where to insert 
*          length                - length to replace     
*          insertString/s/buffer - string to insert
*          ch                    - character to insert
*          bufferLength          - buffer length
* Output : -
* Return : string
* Notes  : -
\***********************************************************************/

String String_replace(String string, ulong index, ulong length, const String insertString);
String String_replaceBuffer(String string, ulong index, ulong length, const char *buffer, ulong bufferLength);
String String_replaceCString(String string, ulong index, ulong length, const char *s);
String String_replaceChar(String string, ulong index, ulong length, char ch);

/***********************************************************************\
* Name   : String_length
* Purpose: get string length
* Input  : string - string
* Output : -
* Return : length of string (0..n)
* Notes  : -
\***********************************************************************/

ulong String_length(const String string);

/***********************************************************************\
* Name   : String_empty
* Purpose: check if string is empty
* Input  : string - string
* Output : -
* Return : TRUE iff string is empty, FALSE otherwise
* Notes  : -
\***********************************************************************/

bool String_empty(const String string);

/***********************************************************************\
* Name   : String_index
* Purpose: get char at index
* Input  : string - string
*          index  - index (0..n-1) or STRING_END to get last
*                   character
* Output : -
* Return : character at position "index"
* Notes  : -
\***********************************************************************/

char String_index(const String string, ulong index);

/***********************************************************************\
* Name   : String_cString
* Purpose: get C-string from string
* Input  : string - string
* Output : -
* Return : C-string
* Notes  : -
\***********************************************************************/

const char *String_cString(const String string);

/***********************************************************************\
* Name   : String_compare
* Purpose: compare two strings
* Input  : string1,string2       - strings to compare
*          stringCompareFunction - string compare function (can be NULL)
*          stringCompareUserData - user data for compare function
* Output : -
* Return : -1 if string1 <  string2
*           0 if string1 == string2
*          +1 if string1 >  string2
* Notes  : -
\***********************************************************************/

int String_compare(const String          string1,
                   const String          string2,
                   StringCompareFunction stringCompareFunction,
                   void                  *stringCompareUserData
                  );

/***********************************************************************\
* Name   : String_equals, String_equalsBuffer, String_equalsCString
*          String_equalsChar
* Purpose: check if strings are equal
* Input  : string1,string2            - strings to compare
*          string/buffer/bufferLength - string/buffer to compare
*          string/s                   - string/C-string to compare
*          string,ch                  - strings/character to compare
* Output : -
* Return : TRUE if strings equal, FALSE otherwise
* Notes  : -
\***********************************************************************/

bool String_equals(const String string1, const String string2);
bool String_equalsBuffer(const String string, const char *buffer, ulong bufferLength);
bool String_equalsCString(const String string, const char *s);
bool String_equalsChar(const String string, char ch);

/***********************************************************************\
* Name   : String_subEquals, String_subEqualsBuffer,
*          String_subEqualsCString, String_subEqualsChar
* Purpose: check if string2 is equal to string1 at some position
* Input  : string1,string2            - strings to compare
*          string/buffer/bufferLength - string/buffer to compare
*          string/s                   - string/C-string to compare
*          string,ch                  - strings/character to compare
*          index                      - position in string1
*          length                     - length to compare
* Output : -
* Return : TRUE if strings equal, FALSE otherwise
* Notes  : -
\***********************************************************************/

bool String_subEquals(const String string1, const String string2, long index, ulong length);
bool String_subEqualsBuffer(const String string, const char *buffer, ulong bufferLength, long index, ulong length);
bool String_subEqualsCString(const String string, const char *s, long index, ulong length);
bool String_subEqualsChar(const String string, char ch, long index);

/***********************************************************************\
* Name   : String_find, String_findCString, String_findChar
* Purpose: find string in string
* Input  : string - string
*          index - index to start search (0..n-1)
*          findString - string to find/s - C-string to find/
*          ch - character to find
* Output : -
* Return : index in string or -1 if not found
* Notes  : -
\***********************************************************************/

long String_find(const String string, ulong index, const String findString);
long String_findCString(const String string, ulong index, const char *s);
long String_findChar(const String string, ulong index, char ch);
long String_findLast(const String string, long index, const String findString);
long String_findLastCString(const String string, long index, const char *s);
long String_findLastChar(const String string, long index, char ch);

/***********************************************************************\
* Name   : String_iterate
* Purpose: iterate over string
* Input  : string                - string
*          stringIterateFunction - iterator function
*          stringIterateUserData - user data for iterator function
* Output : -
* Return : string
* Notes  : Note: returned string of iterate function replaces character
*          in string
\***********************************************************************/

String String_iterate(String                string,
                      StringIterateFunction stringIterateFunction,
                      void                  *stringIterateUserData
                     );

/***********************************************************************\
* Name   : String_toLower, String_toUpper
* Purpose: convert string to lower/upper case
* Input  : string - string
* Output : -
* Return : converted string
* Notes  : -
\***********************************************************************/

String String_toLower(String string);
String String_toUpper(String string);

/***********************************************************************\
* Name   : String_trim, String_trimRight, String_trimLeft
* Purpose: trim string right/left
* Input  : string - string
*          chars  - chars to trim
* Output : -
* Return : string
* Notes  : -
\***********************************************************************/

String String_trim(String string, const char *chars);
String String_trimRight(String string, const char *chars);
String String_trimLeft(String string, const char *chars);

/***********************************************************************\
* Name   : String_escape
* Purpose: escape string
* Input  : string     - string variable
*          chars      - characters to escape
*          escapeChar - escape char
* Output : -
* Return : escaped string
* Notes  : -
\***********************************************************************/

String String_escape(String string, const char *chars, char escapeChar);

/***********************************************************************\
* Name   : String_unescape
* Purpose: unescape string
* Input  : string     - string variable
*          escapeChar - escape char
* Output : -
* Return : unescaped string
* Notes  : -
\***********************************************************************/

String String_unescape(String string, char escapeChar);

/***********************************************************************\
* Name   : String_quote, String_unquote
* Purpose: quote/unquote string
* Input  : string     - string
*          quoteChar  - quote character to add
*          quoteChars - quote characters to remove
* Output : -
* Return : quoted/unquoted string
* Notes  : add quote character and escape enclosed quote characters if
*          needed
\***********************************************************************/

String String_quote(String string, char quoteChar);
String String_unquote(String string, const char *quoteChars);

/***********************************************************************\
* Name   : String_rightPad, String_rightPad
* Purpose: pad string right/left
* Input  : string - string
*          length - length to pad
*          ch     - padding char
* Output : -
* Return : string
* Notes  : -
\***********************************************************************/

String String_padRight(String string, ulong length, char ch);
String String_padLeft(String string, ulong length, char ch);

/***********************************************************************\
* Name   : String_fill
* Purpose: fill string with character
* Input  : string - string
*          length - length to fill
*          ch     - fill char
* Output : -
* Return : string
* Notes  : -
\***********************************************************************/

//String String_fillString(String string, ulong length, char ch);
String String_fillChar(String string, ulong length, char ch);

/***********************************************************************\
* Name   : String_format, String String_vformat
* Purpose: format string and append
* Input  : string - string
*          format - printf-like format string
*          ...    - arguments
* Output : -
* Return : format string
* Notes  : additional format characters
*           %S   String
*           %cS  String with quoting char c
\***********************************************************************/

String String_format(String string, const char *format, ...);
String String_vformat(String string, const char *format, va_list arguments);

/***********************************************************************\
* Name   : String_initTokenizer, String_initTokenizerCString,
*          String_doneTokenizer
* Purpose: initialise/deinitialise string tokenizer
* Input  : stringTokenizer - string tokenizer
*          string          - string
*          separatorChars  - token seperator characters, e. g. " "
*          stringQuotes    - token string quote characters, e. g. ",'
*          skipEmptyTokens - TRUE to skip empty tokens, FALSE to get
*                            also empty tokens
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void String_initTokenizer(StringTokenizer *stringTokenizer,
                          const String    string,
                          ulong           index,
                          const char      *separatorChars,
                          const char      *stringQuotes,
                          bool            skipEmptyTokens
                         );
void String_initTokenizerCString(StringTokenizer *stringTokenizer,
                                 const char      *s,
                                 const char      *separatorChars,
                                 const char      *stringQuotes,
                                 bool            skipEmptyTokens
                                );
void String_doneTokenizer(StringTokenizer *stringTokenizer);

/***********************************************************************\
* Name   : String_getNextToken
* Purpose: find next token
* Input  : stringTokenizer - string tokenizer
* Output : token      - token (internal reference; do not delete!)
*          tokenIndex - token index (could be NULL)
* Return : TRUE if token found, FALSE otherwise
* Notes  : -
\***********************************************************************/

bool String_getNextToken(StringTokenizer *stringTokenizer,
                         String *const   token,
                         long            *tokenIndex
                        );

/***********************************************************************\
* Name   : String_scan
* Purpose: scan string
* Input  : string - string
*          format - format (like scanf)
*          ...    - optional variables
* Output : -
* Return : TRUE is scanned, FALSE on error
* Notes  : for C-strings the max. length have to be specified with %<n>s
\***********************************************************************/

bool String_scan(const String string, ulong index, const char *format, ...);

/***********************************************************************\
* Name   : String_parse
* Purpose: parse string
* Input  : string - string
*          format - format (like scanf)
*          ...    - optional variables
* Output : nextIndex - index of next character in string not parsed (can
*                      be NULL)
* Return : TRUE is parsed, FALSE on error
* Notes  : extended scan-function:
*            - match also specified text
*            - %<n>s will return max. <n-1> characters and always add
*              all \0 at the end of the string
*            - %s and %S are parsed as strings which could be enclosed
*              in "..." or '...'
*            - % s and % S parse rest of string (including spaces)
*            - if a value is NULL, skip value
\***********************************************************************/

bool String_parse(const String string, ulong index, const char *format, ulong *nextIndex, ...);

/***********************************************************************\
* Name   : String_match, String_matchCString
* Purpose: match string pattern
* Input  : string      - string
*          index       - start index in string
*          pattern     - pattern
*          matchString - regular expression match string
*          ...         - optional sub-patterns (strings), last value
*                        have to be NULL!
* Output : -
* Return : TRUE if pattern is matching, FALSE otherwise
* Notes  : -
\***********************************************************************/

bool String_match(const String string, ulong index, const String pattern, String matchString, ...);
bool String_matchCString(const String string, ulong index, const char *pattern, String matchString, ...);

/***********************************************************************\
* Name   : String_toInteger, String_toInteger64, String_toDouble,
*          String_toBoolean
* Purpose: convert string into integer, integer64, double, boolean or
*          string (string without enclosing quotes ' or ")
* Input  : string                   - string variable (for string)
*          convertString            - string to convert
*          index                    - start index in convertString
*          stringUnits              - string units (for integer, double)
*          stringUnitCount          - number of string units (for
*                                     integer, double)
*          trueStrings,falseStrings - string false texts (for boolean)
*          stringQuotes             - string quotes (for string)
* Output : nextIndex - index of next character in string not parsed or
*                      STRING_END if string completely parsed (can be
*                      NULL)
*
* Return : integer/integer64/double/boolean/string value
* Notes  : -
\***********************************************************************/

int String_toInteger(const String convertString, ulong index, long *nextIndex, const StringUnit stringUnits[], uint stringUnitCount);
int64 String_toInteger64(const String string, ulong index, long *nextIndex, const StringUnit stringUnits[], uint stringUnitCount);
double String_toDouble(const String convertString, ulong index, long *nextIndex, const StringUnit stringUnits[], uint stringUnitCount);
bool String_toBoolean(const String convertString, ulong index, long *nextIndex, const char *trueStrings[], uint trueStringCount, const char *falseStrings[], uint falseStringCount);
String String_toString(String string, const String convertString, ulong index, long *nextIndex, const char *stringQuotes);

/***********************************************************************\
* Name   : String_toCString
* Purpose: allocate memory and convert to C-string
* Input  : string - string
* Output : -
* Return : C-string or NULL on insufficient memory
* Notes  : memory have to be deallocated by free()!
\***********************************************************************/

char* String_toCString(const String string);

#ifndef NDEBUG
/***********************************************************************\
* Name   : String_debug
* Purpose: string debug function: output not deallocated strings
* Input  : -
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void String_debug(void);

/***********************************************************************\
* Name   : String_debugDone
* Purpose: done string debug functions
* Input  : -
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void String_debugDone(void);
#endif /* not NDEBUG */

#ifdef __cplusplus
  }
#endif

#endif /* __STRINGS__ */

/* end of file */
