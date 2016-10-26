/***********************************************************************\
*
* File    : global.h
* Author  : Torsten Rupp
* Contents: global definitions
* Systems : Linux
*
************************************************************************/

#ifndef _GLOBAL_H
 #define _GLOBAL_H

// #if (defined DEBUG)
//  #warning DEBUG option set - no LOCAL and no -O2 (optimizer) will be used!
// #endif

/****************************** Includes *******************************/
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <stdbool.h>

/****************** Conditional compilation switches *******************/

/***************************** Constants *******************************/
#define DEBUG_LEVEL 8                          // debug level

/* definition of boolean values */
#ifdef __cplusplus
  #ifndef FALSE
    #define FALSE false
  #endif
  #ifndef TRUE
    #define TRUE true
  #endif
#else
  #ifndef FALSE
    #define FALSE 0
  #endif
  #ifndef TRUE
    #define TRUE  1
  #endif
#endif
#define NO  FALSE
#define YES TRUE
#define OFF FALSE
#define ON  TRUE

#ifndef PI
  #define PI 3.14159265358979323846
#endif

#define MILLI_PER_SECOND 1000

/* minimal and maximal values for some scalar data types */
#define MIN_CHAR           CHAR_MIN
#define MAX_CHAR           CHAR_MAX
#define MIN_SHORTINT       SHRT_MIN
#define MAX_SHORTINT       SHRT_MAX
#define MIN_INT            INT_MIN
#define MAX_INT            INT_MAX
#define MIN_LONG           LONG_MIN
#define MAX_LONG           LONG_MAX

#define MIN_UCHAR          0
#define MAX_UCHAR          UCHAR_MAX
#define MIN_USHORT         0
#define MAX_USHORT         USHRT_MAX
#define MIN_UINT           0
#define MAX_UINT           UINT_MAX
#define MIN_ULONG          0
#define MAX_ULONG          ULONG_MAX

#define MIN_FLOAT          FLT_MIN
#define MAX_FLOAT          FLT_MAX
#define EPSILON_FLOAT      FLT_EPSILON
#define MIN_DOUBLE         DBL_MIN
#define MAX_DOUBLE         DBL_MAX
#define EPSILON_DOUBLE     DBL_EPSILON
#define MIN_LONGDOUBLE     LDBL_MIN
#define MAX_LONGDOUBLE     LDBL_MAX
#define EPSILON_LONGDOUBLE LDBL_EPSILON

/* special constants */
#define NO_WAIT      0
#define WAIT_FOREVER (-1)

/* exit codes */
#define EXITCODE_INTERNAL_ERROR 128

/**************************** Datatypes ********************************/
//#ifndef __cplusplus
//typedef int                bool;
//#endif
typedef unsigned char      uchar;
typedef short int          shortint;
typedef unsigned short int ushortint;

typedef unsigned long int ulong;
typedef unsigned int uint;

typedef enum
{
  CMP_LOWER=-1,
  CMP_EQUAL=0,
  CMP_GREATER=+1
} TCmpResults;

// ??? configure
typedef unsigned char      byte;

typedef char               bool8;
typedef int                bool32;
typedef char               char8;
typedef char               uchar8;
typedef char               int8;
typedef short int          int16;
typedef int                int32;
typedef long long          int64;
typedef unsigned char      uint8;
typedef unsigned short int uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;
typedef void               void32;

/**************************** Variables ********************************/

/****************************** Macros *********************************/
#define GLOBAL extern
#undef LOCAL

#ifndef DEBUG
  #define LOCAL static
#else
  #define LOCAL
#endif

#define UNUSED_VARIABLE(s) (void)s

#define SIZE_OF_ARRAY(a) (sizeof(a)/sizeof(a[0]))

#define ALIGN(n,alignment) (((alignment)>0)?(((n)+(alignment)-1) & ~((alignment)-1)):(n))

/* mathematicl macros */
#ifndef __cplusplus
  #define IS_NAN(d) (!((d)==(d)))                        // check if Not-A-Number (NaN)
  #define IS_INF(d) ((d<-MAX_DOUBLE) || (d>MAX_DOUBLE))  // check if infinit-number
#endif
#ifdef __GNUG__
  #define MIN(x,y) ((x)<?(y))
  #define MAX(x,y) ((x)>?(y))
  #define IN_RANGE(l,x,h) (( ((x)<?(h) )>?(l))
#else
  #define MIN(x,y) (((x)<(y)) ? (x) : (y))
  #define MAX(x,y) (((x)>(y)) ? (x) : (y))
  #define IN_RANGE(l,x,h) ((x)<(l) ? (l) : ((x)>(h) ? (h) : (x)))
#endif
#define SQUARE(x) ((x)*(x))
#define CHECK_RANGE(l,x,u) (( ((l)<=(x)) && ((x)<=(u)) ) || \
                            ( ((u)<=(x)) && ((x)<=(l)) )    \
                           )
#define CHECK_ANGLE_RANGE(l,a,u) (((NormRad(l))<=(NormRad(u)))?CHECK_RANGE(NormRad(l),NormRad(a),NormRad(u)):(CHECK_RANGE(l,a,2*PI) || CHECK_RANGE(0,NormRad(a),NormRad(u))))
#ifndef __cplusplus
 #define IndexMod(l,i,u) ( l+(((i)<0)?( ((u)-(l)+1)-((-(i)%((u)-(l)+1)))%((u)-(l)+1) ):( ((i)>((u)-(l)+1))?( (i)%((u)-(l)+1) ):( (i) ) )) )
#endif

/* used in constant declarations */
#define RAD_TO_DEGREE(n) ((n)>=0?\
                          (((n)<=2*PI)?\
                           (n)*180.0/PI:\
                           ((n)-2*PI)*180.0/PI\
                          ):\
                          (((n)+2*PI)*180.0/PI)\
                         )
#define DEGREE_TO_RAD(n) ((n)>=0?\
                          (((n)<=360)?\
                           (n)*PI/180.0:\
                           ((n)-360)*PI/180.0\
                          ):\
                          (((n)+360)*PI/180.0)\
                         )

#ifndef __cplusplus
  #define RadToDegree(n) RAD_TO_DEGREE(n)
  #define DegreeToRad(n) DEGREE_TO_RAD(n)
#endif
#define NORM_RAD(n)       ((n)<0?((n)+2*PI):(n)>(2*PI)?((n)-2*PI):(n)) /* used in constant declarations */
#define NORM_DEGREE(n)    ((n)<0?((n)+360):(n)>360?((n)-360):(n))      /* used in constant declarations */
#define NORM_RAD90(n)     ((n)>3*PI/2?\
                            ((n)-2*PI):\
                            ((n)>PI/2?\
                             ((n)-PI):\
                             ((n)<-3*PI/2?\
                              ((n)+2*PI):\
                              ((n)<-PI/2?\
                               ((n)+PI):\
                               (n)\
                              )\
                             )\
                            )\
                          )                                            /* used in constant declarations */
#define NORM_RAD180(n)    ((n)>PI?((n)-2*PI):(n)<-PI?((n)+2*PI):(n))   /* used in constant declarations */
#define NORM_DEGREE180(n) ((n)>180?((n)-360):(n)<-180?((n)+360):(n))   /* used in constant declarations */
#define NORM_RAD360(n)    (fmod(n,2*PI))                               /* used in constant declarations */
#define NORM_DEGREE360(n) (fmod(n,360))                                /* used in constant declarations */
/* used in code */
#ifndef __cplusplus
  #define SwapWORD(n) ( ((n & 0xFF00) >> 8) | \
                        ((n & 0x00FF) << 8)   \
                      )
  #define SwapLONG(n) ( ((n & 0xFF000000) >> 24) | \
                        ((n & 0x00FF0000) >>  8) | \
                        ((n & 0x0000FF00) <<  8) | \
                        ((n & 0x000000FF) << 24)   \
                      )
#endif

/* 2 macros necessary, because of "string"-construction */
#define _HALT_STRING1(z) _HALT_STRING2(z) 
#define _HALT_STRING2(z) #z
#undef HALT
#define HALT(errorLevel, format, ...) \
 do \
  { \
   __halt(__FILE__,__LINE__,errorLevel,format, __VA_ARGS__); \
  } \
 while (0)

#define HALT_INSUFFICIENT_MEMORY(...) \
  do \
  { \
     __abort(__FILE__,__LINE__,"FATAL ERROR: ","insufficient memory", __VA_ARGS__); \
  } \
 while (0)

#define HALT_FATAL_ERROR(format, ...) \
  do \
  { \
     __abort(__FILE__, __LINE__,"FATAL ERROR: ",format, __VA_ARGS__); \
  } \
 while (0)

#define HALT_INTERNAL_ERROR(format, ...) \
  do \
  { \
     __abort(__FILE__, __LINE__, "INTERNAL ERROR: ", format, __VA_ARGS__); \
  } \
  while (0)
#define HALT_INTERNAL_ERROR_STILL_NOT_IMPLEMENTED(...) \
  do \
  { \
     HALT_INTERNAL_ERROR("still not implemented", __VA_ARGS__); \
  } \
  while (0)
#define HALT_INTERNAL_ERROR_UNHANDLED_SWITCH_CASE(...) \
  do \
  { \
     HALT_INTERNAL_ERROR("unhandled switch case", __VA_ARGS__); \
  } \
  while (0)
#define HALT_INTERNAL_ERROR_UNREACHABLE() \
  do \
  { \
     HALT_INTERNAL_ERROR("unreachable code"); \
  } \
  while (0)


/* 2 macros necessary, because of "string"-construction */
#define _FAIL_STRING1(z) _FAIL_STRING2(z) 
#define _FAIL_STRING2(z) #z
#define FAIL(errorLevel, format, ...) \
 do \
  { \
   fprintf(stderr, format " - fail in file " __FILE__ ", line " _FAIL_STRING1(__LINE__) "\n" , __VA_ARGS__); \
   exit(errorLevel);\
  } \
 while (0)

#define MEMSET(p,value,size) memset(p,value,size)
#define MEMCLEAR(p,size) memset(p,0,size)

/**************************** Functions ********************************/

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NDEBUG
/***********************************************************************\
* Name       : __dprintf
* Purpose    : debug printf
* Input      : filename - filename
*              lineNb   - line number
*              format   - format string (like printf)
*              ...      - optional arguments
* Output     : -
* Return     : -
* Side-effect: unknown
* Notes      : -
\***********************************************************************/

void __dprintf(const char *filename, unsigned int lineNb, const char *format, ...);
#endif /* NDEBUG */

#ifdef __cplusplus

/***********************************************************************\
* Name       : isNaN
* Purpose    : check if not-a-number (NaN)
* Input      : n - numer
* Output     : -
* Return     : TRUE if n is NaN, FALSE otherwise
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline bool isNaN(double n)
{
  return(n!=n);
}

/***********************************************************************\
* Name       : isInf
* Purpose    : check if number is infinit
* Input      : n - number
* Output     : -
* Return     : TRUE if n is infinit, FALSE otherwise
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline bool isInf(double n)
{
  return((n<-MAX_DOUBLE) || (n>MAX_DOUBLE));
}

#endif

#ifdef __cplusplus

/***********************************************************************\
* Name       : radToDegree
* Purpose    : convert rad to degree
* Input      : n - angle in rad
* Output     : -
* Return     : angle in degree
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline double radToDegree(double n)
{
//???  ASSERT_NaN(n);
//  n=fmod(n,2*PI);
//  if (n<0) n+=2*PI;
  return(n*180/PI);
}

/***********************************************************************\
* Name       : degreeToRad
* Purpose    : convert degree to rad
* Input      : n - angle in degree
* Output     : -
* Return     : angle in rad
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline double degreeToRad(double n)
{
//???  ASSERT_NaN(n);
//  n=fmod(n,360);
//  if (n<0) n+=360;
  return(n*PI/180);
}

#endif

#ifdef __cplusplus

/***********************************************************************\
* Name       : normRad
* Purpose    : normalize angle in rad (0..2PI)
* Input      : n - angle in rad
* Output     : -
* Return     : normalized angle (0..2PI)
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline double normRad(double n)
{
//???  ASSERT(!IsNaN(n));
  n=fmod(n,2*PI);
  if (n<0) n+=2*PI;
  return(n);
}

/***********************************************************************\
* Name       : normDegree
* Purpose    : normalize angle in degree (0..360)
* Input      : n - angle in degree
* Output     : -
* Return     : normalize angle (0..360)
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline double normDegree(double n)
{
//???  ASSERT(!IsNaN(n));
  n=fmod(n,360);
  if (n<0) n+=360;
  return(n);
}

/***********************************************************************\
* Name       : normRad90
* Purpose    : normalize angle in rad (-PI/2..PI/2)
* Input      : n - angle in rad
* Output     : -
* Return     : normalized angle (-PI/2..PI/2)
* Side-Effect: unknown
* Notes      : PI/2..3PI/2   = -PI/2..PI/2
*              3PI/2..2PI    = -PI/2..0
*              -PI/2..-3PI/2 = PI/2..-PI/2
*              -3PI/2..-2PI  = PI/2..0
\***********************************************************************/

inline double normRad90(double n)
{
//???  ASSERT(!IsNaN(n));
  if (n> 3*PI/2) n-=2*PI;
  if (n<-3*PI/2) n+=2*PI;
  if (n>   PI/2) n-=PI;
  if (n<  -PI/2) n+=PI;
  return(n);
}

/***********************************************************************\
* Name       : normRad180
* Purpose    : normalize angle in rad (-PI..PI)
* Input      : n - angle in rad
* Output     : -
* Return     : normalized angle (-PI..PI)
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline double normRad180(double n)
{
//???  ASSERT(!IsNaN(n));
  return(n>PI?(n-2*PI):n<-PI?(n+2*PI):n);
}

/***********************************************************************\
* Name       : normDegree90
* Purpose    : normalize angle in rad (-90..90)
* Input      : n - angle in degree
* Output     : -
* Return     : normalized angle (-90..90)
* Side-Effect: unknown
* Notes      : 90..270    = -90..
*              270..360   = -90..0
*              -90..-270  = 90..-90
*              -270..-360 = 90..0
\***********************************************************************/

inline double normDegree90(double n)
{
//???  ASSERT(!IsNaN(n));
  if (n> 270) n-=360;
  if (n<-270) n+=360;
  if (n>  90) n-=180;
  if (n< -90) n+=180;
  return(n);
}

/***********************************************************************\
* Name       : normDegree180
* Purpose    : normalize angle in degree (-180..180)
* Input      : n - angle in degree
* Output     : -
* Return     : normalize angle (-180..180)
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline double normDegree180(double n)
{
//???  ASSERT(!IsNaN(n));
  return(n>180?(n-360):n<-180?(n+360):n);
}

/***********************************************************************\
* Name       : normRad360
* Purpose    : normalize angle in rad (-2PI...2PI)
* Input      : n - angle in rad
* Output     : -
* Return     : normalized angle (-2PI..2PI)
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline double normRad360(double n)
{
//???  ASSERT(!IsNaN(n));
  return(fmod(n,2*PI));
}

/***********************************************************************\
* Name       : normDegree360
* Purpose    : normalize angle in degree (-360..360)
* Input      : n - angle in degree
* Output     : -
* Return     : normalize angle (-360..360)
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline double normDegree360(double n)
{
//???  ASSERT(!IsNaN(n));
  return(fmod(n,360));
}

#endif

/*---------------------------------------------------------------------*/

#ifdef __cplusplus

/***********************************************************************\
* Name       : swapWORD
* Purpose    : swap low/high byte of word (2 bytes)
* Input      : n - word (a:b)
* Output     : -
* Return     : swapped word (b:a)
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline ushort swapWORD(ushort n)
{
  return( ((n & 0xFF00) >> 8) |
          ((n & 0x00FF) << 8)
        );
}

/***********************************************************************\
* Name       : swapLONG
* Purpose    : swap bytes of long (4 bytes)
* Input      : n - long (a:b:c:d)
* Output     : -
* Return     : swapped long (d:c:b:a)
* Side-Effect: unknown
* Notes      : -
\***********************************************************************/

inline ulong swapLONG(ulong n)
{
  return( ((n & 0xFF000000) >> 24) | 
          ((n & 0x00FF0000) >>  8) |
          ((n & 0x0000FF00) <<  8) |
          ((n & 0x000000FF) << 24)
        );
}

#endif

/***********************************************************************\
* Name   : __halt
* Purpose: halt program
* Input  : filename - filename
*          lineNb   - line number
*          exitcode - exitcode
*          format   - format string (like printf)
*          ...      - optional arguments
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void __halt(const char   *filename,
            unsigned int lineNb,
            int          exitcode,
            const char   *format,
            ...
           );

/***********************************************************************\
* Name   : __abort
* Purpose: abort program
* Input  : filename - filename
*          lineNb   - line number
*          prefix   - prefix text
*          format   - format string (like printf)
*          ...      - optional arguments
* Output : -
* Return : -
* Notes  : -
\***********************************************************************/

void __abort(const char   *filename,
             unsigned int lineNb,
             const char   *prefix,
             const char   *format,
             ...
            );

#ifdef __cplusplus
}
#endif

#endif /* _GLOBAL_H */

/* end of file */
