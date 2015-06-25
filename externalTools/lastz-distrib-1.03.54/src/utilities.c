//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: utilities.c
//
//----------
//
// utilities--
//	Miscellaneous utility functions.
//
//----------

//----------
//
// other files
//
//----------

#include <stdlib.h>				// standard C stuff
#include <stdio.h>
#define  true  1
#define  false 0
#include <stdio.h>				// standard C i/o stuff
#include <string.h>				// standard C string stuff
#include <ctype.h>				// standard C upper/lower stuff
#include <stdarg.h>				// standard C variable arg list stuff
#include <limits.h>				// standard C value limit stuff
#include <errno.h>				// standard C error number stuff
#include "build_options.h"		// build options

#define  utilities_owner		// (make this the owner of its globals)
#include "utilities.h"			// interface to this module

//----------
//
// memory allocation tracking macros
//
// if trackMemoryUsage is #defined (e.g. in the makefile), reportAlloc,
// reportRealloc, and reportFree report memory usage to stderr.  If it is not
// #defined, the functions compile to nothing.
//
//----------

#ifdef trackMemoryUsage

#define reportAlloc(id,p,sz)        fprintf (stderr, "[[ memory alloc %s %08lX #%s ]]\n",          (((id)==NULL)? "" : (id)),               (long)(p), commatize(sz));
#define reportRealloc(id,oldp,p,sz) fprintf (stderr, "[[ memory realloc %s %08lX->%08lX #%s ]]\n", (((id)==NULL)? "" : (id)), (long)(oldp), (long)(p), commatize(sz));
#define reportFree(id,p)            fprintf (stderr, "[[ memory free  %s %08lX ]]\n",              (((id)==NULL)? "" : (id)),               (long)(p));

#else

#define reportAlloc(id,p,sz)        ;
#define reportRealloc(id,oldp,p,sz) ;
#define reportFree(id,p)            ;

#endif // trackMemoryUsage

//----------
//
// fopen_or_die--
//	Open a file.
//
//----------
//
// Arguments:
//	(same as for fopen())
//
// Returns:
//	A pointer to file;  failures result in fatality.
//
//----------

FILE* fopen_or_die
   (const char*	name,
	const char*	mode)
	{
	FILE*		f;

	f = fopen (name, mode);
	if (f == NULL)
		suicidef ("fopen_or_die failed to open \"%s\" for \"%s\"", name, mode);

	if (utilities_dbgDumpFilePointers)
		fprintf (stderr, "fopen_or_die(\"%s\",\"%s\") returns %p\n", name, mode, f);

	return f;
	}

//----------
//
// fclose_if_valid--
//	Close a file previously opened with fopen().
//
//----------
//
// Arguments:
//	(same as for fclose())
//
// Returns:
//	(same as for fclose())
//
//----------

int fclose_if_valid
   (FILE*	f)
	{
	if ((f == NULL) || (f == stdin) || (f == stdout) || (f == stderr))
		return 0;

	if (utilities_dbgDumpFilePointers)
		fprintf (stderr, "fclose_if_valid(%p)\n", f);

	return fclose (f);
	}

//----------
//
// getc_or_die--
//	Read a character from a file.
//
//----------
//
// Arguments:
//	FILE*	f:			(same as for getc())
//	char*	filename:	The name of the file associated with f.  This is used
//						.. only for error reporting, and may be NULL.
//
// Returns:
//	(same as for getc();  except that errors result in program fatality)
//
//----------

int getc_or_die
   (FILE*	f,
	char*	filename)
	{
	int		ch;

	ch = getc (f);
	if (ch != EOF) return ch & 0xFF;

	if (ferror (f))
		{
		if (filename == NULL) filename = "(unnamed file)";
		if (utilities_dbgDumpFilePointers)
			fprintf (stderr, "getc_or_die(%p) (filename \"%s\") reported errno=%d\n", f, filename, errno);
		suicidef_with_perror ("I/O failure for %s (getc reported errno=%d)",
		                      filename, errno);
		}

	return EOF;
	}

//----------
//
// print_prefix--
//	Print a prefix of a string.
//
//----------
//
// Arguments:
//	FILE*		f:	The file to print to.
//	const char*	s:	The string to print the prefix of.
//	int			n:	The number of characters to print.
//
// Returns:
//	The number of characters printed.  This is expected to be n.  However, if
//	the string is shorter than n, then printing terminates when the string
//	does, and the string length is returned.  Also, if n is less than 1,
//	nothing is printed and zero is returned.
//
//----------

int print_prefix
   (FILE*		f,
	const char*	s,
	int			n)
	{
	int			ix;

	if (n < 1) return 0;

	for (ix=0 ; ix<n ; ix++)
		{
		if (s[ix] == 0) break;
		fprintf (f, "%c", s[ix]);
		}

	return ix;
	}

//----------
//
// malloc_or_die, zalloc_or_die, realloc_or_die--
//	Allocate a block of memory from the heap.
//
//----------
//
// Arguments:
//	(same as for malloc() or realloc(), except for the extra id argument)
//	char* id:	an identifying string to be used when trackMemoryUsage is
//				.. turned on;  this can be NULL.
//
// Returns:
//	A pointer to new memory;  failures result in fatality.
//
//----------
//
// notes:
//
// (1)	zalloc_or_die is a malloc() replacement that fills the block of memory
//		.. with zeros.
//
//----------

#ifndef noMemoryWrappers

void* malloc_or_die
   (char*	id,
	size_t	size)
	{
	void*	p;

	// make sure size is legit

#if (SIZE_MAX > mallocLimit)
	if (size > mallocLimit)
		{
		if (id == NULL)
			suicidef ("malloc_or_die blocked large request, for %s bytes (max is %s)",
			          ucommatize(size), ucommatize(mallocLimit));
		else
			suicidef ("malloc_or_die blocked large request, for %s bytes (max is %s), for %s",
			          ucommatize(size), ucommatize(mallocLimit), id);
		}
#endif // overflow possible

	if (size == 0) size = 1;

	// allocate the memory

	p = malloc (size);
	if (p == NULL)
		{
		if (id == NULL)
			suicidef ("call to malloc failed to allocate %s bytes",
			          ucommatize(size));
		else
			suicidef ("call to malloc failed to allocate %s bytes, for %s",
			          ucommatize(size), id);
		}

	reportAlloc (id, p, size);

	return p;
	}

void* zalloc_or_die
   (char*	id,
	size_t	size)
	{
	void*	p;

	// make sure size is legit

	if (size == 0) size = 1;

	// allocate the memory and clear it

	p = malloc_or_die (id, size);
	memset (p, 0, size);

	return p;
	}

void* realloc_or_die
   (char*	id,
	void*	_p,
	size_t	size)
	{
	void*	p;

	// make sure size is legit

#if (SIZE_MAX > mallocLimit)
	if (size > mallocLimit)
		{
		if (id == NULL)
			suicidef ("realloc_or_die blocked large request, for %s bytes (max is %s)",
			          ucommatize(size), ucommatize(mallocLimit));
		else
			suicidef ("realloc_or_die blocked large request, for %s bytes (max is %s), for %s",
			          ucommatize(size), ucommatize(mallocLimit), id);
		}
#endif // overflow possible

	if (size == 0) size = 1;

	// allocate the memory

	p = realloc (_p, size);
	if (p == NULL)
		{
		if (id == NULL)
			suicidef ("call to realloc failed to allocate %s bytes",
			          ucommatize(size));
		else
			suicidef ("call to realloc failed to allocate %s bytes, for %s",
			          ucommatize(size), id);
		}

	reportRealloc (id, _p, p, size);

	return p;
	}

#endif // not noMemoryWrappers

//----------
//
// free_if_valid--
//	De-allocate a block of memory previously allocated from the heap, but
//	first checking to make sure the pointer is not NULL.
//
//----------
//
// Arguments:
//	(same as for free(), except for the extra id argument)
//	char* id:	an identifying string to be used when trackMemoryUsage is
//				.. turned on;  this can be NULL.
//
// Returns:
//	(nothing)
//
//----------

#ifndef noMemoryWrappers

void free_if_valid (arg_dont_complain(char* id), void* p)
	{
	if (p == NULL) return;

	free (p);
	reportFree (id, p);
	}

#endif // not noMemoryWrappers

//----------
//
// copy_string, copy_prefix--
//	Create (in the heap) a copy of a string or a prefix of a string.
//
//----------
//
// Arguments:
//	const char*	s:	The string to copy.
//	int			n:	(copy_prefix only) the number of characters to copy.
//
// Returns:
//	A pointer to new string;  failures result in fatality.
//
//----------

char* copy_string
   (const char*	s)
	{
	char*		ss;

	if (s == NULL) return NULL;

	ss = malloc_or_die ("copy_string", strlen(s) + 1);
	return strcpy (/*to*/ ss, /*from*/ s);
	}

char* copy_prefix
   (const char*	s,
	int			n)
	{
	char*		ss;

	if (s == NULL) return NULL;

	ss = malloc_or_die ("copy_prefix", n + 1);
	memcpy (/*to*/ ss, /*from*/ s, /*how much*/ n);
	ss[n] = 0;
	return ss;
	}

//----------
//
// concatenate_strings--
//	Create (in the heap) concatenation of two strings.
//
//----------
//
// Arguments:
//	const char*	s1, s2:	The strings to copy.
//
// Returns:
//	A pointer to new string;  failures result in fatality.
//
//----------

char* concatenate_strings
   (const char*	s1,
	const char*	s2)
	{
	char*		s, *scan;
	size_t		len = 0;

	if (s1 != NULL) len += strlen (s1);
	if (s2 != NULL) len += strlen (s2);

	s = malloc_or_die ("concatenate_strings", len + 1);

	scan = s;
	if (s1 != NULL) { strcpy (scan, s1);  scan += strlen (s1); }
	if (s2 != NULL) { strcpy (scan, s2);  scan += strlen (s2); }
	*scan = 0;

	return s;
	}

//----------
//
// append_char, append_u8--
//	Append a character to a growable string.
//
//----------
//
// Arguments:
//	type**	s:		(Pointer to) The string to append to.  The string may be
//					.. NULL.  type is either char or u8.
//	u32*	size:	(Pointer to) The number of bytes currently allocated for
//					.. s[].  This may be zero.
//	u32*	len:	(Pointer to) The number of characters currently in the
//					.. string.  Note that this routine does not consider
//					.. terminating zeros.  If the string has one it counts,
//					.. if it doesn't have one, it doesn't get counted.
//	type	ch:		The character to append.  type is either char or u8.
//
// Returns:
//	Nothing;  failures result in fatality.  However, the locations pointed to
//	by s, size, and len may be altered.
//
//----------
//
// Note:	We have here two routines that are identical except for the type of
//			.. the characters in the string (char or u8).  The reson for two
//			.. such routines is to satisfy certain compilers that would require
//			.. a cast from u8** to char** (if we only had a char version of this
//			.. function) but then complain about type-punning when such cast is
//			.. made.
//
//----------

#define create_append_function(function_name,function_string,char_type)      \
void function_name                                                           \
   (char_type**	s,                                                           \
	u32*		size,                                                        \
	u32*		len,                                                         \
	char_type	ch)                                                          \
	{                                                                        \
	/* if we don't have enough room, try to grow */                          \
                                                                             \
	if (*len >= *size)                                                       \
		{                                                                    \
		*size = *size + (*size >> 3) + 30;                                   \
		*s = realloc_or_die (function_string, *s, *size);                    \
		}                                                                    \
                                                                             \
	/* deposit the character */                                              \
                                                                             \
	(*s)[(*len)++] = ch;                                                     \
	}

create_append_function(append_char,"append_char",char)
create_append_function(append_u8,  "append_u8",  u8)

//----------
//
// strcmp_prefix--
//	Determine if a string contains another as a prefix.
//
//----------
//
// Arguments:
//	const char*	str1:	The string.
//	const char*	str2:	The prefix string.
//
// Returns:
//	The same as strcmp(prefix1,str2) would, where prefix1 is str1 truncated
//	to be no longer than str2.
//
//----------

int strcmp_prefix
   (const char*	str1,
	const char*	str2)
	{
	return strncmp (str1, str2, strlen (str2));
	}

//----------
//
// strcmp_suffix, strncmp_suffix--
//	Determine if a string contains another as a suffix.
//
//----------
//
// Arguments:
//	const char*	str1:	The string.
//	const char*	str2:	The suffix string.
//	size_t		n:		(strncmp_suffix only) The max length of str1.
//
// Returns:
//	The same as strcmp(suffix1,str2) or strncmp(suffix1,str2,n) would, where
//	suffix1 is the last N characters of str1, and N is the length of str2.  If
//	str2 is longer than str1, it cannot be a suffix (in this case we compare to
//	the entirety of str1).
//
//----------

int strcmp_suffix
   (const char*	str1,
	const char*	str2)
	{
	size_t len1 = strlen(str1);
	size_t len2 = strlen(str2);

	if (len2 <= len1) return strcmp (str1+len1-len2, str2);
	             else return strcmp (str1,           str2);
	}

int strncmp_suffix
   (const char*	str1,
	const char*	str2,
	size_t		n)
	{
	size_t len1 = strlen(str1);
	size_t len2 = strlen(str2);

	if (len1 > n) len1 = n;

	if (len2 <= len1) return strcmp (str1+len1-len2, str2);
	             else return strcmp (str1,           str2);
	}

//----------
//
// is_blank_string--
//	Determine if a string contains only blank characters.
//
//----------
//
// Arguments:
//
//	const char*	s:	The string.
//
// Returns:
//	true if the string contains only characters for which isspace is true (e.g.
//	spaces, tabs, line feeds);  false otherwise
//
//----------

int is_blank_string
   (const char*	s)
	{
	char*		ss = (char*) s;

	while (*ss != 0)
		{ if (!isspace(*(ss++))) return false; }

	return true;
	}

//----------
//
// string_to_int--
//	Parse a string for the integer value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything other than a valid integer-- failures result in fatality.
//
//----------

int string_to_int
   (const char*	s)
	{
	char*		ss;
	int			v;
	char		extra;

	// skip to first non-blank

	ss = (char*) s;
	while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
		ss++;
	if (*ss == 0) goto empty_string;

	// convert to number

	if (sscanf (ss, "%d%c", &v, &extra) != 1) goto not_an_integer;

	// make sure signs match

	if ((v < 0) && (*ss != '-')) goto out_of_range;
	if ((v > 0) && (*ss == '-')) goto out_of_range;

	return v;

	//////////
	// failure exits
	//////////

empty_string:
	suicidef ("an empty string is not an integer");

not_an_integer:
	suicidef ("\"%s\" is not an integer", s);

out_of_range:
	suicidef ("\"%s\" is outside the range of a signed integer", s);

	return 0;
	}

//----------
//
// string_to_unitized_int, string_to_unitized_int64--
//	Parse a string for the integer value it contains, allowing K, M, and G
//	suffixes.
//
//----------
//
// Arguments:
//	const char*	s:		The string to parse.
//	int byThousands:	true  => K means one thousand
//						false => K means 1,024.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything (except for an opptional suffix) other than a valid integer--
//	failures result in fatality.
//
//----------

int string_to_unitized_int
   (const char*	s,
	int			byThousands)
	{
	char		ss[20];
	int			len = strlen (s);
	char*		parseMe;
	int			v;
	float		vf;
	char		extra;
	int			mult;
	int			isFloat;

	mult = 1;

	if (len >= (int) sizeof (ss))
		parseMe = (char*) s;
	else
		{
		parseMe = ss;
		strcpy (ss, s);

		if (len > 0)
			{
			switch (ss[len-1])
				{
				case 'K': case 'k':
					mult = (byThousands)? 1000 : 1024;
					break;
				case 'M': case 'm':
					mult = (byThousands)? 1000000 : 1024L * 1024L;
					break;
				case 'G': case 'g':
					mult = (byThousands)? 1000000000 : 1024L * 1024L * 1024L;
					break;
				}

			if (mult != 1)
				ss[len-1] = 0;
			}
		}

	isFloat = false;
	if (sscanf (parseMe, "%d%c", &v, &extra) != 1)
		{
		if (sscanf (parseMe, "%f%c", &vf, &extra) != 1) goto bad;
		isFloat = true;
		}

	if (isFloat)
		{
		if ((vf > 0) && ( vf*mult > INT_MAX)) goto overflow;
		if ((vf < 0) && (-vf*mult > INT_MAX)) goto overflow;
		v = (vf * mult) + .5;
		}
	else if (mult != 1)
		{
		if ((v > 0) && ( v > INT_MAX / mult)) goto overflow;
		if ((v < 0) && (-v > INT_MAX / mult)) goto overflow;
		v *= mult;
		}

	return v;

bad:
	suicidef ("\"%s\" is not an integer", s);
	return 0;

overflow:
	suicidef ("\"%s\" is out of range for an integer", s);
	return 0;
	}


int64 string_to_unitized_int64
   (const char*	s,
	int			byThousands)
	{
	char		ss[20];
	int			len = strlen (s);
	char*		parseMe;
	int64		v;
	float		vf;
	char		extra;
	int64		mult;
	int			isFloat;

	mult = 1;

	if (len >= (int) sizeof (ss))
		parseMe = (char*) s;
	else
		{
		parseMe = ss;
		strcpy (ss, s);

		if (len > 0)
			{
			switch (ss[len-1])
				{
				case 'K': case 'k':
					mult = (byThousands)? 1000 : 1024;
					break;
				case 'M': case 'm':
					mult = (byThousands)? 1000000 : 1024L * 1024L;
					break;
				case 'G': case 'g':
					mult = (byThousands)? 1000000000 : 1024L * 1024L * 1024L;
					break;
				}

			if (mult != 1)
				ss[len-1] = 0;
			}
		}

	isFloat = false;
	if (sscanf (parseMe, s64Fmt "%c", &v, &extra) != 1)
		{
		if (sscanf (parseMe, "%f%c", &vf, &extra) != 1) goto bad;
		isFloat = true;
		}

	if (isFloat)
		{
		if ((vf > 0) && ( vf*mult > s64max)) goto overflow;
		if ((vf < 0) && (-vf*mult > s64max)) goto overflow;
		v = (vf * mult) + .5;
		}
	else if (mult != 1)
		{
		if ((v > 0) && ( v > s64max / mult)) goto overflow;
		if ((v < 0) && (-v > s64max / mult)) goto overflow;
		v *= mult;
		}

	return v;

bad:
	suicidef ("\"%s\" is not an integer", s);
	return 0;

overflow:
	suicidef ("\"%s\" is out of range for an integer", s);
	return 0;
	}

//----------
//
// hex_string_to_int--
//	Parse a string for the hexadecimal integer value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The integer value of the string.  Note that the string *must not* contain
//	anything other than a valid integer-- failures result in fatality.
//
//----------

int hex_string_to_int
   (const char*	s)
	{
	int			v;
	char		extra;

	if (sscanf (s, "%X%c", &v, &extra) != 1)
		suicidef ("\"%s\" is not an integer", s);

	return v;
	}

//----------
//
// string_to_double--
//	Parse a string for the double floating point value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The value of the string.  Note that the string *must not* contain anything
//	other than a valid number-- failures result in fatality.
//
//----------

double string_to_double
   (const char*	s)
	{
	double		v;
	char		extra;

	if (sscanf (s, "%lf%c", &v, &extra) != 1)
		suicidef ("\"%s\" is not a number", s);

	return v;
	}

//----------
//
// string_to_unitized_double--
//	Parse a string for the floating point value it contains, allowing K, M, and
//	G suffixes.
//
//----------
//
// Arguments:
//	const char*	s:		The string to parse.
//	int byThousands:	true  => K means one thousand
//						false => K means 1,024.
//
// Returns:
//	The value of the string.  Note that the string *must not* contain anything
//	other than a valid number-- failures result in fatality.
//
//----------

double string_to_unitized_double
   (const char*	s,
	int			byThousands)
	{
	char		ss[20];
	int			len = strlen (s);
	char*		parseMe;
	double		v;
	char		extra;
	int			mult;

	mult = 1;

	if (len >= (int) sizeof (ss))
		parseMe = (char*) s;
	else
		{
		parseMe = ss;
		strcpy (ss, s);

		if (len > 0)
			{
			switch (ss[len-1])
				{
				case 'K': case 'k':
					mult = (byThousands)? 1000 : 1024;
					break;
				case 'M': case 'm':
					mult = (byThousands)? 1000000 : 1024L * 1024L;
					break;
				case 'G': case 'g':
					mult = (byThousands)? 1000000000 : 1024L * 1024L * 1024L;
					break;
				}

			if (mult != 1)
				ss[len-1] = 0;
			}
		}

	if (sscanf (parseMe, "%lf%c", &v, &extra) != 1)
		suicidef ("\"%s\" is not a number", s);

	return v * mult;
	}

//----------
//
// pct_string_to_double--
//	Parse a percentage string for the double floating point value it contains.
//
//----------
//
// Arguments:
//	const char*	s:	The string to parse.
//
// Returns:
//	The value of the string.  Note that the string *must not* contain anything
//	other than a valid percentage-- failures result in fatality.
//
//----------

double pct_string_to_double
   (const char*	s)
	{
	double		v;
	char		pct, extra;

	if ((sscanf (s, "%lf%c%c", &v, &pct, &extra) != 2)
	 || (pct != '%'))
		suicidef ("\"%s\" is not a percentage", s);

	return v / 100.0;
	}

//----------
//
// is_valid_lastz_version--
//	Determine if a string is in proper format to be a lastz version number.
//
//----------
//
// Arguments:
//	char*	s:	The string to test.
//
// Returns:
//	true if the string is valid;  false if not.
//
//----------
//
// notes:
//
// (1)	We don't check whether the version ever existed, only that it is a
//		properly formatted version string.  This means that it is of the form
//		<number>.<number>[.<number>].
//
//----------

int is_valid_lastz_version
   (char*	s)
	{
	int		major, minor, subminor;
	char	extra;

	if (sscanf (s, "%d.%d%c", &major, &minor, &extra) == 2)
		{
		if (major < 0) return false;
		if (minor < 0) return false;
		return true;
		}

	if (sscanf (s, "%d.%d.%d%c", &major, &minor, &subminor, &extra) == 3)
		{
		if (major    < 0) return false;
		if (minor    < 0) return false;
		if (subminor < 0) return false;
		return true;
		}

	return false;
	}

//----------
//
// is_later_lastz_version--
//	Determine if a string is in proper format to be a lastz version number.
//
//----------
//
// Arguments:
//	char*	s1, s2:	The strings to test.  It is assumed that these are both
//					.. valid version strings (as per is_valid_lastz_version).
//
// Returns:
//	true if s1 is a later version than s2;  false if not.
//
//----------

int is_later_lastz_version
   (char*	s1,
	char*	s2)
	{
	int		major1, minor1, subminor1;
	int		major2, minor2, subminor2;
	char	extra;

	if (sscanf (s1, "%d.%d%c", &major1, &minor1, &extra) == 2)
		subminor1 = 0;
	else
		sscanf (s1, "%d.%d.%d", &major1, &minor1, &subminor1);

	if (sscanf (s2, "%d.%d%c", &major2, &minor2, &extra) == 2)
		subminor2 = 0;
	else
		sscanf (s2, "%d.%d.%d", &major2, &minor2, &subminor2);

	if      (major1    > major2)    return true;
	else if (major1    < major2)    return false;

	if      (minor1    > minor2)    return true;
	else if (minor1    < minor2)    return false;

	return (subminor1 > subminor2);
	}

//----------
//
// commatize, ucommatize--
//	Convert an integer to a string, including commas.
//
//----------
//
// Arguments:
//	const int64 v:	The number to convert.
//
// Returns:
//	A string representing that number, including commas.  (see note 1)
//
//----------
//
// notes:
//
// (1)	The memory containing the returned string belongs to this routine, as
//		static memory.  There are only five such memory blocks, and they are
//		used on alternate calls.  So when you make more than five calls, the
//		results of previous calls are clobbered.
//
//----------

char* commatize
   (const int64	v)
	{
	static char	 s1[53];// (big enough for 128-bit decimal value with commas,
	static char	 s2[53];//  .. the biggest being
	static char	 s3[53];//  .. -170,141,183,460,469,231,731,687,303,715,884,105,728)
	static char	 s4[53];
	static char	 s5[53];
	static char* s = s5;
	int		len, commas;
	char*	src, *dst;

	if      (s == s1) s = s2;	// (ping pong)
	else if (s == s2) s = s3;
	else if (s == s3) s = s4;
	else if (s == s4) s = s5;
	else              s = s1;

	sprintf (s, "%jd", (intmax_t) v);	// $$$ this could overflow the buffer
										// $$$ .. if int_max_t > 128 bits

	len = strlen (s);

	if (s[0] == '-') commas = (len-2) / 3;
	            else commas = (len-1) / 3;

	if (commas != 0)
		{
		src = s + len - 1;
		dst = s + len + commas;  *(dst--) = 0;

		while (dst > src)
			{
			*(dst--) = *(src--);
			*(dst--) = *(src--);
			*(dst--) = *(src--);
			*(dst--) = ',';
			}

		}

	return s;
	}


char* ucommatize
   (const u64	v)
	{
	static char	 s1[52];// (big enough for 128-bit decimal value with commas,
	static char	 s2[52];//  .. the biggest being
	static char	 s3[52];//  .. 340,282,366,920,938,463,463,374,607,431,768,211,455)
	static char	 s4[52];
	static char	 s5[52];
	static char* s = s5;
	int		len, commas;
	char*	src, *dst;

	if      (s == s1) s = s2;	// (ping pong)
	else if (s == s2) s = s3;
	else if (s == s3) s = s4;
	else if (s == s4) s = s5;
	else              s = s1;

	sprintf (s, "%jd", (intmax_t) v);	// $$$ this could overflow the buffer
										// $$$ .. if int_max_t > 128 bits

	len = strlen (s);
	commas = (len-1) / 3;

	if (commas != 0)
		{
		src = s + len - 1;
		dst = s + len + commas;  *(dst--) = 0;

		while (dst > src)
			{
			*(dst--) = *(src--);
			*(dst--) = *(src--);
			*(dst--) = *(src--);
			*(dst--) = ',';
			}

		}

	return s;
	}

//----------
//
// unitize--
//	Convert an integer to a string, in units of K, M, or G.
//
//----------
//
// Arguments:
//	const int64 v:		The number to convert.
//	int byThousands:	true  => K means one thousand
//						false => K means 1,024.
//
// Returns:
//	A string representing that number, expressed as units.  (see note 1)
//
//----------
//
// notes:
//
// (1)	The memory containing the returned string belongs to this routine, as
//		static memory.  There are only two such memory blocks, and they are
//		used on alternate calls.  So when you make more than two calls, the
//		results of previous calls are clobbered.
//
//----------

// SI unit prefixes (see, e.g., http://en.wikipedia.org/wiki/SI_prefix)

char* unitName[] = { "", "K", "M", "G", "T", "P", "E", "Z" };


char* unitize
   (const int64	v,
	int			byThousands)
	{
	static char	 s1[10];
	static char	 s2[10];
	static char* s = s2;
	int		sign, unit;
	int64	vv, divisor;
	float	rep;

	s = (s == s1)? s2 : s1;	// (ping pong)

	if (byThousands) divisor = 1000;
	            else divisor = 1024;


	if (v >= 0) { sign = '\0';  vv = v;  }
	       else { sign = '-';   vv = -v; }

	unit = 0;
	for (rep=vv ; vv>1023 ; vv/=divisor,rep/=divisor)
		unit++;

	if (rep > 99) { rep /= divisor;  unit++; }

	if (sign < 0) sprintf (s, "-%.1f%s", rep, unitName[unit]);
	         else sprintf (s, "%.1f%s",  rep, unitName[unit]);

	return s;
	}

//----------
//
// hex_64_string--
//	Convert an integer to a 64-bit hexadecimal string.
//
//----------
//
// Arguments:
//	const int64 v:	The number to convert.
//
// Returns:
//	A string representing that number.  (see note 1)
//
//----------
//
// notes:
//
// (1)	The memory containing the returned string belongs to this routine, as
//		static memory.  There are only two such memory blocks, and they are
//		used on alternate calls.  So when you make more than two calls, the
//		results of previous calls are clobbered.
//
//----------

char* hex_64_string
   (const int64	v)
	{
	static char	 s1[17];
	static char	 s2[17];
	static char* s = s2;
	u64		vv = (u64) v;
	char*	dst;

	s = (s == s1)? s2 : s1;	// (ping pong)

	dst = s + sizeof(s1);
	*(--dst) = 0;
	while (dst > s)
		{
		*(--dst) = "0123456789ABCDEF"[((u8) vv) & 0xF];
		vv >>= 4;
		}

	return s;
	}

//----------
//
// prob_to_string--
//	Convert a proability to a 3-digit string.
//
//----------
//
// Arguments:
//	double	p:	The probability value.
//
// Returns:
//	The string, which always contains 3 characters (plus a terminating zero).
//
//----------

char3 prob_to_string
   (double	p)
	{
	char3	s;
	char	field[5];	// "0.xx" plus a terminator

	if      (p >  1.0)   strcpy (/*to*/ s.s, /*from*/ ">??");
	else if (p >= 0.995) strcpy (/*to*/ s.s, /*from*/ " 1 ");
	else if (p <  0.005) strcpy (/*to*/ s.s, /*from*/ " ~~");
	else if (p <  0.0)   strcpy (/*to*/ s.s, /*from*/ "<??");
	else
		{
		sprintf (field, "%.2f", p);
		strcpy (/*to*/ s.s, /*from*/ field+1);
		}

	return s;
	}

//----------
//
// string_replace--
//	In a string, replace the first instance of one substring with some other
//	string.
//
//----------
//
// Arguments:
//	char*	s:		The working string.
//	int		len:	The size of the buffer allocated for the string, including
//					.. space for a terminating zero.
//	char*	sub:	The substring to replace.
//	char*	rep:	The replacement string (this may be NULL).
//
// Returns:
//	true if the replacement has occurred, false if not.  Note that failure may
//	occur because the substring is not found *or* because there is not room to
//	make the substitution.
//
//----------

int string_replace
   (char*	s,
	int		len,
	char*	sub,
	char*	rep)
	{
	int		sLen, subLen, repLen;
	char*	pos, *src, *dst;

	// locate the substring

	pos = strstr (s, sub);
	if (pos == NULL) return false;

	// see if we have enough room

	if (rep == NULL) rep = "";

	sLen   = strlen(s);
	subLen = strlen(sub);
	repLen = strlen(rep);

	if (sLen + subLen - repLen >= len-1)
		return false;

	// if the replacement is smaller, move the tail toward the start

	if (repLen < subLen)
		{
		src = pos + subLen;
		dst = pos + repLen;
		while (*src != 0) *(dst++) = *(src++);
		*dst = 0;
		}

	// if the replacement is larger, move the tail toward the end

	else if (repLen > subLen)
		{
		src = s   + sLen;
		dst = src + (repLen - subLen);
		*dst = 0;
		while (src >= pos+subLen)
			*(--dst) = *(--src);
		}

	// copy the replacement (note that if the replacement is the same size,
	// we haven't moved the tail at all)

	memcpy (/*to*/ pos, /*from*/ rep, /*how much*/ repLen);

	return true;
	}

//----------
//
// trim_string--
//	Remove blanks (and end-of-line) from both ends of a string.
//
//----------
//
// Arguments:
//	char*	s:	The string.
//
// Returns:
//	The string (the same as s).  Leading blanks are removed by copying
//	characters forward.  Trailing blanks are removed by depositing a
//	terminating zero.
//
//----------

char* trim_string
   (char*	s)
	{
	char*	ss, *dd, *lastInk;

	// skip to first non-blank

	ss = s;
	while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
		ss++;

	if 	(*ss == 0) // (string has nothing but blanks)
		{ *s = 0;  return s; }

	// copy the rest of the string (except the terminating zero)

	dd = lastInk = s;
	while (*ss != 0)
		{
		*(dd++) = *(ss++);

		if ((*ss != 0) && (*ss != ' ') && (*ss != '\t') && (*ss != '\n'))
			lastInk = dd;
		}

	// poke a terminating zero just past the last non-blank

	lastInk[1] = 0;

	return s;
	}

//----------
//
// skip_whitespace--
//	Skip characters until we get something that ain't whitespace.
// skip_darkspace--
//	Skip characters until we get something that ain't darkspace.
// skip_til--
//	Skip characters until we get something in a specified set of characters.
// skip_while--
//	Skip characters until we get something that ain't in a specified set of
//	characters.
//
//----------
//
// Arguments:
//	char*	s:		The sequence to read.
//	char*	chars:	(if needed) The set of characters.
//
// Returns:
//	Pointer to the first character at or beyond s that meets the stopping
//	criteria.  Note that we never scan beyond the end of the string.
//
//----------

char* skip_whitespace (char* s)
	{ while ((*s != 0) && (isspace (*s))) s++;  return s; }

char* skip_darkspace (char* s)
	{ while ((*s != 0) && (!isspace (*s))) s++;  return s; }

char* skip_til (char* s, char* chars)
	{ while ((*s != 0) && (strchr (chars, *s) == NULL)) s++;  return s; }

char* skip_while (char* s, char* chars)
	{ while ((*s != 0) && (strchr (chars, *s) != NULL)) s++;  return s; }

//----------
//
// find_tabbed_tag--
//	Locate a tag in a tab-delimited tagged string.
//
// An example of a tagged string is "ID:TRWFT\tSM:BGDNCSA32".  The two tags are
// ID and SM and their respective values are TRWFT and BGDNCSA32.
//
//----------
//
// Arguments:
//	char*	s:		The tab-delimited tagged string to search.
//	char*	tag:	The tag to search for.  For example, "ID".  Note that the
//					.. colon should NOT be included in the tag.
//
// Returns:
//	Pointer to the first character of the tag, for example to the "I" or "S" in
//	the example given above.  If the tag is not properly found, NULL is returned.
//
//----------
//
// notes:
//
// (1)	If the same tag occurs more than once, we only find the first instance
//		.. of the tag.
//
//----------

char* find_tabbed_tag
   (char*	s,
	char*	tag)
	{
	char*	t;

	t = s;
	while (true)
		{
		if (*t == 0) return NULL;

		t = strstr (t, tag);
		if (t == NULL) return NULL;

		if (t[2] != ':')                 { t++;  continue; }
		if ((t != s) && (t[-1] != '\t')) { t++;  continue; }
		break;
		}

	return t;
	}

//----------
//
// tabbed_tag_length--
//	Determine the length of a tag in a tab-delimited tagged string.
//
//----------
//
// Arguments:
//	char*	tag:	The tag.
//
// Returns:
//	The number of characters in the tag, including the tag's name, colon and
//	value, but not any terminating character(s).
//
//----------

int tabbed_tag_length
   (char*	tag)
	{
	char*	t;

	t = strchr (tag, '\t');
	if (t != NULL) return t - tag;
	          else return strlen(tag);
	}

//----------
//
// swap_64_halves, swap_two32_endian--
//	Perform endian-type shuffling on 64-bit or 32-bit values.
//
// swap_64_halves:     ABCDEFGH IJKLMNOP --> IJKLMNOP ABCDEFGH
// swap_two32_endian:  ABCDEFGH IJKLMNOP --> GHEFCDAB OPMNKLIJ
// swap_32_endian:     ABCDEFGH          --> GHEFCDAB
//
//----------
//
// Arguments:
//	u64/u32 v: The value to shuffle.
//
// Returns:
//	The shuffled value.
//
//----------

u64 swap_64_halves (const u64 v)
	{
	u32 a = (u32) (v >> 32);
	u32 b = (u32)  v;

	return (((u64) b) << 32) + a;
	}

u64 swap_two32_endian (const u64 v)
	{
	u32 a = (u32) (v >> 32);
	u32 b = (u32)  v;

	return (((u64) swap_32_endian(a)) << 32) + swap_32_endian(b);
	}

u32 swap_32_endian (const u32 v)
	{
	return (( v        & 0x000000FF) << 24)
	     + (((v >>  8) & 0x000000FF) << 16)
	     + (((v >> 16) & 0x000000FF) <<  8)
	     +  ((v >> 24) & 0x000000FF);
	}

//----------
//
// bit_count, bit_count64, bit_count16, bit_count8--
//	Count the '1' bits in a 32-bit value.
//
//----------
//
// Arguments:
//	<type>	bits:	The value to count the '1' bits of.
//
// Returns:
//	The number of bits that are '1'.
//
//----------
//
// Notes:
//
//	(1)	This algorithm was adapated from one written by Glenn C. Rhoads
//		<rhoads470@my-deja.com> of the Computer Science Deptartment at Rutgers.
//
//----------

int bit_count
   (u32			bits)
	{
	const u32	allBits      = ~0L;
	const u32	mask10       = (allBits/ 3) << 1;
	const u32	mask0011     =  allBits/ 5;
	const u32	mask00001111 =  allBits/17;

	// convert each pair to a count in the range 0..2
	//  00 => 00
	//  01 => 01
	//  10 => 01
	//  11 => 10

	bits -= (bits & mask10) >> 1;

	// convert each nybble to a count in the range 0..4

	bits = (bits & mask0011) + ((bits>>2) & mask0011);

	// convert each byte to a count in the range 0..8

	bits = (bits + (bits >> 4)) & mask00001111;

	// sum counts over bytes, then 16-bit words

	bits += bits >> 8;
	bits += bits >> 16;

	return bits & 0x000000FF;
	}


int bit_count_64
   (u64			bits)
	{
	const u64	allBits      = ~0L;
	const u64	mask10       = (allBits/ 3) << 1;
	const u64	mask0011     =  allBits/ 5;
	const u64	mask00001111 =  allBits/17;

	// convert each pair to a count in the range 0..2
	//  00 => 00
	//  01 => 01
	//  10 => 01
	//  11 => 10

	bits -= (bits & mask10) >> 1;

	// convert each nybble to a count in the range 0..4

	bits = (bits & mask0011) + ((bits>>2) & mask0011);

	// convert each byte to a count in the range 0..8

	bits = (bits + (bits >> 4)) & mask00001111;

	// sum counts over bytes, then 16-bit words, then 32-bit words

	bits += bits >> 8;
	bits += bits >> 16;
	bits += bits >> 32;

	return bits & 0x000000FF;
	}


int bit_count_16
   (u32			bits)
	{
	const u32	allBits      = ~0;
	const u32	mask10       = (allBits/ 3) << 1;
	const u32	mask0011     =  allBits/ 5;
	const u32	mask00001111 =  allBits/17;

	// convert each pair to a count in the range 0..2
	//  00 => 00
	//  01 => 01
	//  10 => 01
	//  11 => 10

	bits -= (bits & mask10) >> 1;

	// convert each nybble to a count in the range 0..4

	bits = (bits & mask0011) + ((bits>>2) & mask0011);

	// convert each byte to a count in the range 0..8

	bits = (bits + (bits >> 4)) & mask00001111;

	// sum counts over the two bytes

	bits += bits >> 8;

	return bits & 0x00FF;
	}


int bit_count_8
   (u8			bits)
	{
	const u8	allBits      = ~0;
	const u8	mask10       = (allBits/ 3) << 1;
	const u8	mask0011     =  allBits/ 5;
	const u8	mask00001111 =  allBits/17;

	// convert each pair to a count in the range 0..2
	//  00 => 00
	//  01 => 01
	//  10 => 01
	//  11 => 10

	bits -= (bits & mask10) >> 1;

	// convert each nybble to a count in the range 0..4

	bits = (bits & mask0011) + ((bits>>2) & mask0011);

	// convert each byte to a count in the range 0..8

	bits = (bits + (bits >> 4)) & mask00001111;

	return bits;
	}

//----------
//
// hassock_hash--
//	Compute a variant of Austin Appleby's MurmurHash2.
//
//----------
//
// Arguments:
//	const void*	key:	The data block to hash.
//	u32			len:	The length of that block.
//
// Returns:
//	A hash of the block.
//
//----------
//
// Notes:
//
//	(1)	As of Apr/2009, information about this hash function can be found at
//		  murmurhash.googlepages.com
//	(2) This implementation is based on an implementation found at
//	      murmurhash.googlepages.com/MurmurHashNeutral2.cpp
//	    It differs in the following ways:
//	      (a) The "seed" is hardwired.
//		  (b) We parse the data block in reverse;  this allows the caller to
//            prepend an additional seed pattern to his buffer, potentially
//	          getting better mixing for the bits in the final incorporated
//	          bytes.
//		  (c) The last three bytes are incorporated in a different order than
//	          they were in MurmurHash2, because the code just works out better
//	          this way.
//
//----------

u32 hassock_hash
   (const void*	key,
	u32			len)
	{
	const u32	seed = 0x5C3FC4D3;
	const u32	m    = 0x87C10417;
	const int	r    = 24;
	const u8*	data = ((const u8*) key) + len;
	const u8*	stop = ((const u8*) key) + 4;
	u32			h, k;

	h = seed ^ len;
	while (data >= stop)
		{
		k  = *(--data);
		k |= *(--data) << 8;
		k |= *(--data) << 16;
		k |= *(--data) << 24;

		k *= m;
		k ^= k >> r;
		k *= m;

		h *= m;
		h ^= k;

		len  -= 4;
		}

	switch (len)
		{
		case 3: h ^= *(--data) << 16;
		case 2: h ^= *(--data) << 8;
		case 1: h ^= *(--data);
		        h *= m;
		};

	h ^= h >> 13;
	h *= m;
	h ^= h >> 15;

	//printf ("%08X %s\n", h, (char*) key);

	return h;
	}

//----------
//
// suicide, suicidef, suicide_with_perror, suicidef_with_perror--
//	Cause program fatality, after pushing a message out to the user.
//
//----------
//
// Arguments for suicide():
//	const char*	message:	The message to write to stderr before death.  This
//							.. may be NULL.
//
// Arguments for suicidef():
//	const char*	format:		A format string, as per printf.  This may be NULL.
//	...:					(same as for printf)
//
// Returns:
//	(nothing;  it does not return).
//
//----------

void suicide
   (const char*	message)
	{
	if (message == NULL) suicidef (NULL, NULL);
	                else suicidef ("%s", message);
	}

void suicidef
   (const char*	format,
	...)
	{
	va_list	args;

	va_start (args, format);

	fflush  (stdout);
	fprintf (stderr, "FAILURE: ");
	if (format != NULL)
		{
		vfprintf (stderr, format, args);
		fprintf  (stderr, "\n");
		}

	va_end (args);

	exit (EXIT_FAILURE);
	}

// _with_perror adds a call to the system routine perror()

void suicide_with_perror
   (const char*	message)
	{
	if (message == NULL) suicidef_with_perror (NULL, NULL);
	                else suicidef_with_perror ("%s", message);
	}

void suicidef_with_perror
   (const char*	format,
	...)
	{
	va_list	args;

	va_start (args, format);

	fflush  (stdout);
	fprintf (stderr, "FAILURE: ");
	if (format != NULL)
		{
		vfprintf (stderr, format, args);
		fprintf  (stderr, "\n");
		}

	va_end (args);

	perror ("file I/O error");

	exit (EXIT_FAILURE);
	}

