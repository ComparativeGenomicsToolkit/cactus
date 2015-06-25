//-------+---------+---------+---------+---------+---------+---------+--------=
//
// covered_intervals.c-- read a list of nearly-sorted alignment intervals and
//                       report intervals that are covered by at least some
//                       specified number of alignments
//
//----------

#include <stdlib.h>
#define  true  1
#define  false 0
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <limits.h>

#include <inttypes.h>
#include <stdint.h>
typedef int8_t   s8;
typedef uint8_t  u8;
typedef int32_t  s32;
typedef uint32_t u32;

// program revision vitals (not the best way to do this!))

#define programVersionMajor    "0"
#define programVersionMinor    "0"
#define programVersionSubMinor "3"
#define programRevisionDate    "20131202"

//----------
//
// global data and types--
//
//----------

// linked list for chromosomes seen

typedef struct info
    {
    struct info* next;          // next item in a linked list
    char*        chrom;         // chromosome name
    u32          lineNumber;    // line number where this chromosome first seen
    } info;

// command line options

info* chromsSeen      = NULL;
u32   windowSize      = 1*1000*1000;
int   inputHasOffsets = false;
int   originOne       = false;
int   endComment      = false;
int   reportChroms    = false;
u8    depthThreshold  = 1;

#define maxDepth 255

int   debugReportInputIntervals  = false;
int   debugReportParsedIntervals = false;
int   debugWindowSlide           = false;

//----------
//
// prototypes--
//
//----------

int main (int argc, char** argv);

// private functions

static void  parse_options       (int _argc, char** _argv);
static void  emit_intervals      (FILE* f, u8 minDepth,
                                  u8* window, char* chrom,
                                  u32 pendingRun, u32 windowStart, u32 windowEnd);
static u32   emit_some_intervals (FILE* f, u8 minDepth,
                                  u8* window, char* chrom,
                                  u32 pendingRun, u32 windowStart, u32 windowEnd);
static info* find_chromosome     (char* chrom);
static int   read_alignment      (FILE* f,
                                  char* buffer, int bufferLen,
                                  u32* lineNumber,
                                  char** rChrom, u32* rStart, u32* rEnd,
                                  char** qChrom, u32* qStart, u32* qEnd);

static char*  copy_string            (const char* s);
static int    strcmp_prefix          (const char* str1, const char* str2);
static int    string_to_u32          (const char* s);
static int    string_to_unitized_int (const char* s, int byThousands);
static char*  skip_whitespace        (char* s);
static char*  skip_darkspace         (char* s);

//----------
//
// option parsing--
//
//----------

static void  usage    (char* message);
static void  chastise (const char* format, ...);

char* programName = "covered_intervals";


static void usage
   (char*   message)
    {
    if (message != NULL) fprintf (stderr, "%s\n", message);
    fprintf (stderr, "usage: %s [options]\n", programName);
    fprintf (stderr, "\n");
    fprintf (stderr, "Read a list of nearly-sorted alignment intervals and report intervals that are\n");
    fprintf (stderr, "covered by at least some specified number of alignments.\n");
    fprintf (stderr, "\n");
    //                123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789
    fprintf (stderr, "  M=<depth>              report any position that is covered by at least this\n");
    fprintf (stderr, "                         many alignments; the maximum allowed depth is 255\n");
    fprintf (stderr, "                         (by default this is 1)\n");
    fprintf (stderr, "  W=<length>             size of internal bitmap \"window\", in bases;  this\n");
    fprintf (stderr, "                         should be at least twice the size of the query\n");
    fprintf (stderr, "                         fragments being aligned/reported\n");
    fprintf (stderr, "                         (by default this is 1M)\n");
    fprintf (stderr, "  --queryoffsets         input query names contain offsets, as described below\n");
    fprintf (stderr, "                         (by default input query names do not contain offsets)\n");
    fprintf (stderr, "  --origin=zero          *output* intervals are origin-zero, half-open\n");
    fprintf (stderr, "                         (this is the default)\n");
    fprintf (stderr, "  --origin=one           *output* intervals are origin-one, closed\n");
    fprintf (stderr, "                         (*input* intervals are *always* origin-zero)\n");
    fprintf (stderr, "  --markend              write a comment at the end of the output file\n");
    fprintf (stderr, "  --progress=chromosome  report each chromosome as we encounter it\n");
    fprintf (stderr, "  --version              report the program version and quit\n");
    fprintf (stderr, "\n");
    fprintf (stderr, "We expect as input alignments of the form\n");
    fprintf (stderr, "  <refchrom> <refstart> <refend> <qchrom>[_<offset>] <qstart+> <qend+>\n");
    fprintf (stderr, "where\n");
    fprintf (stderr, "  <refchrom> <refstart> <refend>  is the alignment interval on the reference\n");
    fprintf (stderr, "  <qchrom>   <qstart+>  <qend+>   is the alignment interval on the query,\n");
    fprintf (stderr, "                                  with qstart and qend relative to the\n");
    fprintf (stderr, "                                  positive strand\n");
    fprintf (stderr, "  <offset>                        is an (optional) offset to be added to\n");
    fprintf (stderr, "                                  <qstart+> and <qend+>;  usually this is\n");
    fprintf (stderr, "                                  the start of a fragment given to the\n");
    fprintf (stderr, "                                  aligner\n");
    exit (EXIT_FAILURE);
    }


static void chastise (const char* format, ...)
    {
    va_list args;

    va_start (args, format);
    if (format != NULL)
        vfprintf (stderr, format, args);
    va_end (args);

    usage (NULL);
    }


static void parse_options
   (int         _argc,
    char**      _argv)
    {
    int         argc;
    char**      argv;
    char*       arg, *argVal;
    int         tempInt;

    // skip program name

    programName = _argv[0];
    argv = _argv+1;  argc = _argc - 1;

    //////////
    // scan arguments
    //////////

    while (argc > 0)
        {
        arg    = argv[0];
        argVal = strchr(arg,'=');
        if (argVal != NULL) argVal++;

        // M=<depth>

        if ((strcmp_prefix (arg, "M=")       == 0)
         || (strcmp_prefix (arg, "--M=")     == 0)
         || (strcmp_prefix (arg, "--depth=") == 0))
            {
            tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
            if (tempInt == 0)
                chastise ("depth threshold can't be 0 (\"%s\")\n", arg);
            if (tempInt < 0)
                chastise ("depth threshold can't be negative (\"%s\")\n", arg);
            if (tempInt > maxDepth)
                chastise ("depth threshold can't be more than %d (\"%s\")\n", maxDepth, arg);
            depthThreshold = (u8) tempInt;
            goto next_arg;
            }

        // W=<length>

        if ((strcmp_prefix (arg, "W=")        == 0)
         || (strcmp_prefix (arg, "--W=")      == 0)
         || (strcmp_prefix (arg, "--window=") == 0))
            {
            tempInt = string_to_unitized_int (argVal, /*thousands*/ true);
            if (tempInt == 0)
                chastise ("chromosome length can't be 0 (\"%s\")\n", arg);
            if (tempInt < 0)
                chastise ("chromosome length can't be negative (\"%s\")\n", arg);
            windowSize = (u32) tempInt;
            goto next_arg;
            }

        // --queryoffsets

        if ((strcmp (arg, "--queryoffsets") == 0)
         || (strcmp (arg, "--hasoffsets") == 0))
            { inputHasOffsets = true;  goto next_arg; }

        // --origin=one, --origin=zero

        if ((strcmp (arg, "--origin=one") == 0)
         || (strcmp (arg, "--origin=1")   == 0))
            { originOne = true;  goto next_arg; }

        if ((strcmp (arg, "--origin=zero") == 0)
         || (strcmp (arg, "--origin=0")   == 0))
            { originOne = false;  goto next_arg; }

        // --markend

        if (strcmp (arg, "--markend") == 0)
            { endComment = true;  goto next_arg; }

        // --progress=chromosome

        if ((strcmp (arg, "--progress=chromosome")  == 0)
         || (strcmp (arg, "--progress=chromosomes") == 0))
            { reportChroms = true;  goto next_arg; }

        // --version

        if (strcmp (arg, "--version")  == 0)
            {
            fprintf (stderr, "%s (version %s.%s.%s released %s)\n",
                             programName,
                             programVersionMajor, programVersionMinor, programVersionSubMinor, programRevisionDate);
            exit (EXIT_SUCCESS);
            }

        // --debug options (unadvertised)

        if (strcmp (arg, "--debug=report:input") == 0)
            { debugReportInputIntervals = true;  goto next_arg; }

        if (strcmp (arg, "--debug=report:parsed") == 0)
            { debugReportParsedIntervals = true;  goto next_arg; }

        if (strcmp (arg, "--debug=slide") == 0)
            { debugWindowSlide = true;  goto next_arg; }

        // unknown -- argument

        if (strcmp_prefix (arg, "--") == 0)
            chastise ("Can't understand \"%s\"\n", arg);

    next_arg:
        argv++;  argc--;
        continue;
        }

    }

//----------
//
// main program--
//
//----------

int main
   (int     argc,
    char**  argv)
    {
    char    lineBuffer[1000];
    char    prevChrom[1000];
    u8*     window = NULL;
    u32     windowStart, pendingRun, newWindowStart, prefixSize, suffixSize;
    u32     lineNumber;
    char*   rChrom, *qChrom;
    info*   chromInfo, *nextInfo;
    u32     rStart, rEnd, qStart, qEnd, qStartOriginal, qEndOriginal;
    u32     ix;
    int     ok;

    parse_options (argc, argv);

    //////////
    // allocate memory
    //////////

    window = (u8*) malloc (windowSize);
    if (window == NULL) goto cant_allocate_window;

    //////////
    // process intervals
    //////////

    // read intervals and accumulate depth

    prevChrom[0] = 0;
    windowStart = 0;  pendingRun = 0;
    //memset (window, 0, windowSize); // (not necessary, happens later)

    while (true)
        {
        ok = read_alignment (stdin, lineBuffer, sizeof(lineBuffer), &lineNumber,
                             &rChrom, &rStart, &rEnd, &qChrom, &qStart, &qEnd);
        if (!ok) break;

        qStartOriginal = qStart;
        qEndOriginal   = qEnd;

        if (debugReportParsedIntervals)
            fprintf (stderr, "%s %u %u %s %u %u\n",
                             rChrom, rStart, rEnd, qChrom, qStart, qEnd);

        // if this is a new chromosome, emit pending intervals for the previous
        // chromosome and reset the window;  also make sure that we don't see
        // a chromsome in non-consecutive batches

        if (strcmp (qChrom, prevChrom) != 0)
            {
            if (prevChrom[0] != 0)
                emit_intervals (stdout, depthThreshold,
                                window, prevChrom, pendingRun,
                                windowStart, windowStart + windowSize);

            chromInfo = find_chromosome (qChrom);
            if (chromInfo != NULL) goto chrom_not_together;

            chromInfo = (info*) malloc (sizeof(info));
            if (chromInfo == NULL) goto cant_allocate_info;
            chromInfo->next       = chromsSeen;
            chromInfo->chrom      = copy_string (qChrom);
            chromInfo->lineNumber = lineNumber;

            if (reportChroms)
                fprintf (stderr, "progress: reading %s (line %u)\n", qChrom, lineNumber);
            strncpy (prevChrom, qChrom, sizeof(prevChrom));
            windowStart = 0;  pendingRun = 0;
            memset (window, 0, windowSize);
            }

        // ignore trivial self-alignments

        if ((strcmp (qChrom, rChrom) == 0) && (qStart == rStart) && (qEnd == rEnd))
            continue;

        // if this interval won't fit in the window, emit some intervals and
        // slide the window

        if (qStart < windowStart) goto window_slide_problem;

        if (qEnd > windowStart + windowSize)
            {
            newWindowStart = qEnd - windowSize/2;

            if (newWindowStart > windowStart + windowSize)
                {
                // there is no overlap between old window and new
                emit_intervals (stdout, depthThreshold,
                                window, qChrom, pendingRun,
                                windowStart, windowStart + windowSize);
                windowStart = newWindowStart;
                pendingRun  = 0;
                memset (window, 0, windowSize);
                if (debugWindowSlide)
                    fprintf (stderr, "moving window to %u\n", windowStart);
                }
            else
                {
                // there is some overlap between old window and new
                prefixSize = newWindowStart - windowStart;
                suffixSize = windowSize-prefixSize;
                pendingRun = emit_some_intervals (stdout, depthThreshold,
                                                  window, qChrom, pendingRun,
                                                  windowStart, newWindowStart);
                memcpy (/*to*/ window, /*from*/ window+prefixSize, suffixSize);
                memset (window+suffixSize, 0, prefixSize);
                windowStart = newWindowStart;
                if (debugWindowSlide)
                    fprintf (stderr, "sliding window to %u\n", windowStart);
                }
            }

        // mark this interval in the window
        // $$$ note that we're tracking depth;  someday we may want to add a
        // $$$ .. depth threshold instead of 0/1

        qStart -= windowStart;
        qEnd   -= windowStart;
        if (qEnd > windowSize) goto window_too_short;

        for (ix=qStart ; ix<qEnd ; ix++)
            { if (window[ix] < maxDepth) window[ix]++; }
        }

    // emit pending intervals for the final chromosome

    if (prevChrom[0] != 0)
        emit_intervals (stdout, depthThreshold,
                        window, prevChrom, pendingRun,
                        windowStart, windowStart + windowSize);

    //////////
    // success
    //////////

    free (window);

    for (chromInfo=chromsSeen ; chromInfo!=NULL ; chromInfo=nextInfo)
        {
        nextInfo = chromInfo->next;
        if (chromInfo->chrom  != NULL) free (chromInfo->chrom);
        free (chromInfo);
        }
    chromsSeen = NULL;

    if (endComment)
        printf ("# covered_intervals end-of-file\n");

    return EXIT_SUCCESS;

    //////////
    // failure exits
    //////////

cant_allocate_window:
    fprintf (stderr, "failed to allocate %d-entry counting window\n",
                     windowSize);
    return EXIT_FAILURE;

cant_allocate_info:
    fprintf (stderr, "failed to allocate %d-entry info record for %s\n",
                     (int) sizeof(info), qChrom);
    return EXIT_FAILURE;

window_slide_problem:
    fprintf (stderr, "%s %u %u is behind the sliding window (window start = %u)\n",
                     qChrom, qStartOriginal, qEndOriginal, windowStart);
    return EXIT_FAILURE;

window_too_short:
    fprintf (stderr, "W=%u is unable to accomodate %s %u %u (window start = %u)\n",
                     windowSize, qChrom, qStartOriginal, qEndOriginal, windowStart);
    return EXIT_FAILURE;

chrom_not_together:
    fprintf (stderr, "alignments for \"%s\" are not together in the input (lines %u and %u)\n",
                     qChrom, chromInfo->lineNumber, lineNumber);
    return EXIT_FAILURE;
    }

//----------
//
// emit_intervals--
//  Emit covered intervals from the sliding window.
// emit_some_intervals--
//  Emit covered intervals from a prefix of the sliding window.
//
//----------
//
// Arguments:
//  FILE*   f:              file to write to.
//  u8      minDepth:       minimum depth a position must have, to be
//                          .. considered "covered"
//  u8*     window:         depth-of-coverage vector.
//  char*   chrom:          name of the chromosome.
//  u32     pendingRun:     length of run preceding the first entry in the
//                          .. vector.
//  u32     windowStart:    position (on the chromosome) of the first entry in
//                          .. the vector.
//  u32     windowEnd:      position (on the chromosome) beyond the last entry
//                          .. in the vector;  for emit_some_intervals this is
//                          .. at the end of the prefix.
//
// Returns:
//  (emit_intervals)      nothing
//  (emit_some_intervals) length of the unreported run at the end of the prefix
//
//----------

//=== emit_intervals ===

static void emit_intervals
   (FILE*   f,
    u8      minDepth,
    u8*     window,
    char*   chrom,
    u32     pendingRun,
    u32     windowStart,
    u32     windowEnd)
    {
    u32     run;
    u32     o = (originOne)? 1:0;

    run = emit_some_intervals (f, minDepth,
                               window, chrom, pendingRun, windowStart, windowEnd);
    if (run > 0)
        fprintf (f, "%s\t%d\t%d\n", chrom, (windowEnd-run)+o, windowEnd);
    }


//=== emit_some_intervals ===

static u32 emit_some_intervals
   (FILE*   f,
    u8      minDepth,
    u8*     window,
    char*   chrom,
    u32     pendingRun,
    u32     windowStart,
    u32     windowEnd)
    {
    u32     run = pendingRun;
    u32     o = (originOne)? 1:0;
    u32     ix, pos;

    for (ix=0,pos=windowStart ; pos<windowEnd ; ix++,pos++)
        {
        if (window[ix] >= minDepth)
            run++;
        else if (run > 0)
            {
            fprintf (f, "%s\t%d\t%d\n", chrom, (pos-run)+o, pos);
            run = 0;
            }
        }

    return run;
    }

//----------
//
// find_chromosome--
//  Locate a specific chromosome name.
//
//----------
//
// Arguments:
//  char*   chrom:  name of the chromosome to look for.
//
// Returns:
//  a pointer to the record for the chromosome;  NULL if the chromosome is not
//  in our list.
//
//----------

static info* find_chromosome
   (char*   chrom)
    {
    info*   scanInfo;

    for (scanInfo=chromsSeen ; scanInfo!=NULL ; scanInfo=scanInfo->next)
        { if (strcmp (chrom, scanInfo->chrom) == 0) return scanInfo; }

    return NULL;
    }

//----------
//
// read_alignment--
//  Read the next alignment from a file.
//
// We expect alignments to be of the form
//  <refchrom> <refstart> <refend> <qchrom>[_<offset>] <qstart+> <qend+>
// where
//  <refchrom> <refstart> <refend>  is the alignment interval on the reference
//  <qchrom>   <qstart+>  <qend+>   is the alignment interval on the query,
//                                  .. with qstart and qend relative to the
//                                  .. positive strand
//  <offset>                        is an (optional) offset;  usually this is
//                                  .. the start of the fragment being aligned
//
//----------
//
// Arguments:
//  FILE*   f:          File to read from.
//  char*   buffer:     Buffer to read the line into.  Note that the caller
//                      .. should not expect anything about the contents of
//                      .. this buffer upon return.
//  int     bufferLen:  Number of bytes allocated for the buffer.
//  u32*    rStart:     Place to return the line number.
//  char**  rChrom:     Place to return a pointer to the reference chromosome.
//                      .. The returned value will point into the line buffer,
//                      .. and to a zero-terminated string.
//  u32*    rStart:     Place to return the reference start.
//  u32*    rEnd:       Place to return the reference end.
//  char**  qChrom:     Place to return a pointer to the query chromosome.  The
//                      .. returned value will point into the line buffer, and
//                      .. to a zero-terminated string.
//  u32*    qStart:     Place to return the query start.
//  u32*    qEnd:       Place to return the query end.
//
// Returns:
//  true if we were successful;  false if there are no more lines in the file.
//
//----------

static int read_alignment
   (FILE*       f,
    char*       buffer,
    int         bufferLen,
    u32*        _lineNumber,
    char**      _rChrom,
    u32*        _rStart,
    u32*        _rEnd,
    char**      _qChrom,
    u32*        _qStart,
    u32*        _qEnd)
    {
    static u32  lineNumber = 0;
    static int  missingEol = false;
    int         lineLen;
    char*       scan, *mark, *field;
    char*       rChrom, *qChrom;
    u32         rStart, rEnd, qStart, qEnd;
    u32         qOffset;

    // read the next line

try_again:

    if (fgets (buffer, bufferLen, f) == NULL)
        return false;

    lineNumber++;

    // check for lines getting split by fgets (the final line in the file might
    // not have a newline, but no internal lines can be that way)

    if (missingEol) goto missing_eol;

    lineLen = strlen(buffer);
    if (lineLen != 0)
        missingEol = (buffer[lineLen-1] != '\n');

    if (debugReportInputIntervals)
        fprintf (stderr, "line %u: %s", lineNumber, buffer);

    // parse the line

    scan = skip_whitespace(buffer);
    if (*scan == 0)   goto try_again;  // empty line
    if (*scan == '#') goto try_again;  // comment line

    rChrom = scan = buffer;
    if (*scan == ' ') goto no_ref_chrom;
    mark = skip_darkspace(scan);
    scan = skip_whitespace(mark);
    if (*mark != 0) *mark = 0;

    if (*scan == 0) goto no_ref_start;
    field = scan;
    mark = skip_darkspace(scan);
    scan = skip_whitespace(mark);
    if (*mark != 0) *mark = 0;
    rStart = string_to_u32 (field);

    if (*scan == 0) goto no_ref_end;
    field = scan;
    mark = skip_darkspace(scan);
    scan = skip_whitespace(mark);
    if (*mark != 0) *mark = 0;
    rEnd = string_to_u32 (field);

    if (*scan == 0) goto no_query_chrom;
    qChrom = scan;
    mark = skip_darkspace(scan);
    scan = skip_whitespace(mark);
    if (*mark != 0) *mark = 0;

    if (*scan == 0) goto no_query_start;
    field = scan;
    mark = skip_darkspace(scan);
    scan = skip_whitespace(mark);
    if (*mark != 0) *mark = 0;
    qStart = string_to_u32 (field);

    if (*scan == 0) goto no_query_end;
    field = scan;
    mark = skip_darkspace(scan);
    scan = skip_whitespace(mark);
    if (*mark != 0) *mark = 0;
    qEnd = string_to_u32 (field);

    // split offset off of query chromosome

    if (inputHasOffsets)
        {
        qOffset = 0;

        field = strrchr (qChrom, '_');
        if (field == NULL) goto no_offset;
        *(field++) = 0;
        qOffset = string_to_u32 (field);

        qStart += qOffset;
        qEnd   += qOffset;
        }

    //////////
    // success
    //////////

    if (_lineNumber != NULL) *_lineNumber = lineNumber;
    if (_rChrom     != NULL) *_rChrom     = rChrom;
    if (_rStart     != NULL) *_rStart     = rStart;
    if (_rEnd       != NULL) *_rEnd       = rEnd;
    if (_qChrom     != NULL) *_qChrom     = qChrom;
    if (_qStart     != NULL) *_qStart     = qStart;
    if (_qEnd       != NULL) *_qEnd       = qEnd;

    return true;

    //////////
    // failure exits
    //////////

missing_eol:
    fprintf (stderr, "problem at line %u, line is longer than internal buffer\n",
             lineNumber-1);
    exit (EXIT_FAILURE);

no_ref_chrom:
    fprintf (stderr, "problem at line %u, line contains no reference chromosome or begins with whitespace\n",
             lineNumber-1);
    exit (EXIT_FAILURE);

no_ref_start:
    fprintf (stderr, "problem at line %u, line contains no reference interval start\n",
             lineNumber-1);
    exit (EXIT_FAILURE);

no_ref_end:
    fprintf (stderr, "problem at line %u, line contains no reference interval end\n",
             lineNumber-1);
    exit (EXIT_FAILURE);

no_query_chrom:
    fprintf (stderr, "problem at line %u, line contains no query chromosome or begins with whitespace\n",
             lineNumber-1);
    exit (EXIT_FAILURE);

no_query_start:
    fprintf (stderr, "problem at line %u, line contains no query interval start\n",
             lineNumber-1);
    exit (EXIT_FAILURE);

no_query_end:
    fprintf (stderr, "problem at line %u, line contains no query interval end\n",
             lineNumber-1);
    exit (EXIT_FAILURE);

no_offset:
    fprintf (stderr, "problem at line %u, line contains no query offset\n",
             lineNumber-1);
    exit (EXIT_FAILURE);
    }

//----------
//
// copy_string--
//  Create (in the heap) a copy of a string or a prefix of a string.
//
//----------
//
// Arguments:
//  const char* s:  The string to copy.
//  int         n:  (copy_prefix only) the number of characters to copy.
//
// Returns:
//  A pointer to new string;  failures result in program termination.
//
//----------

static char* copy_string
   (const char* s)
    {
    char*       ss;

    if (s == NULL) return NULL;

    ss = malloc (strlen(s) + 1);
    if (ss == NULL)
        {
        fprintf (stderr, "failed to allocate %lld bytes to copy \"%s\"\n",
                         (long long) (strlen(s)+1), s);
        exit (EXIT_FAILURE);
        }

    return strcpy (/*to*/ ss, /*from*/ s);
    }

//----------
//
// strcmp_prefix--
//  Determine if a string contains another as a prefix.
//
//----------
//
// Arguments:
//  const char* str1:   The string.
//  const char* str2:   The prefix string.
//
// Returns:
//  The same as strcmp(prefix1,str2) would, where prefix1 is str1 truncated
//  to be no longer than str2.
//
//----------

static int strcmp_prefix
   (const char* str1,
    const char* str2)
    {
    return strncmp (str1, str2, strlen(str2));
    }

//----------
//
//  string_to_u32--
//  Parse a string for the integer value it contains.
//
//----------
//
// Arguments:
//  const char* s:  The string to parse.
//
// Returns:
//  The integer value of the string.  Note that the string *must not* contain
//  anything other than a valid integer-- failures result in program
//  termination.
//
//----------

int string_to_u32
   (const char* s)
    {
    char*       ss;
    u32         v;
    char        extra;

    // skip to first non-blank

    ss = (char*) s;
    while ((*ss == ' ') || (*ss == '\t') || (*ss == '\n'))
        ss++;
    if (*ss == 0) goto empty_string;

    // convert to number

    if (sscanf (ss, "%u%c", &v, &extra) != 1) goto not_an_integer;

    return v;

    //////////
    // failure exits
    //////////

empty_string:
    fprintf (stderr, "an empty string is not an unsigned integer\n");
    exit (EXIT_FAILURE);

not_an_integer:
    fprintf (stderr, "\"%s\" is not an unsigned integer\n", s);
    exit (EXIT_FAILURE);

    return 0;
    }

//----------
//
// string_to_unitized_int, string_to_unitized_int64--
//  Parse a string for the integer value it contains, allowing K, M, and G
//  suffixes.
//
//----------
//
// Arguments:
//  const char* s:      The string to parse.
//  int byThousands:    true  => K means one thousand
//                      false => K means 1,024.
//
// Returns:
//  The integer value of the string.  Note that the string *must not* contain
//  anything (except for an opptional suffix) other than a valid integer--
//  failures result in fatality.
//
//----------

static int string_to_unitized_int
   (const char* s,
    int         byThousands)
    {
    char        ss[20];
    int         len = strlen (s);
    char*       parseMe;
    int         v;
    float       vf;
    char        extra;
    int         mult;
    int         isFloat;

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
    fprintf (stderr, "\"%s\" is not an integer\n", s);
    exit (EXIT_FAILURE);

overflow:
    fprintf (stderr, "\"%s\" is out of range for an integer\n", s);
    exit (EXIT_FAILURE);

    return 0;
    }

//----------
//
// skip_whitespace--
//  Skip characters until we get something that ain't whitespace.
// skip_darkspace--
//  Skip characters until we get something that ain't darkspace.
//
//----------
//
// Arguments:
//  char*   s:  The string to read.
//
// Returns:
//  Pointer to the first character at or beyond s that meets the stopping
//  criteria.  Note that we never scan beyond the end of the string.
//
//----------

static char* skip_whitespace (char* s)
    { while ((*s != 0) && (isspace (*s))) s++;  return s; }

static char* skip_darkspace (char* s)
    { while ((*s != 0) && (!isspace (*s))) s++;  return s; }

