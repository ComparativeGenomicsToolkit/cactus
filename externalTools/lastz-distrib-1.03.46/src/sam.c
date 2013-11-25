//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: sam.c
//
//----------
//
// sam--
//	Support for printing alignments in SAM format.
//
// SAM format is a pairwise alignment format designed for short read alignments.
// A spec can be found at  samtools.sourceforge.net/SAM1.pdf.
//----------

//----------
//
// other files
//
//----------

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <string.h>				// standard C string stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff
#include "identity_dist.h"		// identity distribution "format" stuff
#include "coverage_dist.h"		// query coverage distribution stuff

#define  sam_owner				// (make this the owner of its globals)
#include "sam.h"				// interface to this module

// SAM bit-encoded flags

#define BAM_FPAIRED         1	// the read is paired in sequencing, no matter whether it is mapped in a pair
#define BAM_FPROPER_PAIR    2	// the read is mapped in a proper pair
#define BAM_FUNMAP          4	// the read itself is unmapped; conflictive with BAM_FPROPER_PAIR
#define BAM_FMUNMAP         8	// the mate is unmapped
#define BAM_FREVERSE       16	// the read is mapped to the reverse strand
#define BAM_FMREVERSE      32	// the mate is mapped to the reverse strand
#define BAM_FREAD1         64	// this is read1
#define BAM_FREAD2        128	// this is read2
#define BAM_FSECONDARY    256	// not primary alignment
#define BAM_FQCFAIL       512	// QC failure
#define BAM_FDUP         1024	// optical or PCR duplicate

//----------
//
// prototypes for private functions
//
//----------

static void print_query_bases (FILE* f, seq* seq2, unspos pos2, unspos length,
                               int softMasked);
static void print_query_quals (FILE* f, seq* seq2, unspos pos2, unspos length,
                               int softMasked);

//----------
//
// sam_rg_tags--
//	Parse a readgroup string and extract the tags that should be attached to
//	each read.
//
//----------
//
// Arguments:
//	char*	readGroup:	The readgroup string (e.g. "ID:TRWFT SM:BGDNCSA32").
//	char**	errorText:	Place to return any error text.  The string pointer
//						.. returned here is local to THIS routine, and does
//						.. not need to be deallocated.  This pointer may be
//						.. NULL, in which case no error text is returned.
//
// Returns:
//	A pointer to the RG tags, allocated from the heap.  The caller is
//	responsible for deallocating that string.  If there is any error in parsing
//	the string, NULLL is returned and the explaination of the error is copied
//	to errorText (unless it is NULL).
//
//----------

char* sam_rg_tags
   (char*	readGroup,
	char**	errorText)
	{
	char*	idTag;
	int		idLen;
	u32		bytesNeeded;
	char*	s, *ss;

	//////////
	// parse and validate the important fields
	//////////

	idTag = find_tabbed_tag (readGroup, "ID");
	if (idTag == NULL)
		{
		if (errorText != NULL) *errorText = "ID is a required field";
		return NULL;
		}
	idLen = tabbed_tag_length (idTag);
	if (idLen <= 3)
		{
		if (errorText != NULL) *errorText = "ID field cannot be empty";
		return NULL;
		}

	// nota bene: SAM spec 0.1.2-draft 20090820 specified SM as a required
	//            tag.  But SAM spec 1.4 shows it as an optional field.

	//smTag = find_tabbed_tag (readGroup, "SM");
	//if (smTag == NULL)
	//	{
	//	if (errorText != NULL) *errorText = "SM is a required field";
	//	return NULL;
	//	}
	//smLen = tabbed_tag_length (smTag);
	//if (smLen <= 3)
	//	{
	//	if (errorText != NULL) *errorText = "SM field cannot be empty";
	//	return NULL;
	//	}

	// nota bene: SAM spec 1.4 appears to indicate that LB and PU tags don't
	//            have to be propagated from the RG header line to alignment
	//            records

	//lbTag = find_tabbed_tag (readGroup, "LB");
	//lbLen = 0; // (placate complier)
	//if (lbTag != NULL)
	//	{
	//	lbLen = tabbed_tag_length (lbTag);
	//	if (lbLen <= 3)
	//		{
	//		if (errorText != NULL) *errorText = "LB field cannot be empty";
	//		return NULL;
	//		}
	//	}

	//puTag = find_tabbed_tag (readGroup, "PU");
	//puLen = 0; // (placate complier)
	//if (puTag != NULL)
	//	{
	//	puLen = tabbed_tag_length (puTag);
	//	if (puLen <= 3)
	//		{
	//		if (errorText != NULL) *errorText = "PU field cannot be empty";
	//		return NULL;
	//		}
	//	}

	//////////
	// collect the ID, LB and PU fields into a single string
	//////////

	bytesNeeded = idLen + 1;
	//if (lbTag != NULL) bytesNeeded += 1 + lbLen;
	//if (puTag != NULL) bytesNeeded += 1 + puLen;

	s = ss = malloc_or_die ("sam_rg_tags", bytesNeeded);
	strncpy (/*to*/ ss, /*from*/ idTag, idLen);
	ss += idLen;

	//if (lbTag != NULL)
	//	{
	//	*(ss++) = '\t';
	//	strncpy (/*to*/ ss, /*from*/ lbTag, lbLen);
	//	ss += lbLen;
	//	}

	//if (puTag != NULL)
	//	{
	//	*(ss++) = '\t';
	//	strncpy (/*to*/ ss, /*from*/ puTag, puLen);
	//	ss += puLen;
	//	}

	*ss = 0;

	if (errorText != NULL) *errorText = NULL;
	return s;
	}

//----------
//
// print_sam_job_header--
//	Print sam format job header.
//
//----------

static int headerPrinted = false;

void print_sam_job_header
   (FILE*		f,
	char*		readGroup)
	{
	fprintf (f, "@HD\tVN:1.0\tSO:unsorted\n");
	if (readGroup != NULL)
		fprintf (f, "@RG\t%s\n", readGroup);
	headerPrinted = false;
	}

//----------
//
// print_sam_header--
//	Print sam format query header.
//
//----------

void print_sam_header
   (FILE* f,
	seq*  seq1,
	arg_dont_complain(seq*  seq2))
	{
	seqpartition*	sp1 = &seq1->partition;
	partition*		p;
	u32				ix;
	char*			name1;

	if (headerPrinted) return;

	// if seq1 is not partitioned, just print the (single) header name

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		fprintf (f, "@SQ\tSN:%s\tLN:" unsposFmt "\n",
		            name1, seq1->trueLen);
		}

	// otherwise, seq1 is partitioned, so print the header name for each
	// partition

	else
		{
		p = sp1->p;

		for (ix=0 ; ix<sp1->len ; ix++)
			fprintf (f, "@SQ\tSN:%s\tLN:" unsposFmt "\n",
			            &sp1->pool[p[ix].header], p[ix].trueLen);
		}

	headerPrinted = true;
	}

//----------
//
// print_sam_align_list--
//	Print a list of gapped alignments in sam format.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	alignel*	alignList:	The list of alignments to print.
//	seq*		seq1:		One sequence.
//	seq*		seq2:		Another sequence.
//	int			softMasked:	true  => sequence ends should be soft masked
//							false => sequence ends should be hard masked
//	char*		rgTags:		Additional tags for RG (read group).  This is a
//							.. valid string per the format described in the SAM
//							.. spec.  This may be NULL.
//
// Returns:
//	(nothing)
//
//----------

void print_sam_align_list
   (FILE*		f,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2,
	int			softMasked,
	char*		rgTags)
	{
	alignel*	a;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		print_sam_align (f,
		                 seq1, a->beg1-1, a->end1,
		                 seq2, a->beg2-1, a->end2,
		                 a->script, a->s, softMasked, rgTags);
		}

	}

//----------
//
// print_sam_align--
//	Print a single gapped alignment in sam format.
//
//----------
//
// Arguments:
//	FILE*		f:			The file to print to.
//	seq*		seq1:		One sequence.
//	unspos		beg1, end1:	Range of positions in sequence 1 (origin 0).
//	seq*		seq2:		Another sequence.
//	unspos		beg2, end2:	Range of positions in sequence 2 (origin 0).
//	editscript*	script:		The script describing the path the alignment takes
//							.. in the DP matrix.
//	score		s:			The alignment's score.
//	int			softMasked:	true  => sequence ends should be soft masked
//							false => sequence ends should be hard masked
//	char*		rgTags:		Additional tags for RG (read group).  This is a
//							.. valid string per the format described in the SAM
//							.. spec.  This may be NULL.
//
// Returns:
//	(nothing)
//
//----------

void print_sam_align
   (FILE*			f,
	seq*			seq1,
	unspos			beg1,
	unspos			end1,
	seq*			seq2,
	unspos			beg2,
	unspos			end2,
	editscript*		script,
	arg_dont_complain(score s),
	int				softMasked,
	char*			rgTags)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	unspos			height, width, i, j, prevI, prevJ, run;
	u32				opIx;
	unspos			len2;
	char*			name1, *name2;
	unspos			offset1, offset2, start1, start2;
	unspos			startLoc1, startLoc2;
	unspos			seq2Len, seq2True;
	int				flag;
	char			maskCh;
	unspos			preMask, postMask, tmp;

	if (seq1->revCompFlags != rcf_forward)
		suicide ("attempt to print - strand or complement for sequence 1 in print_sam_align");

	beg1++; // (internally, we want origin 1, inclusive)
	beg2++;

	       height = end1 - beg1 + 1;
	len2 = width  = end2 - beg2 + 1;

	//////////
	// figure out position offsets and names
	//////////

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, beg1-1);
		name1     = &sp1->pool[part->header];
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		offset2   = 0;
		seq2Len   = seq2->len;
		seq2True  = seq2->trueLen;
		startLoc2 = seq2->startLoc;
		}
	else					// sequence 2 is partitioned
	 	{
		part = lookup_partition (seq2, beg2-1);
		name2     = &sp2->pool[part->header];
		offset2   = part->sepBefore + 1;
		seq2Len   = part->sepAfter - offset2;
		seq2True  = part->trueLen;
		startLoc2 = part->startLoc;
		}

	//////////
	// print sam line (field names indicate below are per sam spec)
	//////////

	start1 = beg1-1 - offset1 + startLoc1;

	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2 = beg2-1 - offset2 + startLoc2;
		end2   = start2-1 + len2;
		flag   = 0;
		}
	else
		{
		start2 = startLoc2 + offset2 + (seq2Len - beg2) - (len2-1);
		end2   = startLoc2 + offset2 + (seq2Len - beg2);
		flag   = BAM_FREVERSE;
		}

	// print qname, flag, rname, pos and mapq

	fprintf (f, "%s\t%d\t%s\t" unsposFmt "\t%d\t",
	            name2, flag, name1, start1, 255);

	// print cigar

	maskCh = (softMasked)? 'S' : 'H';

	preMask = postMask = 0;
	if (start2 > 1) preMask = start2 - 1;
	if (end2 < seq2True) postMask = seq2True - end2;
	if ((seq2->revCompFlags & rcf_rev) != 0)
		{ tmp = preMask;  preMask = postMask;  postMask = tmp; }

	if (preMask != 0) fprintf (f, unsposFmt "%c", preMask, maskCh);

	opIx = 0;
	for (i=j=0 ; (i< height)||(j<width) ; )
		{
		run = edit_script_run_of_subs (script, &opIx);
		if (run > 0)
			{
			fprintf (f, unsposFmt "M", run);
			i += run; j += run;
			}

		if ((i < height) || (j < width))
			{
			prevI = i;  prevJ = j;
			edit_script_indel_len (script, &opIx, &i, &j);
			if (i > prevI)
				fprintf (f, unsposFmt "D", i - prevI);
			if (j > prevJ)
				fprintf (f, unsposFmt "I", j - prevJ);
			}
		}

	if (postMask != 0) fprintf (f, unsposFmt "%c", postMask, maskCh);

	// print mrnm, mpos, and isize

	fprintf (f, "\t%s\t%d\t%d\t", "*", 0, 0);

	// print seq (data from sequence 2)

	print_query_bases (f, seq2, beg2-1, len2, softMasked);

	// print qual (if we have no qual data, we print "*")

	if (seq2->vq == NULL)
		fprintf (f, "\t%s", "*");
	else
		{
		fprintf (f, "\t");
		print_query_quals (f, seq2, beg2-1, len2, softMasked);
		}

	// print tags

	if (rgTags != NULL)
		fprintf (f, "\t%s", rgTags);

	fprintf (f, "\n");
	}

//----------
//
// print_sam_match--
//	Print an hsp in sam format.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to print to.
//	seq*	seq1:		One sequence.
//	unspos	pos1:		The position, in seq1, of first character in the
//						.. match (origin-0).
//	seq*	seq2:		Another sequence.
//	unspos	pos1:		The position, in seq2, of first character in the
//						.. match (origin-0).
//	unspos	length:		The number of nucleotides in the HSP.
//	score	s:			The HSP's score.
//	int		softMasked:	true  => sequence ends should be soft masked
//						false => sequence ends should be hard masked
//	char*	rgTags:		Additional tags for RG (read group).  This is a valid
//						.. string per the format described in the SAM spec.
//						.. This may be NULL.
//
// Returns:
//	(nothing)
//
//----------

void print_sam_match
   (FILE*			f,
	seq*			seq1,
	unspos			pos1,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	arg_dont_complain(score s),
	int				softMasked,
	char*			rgTags)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	char*			name1, *name2;
	unspos			offset1, offset2, start1, start2, end2;
	unspos			startLoc1, startLoc2;
	unspos			seq2Len, seq2True;
	int				flag;
	char			maskCh;
	unspos			preMask, postMask, tmp;

	if (seq1->revCompFlags != rcf_forward)
		suicide ("attempt to print - strand or complement for sequence 1 in print_sam_match");

	//////////
	// figure out position offsets and names
	//////////

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		}
	else					// sequence 1 is partitioned
	 	{
		part     = lookup_partition (seq1, pos1);
		name1    = &sp1->pool[part->header];
		offset1  = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		offset2   = 0;
		seq2Len   = seq2->len;
		seq2True  = seq2->trueLen;
		startLoc2 = seq2->startLoc;
		}
	else					// sequence 2 is partitioned
	 	{
		part      = lookup_partition (seq2, pos2);
		name2     = &sp2->pool[part->header];
		offset2   = part->sepBefore + 1;
		seq2Len   = part->sepAfter - offset2;
		seq2True  = part->trueLen;
		startLoc2 = part->startLoc;
		}

	//////////
	// print sam line (field names indicate below are per sam spec)
	//////////

	start1 = pos1 - offset1 + startLoc1;

	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2 = pos2 - offset2 + startLoc2;
		end2   = start2-1 + length;
		flag   = 0;
		}
	else
		{
		start2 = startLoc2 + offset2 + (seq2Len - pos2) - length;
		end2   = startLoc2 + offset2 + (seq2Len - pos2) - 1;
		flag   = BAM_FREVERSE;
		}

	// print qname, flag, rname, pos and mapq (field names per sam spec)

	fprintf (f, "%s\t%d\t%s\t" unsposFmt "\t%d\t",
	            name2, flag, name1, start1, 255);

	// print cigar

	maskCh = (softMasked)? 'S' : 'H';

	preMask = postMask = 0;
	if (start2 > 1) preMask = start2 - 1;
	if (end2 < seq2True) postMask = seq2True - end2;
	if ((seq2->revCompFlags & rcf_rev) != 0)
		{ tmp = preMask;  preMask = postMask;  postMask = tmp; }

	if (preMask  != 0) fprintf (f, unsposFmt "%c", preMask, maskCh);
	                   fprintf (f, unsposFmt "M",  length);
	if (postMask != 0) fprintf (f, unsposFmt "%c", postMask, maskCh);

	// print mrnm, mpos, and isize

	fprintf (f, "\t%s\t%d\t%d\t", "*", 0, 0);

	// print seq (data from sequence 2)

	print_query_bases (f, seq2, pos2, length, softMasked);

	// print qual (if we have no qual data, we print "*")

	if (seq2->vq == NULL)
		fprintf (f, "\t%s", "*");
	else
		{
		fprintf (f, "\t");
		print_query_quals (f, seq2, pos2, length, softMasked);
		}

	// print rgTags

	if (rgTags != NULL)
		fprintf (f, "\t%s", rgTags);

	fprintf (f, "\n");
	}

//----------
//
// print_query_bases--
//	Print the "seq" field for sam format.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to print to.
//	seq*	seq2:		The query sequence.
//	unspos	pos1:		The position, in seq2, of first character in the
//						.. match (origin-0).
//	unspos	length:		The number of nucleotides in the HSP.
//	int		softMasked:	true  => sequence ends should be soft masked
//						false => sequence ends should be hard masked
//
// Returns:
//	(nothing)
//
//----------

static void print_query_bases
   (FILE*			f,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	int				softMasked)
	{
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	u8*				s2 = seq2->v + pos2;
	u8*				ss2;
	unspos			offset2, start2, end2;
	unspos			startLoc2;
	unspos			seq2True;
	unspos			ix;

	//////////
	// figure out position offsets, etc.
	//////////

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		offset2   = 0;
		startLoc2 = seq2->startLoc;
		seq2True  = seq2->trueLen;
		}
	else					// sequence 2 is partitioned
	 	{
		part      = lookup_partition (seq2, pos2);
		offset2   = part->sepBefore + 1;
		startLoc2 = part->startLoc;
		seq2True  = part->trueLen;
		}

	//////////
	// print seq (data from sequence 2)
	//////////

	start2 = pos2 - offset2 + startLoc2;
	end2   = start2-1 + length;

	if ((softMasked) && (start2 > 1))
		{
		ss2 = seq2->v + pos2 - (start2-1);
		for (ix=0 ; ix<start2-1 ; ix++)
			fprintf (f, "%c", dna_tolower(ss2[ix]));
		}

	for (ix=0 ; ix<length ; ix++)
		fprintf (f, "%c", dna_toupper(s2[ix]));

	if ((softMasked) && (end2 < seq2True))
		{
		for (ix=length ; ix<seq2True - (start2-1) ; ix++)
			fprintf (f, "%c", dna_tolower(s2[ix]));
		}

	}

//----------
//
// print_query_quals--
//	Print the "qual" field for sam format.
//
//----------
//
// Arguments:
//	FILE*	f:			The file to print to.
//	seq*	seq2:		The query sequence.
//	unspos	pos1:		The position, in seq2, of first character in the
//						.. match (origin-0).
//	unspos	length:		The number of nucleotides in the HSP.
//	int		softMasked:	true  => sequence ends should be soft masked
//						false => sequence ends should be hard masked
//
// Returns:
//	(nothing)
//
//----------

static void print_query_quals
   (FILE*			f,
	seq*			seq2,
	unspos			pos2,
	unspos			length,
	int				softMasked)
	{
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	u8*				s2 = seq2->vq + pos2;
	u8*				ss2;
	unspos			offset2, start2, end2;
	unspos			startLoc2;
	unspos			seq2True;
	unspos			ix;

	//////////
	// figure out position offsets, etc.
	//////////

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		offset2   = 0;
		startLoc2 = seq2->startLoc;
		seq2True  = seq2->trueLen;
		}
	else					// sequence 2 is partitioned
	 	{
		part      = lookup_partition (seq2, pos2);
		offset2   = part->sepBefore + 1;
		startLoc2 = part->startLoc;
		seq2True  = part->trueLen;
		}

	//////////
	// print qual (qualities from sequence 2)
	//////////

	start2 = pos2 - offset2 + startLoc2;
	end2   = start2-1 + length;

	if ((softMasked) && (start2 > 1))
		{
		ss2 = seq2->vq + pos2 - (start2-1);
		for (ix=0 ; ix<start2-1 ; ix++)
			fprintf (f, "%c", ss2[ix]);
		}

	for (ix=0 ; ix<length ; ix++)
		fprintf (f, "%c", s2[ix]);

	if ((softMasked) && (end2 < seq2True))
		{
		for (ix=length ; ix<seq2True - (start2-1) ; ix++)
			fprintf (f, "%c", s2[ix]);
		}

	}

