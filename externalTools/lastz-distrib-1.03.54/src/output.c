//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: output.c
//
//----------
//
// output--
//	Interface between the main program and he various output formats.
//
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
#include "sequences.h"			// sequence stuff
#include "infer_scores.h"		// scoring inference stuff
#include "identity_dist.h"		// identity distribution stuff
#include "gfa.h"				// gfa alignment format stuff
#include "lav.h"				// lav alignment format stuff
#include "axt.h"				// axt alignment format stuff
#include "maf.h"				// maf alignment format stuff
#include "sam.h"				// sam alignment format stuff
#include "cigar.h"				// cigar alignment format stuff
#include "genpaf.h"				// genpaf alignment format stuff
#include "text_align.h"			// textual alignment format stuff
#include "align_diffs.h"		// alignment differences format stuff
#include "lastz.h"				// lastz program-wide stuff

#define  output_owner			// (make this the owner of its globals)
#include "output.h"				// interface to this module

extern char* programName;				// (from lastz.c)
extern char* programVersionMajor;		// (from lastz.c)
extern char* programVersionMinor;		// (from lastz.c)
extern char* programVersionSubMinor;	// (from lastz.c)
extern char* programRevisionDate;		// (from lastz.c)

extern control* currParams;				// (from lastz.c)

// debugging defines

//#define snoopGenpaf			// if this is defined, extra code is added to
								// .. track calls to print_genpaf_align() and
								// .. print_genpaf_align();  note that another
								// .. instance of this define is in genpaf.c

//----------
//
// private global data
//
//----------

u32 printedForQuery;		// the number of alignments that have been printed
							// .. for the current query (both strands together)

int strandHeaderPrinted;	// false => we have yet to print a header for the
							//          .. current strand-to-strand alignment

// how often shall we flush the output?

#define matchFlushFrequency 1000

//----------
//
// prototypes for private functions
//
//----------

static void  print_match_composition (FILE* f,
                                      seq* seq1, unspos pos1,
                                      seq* seq2, unspos pos2,
                                      unspos length,
                                      score s,
                                      seed* hitSeed, u32 step);
static void  dump_match              (FILE* f,
                                      seq* seq1, unspos pos1,
                                      seq* seq2, unspos pos2,
                                      unspos length);
static char* program_name            (void);

//----------
//
// init_output_for_query, init_output_for_strand--
//
//----------

void init_output_for_query (void)
	{ printedForQuery = 0; }

void init_output_for_strand (void)
	{ strandHeaderPrinted = false; }

//----------
//
// print_align_list_segments--
//
//----------

//=== stuff for snoopGenpaf ===

#ifndef snoopGenpaf
#define snoopGenpaf_1 ;
#endif // not snoopGenpaf

#ifdef snoopGenpaf

#define snoopGenpaf_1                                                         \
	fprintf (stderr, "segmenting("                                            \
	                 "seq1:" unsposDotsFmt " seq2:" unsposDotsFmt ")\n",      \
		             beg1, end1, beg2, end2);

#endif // snoopGenpaf


//=== print_align_list_segments ===

void print_align_list_segments (alignel* alignList)
	{
	alignel*	a;
	unspos		beg1, end1, beg2, end2;
	unspos		height, width, i, j, prevI, prevJ, run;
	u32			opIx;
	score		s;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		beg1   = a->beg1;
		end1   = a->end1;
		beg2   = a->beg2;
		end2   = a->end2;
		height = end1 - beg1 + 1;
		width  = end2 - beg2 + 1;

		snoopGenpaf_1;

		// print the alignment's segments

		opIx = 0;
		for (i=j=0 ; (i< height)||(j<width) ; )
			{
			prevI = i;  prevJ = j;
			run = edit_script_run_of_subs (a->script, &opIx);
			i += run; j += run;
			if ((i < height) || (j < width))
				edit_script_indel_len (a->script, &opIx, &i, &j);

			s = score_match (currParams->scoring,
			                 currParams->seq1, beg1-1+prevI,
			                 currParams->seq2, beg2-1+prevJ,
			                 run);
			print_match (beg1-1+prevI, beg2-1+prevJ, run, s);
			}
		}

	}

//----------
//
// print_job_header, print_job_footer, print_header, print_align_list,
// print_match, print_comment_open, print_comment_close, print_eof_comment--
//
//----------

void print_job_header (void)
	{
	int outputFormat = currParams->outputFormat;

	switch (outputFormat)
		{
		case fmtGfa:
		case fmtGfaNoScore:
			print_gfa_job_header
			   (currParams->outputFile, program_name(),
			    currParams->seq1->filename, currParams->seq2->filename);
			break;
		case fmtLav:
		case fmtLavComment:
		case fmtLavScore:
		case fmtLavText:
		case fmtLavInfScores:
			print_lav_job_header
			   (currParams->outputFile, program_name(),
			    currParams->seq1->filename, currParams->seq2->filename, currParams->args,
		        currParams->scoring, &currParams->hspThreshold, &currParams->gappedThreshold,
		        currParams->dynamicMasking,
		        /*withExtras*/ (outputFormat==fmtLavComment),
		        currParams->xDrop, currParams->yDrop);
			if (outputFormat == fmtLavText)
				goto text_format;
			if (outputFormat == fmtLavInfScores)
				goto inf_scores_format;
			break;
		case fmtAxt:
		case fmtAxtComment:
		case fmtAxtGeneral:
			print_axt_job_header
			   (currParams->outputFile, program_name(), currParams->args,
		        currParams->scoring,
		        &currParams->hspThreshold, &currParams->gappedThreshold,
		        currParams->xDrop, currParams->yDrop);
			break;
		case fmtMaf:
		case fmtMafComment:
		case fmtMafNoComment:
			print_maf_job_header
			   (currParams->outputFile, program_name(), currParams->args,
		        currParams->scoring,
		        &currParams->hspThreshold, &currParams->gappedThreshold,
		        currParams->xDrop, currParams->yDrop,
		        (outputFormat != fmtMafNoComment));
			break;
		case fmtSoftSam:
		case fmtHardSam:
			print_sam_job_header (currParams->outputFile,currParams->readGroup);
			break;
		case fmtSoftSamNoHeader:
		case fmtHardSamNoHeader:
			; // (do nothing)
			break;
		case fmtCigar:
			print_cigar_job_header (currParams->outputFile);
			break;
		case fmtGenpaf:
			print_genpaf_job_header (currParams->outputFile, currParams->outputInfo);
			break;
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
			print_genpaf_job_header (NULL, NULL);
			break;
		case fmtGenpafBlast:
			print_blast_job_header (currParams->outputFile);
			break;
		case fmtGenpafBlastNoHeader:
			; // (do nothing)
			break;
		case fmtText:
		case fmtZeroText:
		text_format:
			print_text_align_job_header
			   (currParams->outputFile, program_name(),
			    currParams->seq1->filename, currParams->seq2->filename,
			    (outputFormat!=fmtZeroText));
			break;
		case fmtDiffs:
		case fmtDiffsNoBlocks:
			print_align_diffs_job_header
			   (currParams->outputFile, program_name(),
			    currParams->seq1->filename, currParams->seq2->filename);
			break;
		case fmtInfStats:
			init_inference_stats_job (currParams->seq1, currParams->seq2);
			break;
		case fmtIdDist:
			init_identity_dist_job (currParams->seq1, currParams->seq2);
			break;
		case fmtInfScores:
		inf_scores_format:
			; // (do nothing)
			break;
		case fmtHspComp:
		case fmtDeseed:
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_job_header, outputFormat=%d", outputFormat);
		}

//	if (currParams->dotplotFile != NULL)
//		; // (do nothing)
	}

void print_job_footer (void)
	{
	int outputFormat = currParams->outputFormat;

	switch (outputFormat)
		{
		case fmtGfa:
		case fmtGfaNoScore:
			print_gfa_job_footer (currParams->outputFile);
			break;
		case fmtLav:
		case fmtLavComment:
		case fmtLavScore:
		case fmtLavText:
		case fmtLavInfScores:
			print_lav_job_footer (currParams->outputFile);
			if (outputFormat == fmtLavText)
				goto text_format;
			if (outputFormat == fmtLavInfScores)
				goto inf_scores_format;
			break;
		case fmtAxt:
		case fmtAxtComment:
		case fmtAxtGeneral:
			print_axt_job_footer (currParams->outputFile);
			break;
		case fmtMaf:
		case fmtMafComment:
		case fmtMafNoComment:
			print_maf_job_footer (currParams->outputFile);
			break;
		case fmtSoftSam:
		case fmtSoftSamNoHeader:
		case fmtHardSam:
		case fmtHardSamNoHeader:
			; // (do nothing)
			break;
		case fmtCigar:
			print_cigar_job_footer (currParams->outputFile);
			break;
		case fmtGenpaf:
			print_genpaf_job_footer (currParams->outputFile);
			break;
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
			; // (do nothing)
			break;
		case fmtGenpafBlast:
			print_blast_job_footer (currParams->outputFile);
			break;
		case fmtGenpafBlastNoHeader:
			; // (do nothing)
			break;
		case fmtText:
		case fmtZeroText:
		text_format:
			print_text_align_job_footer (currParams->outputFile);
			break;
		case fmtDiffs:
		case fmtDiffsNoBlocks:
			print_align_diffs_job_footer (currParams->outputFile);
			break;
		case fmtInfStats:
			print_inference_stats_job (currParams->outputFile);
			break;
		case fmtIdDist:
			print_identity_dist_job (currParams->outputFile);
			break;
		case fmtInfScores:
		inf_scores_format:
			; // (do nothing)
			break;
		case fmtHspComp:
		case fmtDeseed:
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_job_footer, outputFormat=%d", outputFormat);
		}

//	if (currParams->dotplotFile != NULL)
//		; // (do nothing)
	}

void print_header (void)
	{
	static char*	prevName1 = NULL;
	static char*	prevName2 = NULL;
	static char		prevNameBuff1[maxSequenceName+1];
	static char		prevNameBuff2[maxSequenceName+1];
	char*			name1, *name2;
	int				outputFormat = currParams->outputFormat;

	if (prevName1 == NULL)
		{ prevName1 = prevNameBuff1;  prevNameBuff1[0] = 0; }
	if (prevName2 == NULL)
		{ prevName2 = prevNameBuff2;  prevNameBuff2[0] = 0; }

	switch (outputFormat)
		{
		case fmtGfa:
		case fmtGfaNoScore:
			print_gfa_header (currParams->outputFile,
			                  currParams->seq1, currParams->seq2);
			break;
		case fmtLav:
		case fmtLavComment:
		case fmtLavScore:
		case fmtLavText:
		case fmtLavInfScores:
			print_lav_header (currParams->outputFile,
			                  currParams->seq1, currParams->seq2);
			if (outputFormat == fmtLavText)
				goto text_format;
			if (outputFormat == fmtLavInfScores)
				goto inf_scores_format;
			break;
		case fmtAxt:
		case fmtAxtComment:
		case fmtAxtGeneral:
			print_axt_header (currParams->outputFile,
			                  currParams->seq1, currParams->seq2);
			break;
		case fmtMaf:
		case fmtMafComment:
		case fmtMafNoComment:
			print_maf_header (currParams->outputFile,
			                  currParams->seq1, currParams->seq2);
			break;
		case fmtSoftSam:
		case fmtHardSam:
			print_sam_header (currParams->outputFile,
			                  currParams->seq1, currParams->seq2);
			break;
		case fmtSoftSamNoHeader:
		case fmtHardSamNoHeader:
			; // (do nothing)
			break;
		case fmtCigar:
			print_cigar_header (currParams->outputFile,
			                    currParams->seq1, currParams->seq2);
			break;
		case fmtGenpaf:
			print_genpaf_header (currParams->outputFile,
			                     currParams->seq1, currParams->seq2);
			break;
		case fmtGenpafNoHeader:
			; // (do nothing)
			break;
		case fmtGenpafNameHeader:
			{
			name1 = name2 = NULL;
			if (currParams->seq1->partition.p == NULL) // sequence 1 is not partitioned
				name1 = (currParams->seq1->useFullNames)? currParams->seq1->header
														: currParams->seq1->shortHeader;
			if (currParams->seq2->partition.p == NULL) // sequence 2 is not partitioned
				name2 = (currParams->seq1->useFullNames)? currParams->seq2->header
														: currParams->seq2->shortHeader;
			if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
			if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
			if ((strcmp (name1, prevName1) != 0)
			 || (strcmp (name2, prevName2) != 0))
				{
				if (strcmp(currParams->outputInfo,genpafRDotplotScoreKeys) == 0)
					fprintf (currParams->outputFile, "%s\t%s\tscore\n", name1, name2);
				else
					fprintf (currParams->outputFile, "%s\t%s\n", name1, name2);
				strncpy (/*to*/ prevName1, /*from*/ name1, sizeof(prevNameBuff1));
				strncpy (/*to*/ prevName2, /*from*/ name2, sizeof(prevNameBuff2));
				}
			}
			break;
		case fmtGenpafBlast:
			print_blast_header
			   (currParams->outputFile, program_name(), currParams->args,
			    currParams->seq1, currParams->seq2);
			break;
		case fmtGenpafBlastNoHeader:
			; // (do nothing)
			break;
		case fmtText:
		case fmtZeroText:
		text_format:
			print_text_align_header (currParams->outputFile,
			                         currParams->seq1, currParams->seq2,
			                         (outputFormat!=fmtZeroText));
			break;
		case fmtDiffs:
		case fmtDiffsNoBlocks:
			print_align_diffs_header (currParams->outputFile,
			                          currParams->seq1, currParams->seq2);
			break;
		case fmtHspComp:
		case fmtInfStats:
		case fmtInfScores:
		inf_scores_format:
			; // (do nothing)
			break;
		case fmtIdDist:
		case fmtDeseed:
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_header, outputFormat=%d", outputFormat);
		}

	if (currParams->dotplotFile != NULL)
		{
		name1 = name2 = NULL;
		if (currParams->seq1->partition.p == NULL) // sequence 1 is not partitioned
			name1 = (currParams->seq1->useFullNames)? currParams->seq1->header
													: currParams->seq1->shortHeader;
		if (currParams->seq2->partition.p == NULL) // sequence 2 is not partitioned
			name2 = (currParams->seq1->useFullNames)? currParams->seq2->header
													: currParams->seq2->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		if ((strcmp (name1, prevName1) != 0)
		 || (strcmp (name2, prevName2) != 0))
			{
			if (strcmp(currParams->dotplotKeys,genpafRDotplotScoreKeys) == 0)
				fprintf (currParams->dotplotFile, "%s\t%s\tscore\n", name1, name2);
			else
				fprintf (currParams->dotplotFile, "%s\t%s\n", name1, name2);
			strncpy (/*to*/ prevName1, /*from*/ name1, sizeof(prevNameBuff1));
			strncpy (/*to*/ prevName2, /*from*/ name2, sizeof(prevNameBuff2));
			}
		}
	}

void print_align_list (alignel* alignList)
	{
	int outputFormat = currParams->outputFormat;
	alignel* a;

	if ((currParams->searchLimit > 0)
	 && (printedForQuery >= currParams->searchLimit))
		return;
	printedForQuery++;

	if (!strandHeaderPrinted)
		{ print_header ();  strandHeaderPrinted = true; }

	if (infer_scores_dbgShowIdentity)
		{
		unspos numer, denom;
		u32    bin;

		for (a=alignList ; a!=NULL ; a=a->next)
			{
			alignment_identity (currParams->seq1, currParams->seq2, a,
			                    &numer, &denom);
			bin = identity_bin (numer, denom);
			// nota bene: positions written as 1-based
			print_generic (currParams->outputFile,
			               unsposSlashFmt
			               " pct_identity=" unsposSlashFmt
			               " (bin as " identityBinFormat ")",
			               a->beg1, a->beg2,
			               numer, denom,
			               bin_to_identity (bin));
			}
		}

	switch (outputFormat)
		{
		case fmtGfa:
		case fmtGfaNoScore:
			print_gfa_align_list (currParams->outputFile,
			                      (outputFormat == fmtGfa)? currParams->scoring
			                                              : NULL,
			                      alignList,
								  currParams->seq1, currParams->seq2);
			break;
		case fmtLav:
		case fmtLavComment:
		case fmtLavScore:
		case fmtLavInfScores:
			print_lav_align_list (currParams->outputFile,
			                      alignList,
			                      currParams->seq1, currParams->seq2);
			if (outputFormat == fmtLavInfScores)
				goto inf_scores_format;
			break;
		case fmtLavText:
			for (a=alignList ; a!=NULL ; a=a->next)
				{
				print_lav_align        (currParams->outputFile,
										a->seq1, a->beg1-1, a->end1,
										a->seq2, a->beg2-1, a->end2,
										a->script, a->s);
				print_text_align_align (currParams->outputFile,
										currParams->seq1, a->beg1-1, a->end1,
										currParams->seq2, a->beg2-1, a->end2,
										a->script, a->s,
										false, currParams->textContext);
				}
			break;
		case fmtAxt:
		case fmtAxtComment:
			print_axt_align_list (currParams->outputFile, alignList,
			                      currParams->seq1, currParams->seq2,
			                      /* comments */ outputFormat==fmtAxtComment,
			                      /* extras   */ NULL);
			break;
		case fmtAxtGeneral:
			print_axt_align_list (currParams->outputFile, alignList,
			                      currParams->seq1, currParams->seq2,
			                      /* comments */ false,
			                      /* extras   */ currParams->outputInfo);
			break;
		case fmtMaf:
		case fmtMafNoComment:
			print_maf_align_list (currParams->outputFile,
			                      alignList, currParams->seq1, currParams->seq2,
			                      /* comments */ false);
			break;
		case fmtMafComment:
			print_maf_align_list (currParams->outputFile,
			                      alignList, currParams->seq1, currParams->seq2,
			                      /* comments */ true);
			break;
		case fmtSoftSam:
		case fmtSoftSamNoHeader:
			print_sam_align_list (currParams->outputFile,
			                      alignList, currParams->seq1, currParams->seq2,
			                      /* softMasking */ true,
			                      currParams->samRGTags);
			break;
		case fmtHardSam:
		case fmtHardSamNoHeader:
			print_sam_align_list (currParams->outputFile,
			                      alignList, currParams->seq1, currParams->seq2,
			                      /* softMasking */ false,
			                      currParams->samRGTags);
			break;
		case fmtCigar:
			print_cigar_align_list (currParams->outputFile,
			                        alignList, currParams->seq1, currParams->seq2,
			                        /* withInfo       */ true,
			                        /* markMismatches */ false,
			                        /* letterAfter    */ false,
			                        /* hideSingles    */ false,
			                        /* lowerCase      */ false,
			                        /* withNewLine    */ true);
			break;
		case fmtGenpaf:
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
		case fmtGenpafBlast:
		case fmtGenpafBlastNoHeader:
			print_genpaf_align_list (currParams->outputFile,
			                         alignList, currParams->seq1, currParams->seq2,
			                         currParams->outputInfo);
			break;
		case fmtText:
		case fmtZeroText:
			print_text_align_align_list (currParams->outputFile,
			                             alignList, currParams->seq1, currParams->seq2,
			                             (outputFormat!=fmtZeroText),
			                             currParams->textContext);
			break;
		case fmtDiffs:
		case fmtDiffsNoBlocks:
			print_align_diffs_align_list (currParams->outputFile,
			                              alignList, currParams->seq1, currParams->seq2,
			                              (outputFormat == fmtDiffs),
			                              currParams->nIsAmbiguous);
			break;
		case fmtInfStats:
			infer_stats_from_align_list (alignList, currParams->seq1, currParams->seq2);
			break;
		case fmtInfScores:
		inf_scores_format:
			gather_stats_from_align_list (alignList, currParams->seq1, currParams->seq2);
			break;
		case fmtIdDist:
			identity_dist_from_align_list (alignList, currParams->seq1, currParams->seq2);
			break;
		case fmtHspComp:
		case fmtDeseed:
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_align_list, outputFormat=%d", outputFormat);
		}

	if (currParams->dotplotFile != NULL)
		print_genpaf_align_list_segments (currParams->dotplotFile,
		                                  alignList, currParams->seq1, currParams->seq2,
		                                  currParams->dotplotKeys,
		                                  currParams->scoring);
	}

void print_match (unspos pos1, unspos pos2, unspos length, score s)
	// pos1 and pos2 are the positions of first character in the match,
	// .. (origin-0).
	{
	static u32 printsUntilFlush = matchFlushFrequency;
	int outputFormat = currParams->outputFormat;

	if ((currParams->searchLimit > 0)
	 && (printedForQuery >= currParams->searchLimit))
		return;
	printedForQuery++;

	if (!strandHeaderPrinted)
		{ print_header ();  strandHeaderPrinted = true; }

	if (infer_scores_dbgShowIdentity)
		{
		unspos numer, denom;
		u32    bin;

		segment_identity (currParams->seq1, pos1, currParams->seq2, pos2, length,
		                  &numer, &denom);
		bin = identity_bin (numer, denom);
		// nota bene: positions written as 1-based
		print_generic (currParams->outputFile,
		               unsposSlashFmt
		               " pct_identity=" unsposSlashFmt
		               " (bin as " identityBinFormat ")",
		               pos1+1, pos2+1,
		               numer, denom,
		               bin_to_identity (bin));
		}

	switch (outputFormat)
		{
		case fmtGfa:
		case fmtGfaNoScore:
			print_gfa_match (currParams->outputFile,
			                 currParams->seq1, pos1,
			                 currParams->seq2, pos2, length,
		                     (outputFormat == fmtGfa)? s : 0);
			break;
		case fmtLav:
		case fmtLavComment:
		case fmtLavText:
		case fmtLavInfScores:
			print_lav_match (currParams->outputFile,
			                 currParams->seq1, pos1,
			                 currParams->seq2, pos2, length,
			                 s);
			if (outputFormat == fmtLavText)
				goto text_format;
			if (outputFormat == fmtLavInfScores)
				goto inf_scores_format;
			break;
		case fmtLavScore:
			print_lavscore_match (currParams->outputFile,
			                      currParams->seq1, pos1,
			                      currParams->seq2, pos2, length,
			                      s);
			break;
		case fmtAxt:
		case fmtAxtComment:
			print_axt_match (currParams->outputFile,
			                 currParams->seq1, pos1,
			                 currParams->seq2, pos2, length,
			                 s,
			                 /* comments */ outputFormat==fmtAxtComment,
			                 /* extras   */ NULL);
			break;
		case fmtAxtGeneral:
			print_axt_match (currParams->outputFile,
			                 currParams->seq1, pos1,
			                 currParams->seq2, pos2, length,
			                 s,
			                 /* comments */ false,
			                 /* extras   */ currParams->outputInfo);
			break;
		case fmtMaf:
		case fmtMafNoComment:
			print_maf_match (currParams->outputFile,
			                 currParams->seq1, pos1,
			                 currParams->seq2, pos2, length,
			                 s, /* comments */ false);
			break;
		case fmtMafComment:
			print_maf_match (currParams->outputFile,
			                 currParams->seq1, pos1,
			                 currParams->seq2, pos2, length,
			                 s, /* comments */ true);
			break;
		case fmtSoftSam:
		case fmtSoftSamNoHeader:
			print_sam_match (currParams->outputFile,
			                 currParams->seq1, pos1,
			                 currParams->seq2, pos2, length,
			                 s,
			                 /* softMasking */ true,
			                 currParams->samRGTags);
			break;
		case fmtHardSam:
		case fmtHardSamNoHeader:
			print_sam_match (currParams->outputFile,
			                 currParams->seq1, pos1,
			                 currParams->seq2, pos2, length,
			                 s,
			                 /* softMasking */ false,
			                 currParams->samRGTags);
			break;
		case fmtCigar:
			print_cigar_match (currParams->outputFile,
			                   currParams->seq1, pos1,
			                   currParams->seq2, pos2, length,
			                   s,
			                   /* withInfo       */ true,
			                   /* markMismatches */ false,
			                   /* letterAfter    */ false,
			                   /* hideSingles    */ false,
			                   /* lowerCase      */ false,
			                   /* withNewLine    */ true);
			break;
		case fmtGenpaf:
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
		case fmtGenpafBlast:
		case fmtGenpafBlastNoHeader:
			print_genpaf_match (currParams->outputFile,
			                    currParams->seq1, pos1,
			                    currParams->seq2, pos2, length,
			                    s, currParams->outputInfo);
			break;
		case fmtText:
		case fmtZeroText:
		text_format:
			print_text_align_match (currParams->outputFile,
			                        currParams->seq1, pos1,
			                        currParams->seq2, pos2, length,
			                        s,
			                        (outputFormat!=fmtZeroText),
			                        currParams->textContext);
			break;
		case fmtDiffs:
		case fmtDiffsNoBlocks:
			print_align_diffs_match (currParams->outputFile,
			                         currParams->seq1, pos1,
			                         currParams->seq2, pos2, length,
			                         (outputFormat == fmtDiffs),
			                         currParams->nIsAmbiguous);
			break;
		case fmtHspComp:
			print_match_composition (currParams->outputFile,
			                         currParams->seq1, pos1,
			                         currParams->seq2, pos2, length,
			                         s, currParams->hitSeed, currParams->step);
			break;
		case fmtInfStats:
			infer_stats_from_match (currParams->seq1, pos1,
			                        currParams->seq2, pos2, length);
			break;
		case fmtInfScores:
		inf_scores_format:
			gather_stats_from_match (currParams->seq1, pos1,
			                         currParams->seq2, pos2, length);
			break;
		case fmtIdDist:
			identity_dist_from_match (currParams->seq1, pos1,
			                          currParams->seq2, pos2, length);
			break;
		case fmtDeseed:
			dump_match (currParams->outputFile,
			            currParams->seq1, pos1,
			            currParams->seq2, pos2, length);
			printf ("\n");
			break;
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_match, outputFormat=%d", outputFormat);
		}

	if (currParams->dotplotFile != NULL)
		print_genpaf_match (currParams->dotplotFile,
		                    currParams->seq1, pos1,
		                    currParams->seq2, pos2, length,
		                    s, currParams->dotplotKeys);

	if (--printsUntilFlush == 0)
		{
		fflush (currParams->outputFile);
		printsUntilFlush = matchFlushFrequency;
		}
	}


char* print_comment_open (void)
	{
	int   outputFormat = currParams->outputFormat;
	char* commentPrefix = NULL;

	switch (outputFormat)
		{
		case fmtLav:
		case fmtLavComment:
		case fmtLavScore:
		case fmtLavText:
		case fmtLavInfScores:
			print_lav_comment_open (currParams->outputFile);
			break;
		case fmtGfa:
		case fmtGfaNoScore:
			commentPrefix = "#";
			break;
		case fmtAxt:
		case fmtAxtComment:
		case fmtAxtGeneral:
			commentPrefix = "#";
			break;
		case fmtMaf:
		case fmtMafComment:
		case fmtMafNoComment:
			fprintf (stderr, "WARNING.  Output is not properly MAF format\n");
			commentPrefix = "#";
			break;
		case fmtSoftSam:
		case fmtSoftSamNoHeader:
		case fmtHardSam:
		case fmtHardSamNoHeader:
			fprintf (stderr, "WARNING.  Output is not properly SAM format\n");
			commentPrefix = "#";
			break;
		case fmtCigar:
			fprintf (stderr, "WARNING.  Output is not properly CIGAR format\n");
			commentPrefix = "#";
			break;
		case fmtGenpaf:
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
		case fmtGenpafBlast:
		case fmtGenpafBlastNoHeader:
			commentPrefix = "#";
			break;
		case fmtText:
		case fmtZeroText:
			; // (do nothing)
			break;
		case fmtHspComp:
			commentPrefix = "#";
			break;
		case fmtDiffs:
		case fmtDiffsNoBlocks:
			; // (do nothing)
			break;
		case fmtInfStats:
			commentPrefix = "#";
			break;
		case fmtInfScores:
			commentPrefix = "#";
			break;
		case fmtIdDist:
			commentPrefix = "#";
			break;
		case fmtDeseed:
			commentPrefix = "#";
			break;
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_comment_open, outputFormat=%d", outputFormat);
		}

	return commentPrefix;
	}


void print_comment_close (void)
	{
	int outputFormat = currParams->outputFormat;

	switch (outputFormat)
		{
		case fmtLav:
		case fmtLavComment:
		case fmtLavScore:
		case fmtLavText:
		case fmtLavInfScores:
			print_lav_comment_close (currParams->outputFile);
			break;
		case fmtGfa:
		case fmtGfaNoScore:
		case fmtAxt:
		case fmtAxtComment:
		case fmtAxtGeneral:
		case fmtMaf:
		case fmtMafComment:
		case fmtMafNoComment:
		case fmtSoftSam:
		case fmtSoftSamNoHeader:
		case fmtHardSam:
		case fmtHardSamNoHeader:
		case fmtCigar:
		case fmtGenpaf:
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
		case fmtGenpafBlast:
		case fmtGenpafBlastNoHeader:
		case fmtText:
		case fmtZeroText:
		case fmtHspComp:
		case fmtDiffs:
		case fmtDiffsNoBlocks:
		case fmtInfStats:
		case fmtInfScores:
		case fmtIdDist:
		case fmtDeseed:
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_comment_close, outputFormat=%d", outputFormat);
		}
	}


void print_eof_comment (void)
	{
	int outputFormat = currParams->outputFormat;

	if (outputFormat != fmtNone)
		fprintf (currParams->outputFile, "# lastz end-of-file\n");
	}


void print_m_stanza (census* cen)
	{ // note that census might be NULL
	int outputFormat = currParams->outputFormat;

	switch (outputFormat)
		{
		case fmtLav:
		case fmtLavComment:
		case fmtLavScore:
		case fmtLavText:
		case fmtLavInfScores:
			print_lav_m_stanza (currParams->outputFile, cen);
			break;
		case fmtGfa:
		case fmtGfaNoScore:
		case fmtAxt:
		case fmtAxtComment:
		case fmtAxtGeneral:
		case fmtMaf:
		case fmtMafComment:
		case fmtMafNoComment:
		case fmtSoftSam:
		case fmtSoftSamNoHeader:
		case fmtHardSam:
		case fmtHardSamNoHeader:
		case fmtCigar:
		case fmtGenpaf:
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
		case fmtGenpafBlast:
		case fmtGenpafBlastNoHeader:
		case fmtText:
		case fmtZeroText:
		case fmtHspComp:
		case fmtDiffs:
		case fmtDiffsNoBlocks:
		case fmtInfStats:
		case fmtInfScores:
		case fmtIdDist:
		case fmtDeseed:
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_m_stanza, outputFormat=%d", outputFormat);
		}

//	if (currParams->dotplotFile != NULL)
//		; // (do nothing)
	}

void print_census_stanza (census* cen)
	{
	int outputFormat = currParams->outputFormat;

	switch (outputFormat)
		{
		case fmtLav:
		case fmtLavComment:
		case fmtLavScore:
		case fmtLavText:
		case fmtLavInfScores:
			print_lav_census_stanza (currParams->outputFile, cen);
			break;
		case fmtGfa:
		case fmtGfaNoScore:
		case fmtAxt:
		case fmtAxtComment:
		case fmtAxtGeneral:
		case fmtMaf:
		case fmtMafComment:
		case fmtMafNoComment:
		case fmtSoftSam:
		case fmtSoftSamNoHeader:
		case fmtHardSam:
		case fmtHardSamNoHeader:
		case fmtCigar:
		case fmtGenpaf:
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
		case fmtGenpafBlast:
		case fmtGenpafBlastNoHeader:
		case fmtText:
		case fmtZeroText:
		case fmtHspComp:
		case fmtDiffs:
		case fmtDiffsNoBlocks:
		case fmtInfStats:
		case fmtInfScores:
		case fmtIdDist:
		case fmtDeseed:
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_census_stanza, outputFormat=%d", outputFormat);
		}

//	if (currParams->dotplotFile != NULL)
//		; // (do nothing)
	}

void print_x_stanza (unspos numMasked)
	{
	int outputFormat = currParams->outputFormat;

	switch (outputFormat)
		{
		case fmtGfa:
		case fmtGfaNoScore:
			print_gfa_generic (currParams->outputFile,
			                   'x', "num_masked=" unsposFmt, numMasked);
			break;
		case fmtLav:
		case fmtLavComment:
		case fmtLavScore:
		case fmtLavText:
		case fmtLavInfScores:
			print_lav_x_stanza (currParams->outputFile, numMasked);
			break;
		case fmtAxt:
		case fmtAxtComment:
		case fmtAxtGeneral:
		case fmtMaf:
		case fmtMafComment:
		case fmtMafNoComment:
		case fmtSoftSam:
		case fmtSoftSamNoHeader:
		case fmtHardSam:
		case fmtHardSamNoHeader:
		case fmtCigar:
		case fmtGenpaf:
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
		case fmtGenpafBlast:
		case fmtGenpafBlastNoHeader:
		case fmtText:
		case fmtZeroText:
		case fmtHspComp:
		case fmtDiffs:
		case fmtDiffsNoBlocks:
		case fmtInfStats:
		case fmtIdDist:
		case fmtDeseed:
			print_generic (currParams->outputFile,
			               "num_masked=" unsposFmt, numMasked);
			break;
		case fmtInfScores:
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_x_stanza, outputFormat=%d", outputFormat);
		}

//	if (currParams->dotplotFile != NULL)
//		; // (do nothing)
	}

void print_generic
   (FILE*		f,
	const char*	format,
	...)
	{
	int outputFormat;
	va_list	args;

	va_start (args, format);

	outputFormat = currParams->outputFormat;

	switch (outputFormat)
		{
		case fmtGfa:
		case fmtGfaNoScore:
			vprint_gfa_generic (f, 'z', format, args);
			break;
		case fmtLavComment:
			vprint_lav_comment (f, format, args);
			break;
		case fmtLavText:
			vprint_lav_comment (f, format, args);
			if (format != NULL)
				{
				va_end (args);
				va_start (args, format);
				vfprintf (f, format, args);
				fprintf  (f, "\n");
				}
			break;
		case fmtAxtComment:
			vprint_axt_comment (f, format, args);
			break;
		case fmtMafComment:
			vprint_maf_comment (f, format, args);
			break;
		case fmtText:
		case fmtZeroText:
			if (format != NULL)
				{
				vfprintf (f, format, args);
				fprintf  (f, "\n");
				}
			break;
		case fmtLav:
		case fmtLavScore:
		case fmtLavInfScores:
		case fmtAxt:
		case fmtAxtGeneral:
		case fmtMaf:
		case fmtMafNoComment:
		case fmtSoftSam:
		case fmtSoftSamNoHeader:
		case fmtHardSam:
		case fmtHardSamNoHeader:
		case fmtCigar:
		case fmtGenpaf:
		case fmtGenpafNoHeader:
		case fmtGenpafNameHeader:
		case fmtGenpafBlast:
		case fmtGenpafBlastNoHeader:
		case fmtHspComp:
		case fmtDiffs:
		case fmtDiffsNoBlocks:
		case fmtInfStats:
		case fmtInfScores:
		case fmtIdDist:
		case fmtDeseed:
		case fmtNone:
			; // (do nothing)
			break;
		default:
			suicidef ("internal error, in print_generic, outputFormat=%d", outputFormat);
		}

//	if (currParams->dotplotFile != NULL)
//		; // (do nothing)

	va_end (args);
	}

//----------
//
// print_match_composition--
//	Print a gap-free alignment including position and composition (counts of
//	matched dna letter pairs).
//
// Typical output is shown below, with a header added.  The first letter of the
// pairs is from sequence 1, the second from sequence 2.  P is the 'discovery
// probability'-- the probability that this HSP would be discovered for this
// (seed,Z) combination, over random sequence positions.
//
//	id score pos1/pos2  len p    AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
//	92 121  1475+/1395- 145 .750 27  0  4  0  1 38  2  0  3  0 45  0  1  1  0 23
//	88  28  3374+/4837-  42 .200 18  0  0  0  1  7  0  0  2  1  7  1  0  0  0  5
//	 ...
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The first aligned position in sequence 1.
//	seq*	seq2:	The other sequence.
//	unspos	pos2:	The first aligned position in sequence 2.
//	unspos	length:	The length of the alignment.
//	seed*	hitSeed: Seeding strategy for the hits that found this match.
//	u32		step:	Positional step size in the search for those hits.
//
// Returns:
//	(nothing)
//
//----------

static void print_match_composition
   (FILE*	f,
	seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length,
	score	s,
	seed*	hitSeed,
	u32		step)
	{
	int		pctId;
	unspos	count[4][4];
	float	p;
	char	pstr[6];
	int		ix, iy;

	// compute percent identity, match compostion, and discovery probability

	pctId = percent_identical (seq1, pos1, seq2, pos2, length);
	match_composition (seq1, pos1, seq2, pos2, length, count);
	p = discovery_probability (seq1, pos1+length, seq2, pos2+length, length,
	                           hitSeed, step);

	// convert discovery probability to a string

	if      (p < 0.0) p = 0.0;
	else if (p > 1.0) p = 1.0;

	snprintf (pstr, sizeof(pstr), "%.3f", p);
	if (pstr[0] == '1') // (1.000 -> 1.00)
		pstr[4] = 0;
	else				// (0.XXX -> .XXX)
		{
		pstr[0] = pstr[1];
		pstr[1] = pstr[2];
		pstr[2] = pstr[3];
		pstr[3] = pstr[4];
		pstr[4] = 0;
		}

	// print it

	fprintf (f, "%d " scoreFmtSimple " " unsposSlashSFmt " " unsposFmt " %s",
	            pctId, s,
	            pos1+1, ((seq1->revCompFlags & rcf_rev) != 0)? "-" : "+",
	            pos2+1, ((seq2->revCompFlags & rcf_rev) != 0)? "-" : "+",
	            length, pstr);

	for (ix=0 ; ix<4 ; ix++)
		for (iy=0 ; iy<4 ; iy++)
			fprintf (f, " " unsposFmt, count[ix][iy]);

	fprintf (f, "\n");
	}

//----------
//
// dump_match--
//	Dump the nucleotides (from each sequence) for a gap-free alignment.
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	seq*	seq1:	One sequence.
//	unspos	pos1:	The first aligned position in sequence 1.
//	seq*	seq2:	The other sequence.
//	unspos	pos2:	The first aligned position in sequence 2.
//	unspos	length:	The length of the alignment.
//
// Returns:
//	(nothing)
//
//----------

static void dump_match
   (FILE*	f,
	seq*	seq1,
	unspos	pos1,
	seq*	seq2,
	unspos	pos2,
	unspos	length)
	{
	print_prefix (f, (char*) seq1->v + pos1, length);
	fprintf      (f, "\n");
	print_prefix (f, (char*) seq2->v + pos2, length);
	fprintf      (f, "\n");
	}

//----------
//
// program_name--
//	Determnine the name of this program.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	A string describing the name of this program.  This may point to static
//	memory belonging to this routine, or it may point to memory in global
//	memory space.  But in any case, the caller should *not* deallocate the
//	returned pointer.
//
//----------

static char* program_name
   (void)
	{
	static char	_progName[101];
	int			n;

	n = snprintf (NULL, 0, "%s.v%s.%s.%s",
	              programName, programVersionMajor, programVersionMinor, programVersionSubMinor);
	if (((unsigned) n) < sizeof(_progName))
		{
		sprintf (_progName,
		         "%s.v%s.%s.%s",
		         programName, programVersionMajor, programVersionMinor, programVersionSubMinor);
		return _progName;
		}
	else
		return programName;
	}

