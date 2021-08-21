/*
 * stPaf.h
 *
 * Functions for manipulating paf alignments.
 */

#ifndef ST_PAF_H_
#define ST_PAF_H_

#include "sonLib.h"

/*
 * The structure of a paf alignment is as follows:
 *
 * PAF (from: https://lh3.github.io/minimap2/minimap2.html) and
 *  (https://github.com/lh3/miniasm/blob/master/PAF.md)
 * 1    string	Query sequence name
 * 2	int	Query sequence length
 * 3	int	Query start coordinate (0-based, inclusive)
 * 4	int	Query end coordinate (0-based, exclusive)
 * 5	char	‘+’ if query/target on the same strand; ‘-’ if opposite
 * 6	string	Target sequence name
 * 7	int	Target sequence length
 * 8	int	Target start coordinate on the original strand (0-based, inclusive)
 * 9	int	Target end coordinate on the original strand (0-based, inclusive)
 * 10	int	Number of matching bases in the mapping
 * 11	int	Number bases, including gaps, in the mapping
 * 12	int	Mapping quality (0-255 with 255 for missing)
 *
 * Tags: (only a subset are supported as indicated below)
 *
 * tp	A	Type of aln: P/primary, S/secondary and I,i/inversion
 * cm	i	Number of minimizers on the chain [NOT SUPPORTED/IGNORED]
 * s1	i	Chaining score [NOT SUPPORTED/IGNORED]
 * s2	i	Chaining score of the best secondary chain [NOT SUPPORTED/IGNORED]
 * NM	i	Total number of mismatches and gaps in the alignment [NOT SUPPORTED/IGNORED]
 * MD	Z	To generate the ref sequence in the alignment [NOT SUPPORTED/IGNORED]
 * AS	i	DP alignment score
 * SA	Z	List of other supplementary alignments [NOT SUPPORTED/IGNORED]
 * ms	i	DP score of the max scoring segment in the alignment [NOT SUPPORTED/IGNORED]
 * nn	i	Number of ambiguous bases in the alignment [NOT SUPPORTED/IGNORED]
 * ts	A	Transcript strand (splice mode only) [NOT SUPPORTED/IGNORED]
 * cg	Z	CIGAR string (only in PAF)
 * cs	Z	Difference string [NOT SUPPORTED/IGNORED]
 * dv	f	Approximate per-base sequence divergence [NOT SUPPORTED/IGNORED]
 * de	f	Gap-compressed per-base sequence divergence [NOT SUPPORTED/IGNORED]
 * rl	i	Length of query regions harboring repetitive seeds [NOT SUPPORTED/IGNORED]
 *
 */

typedef enum _cigarOp {
    match = 0,
    query_insert = 1,
    query_delete = 2
} CigarOp;

typedef struct _cigar Cigar;
struct _cigar {
    CigarOp op;
    int64_t length;
    Cigar *next;
};

typedef struct _paf {
    char *query_name;
    int64_t query_length;
    int64_t query_start; // Zero-based
    int64_t query_end; // Zero-based
    char *target_name;
    int64_t target_length;
    int64_t target_start; // Zero-based
    int64_t target_end; // Zero-based
    bool same_strand;
    Cigar *cigar;
    int64_t score; // the dp alignment score
    int64_t mapping_quality;
    int64_t num_matches;
    int64_t num_bases;
    char type; // is 'P' primary / 'S' secondary / 'I' inversion
} Paf;

/*
 * Cleanup the paf record.
 */
void paf_destruct(Paf *paf);

/*
 * Parse a paf from a string.
 */
Paf *paf_parse(char *paf_string);

/*
 * Read a PAF alignment record from the given file. Returns NULL if no record available.
 */
Paf *paf_read(FILE *fh);

/*
 * Prints a paf record
 */
char *paf_print(Paf *paf);

/*
 * Writes a PAF alignment to the given file.
 */
void paf_write(Paf *paf, FILE *fh);

/*
 * Checks a paf alignment coordinates and cigar are valid, error aborts if not.
 */
void paf_check(Paf *paf);

/*
 * Inverts the query and target in the paf record, including inverting any cigar string.
 */
void paf_invert(Paf *paf);

/*
 * Read all the pafs from a file in order.
 */
stList *read_pafs(FILE *paf_file);

/*
 * Write a list of pafs to a file in order.
 */
void write_pafs(FILE *paf_file, stList *pafs);

/*
 * Chain a set of pafs into larger alignments
 */
stList *paf_chain(stList *pafs, int64_t (*gap_cost)(int64_t, void *), void *gap_cost_params, int64_t max_gap_length);

#endif /* ST_PAF_H_ */

