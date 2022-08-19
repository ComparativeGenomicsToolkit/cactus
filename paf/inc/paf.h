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
 * cg	Z	CIGAR string (only in PAF - this is ordered by the target sequence)
 * cs	Z	Difference string [NOT SUPPORTED/IGNORED]
 * dv	f	Approximate per-base sequence divergence [NOT SUPPORTED/IGNORED]
 * de	f	Gap-compressed per-base sequence divergence [NOT SUPPORTED/IGNORED]
 * rl	i	Length of query regions harboring repetitive seeds [NOT SUPPORTED/IGNORED]
 * tl	i	Tile level of paf in the cactus chainining procedure, as output by paf_tile [SUPPORTED BY CACTUS ONLY]
 * cn   i   Chain id, indicating which chain a paf belongs to, as output by paf_chain [SUPPORTED BY CACTUS ONLY]
 */

typedef enum _cigarOp {
    match = 0,
    query_insert = 1, // substring in the query and not the target
    query_delete = 2 // substring in the target and not the query
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
    bool same_strand; // If 0 then query substring is reverse complement with respect to the target
    Cigar *cigar; // Ordered by the target sequence
    int64_t score; // the dp alignment score
    int64_t mapping_quality;
    int64_t num_matches;
    int64_t num_bases;
    int64_t tile_level; // this is a special tag created by paf_tile that indicates the "level" of the alignment in the
    // chaining (somewhat like the nesting of chains within nets in ucsc chains and nets).
    int64_t chain_id; // a tag to indicate which chain a paf belongs to
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
 * Prints a human readable alignment view of the paf
 */
void paf_pretty_print(Paf *paf, char *query_seq, char *target_seq, FILE *fh, bool include_alignment);

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
stList *paf_chain(stList *pafs, int64_t (*gap_cost)(int64_t, int64_t, void *), void *gap_cost_params,
                  int64_t max_gap_length, float percentage_to_trim);

/*
 * Gets the number of aligned bases in the alignment between the query
 * and the target according to the cigar alignment.
 */
int64_t paf_get_number_of_aligned_bases(Paf *paf);

/*
 * Removes prefix and suffix aligned bases.
 */
void paf_trim_ends(Paf *paf, int64_t end_bases_to_trim);

/*
 * Removes a given percentage of prefix and suffix aligned bases.
 */
void paf_trim_end_fraction(Paf *paf, float percentage);

/*
 * Breaks up a paf into its set of constituent matches
 */
stList *paf_shatter(Paf *paf);

/*
 * Calculate stats on the alignment
 */
void paf_stats_calc(Paf *paf, char *query_seq, char *target_seq,
                       int64_t *matches, int64_t *mismatches, int64_t *query_inserts, int64_t *query_deletes);


/*
 * Structure used to represent alignment coverage along a sequence
 */
typedef struct _sequenceCountArray {
    char *name; // Sequence name
    int64_t length; // Sequence length
    uint16_t *counts; // Array of counts, one for each base
} SequenceCountArray;

/*
 * Get a count array for the query sequence of a paf record, creating one if it doesn't exist.
 */
SequenceCountArray *get_alignment_count_array(stHash *seq_names_to_alignment_count_arrays, Paf *paf);

/*
 * Increase the count of alignment coverages for the query bases covered by a paf record.
 */
void increase_alignment_level_counts(SequenceCountArray *seq_count_array, Paf *paf);

typedef struct _interval {
    char *name;
    int64_t start, end, length;
} Interval;

/*
 * Decodes a fasta header into an interval.
 */
Interval *decode_fasta_header(char *fasta_header);

/*
 * Compare intervals;
 */
int cmp_intervals(const void *i, const void *j);

#endif /* ST_PAF_H_ */

