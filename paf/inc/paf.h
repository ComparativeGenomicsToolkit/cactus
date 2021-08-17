/*
 * stPaf.h
 *
 * Functions for manipulating paf alignments.
 */

#ifndef ST_PAF_H_
#define ST_PAF_H_

#include "sonLib.h"

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
    double score;
    double mapping_quality;
    int64_t num_matches;
    int64_t num_bases;
    bool is_primary;
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
 * Write a paf record as a string.
 */
char *paf_print(Paf *paf);

/*
 * Writes a PAF alignment to the given file.
 */
void paf_write(Paf *paf, FILE *fh);

#endif /* ST_PAF_H_ */

