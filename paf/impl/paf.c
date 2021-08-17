#include "paf.h"

/*
 * Library functions for manipulating paf files.
 */

void paf_destruct(Paf *paf) {
    Cigar *c = paf->cigar;
    while(c != NULL) {
        Cigar *c2 = c->next;
        free(c);
        c = c2;
    }
    free(paf);
}

Paf *paf_parse(char *paf_string) {
    Paf *paf = st_calloc(1, sizeof(Paf));

    return paf;
}

Paf *paf_read(FILE *fh) {
    char *c = stFile_getLineFromFile(fh);
    if(c == NULL) {
        return NULL;
    }
    Paf *paf = paf_parse(c);
    free(c);
    return paf;
}



char *paf_print(Paf *paf) {
    /*
     * PAF (from: https://lh3.github.io/minimap2/minimap2.html)
     * 1    string	Query sequence name
     * 2	int	Query sequence length
     * 3	int	Query start coordinate (0-based)
     * 4	int	Query end coordinate (0-based)
     * 5	char	‘+’ if query/target on the same strand; ‘-’ if opposite
     * 6	string	Target sequence name
     * 7	int	Target sequence length
     * 8	int	Target start coordinate on the original strand
     * 9	int	Target end coordinate on the original strand
     * 10	int	Number of matching bases in the mapping
     * 11	int	Number bases, including gaps, in the mapping
     * 12	int	Mapping quality (0-255 with 255 for missing)
     *
     * Tags:
     *
     * tp	A	Type of aln: P/primary, S/secondary and I,i/inversion
     * cm	i	Number of minimizers on the chain
     * s1	i	Chaining score
     * s2	i	Chaining score of the best secondary chain
     * NM	i	Total number of mismatches and gaps in the alignment
     * MD	Z	To generate the ref sequence in the alignment
     * AS	i	DP alignment score
     * SA	Z	List of other supplementary alignments
     * ms	i	DP score of the max scoring segment in the alignment
     * nn	i	Number of ambiguous bases in the alignment
     * ts	A	Transcript strand (splice mode only)
     * cg	Z	CIGAR string (only in PAF)
     * cs	Z	Difference string
     * dv	f	Approximate per-base sequence divergence
     * de	f	Gap-compressed per-base sequence divergence
     * rl	i	Length of query regions harboring repetitive seeds
     *
     */
    char *c = stString_print("%s %i %i %i %c %s %i %i %i %i %i %i\n",
                             paf->query_name, paf->query_length, paf->query_start, paf->query_end, paf->same_strand,
                             paf->target_name, paf->target_length, paf->target_start, paf->target_end,
                             paf->num_matches, paf->num_bases, paf->mapping_quality);

    return c;
}

void paf_write(Paf *paf, FILE *fh) {
    char *c = paf_print(paf);
    fprintf(fh, c);
    free(c);
}


