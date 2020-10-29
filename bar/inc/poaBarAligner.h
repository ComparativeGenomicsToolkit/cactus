/*
 * Copyright (C) 2009-2020 by Benedict Paten, Glenn Hickey
 *
 * Released under the MIT license, see LICENSE.txt
 */

#ifndef POA_END_ALIGNER_H_
#define POA_END_ALIGNER_H_

#include "sonLib.h"
#include "cactus.h"
#include "stPinchIterator.h"

/**
 * Object representing a multiple sequence alignment
 */
typedef struct _Msa {
    int64_t seq_no; // number of sequences
    int *seq_lens; // length of the sequences
    char **seqs; // sequences as ASCII characters
    int column_no; // number of columns in the msa
    uint8_t **msa_seq; // the msa matrix of the aligned sequences
} Msa;

/**
 * Convert a base in the POA alphabet to an ASCII base.
 */
char msa_to_base(uint8_t n);

/**
 * Convert an ASCII base into the corresponding POA alphabet integer
 */
uint8_t msa_to_byte(char c);

/**
 * Clean up an MSA
 */
void msa_destruct(Msa *msa);

/**
 * Pretty prints an MSA object
 */
void msa_print(Msa *msa, FILE *f);

/**
 * Creates a partial order alignment
 * @param seqs An array of DNA string
 * @param seq_lens An array giving the string lengths
 * @param seq_no The number of strings
 * @return An msa of the strings.
 */
Msa *msa_make_partial_order_alignment(char **seqs, int *seq_lens, int64_t seq_no);

/**
 * Takes a set of ends and returns a set of consistent multiple alignments,
 * one for each of them.
 *
 * Each end is a an array of strings, oriented left-to-right from the end.
 * The strings connect the ends.
 * Each string therefore has a reverse complement in the set of strings.
 *
 * The MSAs are consistent with one another if each base in each string and its reverse complement is only
 * aligned to other bases in one MSA.
 *
 * @param end_no The number of ends
 * @param end_lengths The number of strings incident with each the end
 * @param end_strings The strings connecting the ends
 * @param end_string_lengths Length of the strings connecting the ends
 * @param right_end_indexes For each string, the index of the right end that it is connecting
 * @param right_end_row_indexes For each string, the index of the row of its reverse complement
 * @return A consistent Msa for each end
 */
Msa **make_consistent_partial_order_alignments(int64_t end_no, int64_t *end_lengths, char ***end_strings,
        int **end_string_lengths, int64_t **right_end_indexes, int64_t **right_end_row_indexes);

/**
 * Represents a gapless alignment of a set of sequences.
 */
typedef struct _AlignmentBlock {
    int64_t subsequenceIdentifier; // The name of the sequence
    int64_t position; // The start position
    bool strand;
    int64_t length;
    struct _AlignmentBlock *next;
} AlignmentBlock;

void alignmentBlock_destruct(AlignmentBlock *alignmentBlock);

/**
 * Makes alignments of the the unaligned sequence using the bar algorithm.
 *
 * Returns a list of AlignmentBlock ojects
 */
stList *make_flower_alignment_poa(Flower *flower, bool pruneOutStubAlignments);

/**
 * Create a pinch iterator for a list of alignment blocks.
 * @param alignment_blocks A list of AlignmentBlock instances
 * @return A pinch iterator for all the alignments in alignment_blocks
 */
stPinchIterator *stPinchIterator_constructFromAlignedBlocks(stList *alignment_blocks);

#endif
