#include "CuTest.h"
#include "ond.h"

//// Here is a an implementation of NeedlemanWunsch to compare the the O(ND) algorithm
typedef struct _NeedlemanWunsch {
    stList *string1;
    stList *string2;
    int64_t gap_score;
    int64_t match_score;
    int64_t mismatch_score;
    int64_t *edit_matrix;
    bool (*elements_equal)(void *, void *);
} NeedlemanWunsch;

static int64_t *get_cell(NeedlemanWunsch *nw, int64_t x, int64_t y) {
    /*
     * Get a cell from the edit matrix.
     */
    return &(nw->edit_matrix[x * (stList_length(nw->string2) + 1) + y]);
}

static int64_t get_score(NeedlemanWunsch *nw, int64_t i, int64_t j, int64_t *index) {
    /*
     * Gets the highest scoring alignment up to points i, j. index is set to 0 if a match/mismatch, 1 if a gap in i
     * and 2 if a gap in j;
     */
    int64_t a = *get_cell(nw, i-1, j-1) +
                (nw->elements_equal(stList_get(nw->string1, i-1), stList_get(nw->string2, j-1)) ? nw->match_score : nw->mismatch_score);
    int64_t b = *get_cell(nw, i-1, j) + nw->gap_score;
    int64_t c = *get_cell(nw, i, j-1) + nw->gap_score;
    if(a >= b && a >= c) {
        *index = 0;
        return a;
    }
    if(b > c) {
        *index = 1;
        return b;
    }
    *index = 2;
    return c;
}

NeedlemanWunsch *NeedlemanWunsch_construct(stList *string1, stList *string2,
                                           int64_t gap_score, int64_t match_score, int64_t mismatch_score,
                                           bool (*elements_equal)(void *, void *)) {
    /*
     * Finds an optimal global alignment of two strings.
    */
    NeedlemanWunsch *nw = st_calloc(1, sizeof(NeedlemanWunsch));

    nw->edit_matrix = st_calloc((stList_length(string1)+1) * (stList_length(string2)+1), sizeof(int64_t));
    // Numpy matrix representing edit matrix
    // Preinitialized to have zero values

    nw->elements_equal = elements_equal;

    // Initialize prefix indels
    for (int64_t i = 1; i <= stList_length(string1); i++) {
        *get_cell(nw, i, 0) = gap_score * i;
    }

    for (int64_t j = 1; j <= stList_length(string2); j++) {
        *get_cell(nw, 0, j) = gap_score * j;
    }

    // Fill in remaining edit matrix
    for (int64_t i = 1; i <= stList_length(string1); i++) {
        for (int64_t j = 1; j <= stList_length(string2); j++) {
            int64_t k;
            *get_cell(nw, i, j) = get_score(nw, i, j, &k);
        }
    }

    return nw;
}

int64_t NeedlemanWunsch_get_alignment_score(NeedlemanWunsch *nw) {
    /*
     * Return the alignment score
    */
    return *get_cell(nw, stList_length(nw->string1), stList_length(nw->string2));
}

stList *NeedlemanWunsch_get_alignment(NeedlemanWunsch *nw) {
    /* Returns an optimal global alignment of two strings. Aligned
     * is returned as an ordered list of aligned pairs.
     *
     * e.g. For the two strings GATTACA and TACA an global alignment is
     *       is GATTACA
     *          ---TACA
     * This alignment would be returned as:
     * [(3, 0), (4, 1), (5, 2), (6, 3)]
     */
    int64_t i = stList_length(nw->string1);
    int64_t j = stList_length(nw->string2);

    stList *aligned_pairs = stList_construct();
    while (i > 0 && j > 0) {
        int64_t m;
        get_score(nw, i, j, &m);
        if (m == 0) {
            i--; j--;
            if(nw->elements_equal(stList_get(nw->string2, j), stList_get(nw->string1, i))) { // Only include actual matches
                // (not mismatches) in the alignment
                stList_append(aligned_pairs, stList_get(nw->string2, j));
                stList_append(aligned_pairs, stList_get(nw->string1, i));
            }
        } else if (m == 1) {
            i -= 1;
        } else {
            assert(m == 2);
            j -= 1;
        }
    }
    // Put in the right order
    stList_reverse(aligned_pairs);

    return aligned_pairs;
}

static stList *get_random_string() {
    /*
     * Generate a random list of As and Ts to align.
     */
    stList *random_string = stList_construct();
    int64_t length = st_randomInt(1, 10);
    for (int64_t i = 0; i < length; i++) {
        stList_append(random_string, stString_copy(st_random() > 0.5 ? "A" : "T"));
    }
    return random_string;
}

static bool elements_equal(void *a, void *b) {
    // Do the rows match
    return strcmp((char *)a, (char *)b) == 0;
}

void test_ond(CuTest *testCase) {
    // Here we test that the WFA alignment score agrees with Needleman Wunsch for
    // 100 randomly chosen small test examples
    for (int64_t test = 0; test < 100; test++) {
        stList *x = get_random_string();
        stList *y = get_random_string();

        int64_t mismatch_score = st_randomInt(1, 10);
        int64_t gap_score = st_randomInt(1, 10);

        NeedlemanWunsch *nw = NeedlemanWunsch_construct(x, y, 0, -mismatch_score, -gap_score, elements_equal);
        WFA *wfa = WFA_construct(x, y, elements_equal, mismatch_score, gap_score);

        // print("WFA Alignment", wfa.get_alignment())
        // print("NW Alignment", nw.get_alignment())

        // The scores of the alignment should be the same
        CuAssertIntEquals(testCase, -NeedlemanWunsch_get_alignment_score(nw), WFA_get_alignment_score(wfa));

        // The actual alignments should be the same
        stList *nw_alignment = NeedlemanWunsch_get_alignment(nw);
        stList *wfa_alignment = WFA_get_alignment(wfa);
        CuAssertTrue(testCase, stList_length(nw_alignment) % 2 == 0);
        CuAssertIntEquals(testCase, stList_length(nw_alignment), stList_length(wfa_alignment));
        for(int64_t i=0; i<stList_length(nw_alignment); i++) {
            CuAssertTrue(testCase, stList_get(nw_alignment, i) == stList_get(wfa_alignment, i));
        }
    }
}

CuSuite* ond_test_suite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_ond);
    return suite;
}
