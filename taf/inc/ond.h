#ifndef STOND_H_
#define STOND_H_

#include "sonLib.h"

typedef struct _WFA WFA;

WFA *WFA_construct(stList *string1, stList *string2, bool (*elements_equal)(void *, void *),
                   int64_t gap_score, int64_t mismatch_score);

void WFA_destruct(WFA *wfa);

int64_t WFA_get_alignment_score(WFA *wfa);

stList *WFA_get_alignment(WFA *wfa);

#endif /* STOND_H_ */

