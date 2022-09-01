#include "line_iterator.h"

struct _LI {
    FILE *fh;
    char *line;
};

LI *LI_construct(FILE *fh) {
    LI *li = st_calloc(1, sizeof(LI));
    li->fh = fh;
    li->line = li->line = stFile_getLineFromFile(fh);
    return li;
}

void LI_destruct(LI *li) {
    free(li);
}

char *LI_get_next_line(LI *li) {
    char *l = li->line;
    li->line = stFile_getLineFromFile(li->fh);
    return l;
}

char *LI_peek_at_next_line(LI *li) {
    return li->line;
}

