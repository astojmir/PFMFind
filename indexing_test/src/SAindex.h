#include "FSindex.h"


/* Uses suffix array routines for construction,
   otherwise a straight copy of FSindex */

typedef FSINDX SAINDX;

/* Functions */

SAINDX *SA_INDEX_create(const char *database, int print_flag);
void SA_INDEX_destroy(SAINDX *SAI);

int SA_INDEX_save(SAINDX *SAI, const char *filename);
SAINDX *SA_INDEX_load(const char *filename);

char *SA_INDEX_sprint_stats(FSINDX *SAI, int options);
void SA_INDEX_fprint_stats(FSINDX *SAI, FILE *fp, int options);

HIT_LIST_t *SAINDX_rng_srch(SAINDX *SAI, BIOSEQ *query, SCORE_MATRIX *M,
			    int range, int conv_type, HIT_LIST_t *HL);

HIT_LIST_t *SAINDX_kNN_srch(SAINDX *SAI, BIOSEQ *query, SCORE_MATRIX *M,
			    int kNN, HIT_LIST_t *HL);

