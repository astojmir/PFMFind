#include "SAindex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "sarray.h"

extern void alloc_FSSRCH(FSSRCH *FSS, int m, int *nbpt);

static
int *FS_TABLE_order_convert_heap(FS_TABLE *ptable, const char *s, int len)
{
  int i, k;
  int *r;
  int amap[A_SIZE];
  int freq[A_SIZE];
  int ainv[A_SIZE];
  const char *t = s;

  r= mallocec(len*sizeof(int));
  memset(freq, 0, A_SIZE*sizeof(int));

  /* Must scan the whole dataset to determine what letters are used
     so we don't have any 'holes' in order. This is a requirement
     of sarray function */
  for (i=0; i < A_SIZE; i++) 
    ainv[ptable->rank[0][i]] = i;
  
  for (i=0; i < len; i++, t++)
    freq[*t & A_SIZE_MASK]++;

  for (k=0, i=0; i < A_SIZE; i++) {
    if (freq[ainv[i]]) amap[ainv[i]] = k++;
    else amap[ainv[i]] = -1;
  }

  for (i=0; i < len; i++, s++)
    r[i] = amap[*s & A_SIZE_MASK];
  return r;
} 


SAINDX *SA_INDEX_create(const char *database, int print_flag)
{
  EXCEPTION *except;
  SAINDX *FSI = NULL;
  const char *alpha = "STANILVMKRDEQWFYHGPC";
  const char **sepn = &alpha;

  FS_TABLE *ptable;

  char *_base_;
  char *_end_;

  time_t time0 = time(NULL);
  time_t time1;
  double dt;
  char *db_dir = NULL;
  char *db_base = NULL;

  ULINT i;

  uint32_t *tmp_size = NULL;
  int *tmp_lcpx = NULL;

  Try {
    /* Create fake partition table */
    ptable = FS_TABLE_init(sepn, 1);
    if (ptable == NULL)
      Throw FSexcept(BAD_ARGS, "FS_INDEX_create():"
		     " Could not create partition table\n");

    /* Initialise FS_index */
    FSI = callocec(1, sizeof(FSINDX));
    FSI->ptable = ptable;
    path_split(database, &db_dir, &db_base);
    FSI->db_name = strdup(db_base);
    FSI->index_name = strdup("Suffix Array");
    FSI->m = 1; /* The same as the length of ptable */

    /* Load database */
    FSI->s_db = fastadb_open(database); 
 
    /* Global variables */
    _base_ = fastadb_data_pter(FSI->s_db, 0);  
    _end_ = fastadb_end_heap(FSI->s_db);
    
    *_end_ = '_';
    FSI->no_seqs = _end_ - _base_ + 1;
    FSI->db_no_frags = FSI->no_seqs;

    /* Just to make it fully compatible with FSINDX */
    FSI->no_bins = 1;
    FSI->bins = callocec(FSI->no_bins+1, sizeof(uint32_t));
    FSI->bins[0] = 0;
    FSI->bins[1] = FSI->no_seqs;
    
    if (print_flag) {
      printf("***** Suffix Array *****\n");
      printf("Converting sequences ...\n");
    }
    FSI->oa = FS_TABLE_order_convert_heap(FSI->ptable, _base_,
					  FSI->no_seqs);
    if (print_flag) {
      printf("Creating suffix array ...\n");
    }
    if (sarray(FSI->oa, FSI->no_seqs) == -1)
      Throw FSexcept(RUNTIME_ERR, "FS_INDEX_create():"
		     " Could not create suffix array.\n");
    if (print_flag) {
      printf("Creating lcp array ...\n");
    }
    tmp_lcpx = lcp(FSI->oa, _base_, FSI->no_seqs); 
    if (!tmp_lcpx)
      Throw FSexcept(RUNTIME_ERR, "FS_INDEX_create():"
		     " Could not create lcp array.\n");
    FSI->lcpb = mallocec(FSI->no_seqs);
    for (i=0; i < FSI->no_seqs; i++)
      FSI->lcpb[i] = min(tmp_lcpx[i], 255);
    free(tmp_lcpx);
    *_end_='\0';

    /* Take time */
    time1 = time(NULL);
    dt = difftime(time1, time0)/60.0;
    FSI->dtime = dt;

    /* Allocate search specific variables */
    for (i=0; i < THREADS; i++)
      alloc_FSSRCH(FSI->FSS+i, FSI->m, FSI->ptable->nbpt);
  }
  Catch (except) {
    free(tmp_lcpx);
    free(tmp_size);
    SA_INDEX_destroy(FSI);
    Throw except;
  }
  free(db_dir);
  free(db_base);
  return FSI; 
}

void SA_INDEX_destroy(SAINDX *SAI)
{
  FS_INDEX_destroy(SAI);
}

int SA_INDEX_save(SAINDX *SAI, const char *filename)
{
  return FS_INDEX_save(SAI, filename);
}

SAINDX *SA_INDEX_load(const char *filename)
{
  return FS_INDEX_load(filename);
}

char *SA_INDEX_sprint_stats(FSINDX *SAI, int options)
{
  char *buf = NULL;
  int size = 0;
  int len = 0;

  absprintf(&buf, &size, &len, "\n*** SA_INDEX Statistics ***\n");
  absprintf(&buf, &size, &len, "Database: %s\n", SAI->db_name);
  absprintf(&buf, &size, &len, "Database size: %ld\n", SAI->s_db->length); 
  absprintf(&buf, &size, &len, "Number of sequences: %ld\n\n", SAI->s_db->no_seq);  
  absprintf(&buf, &size, &len, "Alphabet: %s\n", SAI->ptable->alphabet);
  absprintf(&buf, &size, &len, "Creation Time: %.2f mins\n\n", SAI->dtime);
  absprintf(&buf, &size, &len, "Suffix array size: %ld\n", SAI->no_seqs);
  return buf;
}

void SA_INDEX_fprint_stats(FSINDX *SAI, FILE *fp, int options)
{
  char *buf = SA_INDEX_sprint_stats(SAI, options);
  fprintf(fp, "%s", buf);
  free(buf);
}

HIT_LIST_t *SAINDX_rng_srch(SAINDX *SAI, BIOSEQ *query, SCORE_MATRIX *M,
			    int range, int conv_type, HIT_LIST_t *HL)
{
  return FSINDX_rng_srch(SAI, query, M, range, conv_type, HL, SUFFIX_ARRAY, SARRAY);
}

HIT_LIST_t *SAINDX_kNN_srch(SAINDX *SAI, BIOSEQ *query, SCORE_MATRIX *M,
			    int kNN, HIT_LIST_t *HL)
{
  return FSINDX_kNN_srch(SAI, query, M, kNN, HL, SUFFIX_ARRAY, SARRAY);
}
