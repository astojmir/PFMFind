#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "FSindex.h"
#include "hit_list.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;

int main(int argc, char **argv)
{
  const char *query_seq;
  int cutoff;
  HIT_LIST_t *hit_list;
  BIOSEQ query;
  FS_INDEX_t *FS_index;
  const char *filename;
  char *matrix_full;
  char *matrix_base;
  char *matrix_dir;
  int no_args = 4;
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix query_seq"
	      " cutoff\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  matrix_full = argv[2]; 
  query_seq = argv[3];
  cutoff = atoi(argv[4]);

  query.len = strlen(query_seq);
  query.start = strdup(query_seq);

  FS_index = FS_INDEX_load(filename);

  ptable = FS_INDEX_get_ptable(FS_index);
  S = SCORE_MATRIX_create(matrix_full, ptable); 
  D = SCORE_MATRIX_S_2_Dquasi(S);
  split_base_dir(matrix_full, &matrix_base, &matrix_dir);  
  FS_INDEX_set_matrix(FS_index, matrix_base, S, D);

  hit_list = FS_INDEX_QD_search(FS_index, &query, cutoff);
  HIT_LIST_sort_by_sequence(hit_list);
  HIT_LIST_sort_incr(hit_list);
  HIT_LIST_print(hit_list, stdout, 0); 

 return EXIT_SUCCESS;
}
