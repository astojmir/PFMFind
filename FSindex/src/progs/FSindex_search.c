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
  enum {T_SIMILARITY, T_QUASI_METRIC, T_METRIC} search_type 
    = T_SIMILARITY; 

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix query_seq"
	      " cutoff [-s | -d | -q];\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  matrix_full = argv[2]; 
  query_seq = argv[3];
  cutoff = atoi(argv[4]);

  if (argc >= no_args+2)
    {
      if (argv[5][0] == '-')
	switch (argv[5][1])
	  {
	  case 's':
	    search_type = T_SIMILARITY;
	    break;
	  case 'm':
	    search_type = T_METRIC;
	    break;
	  case 'q':
	    search_type = T_QUASI_METRIC;
	    break;
	  default:
	    search_type = T_SIMILARITY;
	    break;
	  }
    }
  query.len = strlen(query_seq);
  query.start = strdup(query_seq);

  FS_index = FS_INDEX_load(filename);
  ptable = FS_INDEX_get_ptable(FS_index);
  split_base_dir(matrix_full, &matrix_base, &matrix_dir);  
  S = SCORE_MATRIX_create(matrix_full, ptable); 
  hit_list = HIT_LIST_create(&query, FS_INDEX_get_database(FS_index), 
			     matrix_base, cutoff);

  switch (search_type)
    {
    case T_SIMILARITY:
      D = SCORE_MATRIX_S_2_Dquasi(S);
      FS_INDEX_set_matrix(FS_index, matrix_base, S, D);
      FS_INDEX_search(FS_index, hit_list, &query, cutoff,
		      FS_INDEX_S_process_bin,
		      FS_INDEX_S2QD_convert);
      HIT_LIST_sort_by_sequence(hit_list);
      HIT_LIST_sort_decr(hit_list);
       break;
    case T_QUASI_METRIC:
      D = SCORE_MATRIX_S_2_Dquasi(S);
      FS_INDEX_set_matrix(FS_index, matrix_base, S, D);
      FS_INDEX_search(FS_index, hit_list, &query, cutoff,
		      FS_INDEX_QD_process_bin, 
		      FS_INDEX_identity_convert);
      HIT_LIST_sort_by_sequence(hit_list);
      HIT_LIST_sort_incr(hit_list);
     break;
    case T_METRIC:
      D = SCORE_MATRIX_S_2_Dmax(S);
      FS_INDEX_set_matrix(FS_index, matrix_base, S, D);
      FS_INDEX_search(FS_index, hit_list, &query, cutoff,
		      FS_INDEX_QD_process_bin, 
		      FS_INDEX_identity_convert);
      HIT_LIST_sort_by_sequence(hit_list);
      HIT_LIST_sort_incr(hit_list);
      break;
    default:
      break;
    }

  HIT_LIST_print(hit_list, stdout, 0); 

  return EXIT_SUCCESS;
}
