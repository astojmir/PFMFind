#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "randseq.h"
#include "FSindex.h"
#include "hit_list.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;

int main(int argc, char **argv)
{
  const char *filename;
  const char *count_file;
  char *matrix_full;
  char *matrix_base;
  char *matrix_dir;
  int no_args = 4;


  int cutoff;
  int self_sim;
  int no_runs;
  HIT_LIST_t *hit_list;
  BIOSEQ *query;
  FS_INDEX_t *FS_index;
  SEQ_GENERATOR_t *seq_generator;
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  char *seq_heap;
  ULINT len;
  ULINT i;

  clock_t t0;
  clock_t t_loading;
  clock_t t_search;
  clock_t t_recursive;
  clock_t t_recursive1;

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix count_file"
	      "no_runs cutoff\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  matrix_full = argv[2]; 
  count_file = argv[3];
  no_runs = atoi(argv[4]);
  cutoff = atoi(argv[5]);

  t0 = clock();

  FS_index = FS_INDEX_load(filename);

  len = FS_index->frag_len;
  seq_heap = mallocec(len+1);
  ptable = FS_INDEX_get_ptable(FS_index);
  S = SCORE_MATRIX_create(matrix_full, ptable); 
  D = SCORE_MATRIX_S_2_Dquasi(S);
  split_base_dir(matrix_full, &matrix_base, &matrix_dir);  
  FS_INDEX_set_matrix(FS_index, matrix_base, S, D);
  seq_generator = SEQ_GENERATOR_create(count_file, ptable);

  hit_list = HIT_LIST_create(&query, FS_INDEX_get_database(FS_index), 
			     matrix_base, 0);

  query = mallocec( no_runs * sizeof(BIOSEQ));
  for (i=0; i < no_runs; i++)
    {
      seq_heap = mallocec(len+1);
      SEQ_GENERATOR_rand_seq(seq_generator, query+i, len,
			     seq_heap);
    }
  t_loading = clock() - t0;
  t0 = clock();

  for (i=0; i < no_runs; i++)
    {
      /* Run search */
      HIT_LIST_reset(hit_list, query+i, FS_index->s_db,
		     matrix_base, cutoff);
      FS_INDEX_search(FS_index, hit_list, query+i, cutoff,
		      FS_INDEX_S_process_bin,
		      FS_INDEX_S2QD_convert);
    }
  t_search = clock() - t0;
 
  t0 = clock();
  for (i=0; i < no_runs; i++)
    {
       /* Run search */
      HIT_LIST_reset(hit_list, query+i, FS_index->s_db,
		     matrix_base, cutoff);
      FS_INDEX_recursive_search(FS_index, hit_list, query+i, cutoff,
				FS_INDEX_S_process_bin,
				FS_INDEX_S2QD_convert);
    }
  t_recursive = clock() - t0;

  t0 = clock();
  for (i=0; i < no_runs; i++)
    {
       /* Run search */
      HIT_LIST_reset(hit_list, query+i, FS_index->s_db,
		     matrix_base, cutoff);
      FS_INDEX_recursive_search1(FS_index, hit_list, query+i, cutoff,
				FS_INDEX_S_process_bin,
				FS_INDEX_S2QD_convert);
    }
  t_recursive1 = clock() - t0;

  fprintf(stdout, "Index file: %s\n", filename);
  fprintf(stdout, "Matrix: %s\n", matrix_base);
  fprintf(stdout, "Similarity cutoff: %d\n", cutoff);
  fprintf(stdout, "Number of runs: %d\n", no_runs);
  fprintf(stdout, "Loading time: %.6f s\n", 
	  (double) t_loading / CLOCKS_PER_SEC);
  fprintf(stdout, "Non-recursive search time: %.6f s\n", 
	  (double) t_search / CLOCKS_PER_SEC);
  fprintf(stdout, "Recursive search time: %.6f s\n", 
	  (double) t_recursive / CLOCKS_PER_SEC);
  fprintf(stdout, "New Recursive search time: %.6f s\n", 
	  (double) t_recursive1 / CLOCKS_PER_SEC);
return EXIT_SUCCESS;
}
