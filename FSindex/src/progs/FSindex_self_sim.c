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
  BIOSEQ query;
  FS_INDEX_t *FS_index;
  SEQ_GENERATOR_t *seq_generator;
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  char *seq_heap;
  ULINT len;
  ULINT i;

  ULINT qm_size;


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

  fprintf(stdout, "Similarity cutoff: %d\n\n", cutoff);
  fprintf(stdout, "%*.*s %10s %10s %10s\n",
	  (int) len , (int) len, "Sequence", "self-sim", "qm-range",
	  "no_neighbours"); 

  for (i=0; i < no_runs; i++)
    {
      /* Generate sequence */
      SEQ_GENERATOR_rand_seq(seq_generator, &query, len,
			     seq_heap);
      self_sim = SCORE_MATRIX_evaluate(S, &query, &query);
      /* Run search */
      HIT_LIST_reset(hit_list, &query, FS_index->s_db,
		     matrix_base, cutoff);
      FS_INDEX_search(FS_index, hit_list, &query, cutoff,
		      FS_INDEX_S_process_bin,
		      FS_INDEX_S2QD_convert);

      /* Print results */
      qm_size = HIT_LIST_get_seqs_hits(hit_list);
      fprintf(stdout, "%*.*s %10d %10d %10ld\n",
	      (int) query.len, (int) query.len,  query.start,
	      self_sim, hit_list->converted_range, qm_size);
    }
 
 return EXIT_SUCCESS;
}
