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
  const char *index_file;
  const char *count_file;
  SEQ_GENERATOR_t *seq_generator;
  ULINT no_runs;
  ULINT i;
  ULINT len;
  char *seq_heap;
  ULINT total_indexed_frags;
  ULINT no_scanned = 0;
  ULINT D_cutoff = 0;



  if (argc < 4)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file count_file"
	      " no_runs\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  index_file = argv[1];
  count_file = argv[2];
  no_runs = atoi(argv[3]);

  FS_index = FS_INDEX_load(index_file);
  seq_generator = SEQ_GENERATOR_create(count_file, FS_index->ptable);
  len = FS_index->frag_len;
  seq_heap = mallocec(len+1);
#if 0
  
  hit_list = HIT_LIST_create(query, FS_index->s_db, 
			     FS_index->matrix, D_cutoff);

  total_indexed_frags = FSindex->HT->no_seqs;
  while (no_scanned < total_indexed_frags)
    {
      for (i = 0; i < no_runs; i++)
	{
	  SEQ_GENERATOR_rand_seq(seq_generator, &query, len,
				 seq_heap); 
	  if (!FS_INDEX_search(FS_index, hit_list, &query, 
			       D_cutoff, 
			       FS_INDEX_count_only_process_bin))
	    {
	      fprintf(stderr,"Error in FS_INDEX_search\n");
	      exit(EXIT_FAILURE);
	    }
	}
      /* Print */
      D_cutoff++;
      no_scanned;
    }
#endif
#if 0
  hit_list = FS_INDEX_QD_search(FS_index, &query, cutoff);
  HIT_LIST_sort_by_sequence(hit_list);
  HIT_LIST_sort_incr(hit_list);
  HIT_LIST_print(hit_list, stdout, 0); 
#endif

 return EXIT_SUCCESS;
}
