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

static
ULINT get_max_distance_hit(BIOSEQ *query, HIT_LIST_t *hit_list,
			   SCORE_MATRIX_t *D)
{
  ULINT no_hits = HIT_LIST_get_seqs_hits(hit_list);
  ULINT i;
  ULINT delta = 0;
  SEQ_HIT_t *hit;

  for (i=0; i < no_hits; i++)
    {
      hit = HIT_LIST_get_hit(hit_list, i);
      delta = max(delta, SCORE_MATRIX_evaluate(D, query,
					       hit->subject)); 
    }
  return delta;
}

int main(int argc, char **argv)
{
  const char *filename;
  const char *count_file;
  const char *out_file = NULL;
  char *matrix_full;
  char *matrix_base;
  char *matrix_dir;
  int no_args = 7;

  int K0;
  int K1;
  int step;
  FILE *stream = stdout;
  ULINT one_percent;

  int no_runs;
  HIT_LIST_t *hit_list;
  BIOSEQ query;
  SEQ_GENERATOR_t *seq_generator;
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  SCORE_MATRIX_t *MD;
  FS_HASH_TABLE_t *HT;
  char *seq_heap;
  ULINT len;
  ULINT i;

  ULINT eps;
  ULINT qFS_visited;
  ULINT qFS_hits;
  ULINT qS_visited;
  ULINT qm_size;
  ULINT quv;
  ULINT quh;

  ULINT delta;
  ULINT mFS_visited;
  ULINT mFS_hits;
  ULINT mS_visited;
  ULINT m_size;
  ULINT muv;
  ULINT muh;

  double ratio;


  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix_file count_file"
	      " no_runs start_eps step_size end_eps <out_file>\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  matrix_full = argv[2]; 
  count_file = argv[3];
  no_runs = atoi(argv[4]);
  K0 = atoi(argv[5]);
  step = atoi(argv[6]);
  K1 = atoi(argv[7]);

 if (argc >= no_args+2)
    {
      out_file = argv[8];
      if ((stream = fopen(out_file, "w")) == NULL)
	{
	  fprintf(stderr, "Unable to open %s\n ", out_file);
	  exit(EXIT_FAILURE);
	} 
    } 

  FS_INDEX_load(filename);
  len = FS_INDEX_get_frag_len();
  seq_heap = mallocec(len+1);
  ptable = FS_INDEX_get_ptable();
  S = SCORE_MATRIX_create(matrix_full, ptable); 
  D = SCORE_MATRIX_S_2_Dquasi(S);
  MD = SCORE_MATRIX_S_2_Dmax(S);
  split_base_dir(matrix_full, &matrix_base, &matrix_dir);  
  seq_generator = SEQ_GENERATOR_create(count_file, ptable);

  hit_list = HIT_LIST_create(&query, FS_INDEX_get_database(), 
			     matrix_base, 0);
  if (no_runs < 100)
    one_percent = 1;
  else
    one_percent = no_runs/100;

  fprintf(stream, "***** FSsearch random runs: %s"
	  " database ******\n", FS_INDEX_get_db_name());
  fprintf(stream, "No runs: %d\n", no_runs);
  HT = FS_INDEX_get_hash_table();
  fprintf(stream, "Total fragments: %ld\n", FS_HASH_TABLE_get_total_seqs(HT));
  fprintf(stream, "Total bins: %ld\n", FS_HASH_TABLE_get_no_bins(HT));
  fprintf(stream, "Clustering: %s\n", FS_index_get_alphabet());
  fprintf(stream, "Matrix: %s\n", matrix_base);

  if (K0 <= 0)
    K1 = 1;

  while (K0 <= K1)
    {
      fprintf(stdout,"\neps = %d\n", K0);
      fprintf(stream,"\neps = %d\n", K0);
      fprintf(stream, "%*.*s %6s %10s %10s %10s %10s %10s %10s"
	      " %6s %10s %10s %10s %10s %10s %10s %10s\n",
	      (int) len , (int) len , "Sequence", 
	      "qeps", "qFS", "qFSh", "qsv", "qsh", "quv", "quh",
	      "delta","mFS", "mFSh", "msv", "msh", "muv", "muh",
	      "ratio");
      eps = K0;
      for (i=0; i < no_runs; i++)
	{
	  /* Generate sequence */
	  SEQ_GENERATOR_rand_seq(seq_generator, &query, len, seq_heap);  
	  HIT_LIST_reset(hit_list, &query, FS_INDEX_get_database(),
			 matrix_base, eps);

	  /* Find quasi-metric ball */
	  FS_INDEX_search(hit_list, &query, S, D, eps,
			  FS_INDEX_QD_process_bin, 
			  FS_INDEX_identity_convert);
	  
	  qFS_visited = hit_list->FS_seqs_visited;
	  qFS_hits = hit_list->FS_seqs_hits;
	  qS_visited = hit_list->seqs_visited;
	  qm_size = hit_list->seqs_hits;
	  quv = hit_list->useqs_visited;
	  quh = hit_list->useqs_hits;

	  /* Find metric ball */
	  delta = get_max_distance_hit(&query, hit_list, MD);

	  HIT_LIST_reset(hit_list, &query, FS_INDEX_get_database(),
			 matrix_base, delta);
	  FS_INDEX_search(hit_list, &query, S, MD, delta,
			  FS_INDEX_QD_process_bin, 
			  FS_INDEX_identity_convert);
	  m_size = hit_list->seqs_hits;
	  mFS_visited = hit_list->FS_seqs_visited;
	  mFS_hits = hit_list->FS_seqs_hits;
	  mS_visited = hit_list->seqs_visited;
	  muv = hit_list->useqs_visited;
	  muh = hit_list->useqs_hits;

	  /* Print results */
	  ratio = (double) m_size / qm_size;

	  fprintf(stream, "%*.*s %6ld %10ld %10ld %10ld %10ld %10ld %10ld"
		  " %6ld %10ld %10ld %10ld %10ld %10ld %10ld %10.2f\n",
		  (int) len, (int) len, query.start, eps, qFS_visited, qFS_hits,
		  qS_visited, qm_size, quv, quh, delta,
		  mFS_visited, mFS_hits, mS_visited, m_size, muv, muh,
		  ratio);

	  /* Print progress bar */
	  printbar(stdout, i+1, one_percent, 50);  
	}
      K0 += step;       
    }

  fprintf(stdout,"\n");
  if (out_file != NULL)
    fclose(stream);
  return EXIT_SUCCESS;
}
