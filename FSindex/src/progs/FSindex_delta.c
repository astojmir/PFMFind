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
  char *matrix_full;
  char *matrix_base;
  char *matrix_dir;
  int no_args = 5;


  int k_neighbours;
  int no_runs;
  HIT_LIST_t *hit_list;
  BIOSEQ query;
  SEQ_GENERATOR_t *seq_generator;
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  SCORE_MATRIX_t *MD;
  char *seq_heap;
  ULINT len;
  ULINT i;

  ULINT eps;
  ULINT qFS_visited;
  ULINT qFS_hits;
  ULINT qS_visited;
  ULINT qm_size;
  ULINT delta;
  ULINT mFS_visited;
  ULINT mFS_hits;
  ULINT mS_visited;
  ULINT m_size;
  double ratio;
  double qm_size_av = 0.0;
  double m_size_av = 0.0;
  double ratio_av = 0.0;
  double eps_av = 0.0;
  double delta_av = 0.0;


  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix count_file"
	      " no_neighbours no_runs\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  matrix_full = argv[2]; 
  count_file = argv[3];
  k_neighbours = atoi(argv[4]);
  no_runs = atoi(argv[5]);

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

  fprintf(stdout, "%*.*s %5s %12s %12s %12s %12s %5s %12s %12s"
	  " %12s %12s %10s\n", (int) len , (int) len , "Sequence", 
	  "eps", "qFS_visited", "qFS_hits", "qS_visited", "qS_hits",
	  "delta", "mFS_visited", "mFS_hits", "mS_visited", "mS_hits", 
	  "ratio");

  for (i=0; i < no_runs; i++)
    {
      /* Generate sequence */
      SEQ_GENERATOR_rand_seq(seq_generator, &query, len, seq_heap);  

      /* Find quasi-metric ball */
      FS_INDEX_kNN_search(hit_list, &query, S, D, k_neighbours);

      eps = hit_list->range;
      qFS_visited = hit_list->FS_seqs_visited;
      qFS_hits = hit_list->FS_seqs_hits;
      qS_visited = hit_list->seqs_visited;
      qm_size = HIT_LIST_get_seqs_hits(hit_list);


      /* Find metric ball */
      delta = get_max_distance_hit(&query, hit_list, MD);

      FS_INDEX_search(hit_list, &query, S, MD, delta,
		      FS_INDEX_QD_process_bin, 
		      FS_INDEX_identity_convert);
      m_size = HIT_LIST_get_seqs_hits(hit_list);
      mFS_visited = hit_list->FS_seqs_visited;
      mFS_hits = hit_list->FS_seqs_hits;
      mS_visited = hit_list->seqs_visited;


      /* Print results */
      ratio = (double) m_size / qm_size;

      fprintf(stdout, "%*.*s %5ld %12ld %12ld %12ld %12ld %5ld %12ld"
	      " %12ld %12ld %12ld %10.2f\n", (int) query.len, 
	      (int) query.len, query.start,
	      eps, qFS_visited, qFS_hits, qS_visited, qm_size, delta,
	      mFS_visited, mFS_hits, mS_visited, m_size, ratio);
      qm_size_av += qm_size;
      m_size_av += m_size;
      ratio_av += ratio;
      eps_av += eps;
      delta_av += delta;
    }

  qm_size_av /= no_runs;
  m_size_av /= no_runs;
  ratio_av /= no_runs;
  eps_av /= no_runs;
  delta_av /= no_runs;

  /* Print average */
  fprintf(stdout, "\n\n **** AVERAGES ****\n");
  fprintf(stdout, "Average qm radius: %.2f\n", eps_av);
  fprintf(stdout, "Average size of qm-ball: %.2f\n", qm_size_av);
  fprintf(stdout, "Average metric radius: %.2f\n", delta_av);
  fprintf(stdout, "Average size of m-ball: %.2f\n", m_size_av);
  fprintf(stdout, "Average ratio m / qm : %.2f\n", ratio_av);
 
 return EXIT_SUCCESS;
}
