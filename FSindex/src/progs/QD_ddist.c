#include "misclib.h"
#include "fastadb.h"
#include "blastlib.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "randseq.h"
#include "hit_list.h"
#include "partition.h"
#include "smatrix.h"

int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;

int main(int argc, char **argv)
{
  /* HISTOGRAMS */
  SEED_HIST QD_left_hist; 
  SEED_HIST QD_right_hist;
  SEED_HIST MD_hist;
  
  /* DATABASE */
  SEQUENCE_DB *s_db;
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  ULINT no_frags;

  /* RANDOM SEQUENCE GENERATOR */
  SEQ_GENERATOR_t *seq_generator;
  BIOSEQ query;
  BIOSEQ subject;
  char *seq_heap;

  /* SCORING MATRICES */
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *MD;
  SCORE_MATRIX_t *QMD;


  /* OTHERS */
  const char *matrix_file;
  const char *count_file;
  const char *database_file;
  ULINT no_runs;
  ULINT no_samples;
  ULINT frag_len;
  ULINT i, j;
  int value;
  ULINT rand_frag;
  FS_SEQ_t FS_seq;


  /* Arguments:
     - matrix file
     - database file
     - count file
     - fragment length
     - number of runs (random sequences)
     - number of samples from database
  */

  if (argc < 7)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s matrix_file database_file count_file"
	      " frag_length no_runs no_samples"
	      "\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  matrix_file = argv[1];
  database_file = argv[2];
  count_file = argv[3];
  frag_len = atoi(argv[4]);
  no_runs = atoi(argv[5]);
  no_samples = atoi(argv[6]);

  /* Load database */
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = NONE;
  fastadb_argv[0].access_type = MEMORY;
  fastadb_argv[1].retrieve_deflines = NO;
  s_db = fastadb_open(database_file, fastadb_argt, fastadb_argv); 

  /* Initialise fragment database - we use old code! */
  fastadb_init_frags(s_db, frag_len, frag_len);
  fastadb_get_nofrags(s_db, &no_frags, frag_len, frag_len);

  /* Load matrices */
  ptable = FS_PARTITION_create("STAN#ILVM#KR#EDQ#WFYH#GPC", '#'); 
  S = SCORE_MATRIX_create(matrix_file, ptable);
  QMD = SCORE_MATRIX_S_2_Dquasi(S);
  MD = SCORE_MATRIX_S_2_Dmax(S);

  /* Initialise the random sequence generator */
  seq_generator = SEQ_GENERATOR_create(count_file, ptable);
  seq_heap = mallocec(frag_len+1);


  /* Main loop */
  fprintf(stdout,"Number of runs: %ld\n", no_runs);

  for (i=0; i < no_runs; i++)
    {
      /* Initialise histograms */
      memset(&QD_left_hist, 0, sizeof(SEED_HIST));
      memset(&QD_right_hist, 0, sizeof(SEED_HIST));
      memset(&MD_hist, 0, sizeof(SEED_HIST));
      seedhist_init1(&QD_left_hist, 0, 50); 
      seedhist_init1(&QD_right_hist, 0, 50); 
      seedhist_init1(&MD_hist, 0, 50);
 
      /* Generate random sequence */
      SEQ_GENERATOR_rand_seq(seq_generator, &query, frag_len,
			     seq_heap);       
      for (j=0; j < no_samples; j++)
	{
	  /* Get a random database fragment */
	  do 
	    {
	      rand_frag = lrand48()%no_frags;
	      fastadb_get_frag(s_db, &subject, rand_frag,
			       frag_len, frag_len);
	    }
	  while (!BIOSEQ_2_FS_SEQ(&subject, ptable, &FS_seq));

	  /* Calculate distances and add to histograms*/
	  value = SCORE_MATRIX_evaluate(QMD, &query, &subject);
	  seedhist_add_acount(&QD_left_hist, value);
	  value = SCORE_MATRIX_evaluate(QMD, &subject, &query);
	  seedhist_add_acount(&QD_right_hist, value);
	  value = SCORE_MATRIX_evaluate(MD, &query, &subject);
	  seedhist_add_acount(&MD_hist, value);
#if DEBUG > 1
	  if (value == 0)
	    fprintf(stderr,"%.*s %ld\n", (int) frag_len, 
		    subject.start, rand_frag);
#endif
	}
      /* Print histograms */
      fprintf(stdout,"RANDOM POINT: %.*s\n", (int) frag_len, 
	      query.start); 
      fprintf(stdout, "LEFT_QUASI_METRIC_DISTANCES\n");
      seedhist_print(&QD_left_hist, stdout);
      fprintf(stdout, "RIGHT_QUASI_METRIC_DISTANCES\n");
      seedhist_print(&QD_right_hist, stdout);
      fprintf(stdout, "METRIC_DISTANCES\n");
      seedhist_print(&MD_hist, stdout);

      /* Clear histograms */
      seedhist_clear1(&QD_left_hist);
      seedhist_clear1(&QD_right_hist);
      seedhist_clear1(&MD_hist);    
   }
  return EXIT_SUCCESS;
}
