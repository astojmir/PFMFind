#include "misclib.h"
#include "fastadb.h"
#include "blastlib.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "randseq.h"
#include "hit_list.h"
#include "partition.h"
#include "smatrix.h"
#include "FSindex.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;
int POS_MATRIX_VERBOSE = 0;
FILE *POS_MATRIX_STREAM;



static
double Zscore(BIOSEQ *query, BIOSEQ *subject, double value, 
	      SCORE_MATRIX_t *M)
{
  int j;
  int k;
  int r;
  int val;
  char tmp;
  double x;
  double xx;
  double mean;
  double sd;
  int p = 100;
  
  srand48(time(NULL)); 

  x=0.0;
  xx=0.0;     
  for (j=0; j < p; j++)
    {
      for (k=subject->len-1; k > 0; k--)
	{
	  r = lrand48() % subject->len;
	  tmp = subject->start[r];
	  subject->start[r] = subject->start[k];
	  subject->start[k] = tmp;
	}
      val = (double) SCORE_MATRIX_evaluate(M, query, subject);
      x += val;
      xx += val*val;
    }
  mean = x/p;
  sd = sqrt(xx/p - mean*mean);
  return (value - mean)/sd;
}




int main(int argc, char **argv)
{
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
  SCORE_MATRIX_t *QMD;


  /* OTHERS */
  const char *matrix_file;
  const char *count_file;
  const char *database_file;
  ULINT no_samples;
  ULINT frag_len;
  ULINT j;
  double value;
  double Z;
  ULINT rand_frag;
  FS_SEQ_t FS_seq;

  /* Arguments:
     - matrix file
     - database file
     - count file
     - fragment length
     - number of samples from database
  */

  if (argc < 6)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s matrix_file database_file count_file"
	      " frag_length no_samples"
	      "\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  matrix_file = argv[1];
  database_file = argv[2];
  count_file = argv[3];
  frag_len = atoi(argv[4]);
  no_samples = atoi(argv[5]);

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


  /* Initialise the random sequence generator */
  seq_generator = SEQ_GENERATOR_create(count_file, ptable);
  seq_heap = mallocec(frag_len+1);

  /* Generate random sequence */
  SEQ_GENERATOR_rand_seq(seq_generator, &query, frag_len,
			 seq_heap);       


  
  printf("Query: %.*s\n", (int) frag_len, query.start);

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
      value = (double) SCORE_MATRIX_evaluate(QMD, &query, &subject);

      /* Get Z-Score */
      Z = -Zscore(&query, &subject, value, QMD);

      /* Print */
      printf("%.4f\n", Z);
    }
  return EXIT_SUCCESS;
}
