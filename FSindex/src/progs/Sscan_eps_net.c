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
  BIOSEQ query;
  const char *matrix_file;
  const char *count_file;
  const char *database_file;
  SEQ_GENERATOR_t *seq_generator;
  ULINT no_runs;
  ULINT i;
  ULINT len;
  char *seq_heap;
  const char *out_file = NULL;

  ULINT no_neighbours = 0;
  double P_neighbours = 0.0;
  ULINT D_cutoff = 0;
  FILE *stream = stdout;
  ULINT one_percent;
  ULINT step_size;
  SEQUENCE_DB *s_db;
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];

  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;


  if (argc < 8)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s matrix_file database_file count_file"
	      " frag_length no_runs start_D step_size"
	      " <out_file>\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  matrix_file = argv[1];
  database_file = argv[2];
  count_file = argv[3];
  len = atoi(argv[4]);
  no_runs = atoi(argv[5]);
  D_cutoff = atoi(argv[6]);
  step_size = atoi(argv[7]);
  if (argc >= 9)
    {
      out_file = argv[8];
      if ((stream = fopen(out_file, "w")) == NULL)
	{
	  fprintf(stderr, "Unable to open %s\n ", out_file);
	  exit(EXIT_FAILURE);
	} 
    } 
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = NONE;
  fastadb_argv[0].access_type = MEMORY;
  fastadb_argv[1].retrieve_deflines = NO;

  s_db = fastadb_open(database_file, fastadb_argt, fastadb_argv); 

  ptable = FS_PARTITION_create("STAN#ILVM#KR#EDQ#WFYH#GPC", '#'); 
  S = SCORE_MATRIX_create(matrix_file, ptable);
  D = SCORE_MATRIX_S_2_Dquasi(S);

  seq_generator = SEQ_GENERATOR_create(count_file, ptable);
  len = 10;
  seq_heap = mallocec(len+1);

  one_percent = no_runs/100;
  fprintf(stream, "***** Epsilon neighbourhood of the %s"
	  " database ******\n", database_file);
  fprintf(stream, "No runs: %ld\n", no_runs);
  fprintf(stream, "%10.10s %10.10s %10.10s\n", "D_cutoff",
	  "Neighbours", "Proportion"); 
  if (D_cutoff == 0)
    {
      fprintf(stream, "%10ld %10ld %10.4f\n", D_cutoff, no_neighbours,
	      P_neighbours);        
      fflush(stream);
      D_cutoff++;
    }
  while (no_neighbours < no_runs)
    {
      no_neighbours = 0;
      for (i = 0; i < no_runs; i++)
	{
	  SEQ_GENERATOR_rand_seq(seq_generator, &query, len,
				 seq_heap); 
	  if (SSCAN_has_neighbour(s_db, D, ptable, &query, D_cutoff)) 
	    no_neighbours++;
	  /* Print progress bar */
	  printbar(stdout, i+1, one_percent, 50);  
	}
      P_neighbours = (double) no_neighbours / (double) no_runs;
      fprintf(stdout, "\n");
      fprintf(stream, "%10ld %10ld %10.4f\n", D_cutoff, no_neighbours,
	      P_neighbours);       
      fflush(stream);
      D_cutoff += step_size;
    }
  if (out_file != NULL)
    fclose(stream);
  return EXIT_SUCCESS;
}
