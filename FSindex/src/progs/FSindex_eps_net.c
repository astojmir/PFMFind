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
  FS_INDEX_t *FS_index;
  const char *index_file;
  const char *count_file;
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

  if (argc < 4)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file count_file"
	      " no_runs <out_file>\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  index_file = argv[1];
  count_file = argv[2];
  no_runs = atoi(argv[3]);
  if (argc >= 5)
    {
      out_file = argv[4];
      if ((stream = fopen(out_file, "w")) == NULL)
	{
	  fprintf(stderr, "Unable to open %s\n ", out_file);
	  exit(EXIT_FAILURE);
	} 
    } 

  FS_index = FS_INDEX_load(index_file);
  seq_generator = SEQ_GENERATOR_create(count_file, FS_index->ptable);
  len = FS_index->frag_len;
  seq_heap = mallocec(len+1);

  one_percent = no_runs/100;
  fprintf(stream, "***** Epsilon neighbourhood of the %s"
	  " database ******\n", FS_index->db_name);
  fprintf(stream, "No runs: %ld\n", no_runs);
  fprintf(stream, "%10.10s %10.10s %10.10s\n", "D_cutoff",
	  "Neighbours", "Proportion");  
  fprintf(stream, "%10ld %10ld %10.4f\n", D_cutoff, no_neighbours,
	  P_neighbours);        
  fflush(stream);
  D_cutoff++;
  while (no_neighbours < no_runs)
    {
      no_neighbours = 0;
      for (i = 0; i < no_runs; i++)
	{
	  SEQ_GENERATOR_rand_seq(seq_generator, &query, len,
				 seq_heap); 
	  if (FS_INDEX_has_neighbour(FS_index, &query, D_cutoff))
	    no_neighbours++;
	  /* Print progress bar */
	  printbar(stdout, i+1, one_percent, 50);  
	}
      P_neighbours = (double) no_neighbours / (double) no_runs;
      fprintf(stdout, "\n");
      fprintf(stream, "%10ld %10ld %10.4f\n", D_cutoff, no_neighbours,
	      P_neighbours);       
      fflush(stream);
      D_cutoff++;
    }
  if (out_file != NULL)
    fclose(stream);
  return EXIT_SUCCESS;
}
