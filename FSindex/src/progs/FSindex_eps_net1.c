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

  ULINT no_MDneighbours = 0;
  double P_MDneighbours = 0.0;
  ULINT no_QMDneighbours = 0;
  double P_QMDneighbours = 0.0;
  double ratio;
  ULINT D_cutoff = 0;
  FILE *stream = stdout;
  ULINT one_percent;
  ULINT step_size;
  SCORE_MATRIX_t *MD;
  SCORE_MATRIX_t *QMD;
  char *matrix_full;
  char *matrix_base;
  char *matrix_dir;
  int no_args = 6;
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix_file count_file"
	      " no_runs start_D step_size <out_file>\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  index_file = argv[1];
  matrix_full = argv[2]; 
  count_file = argv[3];
  no_runs = atoi(argv[4]);
  D_cutoff = atoi(argv[5]);
  step_size = atoi(argv[6]);
  if (argc >= no_args+2)
    {
      out_file = argv[7];
      if ((stream = fopen(out_file, "w")) == NULL)
	{
	  fprintf(stderr, "Unable to open %s\n ", out_file);
	  exit(EXIT_FAILURE);
	} 
    } 

  FS_index = FS_INDEX_load(index_file);

  ptable = FS_INDEX_get_ptable(FS_index);
  S = SCORE_MATRIX_create(matrix_full, ptable); 
  QMD = SCORE_MATRIX_S_2_Dquasi(S);
  MD = SCORE_MATRIX_S_2_Dmax(S);
  split_base_dir(matrix_full, &matrix_base, &matrix_dir);  



  seq_generator = SEQ_GENERATOR_create(count_file, ptable);
  len = FS_index->frag_len;
  seq_heap = mallocec(len+1);

  one_percent = no_runs/100;
  fprintf(stream, "***** Epsilon neighbourhood of the %s"
	  " database ******\n", FS_index->db_name);
  fprintf(stream, "No runs: %ld\n", no_runs);
  fprintf(stream, "%8.8s %15.15s %15.15s %15.15s %15.15s %15.15s\n",\
	  "D_cutoff", "QD Neighbours", "QD Proportion",
	  "MD Neighbours", "MD Proportion", "QD/MD Ratio");   
  if (D_cutoff == 0)
    {
      fprintf(stream, "%8ld %15ld %15.4f %15ld %15.4f\n", 
	      D_cutoff, no_QMDneighbours, P_QMDneighbours,
	      no_MDneighbours, P_MDneighbours);        
      fflush(stream);
      D_cutoff++;
    }
  while (no_MDneighbours < no_runs || no_QMDneighbours < no_runs)
    {
      fprintf(stdout,"D = %ld\n", D_cutoff);
      no_QMDneighbours = 0;
      no_MDneighbours = 0;
      for (i = 0; i < no_runs; i++)
	{
	  SEQ_GENERATOR_rand_seq(seq_generator, &query, len,
				 seq_heap); 
	  FS_INDEX_set_matrix(FS_index, matrix_base, S, QMD);
	  if (FS_INDEX_has_neighbour(FS_index, &query, D_cutoff))
	    no_QMDneighbours++;
	  FS_INDEX_set_matrix(FS_index, matrix_base, S, MD);
	  if (FS_INDEX_has_neighbour(FS_index, &query, D_cutoff))
	    no_MDneighbours++;
	  /* Print progress bar */
	  printbar(stdout, i+1, one_percent, 50);  
	}
      P_QMDneighbours = (double) no_QMDneighbours / (double) no_runs;
      P_MDneighbours = (double) no_MDneighbours / (double) no_runs;
      fprintf(stdout, "\n");
      if (no_MDneighbours > 0)
	{
	  ratio = P_QMDneighbours / P_MDneighbours;
	  fprintf(stream, "%8ld %15ld %15.4f %15ld %15.4f %15.4f\n", 
		  D_cutoff, no_QMDneighbours, P_QMDneighbours,
		  no_MDneighbours, P_MDneighbours, ratio);
	}
      else
	{
	  fprintf(stream, "%8ld %15ld %15.4f %15ld %15.4f\n", 
		  D_cutoff, no_QMDneighbours, P_QMDneighbours,
		  no_MDneighbours, P_MDneighbours);
	}
      fflush(stream);
      D_cutoff += step_size;
    }
  if (out_file != NULL)
    fclose(stream);
  return EXIT_SUCCESS;
}
