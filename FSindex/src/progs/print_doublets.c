#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include "FSindex.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;

int main(int argc, char **argv)
{
  const char *filename;
  int i;
  FS_HASH_TABLE_t *HT;
  FS_PARTITION_t *ptable;
  if (argc < 2)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file \n",
	      argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  fprintf(stdout, "Loading index... \n");
  FS_INDEX_load(filename);
  HT = FS_INDEX_get_hash_table();
  ptable = FS_INDEX_get_ptable(); 
  printf("Non-overlapping doublet frequencies\n");
  for (i=0; i < HT->no_bins; i++)
    {
      printf("%s %ld %6.4f\n", FS_seq_print(i, ptable, 2), 
	     HT->bin_size[i],
	     (double) HT->bin_size[i]/ HT->no_seqs);
    }
  printf("\n");
  return EXIT_SUCCESS;
}
