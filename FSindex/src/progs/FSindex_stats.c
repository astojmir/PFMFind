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
  FS_INDEX_t *FS_index;
  const char *filename;

  if (argc < 2)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file \n",
	      argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  FS_index = FS_INDEX_load(filename);
  FS_PARTITION_print(FS_index->ptable, stdout);
  FS_INDEX_print_stats(FS_index, stdout, 0, 0); 

 return EXIT_SUCCESS;
}
