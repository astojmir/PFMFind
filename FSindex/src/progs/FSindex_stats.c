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
  FSINDX *FSI;

  if (argc < 2)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file \n",
	      argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  fprintf(stdout, "Loading index... \n");
  FSI = FS_INDEX_load(filename);
  FS_INDEX_print_stats(FSI, stdout, 0, 0); 

 return EXIT_SUCCESS;
}
