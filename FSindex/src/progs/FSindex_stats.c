#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include "FSindex.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;
int POS_MATRIX_VERBOSE = 0;
FILE *POS_MATRIX_STREAM;

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
  Try {
    FSI = FS_INDEX_load(filename);
    FS_INDEX_print_stats(FSI, stdout, 2); 
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }
 return EXIT_SUCCESS;
}
