#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "hit_list.h"
#include "partition.h"
#include "smatrix.h"

int FS_PARTITION_VERBOSE = 1;
int SCORE_MATRIX_VERBOSE = 1;

int main(int argc, char **argv)
{
  const char *alphabet;
  const char *matrix;
  const char separator = '#';
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *MD;
  SCORE_MATRIX_t *QMD;

  if (argc < 3)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s matrix clusters\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  matrix = argv[1];
  alphabet = argv[2];

  ptable = FS_PARTITION_create(alphabet, separator); 
  S = SCORE_MATRIX_create(matrix, ptable); 
  QMD = SCORE_MATRIX_S_2_Dquasi(S);
  MD = SCORE_MATRIX_S_2_Dmax(S);

  return EXIT_SUCCESS;
}
