#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "randseq.h"
#include "FSindex.h"
#include "hit_list.h"


EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;

#define KNN_SIZE 6

int main(int argc, char **argv)
{
  BIOSEQ query;
  const char *index_file;
  const char *count_file;
  SEQ_GENERATOR_t *seq_generator;
  ULINT no_runs;
  ULINT i;
  ULINT j;
  int len;
  char *seq_heap;

  char *matrix_full;

  int no_args = 5;

  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  HIT_LIST_t *HL = NULL;
  FSINDX *FSI;
  SEQUENCE_DB *sdb;


  int KNN_val[KNN_SIZE] = {1, 50, 100, 250, 500, 1000};
  int KNN_eps;
  int KNN_size=0;

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix_file"
	      " count_file no_runs NN ..\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  index_file = argv[1];
  matrix_full = argv[2]; 
  count_file = argv[3];
  no_runs = atoi(argv[4]);

  for (i=5; i < argc && KNN_size < KNN_SIZE; i++) {
    KNN_val[KNN_size++] = atoi(argv[i]) ;
  }


  Try {
    /* Initialise everything */
    FSI = FS_INDEX_load(index_file);
    ptable = FS_INDEX_get_ptable(FSI);
    S = SCORE_MATRIX_create(matrix_full, ptable); 
    D = SCORE_MATRIX_S_2_Dquasi(S);
    seq_generator = SEQ_GENERATOR_create(count_file, ptable);
    len = FS_INDEX_get_frag_len(FSI);
    seq_heap = mallocec(len+1);
    
    /* Header */
    sdb = FS_INDEX_get_database(FSI);
    printf("Dataset: %s\n", sdb->db_name);
    printf("Matrix: %s\n", matrix_full);
    printf("Count file: %s\n", count_file);
    printf("Fragment length: %d\n", len);
    printf("Random sequences: %ld\n", no_runs); 
    printf("kNN instances: %d\n", KNN_size); 
    printf("%*.*s %s\n", len, len, "Sequence",
	   "          ... Range ...");
    printf("%*.*s ", len, len, "");
    for (j=0; j < KNN_size; j++)
      printf("%6d ", KNN_val[j]);
    printf("\n");
    
    /* All runs */
    for (i = 0; i < no_runs; i++)
      {
	/* generate a random sequence */
	SEQ_GENERATOR_rand_seq(seq_generator, &query, len,
			       seq_heap); 
	printf("%*.*s ", len, len, query.start);
	
	/* Run kNN search to find epsilon for all values of K supplied */
	for (j=0; j < KNN_size; j++)
	  {
	    HL = FSINDX_kNN_srch(FSI, &query, D, KNN_val[j], HL);
	    KNN_eps = (int) HL->range;
	    printf("%6d ", KNN_eps);
	  }
	printf("\n");
      }

    fprintf(stdout,"\n");
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
