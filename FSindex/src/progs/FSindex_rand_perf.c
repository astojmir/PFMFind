#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
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
#define BUF_SIZE 1023
int main(int argc, char **argv)
{
  BIOSEQ query;
  const char *index_file;
  const char *query_file;
  FILE *fp;
  int n=0;
  int no_runs;
  int i;
  int len;
  int frag_len;

  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  HIT_LIST_t *HL = NULL;
  FSINDX *FSI;
  SEQUENCE_DB *sdb;
  char buffer[BUF_SIZE+1];
  int K_item = 0;
  int K = 0;
  int K_items = 0;
  char *matrix_full;
  int no_args = 5;
  int eps;

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix_file query_file"
	      " runs K_item\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  index_file = argv[1];
  matrix_full = argv[2]; 
  query_file = argv[3];
  no_runs = atoi(argv[4]);
  K_item = atoi(argv[5])-1;


  Try {
    /* Open file, read header */
    if((fp = fopen(query_file, "r")) == NULL)
      Throw FSexcept(FOPEN_ERR, "Could not open the file %s.",
		     query_file);

    /* Read header - 7 lines */
    fgets(buffer, BUF_SIZE, fp); /* Dataset - skip                  */
    fgets(buffer, BUF_SIZE, fp); /* Matrix - skip                   */
    fgets(buffer, BUF_SIZE, fp); /* Count file - skip               */
    fgets(buffer, BUF_SIZE, fp); /* Fragment length                 */
    sscanf(buffer+16, "%d", &frag_len);
    fgets(buffer, BUF_SIZE, fp); /* Maximum number of runs          */
    sscanf(buffer+17, "%d", &n);
    fgets(buffer, BUF_SIZE, fp); /* Number of kNN instances         */
    sscanf(buffer+15, "%d", &K_items);
    if (K_item >= K_items)
      Throw FSexcept(INDEX_OUT_OF_RANGE, "K_items out of bounds.\n");
    fgets(buffer, BUF_SIZE, fp); /* Column names - skip             */
    fgets(buffer, BUF_SIZE, fp); /* Numbers of nearset neighbours   */
    sscanf(buffer+frag_len+1+7*K_item,"%d", &K);


    /* Initialise everything */
    FSI = FS_INDEX_load(index_file);
    ptable = FS_INDEX_get_ptable(FSI);
    S = SCORE_MATRIX_create(matrix_full, ptable); 
    D = SCORE_MATRIX_S_2_Dquasi(S);
    len = FS_INDEX_get_frag_len(FSI);

    /* Check consistency */
    if (n > no_runs)
      n = no_runs;
    if (len != frag_len)
      Throw FSexcept(INCONSITENT_DATA, "Fragment lengths do not match.\n");

    query.len = len;
  
    /* Header */
    sdb = FS_INDEX_get_database(FSI);
    printf("***** FSsearch random runs: %s"
	   " database ******\n", sdb->db_name);
    printf("No runs: %d\n", n);
    printf("Fragment length: %d\n", len);
    printf("Total fragments: %ld\n", FSI->no_seqs);
    printf("Total unique fragments: %ld\n", FSI->no_useqs);
    printf("Total bins: %ld\n", FS_INDEX_get_no_bins(FSI));
    printf("Clustering: %s\n", FSI->alphabet);
    printf("Matrix: %s\n", matrix_full);
    printf("\n");
    printf("%6.6s %*.*s %6.6s %6.6s %10.10s %10.10s %10.10s"
	   " %10.10s %10.10s %10.10s %10.10s %10.10s %10.10s\n",
	   "run", (int) len , (int) len, "query", "k", "eps", 
	   "FS_visited", "FS_hits", "FS_ratio", "s_visited", 
	   "s_hits", "s_ratio", "u_visited", "u_hits", "u_ratio");   
  
     /* All runs */
    for (i = 0; i < n; i++)
      {
	/* Read the line */
	fgets(buffer, BUF_SIZE, fp);
	query.start = buffer;
	sscanf(buffer+len+1+7*K_item,"%d", &eps);

	HL = FSINDX_rng_srch(FSI, &query, D, eps, HL);
      
	printf("%6.6d %*.*s %6d %6d %10ld %10ld %10.2f %10ld"
	       " %10ld %10.2f %10ld %10ld %10.2f\n", i, 
	       (int) len , (int) len, query.start,
	       K, HL->range, HL->FS_seqs_visited,
	       HL->FS_seqs_hits, 
	       (double) HL->FS_seqs_visited / HL->FS_seqs_hits, 
	       HL->seqs_visited, HL->seqs_hits,
	       (double) HL->seqs_visited / HL->seqs_hits,
	       HL->useqs_visited, HL->useqs_hits,
	       (double) HL->useqs_visited / HL->useqs_hits);
      }
    fclose(fp);
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
