#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "fastadb.h"
#include "FSindex.h"
#include "SAindex.h"
#include "hit_list.h"
#include "misclib.h"
#include "mvptree.h"

#define MAX_RUNS 5000
#define MAX_LINES 10
#define BUF_SIZE 1023

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int main(int argc, char **argv)
{
  BIOSEQ query;

  FILE *fp;
  char buffer[BUF_SIZE+1];
  char **queries = NULL;

  FSINDX *I;
  SAINDX *A;
  UFRAG_DB *udb;
  SCORE_MATRIX *D;
  HIT_LIST_t *HL = NULL;
  MVP_TREE *MVPT2 = NULL;

  int i;
  int no_tests;

  char *mvpt_file;
  char *FSindex_file;
  char *SAindex_file;
  char *query_file;
  int radius= 15;
  int kNN = 1;

  int no_args = 5;
  if (argc < no_args+1) {
    fprintf(stderr,"Insufficient arguments \n");
    fprintf(stderr,"Usage: %s mvpt_file FSindex_file SAindex_file query_file kNN\n", argv[0]);  
    exit(EXIT_FAILURE);
  }

  mvpt_file = argv[1];
  FSindex_file = argv[2];
  SAindex_file = argv[3];
  query_file = argv[4];
  kNN = atoi(argv[5]);


  Try {
    MVPT2 = MVP_TREE_load2(mvpt_file);
    // MVP_TREE_fprint_stats(MVPT2, stdout, 0);
    D = (SCORE_MATRIX *) MVPT2->matrix;
    // SCORE_MATRIX_fprint(D, stdout);
    udb = (UFRAG_DB *) MVPT2->db;

    I = FS_INDEX_load(FSindex_file);
    A = SA_INDEX_load(SAindex_file);

    queries = mallocec(MAX_RUNS * sizeof(char *));
    i=0;
    if((fp = fopen(query_file, "r")) == NULL)
      Throw FSexcept(FOPEN_ERR, "Could not open the file %s.",
		     query_file);
    
    while (fgets(buffer, BUF_SIZE, fp)) {
      if (buffer[0] == '#') continue;
      if (i >= MAX_RUNS) break;
      buffer[strlen(buffer)-1] = '\0'; /* Removes '\n' */
      
      /* We assume we only have sequences on each line -
	 use original query files */
      queries[i] = strdup(buffer);
      i++;
    }
    no_tests = i;
    fclose(fp);


    
    printf("%-*.*s %5.5s %5.5s %9.9s %9.9s %9.9s %9.9s %12.12s %12.12s ",
	   (int) udb->frag_len, (int) udb->frag_len,
	   "# FRAG", "RANGE",
	   "HITS1", "BINSV1", "BINSH1",
	   "SEQSV1", "SEQSH1",
	   "CHARV1", "CHARH1"); 

    printf("%5.5s %9.9s %9.9s %9.9s %9.9s %12.12s %12.12s ",
	   "HITS2", "BINSV2", "BINSH2",
	   "SEQSV2", "SEQSH2",
	   "CHARV2", "CHARH2"); 

    printf("%5.5s %9.9s %9.9s %9.9s %9.9s %12.12s %12.12s ",
	   "HITS3", "BINSV3", "BINSH3",
	   "SEQSV3", "SEQSH3",
	   "CHARV3", "CHARH3"); 

    printf("%5.5s %9.9s %9.9s\n",
	   "HITS4", 
	   "SEQSV4", "SEQSH4");

    for (i=0; i < no_tests; i++) {
      query.start = queries[i];
      query.len = udb->frag_len;

      /* FSindex kNN run */
      HL = FSINDX_kNN_srch(I, &query, D, kNN, HL, FS_BINS, SARRAY);
      HIT_LIST_sort_decr(HL);
      radius = HL->hits[0].dist;

      printf("%s %5d %5ld %9ld %9ld %9ld %9ld %12ld %12ld ",
	     queries[i], radius, HL->actual_seqs_hits,
	     HL->FS_seqs_visited, HL->FS_seqs_hits,
	     HL->seqs_visited, HL->seqs_hits,
	     HL->useqs_visited, HL->useqs_hits);

      /* FSindex range run */
      HL = FSINDX_rng_srch(I, &query, D, radius, MAX, HL, FS_BINS, SARRAY);
      printf("%5ld %9ld %9ld %9ld %9ld %12ld %12ld ",
	     HL->actual_seqs_hits,
	     HL->FS_seqs_visited, HL->FS_seqs_hits,
	     HL->seqs_visited, HL->seqs_hits,
	     HL->useqs_visited, HL->useqs_hits);

      /* FSindex suffix array run */
      HL = SAINDX_rng_srch(A, &query, D, radius, MAX, HL);
      printf("%5ld %9ld %9ld %9ld %9ld %12ld %12ld ",
	     HL->actual_seqs_hits,
	     HL->FS_seqs_visited, HL->FS_seqs_hits,
	     HL->seqs_visited, HL->seqs_hits,
	     HL->useqs_visited, HL->useqs_hits);

      /* MVPT run */
      HL = MVP_TREE_rng_srch2(MVPT2, &query, radius, HL);
      printf("%5ld %9ld %9ld\n",
	     HL->actual_seqs_hits,
	     HL->seqs_visited, HL->seqs_hits);
    }
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}


