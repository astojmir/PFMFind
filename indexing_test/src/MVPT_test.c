#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "fastadb.h"
#include "FSindex.h"
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
  UFRAG_DB *udb;
  SCORE_MATRIX *D;
  HIT_LIST_t *HL1 = NULL;
  HIT_LIST_t *HL2 = NULL;
  MVP_TREE *MVPT2 = NULL;

  int i;
  int no_tests;

  char *mvpt_file;
  char *FSindex_file;
  char *query_file;
  int radius= 15;
  int kNN = 1;

  int no_args = 4;
  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s mvpt_file FSindex_file query_file kNN [out_file]\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  mvpt_file = argv[1];
  FSindex_file = argv[2];
  query_file = argv[3];
  kNN = atoi(argv[4]);


  Try {
    printf ("Loading MVPT2\n");
    MVPT2 = MVP_TREE_load2(mvpt_file);
    MVP_TREE_fprint_stats(MVPT2, stdout, 0);
    D = (SCORE_MATRIX *) MVPT2->matrix;
    SCORE_MATRIX_fprint(D, stdout);
    udb = (UFRAG_DB *) MVPT2->db;

#if 0
    printf ("Loading FSindex\n");
    I = FS_INDEX_load(FSindex_file);
#endif
    printf ("Loading queries\n");
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


    printf ("\n**** Running Tests ****\n");
    for (i=0; i < no_tests; i++) {
      query.start = queries[i];
      query.len = udb->frag_len;

#if 0
      /* FSindex run */
      HL1 = FSINDX_kNN_srch(I, &query, D, kNN, HL1, FS_BINS, SARRAY);
      HIT_LIST_sort_decr(HL1);
      radius = HL1->hits[0].dist;
#endif
      /* MVPT run */
      HL2 = MVP_TREE_rng_srch2(MVPT2, &query, radius, HL2);
      HIT_LIST_sort_by_sequence(HL2);
       
      printf("%.3d %s %5ld %8ld (%ld%%)\n", i, queries[i], 
	     HL2->actual_seqs_hits,
	     HL2->seqs_visited, 
	     HL2->seqs_visited*100/MVPT2->no_pts);
#if 0
      /* Comparison */

      test_passed = 1;
      if (HL1->actual_seqs_hits != HL2->actual_seqs_hits)
	test_passed = 0;
      else
	for (j=0; j<HL1->actual_seqs_hits; j++)
	  if (SEQ_HIT_cmp_by_seq(HL1->hits+j, HL2->hits+j) != 0)
	    test_passed = 0;

      if (test_passed) 
	printf("%.3d %s %5ld %5ld %8ld %8ld (%ld%%) PASSED\n", i, queries[i], 
	       HL1->actual_seqs_hits,
	       HL2->actual_seqs_hits, HL1->seqs_visited, HL2->seqs_visited,
	       HL2->seqs_visited*100/MVPT->no_pts);
      else 
	printf("%.3d %s %5ld %5ld %8ld %8ld (%ld%%) FAILED\n", i, queries[i], 
	       HL1->actual_seqs_hits,
	       HL2->actual_seqs_hits, HL1->seqs_visited, HL2->seqs_visited,
	       HL2->seqs_visited*100/MVPT->no_pts);
#endif
#if 0
      j=0;
      k=0;
      while (j<HL1->actual_seqs_hits && k<HL2->actual_seqs_hits) {
	if (j == HL1->actual_seqs_hits) {
	  printf("- %*.*s %7ld %7ld %7d\n", (int) udb->frag_len, (int) udb->frag_len,
		 HL2->hits[k].subject.start, HL2->hits[k].sequence_id,
		 HL2->hits[k].sequence_from, 
		 fastadb_data_offset(udb->sdb, HL2->hits[k].subject.start));
	  k++;
	}
	else if (k == HL2->actual_seqs_hits) {
	  bad_offsets[b++] = 
	    fastadb_data_offset(udb->sdb, udb->sdb->seq[HL1->hits[j].sequence_id].start) 
	    + HL1->hits[j].sequence_from;
	  
	  printf("+ %*.*s %7ld %7ld %7d\n", (int) udb->frag_len, (int) udb->frag_len,
		 HL1->hits[j].subject.start, HL1->hits[j].sequence_id,
		 HL1->hits[j].sequence_from, bad_offsets[b-1]);
	  j++;
	}
	else if (SEQ_HIT_cmp_by_seq(HL1->hits+j, HL2->hits+k) < 0) {
	  bad_offsets[b++] = 
	    fastadb_data_offset(udb->sdb, udb->sdb->seq[HL1->hits[j].sequence_id].start) 
	    + HL1->hits[j].sequence_from;
	  
	  printf("+ %*.*s %7ld %7ld %7d\n", (int) udb->frag_len, (int) udb->frag_len,
		 HL1->hits[j].subject.start, HL1->hits[j].sequence_id,
		 HL1->hits[j].sequence_from, bad_offsets[b-1]);
	  j++;
	}
	else if (SEQ_HIT_cmp_by_seq(HL1->hits+j, HL2->hits+k) > 0) {
	  printf("- %*.*s %7ld %7ld %7d\n", (int) udb->frag_len, (int) udb->frag_len,
		 HL2->hits[k].subject.start, HL2->hits[k].sequence_id,
		 HL2->hits[k].sequence_from, 
		 fastadb_data_offset(udb->sdb, HL2->hits[k].subject.start));
	  k++;
	}
	else {
	  j++;
	  k++;
	}
	
      }

      printf ("\nSearching missing offsets\n");
      for (j=0; j < udb->pts_size; j++) {
	for (k=0; k < b; k++) 
	  if (bad_offsets[k] == udb->pts[j]) {
	    printf("%ld Found in pts!\n", bad_offsets[k]); 
	    bad_oid[k] = j;
	  }
      }
      for (j=0; j < udb->dup_heap_size; j++) {
	for (k=0; k < b; k++) 
	  if (bad_offsets[k] == udb->dup_heap[j]) {
	    printf("%ld Found in dup_heap!\n", bad_offsets[k]);
	    bad_oid[k] = j;
	  } 
      }

      printf ("\nSearching in MVPT\n");
      for (j=0; j < MVPT->nodes_used; j++) {
	for (k=0; k < b; k++) { 
	  if (bad_oid[k] == MVPT->nodes[j].vp1)
	    printf("%ld Found as vp1 in node %d!\n", bad_oid[k], j);
	  else if (bad_oid[k] == MVPT->nodes[j].vp2)
	    printf("%ld Found as vp2 in node %d!\n", bad_oid[k], j);
	  
	  if (MVPT->nodes[j].node_type == LEAF) {
	    for (l=0; l < MVPT->nodes[j].size - 1; l++)
	      if (bad_oid[k] == MVPT->nodes[j].CD.leaf.datapoint[l])
		printf("%ld Found as datapoint %d in leaf node %d!\n", bad_oid[k], l, j);
	    
	  }
	  
	}
      }

      printf ("\n******* FSINDEX list *******\n");
      HIT_LIST_print(HL1, stdout, 0);
      printf ("\n******* MVPT list *******\n");
      HIT_LIST_print(HL2, stdout, 0);
      printf ("\n\n");
#endif
    }
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

