#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "FSindex.h"
#include "hit_list.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;

int main(int argc, char **argv)
{
  const char *query_seq;
  int cutoff;
  HIT_LIST_t *hit_list;
  BIOSEQ query;
  FS_INDEX_t *FS_index;
  const char *filename;
  ULINT no_seq;
  ULINT i;

  if (argc < 4)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file no_seq"
	      " cutoff\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  no_seq = atoi(argv[2]);
  cutoff = atoi(argv[3]);

  FS_index = FS_INDEX_load(filename);
  
  for (i=0; i < no_seq; i++)
    {
      bioseq_random(&query, 10, "QWERTYIPASDFGHKLCVNM",
		   ULINT a_len)
  hit_list = FS_INDEX_QD_search(FS_index, &query, cutoff);
  HIT_LIST_sort_by_sequence(hit_list);
  HIT_LIST_sort_incr(hit_list);
  HIT_LIST_print(hit_list, stdout, 0); 



 return EXIT_SUCCESS;
}
