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
  const char *db_name;
  const char *matrix;
  const char *query_seq;
  int cutoff;
  HIT_LIST_t *hit_list;
  BIOSEQ query;
  SEQUENCE_DB *s_db;
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];

  if (argc < 5)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s database matrix query_seq"
	      " cutoff\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  db_name = argv[1];
  matrix = argv[2];
  query_seq = argv[3];
  cutoff = atoi(argv[4]);

  query.len = strlen(query_seq);
  query.start = strdup(query_seq);

  /* Load database */
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = NONE;
  fastadb_argv[0].access_type = MEMORY;
  fastadb_argv[1].retrieve_deflines = NO;

  s_db = fastadb_open(db_name, fastadb_argt, fastadb_argv); 

  hit_list = SSCAN_QD_search(s_db, matrix, &query, cutoff);

  HIT_LIST_sort_by_sequence(hit_list);
  HIT_LIST_sort_incr(hit_list);
  HIT_LIST_print(hit_list, stdout, 0); 

  /* Free sequence database */
  fastadb_close(s_db);


 return EXIT_SUCCESS;
}
