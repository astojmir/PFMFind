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

static 
void print_help(const char *progname)
{
  fprintf(stderr, "%s Arguments:\n\n"
	  "-F  f   Use FSindex index file f (Mandatory)\n"
	  "-M  f   Use score matrix file f (Mandatory)\n"
	  "[-s | -q | -d | -p] n (Mandatory)\n"
	  "        Use cutoff value n with\n"
	  "           -s similarity scores\n"
	  "           -q non-symmetric distance scores\n"
	  "           -d symmetric distance scores\n"
	  "           -p PSM scores\n"
	  "-i  f   Read input from file f (Optional)\n"
	  "          If ommitted, read from stdin\n"
	  "-o  f   Write output to file f (Optional)\n"
	  "          If ommitted, write to stdout\n"
	  "-I  n   Run n profile iterations (Optional)\n"
	  "          If n = 0 run until convergence\n"
	  "-L  n   Scaling factor to use when producing PSMs\n"
	  "-d  f   Read background distribution of amino acids from"
	  " file f\n", progname);
}





int main(int argc, char **argv)
{
  const char *query_seq;
  int cutoff;
  HIT_LIST_t *hit_list;
  BIOSEQ query;
  const char *filename;
  char *matrix_full;
  char *matrix_base;
  char *matrix_dir;
  int no_args = 4;
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  POS_MATRIX *PS;
  POS_MATRIX *PD;


  enum {T_SIMILARITY, T_QUASI_METRIC, T_METRIC, T_PROFILE} search_type 
    = T_SIMILARITY; 

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix query_seq"
	      " cutoff [-s | -m | -q];\n", argv[0]);
      print_help(argv[0]);
      exit(EXIT_FAILURE);

    }
  filename = argv[1];
  matrix_full = argv[2]; 
  query_seq = argv[3];
  cutoff = atoi(argv[4]);

  if (argc >= no_args+2)
    {
      if (argv[5][0] == '-')
	switch (argv[5][1])
	  {
	  case 'p':
	    search_type = T_PROFILE;
	    break;
	  case 's':
	    search_type = T_SIMILARITY;
	    break;
	  case 'm':
	    search_type = T_METRIC;
	    break;
	  case 'q':
	    search_type = T_QUASI_METRIC;
	    break;
	  default:
	    search_type = T_SIMILARITY;
	    break;
	  }
    }
  fprintf(stdout,"Loading database and index ... \n");
  query.len = strlen(query_seq);
  query.start = strdup(query_seq);

  FS_INDEX_load(filename);
  ptable = FS_INDEX_get_ptable();
  split_base_dir(matrix_full, &matrix_base, &matrix_dir);  
  S = SCORE_MATRIX_create(matrix_full, ptable); 
  hit_list = HIT_LIST_create(&query, FS_INDEX_get_database(), 
			     matrix_base, cutoff);

  fprintf(stdout,"Searching ... \n");
  switch (search_type)
    {
    case T_PROFILE:
      PS = SCORE_2_POS_MATRIX(S, &query);
      PD = POS_MATRIX_S_2_D(PS, &query);      
      FS_INDEX_profile_search(hit_list, &query, PS, PD, 0,
			      query.len-1, cutoff,
			      FS_INDEX_profile_S_process_bin,
			      FS_INDEX_S2QD_convert);
      HIT_LIST_sort_by_sequence(hit_list);
      HIT_LIST_sort_decr(hit_list);
      break;
    case T_SIMILARITY:
      D = SCORE_MATRIX_S_2_Dquasi(S);
      FS_INDEX_search(hit_list, &query, S, D, cutoff,
		      FS_INDEX_S_process_bin,
		      FS_INDEX_S2QD_convert);
      HIT_LIST_sort_by_sequence(hit_list);
      HIT_LIST_sort_decr(hit_list);
      break;
    case T_QUASI_METRIC:
      D = SCORE_MATRIX_S_2_Dquasi(S);
      FS_INDEX_search(hit_list, &query, S, D, cutoff,
		      FS_INDEX_QD_process_bin, 
		      FS_INDEX_identity_convert);
      HIT_LIST_sort_by_sequence(hit_list);
      HIT_LIST_sort_incr(hit_list);
     break;
    case T_METRIC:
      D = SCORE_MATRIX_S_2_Dmax(S);
      FS_INDEX_search(hit_list, &query, S, D, cutoff,
		      FS_INDEX_QD_process_bin, 
		      FS_INDEX_identity_convert);
      HIT_LIST_sort_by_sequence(hit_list);
      HIT_LIST_sort_incr(hit_list);
      break;
    default:
      break;
    }

  HIT_LIST_print(hit_list, stdout, 0); 

  return EXIT_SUCCESS;
}
