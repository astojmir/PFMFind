#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "FSindex.h"
#include "hit_list.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;

static 
void print_help(const char *progname)
{
  fprintf(stderr, "Usage: %s options\n"
	  "Options:\n"
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
	  "-S  n   Cutoff score to use on second and subsequent"
	  " iterations\n"
	  "-f  f   Read background distribution of amino acids from"
	  " file f\n" 
	  "-A  x   Add x pseudo counted sequences\n", progname);
}





int main(int argc, char **argv)
{
  int c;                              /* Option character          */ 
  extern char *optarg;                /* Option argument           */
  extern int optind;
  extern int optopt;
  extern int opterr;
  int errflg = 0;                     /* Option error flag         */


  int cutoff;
  HIT_LIST_t *hit_list;
  BIOSEQ *query = NULL;

  const char *filename;
  int index_flag = 0;

  const char *in_file = NULL;
  const char *out_file = NULL;
  FILE *out_stream = stdout;
  SEQUENCE_DB *input_db;
  fastadb_arg fastadb_argt[4];
  fastadb_argv_t fastadb_argv[4];

  char *matrix_full;
  char *matrix_base;
  char *matrix_dir;
  int matrix_flag = 0;

  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  POS_MATRIX *PS;
  double lambda = 1.0;
  double A = 20.0;
  const char *freq_filename = NULL;
  int cutoff2;

  enum {T_SIMILARITY, T_QUASI_METRIC, T_METRIC, T_PROFILE} search_type 
    = T_SIMILARITY; 
  int cutoff_flag = 0;

  while ((c = getopt(argc, argv, "F:M:s:d:q:p:i:o:I:L:S:f:A:")) != EOF)
    switch (c) 
      {
      case 'F':
	filename = optarg;
	index_flag++;
	break;
      case 'M':
	matrix_full = optarg;
	matrix_flag++;
	break;
      case 'p':
	search_type = T_PROFILE;
	cutoff = atoi(optarg);
	cutoff_flag++;
	break;
      case 's':
	search_type = T_SIMILARITY;
	cutoff = atoi(optarg);
	cutoff_flag++;
	break;
      case 'd':
	search_type = T_METRIC;
	cutoff = atoi(optarg);
	cutoff_flag++;
	break;
      case 'q':
	search_type = T_QUASI_METRIC;
	cutoff = atoi(optarg);
	cutoff_flag++;
	break;
      case 'i':
	in_file = optarg;
	break;
      case 'o':
	out_file = optarg;
	if((out_stream = fopen(out_file, "w")) == NULL)
	  {
	    fprintf(stderr, 
		    "Could not open output file %s!\n", out_file);
	    exit(EXIT_FAILURE);
	  }
	break;
      case 'L':
	lambda = atof(optarg);
	break;
      case 'f':
	freq_filename = optarg;
	break;
      case 'A':
	A = atof(optarg);
	break;
      case 'S':
	cutoff2 = atoi(optarg);
	break;

      case '?':
      default:
	errflg++;
	break;
      }

  if (!index_flag)
    {
      errflg++;
      fprintf(stderr,"Error: Missing index filename.\n");
    }
  if (!matrix_flag)
    {
      errflg++;
      fprintf(stderr,"Error: Missing score matrix filename.\n");
    }
  if (!cutoff_flag)
    {
      errflg++;
      fprintf(stderr,"Error: Missing search type and cutoff"
	      " value.\n");
    }
  else if (cutoff_flag > 1)
    {
      fprintf(stderr,"Warning: Too many search types and cutoff"
	      " values specified.\n");
    }
  if (errflg)
    {
      fprintf(stderr,"Error: Insufficient or invalid arguments \n");
      print_help(argv[0]);
      exit(EXIT_FAILURE);
    }

  fprintf(stdout,"Loading database and index ... \n");
  FS_INDEX_load(filename);
  ptable = FS_INDEX_get_ptable();
  split_base_dir(matrix_full, &matrix_base, &matrix_dir);  
  S = SCORE_MATRIX_create(matrix_full, ptable); 
  hit_list = HIT_LIST_create(query, FS_INDEX_get_database(), 
			     matrix_base, cutoff);

  fprintf(stdout,"Searching ... \n");
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = MAX_ROW_LENGTH;
  fastadb_argt[3] = KEEP_OLDSEQS;
  fastadb_argv[0].access_type = SEQUENTIAL;
  fastadb_argv[1].retrieve_deflines = YES;
  fastadb_argv[2].max_row_length = 1 << 16; /* 2^16 - 65 kB */ 
  fastadb_argv[3].keep_oldseqs = NO; 

  input_db = fastadb_open(in_file, fastadb_argt, fastadb_argv); 

  switch (search_type)
    {
    case T_SIMILARITY:
    case T_QUASI_METRIC:
      D = SCORE_MATRIX_S_2_Dquasi(S);
     break;
    case T_METRIC:
      D = SCORE_MATRIX_S_2_Dmax(S);
      break;
    default:
      break;
    }

  while (fastadb_get_next_seq(input_db, &query))
    {
      HIT_LIST_reset(hit_list, query, FS_INDEX_get_database(),
		     matrix_base, cutoff);
      
      switch (search_type)
	{
	case T_PROFILE:
	  PS = SCORE_2_POS_MATRIX(S, query);
	  POS_MATRIX_init(PS, lambda, POS_MATRIX_simple_pseudo_counts,
			  freq_filename);
	  POS_MATRIX_simple_pseudo_counts_init(PS, A);

	  
	  FS_INDEX_profile_search(hit_list, PS,0, query->len-1, cutoff, 
				  FS_INDEX_profile_S_process_bin,
				  FS_INDEX_S2QD_convert);
      
	  HIT_LIST_sort_by_sequence(hit_list);
	  HIT_LIST_sort_decr(hit_list);

	  HIT_LIST_get_hit_seqs(hit_list, &PS->seq, cutoff, 
				&PS->no_seqs, &PS->max_no_seqs);
	  POS_MATRIX_equal_weights(PS);
	  POS_MATRIX_update(PS);

	  break;
	case T_SIMILARITY:
	  FS_INDEX_search(hit_list, query, S, D, cutoff,
			  FS_INDEX_S_process_bin,
			  FS_INDEX_S2QD_convert);
	  HIT_LIST_sort_by_sequence(hit_list);
	  HIT_LIST_sort_decr(hit_list);
	  break;
	case T_QUASI_METRIC:
	case T_METRIC:
	  D = SCORE_MATRIX_S_2_Dmax(S);
	  FS_INDEX_search(hit_list, query, S, D, cutoff,
			  FS_INDEX_QD_process_bin, 
			  FS_INDEX_identity_convert);
	  HIT_LIST_sort_by_sequence(hit_list);
	  HIT_LIST_sort_incr(hit_list);
	  break;
	default:
	  break;
	}
      HIT_LIST_print(hit_list, out_stream, 0); 
    }
  return EXIT_SUCCESS;
}
