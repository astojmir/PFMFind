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
int POS_MATRIX_VERBOSE = 0;
FILE *POS_MATRIX_STREAM;

/* TO DO:

- main conversions
- profiles transformed into distances
- back transformation of hit lists prior to output

*/



static 
void print_help(const char *progname)
{
  fprintf(stderr, "Usage: %s [options] index matrix cutoff\n"
	  "\n"
	  "Arguments:\n"
	  "index    - FSindex file name;\n"
	  "matrix   - similarity score matrix file name;\n"
	  "cutoff   - cutoff value.\n\n"
	  "General Options:\n" 
	  "-i  f   Read input from file f (Default = stdin);\n"
	  "-o  f   Write output to file f (Default = stdout);\n"
	  "-t  [s|q|d|p] Search type:\n"
	  "     s similarity;\n"
	  "     q non-symmetric distance (default);\n"
	  "     d symmetric distance;\n"
	  "     p PSM similarity;\n"
	  "-c  [s|k] Cutoff type:\n"
	  "     s score (default);\n"
          "     k number of nearest neighbours.\n"
	  "-v      Verbose output\n"
	  "\n"
	  "Profile General Options (only if PSM search type specified):\n"
	  "-I  n   Run n profile iterations (Default = 5)\n" 
	  "          If n = 0 run until convergence\n"
	  "-L  x   Scaling factor to use when producing PSMs"
	  " (Default = 1.0)\n"
	  "-S  x   Cutoff score to use on second and subsequent"
	  " iterations\n"
	  "         (Default = the same score as specified in -p"
	  " option)\n" 
	  "\n"
	  "Profile Pseudo-count Options:\n"
	  "-A  x   Add x pseudo counted sequences (Default = 20.0)\n"
	  "-f  f   Read background amino acid distribution from" 
	  "file f\n"
	  "        (must be specified for PSM similarity searches)\n"
	  "\n"
	  ,progname);
}





int main(int argc, char **argv)
{
  int c;                              /* Option character          */ 
  extern char *optarg;                /* Option argument           */
  extern int optind;
  extern int optopt;
  extern int opterr;

  FSINDX *FSI = NULL;
  HIT_LIST_t *HL = NULL;
  BIOSEQ *query = NULL;

  const char *filename = NULL;
  const char *in_file = NULL;
  const char *out_file = NULL;
  FILE *out_stream = stdout;
  SEQUENCE_DB *input_db = NULL;
  fastadb_arg fastadb_argt[5];
  fastadb_argv_t fastadb_argv[5];

  char *matrix_full = NULL;

  FS_PARTITION_t *ptable = NULL;
  SCORE_MATRIX_t *S = NULL;
  SCORE_MATRIX_t *D = NULL;

  double lambda = 1.0;
  double A = 20.0;
  const char *freq_filename = NULL;
  int iters = 5;

  int cutoff  = 0;
  int cutoff2 = 0;

  int d0;
  int s0;

  int cutoff2_flag = 0;


  enum {T_SIMILARITY, T_QUASI_METRIC, T_METRIC, T_PROFILE} search_type 
    = T_QUASI_METRIC; 
  enum {SCORE, kNN} cutoff_type = SCORE; 


  POS_MATRIX_STREAM = stdout;

  /* Input database arguments */
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = MAX_ROW_LENGTH;
  fastadb_argt[3] = KEEP_OLDSEQS;
  fastadb_argt[4] = NONE;
  fastadb_argv[0].access_type = SEQUENTIAL;
  fastadb_argv[1].retrieve_deflines = YES;
  fastadb_argv[2].max_row_length = 1 << 16; /* 2^16 - 65 kB */ 
  fastadb_argv[3].keep_oldseqs = NO; 


  /* Process options */
  while ((c = getopt(argc, argv, "i:o:t:c:f:I:L:S:A:v")) != EOF)
    switch (c) 
      {
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
	POS_MATRIX_STREAM = out_stream;
	break;
      case 't':
	switch (*optarg)
	  {
	  case 'p':
	    search_type = T_PROFILE;
	    break;
	  case 's':
	    search_type = T_SIMILARITY;
	    break;
	  case 'q':
	    search_type = T_QUASI_METRIC;
	    break;
	  case 'd':
	    search_type = T_METRIC;
	    break;
	  default:
	    fprintf(stderr, "Unreckognised -t argument.\n");
	    exit(EXIT_FAILURE);
	  }
	break;
      case 'c':
	switch (*optarg)
	  {
	  case 's':
	    cutoff_type = SCORE;
	    break;
	  case 'k':
	    cutoff_type = kNN;
	    break;
	  default:
	    fprintf(stderr, "Unreckognised -c argument.\n");
	    exit(EXIT_FAILURE);
	  }
	break;
      case 'f':
	freq_filename = optarg;
	break;
      case 'I':
	iters = atoi(optarg);
	break;
      case 'L':
	lambda = atof(optarg);
	break;
      case 'S':
	cutoff2 = atoi(optarg);
	cutoff2_flag++;
	break;
      case 'A':
	A = atof(optarg);
	break;
      case 'v':
	POS_MATRIX_VERBOSE = 1;
	break;
      case '?':
      default:
	fprintf(stderr,"Error: Invalid argument. \n");
	exit(EXIT_FAILURE);
	break;
      }

  /* Process mandatory arguments */ 

  if (optind + 3 > argc)
    {
      fprintf(stderr, "Missing arguments.\n");
      print_help(argv[0]);
      exit(EXIT_FAILURE);
    }

  filename = argv[optind];
  matrix_full = argv[optind+1];
  cutoff = atoi(argv[optind+2]);

  if (freq_filename == NULL && search_type == T_PROFILE)
    {
      fprintf(stderr,"Missing amino acid frequency filename.\n");
      exit(EXIT_FAILURE);
    }


  fprintf(stdout,"Loading database and index ... \n");

  if (!cutoff2_flag) cutoff2 = cutoff;

  input_db = fastadb_open(in_file, fastadb_argt, fastadb_argv); 
  FSI = FS_INDEX_load(filename);
  ptable = FS_INDEX_get_ptable(FSI);
  S = SCORE_MATRIX_create(matrix_full, ptable); 

  fprintf(stdout,"Searching ... \n");


  /* Convert matrices for range searches*/
  switch (search_type)
    {
    case T_PROFILE:
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


  /* Iterate over all query sequences */
  while (fastadb_get_next_seq(input_db, &query))
    {
      switch (search_type)
	{
	case T_PROFILE:
	  if (cutoff_type == SCORE)
	    {
	      HL = FSINDX_prof_rng_srch(FSI, query, D, cutoff, lambda,
					A, freq_filename, iters,
					 cutoff2, HL);
	      HIT_LIST_sort_decr(HL);
	    }
	  else /* cutoff_type == kNN */
	    {
	      HL = FSINDX_prof_kNN_srch(FSI, query, D, cutoff, 
					A, freq_filename, iters, HL); 
	    }

	  HIT_LIST_print(HL, out_stream, 0); 	
	  break;
	case T_SIMILARITY:
	  s0 = SCORE_MATRIX_evaluate(S, query, query);
	  if (cutoff_type == SCORE)
	    {
	      d0 = s0 - cutoff;
	      HL = FSINDX_rng_srch(FSI, query, D, d0, HL);
	    }
	  else /* cutoff_type == kNN */
	    {
	      HL = FSINDX_kNN_srch(FSI, query, D, cutoff, HL);
	    }
	  HIT_LIST_sort_incr(HL);
	  SCORE_MATRIX_convert(s0, HL);
	  HIT_LIST_print(HL, out_stream, 0); 
	  break;
	case T_QUASI_METRIC:
	case T_METRIC:
	  if (cutoff_type == SCORE)
	    {
	      d0 = cutoff;
	      HL = FSINDX_rng_srch(FSI, query, D, d0, HL);
	    }
	  else /* cutoff_type == kNN */
	    {
	      HL = FSINDX_kNN_srch(FSI, query, D, cutoff, HL);
	    }
	  HIT_LIST_sort_incr(HL);
	  HIT_LIST_print(HL, out_stream, 0); 
	  break;
	default:
	  break;
	}
    }
  return EXIT_SUCCESS;
}
