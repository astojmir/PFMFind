#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "partition.h"
#include "smatrix.h"
#include "FSindex.h"
#include "randseq.h"

int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 0;

static 
void print_help(const char *progname)
{
  fprintf(stderr, "Usage: %s options\n"
	  "\n"
	  "General Options:\n"
	  "-F  f   Use FSindex index file f (Mandatory)\n"
	  "-M  f   Use score matrix file f (Mandatory)\n"
	  "-k  n   Number of points\n"
	  "-r  n   Number of runs\n"
	  "-o  f   Write output to file f (Default = stdout)\n"
	  "-f  f   Read background amino acid distribution from" 
	  " file f"
	  "          (Mandatory)\n" 
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
  int errflg = 0;                     /* Option error flag         */


  BIOSEQ query;
  BIOSEQ subject;

  const char *filename = NULL;
  int index_flag = 0;
  FS_HASH_TABLE_t *HT;
  const char *freq_filename = NULL;
  SEQ_GENERATOR_t *seq_generator;
  char *seq_heap;

  const char *out_file = NULL;
  FILE *out_stream = stdout;

  char *matrix_full = NULL;
  char *matrix_base = NULL;
  char *matrix_dir = NULL;
  int matrix_flag = 0;

  FS_PARTITION_t *ptable = NULL;
  SCORE_MATRIX_t *S = NULL;
  SCORE_MATRIX_t *D = NULL;

  int i, j, k;
  int no_pts = 0;
  int no_runs = 0;
  int len;
  ULINT one_percent;
  int s;
  int dist;
  int poffset;
  

  while ((c = getopt(argc, argv, 
		     "F:M:k:r:o:f:")) != EOF)
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
      case 'k':
	no_pts = atoi(optarg);
	break;
      case 'r':
	no_runs = atoi(optarg);
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
      case 'f':
	freq_filename = optarg;
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
  if (freq_filename == NULL)
    {
      errflg++;
      fprintf(stderr,"Error: Missing amino acid frequency filename.\n");
    }
  if (!no_runs || !no_pts)
    {
      errflg++;
      fprintf(stderr,"Error: Missing number of samples and runs.\n");
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
  HT = FS_INDEX_get_hash_table();
  len = FS_INDEX_get_frag_len();
  seq_heap = mallocec(len+1);
  split_base_dir(matrix_full, &matrix_base, &matrix_dir);  
  S = SCORE_MATRIX_create(matrix_full, ptable); 
  D = SCORE_MATRIX_S_2_Dquasi(S);
  seq_generator = SEQ_GENERATOR_create(freq_filename, ptable);

  if (no_runs < 100)
    one_percent = 1;
  else
    one_percent = no_runs/100;


  /* calculation */

  fprintf(out_stream, "%s %d %d\n", FS_INDEX_get_db_name(), no_pts, no_runs);
  for (i=0; i < no_pts; i++)
    {
      fprintf(stdout, "pt #%d\n", i);

      SEQ_GENERATOR_rand_seq(seq_generator, &query, len, seq_heap);  
      fprintf(out_stream, "%*.*s\n", (int) query.len, (int) query.len,
	      query.start);
      for (j=0; j < no_runs; j++)
	{
	  s = HT->heap[lrand48()%FS_HASH_TABLE_get_total_seqs(HT)];
	  fastadb_get_Ffrag(FS_INDEX_get_database(), 10, &subject, s);
	  fprintf(out_stream, "%*.*s %d", (int) subject.len,
		  (int) subject.len, subject.start, 
		  SCORE_MATRIX_evaluate(D, &query, &subject));
	  for (dist=0, k=0; k < subject.len; k++)
	    {
	      poffset = FS_PARTITION_get_poffset(ptable,
		FS_PARTITION_get_pttn(ptable, subject.start[k]));	       
	      dist += D->pM[query.start[k] & A_SIZE_MASK][poffset];
	      fprintf(out_stream, " %d", dist);
	    }
	  printbar(stdout, (ULINT) j+1, one_percent, 50);  
	  fprintf(out_stream, "\n");

	}      
      fprintf(out_stream, "\n");
    }
  fprintf(stdout,"\n");
  if (out_file != NULL)
    fclose(out_stream);
  return EXIT_SUCCESS;
}
