#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "fastadb.h"
#include "FSindex.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;

static 
void print_help(const char *progname)
{
  fprintf(stderr, "Usage: %s options\n"
	  "\n"
	  "Mandatory:\n"
	  "-F  f   Use FSindex index file f\n"
	  "-b  s   Bin sequence \n"
	  "Optional:\n"
	  "-s      Show sequences without descriptions\n"
	  "-c      Show sequences with descriptions\n"
	  "-p  f   Show expectation based on product measure.\n"
	  "        f is a file with amino acid frequencies\n"
	  "\n"
	  ,progname);
}





int main(int argc, char **argv)
{
  /* Index data */
  const char *filename = NULL;
  FSINDX *FSI;
  int frag_len;
  SEQUENCE_DB *s_db;
  
  /* Matrices */
  FS_PARTITION_t *ptable = NULL;
#if 0
  SCORE_MATRIX_t *S = NULL;
  SCORE_MATRIX_t *D = NULL;
  char *matrix_full = NULL;
  char *matrix_base = NULL;
  char *matrix_dir = NULL;
#endif

  /* Getopt variables */
  int c;                              /* Option character          */ 
  extern char *optarg;                /* Option argument           */
  extern int optind;
  extern int optopt;
  extern int opterr;  
  int errflg = 0;                     /* Option error flag         */
  
  FS_SEQ_t FS_query;
  int j;
  ULINT frag_offset = 0;
  BIOSEQ subject;
  int m=0;
  char sc[2];
  int C;
  int K;
  int KK;
  char *seq = NULL;
  ULINT id;
  ULINT from;
  char *defline;

  while ((c = getopt(argc, argv, "F:b:")) != EOF)
    switch (c) 
      {
      case 'F':
	filename = optarg;
	break;
      case 'b':
	m = strlen(optarg);
	seq = optarg;
	break;
      case '?':
      default:
	errflg++;
	break;
      }

  if (filename == NULL || seq == NULL)
    errflg++;

  if (errflg)
    {
      fprintf(stderr,"Error: Insufficient or invalid arguments \n");
      print_help(argv[0]);
      exit(EXIT_FAILURE);
    }


  fprintf(stdout, "Loading index... \n");
  FSI = FS_INDEX_load(filename);
  ptable = FS_INDEX_get_ptable(FSI);
  s_db = FS_INDEX_get_database(FSI);
  frag_len = FS_INDEX_get_frag_len(FSI);
  if (m != frag_len)
    {
      fprintf(stderr, "Wrong sequence length\n");
      exit(EXIT_FAILURE);
    }

  FS_query = 0;
  K = FS_PARTITION_get_no_partitions(ptable);
  KK = 1;
  for (KK=1, j=0; j < m; j++, KK*=K)
    {
      sc[0] = seq[j];
      C = atoi(sc);
      if (C > K)
	{
	  fprintf(stderr, "Invalid sequence\n");
	  exit(EXIT_FAILURE);
	}
      FS_query += C*KK;
    }
  fprintf(stdout, "Bin: %s\n", seq);
  fprintf(stdout, "Size: %ld\n", 
	  FS_INDEX_get_no_seqs(FSI, FS_query)); 
  fprintf(stdout, "'Unique' size: %ld\n\n",
	  FS_INDEX_get_no_useqs(FSI, FS_query));
  for (j=0; j < FS_INDEX_get_no_seqs(FSI, FS_query); j++)
    {
      frag_offset = FS_INDEX_retrieve_seq(FSI, FS_query, j);
      fastadb_get_Ffrag(s_db, frag_len, &subject, frag_offset);
      fastadb_find_Ffrag_seq(s_db, &subject, &id, &from);
      defline = s_db->seq[id].id.defline;
      fprintf(stdout, "%*.*s  %s\n", m, m, subject.start, defline);
    }
 

 return EXIT_SUCCESS;
}
