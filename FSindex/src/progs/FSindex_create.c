#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "FSindex.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 0;
int POS_MATRIX_VERBOSE;
FILE *POS_MATRIX_STREAM;

static 
void print_help(const char *progname)
{
  fprintf(stderr, "Usage: %s options\n"
	  "\n"
	  "Mandatory:\n"
	  "-d  f   Fasta sequence database\n"
	  "-p  s   Alphabet partitions (separated by '#')\n"
	  "-m  n   Fragment length\n"
	  "\n"
	  "Optional:\n"
	  "-F  f   Write index to file f\n"
	  "-s  n   Skip n fragments\n"
	  "-v      Verbose: print indexing messages\n"
	  "-V      Verbose: print partitioning messages as well\n"
	  "-b      Print progress bar\n"
	  "\n"
	  ,progname);
}


int main(int argc, char **argv)
{
  int c;                              /* Option character          */ 
  extern char *optarg;                /* Option argument           */
  int errflg = 0;                     /* Option error flag         */

  const char separator = '#';
  ULINT frag_len = 0;
  const char *db_name = NULL;
  const char *alphabet = NULL;
  char *filename = NULL;
  char *index_dir;
  char *db_dir;
  char *index_base;
  char *db_base;
  int skip = 1;

  int filename_flag = 0;
  FSINDX *FSI;

  while ((c = getopt(argc, argv, "d:p:m:F:s:Vvb")) != EOF)
    switch (c) 
      {
      case 'd':
	db_name = optarg;
	break;
      case 'p':
	alphabet = optarg;
	break;
      case 'm':
	frag_len = atoi(optarg);
	break;
      case 'F':
	filename = optarg;
	filename_flag = 1;
	break;
      case 's':
	skip = atoi(optarg);
	if (skip < 1)
	  errflg++;
	break;
      case 'V':
	FS_PARTITION_VERBOSE = 1;
      case 'v':
	FS_INDEX_VERBOSE = 1;
	break;
      case 'b':
	FS_INDEX_PRINT_BAR = 1;
	break;
      case '?':
      default:
	errflg++;
	break;
      }

  if ((db_name == NULL) || (alphabet == NULL) || (frag_len == 0))
    errflg++;

  if (errflg)
    {
      fprintf(stderr,"Error: Insufficient or invalid arguments \n");
      print_help(argv[0]);
      exit(EXIT_FAILURE);
    }

  if (filename_flag)
    {
      split_base_dir(filename, &index_base, &index_dir);
      split_base_dir(db_name, &db_base, &db_dir);
      if (strcmp(db_dir, index_dir) != 0)
	{
	  fprintf(stderr, "Warning! The directories of database file "
		  "and index file must match.\nSetting the index directory"
		  " to the database directory.\n");
	  cat_base_dir(&filename, index_base, db_dir);
	  fprintf(stderr, "The new path is %s\n", filename);      
	}
    }

  printf("Creating index...\n");
  FSI = FS_INDEX_create(db_name, frag_len, alphabet, separator, skip); 

  if (filename_flag)
    {
      printf("Saving index...\n");
      FS_INDEX_save(FSI, filename);
    }

  return EXIT_SUCCESS;
}
