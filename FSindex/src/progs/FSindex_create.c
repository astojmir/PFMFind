#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "FSindex.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 0;

int main(int argc, char **argv)
{
  FS_INDEX_t *FS_index;
  const char separator = '#';
  ULINT frag_len;
  const char *db_name = NULL;
  const char *alphabet = NULL;
  char *filename = NULL;
  char *index_dir;
  char *db_dir;
  char *index_base;
  char *db_base;

  int filename_flag = 0;
  int i;

  if (argc < 4)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s database clusters"
	      " fragment_length index_file [-v -m -b -V]\n",
	      argv[0]);  
      exit(EXIT_FAILURE);
    }
  db_name = argv[1];
  alphabet = argv[2];
  frag_len = atoi(argv[3]);
  i = 4;
  if (argc >= 5 && argv[4][0] != '-')
    {
      filename = argv[4];
      filename_flag = 1;
      i = 5;
    }
  while (i < argc)
    {
      if (argv[i][0] != '-')
	break;
      switch (argv[i][1])
	{
	case 'v':
	  FS_INDEX_VERBOSE = 1;
	  break;
	case 'm':
	  FS_PARTITION_VERBOSE = 1;
	  SCORE_MATRIX_VERBOSE = 1;
	  break;
	case 'b':
	  FS_INDEX_PRINT_BAR = 1;
	  break;
	case 'V':
	  FS_PARTITION_VERBOSE = 1;
	  SCORE_MATRIX_VERBOSE = 1;
	  FS_INDEX_VERBOSE = 1;
	  break;
	default:
	  break;
	}
      i++;
    }

  split_base_dir(filename, &index_base, &index_dir);
  split_base_dir(db_name, &db_base, &db_dir);
  if (strcmp(db_dir, index_dir) != 0)
    {
      fprintf(stderr, "Warning! The directories of database file "
	      "and index file must match.\n Setting the index directory"
	      "to the database directory.\n");
      cat_base_dir(&filename, index_base, db_dir);
      fprintf(stderr, "The new path is %s\n", filename);      
    }

  FS_index = FS_INDEX_create(db_name, frag_len, alphabet, 
			    separator); 

  if (filename_flag)
    {
      printf("Saving database\n");
      FS_INDEX_save(FS_index, filename);
    }

  /* This part is removed as it was used for testing only */
#if 0
  printf("Removing index\n");
  FS_INDEX_destroy(FS_index); 

  printf("Loading database\n");
  FS_index = FS_INDEX_load(filename);
  /* This takes way too long */
  printf("Removing index\n");
  FS_INDEX_destroy(FS_index); 
#endif
  return EXIT_SUCCESS;
}
