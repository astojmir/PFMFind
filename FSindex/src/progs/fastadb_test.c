/********************************************************************
 *                                                                  *
 *  fastadb_test.c - Tests the routines in fastadb library.         *
 *                                                                  *
 *              Aleksandar Stojmirovic, March 4th 2002              *
 *                                                                  *
 ********************************************************************/

#include "fastadb.h"
#include <time.h>
#include <stdlib.h>


static 
void fastadb_test(char *db_name, fastadb_arg *argt, fastadb_argv_t *argv)
{
  SEQUENCE_DB *s_db;
  BIOSEQ *seq;
#ifdef COUNTS
  ULINT length;
  ULINT no_seqs;
#endif
  yn_bool_t defline;
#ifdef RANDOM
  ULINT i;
#endif
  fastadb_arg *argt0 = argt;
  fastadb_argv_t *argv0 = argv;

  while (*argt0 != NONE)
    {
      switch (*argt0)
	{
	case (RETREIVE_DEFLINES):
	  defline = argv0->retrieve_deflines;
	  break;
	default:
	  break;
	}
      argt0++;
      argv0++;
    }
#ifdef COUNTS
  printf("******* TESTING %s *******\n", db_name); 
#endif
  s_db = fastadb_open(db_name, argt, argv);
#ifdef COUNTS
  fastadb_get_length(s_db, &length);
  fastadb_get_noseqs(s_db, &no_seqs);
  printf("Database length: %ld\n", length);
  printf("Database number of sequences: %ld\n", no_seqs);
  printf("******* TESTING SEQUENTIAL ACCESS *******\n"); 
#endif
  fastadb_get_seq(s_db, 0, &seq);
  bioseq_parse(seq, NULL, defline);
  while (fastadb_get_next_seq(s_db, &seq))
    bioseq_parse(seq, NULL, defline);
#ifdef RANDOM
  printf("******* TESTING RANDOM ACCESS *******\n"); 
  srand48(time(NULL)); 
  for (i=0; i<100; i++)
    {
      fastadb_get_seq(s_db, lrand48()%no_seqs, &seq);    
      bioseq_parse(seq, NULL, defline);
    }
#endif
  fastadb_close(s_db);
}

static 
void fastadb_frag_test(char *db_name, fastadb_arg *argt, 
		       fastadb_argv_t *argv)
{
  /* Idea is to get a few fragments to see if this is going to
     work. Not really a detailed test. */
  ULINT i = 0;
  SEQUENCE_DB *s_db;
  BIOSEQ frag;

  s_db = fastadb_open(db_name, argt, argv);
  fastadb_init_frags(s_db, 15, 30);
  for (i=790; i < 1000; i++)
    {
      printf("i= %ld\n",i);
      memset(&frag, 0, sizeof(BIOSEQ));
      fastadb_get_frag(s_db, &frag, i, 15, 30);
      bioseq_parse(&frag, NULL, NO);
    }      
  fastadb_clear_frags(s_db);
  fastadb_close(s_db);
}

int main(int argc, char **argv)
{
  char *db_name;
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  int type = 255;

  if (argc < 3)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s database type\n", argv[0]); 
      exit(EXIT_FAILURE);
    }
  db_name = argv[1];
  type = atoi(argv[2]);
  fastadb_argt[2] = NONE;
  
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  if (type == 0)
    {
      fastadb_argv[0].access_type = MEMORY;
      fastadb_argv[1].retrieve_deflines = NO;
      fastadb_test(db_name, fastadb_argt, fastadb_argv);
    }
  if (type == 1)
    {
      fastadb_argv[0].access_type = RANDOM;
      fastadb_argv[1].retrieve_deflines = NO;
      fastadb_test(db_name, fastadb_argt, fastadb_argv);
    }
 
  /* Now switch order to see if union works as it should */
  fastadb_argt[0] = RETREIVE_DEFLINES;
  fastadb_argt[1] = ACCESS_TYPE;
  if (type == 2)
    {
      fastadb_argv[0].retrieve_deflines = YES;
      fastadb_argv[1].access_type = MEMORY;
      fastadb_test(db_name, fastadb_argt, fastadb_argv);
    }
  if (type == 3)
    {
      fastadb_argv[0].retrieve_deflines = YES;
      fastadb_argv[1].access_type = RANDOM;
      fastadb_test(db_name, fastadb_argt, fastadb_argv);
    }
  if (type == 4)
    {
      fastadb_argv[0].retrieve_deflines = YES;
      fastadb_argv[1].access_type = MEMORY;
      fastadb_frag_test(db_name, fastadb_argt, fastadb_argv);    
    }

  return 0;
}
