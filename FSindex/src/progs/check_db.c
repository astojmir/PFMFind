/********************************************************************

  check_db.c

  counts letters in molecular sequence database

 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "misclib.h"
#include "fastadb.h"
#include "partition.h"

int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;


int main(int argc, char **argv)
{
  ULINT i;
  ULINT j;
  ULINT c4 = 0;
  ULINT c5 = 0;

  BIOSEQ *frag = mallocec(sizeof(BIOSEQ));
  FS_SEQ_t FS_seq;
  SEQUENCE_DB *s_db;
  FS_PARTITION_t *ptable4;
  FS_PARTITION_t *ptable5;
  int flag4;
  int flag5;
  ULINT frag_len = 10;
  ULINT seq_id;
  ULINT rel_offset; 
  const char *database;
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];

  if (argc > 1)
    {
      database = argv[1];
    }
  else
    {
      fprintf(stderr, "Use %s fasta_db\n", argv[0]);
      exit(EXIT_FAILURE);
    }

  /* Load database */
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = NONE;
  fastadb_argv[0].access_type = MEMORY;
  fastadb_argv[1].retrieve_deflines = NO;

  s_db = fastadb_open(database, fastadb_argt, fastadb_argv); 
  fastadb_init_Ffrags(s_db, frag_len);

  /* Create partition table */
  ptable4 = FS_PARTITION_create("STAN#ILVM#KRDEQ#WFYHGPC", '#'); 
  ptable5 = FS_PARTITION_create("STAN#ILVM#KRDEQ#WFYH#GPC", '#'); 

  j = 0;
  while (fastadb_get_next_Ffrag(s_db, frag_len, frag, &i, 1))  
    {
      j++;
      flag4=0;
      flag5=0;
      /* Calculate its FS_seq */
      if (BIOSEQ_2_FS_SEQ(frag, ptable4, &FS_seq))
	{
	  flag4=1;
	  c4++;
	}
      if (BIOSEQ_2_FS_SEQ(frag, ptable5, &FS_seq))
	{
	  flag5=1;
	  c5++;
	}
      else
	{
	  fprintf(stderr, "%ld %.*s\n", j, (int) frag->len, frag->start);
	}

      if (flag4 != flag5)
	{
	  fastadb_find_Ffrag_seq(s_db, frag, &seq_id, &rel_offset); 
	  printf("%d%d SEQ: %ld, OFFSET: %ld, FRAGMENT: %ld %.*s\n",
		 flag4, flag5,
		 seq_id, rel_offset, j, (int) frag_len, frag->start);
	}
    }
  printf("j = %ld; c4 = %ld; c5 = %ld\n",j,c4,c5);   
  exit(EXIT_SUCCESS);
}
