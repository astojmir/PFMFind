#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "hit_list.h"
#include "partition.h"
#include "smatrix.h"
#include "FSindex.h"

/********************************************************************/
/*                                                                  */
/*                     SEQUENTIAL SCAN module                       */
/*                                                                  */
/********************************************************************/

/* For now, all prototypes in FSindex.h  */


HIT_LIST_t *SSCAN_QD_search(SEQUENCE_DB *s_db, const char *matrix, 
			    BIOSEQ *query, ULINT D_cutoff)
{
  ULINT one_percent_fragments;

  ULINT i, j;
  BIOSEQ *frag = mallocec(sizeof(BIOSEQ));


  ULINT frag_len = query->len;
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  int D_value;
  HIT_LIST_t *hit_list;
  ULINT no_frags;
  FS_SEQ_t FS_seq;
  ULINT len;


  no_frags = fastadb_count_Ffrags(s_db, frag_len);
  one_percent_fragments = no_frags / 100;
  fastadb_init_Ffrags(s_db, frag_len);

  /* Create partition table */
  ptable = FS_PARTITION_create("STAN#ILVM#KR#EDQ#WFYH#GPC", '#'); 

  /* Load matrices, create quasi-metric */
  S = SCORE_MATRIX_create(matrix, ptable); 
  D = SCORE_MATRIX_S_2_Dquasi(S); 

  hit_list = HIT_LIST_create(query, s_db, strdup(matrix), D,
			     smatrix_eval, D_cutoff);
  frag_len = query->len;



  /* For each fragment */ 
  j = 0;
  i = 0;
  while (i < s_db->seq_data_len)  
    {
      len = strlen(s_db->seq_data + i);
      if (len < frag_len)
	{
	  i += len + 1;
	  continue;
	}
      
      frag->len = frag_len;
      frag->start = s_db->seq_data + i;
      i++;
      if (!BIOSEQ_2_FS_SEQ(frag, ptable, &FS_seq))
	continue;
      if (SCORE_MATRIX_evaluate_max(D, query, frag, 
				    D_cutoff, &D_value))
	{
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  HIT_LIST_insert_seq_hit(hit_list, frag, 
				  (float) D_value); 
	}
      HIT_LIST_count_seq_visited(hit_list, 1);
      j++;
   }

  /* Take time */
  HIT_LIST_stop_timer(hit_list);

  /* Clean up */
  free(frag);

  /* Free matrices */
  SCORE_MATRIX_destroy(D);
  SCORE_MATRIX_destroy(S);

  /* Free partition table */
  FS_PARTITION_destroy(ptable);  

  return hit_list; 
}

int SSCAN_has_neighbour(SEQUENCE_DB *s_db,  SCORE_MATRIX_t *D, 
			FS_PARTITION_t *ptable,
			BIOSEQ *query, ULINT D_cutoff)
{
  ULINT i;
  BIOSEQ frag;


  ULINT frag_len = query->len;
  int D_value;
  ULINT no_frags;
  FS_SEQ_t FS_seq;
  ULINT len;


  no_frags = fastadb_count_Ffrags(s_db, frag_len);
  fastadb_init_Ffrags(s_db, frag_len);
  frag_len = query->len;



  /* For each fragment */ 
  i = 0;
  while (i < s_db->seq_data_len)  
    {
      len = strlen(s_db->seq_data + i);
      if (len < frag_len)
	{
	  i += len + 1;
	  continue;
	}
      
      frag.len = frag_len;
      frag.start = s_db->seq_data + i;
      i++;
      if (!BIOSEQ_2_FS_SEQ(&frag, ptable, &FS_seq))
	continue;
      if (SCORE_MATRIX_evaluate_max(D, query, &frag, 
				    D_cutoff, &D_value))
	{
#if 0
	  fprintf(stderr, "%d\n%.*s\n%.*s\n\n", D_value,
		  (int) query->len, query->start, 
		  (int) frag_len, frag.start);
#endif
	  return 1;
	}
   }
  return 0; 
}
