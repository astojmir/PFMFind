#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "hit_list.h"


#ifdef DEBUG
#endif

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               LIBRARY FUNCTIONS                              ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

/********************************************************************/    
/*                                                                  */
/*                    SEQ_HIT module                                */ 
/*                                                                  */
/********************************************************************/    


void SEQ_HIT_print(SEQ_HIT_t *hit, FILE *stream, 
		   SEQ_HIT_PRINT_OPT_t options)
{
  fprintf(stream, "%10ld %6ld %6ld %.*s %6.2f \n",
	  hit->sequence_id, 
	  hit->sequence_from, 
	  hit->sequence_from + hit->subject->len,
	  (int) hit->subject->len,
	  hit->subject->start,
	  hit->value);
}

static
int SEQ_HIT_cmp_decr(const void *S1, const void *S2)
{
  SEQ_HIT_t *T1 = (SEQ_HIT_t *)S1;
  SEQ_HIT_t *T2 = (SEQ_HIT_t *)S2;

  /* Want to use for sorting in decreasing order */
  return (int) ((T2)->value - (T1)->value);
}

static
int SEQ_HIT_cmp_incr(const void *S1, const void *S2)
{
  SEQ_HIT_t *T1 = (SEQ_HIT_t *)S1;
  SEQ_HIT_t *T2 = (SEQ_HIT_t *)S2;

  return (int) ((T1)->value - (T2)->value);
}

static
int SEQ_HIT_cmp_by_seq(const void *S1, const void *S2)
{
  SEQ_HIT_t *T1 = (SEQ_HIT_t *)S1;
  SEQ_HIT_t *T2 = (SEQ_HIT_t *)S2;

  if (T1->sequence_id > T2->sequence_id)
    return 1;
  else if (T1->sequence_id < T2->sequence_id)
    return -1;
  else if (T1->sequence_from > T2->sequence_from)
    return 1;
  else if (T1->sequence_from < T2->sequence_from)
    return -1;
  else
    return 0;
}


/********************************************************************/    
/*                                                                  */
/*                    HIT_LIST module                               */ 
/*                                                                  */
/********************************************************************/    
#define MAX_HITS 100

/* Main constructor */
HIT_LIST_t *HIT_LIST_create(BIOSEQ *query, SEQUENCE_DB *s_db, 
			    const char *matrix, int range)
{
  HIT_LIST_t *HIT_list = callocec(1, sizeof(HIT_LIST_t));
  HIT_list->query = query;
  HIT_list->s_db = s_db;
  HIT_list->matrix = matrix;
  HIT_list->range = range;
  HIT_list->start_time = time(NULL);
  HIT_list->max_hits = MAX_HITS;
  HIT_list->hits = mallocec(MAX_HITS * sizeof(SEQ_HIT_t));

  return HIT_list;
}  

/* Reset list */
void HIT_LIST_reset(HIT_LIST_t *HIT_list, BIOSEQ *query, 
		    SEQUENCE_DB *s_db, const char *matrix, int range)
{
  free(HIT_list->hits);
  memset(HIT_list, 0, sizeof(HIT_LIST_t));
  HIT_list->query = query;
  HIT_list->s_db = s_db;
  HIT_list->matrix = matrix;
  HIT_list->range = range;
  HIT_list->start_time = time(NULL);
  HIT_list->max_hits = MAX_HITS;
  HIT_list->hits = mallocec(MAX_HITS * sizeof(SEQ_HIT_t));
}

/* Destructor */

/* We destroy all hits stored as well */
void HIT_LIST_destroy(HIT_LIST_t *HIT_list)
{
  free(HIT_list->hits);
  free(HIT_list);
}

/* Insertion */
void HIT_LIST_insert_seq_hit(HIT_LIST_t *HIT_list, BIOSEQ *subject, 
			     float value) 
{
  BIOSEQ *new_subject = mallocec(sizeof(BIOSEQ));
  memcpy(new_subject, subject, sizeof(BIOSEQ));

  if (HIT_list->actual_seqs_hits >= HIT_list->max_hits)
    {
      HIT_list->max_hits += MAX_HITS;
      HIT_list->hits = reallocec(HIT_list->hits, HIT_list->max_hits *
				 sizeof(SEQ_HIT_t));
    }
  HIT_list->hits[HIT_list->actual_seqs_hits].subject = new_subject;
  HIT_list->hits[HIT_list->actual_seqs_hits].value = value;

  /* TO DO: get sequence_id, sequence_from from database */
  fastadb_find_Ffrag_seq(HIT_list->s_db, new_subject, 
     &(HIT_list->hits[HIT_list->actual_seqs_hits].sequence_id),  
     &(HIT_list->hits[HIT_list->actual_seqs_hits].sequence_from)); 

  HIT_list->actual_seqs_hits++;
}

void HIT_LIST_count_FS_seq_visited(HIT_LIST_t *HIT_list, ULINT count)
{
  HIT_list->FS_seqs_visited += count;
}

void HIT_LIST_count_FS_seq_hit(HIT_LIST_t *HIT_list, ULINT count)
{
  HIT_list->FS_seqs_hits += count;
}

void HIT_LIST_count_seq_visited(HIT_LIST_t *HIT_list, ULINT count)
{
  HIT_list->seqs_visited += count;
}

void HIT_LIST_count_seq_hit(HIT_LIST_t *HIT_list, ULINT count)
{
  HIT_list->seqs_hits += count;
}

void HIT_LIST_stop_timer(HIT_LIST_t *HIT_list)
{
  HIT_list->end_time = time(NULL);
  HIT_list->search_time = HIT_list->end_time - HIT_list->start_time; 
}


/* Printing */
void HIT_LIST_print(HIT_LIST_t *HIT_list, FILE *stream, 
		    HIT_LIST_PRINT_OPT_t *options)
{
  ULINT i;
  ULINT no_frags;

  fastadb_init_Ffrags(HIT_list->s_db, HIT_list->query->len);
  no_frags = fastadb_count_Ffrags(HIT_list->s_db,
				  HIT_list->query->len); 

  fprintf(stream, "**** List of hits ***\n\n");
  fprintf(stream, "Database: %s\n", HIT_list->s_db->db_name);
  fprintf(stream, "Matrix: %s\n", HIT_list->matrix);
  fprintf(stream, "Fragment Length: %ld\n", HIT_list->query->len);
  fprintf(stream, "Total Fragments: %ld\n", no_frags);
  fprintf(stream, "Query: %*s\n", (int) HIT_list->query->len,
	  HIT_list->query->start);
  fprintf(stream, "Cutoff Range: %d\n", HIT_list->range);
  fprintf(stream, "Cutoff Distance Range: %d\n", 
	  HIT_list->converted_range);
  fprintf(stream, "\n");

  fprintf(stream, "Number of bins checked: %ld\n", 
	  HIT_list->FS_seqs_visited);
  fprintf(stream, "Number of bins accepted: %ld\n", 
	  HIT_list->FS_seqs_hits);
  fprintf(stream, "Number of fragments checked: %ld (%.2f %%)\n", 
	  HIT_list->seqs_visited, 
	  (double) HIT_list->seqs_visited / no_frags);
  fprintf(stream, "Number of fragments accepted: %ld (%.2f %%)\n", 
	  HIT_list->seqs_hits,
	  (double) HIT_list->actual_seqs_hits / no_frags);
  fprintf(stream, "Number of hits displayed: %ld\n", 
	  HIT_list->actual_seqs_hits);
  fprintf(stream, "Search time: %.2f sec. \n", HIT_list->search_time);
  fprintf(stream, "\n");
  fprintf(stream, "* Hits *\n\n");
  fprintf(stream, "%6.6s. %10.10s %6.6s %6.6s %.*s %6.6s \n",
	  "Hit #", "Sequence #", "From", "To", 
	  (int) HIT_list->query->len, "Fragment", "Score");

  for (i=0; i <  HIT_list->actual_seqs_hits; i++)
    {
      fprintf(stream, "%6ld. ", i+1);
      SEQ_HIT_print(HIT_list->hits +i, stream, 0); 
    }
}

/* Element Access */
ULINT HIT_LIST_get_seqs_hits(HIT_LIST_t *HIT_list)
{
  return HIT_list->actual_seqs_hits;
}

SEQ_HIT_t *HIT_LIST_get_hit(HIT_LIST_t *HIT_list, ULINT i)
{
  if (i < HIT_list->seqs_hits)
    return HIT_list->hits+i;
  else
    return NULL;
}

void HIT_LIST_set_converted_range(HIT_LIST_t *HIT_list, int crange)
{
  HIT_list->converted_range = crange;  
}

/* Sorting */

void HIT_LIST_sort_decr(HIT_LIST_t *HIT_list)
{
  qsort(HIT_list->hits, HIT_list->actual_seqs_hits, 
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_decr);
}

void HIT_LIST_sort_incr(HIT_LIST_t *HIT_list)
{
  qsort(HIT_list->hits, HIT_list->actual_seqs_hits, 
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_incr);
}

void HIT_LIST_sort_by_sequence(HIT_LIST_t *HIT_list)
{
  qsort(HIT_list->hits, HIT_list->actual_seqs_hits, 
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_by_seq);
}
