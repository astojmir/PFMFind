#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
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
  char dots[] = "...";
  
  if (strlen(hit->subject->id.defline) > 57)
    fprintf(stream, "%-57.57s%3.3s %6.0f %5.2e\n",
	    hit->subject->id.defline,
	    dots, hit->value, hit->pvalue);
  else
    fprintf(stream, "%-57.57s%3.3s %6.0f %5.2e\n",
	    hit->subject->id.defline,
	    "", hit->value, hit->pvalue);
}

static
void SEQ_HIT_print_alignment(SEQ_HIT_t *hit, FILE *stream,
			     BIOSEQ *query,
			     SEQ_HIT_PRINT_OPT_t options)
{
  fprintf(stream, "%s\n", hit->subject->id.defline);
  fprintf(stream, "        (Length = %ld, ID = %ld)\n\n", 
	  hit->subject->len, hit->sequence_id); 
  fprintf(stream, " Score = %5.0f, P-value = %5.2e\n", hit->value,
	  hit->pvalue);
#if 0
  fprintf(stream, "\n");
#endif
  fprintf(stream, "Query: %5ld %*.*s %5ld\n", 
	  1l, (int) hit->subject->len, (int) hit->subject->len,
	  query->start, query->len);
#if 0
  fprintf(stream, "\n");
#endif
  fprintf(stream, "Sbjct: %5ld %*.*s %5ld\n\n", 
	  hit->sequence_from+1, (int) hit->subject->len, 
	  (int) hit->subject->len, hit->subject->start,
	  (int) hit->sequence_from + hit->subject->len);
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

  HIT_list->hits[HIT_list->actual_seqs_hits].subject->id.defline
    = HIT_list->s_db->seq[HIT_list->hits[HIT_list->actual_seqs_hits]
			 .sequence_id].id.defline;
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
  int i;
  ULINT no_frags;

  fastadb_init_Ffrags(HIT_list->s_db, HIT_list->query->len);
  no_frags = fastadb_count_Ffrags(HIT_list->s_db,
				  HIT_list->query->len); 

  fprintf(stream, "**** List of hits ****\n\n");
  fprintf(stream, "Query = %s\n", HIT_list->query->id.defline); 
  fprintf(stream, "        %*s\n", (int) HIT_list->query->len,
	  HIT_list->query->start);
  fprintf(stream, "        (%ld letters)\n\n", HIT_list->query->len); 
  
  fprintf(stream, "Database: %s\n\n", HIT_list->s_db->db_name);
  
  if (HIT_list->actual_seqs_hits == 0)
    fprintf(stream, "%25s\n --- No significant hits found ---\n\n",
	    "");
  else {
  fprintf(stream, "%25s **** Significant hits ****\n\n", "");
  fprintf(stream, "%76s P\n","");
  fprintf(stream, "%67s  Score Value\n", "");

  for (i=0; i <  HIT_list->actual_seqs_hits; i++)
    {
      fprintf(stream, "%5d. ", i+1);
      SEQ_HIT_print(HIT_list->hits +i, stream, 0); 
    }
  fprintf(stream, "\n");
  fprintf(stream, "%25s **** Alignments ****\n\n", "");
  for (i=0; i <  HIT_list->actual_seqs_hits; i++)
    {
      fprintf(stream, "%d. ", i+1);
      SEQ_HIT_print_alignment(HIT_list->hits +i, stream, 
			      HIT_list->query, 0);  
    }
  fprintf(stream, "\n\n");}

  fprintf(stream, "**** Statistics ****\n\n");
  fprintf(stream, "Database: %s\n", HIT_list->s_db->db_name);
  fprintf(stream, "Number of letters in database: %ld\n", 
	  HIT_list->s_db->length);
  fprintf(stream, "Number of sequences in database: %ld\n", 
	  HIT_list->s_db->no_seq);
  fprintf(stream, "\n");
  
  if (HIT_list->index_name != NULL)
    {
      fprintf(stream, "Index: %s\n", HIT_list->index_name);
      fprintf(stream, "Alphabet Partitions: %s\n",
	      HIT_list->alphabet); 
    }
  fprintf(stream, "Fragment Length: %ld\n", HIT_list->query->len);
  fprintf(stream, "Number of fragments in database: %ld\n", no_frags); 
  if (HIT_list->index_name != NULL)
    {
      fprintf(stream, "Number of fragments in index: %ld\n", 
	      HIT_list->index_seqs_total); 
      fprintf(stream, "Number of bins in index: %ld\n", 
	      HIT_list->FS_seqs_total);
    } 
  fprintf(stream, "\n");

  fprintf(stream, "Matrix: %s\n", HIT_list->matrix);
  fprintf(stream, "Score cutoff value: %d\n", HIT_list->range);
  fprintf(stream, "Distance cutoff value: %d\n", 
	  HIT_list->converted_range);
  if (HIT_list->index_name != NULL)
    {
      fprintf(stream, "Number of generated neighbour bins: %ld"
	      " (%.2f %%)\n", HIT_list->FS_seqs_visited,
	      (double) HIT_list->FS_seqs_visited * 100 /
	      HIT_list->FS_seqs_total);
      fprintf(stream, "Number of accepted neighbour bins: %ld "
	      "(%.2f %%)\n", HIT_list->FS_seqs_hits,
	      (double) HIT_list->FS_seqs_hits * 100 /
	      HIT_list->FS_seqs_total);
      fprintf(stream, "Accepted out of generated neighbour bins:"
	      " %.2f %%\n", (double) HIT_list->FS_seqs_hits * 100 /
	      HIT_list->FS_seqs_visited);
    }
  fprintf(stream, "Number of checked fragments: %ld (%.2f %%)\n", 
	  HIT_list->seqs_visited, 
	  (double) HIT_list->seqs_visited * 100/ no_frags);
  fprintf(stream, "Number of accepted fragments: %ld (%.2f %%)\n",  
	  HIT_list->seqs_hits,
	  (double) HIT_list->actual_seqs_hits * 100 / no_frags);
  fprintf(stream, "Number of displayed hits: %ld\n", 
	  HIT_list->actual_seqs_hits);
  fprintf(stream, "Search time: %.2f sec. \n", HIT_list->search_time);
  fprintf(stream, "\n");
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

void HIT_LIST_set_index_data(HIT_LIST_t *HIT_list, 
			     const char *index_name,
			     const char *alphabet,
			     ULINT FS_seqs_total,
			     ULINT index_seqs_total)
{
  HIT_list->index_name = index_name;
  HIT_list->alphabet = alphabet;
  HIT_list->FS_seqs_total = FS_seqs_total;
  HIT_list->index_seqs_total = index_seqs_total; 
}

void HIT_LIST_get_hit_seqs(HIT_LIST_t *HIT_list, BIOSEQ **seqs,
			   int cutoff, int *n, int *max_n)
{
  int i;
  SEQ_HIT_t *hit;

  if (*max_n < HIT_list->actual_seqs_hits)
    {
      *max_n = HIT_list->actual_seqs_hits;
      *seqs = reallocec(*seqs, *max_n * sizeof(BIOSEQ));
    }

  for (i=0, *n=0, hit=HIT_list->hits; i < HIT_list->actual_seqs_hits; 
       i++, hit++)
    {
      if (hit->value >= cutoff)
	{
	  (*seqs)[*n] = *(hit->subject);
	  (*n)++;
	}
      /* Some garbage collection */
      free(hit->subject);
      hit->subject = NULL;
    }
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


/* Add p-values */
void HIT_LIST_Gaussian_pvalues(HIT_LIST_t *HT, double mean, 
			       double var)
{
  int i;
  double sd = sqrt(var);
  double Z;

  for (i=0; i <  HT->actual_seqs_hits; i++)
    {
      Z = (HT->hits[i].value - mean)/sd;
      if (Z < 0) Z = -Z;
      HT->hits[i].pvalue = 0.5 * erfc(Z/M_SQRT2);
    }
}
