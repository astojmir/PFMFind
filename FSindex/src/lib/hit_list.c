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

#define PQDATUM SEQ_HIT_t 
#define PQPRIO(p) (p.value) 

#include "pqueue.h"

/* Main constructor */
HIT_LIST_t *HIT_LIST_create(BIOSEQ *query, SEQUENCE_DB *s_db, 
			    const char *matrix, int range)
{
  HIT_LIST_t *HL = callocec(1, sizeof(HIT_LIST_t));
  HL->query = query;
  HL->s_db = s_db;
  HL->matrix = matrix;
  HL->range = range;
  HL->start_time = time(NULL);
  HL->max_hits = MAX_HITS;
  HL->hits = mallocec(MAX_HITS * sizeof(SEQ_HIT_t));
  HL->kNN = -1;

  return HL;
}  

/* Reset list */
void HIT_LIST_reset(HIT_LIST_t *HL, BIOSEQ *query, 
		    SEQUENCE_DB *s_db, const char *matrix, int range)
{
  ULINT i;
  /* hits and tmp_hits are not reallocated */
  /* Everything is set to 0 except max_hits and max_tmp_hits */

  for (i=0; i <  HL->actual_seqs_hits; i++)
    free(HL->hits[i].subject);

  HL->query = query;
  HL->frag_len = 0;
  HL->s_db = s_db;
  HL->matrix = matrix;
  HL->range = range;
  HL->converted_range = 0;
  HL->kNN = -1;
  HL->index_name = NULL;
  HL->alphabet = NULL;
  HL->FS_seqs_total = 0;
  HL->FS_seqs_visited = 0;
  HL->FS_seqs_hits = 0;
  HL->index_seqs_total = 0;
  HL->seqs_visited = 0;
  HL->seqs_hits = 0;
  HL->useqs_visited = 0;
  HL->useqs_hits = 0;
  HL->start_time = time(NULL);
  HL->end_time = 0;
  HL->search_time = 0;

  HL->actual_seqs_hits = 0;
  if (HL->p_queue != NULL)
    {
      free(HL->p_queue->d);
      free(HL->p_queue);
      HL->p_queue = NULL;
    }
  HL->no_tmp_hits = 0;
}

/* Destructor */

/* We destroy all hits stored as well */
void HIT_LIST_destroy(HIT_LIST_t *HL)
{
  ULINT i;
  for (i=0; i <  HL->actual_seqs_hits; i++)
    free(HL->hits[i].subject);
  if (HL->p_queue != NULL)
    {
      free(HL->p_queue->d);
      free(HL->p_queue);
    }
  free(HL->hits);
  free(HL->tmp_hits);
  free(HL);
}

/* Insertion */
void HIT_LIST_insert_seq_hit(HIT_LIST_t *HL, BIOSEQ *subject, 
			     float value) 
{
  BIOSEQ *new_subject = mallocec(sizeof(BIOSEQ));
  memcpy(new_subject, subject, sizeof(BIOSEQ));

  if (HL->actual_seqs_hits >= HL->max_hits)
    {
      HL->max_hits += MAX_HITS;
      HL->hits = reallocec(HL->hits, HL->max_hits *
				 sizeof(SEQ_HIT_t));
    }
  HL->hits[HL->actual_seqs_hits].subject = new_subject;
  HL->hits[HL->actual_seqs_hits].value = value;

  /* TO DO: get sequence_id, sequence_from from database */
  fastadb_find_Ffrag_seq(HL->s_db, new_subject, 
     &(HL->hits[HL->actual_seqs_hits].sequence_id),  
     &(HL->hits[HL->actual_seqs_hits].sequence_from)); 

  HL->hits[HL->actual_seqs_hits].subject->id.defline
    = HL->s_db->seq[HL->hits[HL->actual_seqs_hits]
			 .sequence_id].id.defline;
  HL->actual_seqs_hits++;
}

int HIT_LIST_insert_seq_hit_queue(HIT_LIST_t *HL, 
				   BIOSEQ *subject, 
				   float value) 
{
  BIOSEQ *new_subject = mallocec(sizeof(BIOSEQ));
  struct pqueue *p_queue = HL->p_queue;
  SEQ_HIT_t top;
  SEQ_HIT_t temp;
  SEQ_HIT_t new_hit;
  int i;

 *new_subject = *subject;

  new_hit.subject = new_subject;
  new_hit.value = value;

  fastadb_find_Ffrag_seq(HL->s_db, new_subject, 
     &(new_hit.sequence_id),  
     &(new_hit.sequence_from)); 

  new_hit.subject->id.defline
    = HL->s_db->seq[new_hit.sequence_id].id.defline; 

  if (HL->actual_seqs_hits < HL->kNN)
    {
      pqinsert(p_queue, new_hit);
      
      if (++HL->actual_seqs_hits == HL->kNN)
	{
	  pqpeek(p_queue, &top);
	  HL->range = (int) top.value;
	}
      else
	HL->range = INT_MAX;
    }
  else if (value < HL->range)
    {
      pqremove(p_queue, &temp);
      pqinsert(p_queue, new_hit);

      pqpeek(p_queue, &top);
      HL->range = (int) top.value;

      if (temp.value == top.value)
	{
	  if (HL->no_tmp_hits >= HL->max_tmp_hits)
	    {
	      HL->max_tmp_hits += HL->kNN;
	      HL->tmp_hits = reallocec(HL->tmp_hits, 
			HL->max_tmp_hits * sizeof(SEQ_HIT_t));
	    }
	  HL->tmp_hits[HL->no_tmp_hits] = temp;
	  HL->no_tmp_hits++;
	  HL->actual_seqs_hits++; 
	}
      else
	{
	  free(temp.subject);
	  for (i=0; i < HL->no_tmp_hits; i++)
	    free(HL->tmp_hits[i].subject);
	  HL->no_tmp_hits = 0;
	  HL->actual_seqs_hits = HL->kNN; 
	}
    }
  else /* if (value == HL->range) */
    {
      if (HL->no_tmp_hits >= HL->max_tmp_hits)
	{
	  HL->max_tmp_hits += HL->kNN;
	  HL->tmp_hits = reallocec(HL->tmp_hits, 
	    HL->max_tmp_hits * sizeof(SEQ_HIT_t));
	}
      HL->tmp_hits[HL->no_tmp_hits] = new_hit;
      HL->no_tmp_hits++;
      HL->actual_seqs_hits++; 
    }

  return HL->range;
}


void HIT_LIST_count_FS_seq_visited(HIT_LIST_t *HL, ULINT count)
{
  HL->FS_seqs_visited += count;
}

void HIT_LIST_count_FS_seq_hit(HIT_LIST_t *HL, ULINT count)
{
  HL->FS_seqs_hits += count;
}

void HIT_LIST_count_seq_visited(HIT_LIST_t *HL, ULINT count)
{
  HL->seqs_visited += count;
}

void HIT_LIST_count_seq_hit(HIT_LIST_t *HL, ULINT count)
{
  HL->seqs_hits += count;
}

void HIT_LIST_count_useq_visited(HIT_LIST_t *HL, ULINT count)
{
  HL->useqs_visited += count;
}

void HIT_LIST_count_useq_hit(HIT_LIST_t *HL, ULINT count)
{
  HL->useqs_hits += count;
}


void HIT_LIST_stop_timer(HIT_LIST_t *HL)
{
  HL->end_time = time(NULL);
  HL->search_time = HL->end_time - HL->start_time; 
}


/* Printing */
void HIT_LIST_print(HIT_LIST_t *HL, FILE *stream, 
		    HIT_LIST_PRINT_OPT_t *options)
{
  int i;
  ULINT no_frags;

  fastadb_init_Ffrags(HL->s_db, HL->query->len);
  no_frags = fastadb_count_Ffrags(HL->s_db,
				  HL->query->len); 

  fprintf(stream, "**** List of hits ****\n\n");
  fprintf(stream, "Query = %s\n", HL->query->id.defline); 
  fprintf(stream, "        %*s\n", (int) HL->query->len,
	  HL->query->start);
  fprintf(stream, "        (%ld letters)\n\n", HL->query->len); 
  
  fprintf(stream, "Database: %s\n\n", HL->s_db->db_name);
  
  if (HL->actual_seqs_hits == 0)
    fprintf(stream, "%25s\n --- No significant hits found ---\n\n",
	    "");
  else {
  fprintf(stream, "%25s **** Significant hits ****\n\n", "");
  fprintf(stream, "%76s P\n","");
  fprintf(stream, "%67s  Score Value\n", "");

  for (i=0; i <  HL->actual_seqs_hits; i++)
    {
      fprintf(stream, "%5d. ", i+1);
      SEQ_HIT_print(HL->hits+i, stream, 0); 
    }
  fprintf(stream, "\n");
  fprintf(stream, "%25s **** Alignments ****\n\n", "");
  for (i=0; i <  HL->actual_seqs_hits; i++)
    {
      fprintf(stream, "%d. ", i+1);
      SEQ_HIT_print_alignment(HL->hits +i, stream, 
			      HL->query, 0);  
    }
  fprintf(stream, "\n\n");}

  fprintf(stream, "**** Statistics ****\n\n");
  fprintf(stream, "Database: %s\n", HL->s_db->db_name);
  fprintf(stream, "Number of letters in database: %ld\n", 
	  HL->s_db->length);
  fprintf(stream, "Number of sequences in database: %ld\n", 
	  HL->s_db->no_seq);
  fprintf(stream, "\n");
  
  if (HL->index_name != NULL)
    {
      fprintf(stream, "Index: %s\n", HL->index_name);
      fprintf(stream, "Alphabet Partitions: %s\n",
	      HL->alphabet); 
    }
  fprintf(stream, "Fragment Length: %ld\n", HL->query->len);
  fprintf(stream, "Number of fragments in database: %ld\n", no_frags); 
  if (HL->index_name != NULL)
    {
      fprintf(stream, "Number of fragments in index: %ld\n", 
	      HL->index_seqs_total); 
      fprintf(stream, "Number of bins in index: %ld\n", 
	      HL->FS_seqs_total);
    } 
  fprintf(stream, "\n");

  fprintf(stream, "Matrix: %s\n", HL->matrix);
  fprintf(stream, "Score cutoff value: %d\n", HL->range);
  fprintf(stream, "Distance cutoff value: %d\n", 
	  HL->converted_range);
  if (HL->index_name != NULL)
    {
      fprintf(stream, "Number of generated neighbour bins: %ld"
	      " (%.2f %%)\n", HL->FS_seqs_visited,
	      (double) HL->FS_seqs_visited * 100 /
	      HL->FS_seqs_total);
      fprintf(stream, "Number of accepted neighbour bins: %ld "
	      "(%.2f %%)\n", HL->FS_seqs_hits,
	      (double) HL->FS_seqs_hits * 100 /
	      HL->FS_seqs_total);
      fprintf(stream, "Accepted out of generated neighbour bins:"
	      " %.2f %%\n", (double) HL->FS_seqs_hits * 100 /
	      HL->FS_seqs_visited);
    }
  fprintf(stream, "Number of checked fragments: %ld (%.2f %%)\n", 
	  HL->seqs_visited, 
	  (double) HL->seqs_visited * 100/ no_frags);
  fprintf(stream, "Number of accepted fragments: %ld (%.2f %%)\n",  
	  HL->seqs_hits,
	  (double) HL->actual_seqs_hits * 100 / no_frags);
  fprintf(stream, "Number of checked distinct fragments: %ld\n", 
	  HL->useqs_visited);
  fprintf(stream, "Number of accepted distinct fragments: %ld\n",  
	  HL->useqs_hits);

  fprintf(stream, "Number of displayed hits: %ld\n", 
	  HL->actual_seqs_hits);
  fprintf(stream, "Search time: %.2f sec. \n", HL->search_time);
  fprintf(stream, "\n");
}


/* Element Access */
ULINT HIT_LIST_get_seqs_hits(HIT_LIST_t *HL)
{
  return HL->actual_seqs_hits;
}

SEQ_HIT_t *HIT_LIST_get_hit(HIT_LIST_t *HL, ULINT i)
{
  if (i < HL->seqs_hits)
    return HL->hits+i;
  else
    return NULL;
}

void HIT_LIST_set_kNN(HIT_LIST_t *HL, ULINT kNN)
{
  HL->kNN = kNN;
  HL->p_queue = pqinit(NULL, kNN);
  HL->max_tmp_hits = kNN;
  HL->no_tmp_hits = 0;
  HL->tmp_hits = 
    reallocec(HL->tmp_hits, kNN * sizeof(SEQ_HIT_t));
}


void HIT_LIST_set_converted_range(HIT_LIST_t *HL, int crange)
{
  HL->converted_range = crange;  
}

void HIT_LIST_set_index_data(HIT_LIST_t *HL, 
			     const char *index_name,
			     const char *alphabet,
			     ULINT FS_seqs_total,
			     ULINT index_seqs_total)
{
  HL->index_name = index_name;
  HL->alphabet = alphabet;
  HL->FS_seqs_total = FS_seqs_total;
  HL->index_seqs_total = index_seqs_total; 
}

void HIT_LIST_get_hit_seqs(HIT_LIST_t *HL, BIOSEQ **seqs,
			   int cutoff, int *n, int *max_n)
{
  int i;
  SEQ_HIT_t *hit;

  if (*max_n < HL->actual_seqs_hits)
    {
      *max_n = HL->actual_seqs_hits;
      *seqs = reallocec(*seqs, *max_n * sizeof(BIOSEQ));
    }

  for (i=0, *n=0, hit=HL->hits; i < HL->actual_seqs_hits; 
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

void HIT_LIST_sort_decr(HIT_LIST_t *HL)
{
  qsort(HL->hits, HL->actual_seqs_hits, 
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_decr);
}

void HIT_LIST_sort_incr(HIT_LIST_t *HL)
{
  qsort(HL->hits, HL->actual_seqs_hits, 
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_incr);
}

void HIT_LIST_sort_by_sequence(HIT_LIST_t *HL)
{
  qsort(HL->hits, HL->actual_seqs_hits, 
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_by_seq);
}

void HIT_LIST_sort_kNN(HIT_LIST_t *HL)
{
  /* Assuming increasing order (distances) */
  int i;
  struct pqueue *p_queue = HL->p_queue;

  if (HL->actual_seqs_hits > HL->max_hits)
    {
      HL->max_hits = HL->actual_seqs_hits;
      HL->hits = reallocec(HL->hits, 
	HL->max_hits * sizeof(SEQ_HIT_t));
    }
  for (i=0; i < HL->kNN; i++)
    {
      pqremove(p_queue, HL->hits + i);
    }
  memcpy(HL->hits + i, HL->tmp_hits, 
	 HL->no_tmp_hits * sizeof(SEQ_HIT_t));
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
