#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include "avl.h"
#include "partition.h"
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

static
int split_line(char *line,       /* line string */ 
	       int n,            /* max number of characters */
	       char c)           /* separator char */
/* returns number of chars split */
{
  int l = strlen(line);
  char *cp;

  if (l <= n) return l;
  for (cp = line + n - 1; cp >= line; cp--)
    if (*cp == c)
      break;

  if (cp > line) 
    return cp - line;
  else return l;
}


/********************************************************************/    
/*                                                                  */
/*                    SEQ_HIT module                                */ 
/*                                                                  */
/********************************************************************/    


void SEQ_HIT_print(SEQ_HIT_t *hit, FILE *stream, 
		   SEQ_HIT_PRINT_OPT_t options)
{
  char dots[] = "...";
#if 0  
  if (strlen(hit->subject->id.defline) > 57)
    fprintf(stream, "%-57.57s%3.3s %6.0f %5.2e\n",
	    hit->subject->id.defline,
	    dots, hit->value, hit->pvalue);
  else
    fprintf(stream, "%-57.57s%3.3s %6.0f %5.2e\n",
	    hit->subject->id.defline,
	    "", hit->value, hit->pvalue);
#endif
  if (strlen(hit->subject->id.defline) > 57)
    fprintf(stream, "%-57.57s%3.3s %6.0f %5.2e\n",
	    hit->subject->id.defline,
	    dots, hit->value, hit->evalue);
  else
    fprintf(stream, "%-57.57s%3.3s %6.0f %5.2e\n",
	    hit->subject->id.defline,
	    "", hit->value, hit->evalue);
}

static
void SEQ_HIT_print_alignment(SEQ_HIT_t *hit, FILE *stream,
			     BIOSEQ *query, KW_INDEX *KWI)
{
  int i;
  int j;
  int c;
  int k;
  int m;

  fprintf(stream, "%s\n", hit->subject->id.defline);
  fprintf(stream, "        (Length = %ld, ID = %ld)\n\n", 
	  hit->subject->len, hit->sequence_id); 
  fprintf(stream, "Score = %.0f, Z-score = %.2f, P-value = %.2e "
	  "E-value = %.2e\n",
	  hit->value, hit->zvalue, hit->pvalue, hit->evalue);

  if (KWI != NULL)
    fprintf(stream, "C-ratio = %.2f, KW-score = %.2f\n", hit->cratio,
	    hit->kw_score);
  fprintf(stream, "\n");

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
  if (KWI != NULL)
    {
      if (hit->oc_cluster > 0)
	{
	  m = split_line(KWI->cld_h + KWI->cld[hit->oc_cluster], 60, ' ');
	  fprintf(stream, "O-cluster: (%d) %.*s\n", hit->oc_cluster, m, 
		  KWI->cld_h + KWI->cld[hit->oc_cluster]);
	  while ((k = split_line(KWI->cld_h + KWI->cld[hit->oc_cluster] + m,
				 60, ' ')) > 0)
	    {
	      fprintf(stream, "%20s %.*s\n", "", k,
		      KWI->cld_h + KWI->cld[hit->oc_cluster] + m);
	      m += k;
	    }
	}

      if (KWI->seq2kw1[hit->sequence_id] > 0)
	{
	  fprintf(stream, "Keywords: ");
	  c = 10;
	  for (i = KWI->seq2kw0[hit->sequence_id]; 
	       i < KWI->seq2kw0[hit->sequence_id] + KWI->seq2kw1[hit->sequence_id];
	       i++)
	    {
	      j = KWI->skw_h[i];
	      if (c + strlen(KWI->kwd_h + KWI->kwd[j]) + 2 > 80)
		{
		  fprintf(stream, "\n          ");
		  c = 10;
		}
	      fprintf(stream, "%s; ", KWI->kwd_h + KWI->kwd[j]);
	      c += strlen(KWI->kwd_h + KWI->kwd[j]) + 2;
	    }
	  fprintf(stream, "\n\n");
	}

    }
  fprintf(stream, "----------------------------------------"
	  "----------------------------------------\n");
}

static
int SEQ_HIT_cmp_decr(const void *S1, const void *S2)
{
  const SEQ_HIT_t *T1 = S1;
  const SEQ_HIT_t *T2 = S2;

  /* Want to use for sorting in decreasing order */
  return (int) ((T2)->value - (T1)->value);
}

static
int SEQ_HIT_cmp_incr(const void *S1, const void *S2)
{
  const SEQ_HIT_t *T1 = S1;
  const SEQ_HIT_t *T2 = S2;

  return (int) ((T1)->value - (T2)->value);
}

static
int SEQ_HIT_cmp_by_seq(const void *S1, const void *S2)
{
  const SEQ_HIT_t *T1 = S1;
  const SEQ_HIT_t *T2 = S2;

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

static
int SEQ_HIT_cmp_by_oc(const void *S1, const void *S2)
{
  const SEQ_HIT_t *T1 = S1;
  const SEQ_HIT_t *T2 = S2;

  if (T1->oc_cluster > T2->oc_cluster)
    return 1;
  else if (T1->oc_cluster < T2->oc_cluster)
    return -1;
  else if (T1->sequence_id > T2->sequence_id)
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

static
int SEQ_HIT_cmp_by_evalue(const void *S1, const void *S2)
{
  const SEQ_HIT_t *T1 = S1;
  const SEQ_HIT_t *T2 = S2;

  if (T1->evalue < T2->evalue)
    return 1;
  else if (T1->evalue > T2->evalue)
    return -1;
  else
    return 0;
}

static
int SEQ_HIT_cmp_by_cratio(const void *S1, const void *S2)
{
  const SEQ_HIT_t *T1 = S1;
  const SEQ_HIT_t *T2 = S2;

  if (T1->cratio > T2->cratio)
    return 1;
  else if (T1->cratio < T2->cratio)
    return -1;
  else
    return 0;
}

static
int SEQ_HIT_cmp_by_kwscore(const void *S1, const void *S2)
{
  const SEQ_HIT_t *T1 = S1;
  const SEQ_HIT_t *T2 = S2;

  if (T1->kw_score > T2->kw_score)
    return 1;
  else if (T1->kw_score < T2->kw_score)
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
			    const char *matrix, void *M,
			    eval_func *efunc, int range)
{
  HIT_LIST_t *HL = callocec(1, sizeof(HIT_LIST_t));
  HL->query = query;
  HL->s_db = s_db;
  HL->matrix = matrix;
  HL->M = M;
  HL->efunc = efunc;
  HL->range = range;
  HL->start_time = time(NULL);
  HL->max_hits = MAX_HITS;
  HL->hits = mallocec(MAX_HITS * sizeof(SEQ_HIT_t));
  HL->kNN = -1;
  return HL;
}  

/* Reset list */
void HIT_LIST_reset(HIT_LIST_t *HL, BIOSEQ *query, 
		    SEQUENCE_DB *s_db, const char *matrix, 
		    void *M, eval_func *efunc, int range)
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
  HL->M = M;
  HL->efunc = efunc;
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

  /* Reset arrays of keywords (but not free) */
  HL->oc_size = 0;
  HL->kw_size = 0;
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
#if 0
  WRD_CNTR_clear(HL->oc);
  WRD_CNTR_clear(HL->kw);
#endif
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
  memset(HL->hits +HL->actual_seqs_hits,0, sizeof(SEQ_HIT_t));
  HL->hits[HL->actual_seqs_hits].subject = new_subject;
  HL->hits[HL->actual_seqs_hits].value = value;
  HL->hits[HL->actual_seqs_hits].rejected = 0;

  /* TO DO: get sequence_id, sequence_from from database */
  fastadb_find_Ffrag_seq(HL->s_db, new_subject, 
     &(HL->hits[HL->actual_seqs_hits].sequence_id),  
     &(HL->hits[HL->actual_seqs_hits].sequence_from)); 

  HL->hits[HL->actual_seqs_hits].subject->id.defline
    = HL->s_db->seq[HL->hits[HL->actual_seqs_hits]
			 .sequence_id].id.defline;
  HL->actual_seqs_hits++;
  HL->accepted = HL->actual_seqs_hits;
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
  int k;
  int m;
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
  fprintf(stream, "%76s E\n","");
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
			      HL->query, HL->KWI);  
    }
  fprintf(stream, "\n\n");

  if (HL->KWI != NULL) {
  fprintf(stream, "%25s **** Orthologous clusters ****\n\n", "");
  fprintf(stream, "%5s %7s %5s %5s %5s %s\n", "", "Cluster", "Hits", "Total",
	  "Ratio", "Cluster Name");
  for (i=0; i <  HL->oc_size; i++)
    {
      m = split_line(HL->KWI->cld_h + HL->KWI->cld[HL->oc[i].item], 48, ' ');
      fprintf(stream, "%4d. %7d %5d %5d %5.2f %.*s\n", i+1, HL->oc[i].item,
	      HL->oc[i].count, HL->KWI->clf[HL->oc[i].item], 
	      HL->oc[i].score, m, 
	      HL->KWI->cld_h + HL->KWI->cld[HL->oc[i].item]);
      while ((k = split_line(HL->KWI->cld_h + HL->KWI->cld[HL->oc[i].item] + m,
			     48, ' ')) > 0)
	{
	  fprintf(stream, "%31s %.*s\n", "", k,
		  HL->KWI->cld_h + HL->KWI->cld[HL->oc[i].item] + m);
	  m += k;
	}
    }

  fprintf(stream, "\n\n");
    
  fprintf(stream, "%25s **** Keywords ****\n\n", "");
  fprintf(stream, "%5s %5s %5s %6s %s\n", "", "ID", "Hits",
	  "Score", "Keyword");
  for (i=0; i <  HL->kw_size; i++)
    {
      fprintf(stream, "%4d. %5d %5d %6.2f %s\n", 
	      i+1, HL->kw[i].item, HL->kw[i].count, 
	      HL->kw[i].score, HL->KWI->kwd_h + HL->KWI->kwd[HL->kw[i].item]);
    }
  fprintf(stream, "\n\n");
  }
  }
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
  fprintf(stream, "\n");

  fprintf(stream, "Gamma distribution shape: %.4f\n", HL->shape);
  fprintf(stream, "Gamma distribution rate: %.4f\n", HL->rate);
  fprintf(stream, "Gamma distribution beta: %.4f\n", HL->Zmin);
  fprintf(stream, "\n");

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
	  (double) HL->seqs_hits * 100 / no_frags);
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
  if (i >= HL->actual_seqs_hits)
    Throw FSexcept(INDEX_OUT_OF_RANGE, 
		   "HIT_LIST_get_hit(): Index out of range.");
  return HL->hits+i;
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

#if 0
  if (FSS->pfunc == profile_process_bin ||
      FSS->pfunc == profile_kNN_process_bin)
    efunc = profile_eval;
  else
    efunc = smatrix_eval;
#endif



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
#if 0
      /* Some garbage collection */
      free(hit->subject);
      hit->subject = NULL;
#endif
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

void HIT_LIST_sort_oc(HIT_LIST_t *HL)
{
  qsort(HL->hits, HL->actual_seqs_hits, 
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_by_oc);
}

void HIT_LIST_sort_evalue(HIT_LIST_t *HL, int offset)
{
  qsort(HL->hits + offset, HL->actual_seqs_hits - offset,
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_by_evalue);
}

void HIT_LIST_sort_cratio(HIT_LIST_t *HL, int offset)
{
  qsort(HL->hits + offset, HL->actual_seqs_hits - offset,
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_by_cratio);
}

void HIT_LIST_sort_kwscore(HIT_LIST_t *HL, int offset)
{
  qsort(HL->hits + offset, HL->actual_seqs_hits - offset,
	sizeof(SEQ_HIT_t), SEQ_HIT_cmp_by_kwscore);
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


/********************************************************************/ 
/*         Z-score functions                                        */
/********************************************************************/ 
static inline
double Zscore(BIOSEQ *query, BIOSEQ *subject, double value, 
	      void *M, eval_func *efunc)
{
  int j;
  int k;
  int r;
  int val;
  char tmp;
  double x;
  double xx;
  double mean;
  double sd;
  const int p = 100;
  
  x=0.0;
  xx=0.0;     
  for (j=0; j < p; j++)
    {
      for (k=subject->len-1; k > 0; k--)
	{
	  r = lrand48() % subject->len;
	  tmp = subject->start[r];
	  subject->start[r] = subject->start[k];
	  subject->start[k] = tmp;
	}
      val = (double) efunc(M, query, subject);
      x += val;
      xx += val*val;
    }
  mean = x/p;
  sd = sqrt(xx/p - mean*mean);
  if (sd == 0.0)
    return 0.0;
  else
    return (value - mean)/sd;
}

static
void Z_gamma_params(double *shape, double *rate, double *Zmin,
		    SEQUENCE_DB *s_db, void *M, eval_func *efunc, 
		    BIOSEQ *query, int no_samples)
{
  ULINT no_frags;
  BIOSEQ subject;
  ULINT frag_len = query->len;
  int j;
  double value;
  double *Z = mallocec(no_samples * sizeof(double));
  ULINT rand_frag;
  FS_SEQ_t FS_seq;
  
  double arith_mean = 0.0;
  double geom_mean = 0.0;
  double MM;
  double num;
  double den;
  FS_PARTITION_t *ptable; 

  /* Any ptable corresponding to the whole protein alphabet is OK */
  ptable = FS_PARTITION_create("STAN#ILVM#KR#EDQ#WFYH#GPC", '#'); 

  *Zmin = 1000000;
  fastadb_get_nofrags(s_db, &no_frags, frag_len, frag_len);
  for (j=0; j < no_samples; j++)
    {
      /* Get a random database fragment */
      do 
	{
	  rand_frag = lrand48()%no_frags;
	  fastadb_get_frag(s_db, &subject, rand_frag,
			   frag_len, frag_len);
	}
      while (!BIOSEQ_2_FS_SEQ(&subject, ptable, &FS_seq));
      
      /* Calculate distances and add to histograms*/
      value = (double) efunc(M, query, &subject);
      
      /* Get Z-Score */
      Z[j] = -Zscore(query, &subject, value, M, efunc);
      if (Z[j] < *Zmin)
	*Zmin = Z[j];
    }

  /* We need to shift Z-scores by their min. value so that they all are
     gt. 0 */
     
  for (j=0; j < no_samples; j++)
    {
      arith_mean += (Z[j]-*Zmin+0.001);
      geom_mean += log(Z[j]-*Zmin+0.001);
    }


  /* gam.apprxFIT
     fitting parameters shape(k) and rate (1/theta) in gamma distribution
     method is taken from the encyclopedia of biostatistics: 
     gamma distibution, page 292

     Adapted from S code found at
     http://www.biostat.wustl.edu/mailinglists/s-news/200106/msg00152.html
  */

  arith_mean /= no_samples;
  geom_mean = exp(geom_mean/no_samples);
  MM = log(arith_mean/geom_mean);

  if(MM >= 0 && MM <= 0.5772) 
    *shape = (0.5000876 + 0.1648852 * MM - 0.0544276 * MM * MM)/MM;
  else if(MM >= 0.5772 && MM <= 17) 
    {
      num = 8.898919 + 9.05995 * MM + 0.9775373 * MM * MM;
      den = MM * (17.79728 + 11.968477 * MM + MM * MM);
      *shape = num/den;
    }
  else if(MM > 17)
    *shape = 1/MM;

  *rate = *shape/arith_mean;
  free(Z);
  FS_PARTITION_destroy(ptable);  
}

#if 0
/* Original S code */

function(x ) 
{ 
#gam.apprxFIT
#fitting parameters shape(k) and rate (1/theta) in gamma distribution
#method is taken from the encyclopedia of biostatistics: gamma distibution, page 292
http://www.biostat.wustl.edu/mailinglists/s-news/200106/msg00152.html
n <- length(x)
arith.mean <- mean(x)
geom.mean <- exp(sum(log(x))/n)


M <- log(arith.mean/geom.mean)
if(M >= 0 & M <= 0.5772) {
est.shape <- (0.5000876 + 0.1648852 * M - 0.0544276 * M^2)/M
}
if(M >= 0.5772 & M <= 17) {
num <- 8.898919 + 9.05995 * M + 0.9775373 * M^2
den <- M * (17.79728 + 11.968477 * M + M^2)
est.shape <- num/den
}
if(M > 17)
est.shape <- 1/M
est.rate <- est.shape/arith.mean
list(est.shape = est.shape, est.rate = est.rate)
} 

#endif 


extern double igamc( double, double );

void HIT_LIST_Zscores(HIT_LIST_t *HL)
{
  int i;
  int hits;
  SEQ_HIT_t *hit;
  double shape;
  double rate;
  double Zmin;
  ULINT no_frags;

  BIOSEQ *query = HL->query;
  BIOSEQ subject;
  eval_func *efunc = HL->efunc;
  
  hits = HIT_LIST_get_seqs_hits(HL);



  if (HL->index_name != NULL)
    no_frags = HL->index_seqs_total; 
  else
    no_frags = fastadb_count_Ffrags(HL->s_db, HL->query->len); 

  srand48(time(NULL)); 

  subject.start = mallocec(query->len+1);
  subject.start[query->len] = '\0';
  subject.len = query->len;

  /* Get the parameters of the Gamma distribution */

  Z_gamma_params(&shape, &rate, &Zmin, HL->s_db, HL->M, efunc, 
		 query, 5000);

  HL->shape = shape;
  HL->rate = rate;
  HL->Zmin = Zmin;
  
  /* Evaluate Z-scores and p-values for hits */
  for (i=0; i < hits; i++)
    {
      hit = HIT_LIST_get_hit(HL, i);
      memcpy(subject.start, hit->subject->start, hit->subject->len);
      hit->zvalue = -Zscore(query, &subject, hit->value, 
			    HL->M, efunc);
      hit->pvalue = igamc(shape, (hit->zvalue-Zmin+0.001)*rate);
      hit->evalue = no_frags * hit->pvalue;
    }
  free(subject.start);
}



