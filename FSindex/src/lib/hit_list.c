#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include "expat.h"
#include "partition.h"
#include "hit_list.h"



/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               LIBRARY FUNCTIONS                              ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

#if 0
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

#endif
/********************************************************************/    
/*                                                                  */
/*                    SEQ_HIT module                                */ 
/*                                                                  */
/********************************************************************/    


void SEQ_HIT_print(SEQ_HIT_t *hit, FILE *stream, 
		   SEQ_HIT_PRINT_OPT_t options)
{
  char dots[] = "...";

  if (strlen(hit->subject.id.defline) > 57)
    fprintf(stream, "%-57.57s%3.3s %6.0f\n",
	    hit->subject.id.defline,
	    dots, hit->value);
  else
    fprintf(stream, "%-57.57s%3.3s %6.0f\n",
	    hit->subject.id.defline,
	    "", hit->value);
}

void 
SEQ_HIT_xml(SEQ_HIT_t *hit, FILE *fp, int i,
	    SEQ_HIT_PRINT_OPT_t options)
{
  /* i is indentation */
  xml_open_tag(fp, "hit", i,1);

  xml_open_tag(fp, "seq_id", i+1,0);
  fprintf(fp,"%ld", hit->sequence_id);
  xml_close_tag(fp, "seq_id", 0);
  xml_open_tag(fp, "seq_def", i+1,0);
  fprintf(fp,"%s", hit->subject.id.defline);
  xml_close_tag(fp, "seq_def", 0);
  xml_open_tag(fp, "seq_from", i+1,0);
  fprintf(fp,"%ld", hit->sequence_from);
  xml_close_tag(fp, "seq_from", 0);
  xml_open_tag(fp, "seq_to", i+1,0);
  fprintf(fp,"%ld", hit->sequence_to);
  xml_close_tag(fp, "seq_to", 0);
  xml_open_tag(fp, "seq_len", i+1,0);
  fprintf(fp,"%ld", hit->subject.len);
  xml_close_tag(fp, "seq_len", 0);
  xml_open_tag(fp, "seq_seq", i+1,0);
  fprintf(fp,"%*.*s", (int) hit->subject.len, 
	  (int) hit->subject.len, hit->subject.start);
  xml_close_tag(fp, "seq_seq", 0);
  xml_open_tag(fp, "score", i+1,0);
  fprintf(fp,"%e", (double) hit->value);
  xml_close_tag(fp, "score", 0);

  xml_close_tag(fp, "hit", i);
}


static
void SEQ_HIT_print_alignment(SEQ_HIT_t *hit, FILE *stream,
			     BIOSEQ *query)
{
  fprintf(stream, "%s\n", hit->subject.id.defline);
  fprintf(stream, "        (Length = %ld, ID = %ld)\n\n", 
	  hit->subject.len, hit->sequence_id); 
  fprintf(stream, "Score = %.0f\n", hit->value);
  fprintf(stream, "\n");

  fprintf(stream, "Query: %5ld %*.*s %5ld\n", 
	  1l, (int) hit->subject.len, (int) hit->subject.len,
	  query->start, query->len);

  fprintf(stream, "Sbjct: %5ld %*.*s %5ld\n\n", 
	  hit->sequence_from+1, (int) hit->subject.len, 
	  (int) hit->subject.len, hit->subject.start,
	  (int) hit->sequence_from + hit->subject.len);
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
static inline
HIT_LIST_t *alloc_HL(void)
{
  HIT_LIST_t *HL = callocec(1, sizeof(HIT_LIST_t));
  HL->max_hits = MAX_HITS;
  HL->hits = mallocec(MAX_HITS * sizeof(SEQ_HIT_t));
  return HL;
}

static inline
void reset_HL(HIT_LIST_t *HL)
{
  SEQ_HIT_t *hits = HL->hits;
  ULINT max_hits = HL->max_hits;
  ULINT max_tmp_hits = HL->max_tmp_hits;
  SEQ_HIT_t *tmp_hits = HL->tmp_hits;
  
  /* Reset all values to 0 except array pointers 
     and the sizes allocated. */

  if (HL->p_queue != NULL) {
    free(HL->p_queue->d);
    free(HL->p_queue);
  }

  memset(HL, 0, sizeof(HIT_LIST_t));
  
  HL->hits = hits;
  HL->max_hits = max_hits;
  HL->max_tmp_hits = max_tmp_hits;
  HL->tmp_hits = tmp_hits;

#if 0
  /* Reset arrays of keywords (but not free) */
  HL->oc_size = 0;
  HL->kw_size = 0;
#endif
}

static inline
void init_HL_data(HIT_LIST_t *HL, BIOSEQ *query, 
		  SEQUENCE_DB *s_db, const char *matrix, 
		  void *M, eval_func *efunc, int range)
{
  /* Set the values as per arguments */
  HL->query = *query;

  HL->s_db = s_db;
  HL->db_name = strdup(s_db->db_name);
  HL->db_path = strdup(s_db->db_name);
  HL->db_length = s_db->length;
  HL->db_no_seq = s_db->no_seq;
  HL->db_min_frag_len = s_db->min_len;
  HL->db_max_frag_len = s_db->max_len;
  HL->db_no_frags = s_db->no_frags;

  HL->matrix = strdup(matrix);
  HL->M = M;

  HL->efunc = efunc;
  HL->range = range;
  HL->kNN = -1;
  HL->start_time = time(NULL);
}

HIT_LIST_t *HIT_LIST_create(BIOSEQ *query, SEQUENCE_DB *s_db, 
			    const char *matrix, void *M,
			    eval_func *efunc, int range)
{

  HIT_LIST_t *HL = alloc_HL();
  init_HL_data(HL, query, s_db, matrix, M, efunc, range);
  return HL;
}  

/* Reset list */
void HIT_LIST_reset(HIT_LIST_t *HL, BIOSEQ *query, 
		    SEQUENCE_DB *s_db, const char *matrix, 
		    void *M, eval_func *efunc, int range)
{
  reset_HL(HL);
  init_HL_data(HL, query, s_db, matrix, M, efunc, range);
}

/* Destructor */

/* We destroy all hits stored as well */
void HIT_LIST_destroy(HIT_LIST_t *HL)
{
  free(HL->db_name);
  free(HL->db_path);
  free(HL->index_name);
  free(HL->alphabet);
  free(HL->matrix);
  
  if (HL->p_queue != NULL)
    {
      free(HL->p_queue->d);
      free(HL->p_queue);
    }
  if (HL->alloc_seqs) {
    int i;
    for (i=0; i < HL->actual_seqs_hits; i++) {
      free(HL->hits[i].subject.id.defline);
      free(HL->hits[i].subject.start);
    }
  }
  free(HL->hits);

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
  if (HL->actual_seqs_hits >= HL->max_hits)
    {
      HL->max_hits += MAX_HITS;
      HL->hits = reallocec(HL->hits, HL->max_hits *
				 sizeof(SEQ_HIT_t));
      memset(HL->hits +HL->actual_seqs_hits,0, 
	     (HL->max_hits - HL->actual_seqs_hits)* 
	     sizeof(SEQ_HIT_t));
    }

  HL->hits[HL->actual_seqs_hits].subject = *subject;
  HL->hits[HL->actual_seqs_hits].value = value;

  fastadb_find_Ffrag_seq(HL->s_db, subject, 
     &(HL->hits[HL->actual_seqs_hits].sequence_id),  
     &(HL->hits[HL->actual_seqs_hits].sequence_from)); 
  HL->hits[HL->actual_seqs_hits].sequence_to = 
    HL->hits[HL->actual_seqs_hits].sequence_from +
    HL->hits[HL->actual_seqs_hits].subject.len;
  HL->hits[HL->actual_seqs_hits].subject.id.defline
    = HL->s_db->seq[HL->hits[HL->actual_seqs_hits]
			 .sequence_id].id.defline;
  HL->actual_seqs_hits++;
  HL->accepted = HL->actual_seqs_hits;
}

int HIT_LIST_insert_seq_hit_queue(HIT_LIST_t *HL, 
				  BIOSEQ *subject, 
				  float value,
				  ULINT bin,
				  ULINT pos) 
{
  struct pqueue *p_queue = HL->p_queue;
  SEQ_HIT_t top;
  SEQ_HIT_t temp;
  SEQ_HIT_t new_hit;

  new_hit.subject = *subject;
  new_hit.value = value;

  fastadb_find_Ffrag_seq(HL->s_db, subject, 
     &(new_hit.sequence_id),  
     &(new_hit.sequence_from)); 

  new_hit.subject.id.defline
    = HL->s_db->seq[new_hit.sequence_id].id.defline; 
  new_hit.bin = bin;
  new_hit.pos = pos;

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

  fastadb_init_Ffrags(HL->s_db, HL->query.len);
  no_frags = fastadb_count_Ffrags(HL->s_db,
				  HL->query.len); 

  fprintf(stream, "**** List of hits ****\n\n");
  fprintf(stream, "Query = %s\n", HL->query.id.defline); 
  fprintf(stream, "        %*s\n", (int) HL->query.len,
	  HL->query.start);
  fprintf(stream, "        (%ld letters)\n\n", HL->query.len); 
  
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
			      &HL->query);  
    }
  fprintf(stream, "\n\n");

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
  fprintf(stream, "Fragment Length: %ld\n", HL->query.len);
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


void HIT_LIST_xml(HIT_LIST_t *HL, FILE *fp, 
		  HIT_LIST_PRINT_OPT_t *options)
{
  int i=0; /* Indentation */
  int j;   /* Counter for hits */

  fprintf(fp,"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n");
  xml_open_tag(fp, "FSOutput", i+0,1);

  /* Header */
  xml_open_tag(fp, "FSOutput_program", i+1,0);
  fprintf(fp,"%s", "FSindex_search");
  xml_close_tag(fp, "FSOutput_program", 0);
  xml_open_tag(fp, "FSOutput_version", i+1,0);
  fprintf(fp,"%s", "0.1");
  xml_close_tag(fp, "FSOutput_version", 0);


  /* Database statistics */
  xml_open_tag(fp, "db", i+1,1);
  xml_open_tag(fp, "db_name", i+2,0);
  fprintf(fp,"%s", HL->db_name);
  xml_close_tag(fp, "db_name", 0);
  xml_open_tag(fp, "db_path", i+2,0);
  fprintf(fp,"%s", HL->db_path);
  xml_close_tag(fp, "db_path", 0);
  xml_open_tag(fp, "db_length", i+2,0);
  fprintf(fp,"%ld", HL->db_length);
  xml_close_tag(fp, "db_length", 0);
  xml_open_tag(fp, "db_no_seq", i+2,0);
  fprintf(fp,"%ld", HL->db_no_seq);
  xml_close_tag(fp, "db_no_seq", 0);
  xml_open_tag(fp, "db_min_frag_len", i+2,0);
  fprintf(fp,"%ld", HL->db_min_frag_len);
  xml_close_tag(fp, "db_min_frag_len", 0);
  xml_open_tag(fp, "db_max_frag_len", i+2,0);
  fprintf(fp,"%ld", HL->db_max_frag_len);
  xml_close_tag(fp, "db_max_frag_len", 0);
  xml_open_tag(fp, "db_no_frags", i+2,0);
  fprintf(fp,"%ld", HL->db_no_frags);
  xml_close_tag(fp, "db_no_frags", 0);
  xml_close_tag(fp, "db", i+1);
  
  /* Query parameters */
  xml_open_tag(fp, "Query", i+1,1);
  xml_open_tag(fp, "query_seq", i+2,0);
  fprintf(fp,"%*.*s", (int) HL->query.len, 
	  (int) HL->query.len, HL->query.start); 
  xml_close_tag(fp, "query_seq", 0);
  xml_open_tag(fp, "query_len", i+2,0);
  fprintf(fp,"%ld", HL->query.len);
  xml_close_tag(fp, "query_len", 0);
  xml_open_tag(fp, "query_def", i+2,0);
  fprintf(fp,"%s",  HL->query.id.defline);
  xml_close_tag(fp, "query_def", 0);
  xml_open_tag(fp, "query_type", i+2,0);
  if (HL->kNN == -1)
    fprintf(fp,"%s", "range");
  else
    fprintf(fp,"%s", "kNN");
  xml_close_tag(fp, "query_type", 0);
  xml_open_tag(fp, "similarity_measure_type", i+2,0);
  fprintf(fp,"%s", "distance");
  xml_close_tag(fp, "similarity_measure_type", 0);
  xml_open_tag(fp, "score_matrix", i+2,0);
  fprintf(fp,"%s", HL->matrix);
  xml_close_tag(fp, "score_matrix", 0);
  if (HL->kNN == -1) {
    xml_open_tag(fp, "range", i+2,0);
    fprintf(fp,"%d", HL->range);
    xml_close_tag(fp, "range", 0);
    xml_open_tag(fp, "distance_range", i+2,0);
    fprintf(fp,"%d", HL->converted_range);
    xml_close_tag(fp, "distance_range", 0);
  }
  else {
    xml_open_tag(fp, "NN", i+2,0);
    fprintf(fp,"%d", HL->kNN);
    xml_close_tag(fp, "NN", 0);
  }
  xml_close_tag(fp, "Query", i+1);

  /* List of hits */
  xml_open_tag(fp, "Hit_List", i+1,1);
  xml_open_tag(fp, "Hits_size", i+2,0);
  fprintf(fp,"%ld", HL->actual_seqs_hits );
  xml_close_tag(fp, "Hits_size", 0);
  xml_open_tag(fp, "Hits", i+2,1);
  for (j=0; j < HL->actual_seqs_hits; j++) {
    SEQ_HIT_xml(HL->hits+j, fp, i+3, 0);
  }
  xml_close_tag(fp, "Hits", i+2);
  xml_close_tag(fp, "Hit_List", i+1);

  
  /* Index stats */
  xml_open_tag(fp, "Index_Statistics", i+1,1);
  xml_open_tag(fp, "frag_len", i+2,0);
  fprintf(fp,"%ld", HL->frag_len);
  xml_close_tag(fp, "frag_len", 0);
  xml_open_tag(fp, "alphabet_partitions", i+2,0);
  fprintf(fp,"%s", HL->alphabet);
  xml_close_tag(fp, "alphabet_partitions", 0);
  xml_open_tag(fp, "bins", i+2,0);
  fprintf(fp,"%ld", HL->FS_seqs_total);
  xml_close_tag(fp, "bins", 0);
  xml_open_tag(fp, "seqs_total", i+2,0);
  fprintf(fp,"%ld", HL->index_seqs_total);
  xml_close_tag(fp, "seqs_total", 0);
  xml_open_tag(fp, "useqs_total", i+2,0);
  fprintf(fp,"%ld", HL->useqs_total);
  xml_close_tag(fp, "useqs_total", 0);
  xml_close_tag(fp, "Index_Statistics", i+1);

  /* Search stats */
  xml_open_tag(fp, "Search_Statistics", i+1,1);
  xml_open_tag(fp, "FS_seqs_total", i+2,0);
  fprintf(fp,"%ld", HL->FS_seqs_total);
  xml_close_tag(fp, "FS_seqs_total", 0);
  xml_open_tag(fp, "FS_seqs_visited", i+2,0);
  fprintf(fp,"%ld", HL->FS_seqs_visited);
  xml_close_tag(fp, "FS_seqs_visited", 0);
  xml_open_tag(fp, "FS_seqs_hits", i+2,0);
  fprintf(fp,"%ld", HL->FS_seqs_hits);
  xml_close_tag(fp, "FS_seqs_hits", 0);
  xml_open_tag(fp, "seqs_total", i+2,0);
  fprintf(fp,"%ld", HL->index_seqs_total);
  xml_close_tag(fp, "seqs_total", 0);
  xml_open_tag(fp, "seqs_visited", i+2,0);
  fprintf(fp,"%ld", HL->seqs_visited);
  xml_close_tag(fp, "seqs_visited", 0);
  xml_open_tag(fp, "seqs_hits", i+2,0);
  fprintf(fp,"%ld", HL->seqs_hits);
  xml_close_tag(fp, "seqs_hits", 0);
  xml_open_tag(fp, "useqs_total", i+2,0);
  fprintf(fp,"%ld", HL->useqs_total);
  xml_close_tag(fp, "useqs_total", 0);
  xml_open_tag(fp, "useqs_visited", i+2,0);
  fprintf(fp,"%ld", HL->useqs_visited);
  xml_close_tag(fp, "useqs_visited", 0);
  xml_open_tag(fp, "useqs_hits", i+2,0);
  fprintf(fp,"%ld", HL->seqs_hits);
  xml_close_tag(fp, "useqs_hits", 0);
  xml_open_tag(fp, "start_time", i+2,0);
  fprintf(fp,"%e", HL->start_time);
  xml_close_tag(fp, "start_time", 0);
  xml_open_tag(fp, "end_time", i+2,0);
  fprintf(fp,"%e", HL->end_time);
  xml_close_tag(fp, "end_time", 0);
  xml_open_tag(fp, "search_time", i+2,0);
  fprintf(fp,"%e", HL->search_time);
  xml_close_tag(fp, "search_time", 0);
  xml_close_tag(fp, "Search_Statistics", i+1);

  xml_close_tag(fp, "FSOutput", i);
}



#define BUFSIZE 128
#define ROOT_TAG 51

const
VALID_TAG vtags[] = {
/* 0  */  {"FSOutput", 7, {1,2,3,11,21,32,38}, 51, NOVAL, NOOBJ, 0}, 

/* 1  */  {"FSOutput_program", 0, {}, 0, NOVAL, NOOBJ, 0}, 
/* 2  */  {"FSOutput_version", 0, {}, 0, NOVAL, NOOBJ, 0}, 
/* Needs to be revisited - i.e. add these members to HIT_LIST */
/* 3  */  {"db", 7, {4,5,6,7,8,9,10}, 0, NOVAL, NOOBJ, 0}, 
/* 4  */  {"db_name", 0, {}, 3, CHAR, HIT_LIST, offsetof(HIT_LIST_t, db_name)}, 
/* 5  */  {"db_path", 0, {}, 3, CHAR, HIT_LIST, offsetof(HIT_LIST_t, db_path)}, 
/* 6  */  {"db_length", 0, {}, 3, LONG, HIT_LIST, offsetof(HIT_LIST_t, db_length)}, 
/* 7  */  {"db_no_seq", 0, {}, 3, LONG, HIT_LIST, offsetof(HIT_LIST_t, db_no_seq)}, 
/* 8  */  {"db_min_frag_len", 0, {}, 3, LONG, HIT_LIST, offsetof(HIT_LIST_t, db_min_frag_len)}, 
/* 9  */  {"db_max_frag_len", 0, {}, 3, LONG, HIT_LIST, offsetof(HIT_LIST_t, db_max_frag_len)}, 
/* 10 */  {"db_no_frags", 0, {}, 3, LONG, HIT_LIST, offsetof(HIT_LIST_t, db_no_frags)}, 

/* 11 */  {"Query", 9, {12, 13, 14, 15, 16, 17, 18, 19, 20}, 0, NOVAL, NOOBJ, 0}, 
/* 12 */  {"query_seq", 0, {}, 11, CHAR, HIT_LIST, offsetof(HIT_LIST_t, query) + offsetof(BIOSEQ, start)}, 
/* 13 */  {"query_len", 0, {}, 11, LONG, HIT_LIST, offsetof(HIT_LIST_t, query) + offsetof(BIOSEQ, len)}, 
/* 14 */  {"query_def", 0, {}, 11, CHAR, HIT_LIST, offsetof(HIT_LIST_t, query) + offsetof(BIOSEQ, id)}, 
/* 15 */  {"query_type", 0, {}, 11, NOVAL, NOOBJ, 0}, 
/* 16 */  {"similarity_measure_type", 0, {}, 11, NOVAL, NOOBJ, 0}, 
/* 17 */  {"score_matrix", 0, {}, 11, CHAR, HIT_LIST, offsetof(HIT_LIST_t, matrix)}, 
/* 18 */  {"range", 0, {}, 11, LONG, HIT_LIST, offsetof(HIT_LIST_t, range)}, 
/* 19 */  {"distance_range", 0, {}, 11, LONG, HIT_LIST, offsetof(HIT_LIST_t, converted_range)}, 
/* 20 */  {"NN", 0, {}, 11, LONG, HIT_LIST, offsetof(HIT_LIST_t, kNN)}, 

/* 21 */  {"Hit_List", 2, {22, 23}, 0, NOVAL, NOOBJ, 0},
/* "other version"  {"Hits_size", 0, {}, 21, LONG, HIT_LIST, offsetof(HIT_LIST_t, actual_seqs_hits)} ,*/ 
/* 22 */  {"Hits_size", 0, {}, 21, NOVAL, NOOBJ, 0}, 
/* 23 */  {"Hits", 1, {24}, 21, NOVAL, NOOBJ, 0}, 
/* 24 */  {"hit", 7, {25, 26, 27, 28, 29, 30, 31}, 23, NOVAL, NOOBJ, 0}, 
/* 25 */  {"seq_id", 0, {}, 24, LONG, LIST, offsetof(SEQ_HIT_t, sequence_id)}, 
/* 26 */  {"seq_def", 0, {}, 24, CHAR, LIST, offsetof(SEQ_HIT_t, subject) + offsetof(BIOSEQ, id)}, 
/* 27 */  {"seq_from", 0, {}, 24, LONG, LIST, offsetof(SEQ_HIT_t, sequence_from)}, 
/* 28 */  {"seq_to", 0, {}, 24, LONG, LIST, offsetof(SEQ_HIT_t, sequence_to)}, 
/* 29 */  {"seq_len", 0, {}, 24, LONG, LIST, offsetof(SEQ_HIT_t, subject) + offsetof(BIOSEQ, len)}, 
/* 30 */  {"seq_seq", 0, {}, 24, CHAR, LIST, offsetof(SEQ_HIT_t, subject) + offsetof(BIOSEQ, start)}, 
/* 31 */  {"score", 0, {}, 24, FLOAT, LIST, offsetof(SEQ_HIT_t, value)}, 

/* 32 */  {"Index_Statistics", 5, {33, 34, 35, 36, 37}, 0, NOVAL, NOOBJ, 0},
/* 33 */  {"frag_len", 0, {}, 32, LONG, HIT_LIST, offsetof(HIT_LIST_t, frag_len)},
/* 34 */  {"alphabet_partitions", 0, {}, 32, CHAR, HIT_LIST, offsetof(HIT_LIST_t, alphabet)},
/* 35 */  {"bins", 0, {}, 32, LONG, HIT_LIST, offsetof(HIT_LIST_t, FS_seqs_total)},
/* 36 */  {"seqs_total", 0, {}, 32, LONG, HIT_LIST, offsetof(HIT_LIST_t, index_seqs_total)},
/* 37 */  {"useqs_total", 0, {}, 32, NOVAL, NOOBJ, 0},

/* 38 */  {"Search_Statistics", 12, {39, 40, 41, 42, 43, 44, 45, 46 ,47, 48, 49, 50}, 0, NOVAL, NOOBJ, 0},
/* 39 */  {"FS_seqs_total", 0, {}, 38, LONG, HIT_LIST, offsetof(HIT_LIST_t, FS_seqs_total)},
/* 40 */  {"FS_seqs_visited", 0, {}, 38, LONG, HIT_LIST, offsetof(HIT_LIST_t, FS_seqs_visited)},
/* 41 */  {"FS_seqs_hits", 0, {}, 38, LONG, HIT_LIST, offsetof(HIT_LIST_t, FS_seqs_hits)},
/* 42 */  {"seqs_total", 0, {}, 38, LONG, HIT_LIST, offsetof(HIT_LIST_t, index_seqs_total)},
/* 43 */  {"seqs_visited", 0, {}, 38, LONG, HIT_LIST, offsetof(HIT_LIST_t, seqs_visited)},
/* 44 */  {"seqs_hits", 0, {}, 38, LONG, HIT_LIST, offsetof(HIT_LIST_t, seqs_hits)},
/* 45 */  {"useqs_total", 0, {}, 38, LONG, HIT_LIST, offsetof(HIT_LIST_t, useqs_total)},
/* 46 */  {"useqs_visited", 0, {}, 38, LONG, HIT_LIST, offsetof(HIT_LIST_t, useqs_visited)},
/* 47 */  {"useqs_hits", 0, {}, 38, LONG, HIT_LIST, offsetof(HIT_LIST_t, useqs_hits)},
/* 48 */  {"start_time", 0, {}, 38, DOUBLE, HIT_LIST, offsetof(HIT_LIST_t, start_time)},
/* 49 */  {"end_time", 0, {}, 38, DOUBLE, HIT_LIST, offsetof(HIT_LIST_t, end_time)},
/* 50 */  {"search_time", 0, {}, 38, DOUBLE, HIT_LIST, offsetof(HIT_LIST_t, search_time)},

/* This is root ! - points to 0 */
/* 51 */  {"", 1, {0}, 51, NOVAL, NOOBJ, 0}

};

XML_PARSE_DATA *init_parse_data(void)
{
  XML_PARSE_DATA *PD = callocec(1, sizeof(XML_PARSE_DATA));
  PD->buffer = mallocec(BUFSIZE);
  PD->max_bufsize = BUFSIZE;
  
  /* Allocate HL */
  PD->HL = alloc_HL();
  PD->HL->kNN = -1;

  PD->tags = &vtags[0];
  PD->top = &vtags[ROOT_TAG];
  
  return PD;
}

void destroy_parse_data(XML_PARSE_DATA *PD)
{
  free(PD->buffer);
  free(PD);
}





void
HIT_LIST_XML_StartElementHandler(void *userData,
				 const XML_Char *name,
				 const XML_Char **atts)
{
  int i;
  XML_PARSE_DATA *PD = userData;
  const VALID_TAG *ct = NULL;
  const VALID_TAG *ot = PD->top;

  /* printf("open: %s\n", name); */
  /* Skip if skip flag is set */
  if (PD->skip) {
    PD->skip++;
    return;
  }
  /* Check if we have a valid tag.
     If we do, point to it by top,
     otherwise increment skip flag. */

  for (i=0; i < PD->top->children; i++) {
    ct = PD->tags + PD->top->child[i];
    /* printf("compare with: %s\n", ct->tag); */
    if (!strcmp(name, ct->tag)) {
      PD->top = ct;
      break;
    }
  }
  if (PD->top == ot) {
      PD->skip++;
      return;
  }
  
  PD->opentag = 1;

  if (!strcmp(name, "hit")) {
    if (PD->HL->actual_seqs_hits >= PD->HL->max_hits) {
      PD->HL->max_hits += MAX_HITS;
      PD->HL->hits = reallocec(PD->HL->hits, PD->HL->max_hits *
			   sizeof(SEQ_HIT_t));
      memset(PD->HL->hits +PD->HL->actual_seqs_hits,0, 
	     (PD->HL->max_hits - PD->HL->actual_seqs_hits)* 
	     sizeof(SEQ_HIT_t));
    }
  }
}

void
HIT_LIST_XML_EndElementHandler(void *userData,
			       const XML_Char *name)
{
  XML_PARSE_DATA *PD = userData;
  void *p;
  
  /* printf("close: %s\n", name); */
  /* If we skipped now we go back */
  if (PD->skip) {
    PD->skip--;
    return;
  }
  PD->opentag = 0;
  /* Convert and assign if CDATA */
  if (!PD->top->children) {
    switch (PD->top->otype) 
      {
      case (NOOBJ):
	break;
      case (HIT_LIST):
	p = (char *)PD->HL + PD->top->offset;
	break;
      case (LIST):
	p = (char *)(PD->HL->hits + PD->HL->actual_seqs_hits) + 
	  PD->top->offset;
	break;
      }
    switch (PD->top->vtype)
      {
      case (NOVAL):
	break;
      case (CHAR):
	*((char **)p) = strdup(PD->buffer);
	break;
      case (SHORT):
	*((short *)p) = (short)atol(PD->buffer);
	break;
      case (LONG):
	*((long *)p) = atol(PD->buffer);
	break;
      case (FLOAT):
	*((float *)p) = (float)atof(PD->buffer);
	break;
      case (DOUBLE):
	*((double *)p) = atof(PD->buffer);
	break;
     }
  }

  /* Go up the tree */
  PD->top = PD->tags + PD->top->parent;

  /* Reset buffer */
  PD->bufsize = 0;

  /* If we encounter </hit> */
  if (!strcmp(name, "hit")) {
    PD->HL->actual_seqs_hits++;
  }    
}

void HIT_LIST_XML_CharacterDataHandler(void *userData,
				       const XML_Char *s,
				       int len)
{
  XML_PARSE_DATA *PD = userData;
  /* printf("opentag: %d char <<%*.*s>>\n", PD->opentag, len, len, s); */
  if (PD->skip || !PD->opentag || PD->top->children) return;
  /* Grow buffer */
  if (PD->bufsize + len + 1 >= PD->max_bufsize) {
    PD->max_bufsize += len+1;
    PD->buffer = reallocec(PD->buffer, PD->max_bufsize);
  }
  memcpy(PD->buffer+PD->bufsize,s,len);
  PD->bufsize += len;
  PD->buffer[PD->bufsize]='\0';
}

#define BUFFSIZE        8192
char Buff[BUFFSIZE];

HIT_LIST_t *HIT_LIST_parse_xml(FILE *stream)
{
  HIT_LIST_t *new_HL;
  XML_PARSE_DATA *PD = init_parse_data();
  XML_Parser p = XML_ParserCreate(NULL);

  /* Throw runtime exception */
  if (! p) {
    Throw FSexcept(NO_MEM, "HIT_LIST_parse_xml():"
		   " Couldn't allocate memory for parser.");
  }

  XML_SetElementHandler(p, HIT_LIST_XML_StartElementHandler, 
			HIT_LIST_XML_EndElementHandler);
  XML_SetCharacterDataHandler(p, HIT_LIST_XML_CharacterDataHandler);
  XML_SetUserData(p, PD);

  for (;;) {
    int done;
    int len;

    len = fread(Buff, 1, BUFFSIZE, stream);
    if (ferror(stream)) {
      Throw FSexcept(IO_ERR, 
		     "HIT_LIST_parse_xml(): Read error.");
    }
    done = feof(stream);

    if (XML_Parse(p, Buff, len, done) == XML_STATUS_ERROR) {
      Throw FSexcept(RUNTIME_ERR, "Parse error at line %d:\n%s.",
		     XML_GetCurrentLineNumber(p),
		     XML_ErrorString(XML_GetErrorCode(p)));
    } 

    if (done)
      break;
  }
  new_HL = PD->HL;
  destroy_parse_data(PD);
  XML_ParserFree(p);
  return new_HL;
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
			     ULINT frag_len,
			     ULINT FS_seqs_total,
			     ULINT index_seqs_total,
			     ULINT useqs_total)
{
  HL->index_name = strdup(index_name);
  HL->alphabet = strdup(alphabet);
  HL->frag_len = frag_len;
  HL->FS_seqs_total = FS_seqs_total;
  HL->index_seqs_total = index_seqs_total; 
  HL->useqs_total = useqs_total;
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
	  (*seqs)[*n] = hit->subject;
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
  /* Assuming decreasing order (distances) */
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


#if 0

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

  BIOSEQ *query = &HL->query;
  BIOSEQ subject;
  eval_func *efunc = HL->efunc;
  
  hits = HIT_LIST_get_seqs_hits(HL);



  if (HL->index_name != NULL)
    no_frags = HL->index_seqs_total; 
  else
    no_frags = fastadb_count_Ffrags(HL->s_db, HL->query.len); 

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
      memcpy(subject.start, hit->subject.start, hit->subject.len);
      hit->zvalue = -Zscore(query, &subject, hit->value, 
			    HL->M, efunc);
      hit->pvalue = igamc(shape, (hit->zvalue-Zmin+0.001)*rate);
      hit->evalue = no_frags * hit->pvalue;
    }
  free(subject.start);
}

#endif


