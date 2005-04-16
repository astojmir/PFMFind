#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include "partition.h"
#include "hit_list.h"
#include "misclib.h"

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

  if (strlen(hit->subject.id.defline) > 57)
    fprintf(stream, "%-57.57s%3.3s %6d\n",
	    hit->subject.id.defline,
	    dots, hit->dist);
  else
    fprintf(stream, "%-57.57s%3.3s %6d\n",
	    hit->subject.id.defline,
	    "", hit->dist);
}

static
void SEQ_HIT_print_alignment(SEQ_HIT_t *hit, FILE *stream,
			     BIOSEQ *query)
{
  fprintf(stream, "%s\n", hit->subject.id.defline);
  fprintf(stream, "        (Length = %ld, ID = %ld)\n\n", 
	  hit->subject.len, hit->sequence_id); 
  fprintf(stream, "Distance = %d\n", hit->dist);
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
  return (int) ((T2)->dist - (T1)->dist);
}

static
int SEQ_HIT_cmp_incr(const void *S1, const void *S2)
{
  const SEQ_HIT_t *T1 = S1;
  const SEQ_HIT_t *T2 = S2;

  return (int) ((T1)->dist - (T2)->dist);
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
  int i;

  free(HL->db_name);
  free(HL->db_path);
  free(HL->index_name);

  free((char *)HL->query.id.defline);
  free(HL->query.start);

  if (HL->sepn != NULL) { 
    for (i=0; i < HL->len; i++)
      free(HL->sepn[i]);
    free(HL->sepn);
  }
  if (HL->alloc_seqs) {
    for (i=0; i < HL->actual_seqs_hits; i++) {
      free((char *)HL->hits[i].subject.id.defline);
      free(HL->hits[i].subject.start);
    }
  }

  /* Reset all values to 0 except array pointers 
     and the sizes allocated. */

  memset(HL, 0, sizeof(HIT_LIST_t));
  
  HL->hits = hits;
  HL->max_hits = max_hits;
}

static inline
void init_HL_data(HIT_LIST_t *HL, BIOSEQ *query, 
		  SEQUENCE_DB *s_db, SCORE_MATRIX *M, 
		  int range, int kNN, int conv_type)
{
  HL->s_db = s_db;
  HL->db_name = strdup(s_db->db_name);
  HL->db_path = strdup(s_db->db_name);
  HL->db_length = s_db->length;
  HL->db_no_seq = s_db->no_seq;

  /* TO DO: remove these */
  HL->db_min_frag_len = 0;
  HL->db_max_frag_len = 0;
  HL->db_no_frags = s_db->length;

  HL->query.len = query->len;
  HL->query.id.defline = strdup(query->id.defline);
  HL->query.start = mallocec(query->len+1);
  memcpy(HL->query.start, query->start, query->len);
  HL->query.start[query->len] = '\0';

  HL->M = M;
  if (kNN == -1) {
    HL->conv_type = conv_type;
    M->set_conv_type(M, conv_type);
    HL->dist_range = M->range_conv(M, query->start, query->len, range);
    if (range != HL->dist_range)
      HL->sim_range = range;
    else
      HL->sim_range = 0;
  }
  HL->kNN = kNN;
  HL->start_time = time(NULL);
}

HIT_LIST_t *HIT_LIST_create(BIOSEQ *query, SEQUENCE_DB *s_db,
			    SCORE_MATRIX *M, int range, int kNN,
			    int conv_type)
{

  HIT_LIST_t *HL;
  HL = alloc_HL();
  init_HL_data(HL, query, s_db, M, range, kNN, conv_type);
  return HL;
}  

/* Reset list */
void HIT_LIST_reset(HIT_LIST_t *HL, BIOSEQ *query, 
		    SEQUENCE_DB *s_db, SCORE_MATRIX *M, 
		    int range, int kNN, int conv_type)
{
  reset_HL(HL);
  init_HL_data(HL, query, s_db, M, range, kNN, conv_type);
}

void HIT_LIST_index_data(HIT_LIST_t *HL, char *index_name,
			 char **sepn, int len,
			 ULINT FS_seqs_total,
			 ULINT index_seqs_total,
			 ULINT useqs_total)
{
  int i;
  HL->index_name = strdup(index_name);
  HL->len = len;
  HL->sepn = mallocec(len*sizeof(char *));
  for (i=0; i < len; i++)
    HL->sepn[i] = strdup(sepn[i]);
  HL->FS_seqs_total = FS_seqs_total;
  HL->index_seqs_total = index_seqs_total; 
  HL->useqs_total = useqs_total;
}

/* Destructors */
/* We destroy all hits stored as well */
void HIT_LIST_cleanup(HIT_LIST_t *HL)
{
  int i;
  free(HL->db_name);
  free(HL->db_path);
  free(HL->index_name);

  free((char *)HL->query.id.defline);
  free(HL->query.start);

  if (HL->sepn != NULL) { 
    for (i=0; i < HL->len; i++)
      free(HL->sepn[i]);
    free(HL->sepn);
  }
  if (HL->alloc_seqs) {
    for (i=0; i < HL->actual_seqs_hits; i++) {
      free((char *)HL->hits[i].subject.id.defline);
      free(HL->hits[i].subject.start);
    }
  }
  free(HL->hits);
}

void HIT_LIST_destroy(HIT_LIST_t *HL)
{
  HIT_LIST_cleanup(HL);
  free(HL);
}

/* Insertion */
void HIT_LIST_insert_hit(HIT_LIST_t *HL, BIOSEQ *subject, 
			     int dist, int sim) 
{
  ULINT *n = &(HL->actual_seqs_hits);
  ULINT *N = &(HL->max_hits);
  SEQ_HIT_t *h;

  if (*n >= *N) {
    *N += MAX_HITS;
    HL->hits = reallocec(HL->hits, *N*sizeof(SEQ_HIT_t));
    memset(HL->hits+*n, 0, (*N - *n)* sizeof(SEQ_HIT_t));
  }
  h = HL->hits+(*n)++;
  
  h->subject = *subject;
  h->dist = dist;
  h->sim = sim;

  fastadb_find_Ffrag_seq(HL->s_db, subject, &(h->sequence_id),  
			 &(h->sequence_from)); 
  h->sequence_to = h->sequence_from + h->subject.len;
  h->subject.id.defline = 
    HL->s_db->seq[h->sequence_id].id.defline;
}


/* Printing */
void HIT_LIST_print(HIT_LIST_t *HL, FILE *stream, 
		    HIT_LIST_PRINT_OPT_t *options)
{
  int i;
  ULINT no_frags;

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
      fprintf(stream, "Alphabet Partitions:\n");
      for (i=0; i < HL->len; i++) 
	fprintf(stream,"%s\n", HL->sepn[i]); 
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

  /* fprintf(stream, "Matrix: %s\n", HL->matrix); */
  fprintf(stream, "Similarity cutoff value: %d\n", HL->sim_range);
  fprintf(stream, "Distance cutoff value: %d\n", 
	  HL->dist_range);
  fprintf(stream, "\n");

  if (HL->FS_seqs_visited)
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

void HIT_LIST_stop_timer(HIT_LIST_t *HL)
{
  HL->end_time = time(NULL);
  HL->search_time = HL->end_time - HL->start_time; 
}

