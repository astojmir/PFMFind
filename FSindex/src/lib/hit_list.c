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
    fprintf(stream, "%-57.57s%3.3s %6d\n",
	    hit->subject.id.defline,
	    dots, hit->dist);
  else
    fprintf(stream, "%-57.57s%3.3s %6d\n",
	    hit->subject.id.defline,
	    "", hit->dist);
}

#if 0
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
#endif

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
  HL->db_min_frag_len = s_db->min_len;
  HL->db_max_frag_len = s_db->max_len;
  HL->db_no_frags = s_db->no_frags;

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

#if 0
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
#endif



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

