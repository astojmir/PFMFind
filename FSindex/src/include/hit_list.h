#ifndef _HIT_LIST_H
#define _HIT_LIST_H

#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fastadb.h"
#ifdef USE_MPATROL
#include <mpatrol.h>
#endif

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               PROTOTYPES                                     ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    


/********************************************************************/    
/*                                                                  */
/*                    SEQ_HIT module                                */ 
/*                                                                  */
/********************************************************************/    

typedef struct 
{
  BIOSEQ *subject;
  ULINT sequence_id;
  ULINT sequence_from;
  float value;
  double pvalue;
} SEQ_HIT_t;

/* This type is presently not used */
typedef int SEQ_HIT_PRINT_OPT_t;

void SEQ_HIT_print(SEQ_HIT_t *hit, FILE *stream, 
		   SEQ_HIT_PRINT_OPT_t options);


/********************************************************************/    
/*                                                                  */
/*                    HIT_LIST module                               */ 
/*                                                                  */
/********************************************************************/    


  
typedef enum {SIMILARITY, QUASI_METRIC} SIMILARITY_MEASURE_t;
typedef enum {RANGE, K_NN, NET, FS_RANGE} SEARCH_TYPE_t;

typedef struct
{
  BIOSEQ *query;
  ULINT frag_len;
  SEQUENCE_DB *s_db;
  const char *matrix;
  int range;
  int converted_range;

  const char *index_name;
  const char *alphabet;
  ULINT FS_seqs_total;
  ULINT FS_seqs_visited;
  ULINT FS_seqs_hits;
  ULINT index_seqs_total;
  ULINT seqs_visited;
  ULINT seqs_hits;
  double start_time;
  double end_time;
  double search_time;

  ULINT max_hits;
  ULINT actual_seqs_hits;
  SEQ_HIT_t *hits;  
} HIT_LIST_t;

/* This type is presently not used */
typedef int HIT_LIST_PRINT_OPT_t;

/* Main constructor */
HIT_LIST_t *HIT_LIST_create(BIOSEQ *query, SEQUENCE_DB *s_db, 
			    const char *matrix, int range);
/* Reset list */
void HIT_LIST_reset(HIT_LIST_t *HIT_list, BIOSEQ *query, 
		    SEQUENCE_DB *s_db, const char *matrix, int range); 

/* Destructor */
void HIT_LIST_destroy(HIT_LIST_t *HIT_list);

/* Insertion */
void HIT_LIST_count_FS_seq_visited(HIT_LIST_t *HIT_list, ULINT count);
void HIT_LIST_count_FS_seq_hit(HIT_LIST_t *HIT_list, ULINT count);
void HIT_LIST_count_seq_visited(HIT_LIST_t *HIT_list, ULINT count);
void HIT_LIST_count_seq_hit(HIT_LIST_t *HIT_list, ULINT count);

void HIT_LIST_insert_seq_hit(HIT_LIST_t *HIT_list, BIOSEQ *subject, 
			     float value); 
void HIT_LIST_stop_timer(HIT_LIST_t *HIT_list); 

/* Printing */
void HIT_LIST_print(HIT_LIST_t *HIT_list, FILE *stream, 
		    HIT_LIST_PRINT_OPT_t *options);

/* Element Access */
ULINT HIT_LIST_get_seqs_hits(HIT_LIST_t *HIT_list);
SEQ_HIT_t *HIT_LIST_get_hit(HIT_LIST_t *HIT_list, ULINT i);
void HIT_LIST_set_converted_range(HIT_LIST_t *HIT_list, int crange);
void HIT_LIST_set_index_data(HIT_LIST_t *HIT_list, 
			     const char *index_name,
			     const char *alphabet,
			     ULINT FS_seqs_total,
			     ULINT index_seqs_total);
void HIT_LIST_get_hit_seqs(HIT_LIST_t *HIT_list, BIOSEQ **seqs,
			   int cutoff, int *n, int *max_n);

/* Sorting */
void HIT_LIST_sort_decr(HIT_LIST_t *HIT_list);
void HIT_LIST_sort_incr(HIT_LIST_t *HIT_list);
void HIT_LIST_sort_by_sequence(HIT_LIST_t *HIT_list);

/* Add p-values */
void HIT_LIST_Gaussian_pvalues(HIT_LIST_t *HT, double mean, 
			       double var);


#endif
