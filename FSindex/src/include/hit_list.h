#ifndef _HIT_LIST_H
#define _HIT_LIST_H

#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "expat.h"
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
  BIOSEQ subject;
  ULINT sequence_id;
  ULINT sequence_from;
  ULINT sequence_to;
  ULINT bin;
  ULINT pos;
  float value;
} SEQ_HIT_t;

/* This type is presently not used */
typedef int SEQ_HIT_PRINT_OPT_t;

void SEQ_HIT_print(SEQ_HIT_t *hit, FILE *stream, 
		   SEQ_HIT_PRINT_OPT_t options);

void 
SEQ_HIT_xml(SEQ_HIT_t *hit, FILE *fp, int i,
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

  /* db_stats */
  SEQUENCE_DB *s_db;
  char *db_name;
  char *db_path;
  ULINT db_length;
  ULINT db_no_seq;
  ULINT db_min_frag_len;
  ULINT db_max_frag_len;
  ULINT db_no_frags;

  /* index_stats */
  char *index_name;
  char *alphabet;

  /* FS_query */
  ULINT frag_len;
  SEARCH_TYPE_t srch_type;
  BIOSEQ query;
  char *matrix;
  void *M;
  eval_func *efunc;
  int range;
  int converted_range;
  int kNN;
 
  /* FS_search_stats */
  ULINT FS_seqs_total;
  ULINT FS_seqs_visited;
  ULINT FS_seqs_hits;
  ULINT index_seqs_total;
  ULINT seqs_visited;
  ULINT seqs_hits;
  ULINT useqs_total;
  ULINT useqs_visited;
  ULINT useqs_hits;
  double start_time;
  double end_time;
  double search_time;

  /* FS_hits */
  ULINT max_hits;
  ULINT actual_seqs_hits;
  ULINT accepted;
  SEQ_HIT_t *hits;

  /* Priority queue */
  struct pqueue *p_queue;
  ULINT max_tmp_hits;
  ULINT no_tmp_hits;
  SEQ_HIT_t *tmp_hits;

  
  int alloc_seqs;
} HIT_LIST_t;

/* This type is presently not used */
typedef int HIT_LIST_PRINT_OPT_t;

/* Main constructor */
HIT_LIST_t *HIT_LIST_create(BIOSEQ *query, SEQUENCE_DB *s_db, 
			    const char *matrix, void *M,
			    eval_func *efunc, int range);
/* Reset list */
void HIT_LIST_reset(HIT_LIST_t *HL, BIOSEQ *query, 
		    SEQUENCE_DB *s_db, const char *matrix, 
		    void *M, eval_func *efunc, int range);

/* Destructor */
void HIT_LIST_destroy(HIT_LIST_t *HIT_list);

/* Insertion */
void HIT_LIST_count_FS_seq_visited(HIT_LIST_t *HIT_list, ULINT count);
void HIT_LIST_count_FS_seq_hit(HIT_LIST_t *HIT_list, ULINT count);
void HIT_LIST_count_seq_visited(HIT_LIST_t *HIT_list, ULINT count);
void HIT_LIST_count_seq_hit(HIT_LIST_t *HIT_list, ULINT count);
void HIT_LIST_count_useq_visited(HIT_LIST_t *HIT_list, ULINT count);
void HIT_LIST_count_useq_hit(HIT_LIST_t *HIT_list, ULINT count);

void HIT_LIST_insert_seq_hit(HIT_LIST_t *HIT_list, BIOSEQ *subject, 
			     float value); 
int HIT_LIST_insert_seq_hit_queue(HIT_LIST_t *HL, BIOSEQ *subject, 
				  float value, ULINT bin, ULINT pos); 
void HIT_LIST_stop_timer(HIT_LIST_t *HIT_list); 

/* Printing */
void HIT_LIST_print(HIT_LIST_t *HIT_list, FILE *stream, 
		    HIT_LIST_PRINT_OPT_t *options);

/* XML routines */
#define MAX_CHILDREN 15
typedef enum {NOVAL, CHAR, SHORT, LONG, FLOAT, DOUBLE} val_type;
typedef enum {NOOBJ, HIT_LIST, LIST} obj_type;

struct valid_tag
{
  const char *tag;
  int children;
  int child[MAX_CHILDREN];
  int parent;
  val_type vtype;
  obj_type otype;
  size_t offset;
};

typedef struct valid_tag VALID_TAG;

typedef struct
{
  const VALID_TAG *tags;
  const VALID_TAG *top;
  char *buffer;
  int bufsize;
  int max_bufsize;
  HIT_LIST_t *HL;
  int skip;
  int opentag;
} XML_PARSE_DATA;

void HIT_LIST_xml(HIT_LIST_t *HL, FILE *fp, 
		  HIT_LIST_PRINT_OPT_t *options);

XML_PARSE_DATA *init_parse_data(void);
void destroy_parse_data(XML_PARSE_DATA *PD);

void
HIT_LIST_XML_StartElementHandler(void *userData,
				 const XML_Char *name,
				 const XML_Char **atts);
void
HIT_LIST_XML_EndElementHandler(void *userData,
			       const XML_Char *name);
void HIT_LIST_XML_CharacterDataHandler(void *userData,
				       const XML_Char *s,
				       int len);

HIT_LIST_t *HIT_LIST_parse_xml(FILE *stream);

/* Element Access */
ULINT HIT_LIST_get_seqs_hits(HIT_LIST_t *HIT_list);
SEQ_HIT_t *HIT_LIST_get_hit(HIT_LIST_t *HIT_list, ULINT i);
void HIT_LIST_set_kNN(HIT_LIST_t *HIT_list, ULINT kNN);
void HIT_LIST_set_converted_range(HIT_LIST_t *HIT_list, int crange);
void HIT_LIST_set_index_data(HIT_LIST_t *HIT_list, 
			     const char *index_name,
			     const char *alphabet,
			     ULINT frag_len,
			     ULINT FS_seqs_total,
			     ULINT index_seqs_total,
			     ULINT useqs_total);

void HIT_LIST_get_hit_seqs(HIT_LIST_t *HIT_list, BIOSEQ **seqs,
			   int cutoff, int *n, int *max_n);

/* Sorting */
void HIT_LIST_sort_decr(HIT_LIST_t *HIT_list);
void HIT_LIST_sort_incr(HIT_LIST_t *HIT_list);
void HIT_LIST_sort_by_sequence(HIT_LIST_t *HIT_list);
void HIT_LIST_sort_kNN(HIT_LIST_t *HL);

#if 0 /* Z-scores */
void HIT_LIST_Zscores(HIT_LIST_t *HL);
#endif

#endif
