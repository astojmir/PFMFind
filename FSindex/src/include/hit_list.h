#ifndef _HIT_LIST_H
#define _HIT_LIST_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fastadb.h"
#include "smatrix.h"
#include "misclib.h"

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
  int dist;
  int sim;
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
  char **sepn;
  int len;

  /* FS_query */
  BIOSEQ query;
  SCORE_MATRIX *M;
  int conv_type;
  int sim_range;
  int dist_range;
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
  SEQ_HIT_t *hits;

  int alloc_seqs;
} HIT_LIST_t;

/* This type is presently not used */
typedef int HIT_LIST_PRINT_OPT_t;

/* Main constructor */
HIT_LIST_t *HIT_LIST_create(BIOSEQ *query, SEQUENCE_DB *s_db,
			    SCORE_MATRIX *M, int range, int kNN,
			    int conv_type);

/* Reset list */
void HIT_LIST_reset(HIT_LIST_t *HL, BIOSEQ *query, 
		    SEQUENCE_DB *s_db, SCORE_MATRIX *M, 
		    int range, int kNN, int conv_type);

void HIT_LIST_index_data(HIT_LIST_t *HL, char *index_name,
			 char **sepn, int len,
			 ULINT FS_seqs_total,
			 ULINT index_seqs_total,
			 ULINT useqs_total);

/* Destructor */
void HIT_LIST_destroy(HIT_LIST_t *HL);
void HIT_LIST_cleanup(HIT_LIST_t *HL);

/* Counts */
void HIT_LIST_count_FS_seq_visited(HIT_LIST_t *HL, ULINT count);
void HIT_LIST_count_FS_seq_hit(HIT_LIST_t *HL, ULINT count);
void HIT_LIST_count_seq_visited(HIT_LIST_t *HL, ULINT count);
void HIT_LIST_count_seq_hit(HIT_LIST_t *HL, ULINT count);
void HIT_LIST_count_useq_visited(HIT_LIST_t *HL, ULINT count);
void HIT_LIST_count_useq_hit(HIT_LIST_t *HL, ULINT count);

/* Insertion */
void HIT_LIST_insert_hit(HIT_LIST_t *HL, BIOSEQ *subject, 
			 int dist, int sim); 
void HIT_LIST_stop_timer(HIT_LIST_t *HL); 

/* Printing */
void HIT_LIST_print(HIT_LIST_t *HL, FILE *stream, 
		    HIT_LIST_PRINT_OPT_t *options);

#if 0
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
#endif

/* Element Access */
ULINT HIT_LIST_get_seqs_hits(HIT_LIST_t *HL);
SEQ_HIT_t *HIT_LIST_get_hit(HIT_LIST_t *HL, ULINT i);

/* Sorting */
void HIT_LIST_sort_decr(HIT_LIST_t *HL);
void HIT_LIST_sort_incr(HIT_LIST_t *HL);
void HIT_LIST_sort_by_sequence(HIT_LIST_t *HL);


/********************************************************************/    
/*                    Counting Macros                               */ 
/********************************************************************/    

#define HIT_LIST_count_FS_seq_visited(HL, count) \
  (HL)->FS_seqs_visited += (count)

#define HIT_LIST_count_FS_seq_hit(HL, count)     \
  (HL)->FS_seqs_hits += (count)

#define HIT_LIST_count_seq_visited(HL, count)    \
  (HL)->seqs_visited += (count)

#define HIT_LIST_count_seq_hit(HL, count)        \
  (HL)->seqs_hits += (count)

#define HIT_LIST_count_useq_visited(HL, count)   \
  (HL)->useqs_visited += (count);

#define HIT_LIST_count_useq_hit(HL, count)       \
  (HL)->useqs_hits += (count);

#endif
