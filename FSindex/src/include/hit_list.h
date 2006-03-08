/*
 * Copyright (C) 2004-2006 Victoria University of Wellington
 *
 * This file is part of the PFMFind module.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2,
 * or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


#ifndef _HIT_LIST_H
#define _HIT_LIST_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fastadb.h"
#include "smatrix.h"
#include "misclib.h"

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
