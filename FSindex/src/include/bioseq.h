#ifndef _BIOSEQ_H
#define _BIOSEQ_H

#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef USE_MPATROL
#include <mpatrol.h>
#endif

typedef union 
{
  ULINT id_num;
  char *defline;
} seq_id_t;

typedef enum {NO, YES} yn_bool_t; 

typedef struct
  {
    seq_id_t id;
    ULINT len;
    char *start;
  } BIOSEQ; 

typedef int eval_func(void *M, BIOSEQ *query, BIOSEQ *subject);

/* ************* bioseq functions ********************* */
int cmp_bioseq(const void *S1, const void *S2);
int bioseq_parse(BIOSEQ *seq, char *filename, yn_bool_t defline);
int bioseq_get_frag(BIOSEQ *seq, BIOSEQ *frag, ULINT start, ULINT
		    length, ULINT id_no);
void string2seq(BIOSEQ *seq, char *string);
void bioseq_random(BIOSEQ *seq, ULINT seq_len, char *alphabet,
		   ULINT a_len);
void bioseq_seq2string(BIOSEQ *seq, char *string, ULINT from, 
		       ULINT to);
BIOSEQ *bioseq_copy(BIOSEQ *seq);

#endif /* #ifndef _BIOSEQ_H */ 
