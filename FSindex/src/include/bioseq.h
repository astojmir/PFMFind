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


#ifndef _BIOSEQ_H
#define _BIOSEQ_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "misclib.h"

typedef union 
{
  ULINT id_num;
  const char *defline;
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
