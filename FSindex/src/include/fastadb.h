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


#ifndef _FASTADB_H
#define _FASTADB_H

#include "bioseq.h"
#include "misclib.h"

typedef struct
{
  /* General */
  const char *db_name;  /* Name of the database */
  ULINT length;         /* Total Length of Database in residues*/

  /* Sequences */
  ULINT no_seq;         /* Number of sequences in the database */
  BIOSEQ *seq;          /* Array of sequences                  */

  /* Storage Heap */
  char *seq_data;       /* Data heap */
  ULINT seq_data_len;   /* Length of data heap */
  char *deflines;       /* Defline heap */
  ULINT deflines_len;   /* Length of defline heap */

} SEQUENCE_DB;


/* General Purpose Functions */
SEQUENCE_DB *fastadb_open(const char *db_name);
int fastadb_close(SEQUENCE_DB *s_db);

/* Offset <--> char pter conversion functions (as macros) */
char *fastadb_data_pter(SEQUENCE_DB *s_db, ULINT offset);
ULINT fastadb_data_offset(SEQUENCE_DB *s_db, char *s);

/* Ffragment functions */
ULINT fastadb_count_Ffrags(SEQUENCE_DB *s_db, ULINT len);
int fastadb_find_Ffrag_seq(SEQUENCE_DB *s_db, BIOSEQ *frag, 
			   ULINT *seq_id, ULINT *rel_offset);

/* Macro definitions */

#define fastadb_data_pter(db, i) ((db)->seq_data + (i))
#define fastadb_data_offset(db, s) ((s) - (db)->seq_data)
#define fastadb_end_heap(db) ((db)->seq[(db)->no_seq-1].start        \
        + (db)->seq[(db)->no_seq-1].len)

#endif /* #ifndef _FASTADB_H */ 

