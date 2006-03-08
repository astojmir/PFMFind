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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fastadb.h"
#include "misclib.h"

/* Allocate memory in 8 Mb blocks */
#ifndef ALLOC_BLOCK_SIZE
#define ALLOC_BLOCK_SIZE (1 << 23)
#endif
#define BUF_SIZE 512

static inline
int fastadb_put_line(char *buffer, char **heap, ULINT *len, 
		     ULINT *max_len)
{
  int l;
  l = strlen(buffer);

  if (*len + l + 2 >= *max_len) {
    *max_len += ALLOC_BLOCK_SIZE;
    *heap = reallocec(*heap, *max_len); 
  }
  if (buffer[l-1] == '\n')
    l--;
  
  if (l > 0)
    memcpy(*heap+*len, buffer, l);

  *len += l;
  return l;  
}

SEQUENCE_DB *fastadb_open(const char *db_name)
{
  FILE *fp;
  SEQUENCE_DB *s_db;  
  int c, l, i;
  static char buffer[BUF_SIZE+1];
  char *s1, *s2;
  ULINT data_size = ALLOC_BLOCK_SIZE;
  ULINT desc_size = ALLOC_BLOCK_SIZE;

  buffer[BUF_SIZE] = '\0'; /* just in case */

  if (db_name == NULL) {
    fp = stdin;
    db_name = "stdin";
  }
  else if ((fp = fopen(db_name, "r")) == NULL)
    Throw FSexcept(FOPEN_ERR, 
		   "fastadb_open(): Could not open fasta database %s",
		   db_name);


  s_db = callocec(1, sizeof(SEQUENCE_DB));

  s_db->db_name = strdup(db_name);

  /* Allocate the storage heap */
  s_db->seq_data = mallocec(data_size);
  s_db->deflines = mallocec(desc_size);

  /* Load sequences into memory */
  while(!feof(fp)) {
    /* First character must be '>' because this is a new sequence. */

    /* Description */
    while (1) {
      fgets(buffer, BUF_SIZE, fp);
      l = fastadb_put_line(buffer+1, &(s_db->deflines), 
			   &(s_db->deflines_len),
			   &desc_size);
      if (buffer[l+1] == '\n') break;
    }
    s_db->deflines[s_db->deflines_len++] = '\0';

    /* Sequence */
    while (1) {
      if ((c = fgetc(fp)) == EOF) break;
      if (c == '>') { 
	ungetc(c, fp);
	break;
      }

      buffer[0] = c;
      fgets(buffer+1, BUF_SIZE, fp);
      l = fastadb_put_line(buffer, &(s_db->seq_data), 
			   &(s_db->seq_data_len),
			   &data_size); 
    }
    s_db->seq_data[s_db->seq_data_len++] = '\0';

    s_db->no_seq++;
  }

  if (db_name != NULL) fclose(fp);

  /* Trim memory segments - we add 30 chars more to seq_data
     so that we can have fast fragment comparison near the end
  */
  s_db->seq_data = reallocec(s_db->seq_data, s_db->seq_data_len+30);
  s_db->deflines = reallocec(s_db->deflines, s_db->deflines_len);

  /* Allocate sequences */
  s_db->seq = mallocec(s_db->no_seq * sizeof(BIOSEQ));
  
  /* Write offsets */
  s1 = s_db->seq_data;
  s2 = s_db->deflines;
  for (i=0; i < s_db->no_seq; i++) {
    s_db->seq[i].start = s1;
    s_db->seq[i].len = strlen(s1);
    s1 += s_db->seq[i].len + 1;
    s_db->length += s_db->seq[i].len;
    
    s_db->seq[i].id.defline = s2;
    s2 += strlen(s2) + 1;
  }

  return s_db;
}

int fastadb_close(SEQUENCE_DB *s_db)
{
  if (s_db == NULL) return 1;
  free((char *)s_db->db_name);
  free(s_db->seq_data);
  free(s_db->deflines);
  free(s_db->seq);
  free(s_db);
  return 1;
}

ULINT fastadb_count_Ffrags(SEQUENCE_DB *s_db, ULINT len)
{
  ULINT i;
  ULINT S = 0;
  for (i=0; i < s_db->no_seq; i++)
    if (s_db->seq[i].len >= len) 
      S += s_db->seq[i].len - len + 1;
  return S;
}

int fastadb_find_Ffrag_seq(SEQUENCE_DB *s_db, BIOSEQ *frag, 
			   ULINT *seq_id, ULINT *rel_offset) 
{
  ULINT a = 0;
  ULINT b = s_db->no_seq - 1;
  ULINT k = b/2; 

  if (frag->start > s_db->seq[s_db->no_seq - 1].start) {
    *seq_id = s_db->no_seq - 1;
    *rel_offset = frag->start - s_db->seq[s_db->no_seq - 1].start;
    return 1;
  }

  while (a<=b) {
    k=(a+b)/2;
    if (s_db->seq[k].start <= frag->start) {
      if (s_db->seq[k+1].start > frag->start) {
	*seq_id = k;
	*rel_offset = frag->start - s_db->seq[k].start;
	return 1;
      }
      else a=k+1;
    }
    else b=k;
  }
  *seq_id = k;
  *rel_offset = frag->start - s_db->seq[k].start;
  return 1;
}

