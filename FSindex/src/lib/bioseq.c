/*
 * Copyright (C) 2005-2006 Victoria University of Wellington
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
#include <time.h>
#include "bioseq.h"
#include "misclib.h"

int cmp_bioseq(const void *S1, const void *S2)
{
  BIOSEQ *T1 = (BIOSEQ *)S1;
  BIOSEQ *T2 = (BIOSEQ *)S2;

  
  if (T1->len > T2->len)
    return 1;
  else if (T1->len < T2->len)
    return -1;
  else
    return strncmp(T1->start, T2->start, T1->len);
}

/* int bioseq_parse(BIOSEQ *seq, char *filename, yn_bool_t defline)

   Writes the sequence seq into file filename (or stdout if filename
   is NULL. 
   Returns 1 if successful, otherwise 0. 
                      

   Parameters: 

   BIOSEQ *seq - Sequence to print.

   char *filename - The name of file to write the sequence to. If NULL
                    the sequence is written to stdout.

   yn_bool_t defline - Whether or not the comment line is stored.

*/
int bioseq_parse(BIOSEQ *seq, char *filename, yn_bool_t defline)
{
  ULINT i;
  const USINT char_per_row = 80; 
  FILE *out_stream;

  if (filename == NULL)
    out_stream = stdout;
  else if ((out_stream = fopen(filename, "w")) == NULL)
    return 0;

  if (defline)
    fprintf(out_stream, ">%s\n", seq->id.defline);
  else
    fprintf(out_stream, ">%ld\n", seq->id.id_num);

  for (i=0; i < seq->len; i++)
    {
      fputc(seq->start[i], out_stream);
      if ((i+1)%char_per_row == 0)
	fputc('\n', out_stream);
    }
  if (i%char_per_row != 0)
    fputc('\n', out_stream);

  if (filename != NULL)
    fclose(out_stream);
  
  return 1;
} 

/* We assume frag is allocated ! */
/* We always put id_no in id field - no text */
int bioseq_get_frag(BIOSEQ *seq, BIOSEQ *frag, ULINT start, ULINT
		    length, ULINT id_no)
{
  if (start + length > seq->len)
    return BAD_ARGS;

  frag->id.id_num = id_no;
  frag->len = length;
  frag->start = seq->start+start;
  return NO_ERR;
}

void string2seq(BIOSEQ *seq, char *string)
{
  /* seq must be allocated */
  seq->len = strlen(string);
  seq->start = string;
  seq->id.defline = NULL;
}

void bioseq_random(BIOSEQ *seq, ULINT seq_len, char *alphabet,
		   ULINT a_len) 
{
  ULINT i;
  long seed = time(NULL);
  srand48(seed);
  seq->len = seq_len;
  seq->id.defline = NULL;
  seq->start = mallocec(seq_len);
  for (i=0; i < seq_len; i++)
    {
      seq->start[i] = alphabet[lrand48()%a_len];
    }
}

void bioseq_seq2string(BIOSEQ *seq, char *string, ULINT from, 
		       ULINT to)
{
  /* string must be allocated to to-from+2 bytes! */
  if (from > to || from >= seq->len) return;
  memcpy(string, seq->start+from, to-from+1);
  string[to-from+1] = '\0';
}

BIOSEQ *bioseq_copy(BIOSEQ *seq)
{
  BIOSEQ *newseq = mallocec(sizeof(BIOSEQ));
  
  newseq->id = seq->id;
  newseq->len = seq->len;
  newseq->start = mallocec(newseq->len+1);
  memcpy(newseq->start, seq->start, newseq->len);
  newseq->start[newseq->len]='\0';
  return newseq;
}
