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


/********************************************************************/    
/*                                                                  */
/*                     SCORE_MATRIX module                          */ 
/*                                                                  */
/********************************************************************/    

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "partition.h"
#include "smatrix.h"
#include "misclib.h"

/* Evaluation functions */
static
int score_eval(SCORE_MATRIX *S, const char *s1, const char *s2, int len1)
{
  int **MM = S->M;
  int i;
  int H = 0;

  for(i=0; i < len1; s1++, s2++, i++)    
    H += MM[*s1 & A_SIZE_MASK][*s2 & A_SIZE_MASK]; 
  return H;
}

static
int profile_eval(SCORE_MATRIX *S, const char *s1, const char *s2, int len1)
{
  int **MM = S->M;
  int i;
  int H = 0;
  int l = min(S->len, len1);

  for(i=0; i < l; s1++, MM++, i++)    
    H += (*MM)[*s1 & A_SIZE_MASK];
  return H;
}

/* Item conversion functions */
static
int PD_item_conv(SCORE_MATRIX *S, char c, int i, int j)
{
  char cj = 64+j;
  if (!bsearch(&cj, S->alphabet, S->alen, 1, compare_char))
    return VERY_LARGE_DISTANCE;
  else
    return S->M[i][j];
}

static 
int SD_item_conv(SCORE_MATRIX *S, char c, int i, int j)
{
  char cj = 64+j;
  if (!bsearch(&cj, S->alphabet, S->alen, 1, compare_char))
    return VERY_LARGE_DISTANCE;
  else
    return S->M[c & A_SIZE_MASK][j];
}

static 
int PS_item_conv(SCORE_MATRIX *S, char c, int i, int j)
{
  char cj = 64+j;
  if (!bsearch(&cj, S->alphabet, S->alen, 1, compare_char))
    return VERY_LARGE_DISTANCE;
  else 
    return S->M[i][S->qseq[i] & A_SIZE_MASK] - S->M[i][j];
}

static 
int SS_item_conv_quasi(SCORE_MATRIX *S, char c, int i, int j)
{
  int k = c & A_SIZE_MASK;
  char cj = 64+j;
  if (!bsearch(&cj, S->alphabet, S->alen, 1, compare_char) ||
      !bsearch(&c, S->alphabet, S->alen, 1, compare_char))
    return VERY_LARGE_DISTANCE;
  else 
    return S->M[k][k] - S->M[k][j];
}

static 
int SS_item_conv_max(SCORE_MATRIX *S, char c, int i, int j)
{
  int k = c & A_SIZE_MASK;
  char cj = 64+j;
  if (!bsearch(&cj, S->alphabet, S->alen, 1, compare_char) ||
      !bsearch(&c, S->alphabet, S->alen, 1, compare_char))
    return VERY_LARGE_DISTANCE;
  else 
    return max(S->M[k][k] - S->M[k][j], S->M[j][j] - S->M[j][k]);
}

static 
int SS_item_conv_avg(SCORE_MATRIX *S, char c, int i, int j)
{
  int k = c & A_SIZE_MASK;
  char cj = 64+j;
  if (!bsearch(&cj, S->alphabet, S->alen, 1, compare_char) ||
      !bsearch(&c, S->alphabet, S->alen, 1, compare_char))
    return VERY_LARGE_DISTANCE;
  else 
    return S->M[k][k] - S->M[k][j] + S->M[j][j] - S->M[j][k];
}

/* Range conversion functions */
static
int id_range_conv(SCORE_MATRIX *S, const char *s, int len1, int r)
{
  return r;
}

static
int SS_range_conv(SCORE_MATRIX *S, const char *s, int len1, int r)
{
  return score_eval(S, s, s, len1) - r;
}

static
int PS_range_conv(SCORE_MATRIX *S, const char *s, int len1, int r)
{
  return profile_eval(S, S->qseq, NULL, S->len) - r;
}
 
/* Matrix conversion functions */
static
SCORE_MATRIX *S_matrix_conv(SCORE_MATRIX *S, const char *s, int len1)
{
  int i, j;
  int **NM;
  if (S->conv_type & POS) {
    NM = mallocec(len1 * sizeof(int *));
    for (i=0; i < len1; i++) {
      NM[i] = mallocec(A_SIZE * sizeof(int));
      for (j=0; j < A_SIZE; j++) 
	NM[i][j] = S->item_conv(S, s[i], i, j);
    }
    return SCORE_MATRIX_init(NM, len1, POSITIONAL, DISTANCE,
			     strdup(S->alphabet));
  }
  else {
    NM = mallocec(A_SIZE * sizeof(int *));
    for (i=0; i < A_SIZE; i++) {
      NM[i] = mallocec(A_SIZE * sizeof(int));
      for (j=0; j < A_SIZE; j++) 
	NM[i][j] = S->item_conv(S, i+64, i, j);
    }
    return SCORE_MATRIX_init(NM, len1, SCORE, DISTANCE,
			     strdup(S->alphabet));
  }
}

static
SCORE_MATRIX *P_matrix_conv(SCORE_MATRIX *S, const char *s, int len1)
{
  int i, j;
  int **NM;
  NM = mallocec(S->len * sizeof(int *));
  for (i=0; i < S->len; i++) {
    NM[i] = mallocec(A_SIZE * sizeof(int));
    for (j=0; j < A_SIZE; j++) 
      NM[i][j] = S->item_conv(S, 0, i, j);
  }
  return SCORE_MATRIX_init(NM, S->len, POSITIONAL, DISTANCE,
			   strdup(S->alphabet));
}

/* Set conversion type functions */
static
void id_set_conv_type(SCORE_MATRIX *S, int conv_type)
{
  return;
}

static
void SS_set_conv_type(SCORE_MATRIX *S, int conv_type)
{
  switch (conv_type & 6) {
  case QUASI:
    S->item_conv = SS_item_conv_quasi;
    S->range_conv = SS_range_conv;
    break;
  case MAX:
    S->item_conv = SS_item_conv_max;
    S->range_conv = id_range_conv;
    break;
  case AVG:
    S->item_conv = SS_item_conv_avg;
    S->range_conv = id_range_conv;
    break;
  default:
    conv_type = (S->conv_type & 6) | (conv_type & 1);
  }
  S->conv_type = conv_type & 7;
}

/* Constructors and destructor */
SCORE_MATRIX *SCORE_MATRIX_init(int **M, int len, MATRIX_TYPE Mtype,
				SCORE_TYPE Stype, char *a)
{
  /* M and alphabet are NOT hard copied - their pointers are just
     assigned */

  int i,j;
  int k;
  int l = strlen(a);
  char *c;
  SCORE_MATRIX *S;
  
  /* Check that each letter in alphabet appears only once */
  qsort(a, l, 1, compare_char);
  for (j=1; j < l; j++) 
    if (a[j] == a[j-1]) 
      return NULL;
  
  S = callocec(1, sizeof(SCORE_MATRIX));
  S->M = M;
  S->Mtype = Mtype;
  S->Stype = Stype;
  S->alphabet = a;
  S->alen = l;
  if (Mtype == SCORE) {
    S->matrix_conv = S_matrix_conv;
    S->eval_score = score_eval;
    if (Stype == DISTANCE) {
      S->conv_type = POS;
      S->item_conv = SD_item_conv;
      S->range_conv = id_range_conv;
      S->set_conv_type = id_set_conv_type;
    }
    else { /* Stype == SIMILARITY */
      S->conv_type = QUASI;
      S->item_conv = SS_item_conv_quasi;
      S->range_conv = SS_range_conv;
      S->set_conv_type = SS_set_conv_type;
    }
  }
  else { /* Mtype == POSITIONAL */
    S->len = len;
    S->eval_score = profile_eval;
    S->qseq = mallocec(len+1);
    S->qseq[len] = '\0';
    if (Stype == DISTANCE) {
      S->matrix_conv = NULL;
      S->conv_type = 0;
      S->item_conv = PD_item_conv;
      S->range_conv = id_range_conv;
      S->set_conv_type = id_set_conv_type;
      /* Calculate qseq - should have letters with minimal dist (0)*/
      for (i=0; i < len; i++) {
	for (c=a, k=INT_MAX; *c; c++) {
	  if (M[i][*c & A_SIZE_MASK] < k) {
	    k = M[i][*c & A_SIZE_MASK];
	    S->qseq[i] = *c;
	  }
	}
      }
    }
    else { /* Stype == SIMILARITY */
      S->matrix_conv = P_matrix_conv;
      S->conv_type = QUASI;
      S->item_conv = PS_item_conv;
      S->range_conv = PS_range_conv;
      S->set_conv_type = id_set_conv_type;
      /* Calculate qseq */
      for (i=0; i < len; i++) {
	for (c=a, k=INT_MIN; *c; c++) {
	  if (M[i][*c & A_SIZE_MASK] > k) {
	    k = M[i][*c & A_SIZE_MASK];
	    S->qseq[i] = *c;
	  }
	}
      }
    }
  }
  return S;
}

SCORE_MATRIX *SCORE_MATRIX_copy(SCORE_MATRIX *S)
{
  int i;
  SCORE_MATRIX *T = callocec(1, sizeof(SCORE_MATRIX));
  int ldim;

  T->Mtype = S->Mtype;
  T->Stype = S->Stype;
  T->len = S->len;
  T->alen = S->alen;
  T->conv_type = S->conv_type;
  T->eval_score = S->eval_score;
  T->item_conv = S->item_conv;
  T->range_conv = S->range_conv;
  T->set_conv_type = S->set_conv_type;
  T->matrix_conv = S->matrix_conv;

  if (S->len == 0) {
    ldim = A_SIZE;
    T->qseq = NULL;
  }
  else {
    ldim = S->len;
    T->qseq = mallocec(ldim);
    memcpy(T->qseq, S->qseq, ldim);
  }

  T->M = mallocec(ldim * sizeof(int *));
  for (i=0; i < ldim; i++) {
    T->M[i] = mallocec(A_SIZE * sizeof(int));
    memcpy(T->M[i], S->M[i], A_SIZE * sizeof(int));
  }

  T->alphabet = strdup(S->alphabet);

  return T;
}

SCORE_MATRIX *SCORE_MATRIX_from_file(const char *filename)
{
  /* Warning: as it is, this routine is only guaranteed to work with
     BLOSUM62 - it can possibly be generalised */


  char linebuf[512]; /* line buffer - hopefully 512 bytes should be
		       enough for now */
  FILE *stream = fopen(filename, "r");
  int header_read = 0;
  int i, j;
  int row, col;
  int value;
  int field_width = 0;
  int **M;
  char *a;

  if (stream == NULL)
    Throw FSexcept(FOPEN_ERR,
		   "SCORE_MATRIX_from_file():"
		   "Could not open file %s", filename);

  M = mallocec(A_SIZE * sizeof(int *));
  for (i=0; i < A_SIZE; i++)
    M[i] = callocec(A_SIZE, sizeof(int));
  a = mallocec(A_SIZE + 1);

  /* Load similarity matrix */
  j = 0;
  while (fgets(linebuf,512,stream) != NULL) {
    if (linebuf[0] == '#') continue;
    if (linebuf[0] == ' ' && !header_read) {
      /* Determine field width */
      for(; linebuf[field_width] == ' '; field_width++);
      /* Load header */
      for (i=0; i < 23; i++)
	a[i] = linebuf[field_width+field_width*i]; 
      header_read = 1;
    }
    else {
      /* Load data */
      row = a[j] & A_SIZE_MASK;
      for (i=0; i < 23; i++) {
	col = a[i] & A_SIZE_MASK;
	sscanf(linebuf+2+field_width*i, "%d", &value); 
	M[row][col] = value;
      }
      j++;
      if (j >= 23) break;
    }
  }
  fclose(stream);
  a[23] = '\0';
  return SCORE_MATRIX_init(M, 0, SCORE, SIMILARITY, a);
}

void SCORE_MATRIX_del(SCORE_MATRIX *S)
{
  int i;
  int m = S->Mtype == SCORE ? A_SIZE : S->len;
  for (i = 0; i < m; i++)
    free(S->M[i]);
  free(S->M);
  if (S->qseq != NULL)
    free(S->qseq);
  free(S->alphabet);
  free(S);
}


/* Printing */
char *SCORE_MATRIX_sprint(SCORE_MATRIX *S)
{
  char *buf = NULL;
  int size = 0;
  int len = 0;
  int i;
  char *a;
  char *b;
  
  if (S->Stype == DISTANCE)
    absprintf(&buf, &size, &len, "DISTANCE MATRIX\n");
  else
    absprintf(&buf, &size, &len, "SIMILARITY MATRIX\n");

  if (S->Mtype == POSITIONAL) {
    absprintf(&buf, &size, &len, "%2.2s ", "");
    for (i=0; i < S->len; i++) 
      absprintf(&buf, &size, &len, "%4d ", i);
    absprintf(&buf, &size, &len, "\n");
    absprintf(&buf, &size, &len, "%2.2s ", "");
    for (i=0; i < S->len; i++) 
      absprintf(&buf, &size, &len, "   %c ", S->qseq[i]);
    absprintf(&buf, &size, &len, "\n");
    for (a=S->alphabet; *a; a++) {
      absprintf(&buf, &size, &len, " %c ", *a);
      for (i=0; i < S->len; i++) 
	absprintf(&buf, &size, &len, "%4d ", S->M[i][*a & A_SIZE_MASK]);
      absprintf(&buf, &size, &len, "\n");
    }
  }
  else {
    absprintf(&buf, &size, &len, "%3.3s ", "");
    for (a=S->alphabet; *a; a++) 
      absprintf(&buf, &size, &len, "  %c ", *a);
    absprintf(&buf, &size, &len, "\n");
    for (b=S->alphabet; *b; b++) {
      absprintf(&buf, &size, &len, "  %c ", *b);
      for (a=S->alphabet; *a; a++)
	absprintf(&buf, &size, &len, "%3d ", S->M[*b & A_SIZE_MASK][*a & A_SIZE_MASK]);
      absprintf(&buf, &size, &len, "\n");
    }
  }
  absprintf(&buf, &size, &len, "\n");
  return buf;
}

void SCORE_MATRIX_fprint(SCORE_MATRIX *S, FILE *fp)
{
  char *buf = SCORE_MATRIX_sprint(S);
  fprintf(fp, "%s", buf);
  free(buf);
}


void SCORE_MATRIX_write(SCORE_MATRIX *S, FILE *fp)
{
  int i, m;

  fwrite(&S->Mtype, sizeof(MATRIX_TYPE), 1, fp);
  fwrite(&S->Stype, sizeof(SCORE_TYPE), 1, fp);
  fwrite(&S->len, sizeof(int), 1, fp);
  fwrite(&S->alen, sizeof(int), 1, fp);
  fwrite(S->alphabet, sizeof(char), S->alen, fp);
  
  m = S->Mtype == SCORE ? A_SIZE : S->len;
  for (i=0; i < m; i++)
    fwrite(S->M[i], sizeof(int), A_SIZE, fp);
}

SCORE_MATRIX *SCORE_MATRIX_read(FILE *fp)
{
  int i, m;
  MATRIX_TYPE Mtype;
  SCORE_TYPE Stype;
  int len;
  int alen;
  int **M;
  char *alphabet;


  fread(&Mtype, sizeof(MATRIX_TYPE), 1, fp);
  fread(&Stype, sizeof(SCORE_TYPE), 1, fp);
  fread(&len, sizeof(int), 1, fp);
  fread(&alen, sizeof(int), 1, fp);
  alphabet = callocec(alen+1, sizeof(char)); 
  fread(alphabet, sizeof(char), alen, fp);
  
  m = Mtype == SCORE ? A_SIZE : len;
  M = mallocec(m * sizeof(int *));
  for (i=0; i < m; i++) {
    M[i] = mallocec(A_SIZE * sizeof(int));
    fread(M[i], sizeof(int), A_SIZE, fp);
  }

  return SCORE_MATRIX_init(M, len, Mtype, Stype, alphabet);
}


unsigned char *SCORE_MATRIX_swrite(SCORE_MATRIX *S, int *size)
{
  /* Use byte-order-independent functions to write the matrix as a
     byte stream. All ints are 32 bits. */ 

  int i, msize, x, k;
  unsigned char *s;
  char *c, *d;

  msize = S->Mtype == SCORE ? S->alen*S->alen : S->len*S->alen;
 
  *size = 4 + 4 + 4 + msize*4 + S->alen+1;
  s = mallocec(*size);

  k = 0;
  x = (int) S->Mtype;
  swrite_int32(&x, 1, s+k);
  k += 4;
  x = (int) S->Stype;
  swrite_int32(&x, 1, s+k);
  k += 4;
  swrite_int32(&S->len, 1, s+k);
  k += 4;
  strncpy(s+k, S->alphabet, S->alen+1);
  k += S->alen+1;


  if (S->Mtype == SCORE) {
    for (d=S->alphabet; *d; d++)
      for (c=S->alphabet; *c; c++) {
	swrite_int32(S->M[*d & A_SIZE_MASK] + (*c & A_SIZE_MASK), 1, s+k); 
	k += 4;
      }
  }
  else {
    for (i=0; i < S->len; i++)
      for (c=S->alphabet; *c; c++) {
	swrite_int32(S->M[i] + (*c & A_SIZE_MASK), 1, s+k);
	k += 4;
      }
  }

  printf("%d %d\n", *size, k);
  return s;
}

SCORE_MATRIX *SCORE_MATRIX_sread(unsigned char *s)
{
  /* Use byte-order-independent functions to write the matrix as a
     byte stream. All ints are 32 bits. */ 

  int i, x, k;
  char *c, *d;

  MATRIX_TYPE Mtype;
  SCORE_TYPE Stype;
  int len;
  int **NM;
  char *alphabet;


  k = 0;

  sread_int32(&x, 1, s+k);
  Mtype = x;
  k += 4;
  sread_int32(&x, 1, s+k);
  Stype = x;
  k += 4;
  sread_int32(&len, 1, s+k);
  k += 4;
  alphabet = strdup(s+k);
  k += strlen(alphabet)+1;

  if (Mtype == SCORE) {
    NM = mallocec(A_SIZE * sizeof(int *));
    for (i=0; i < A_SIZE; i++)
      NM[i] = mallocec(A_SIZE * sizeof(int));
    for (d=alphabet; *d; d++)
      for (c=alphabet; *c; c++) {
	sread_int32(NM[*d & A_SIZE_MASK] + (*c & A_SIZE_MASK), 1, s+k); 
	k += 4;
      }
  }
  else {
    NM = mallocec(len * sizeof(int *));
    for (i=0; i < len; i++)
      NM[i] = mallocec(A_SIZE * sizeof(int));
      for (c=alphabet; *c; c++) {
	sread_int32(NM[i] + (*c & A_SIZE_MASK), 1, s+k);
	k += 4;
      }
  }

  return SCORE_MATRIX_init(NM, len, Mtype, Stype, alphabet);
}
