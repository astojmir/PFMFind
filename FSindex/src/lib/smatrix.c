/********************************************************************/    
/*                                                                  */
/*                     SCORE_MATRIX module                          */ 
/*                                                                  */
/********************************************************************/    

#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "partition.h"

#ifdef DEBUG
#define SMATRIX_INLINE
#endif

#include "smatrix.h"

#ifdef DEBUG
#undef SMATRIX_INLINE
#endif


/* Routines for computing distances to pattern letters (i.e. subsets
   of alphabet partitions. */

static int SM_psize;          /* Current partition size */
static int SM_poffset;        /* Current partition offset */
static SCORE_MATRIX_t *SM_S;  /* Score matrix    */
static int SM_i;              /* Current Alphabet letter */
static int SM_j;              /* Current Partition */

static
void pattern_similarity(int Score, int P, int pos)
{
  int letter;
  if (pos == SM_psize)
    {
      P += SM_poffset;
      SM_S->pM[SM_i][P] = Score;
      return;
    }
  pos++;
  pattern_similarity(Score, P, pos);

  P = P | (1 << (pos-1));
  letter = FS_PARTITION_get_letter(SM_S->ptable, SM_j, pos-1); 
  letter = letter & A_SIZE_MASK;
  if (Score < SM_S->M[SM_i][letter])
    Score = SM_S->M[SM_i][letter];
  pattern_similarity(Score, P, pos);
  return;
}

static
void pattern_distance(int Score, int P, int pos)
{
  int letter;
  if (pos == SM_psize)
    {
      P += SM_poffset;
      SM_S->pM[SM_i][P] = Score;
      return;
    }
  pos++;
  pattern_distance(Score, P, pos);

  P = P | (1 << (pos-1));
  letter = FS_PARTITION_get_letter(SM_S->ptable, SM_j, pos-1); 
  letter = letter & A_SIZE_MASK;
  if (Score > SM_S->M[SM_i][letter])
    Score = SM_S->M[SM_i][letter];
  pattern_distance(Score, P, pos);
  return;
}


/* Main constructor */
SCORE_MATRIX_t *SCORE_MATRIX_create(const char *filename,
				    FS_PARTITION_t *ptable)
{
  /* Warning: as it is, this routine is only guranteed to work with
     BLOSUM62 */


  char linebuf[512]; /* line buffer - hopefully 512 bytes should be
		       enough for now */
  FILE *stream = fopen(filename, "r");
  int A[A_SIZE]; /* Alphabet table */
  int header_read = 0;
  UINT_t i, j;
  UINT_t row, col;
  int value;
  int field_width = 0;
  int j_offset;

  SCORE_MATRIX_t *S;

  if (stream == NULL)
    return NULL;

  S = callocec(1, sizeof(SCORE_MATRIX_t));
  S->similarity_flag = 1;
  S->ptable = ptable;


  /* Load similarity matrix */
  j = 0;
  while (fgets(linebuf,512,stream) != NULL)
    {
      if (linebuf[0] == '#') continue;
      if (linebuf[0] == ' ' && !header_read)
	{
	  /* Determine field width */
	  for(; linebuf[field_width] == ' '; field_width++);

	  /* Load header */
	  for (i=0; i < 23; i++)
	    A[i] = linebuf[field_width+field_width*i]; 
	  header_read = 1;
	}
      else
	{
	  /* Load data */
	  row = A[j] & A_SIZE_MASK;
	  for (i=0; i < 23; i++)
	    {
	      col = A[i] & A_SIZE_MASK;
	      sscanf(linebuf+2+field_width*i, "%d", &value); 
	      S->M[row][col] = value;
	    }
	  j++;
	  if (j >= 23) break;
	}
    }
  fclose(stream);

  /* Compute maximum similarities to pattern letters (subsets of
     alphabet partitions) */
  for (SM_i=0; SM_i < A_SIZE; SM_i++)
    {      
      if (ptable->partition_table[SM_i] == -1)
	continue;
      for (SM_j=0; SM_j < ptable->no_partitions; SM_j++)
	{
	  SM_S = S;
	  SM_psize = FS_PARTITION_get_size(S->ptable, SM_j);
	  SM_poffset =  FS_PARTITION_get_poffset(S->ptable, SM_j);

	  /* Now traverse the binary tree representing the subsets and
	     compute maximum similarities */
	  pattern_similarity(SHRT_MIN, 0, 0);
	  SM_S->pM[SM_i][SM_poffset] =
	    SM_S->pM[SM_i][SM_poffset + (1 << SM_psize) - 1];
	}
    }

  /* Compute maximum similarities to any other cluster */
  for (i=0; i < A_SIZE; i++)
    {
      if (ptable->partition_table[i] == -1)
	continue;
      S->pMclosest[i] = SHRT_MIN;
      for (j=0; j < ptable->no_partitions; j++)
	{
	  if (ptable->partition_table[i] != j)
	    {
	      j_offset = FS_PARTITION_get_poffset(ptable, j); 
	      S->pMclosest[i] = max(S->pMclosest[i], 
				    S->pM[i][j_offset]); 
	    }
	}
    }


  if (SCORE_MATRIX_VERBOSE > 0)
    SCORE_MATRIX_print(S, stdout, "*** SIMILARITY MATRIX ***");
    
  return S;
}
/* Destructor */
void SCORE_MATRIX_destroy(SCORE_MATRIX_t *Score_matrix)
{
  free(Score_matrix);
}

/* Read, Write */
int SCORE_MATRIX_write(SCORE_MATRIX_t *Score_matrix, FILE *stream)
{
  fwrite(&(Score_matrix->M[0]), sizeof(SSINT), 
	 A_SIZE * A_SIZE, stream);   
  fwrite(&(Score_matrix->pM[0]), sizeof(SSINT), 
	 A_SIZE * P_SIZE, stream);
  fwrite(&(Score_matrix->pMclosest[0]), sizeof(SSINT), 
	 A_SIZE, stream);
  fwrite(&(Score_matrix->similarity_flag), sizeof(char), 1, stream);   
  return 1;
} 

SCORE_MATRIX_t *SCORE_MATRIX_read(FILE *stream)
{
  SCORE_MATRIX_t *S = mallocec(sizeof(SCORE_MATRIX_t));
  fread(&(S->M[0]), sizeof(SSINT), 
	 A_SIZE * A_SIZE, stream);   
  fread(&(S->pM[0]), sizeof(SSINT), 
	 A_SIZE * P_SIZE, stream);
  fread(&(S->pMclosest[0]), sizeof(SSINT), 
	 A_SIZE, stream);
  fread(&(S->similarity_flag), sizeof(char), 1, stream);   

  if (SCORE_MATRIX_VERBOSE > 0)
    SCORE_MATRIX_print(S, stdout, "*** LOADED MATRIX ***");
  return S;
} 

/* Element Access + Properties */
int SCORE_MATRIX_max_entry(SCORE_MATRIX_t *S)
{
  int M = INT_MIN;
  UINT_t i, j;
  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      for (j=0; j < A_SIZE; j++)
	{
	  if (S->ptable->partition_table[j] == -1)
	    continue;
	  M = max(M, S->M[i][j]);
	}
    }
  return M;
}

int SCORE_MATRIX_min_entry(SCORE_MATRIX_t *S)
{
  int m = INT_MAX;
  UINT_t i, j;
  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      for (j=0; j < A_SIZE; j++)
	{
	  if (S->ptable->partition_table[j] == -1)
	    continue;
	  m = min(m, S->M[i][j]);
	}
    }
  return m;
}

int SCORE_MATRIX_max_entry_col(SCORE_MATRIX_t *S, int col)
{
  int M = INT_MIN;
  UINT_t i;
  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      M = max(M, S->M[i][col]);
    }
  return M;
}

int SCORE_MATRIX_max_entry_row(SCORE_MATRIX_t *S, int row)
{
  int M = INT_MIN;
  UINT_t i;
  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      M = max(M, S->M[row][i]);
    }
  return M;
}

int SCORE_MATRIX_min_entry_col(SCORE_MATRIX_t *S, int col)
{
  int m = INT_MAX;
  ULINT i;
  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      m = min(m, S->M[i][col]);
    }
  return m;
}

int SCORE_MATRIX_min_entry_row(SCORE_MATRIX_t *S, int row)
{
  int m = INT_MAX;
  UINT_t i;
  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      m = min(m, S->M[row][i]);
    }
  return m;
}

/* Similarities to Distances, Quasi-metrics ... */
SCORE_MATRIX_t *SCORE_MATRIX_S_2_Dquasi(SCORE_MATRIX_t *S)
{
  UINT_t i, j;
  SCORE_MATRIX_t *D = callocec(1, sizeof(SCORE_MATRIX_t));
  int j_offset;

  D->similarity_flag = 0;
  D->ptable = S->ptable;

  for (i=0; i < A_SIZE; i++)
    {
      if (D->ptable->partition_table[i] == -1)
	continue;
      for (j=0; j < A_SIZE; j++)
	{
	  if (D->ptable->partition_table[j] == -1)
	    continue;
	  D->M[i][j] = S->M[i][i] - S->M[i][j];
	}
    }

  /* TO DO: */
  /* Here use some of the graph shortest path algorithms */

  /* Compute maximum distances to pattern letters (subsets of
     alphabet partitions) */
  for (SM_i=0; SM_i < A_SIZE; SM_i++)
    {      
      if (D->ptable->partition_table[SM_i] == -1)
	continue;
      for (SM_j=0; SM_j < D->ptable->no_partitions; SM_j++)
	{
	  SM_S = D;
	  SM_psize = FS_PARTITION_get_size(D->ptable, SM_j);
	  SM_poffset =  FS_PARTITION_get_poffset(D->ptable, SM_j);

	  /* Now traverse the binary tree representing the subsets and
	     compute maximum similarities */
	  pattern_distance(SHRT_MAX, 0, 0);
	  SM_S->pM[SM_i][SM_poffset] =
	    SM_S->pM[SM_i][SM_poffset + (1 << SM_psize) - 1];
	}
    }

  /* Compute minimum distances to any cluster */
  for (i=0; i < A_SIZE; i++)
    {
      if (D->ptable->partition_table[i] == -1)
	continue;
      D->pMclosest[i] = SHRT_MAX;
      for (j=0; j < D->ptable->no_partitions; j++)
	{
	  if (D->ptable->partition_table[i] != j)
	    {
	      j_offset = 
		FS_PARTITION_get_poffset(D->ptable, j);
	      D->pMclosest[i] = min(D->pMclosest[i], 
				    D->pM[i][j_offset]);
	    }
	}
    }

  if (SCORE_MATRIX_VERBOSE > 0)
    SCORE_MATRIX_print(D, stdout, "*** QUASI-METRIC MATRIX ***");
  return D;
}

SCORE_MATRIX_t *SCORE_MATRIX_S_2_Dmax(SCORE_MATRIX_t *S)
{
  SCORE_MATRIX_t *D = SCORE_MATRIX_S_2_Dquasi(S);
  UINT_t i, j;
  int j_offset;

  /* Symmetrise */
  for (i=0; i < A_SIZE-1; i++)
    {
      if (D->ptable->partition_table[i] == -1)
	continue;
      for (j=i+1; j < A_SIZE; j++)
	{
	  if (D->ptable->partition_table[j] == -1)
	    continue;
	  D->M[i][j] = max(D->M[i][j], D->M[j][i]);
	  D->M[j][i] = D->M[i][j];
	}
    }


  /* Compute maximum distances to pattern letters (subsets of
     alphabet partitions) */
  for (SM_i=0; SM_i < A_SIZE; SM_i++)
    {      
      if (D->ptable->partition_table[SM_i] == -1)
	continue;
      for (SM_j=0; SM_j < D->ptable->no_partitions; SM_j++)
	{
	  SM_S = D;
	  SM_psize = FS_PARTITION_get_size(D->ptable, SM_j);
	  SM_poffset =  FS_PARTITION_get_poffset(D->ptable, SM_j);

	  /* Now traverse the binary tree representing the subsets and
	     compute maximum similarities */
	  pattern_distance(SHRT_MAX, 0, 0);
	  SM_S->pM[SM_i][SM_poffset] =
	    SM_S->pM[SM_i][SM_poffset + (1 << SM_psize) - 1];
	}
    }

  /* Compute minimum distances to any cluster */
  for (i=0; i < A_SIZE; i++)
    {
      if (D->ptable->partition_table[i] == -1)
	continue;
      D->pMclosest[i] = SHRT_MAX;
      for (j=0; j < D->ptable->no_partitions; j++)
	{
	  if (D->ptable->partition_table[i] != j)
	    {
	      j_offset = 
		FS_PARTITION_get_poffset(D->ptable, j);
	      D->pMclosest[i] = min(D->pMclosest[i], 
				    D->pM[i][j_offset]);
	    }
	}
    }

  if (SCORE_MATRIX_VERBOSE > 0)
    SCORE_MATRIX_print(D, stdout, "*** Dmax MATRIX ***");
  return D;
}

SCORE_MATRIX_t *SCORE_MATRIX_S_2_Davg(SCORE_MATRIX_t *S)
{
  SCORE_MATRIX_t *D = SCORE_MATRIX_S_2_Dquasi(S);
  UINT_t i, j;
  int j_offset;

  /* Symmetrise */
  for (i=0; i < A_SIZE-1; i++)
    {
      if (D->ptable->partition_table[i] == -1)
	continue;
      for (j=i+1; j < A_SIZE; j++)
	{
	  if (D->ptable->partition_table[j] == -1)
	    continue;
	  D->M[i][j] = D->M[i][j] + D->M[j][i];
	  D->M[j][i] = D->M[i][j];
	}
    }

  /* Compute maximum distances to pattern letters (subsets of
     alphabet partitions) */
  for (SM_i=0; SM_i < A_SIZE; SM_i++)
    {      
      if (D->ptable->partition_table[SM_i] == -1)
	continue;
      for (SM_j=0; SM_j < D->ptable->no_partitions; SM_j++)
	{
	  SM_S = D;
	  SM_psize = FS_PARTITION_get_size(D->ptable, SM_j);
	  SM_poffset =  FS_PARTITION_get_poffset(D->ptable, SM_j);

	  /* Now traverse the binary tree representing the subsets and
	     compute maximum similarities */
	  pattern_distance(SHRT_MAX, 0, 0);
	  SM_S->pM[SM_i][SM_poffset] =
	    SM_S->pM[SM_i][SM_poffset + (1 << SM_psize) - 1];
	}
    }

  /* Compute minimum distances to any cluster */
  for (i=0; i < A_SIZE; i++)
    {
      if (D->ptable->partition_table[i] == -1)
	continue;
      D->pMclosest[i] = SHRT_MAX;
      for (j=0; j < D->ptable->no_partitions; j++)
	{
	  if (D->ptable->partition_table[i] != j)
	    {
	      j_offset = 
		FS_PARTITION_get_poffset(D->ptable, j);
	      D->pMclosest[i] = min(D->pMclosest[i], 
				    D->pM[i][j_offset]);
	    }
	}
    }

  if (SCORE_MATRIX_VERBOSE > 0)
    SCORE_MATRIX_print(D, stdout, "*** Davg MATRIX ***");
  return D;
}

/* Printing */
void SCORE_MATRIX_print(SCORE_MATRIX_t *S, FILE *stream, 
			const char *title)
{
  UINT_t i, j;
  int p=0;
  int lptl;
  int lp;
  char s[C_SIZE+1];

  fprintf(stream, "%s\n\n", title);

  /* Print M */
  fprintf(stream, "%s\n", "Scoring matrix");
  fprintf(stream, "    ");
  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      fprintf(stream, " %c  ", i+64);
    }
  fprintf(stream, "\n");

  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      fprintf(stream, " %c ", i+64);
      for (j=0; j < A_SIZE; j++)
	{
	  if (S->ptable->partition_table[j] == -1)
	    continue;
	  fprintf(stream, "%3d ", S->M[i][j]);
	}
      fprintf(stream, "\n");
    }

  /* Print pM - transposed for easier viewing */

  fprintf(stream, "\n\n");
  fprintf(stream, "%s\n", "Scores w.r.t. partitions");

  lp = FS_PARTITION_get_no_partitions(S->ptable) -1;
  lptl =  FS_PARTITION_get_poffset(S->ptable, lp)
    + (1 << FS_PARTITION_get_size(S->ptable, lp)); 
  for (j=0; j < lptl; j++)
    {
      if (j == FS_PARTITION_get_poffset(S->ptable, p))
	{
	  p++;
	  fprintf(stream, "\n");
	  fprintf(stream, "%*s  ", C_SIZE, "");
	  for (i=0; i < A_SIZE; i++)
	    {
	      if (S->ptable->partition_table[i] == -1)
		continue;
	      fprintf(stream, " %c  ", i+64);
	    }
	  fprintf(stream, "\n");
	  fprintf(stream, "\n");
#ifdef DEBUG
	  fprintf(stream, "CLSTR #%1d ", p);
	  for (i=0; i < A_SIZE; i++)
	    {
	      if (S->ptable->partition_table[i] == -1)
		continue;
	      fprintf(stream, "%3d ", S->pM[i][j]);
	    }
	  fprintf(stream, "\n");
#endif /* #ifdef DEBUG */
	  continue;
	}
      FS_PARTITION_pletter_2_string(S->ptable, j, &s[0]);
      fprintf(stream, "%-*s ", C_SIZE, s);
      for (i=0; i < A_SIZE; i++)
	{
	  if (S->ptable->partition_table[i] == -1)
	    continue;
	  fprintf(stream, "%3d ", S->pM[i][j]);
	}
      fprintf(stream, "\n");
    }

  /* Print pMclosest - transposed for easier viewing */
  fprintf(stream, "\n\n");
  fprintf(stream, "%s\n", "Best scores w.r.t. other partitions");
  fprintf(stream, "   ");

  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      fprintf(stream, " %2c ", i+64);
    }
  fprintf(stream, "\n");

  fprintf(stream, "%2s ", "");
  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      fprintf(stream, "%3d ", S->pMclosest[i]);
    }
  fprintf(stream, "\n");
  fprintf(stream, "\n");
}
