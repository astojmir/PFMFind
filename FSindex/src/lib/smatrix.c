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
  UINT_t i, j, k;
  UINT_t row, col;
  int value;
  int field_width = 0;

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

  /* Compute maximum similarities to clusters */
  for (i=0; i < A_SIZE; i++)
    {
      if (ptable->partition_table[i] == -1)
	continue;
      for (j=0; j < ptable->no_partitions; j++)
	{
	  S->pM[i][j] = SHRT_MIN;
	  for (k=0; k < A_SIZE; k++)
	    if (ptable->partition_table[k] == j)
	      {
		S->pM[i][j] = max(S->pM[i][j], S->M[i][k]);
	      }
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
	    S->pMclosest[i] = max(S->pMclosest[i], S->pM[i][j]);
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
  UINT_t i, j, k;
  SCORE_MATRIX_t *D = callocec(1, sizeof(SCORE_MATRIX_t));

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

  /* Compute minimum distances to clusters */
  for (i=0; i < A_SIZE; i++)
    {
      if (D->ptable->partition_table[i] == -1)
	continue;
      for (j=0; j < D->ptable->no_partitions; j++)
	{
	  D->pM[i][j] = SHRT_MAX;
	  for (k=0; k < A_SIZE; k++)
	    if (D->ptable->partition_table[k] == j)
	      D->pM[i][j] = min(D->pM[i][j], D->M[i][k]);
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
	    D->pMclosest[i] = min(D->pMclosest[i], D->pM[i][j]);
	}
    }

  if (SCORE_MATRIX_VERBOSE > 0)
    SCORE_MATRIX_print(D, stdout, "*** QUASI-METRIC MATRIX ***");
  return D;
}

SCORE_MATRIX_t *SCORE_MATRIX_S_2_Dmax(SCORE_MATRIX_t *S)
{
  SCORE_MATRIX_t *D = SCORE_MATRIX_S_2_Dquasi(S);
  UINT_t i, j, k;

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

  /* Compute minimum distances to clusters */
  for (i=0; i < A_SIZE; i++)
    {
      if (D->ptable->partition_table[i] == -1)
	continue;
      for (j=0; j < D->ptable->no_partitions; j++)
	{
	  D->pM[i][j] = SHRT_MAX;
	  for (k=0; k < A_SIZE; k++)
	    if (D->ptable->partition_table[k] == j)
	      D->pM[i][j] = min(D->pM[i][j], D->M[i][k]);
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
	    D->pMclosest[i] = min(D->pMclosest[i], D->pM[i][j]);
	}
    }

  if (SCORE_MATRIX_VERBOSE > 0)
    SCORE_MATRIX_print(D, stdout, "*** Dmax MATRIX ***");
  return D;
}

SCORE_MATRIX_t *SCORE_MATRIX_S_2_Davg(SCORE_MATRIX_t *S)
{
  SCORE_MATRIX_t *D = SCORE_MATRIX_S_2_Dquasi(S);
  UINT_t i, j, k;

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

  /* Compute minimum distances to clusters */
  for (i=0; i < A_SIZE; i++)
    {
      if (D->ptable->partition_table[i] == -1)
	continue;
      for (j=0; j < D->ptable->no_partitions; j++)
	{
	  D->pM[i][j] = SHRT_MAX;
	  for (k=i; k < A_SIZE; k++)
	    if (D->ptable->partition_table[k] == j)
	      D->pM[i][j] = min(D->pM[i][j], D->M[i][k]);
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
	    D->pMclosest[i] = min(D->pMclosest[i], D->pM[i][j]);
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
  fprintf(stream, "   ");

  for (i=0; i < A_SIZE; i++)
    {
      if (S->ptable->partition_table[i] == -1)
	continue;
      fprintf(stream, " %2c ", i+64);
    }
  fprintf(stream, "\n");


  for (j=0; j < S->ptable->no_partitions; j++)
    {
      fprintf(stream, "%2d ", j);
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
