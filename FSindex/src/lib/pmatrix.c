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
#include "smatrix.h"

#ifdef DEBUG
#define POS_MATRIX_INLINE
#endif

#include "pmatrix.h"

#ifdef DEBUG
#undef POS_MATRIX_INLINE
#endif


POS_MATRIX *SCORE_2_POS_MATRIX(SCORE_MATRIX_t *S, BIOSEQ *query)
{
  /* As all this routine does is copying, it can be used for both
     similarity and distance matrices */


  POS_MATRIX *PS = callocec(1, sizeof(POS_MATRIX)); 
  int i;
  int j;
  int k;
  char *q = query->start;

  PS->len = query->len;
  PS->ptable = S->ptable;
  PS->similarity_flag = S->similarity_flag;

  PS->M = mallocec(PS->len * A_SIZE * sizeof(SSINT));
  PS->pM = mallocec(PS->len * P_SIZE * sizeof(SSINT));
  PS->pMclosest = mallocec(PS->len * sizeof(SSINT));



  /* Copy matrices position-wise */

  for (i=0; i < PS->len; i++, q++)
    {
      k = *q & A_SIZE_MASK;
      for (j=0; j < A_SIZE; j++)
	{
	  PS->M[PM_M(i,j)] = 
	    S->M[k][j];
	}
      for (j=0; j < P_SIZE; j++)
	{
	  PS->pM[PM_pM(i,j)] = 
	    S->pM[k][j];
	}
      PS->pMclosest[i] = S->pMclosest[k];
    }

  return PS;
}

/* Destructor */
void POS_MATRIX_destroy(POS_MATRIX *PS)
{
  free(PS->M);
  free(PS->pM);
  free(PS->pMclosest);
  free(PS);
}

/* Read, Write */
int POS_MATRIX_write(POS_MATRIX *PS, FILE *stream) 
{
  fwrite(&(PS->len), sizeof(ULINT), 1, stream);   
  fwrite(PS->M, sizeof(SSINT), PS->len * A_SIZE, stream);   
  fwrite(PS->pM, sizeof(SSINT), PS->len * P_SIZE, stream);   
  fwrite(PS->pMclosest, sizeof(SSINT), PS->len, stream);   
  fwrite(&(PS->similarity_flag), sizeof(char), 1, stream);   
  return 1;
} 

POS_MATRIX *POS_MATRIX_read(FILE *stream)
{
  POS_MATRIX *PS = mallocec(sizeof(POS_MATRIX));  
  fread(&(PS->len), sizeof(ULINT), 1, stream);   
  PS->M = mallocec(PS->len * A_SIZE * sizeof(SSINT));
  PS->pM = mallocec(PS->len * P_SIZE * sizeof(SSINT));
  PS->pMclosest = mallocec(PS->len * sizeof(SSINT));
  fread(PS->M, sizeof(SSINT), PS->len * A_SIZE, stream);   
  fread(PS->pM, sizeof(SSINT), PS->len * P_SIZE, stream);   
  fread(PS->pMclosest, sizeof(SSINT), PS->len, stream);   
  fread(&(PS->similarity_flag), sizeof(char), 1, stream);   
  return PS;
} 

/* Printing */
void POS_MATRIX_print(POS_MATRIX *PS, FILE *stream, 
			const char *title)
{
  int i, j;
  
  fprintf(stream, "%s\n\n", title);

  /* Print M - transposed*/
  fprintf(stream, "%s\n", "Position Scoring matrix");
  fprintf(stream, "    ");
  for (i=0; i <  PS->len; i++)
    fprintf(stream, " %3d ", i);
  fprintf(stream, "\n");

  for (j=0; j < A_SIZE; j++)
    {
      if (PS->ptable->partition_table[j] == -1)
	continue;
      fprintf(stream, " %3c ", j+64);
      for (i=0; i < PS->len; i++)
	fprintf(stream, "%3d ", PS->M[PM_M(i,j)]);
      fprintf(stream, "\n");
    }
      

  /* Print pM - transposed for easier viewing */

  fprintf(stream, "\n\n");
  fprintf(stream, "%s\n", "Scores w.r.t. partitions");
  fprintf(stream, "      ");

  for (i=0; i <  PS->len; i++)
    fprintf(stream, " %3d ", i);
  fprintf(stream, "\n");

  for (j=0; j < P_SIZE; j++)
    {
      fprintf(stream, " %5d ", j);
      for (i=0; i < PS->len; i++)
	fprintf(stream, "%3d ", PS->pM[PM_pM(i,j)]);
      fprintf(stream, "\n");
    }
      
  /* Print pMclosest - transposed for easier viewing */
  fprintf(stream, "\n\n");
  fprintf(stream, "%s\n", "Best scores w.r.t. other partitions");
  fprintf(stream, "   ");

  for (i=0; i <  PS->len; i++)
    fprintf(stream, " %3d ", i);
  fprintf(stream, "\n");

  fprintf(stream, "%2s ", "");
  for (i=0; i < PS->len; i++)
    fprintf(stream, "%3d ", PS->pMclosest[i]);

  fprintf(stream, "\n");
  fprintf(stream, "\n");
}


/* Similarities to Distances */
POS_MATRIX *POS_MATRIX_S_2_D(POS_MATRIX *PS, BIOSEQ *query)
{
  int i = 0;
  int j = 0;
  int k = 0;
  POS_MATRIX *PD = mallocec(sizeof(POS_MATRIX));  
  int maxS;
  
  query->len = PS->len;
  query->start = mallocec(PS->len + 1); 
  query->start[PS->len] = '\0';

  PD->len = PS->len;
  PD->ptable = PS->ptable;
  PD->similarity_flag = 0;
  PD->M = mallocec(PS->len * A_SIZE * sizeof(SSINT));
  PD->pM = mallocec(PS->len * P_SIZE * sizeof(SSINT));
  PD->pMclosest = mallocec(PS->len * sizeof(SSINT));

  for (i=0; i < PS->len; i++)
    {
      maxS = POS_MATRIX_max_entry_pos(PS, i, &k);
      query->start[i] = 'A' + k;
      for (j=0; j < A_SIZE; j++)
	PD->M[PM_M(i,j)] = maxS -  PS->M[PM_M(i,j)];

      for (j=0; j < P_SIZE; j++)
	PD->pM[PM_pM(i,j)] = maxS - PS->pM[PM_pM(i,j)];

      PD->pMclosest[i] = maxS  - PS->pMclosest[i]; 
      
    }      
  return PD;
}
