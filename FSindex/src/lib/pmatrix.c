/********************************************************************/    
/*                                                                  */
/*                     SCORE_MATRIX module                          */ 
/*                                                                  */
/********************************************************************/    

#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "partition.h"
#include "smatrix.h"
#include "fastadb.h"

#ifdef DEBUG
#define POS_MATRIX_INLINE
#endif

#include "pmatrix.h"

#ifdef DEBUG
#undef POS_MATRIX_INLINE
#endif



static inline
void create_counts(POS_MATRIX *PS)
{
  int i;
  int j;


/* We try to use memset instead 
  for (i=0; i < len; i++) 
    for (j=0; j < A_SIZE; j++) 
      count[] = 0.0; 
*/
  memset(PS->count, 0, PS->len * A_SIZE * sizeof(double));

  for (j=0; j < PS->no_seqs; j++)
    for (i=0; i < PS->len; i++)
      {
	PS->count[PM_M(i, *(PS->seq[j].start+i) & A_SIZE_MASK)] 
	  += PS->weight[j];
	PS->count[PM_M(i, 0)] += PS->weight[j];
      }
}

static
void print_counts(POS_MATRIX *PS, FILE *stream)
{
  int i, j;
  

  fprintf(stream, "%s\n", "Amino Acid Counts");
  fprintf(stream, "    ");
  for (i=0; i <  PS->len; i++)
    fprintf(stream, " %5d ", i);
  fprintf(stream, "\n");

  for (j=0; j < A_SIZE; j++)
    {
      if (PS->ptable->partition_table[j] == -1)
	continue;
      fprintf(stream, " %3c ", j+64);
      for (i=0; i < PS->len; i++)
	fprintf(stream, "%6.3g ", PS->count[PM_M(i,j)]);
      fprintf(stream, "\n");
    }
  fprintf(stream, "TOT. ");
  for (i=0; i < PS->len; i++)
    fprintf(stream, "%6.3g ", PS->count[PM_M(i,0)]);
  fprintf(stream, "\n");
  fprintf(stream, "\n");
      

}
static
void print_M(POS_MATRIX *PS, FILE *stream)
{
  int i, j;

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
	fprintf(stream, "%4d ", PS->M[PM_M(i,j)]);
      fprintf(stream, "\n");
    }
}

static inline
void compute_log_odds(POS_MATRIX *PS)
{
  int i;
  int j;
  ULINT len = PS->len;
  double lambda = PS->lambda;
  double *bkgrnd = PS->bkgrnd;
  double *count = PS->count;
  double score;

  for (i=0; i < len; i++) 
    for (j=1; j < A_SIZE; j++) 
      if (bkgrnd[j] != 0.0)
	{
	  score = lambda * log(count[PM_M(i, j)] / 
			       (count[PM_M(i, 0)] * bkgrnd[j]));
	  PS->M[PM_M(i, j)] = (SSINT) rint(score);
	}
}

#define BUF_SIZE 256 


void POS_MATRIX_init(POS_MATRIX *PS, double lambda,
		     pseudo_counts_func *pcnf, 
		     const char *freq_filename)
{
  FILE *stream;
  ULINT freq;
  ULINT total = 0;
  char buffer[BUF_SIZE];
  int i;
  int n; 
  char c;


  PS->lambda = lambda;

  if (PS->count == NULL)
    PS->count = mallocec(PS->len * A_SIZE * sizeof(double));

  if (PS->bkgrnd == NULL)
    PS->bkgrnd = mallocec(A_SIZE * sizeof(double));
 
  if (PS->query == NULL)
    {
      PS->query = mallocec(sizeof(BIOSEQ));
      PS->query->start = mallocec(PS->len + 1);
      PS->query->start[PS->len] = '\0';
      PS->query->id.defline = strdup("PSM");
    }

  PS->pcnf = pcnf;
  
  if (freq_filename == NULL)
    stream = stdin;
  else if ((stream = fopen(freq_filename, "r")) == NULL)
    {
      fprintf(stderr, 
	      "POS_MATRIX_init(): Could not frequency file %s!\n", 
	      freq_filename);
      exit(EXIT_FAILURE);
    }

  /* Skip first lines */
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  sscanf(buffer+15, "%d", &n); 
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);

  /* Now get frequencies */

  for (i=0; i < n; i++)
    {
      fgets(buffer, BUF_SIZE, stream);
      sscanf(buffer+15, "%c", &c);
      if (PS->ptable->partition_table[c & A_SIZE_MASK] == -1)
	PS->bkgrnd[c & A_SIZE_MASK] = 0.0;
      else
	{
	  sscanf(buffer+25, "%ld", &freq);
	  PS->bkgrnd[c & A_SIZE_MASK] = (double) freq;
	  total += freq;
	}
    }

  for (i=1; i < A_SIZE; i++)
    PS->bkgrnd[i] /= total;   

  fclose(stream);
}

/* Routines for computing distances to pattern letters (i.e. subsets
   of alphabet partitions. - same as in smatrix.h */

static int SM_psize;          /* Current partition size */
static int SM_poffset;        /* Current partition offset */
static POS_MATRIX *SM_S;      /* Score matrix    */
static int SM_i;              /* Current sequence position */
static int SM_j;              /* Current Partition */

static
void pattern_similarity(int Score, int P, int pos)
{
  int letter;
  if (pos == SM_psize)
    {
      P += SM_poffset;
      SM_S->pM[PM_pM(SM_i, P)] = Score;
      return;
    }
  pos++;
  pattern_similarity(Score, P, pos);

  P = P | (1 << (pos-1));
  letter = FS_PARTITION_get_letter(SM_S->ptable, SM_j, pos-1); 
  letter = letter & A_SIZE_MASK;
  if (Score < SM_S->M[PM_M(SM_i,letter)])
    Score = SM_S->M[PM_M(SM_i,letter)];
  pattern_similarity(Score, P, pos);
  return;
}


void POS_MATRIX_update(POS_MATRIX *PS)
{
  int j;
  int k;
  int j_offset;
  int maxS;

  /* Update score matrix */

  PS->pcnf(PS);

#ifdef DEBUG
  print_counts(PS, stdout);
#endif
  compute_log_odds(PS);
#ifdef DEBUG
  print_M(PS, stdout);
#endif

  /* Compute maximum similarities to pattern letters (subsets of
     alphabet partitions) */
  for (SM_i=0; SM_i < PS->len; SM_i++)
    {      
      for (SM_j=0; SM_j < PS->ptable->no_partitions; SM_j++)
	{
	  SM_S = PS;
	  SM_psize = FS_PARTITION_get_size(PS->ptable, SM_j);
	  SM_poffset =  FS_PARTITION_get_poffset(PS->ptable, SM_j);

	  /* Now traverse the binary tree representing the subsets and
	     compute maximum similarities */
	  pattern_similarity(SHRT_MIN, 0, 0);
	  SM_S->pM[PM_pM(SM_i, SM_poffset)] =
	    SM_S->pM[PM_pM(SM_i, SM_poffset + (1 << SM_psize) - 1)];
	}
      
      /* Convert similarities to distances */
      maxS = POS_MATRIX_max_entry_pos(PS, SM_i, &k);
      PS->query->start[SM_i]= 64 + k;

      for (j=0; j < P_SIZE; j++)
	PS->pM[PM_pM(SM_i,j)] = maxS - PS->pM[PM_pM(SM_i,j)];

      /* Compute minimum distances to any other cluster */
      PS->pMclosest[SM_i] = 0;
      for (j=0; j < PS->ptable->no_partitions; j++)
	{
	  if (PS->ptable->partition_table[k] != j)
	    {
	      j_offset = FS_PARTITION_get_poffset(PS->ptable, j); 
	      PS->pMclosest[SM_i] = 
		maxS - max(PS->pMclosest[SM_i], 
			   PS->pM[PM_pM(SM_i, j_offset)]); 
			  
	    }
	}
    }     
}



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

  PS->M = mallocec(PS->len * A_SIZE * sizeof(SSINT));
  PS->pM = mallocec(PS->len * P_SIZE * sizeof(SSINT));
  PS->pMclosest = mallocec(PS->len * sizeof(SSINT));

  PS->query = mallocec(sizeof(BIOSEQ));
  PS->query->start = mallocec(PS->len + 1);
  memcpy(PS->query->start, query->start, PS->len + 1);
  PS->query->id.defline = strdup("PSM");
  PS->query->len = query->len;

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

  POS_MATRIX_S_2_D(PS);
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

/* Pseudo-count functions */
void POS_MATRIX_simple_pseudo_counts(POS_MATRIX *PS)
{
  double *A = (double *) PS->params;
  int i;
  int j;

  for (i=0; i < PS->len; i++) 
    for (j=1; j < A_SIZE; j++) 
      if (PS->ptable->partition_table[j] != -1)
	{
	  PS->count[PM_M(i, j)] += *A * PS->bkgrnd[j];
	  PS->count[PM_M(i, 0)] += *A * PS->bkgrnd[j];  
	}
}


void POS_MATRIX_simple_pseudo_counts_init(POS_MATRIX *PS, double A)
{
  double *pA = mallocec(sizeof(double));

  *pA = A;
  PS->params = pA;
}


/* Weight functions */
void POS_MATRIX_equal_weights(POS_MATRIX *PS)
{
  int i;
  PS->weight = reallocec(PS->weight, PS->max_no_seqs * sizeof(double));
  for (i=0; i < PS->no_seqs; i++)
    PS->weight[i] = 1.0;
  create_counts(PS);
#ifdef DEBUG
  print_counts(PS, stdout);
#endif
}

void POS_MATRIX_Henikoff_weights(POS_MATRIX *PS)
{
  int i, j, k;
  int r=0;
  double *bkgrnd = PS->bkgrnd;
  double *count = PS->count;
  double w;
  double w_sum;

  POS_MATRIX_equal_weights(PS);

  memset(PS->weight, 0, PS->no_seqs * sizeof(double));

  w_sum = 0.0;
  for (i=0; i < PS->len; i++)
    {
      for (j=1; j < A_SIZE; j++) 
	if (bkgrnd[j] != 0.0 && count[PM_M(i, j)] > 0.0)
	  r++;
      for (k=0; k < PS->no_seqs; k++)
	{
	  w = 1/(r * PS->count[PM_M(i, *(PS->seq[k].start+i) & A_SIZE_MASK)]);
	  PS->weight[k] += w;
	  w_sum += w;
	}
    }
  w_sum = (double) PS->no_seqs / w_sum;

  for (k=0; k < PS->no_seqs; k++)
    PS->weight[k] *= w_sum;

  create_counts(PS);
#ifdef DEBUG
  /* Print Weights */
  fprintf(stdout, "SEQUENCE WEIGHTS\n");
  for (k=0; k < PS->no_seqs; k++)
    {
      fprintf(stdout, "%.*s %6.2f\n", (int) PS->len, PS->seq[k].start,
	      PS->weight[k]);
    }  
  fprintf(stdout, "\n\n");
  print_counts(PS, stdout);
#endif
}

/* Read, Write */
int POS_MATRIX_write(POS_MATRIX *PS, FILE *stream) 
{
  fwrite(&(PS->len), sizeof(ULINT), 1, stream);   
  fwrite(PS->M, sizeof(SSINT), PS->len * A_SIZE, stream);   
  fwrite(PS->pM, sizeof(SSINT), PS->len * P_SIZE, stream);   
  fwrite(PS->pMclosest, sizeof(SSINT), PS->len, stream);   
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
void POS_MATRIX_S_2_D(POS_MATRIX *PS)
{
  int i = 0;
  int j = 0;
  int k = 0;
  int maxS;
  
  for (i=0; i < PS->len; i++)
    {
      maxS = POS_MATRIX_max_entry_pos(PS, i, &k);
      PS->query->start[i] = 64 + k;

      for (j=0; j < P_SIZE; j++)
	PS->pM[PM_pM(i,j)] = maxS - PS->pM[PM_pM(i,j)];

      PS->pMclosest[i] = maxS  - PS->pMclosest[i]; 
    }      
}
