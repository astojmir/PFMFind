/********************************************************************/    
/*                                                                  */
/*                     SCORE_MATRIX module                          */ 
/*                                                                  */
/********************************************************************/    

#ifndef _SMATRIX_H
#define _SMATRIX_H

#include "misclib.h"
#include <stdio.h>
#include "fastadb.h"
#include "partition.h"
#ifdef USE_MPATROL
#include <mpatrol.h>
#endif


/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               PROTOTYPES                                     ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

typedef struct 
{
  SSINT M[A_SIZE][A_SIZE];
  SSINT pM[A_SIZE][P_SIZE];
  SSINT pMclosest[A_SIZE];
  char similarity_flag;
  FS_PARTITION_t *ptable;  
} SCORE_MATRIX_t;

extern int SCORE_MATRIX_VERBOSE;

/* Main constructor */
SCORE_MATRIX_t *SCORE_MATRIX_create(const char *filename,
				    FS_PARTITION_t *ptable); 
/* Destructor */
void SCORE_MATRIX_destroy(SCORE_MATRIX_t *Score_matrix);

/* Read, Write */
int SCORE_MATRIX_write(SCORE_MATRIX_t *Score_matrix, FILE *stream); 
SCORE_MATRIX_t *SCORE_MATRIX_read(FILE *stream);

/* Set members */
void SCORE_MATRIX_set_sim_flag(SCORE_MATRIX_t *Score_matrix, 
			       char sim_flag);
void SCORE_MATRIX_set_ptable(SCORE_MATRIX_t *Score_matrix, 
			     FS_PARTITION_t *ptable);

/* Element Access + Properties */
int SCORE_MATRIX_max_entry(SCORE_MATRIX_t *Score_matrix);
int SCORE_MATRIX_min_entry(SCORE_MATRIX_t *Score_matrix);
int SCORE_MATRIX_max_entry_col(SCORE_MATRIX_t *Score_matrix, 
			       int col);
int SCORE_MATRIX_max_entry_row(SCORE_MATRIX_t *Score_matrix, 
			       int row);
int SCORE_MATRIX_min_entry_col(SCORE_MATRIX_t *Score_matrix, 
			       int col);
int SCORE_MATRIX_min_entry_row(SCORE_MATRIX_t *Score_matrix, 
			       int row);
int SCORE_MATRIX_entry(SCORE_MATRIX_t *Score_matrix,
		       int row, int col); 
int SCORE_MATRIX_p_entry(SCORE_MATRIX_t *Score_matrix,
			 int row, int col);
 
/* Similarities to Distances, Quasi-metrics ... */
SCORE_MATRIX_t *SCORE_MATRIX_S_2_Dmax(SCORE_MATRIX_t *S);
SCORE_MATRIX_t *SCORE_MATRIX_S_2_Davg(SCORE_MATRIX_t *S);
SCORE_MATRIX_t *SCORE_MATRIX_S_2_Dquasi(SCORE_MATRIX_t *S);
int Davg_2_S(int Davg, int Sxx, int Syy);
int Dquasi_2_S(int Dquasi, int Sxx);

 /* Evaluation of similarities, distances */
int SCORE_MATRIX_evaluate(SCORE_MATRIX_t *S, BIOSEQ *query,
			  BIOSEQ *subject);
int SCORE_MATRIX_evaluate_min(SCORE_MATRIX_t *S, 
			      BIOSEQ *query, BIOSEQ *subject,
				int Tmin, int *value);
int SCORE_MATRIX_evaluate_max(SCORE_MATRIX_t *S, 
			      BIOSEQ *query, BIOSEQ *subject,
			      int Tmax, int *value);
int SCORE_MATRIX_verify_pos(SCORE_MATRIX_t *D, BIOSEQ *query,
			    USINT *TT, int k, int cutoff);

/* Printing */
void SCORE_MATRIX_print(SCORE_MATRIX_t *S, FILE *stream, 
			const char *title);

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS                    ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    
#ifndef DEBUG

#ifndef MY_INLINE
#define MY_INLINE extern inline
#endif

#define SMATRIX_INLINE

#else

#ifndef MY_INLINE
#define MY_INLINE 
#endif

#endif /* #ifndef DEBUG */

#ifdef SMATRIX_INLINE
/* Set members */

MY_INLINE
void SCORE_MATRIX_set_sim_flag(SCORE_MATRIX_t *Score_matrix, 
			       char sim_flag)
{
  Score_matrix->similarity_flag = sim_flag;
}


MY_INLINE
void SCORE_MATRIX_set_ptable(SCORE_MATRIX_t *Score_matrix, 
			     FS_PARTITION_t *ptable)
{
  Score_matrix->ptable = ptable;
}

/* Get entries */ 


MY_INLINE
int SCORE_MATRIX_entry(SCORE_MATRIX_t *Score_matrix,
		       int row, int col) 
{
  return Score_matrix->M[row][col];
}


MY_INLINE
int SCORE_MATRIX_p_entry(SCORE_MATRIX_t *Score_matrix,
			 int row, int col)
{
  return Score_matrix->pM[row][col];
}

/* Similarities to Distances, Quasi-metrics ... */


MY_INLINE
int Davg_2_S(int Davg, int Sxx, int Syy)
{
  return (Sxx + Syy - Davg)/2;
}


MY_INLINE
int Dquasi_2_S(int Dquasi, int Sxx)
{
  return Sxx - Dquasi;
}

 /* Evaluation of similarities, distances */

MY_INLINE
int SCORE_MATRIX_evaluate(SCORE_MATRIX_t *S, BIOSEQ *query,
			  BIOSEQ *subject)
{
  UINT_t i = query->len - 1;
  int H = 0;
  char *q = query->start;
  char *s = subject->start;
#if 0
  for(; i--; )
#endif
  for(i=0; i < query->len; i++)    
    H += S->M[*(q+i) & A_SIZE_MASK][*(s+i) & A_SIZE_MASK]; 
  return H;
}

MY_INLINE
int SCORE_MATRIX_evaluate_min(SCORE_MATRIX_t *S, 
			      BIOSEQ *query, BIOSEQ *subject,
			      int Tmin, int *value)
{
  /* For now no particular optimisation */
  if ((*value = SCORE_MATRIX_evaluate(S, query, subject)) >= Tmin)
    return 1;
  else
    return 0;
}


MY_INLINE
int SCORE_MATRIX_evaluate_max(SCORE_MATRIX_t *S, 
			      BIOSEQ *query, BIOSEQ *subject,
			      int Tmax, int *value)
{
  /* For now no particular optimisation */
  if ((*value = SCORE_MATRIX_evaluate(S, query, subject)) <= Tmax)
    return 1;
  else
    return 0;
}

MY_INLINE
int SCORE_MATRIX_verify_pos(SCORE_MATRIX_t *D, BIOSEQ *query,
			    USINT *TT, int k, int cutoff)
{
  int i;
  int Sum = 0;

  
  for (i = 0; i < k; i++)
    Sum += D->pMclosest[query->start[TT[i]] & A_SIZE_MASK];

  if (Sum <= cutoff)
    return 1;
  else
    return 0;
}

#endif /* #ifdef SMATRIX_INLINE */

#endif /* #ifndef _SMATRIX_H */
