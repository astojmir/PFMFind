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
#include "hit_list.h"
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
  SSINT SS[A_SIZE];
  const char *filename;
  char similarity_flag;
  FS_PARTITION_t *ptable;
  double mean;
  double var;
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
MY_INLINE
void SCORE_MATRIX_get_meanvar(SCORE_MATRIX_t *S, double *mean, 
			      double *var);

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
void SCORE_MATRIX_convert(int s0, HIT_LIST_t *HL);

int Davg_2_S(int Davg, int Sxx, int Syy);
int Dquasi_2_S(int Dquasi, int Sxx);

/* Scores and p-values */
int SCORE_MATRIX_Gaussian_cutoff(SCORE_MATRIX_t *S, BIOSEQ *query,
				 double *bkgrnd, double pcutoff);

 /* Evaluation of similarities, distances */
MY_INLINE
int SCORE_MATRIX_evaluate(SCORE_MATRIX_t *S, BIOSEQ *query,
			  BIOSEQ *subject);
MY_INLINE
int SCORE_MATRIX_evaluate_min(SCORE_MATRIX_t *S, 
			      BIOSEQ *query, BIOSEQ *subject,
				int Tmin, int *value);
MY_INLINE
int SCORE_MATRIX_evaluate_max(SCORE_MATRIX_t *S, 
			      BIOSEQ *query, BIOSEQ *subject,
			      int Tmax, int *value);
MY_INLINE
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

/* Set members */

#define SCORE_MATRIX_set_sim_flag(S, sim_flag) \
        ((S)->similarity_flag = (sim_flag))

#define SCORE_MATRIX_set_ptable(S, ptable) \
        ((S)->ptable = (ptable))


/* Get entries */ 

#define SCORE_MATRIX_entry(S, row, col) \
        ((S)->M[(row)][(col)])

#define SCORE_MATRIX_p_entry(S, row, col) \
        ((S)->pM[(row)][(col)])

MY_INLINE
void SCORE_MATRIX_get_meanvar(SCORE_MATRIX_t *S, double *mean, 
			      double *var)
{
  *mean = S->mean;
  *var = S->var;
}

#define SCORE_MATRIX_filename(S) \
        ((S)->filename)


/* Similarities to Distances, Quasi-metrics ... */

#define Davg_2_S(Davg, Sxx, Syy) \
        (((Sxx) + (Syy) - (Davg))/2)

#define Dquasi_2_S(Dquasi, Sxx) \
        ((Sxx) - (Dquasi))


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

#endif /* #ifndef _SMATRIX_H */
