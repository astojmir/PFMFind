/********************************************************************/    
/*                                                                  */
/*                     POS_MATRIX module                            */ 
/*                                                                  */
/********************************************************************/    

#ifndef _PMATRIX_H
#define _PMATRIX_H

#include "misclib.h"
#include <stdio.h>
#include "fastadb.h"
#include "partition.h"
#include "smatrix.h"
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
  SSINT *M;                 /* Pure alphabet matrix: len * A_SIZE */
  SSINT *pM;                /* Pattern matrix: len * P_SIZE */
  SSINT *pMclosest;         /* Changed pos array */
  ULINT len;
  char similarity_flag;
  FS_PARTITION_t *ptable;  
} POS_MATRIX;


/* Similarity to distance cutoff conversion function type */
typedef int POS_MATRIX_range_convert_func(POS_MATRIX *PS, 
					  BIOSEQ *query, 
					  int cutoff,
					  ULINT Pfrom, ULINT Pto);	


/* Constructors and converters */
POS_MATRIX *POS_MATRIX_load(const char *filename,
			    FS_PARTITION_t *ptable); 

POS_MATRIX *POS_MATRIX_create(ULINT no_seq, BIOSEQ *seq,
			      ULINT *Pfrom, ULINT *Pto,
			      FS_PARTITION_t *ptable); 

POS_MATRIX *SCORE_2_POS_MATRIX(SCORE_MATRIX_t *S, BIOSEQ *query);

/* Destructor */
void POS_MATRIX_destroy(POS_MATRIX *PS);

/* Read, Write */
int POS_MATRIX_write(POS_MATRIX *PS, FILE *stream); 
POS_MATRIX *POS_MATRIX_read(FILE *stream);

/* Printing */
void POS_MATRIX_print(POS_MATRIX *PS, FILE *stream, 
			const char *title);

/* Set members */
void POS_MATRIX_set_sim_flag(POS_MATRIX *PS, char sim_flag);
void POS_MATRIX_set_ptable(POS_MATRIX *PS, 
			   FS_PARTITION_t *ptable);

/* Element Access + Properties */
int POS_MATRIX_max_entry_pos(POS_MATRIX *PS, int pos, int *col);


int PM_M(int i, int j);
int PM_pM(int i, int j);


/* Similarities to Distances */
POS_MATRIX *POS_MATRIX_S_2_D(POS_MATRIX *PS, BIOSEQ *query);

/* Cutoff value conversion */

int POS_MATRIX_id_convert(POS_MATRIX *PS, BIOSEQ *query, int cutoff,
			  ULINT Pfrom, ULINT Pto); 

int POS_MATRIX_S2D_convert(POS_MATRIX *PS, BIOSEQ *query, int cutoff,
			   ULINT Pfrom, ULINT Pto); 
 /* Evaluation of similarities, distances */
int POS_MATRIX_evaluate(POS_MATRIX *PS, BIOSEQ *subject, 
			ULINT Pfrom, ULINT Pto);

int POS_MATRIX_evaluate_min(POS_MATRIX *PS, BIOSEQ *subject, 
			    ULINT Pfrom, ULINT Pto,
			    int Tmin, int *value); 

int POS_MATRIX_evaluate_max(POS_MATRIX *PS, BIOSEQ *subject,
			    ULINT Pfrom, ULINT Pto,
			    int Tmax, int *value);

int POS_MATRIX_verify_pos(POS_MATRIX *PD, ULINT Pfrom, ULINT Pto,
			  USINT *TT, int k, int cutoff);





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

#define POS_MATRIX_INLINE

#else

#ifndef MY_INLINE
#define MY_INLINE 
#endif

#endif /* #ifndef DEBUG */

#ifdef POS_MATRIX_INLINE
/* Set members */

MY_INLINE
void POS_MATRIX_set_sim_flag(POS_MATRIX *PS, char sim_flag)
{
  PS->similarity_flag = sim_flag;
}


MY_INLINE
void POS_MATRIX_set_ptable(POS_MATRIX *PS, FS_PARTITION_t *ptable)
{
  PS->ptable = ptable;
}


/* Element Access + Properties */
MY_INLINE
int PM_M(int i, int j)
{
  /* Remember A_SIZE = 32 (5 bits) */
  return (i << 5) + j;
}

MY_INLINE
int PM_pM(int i, int j)
{
  /* Remember P_SIZE = 256 (8 bits) */
  return (i << 8) + j;
}


MY_INLINE
int POS_MATRIX_max_entry_pos(POS_MATRIX *PS, int pos, int *col)
{
  int j;
  int maxS = 0; /* Must have at least one entry gt 0 */
  
  *col = 0;
  for (j=0; j < A_SIZE; j++)
    {
      if (maxS <  PS->M[PM_M(pos,j)])
	{
	  maxS = PS->M[PM_M(pos,j)];
	  *col = j;
	}
    }  
  return maxS;
}

/* Evaluation of similarities, distances */

MY_INLINE
int POS_MATRIX_evaluate(POS_MATRIX *PS, BIOSEQ *subject, 
			ULINT Pfrom, ULINT Pto)
{
  int i;
  int H = 0;
  char *s = subject->start;

  for(i=Pfrom; i <= Pto; i++, s++)    
    H += PS->M[PM_M(i, *s & A_SIZE_MASK)]; 
  return H;
}

MY_INLINE
int POS_MATRIX_evaluate_min(POS_MATRIX *PS, BIOSEQ *subject, 
			    ULINT Pfrom,ULINT Pto,
			    int Tmin, int *value)
{
  if ((*value = 
       POS_MATRIX_evaluate(PS, subject, Pfrom, Pto)) >= Tmin)  
    return 1;
  else
    return 0;
} 

MY_INLINE
int POS_MATRIX_evaluate_max(POS_MATRIX *PS, BIOSEQ *subject, 
			    ULINT Pfrom, ULINT Pto,
			    int Tmax, int *value)
{
  if ((*value = 
       POS_MATRIX_evaluate(PS, subject, Pfrom, Pto)) <= Tmax)  
    return 1;
  else
    return 0;
}


MY_INLINE
int POS_MATRIX_verify_pos(POS_MATRIX *PD, ULINT Pfrom, ULINT Pto,
			  USINT *TT, int k, int cutoff)
{
  int i;
  int Sum = 0;

  for (i = 0; i < k; i++)
    Sum += PD->pMclosest[Pfrom + TT[i]];

  if (Sum <= cutoff)
    return 1;
  else
    return 0;
}

/* Cutoff value conversion */
MY_INLINE
int POS_MATRIX_id_convert(POS_MATRIX *PS, BIOSEQ *query, int cutoff,
			  ULINT Pfrom, ULINT Pto) 
{
  return cutoff;
}

MY_INLINE
int POS_MATRIX_S2D_convert(POS_MATRIX *PS, BIOSEQ *query, int cutoff,
			   ULINT Pfrom, ULINT Pto) 
{
  return POS_MATRIX_evaluate(PS, query, Pfrom, Pto) - cutoff;
}




#endif /* #ifdef POS_MATRIX_INLINE */

#endif /* #ifndef _PMATRIX_H */
