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
extern int POS_MATRIX_VERBOSE;
extern FILE *POS_MATRIX_STREAM;


struct POS_MATRIX_s 
{
  SSINT *M;                 /* Pure alphabet matrix: len * A_SIZE */
  SSINT *pM;                /* Pattern matrix: len * P_SIZE */
  SSINT *pMclosest;         /* Changed pos array */
  SSINT *SS;                /* Positional self-similarities */
  SSINT qS;                 /* Query self-similarity        */
  ULINT len;
  FS_PARTITION_t *ptable;   /* Partition table */
  const char *filename;     /* Original matrix filename */
  double lambda;            /* Scaling factor */
  double *bkgrnd;           /* Background probabilities */
  double *count;            /* Weighted frequency counts */  
  void (*pcnf) (struct POS_MATRIX_s *);
  /* pseudo_counts_func *pcnf; */
  void *params;             /* Parameters of pseudo count functions */  
  BIOSEQ *query;
  int no_seqs;
  int max_no_seqs;
  BIOSEQ *seq;
  double *weight;
  int iteration;
  double mean;
  double var;
};

typedef struct POS_MATRIX_s POS_MATRIX;

typedef void pseudo_counts_func(POS_MATRIX *PS);


/* Similarity to distance cutoff conversion function type */
typedef int POS_MATRIX_range_convert_func(POS_MATRIX *PS, 
					  BIOSEQ *query, 
					  int cutoff,
					  ULINT Pfrom, ULINT Pto);	


/* Constructors and converters */
POS_MATRIX *POS_MATRIX_load(const char *filename,
			    FS_PARTITION_t *ptable); 

void POS_MATRIX_init(POS_MATRIX *PS, double lambda,
		     pseudo_counts_func *pcnf, 
		     const char *freq_filename);

void POS_MATRIX_update(POS_MATRIX *PS);

POS_MATRIX *SCORE_2_POS_MATRIX(SCORE_MATRIX_t *S, BIOSEQ *query);

void POS_MATRIX_convert(POS_MATRIX *PS, HIT_LIST_t *HL);


/* Destructor */
void POS_MATRIX_destroy(POS_MATRIX *PS);

/* Pseudo-count functions */

pseudo_counts_func POS_MATRIX_simple_pseudo_counts;
void POS_MATRIX_simple_pseudo_counts_init(POS_MATRIX *PS, double A);

/* Weight functions */
void POS_MATRIX_equal_weights(POS_MATRIX *PS);
void POS_MATRIX_Henikoff_weights(POS_MATRIX *PS);

/* Load background frequences */
double *load_bkgrnd_probs(const char *filename);


/* Read, Write */
int POS_MATRIX_write(POS_MATRIX *PS, FILE *stream); 
POS_MATRIX *POS_MATRIX_read(FILE *stream);

/* Scores and p-values */
int POS_MATRIX_Gaussian_cutoff(POS_MATRIX *S, double pcutoff);

/* Printing */
void POS_MATRIX_print(POS_MATRIX *PS, FILE *stream, 
			const char *title);

/* Set members */
void POS_MATRIX_set_ptable(POS_MATRIX *PS, 
			   FS_PARTITION_t *ptable);

/* Element Access + Properties */
MY_INLINE
int POS_MATRIX_max_entry_pos(POS_MATRIX *PS, int pos, int *col);

MY_INLINE
void POS_MATRIX_get_meanvar(POS_MATRIX *PS, double *mean, 
			    double *var); 

int PM_M(int i, int j);
int PM_pM(int i, int j);


/* Similarities to Distances */
void POS_MATRIX_S_2_D(POS_MATRIX *PS);


/* Cutoff value conversion */

MY_INLINE
int POS_MATRIX_id_convert(POS_MATRIX *PS, BIOSEQ *query, int cutoff,
			  ULINT Pfrom, ULINT Pto); 

MY_INLINE
int POS_MATRIX_S2D_convert(POS_MATRIX *PS, BIOSEQ *query, int cutoff,
			   ULINT Pfrom, ULINT Pto);
 
 /* Evaluation of similarities, distances */
MY_INLINE
int POS_MATRIX_evaluate(POS_MATRIX *PS, BIOSEQ *subject, 
			ULINT Pfrom, ULINT Pto);
MY_INLINE
int POS_MATRIX_evaluate_min(POS_MATRIX *PS, BIOSEQ *subject, 
			    ULINT Pfrom, ULINT Pto,
			    int Tmin, int *value); 
MY_INLINE
int POS_MATRIX_evaluate_max(POS_MATRIX *PS, BIOSEQ *subject,
			    ULINT Pfrom, ULINT Pto,
			    int Tmax, int *value);
MY_INLINE
int POS_MATRIX_verify_pos(POS_MATRIX *PD, ULINT Pfrom, ULINT Pto,
			  USINT *TT, int k, int cutoff);





/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS                    ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

/* Set members */
#define POS_MATRIX_set_ptable(PS, ptable) \
        ((PS)->ptable = (ptable))

/* Remember A_SIZE = 32 (5 bits) */ 
#define PM_M(i, j) \
       (((i) << 5) + (j))

/* Remember P_SIZE = 256 (8 bits) */ 
#define PM_pM(i, j) \
        (((i) << 8) + (j))

#define POS_MATRIX_filename(PS) \
        ((PS)->filename)



MY_INLINE
void POS_MATRIX_get_meanvar(POS_MATRIX *PS, double *mean, 
			    double *var)
{
  *mean = PS->mean;
  *var = PS->var;
}

MY_INLINE
int POS_MATRIX_max_entry_pos(POS_MATRIX *PS, int pos, int *col)
{
  int j;
  int maxS = 0; /* Must have at least one entry gt 0 */
  
  *col = 0;
  for (j=0; j < A_SIZE; j++)
    {
      if (PS->ptable->partition_table[j] == -1)
	continue;
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

#endif /* #ifndef _PMATRIX_H */
