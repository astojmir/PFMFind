/********************************************************************/    
/*                                                                  */
/*                     SCORE_MATRIX module                          */ 
/*                                                                  */
/********************************************************************/    

#ifndef _SMATRIX_H
#define _SMATRIX_H

#include "partition.h"
#include <stdio.h>
#include "misclib.h"


/********************************************************************/    
/*                                                                  */
/* This module implements the scoring matrices, both the simple     */
/* (alphabet) and positional (profiles). All functions (except      */
/* constructors, destructors or similar) are implemented as member  */
/* function pointers so that different functions can be assigned    */
/* according to the type.                                           */
/*                                                                  */
/* Members:                                                         */
/*                                                                  */
/* - Mtype - matrix type {SCORE, POSITIONAL}                        */
/* - Stype - score type  {SIMILARITY, DISTANCE}                     */
/* - len   - length (set only if Mtype == POSITIONAL)               */
/* - M     - actual values of the matrix (A_SIZE by A_SIZE if       */
/*           Mtype == SCORE, len by A_SIZE if Mtype == POSITIONAL   */
/* - qseq  - only if Mtype == POSITIONAL: a sequence with maximal   */
/*           similarity or zero distance; may not be unique; length */
/*           is len.                                                */
/* - alphabet - alphabet used.                                      */
/*                                                                  */
/*                                                                  */
/*                                                                  */
/* Functions:                                                       */
/*                                                                  */
/* - eval_score(S, s1, s2, len2) - evaluates the score between      */
/*     sequences s1 and s2, of length len2. If Mtype == POSITIONAL  */
/*     s2 is ignored and the min(len, len2) is used for score       */
/*     evaluation.                                                  */
/* - set_conv_type(S, conv_type) - sets other conversion functions  */
/*     according to the conv_type.                                  */
/* - item_conv(S, c, i, j) - converts the item M[i][j] if Mtype ==  */
/*     POSITIONAL or M[c & A_SIZE_MASK][j] if Mtype == SCORE.       */
/* - range_conv(S, q, len1, r) - converts the range r. q : query    */
/*     sequence, len1 its length - ignored if Mtype == POSITIONAL.  */
/* - matrix_conv(S, q, len) - produces new matrix of desired type.  */
/*                                                                  */
/* Conversion type flag:                                            */
/*                                                                  */
/* Bitwise flag: first bit set (POSITIONAL) gives conversion to     */
/*   positional matrix. Next two bits {QUASI, MAX, AVG} give ways   */
/*   of conversion to distances. The parts that are inappropriate   */
/*   are ignored.                                                   */
/*                                                                  */
/* Conversions and cases:                                           */
/*                                                                  */
/* 1. Mtype == POSITIONAL && Stype == DISTANCE                      */
/*      No conversion performed for any conv_type.                  */
/*                                                                  */
/* 2. Mtype == SCORE && Stype == DISTANCE                           */
/*      No conversion by default. Can convert to POSITIONAL matrix. */
/*                                                                  */
/* 3. Mtype == POSITIONAL && Stype == SIMILARITY                    */
/*      By default, conv_flag is POS        + QUASI. This is the    */
/*      only possible conversion.                                   */
/*                                                                  */
/* 4. Mtype == SCORE && Stype == SIMILARITY                         */
/*      By default, conv_flag == QUASI. All other combinations are  */
/*      allowed but range_conv is identity (no conversion possible) */
/*      if MAX or AVG conversion is set. Thus, such radii are       */
/*      interpreted as distance radii.                              */
/*                                                                  */
/********************************************************************/    

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               PROTOTYPES                                     ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    
/* This is used so that '\0' maps to a very large distance and
   so the sequences of shorter length get rejected */
#define VERY_LARGE_DISTANCE 5000

typedef enum {SCORE, POSITIONAL} MATRIX_TYPE;
typedef enum {SIMILARITY, DISTANCE} SCORE_TYPE;

#define POS 1
#define QUASI 2
#define MAX 4
#define AVG 6

typedef struct SM_s
{
  MATRIX_TYPE Mtype;
  SCORE_TYPE Stype;
  int len;
  int **M;
  char *qseq;
  char *alphabet;
  int alen;
  int conv_type;
  int (*eval_score) (struct SM_s *, const char *, const char *, int);
  int (*item_conv) (struct SM_s *, char, int, int);
  int (*range_conv) (struct SM_s *, const char *, int, int);
  void (*set_conv_type) (struct SM_s *, int);
  struct SM_s * (*matrix_conv) (struct SM_s *, const char *, int);
} SCORE_MATRIX;


/* Constructors and destructor */
SCORE_MATRIX *SCORE_MATRIX_init(int **M, int len, MATRIX_TYPE Mtype,
				SCORE_TYPE Stype, char *alphabet);

SCORE_MATRIX *SCORE_MATRIX_copy(SCORE_MATRIX *S);

SCORE_MATRIX *SCORE_MATRIX_from_file(const char *filename);

void SCORE_MATRIX_del(SCORE_MATRIX *S);

void SCORE_MATRIX_fprint(SCORE_MATRIX *S, FILE *fp);

char *SCORE_MATRIX_sprint(SCORE_MATRIX *S);

#endif /* #ifndef _SMATRIX_H */
