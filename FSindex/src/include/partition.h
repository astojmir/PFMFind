/********************************************************************/    
/*                                                                  */
/*                    FS_TABLE module                               */ 
/*                                                                  */
/********************************************************************/    
#ifndef _FS_TABLE_H
#define _FS_TABLE_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "misclib.h"

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               PROTOTYPES                                     ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

/* Predefined sizes - powers of two:
   A_SIZE - maximum size of alphabet;
*/

#define A_SIZE 32
#define A_SIZE_MASK 31

/* We store SEQ_index list in FS_BIN */
typedef ULINT SEQ_index_t;

typedef int FS_SEQ_t;

typedef struct
{
  int len;                        /* Sequence length */
  char **sepn;                    /* Alphabet divided in partitions by
				     separator. e.g. ABC#DEF#GHIJK for
				     each sequence position */ 
  char *alphabet;                 /* Unseparated alphabet */
  int **pttn;                     /* Maps, for each position, 
				     each letter in the alphabet
				     (0 to A_SIZE_MASK - use &
				     A_SIZE_MASK with any character to
				     convert), to the partition it
				     belongs to  (-1 if none) */
  int **hash;                     /* Maps for each position, 
				     each partition to a bin number. The
				     bin number of a sequence is the sum 
				     of the bin numbers of its residues.
				  */
  int *K;                        /* For each position gives the number
				     of partitions. */
  int *nbpt;                     /* nbpt = K-1            */
  ULINT bins;                    /* Total number of possible reduced
				    sequences. */
  int **rank;                    /* Rank of each alphabet letter at each position */
  int same_order;                /* 1 if the alphabet order at each position is 
				    the same, 0 otherwise */
} FS_TABLE;

/* Main constructor */

FS_TABLE *FS_TABLE_init(const char **sepn, int len); 

/* Destructor */

void FS_TABLE_del(FS_TABLE *ptable);

/* I/O */
void FS_TABLE_write(FS_TABLE *ptable, FILE *fp);
FS_TABLE *FS_TABLE_read(FILE *fp);

/* Sequence conversion and validation */

MY_INLINE
FS_SEQ_t FS_SEQ(FS_TABLE *ptable, const char *seq, int len);

char *FS_SEQ_print(FS_TABLE *ptable, FS_SEQ_t FSseq);

/* Printing */
char *FS_TABLE_sprint(FS_TABLE *ptable);
void FS_TABLE_fprint(FS_TABLE *ptable, FILE *fp);

/* Sequence database (heap) conversion into array of ints */
int *FS_TABLE_order_convert_heap(FS_TABLE *ptable, const char *s, int len);

/* Access macros - we will access it directly for indexing.
   No bounds checking either. */

int FS_TABLE_get_len(FS_TABLE *ptable);

char **FS_TABLE_get_sepn(FS_TABLE *ptable);

char *FS_TABLE_get_alphabet(FS_TABLE *ptable);

int FS_TABLE_get_K(FS_TABLE *ptable, int pos);

int FS_TABLE_get_pttn(FS_TABLE *ptable, int pos, char letter);

int FS_TABLE_get_hash(FS_TABLE *ptable, int pos, int pttn);


/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS                    ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

/********************************************************************/    
/*                                                                  */
/*                    FS_TABLE module                               */ 
/*                                                                  */
/********************************************************************/    

#define FS_TABLE_get_len(ptable) ((ptable)->len)

#define FS_TABLE_get_sepn(ptable) ((ptable)->sepn)

#define FS_TABLE_get_alphabet(ptable) ((ptable)->alphabet)

#define FS_TABLE_get_pttn(ptable, pos, letter) \
        ((ptable)->pttn[(pos)][(letter) & A_SIZE_MASK])

#define FS_TABLE_get_K(ptable, pos) ((ptable)->K[(pos)])

#define FS_TABLE_get_hash(ptable, pos, pttn) \
        ((ptable)->hash[(pos)][(pttn)])


MY_INLINE 
FS_SEQ_t FS_SEQ(FS_TABLE *ptable, const char *seq, int len)
{
  int i;
  int p;
  int l = len;
  ULINT FSseq;

  FSseq = 0;
  if (len > ptable->len)
    l = ptable->len;
  for (i=0; (i < l) && *seq; i++, seq++) {
    p = ptable->pttn[i][*seq & A_SIZE_MASK];
    if (p == -1) {
      FSseq = -1;
      break;
    }
    FSseq += ptable->hash[i][p];
  }
  return FSseq;
}

#endif /* #ifndef _FS_TABLE_H */
