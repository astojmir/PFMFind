/********************************************************************/    
/*                                                                  */
/*                    FS_PARTITION module                           */ 
/*                                                                  */
/********************************************************************/    
#ifndef _FS_PARTITION_H
#define _FS_PARTITION_H

#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fastadb.h"
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

/* Predefined sizes - must be powers of two:
   A_SIZE - maximum size of alphabet;
   P_SIZE - maxumum size of pattern alphabet.
*/
#define A_SIZE 32
#define A_SIZE_MASK 31
#define P_SIZE 256

/* We store SEQ_index list in FS_BIN */
typedef ULINT SEQ_index_t;


typedef ULINT FS_SEQ_t;

typedef struct
{
  ULINT no_partitions;
  SLINT partition_table[A_SIZE];
} FS_PARTITION_t;

extern int FS_PARTITION_VERBOSE;

/* Main constructor */
FS_PARTITION_t *FS_PARTITION_create(const char *alphabet, 
				    char separator);

/* Destructor */
void FS_PARTITION_destroy(FS_PARTITION_t *FS_partition);

/* Read, Write */
int FS_PARTITION_write(FS_PARTITION_t *FS_partition, FILE *stream); 
FS_PARTITION_t *FS_PARTITION_read(FILE *stream);

/* Access and conversion */
ULINT FS_PARTITION_get_no_partitions(FS_PARTITION_t *FS_partition);
int BIOSEQ_2_FS_SEQ(BIOSEQ *seq, FS_PARTITION_t *FS_partition,
		    FS_SEQ_t *FS_seq);
char *FS_seq_print(FS_SEQ_t FS_seq, FS_PARTITION_t *FS_partition,
		   ULINT frag_len);
int FS_PARTITION_check_seq(BIOSEQ *seq, FS_PARTITION_t *FS_partition);

/* Printing */
void FS_PARTITION_print(FS_PARTITION_t *ptable, FILE *stream);


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

#define FS_PARTITION_INLINE

#else

#ifndef MY_INLINE
#define MY_INLINE 
#endif

#endif


#ifdef  FS_PARTITION_INLINE

/********************************************************************/    
/*                                                                  */
/*                    FS_PARTITION module                           */ 
/*                                                                  */
/********************************************************************/    

MY_INLINE
ULINT FS_PARTITION_get_no_partitions(FS_PARTITION_t *FS_partition)
{
  return FS_partition->no_partitions;
}

MY_INLINE
int FS_PARTITION_check_seq(BIOSEQ *seq, FS_PARTITION_t *FS_partition)
{
  char *c = seq->start;
  char *c1 = seq->start + seq->len;

  while(c < c1)
    {
      if (FS_partition->partition_table[*c & A_SIZE_MASK] == -1)
	return 0;
      c++;
    }
  return 1;

} 

MY_INLINE
int BIOSEQ_2_FS_SEQ(BIOSEQ *seq, FS_PARTITION_t *FS_partition,
		    FS_SEQ_t *FS_seq)
{
  register ULINT K1 = 1;
  register ULINT K = FS_partition->no_partitions;
  register char *c = seq->start;
  register char *c1 = seq->start + seq->len;
  register SLINT i;
  register FS_SEQ_t FS_seq0 = 0;

  while(c < c1)
    {
      i = FS_partition->partition_table[*c & A_SIZE_MASK];
      if (i == -1) return 0;
      FS_seq0 += K1 * i;
      K1 *= K;
      c++;
    }
  *FS_seq = FS_seq0;
  return 1;
}

#endif /* #ifdef  FS_PARTITION_INLINE */
#endif /* #ifndef _FS_PARTITION_H */
