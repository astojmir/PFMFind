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
   P_SIZE - maxumum size of pattern alphabet;
   C_SIZE - maximum number of partitions.
*/
#define A_SIZE 32
#define A_SIZE_MASK 31
#define P_SIZE 256
#define C_SIZE 8

/* We store SEQ_index list in FS_BIN */
typedef ULINT SEQ_index_t;


typedef ULINT FS_SEQ_t;

typedef struct
{
  int no_partitions;
  char *alphabet;                 /* Alphabet divided in partitions by
				     separator. e.g. ABC#DEF#GHIJK */
  char separator;                 /* Separator character */
  int partition_table[A_SIZE];    /* For each letter in the alphabet
				     (0 to A_SIZE_MASK - use &
				     A_SIZE_MASK with any character to
				     convert), stores the partition it
				     belongs to  (-1 if none) */ 
  int partition_pos[A_SIZE];      /* For each letter in the alphabet
				     (0 to A_SIZE_MASK), stores the
				     position within the partition or
				     -1 if none */ 
  int partition_size[C_SIZE];     /* For each partition, stores
				     the number of letters in it. */
  int partition_poffset[C_SIZE];  /* For each partition, stores the
				     offset to the first pattern
				     letter of that partition. */
  int partition_loffset[C_SIZE];  /* For each partition, stores the
				     offset of its first letter in the
				     alphabet string. */ 
  
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
int FS_PARTITION_get_no_partitions(FS_PARTITION_t *FS_partition);
MY_INLINE
int BIOSEQ_2_FS_SEQ(BIOSEQ *seq, FS_PARTITION_t *FS_partition,
		    FS_SEQ_t *FS_seq);
char *FS_seq_print(FS_SEQ_t FS_seq, FS_PARTITION_t *FS_partition,
		   ULINT frag_len);
MY_INLINE
int FS_PARTITION_check_seq(BIOSEQ *seq, FS_PARTITION_t *FS_partition);


int FS_PARTITION_get_pttn(FS_PARTITION_t *ptable, int letter);

int FS_PARTITION_get_pos(FS_PARTITION_t *ptable, int letter);

int FS_PARTITION_get_size(FS_PARTITION_t *ptable, int partition);

int FS_PARTITION_get_poffset(FS_PARTITION_t *ptable, int partition);

int FS_PARTITION_get_letter(FS_PARTITION_t *ptable, int partition, 
			    int pos);


/* Printing */
void FS_PARTITION_print(FS_PARTITION_t *ptable, FILE *stream);
MY_INLINE
void FS_PARTITION_pletter_2_string(FS_PARTITION_t *ptable, int p,
				   char *s);

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS                    ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    


/********************************************************************/    
/*                                                                  */
/*                    FS_PARTITION module                           */ 
/*                                                                  */
/********************************************************************/    

#define FS_PARTITION_get_no_partitions(FS_partition) \
        ((FS_partition)->no_partitions)

#define FS_PARTITION_get_pttn(ptable, letter) \
        ((ptable)->partition_table[(letter) & A_SIZE_MASK])

#define FS_PARTITION_get_pos(ptable, letter) \
        ((ptable)->partition_pos[(letter) & A_SIZE_MASK]) 

#define FS_PARTITION_get_size(ptable, partition) \
        ((ptable)->partition_size[(partition)]) 

#define FS_PARTITION_get_poffset(ptable, partition) \
        ((ptable)->partition_poffset[(partition)])

#define FS_PARTITION_get_letter(ptable, partition, pos) \
        ((int) \
        ((ptable)->alphabet[ptable->partition_loffset[(partition)] + (pos)]))

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

MY_INLINE
void FS_PARTITION_pletter_2_string(FS_PARTITION_t *ptable, int p,
				   char *s)
{
  /* Assume that s is pre-allocated to at least C_SIZE + 1*/

  int i = ptable->no_partitions - 1;
  int j;
  int k = 0;

  while (i)
    {
      if (p >= FS_PARTITION_get_poffset(ptable, i))
	{
	  p -= FS_PARTITION_get_poffset(ptable, i);
	  break;
	}
      i--;
    }

  for (j=0; j < FS_PARTITION_get_size(ptable, i); j++)
    {
      if (p & 1)
	{
	  s[k++] = (char) FS_PARTITION_get_letter(ptable, i, j);
	}
      p = p >> 1;
    }
  s[k] = '\0';
}


#endif /* #ifndef _FS_PARTITION_H */
