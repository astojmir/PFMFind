/********************************************************************/    
/*                                                                  */
/*                    FS_PARTITION module                           */ 
/*                                                                  */
/********************************************************************/    
#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>

#ifdef DEBUG
#define FS_PARTITION_INLINE
#endif

#include "partition.h"

#ifdef DEBUG
#undef FS_PARTITION_INLINE
#endif


/* Main constructor */
FS_PARTITION_t *FS_PARTITION_create(const char *alphabet, 
				    char separator)
{
  FS_PARTITION_t *ptable = callocec(1, sizeof(FS_PARTITION_t));
  USINT p=0;
  memset(&(ptable->partition_table[0]), -1, A_SIZE*sizeof(SLINT));

  while (*alphabet)
    {
      if (*alphabet == separator)
	p++;
      else
	ptable->partition_table[*alphabet & A_SIZE_MASK] = p;
      alphabet++;
    }
  ptable->no_partitions = p+1;
  if (FS_PARTITION_VERBOSE > 0)
    FS_PARTITION_print(ptable, stdout);
  return ptable;
}

/* Destructor */
void FS_PARTITION_destroy(FS_PARTITION_t *FS_partition)
{
  free(FS_partition);
}

/* Read, Write */
int FS_PARTITION_write(FS_PARTITION_t *FS_partition, FILE *stream)
{
  USINT table_size = A_SIZE;
  fwrite(&table_size, sizeof(USINT), 1, stream);
  fwrite(FS_partition, sizeof(FS_PARTITION_t), 1, stream);
  return 1;
} 

FS_PARTITION_t *FS_PARTITION_read(FILE *stream)
{
  FS_PARTITION_t *ptable = callocec(1, sizeof(FS_PARTITION_t));
  USINT table_size;
  fread(&table_size, sizeof(USINT), 1, stream);
  if (table_size != A_SIZE)
    {
      free(ptable);
      return NULL;
    }
  fread(ptable, sizeof(FS_PARTITION_t), 1, stream);
  return ptable;
}

/* Access and conversion */

char *FS_seq_print(FS_SEQ_t FS_seq, FS_PARTITION_t *FS_partition,
		   ULINT frag_len)
{
  SLINT i;
  ULINT K = FS_partition->no_partitions;
  char *seq = mallocec(frag_len+1);

  for (i=0; i < frag_len; i++)
    {
      seq[i] = '0' + (FS_seq % K);
      FS_seq /= K;
    }
  seq[i] = '\0';
  return seq;
}

void FS_PARTITION_print(FS_PARTITION_t *ptable, FILE *stream)
{
  UINT_t i;

  fprintf(stream, "\n *** PARTITION TABLE ***\n");
  for (i=0; i < A_SIZE; i++)
    {
      if (ptable->partition_table[i] == -1)
	continue;
      fprintf(stream, "%2c ", 64+i);
    }
  fprintf(stream, "\n");
  for (i=0; i < A_SIZE; i++)
    {
      if (ptable->partition_table[i] == -1)
	continue;
      fprintf(stream, "%2ld ", ptable->partition_table[i]);
    }
  fprintf(stream, "\n");
  fprintf(stream, "\n");
}


