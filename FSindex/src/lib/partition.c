/********************************************************************/    
/*                                                                  */
/*                    FS_PARTITION module                           */ 
/*                                                                  */
/********************************************************************/    
#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef MY_INLINE
#define MY_INLINE
#include "partition.h"
#undef MY_INLINE
#endif

/* Main constructor */
FS_PARTITION_t *FS_PARTITION_create(const char *alphabet, 
				    char separator)
{
  FS_PARTITION_t *ptable = mallocec(sizeof(FS_PARTITION_t));
  int partition = 0;
  int pos = 0;
  int i = 0;
  int poffset = 0;

  /* Make sure that the first letter is not separator */
  if (*alphabet == separator)
    {
      alphabet++;
      i++;
    }

  memset(&(ptable->partition_table[0]), -1, A_SIZE*sizeof(int));
  memset(&(ptable->partition_pos[0]), -1, A_SIZE*sizeof(int));
  memset(&(ptable->partition_size[0]), -1, C_SIZE*sizeof(int));
  memset(&(ptable->partition_poffset[0]), -1, C_SIZE*sizeof(int));
  memset(&(ptable->partition_loffset[0]), -1, C_SIZE*sizeof(int));

  ptable->alphabet = strdup(alphabet);
  ptable->separator = separator;

  ptable->partition_poffset[partition] = 0;
  ptable->partition_loffset[partition] = i;

  while (*alphabet)
    {
      if (*alphabet == separator)
	{
	  ptable->partition_size[partition] = pos;

	  partition++;
	  if (partition >= C_SIZE)
	    {
	      fprintf(stderr, "Too many partitions: maximum %d!\n",
		      C_SIZE);
	      FS_PARTITION_destroy(ptable);      
	      return NULL;
	    }
	  poffset += (1 << pos); 
	  if (poffset > P_SIZE)
	    {
	      fprintf(stderr, "Too many pattern letters: maximum %d!\n",
		      P_SIZE);
	      FS_PARTITION_destroy(ptable);      
	      return NULL;
	    }
	  ptable->partition_poffset[partition] = poffset;
	  ptable->partition_loffset[partition] = i+1;
	  pos = 0;
	}
      else
	{
	  ptable->partition_table[*alphabet & A_SIZE_MASK] =partition;   
	  ptable->partition_pos[*alphabet & A_SIZE_MASK] = pos;
	  pos++;
	}
      alphabet++;
      i++;
    }

  ptable->partition_size[partition] = pos;
  ptable->no_partitions = partition+1;
  if (FS_PARTITION_VERBOSE > 0)
    FS_PARTITION_print(ptable, stdout);
  return ptable;
}

/* Destructor */
void FS_PARTITION_destroy(FS_PARTITION_t *FS_partition)
{
  free(FS_partition->alphabet);
  free(FS_partition);
}

/* Read, Write */
int FS_PARTITION_write(FS_PARTITION_t *ptable, FILE *stream)
{
  int a_size = A_SIZE;
  int c_size = C_SIZE;
  int p_size = P_SIZE;
  int alphabet_size = strlen(ptable->alphabet)+1;

  fwrite(&a_size, sizeof(int), 1, stream);
  fwrite(&c_size, sizeof(int), 1, stream);
  fwrite(&p_size, sizeof(int), 1, stream);
  fwrite(&ptable->no_partitions, sizeof(int), 1, stream);
  fwrite(&alphabet_size, sizeof(int), 1, stream);
  fwrite(ptable->alphabet, sizeof(char), alphabet_size, stream);
  fwrite(&ptable->separator, sizeof(char), 1, stream);
  fwrite(ptable->partition_table, sizeof(int), a_size, stream);
  fwrite(ptable->partition_pos, sizeof(int), a_size, stream);
  fwrite(ptable->partition_size, sizeof(int), c_size, stream);
  fwrite(ptable->partition_poffset, sizeof(int), c_size, stream);
  fwrite(ptable->partition_loffset, sizeof(int), c_size, stream);
  return 1;
} 

FS_PARTITION_t *FS_PARTITION_read(FILE *stream)
{
  FS_PARTITION_t *ptable = mallocec(sizeof(FS_PARTITION_t));

  int a_size;
  int c_size;
  int p_size;
  int alphabet_size;

  fread(&a_size, sizeof(int), 1, stream);
  fread(&c_size, sizeof(int), 1, stream);
  fread(&p_size, sizeof(int), 1, stream);

  if (a_size != A_SIZE || c_size != C_SIZE || p_size != P_SIZE)
    {
      fprintf(stderr, "Cannot load partition - different A_SIZE,"
	      " C_SIZE or P_SIZE!\n");
      free(ptable);      
      return NULL;
    }

  fread(&ptable->no_partitions, sizeof(int), 1, stream);
  fread(&alphabet_size, sizeof(int), 1, stream);
  ptable->alphabet  = mallocec(alphabet_size);
  fread(ptable->alphabet, sizeof(char), alphabet_size, stream);
  fread(&ptable->separator, sizeof(char), 1, stream);
  fread(ptable->partition_table, sizeof(int), a_size, stream);
  fread(ptable->partition_pos, sizeof(int), a_size, stream);
  fread(ptable->partition_size, sizeof(int), c_size, stream);
  fread(ptable->partition_poffset, sizeof(int), c_size, stream);
  fread(ptable->partition_loffset, sizeof(int), c_size, stream);

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
      fprintf(stream, "%2d ", ptable->partition_table[i]);
    }
  fprintf(stream, "\n");
  fprintf(stream, "\n");
}


