/********************************************************************/    
/*                                                                  */
/*                     RANDOM SEQUENCES module                      */ 
/*                                                                  */
/********************************************************************/    
#include <time.h>
#include "randseq.h"
#include "fastadb.h"
#include "partition.h"


#define BUF_SIZE 256 

/* Constructor */
SEQ_GENERATOR_t *SEQ_GENERATOR_create(const char *filename,
				      FS_PARTITION_t *ptable)
{
  FILE *stream;
  char buffer[BUF_SIZE];
  ULINT i;
  ULINT n; 
  char c;
  ULINT freq[A_SIZE];
  SEQ_GENERATOR_t *seq_generator;

  memset(freq, 0, A_SIZE * sizeof(ULINT));
  if ((stream = fopen(filename, "r")) == NULL)
    Throw FSexcept(FOPEN_ERR, "SEQ_GENERATOR_create():"
		   " Could not open the file %s.",
		   filename);

  seq_generator = mallocec(sizeof(SEQ_GENERATOR_t));
  srand48(time(NULL)); 

  /* Skip first lines */
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  sscanf(buffer+15, "%ld", &n); 
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);

  /* Now get frequencies */

  for (i=0; i < n; i++)
    {
      fgets(buffer, BUF_SIZE, stream);
      sscanf(buffer+15, "%c", &c);
      if (ptable->partition_table[c & A_SIZE_MASK] == -1)
	{
	  freq[c & A_SIZE_MASK] = 0;
	  seq_generator->cum_freq[i] = seq_generator->total_residues; 
	}
      else
	{
	  sscanf(buffer+23, "%ld", &freq[c & A_SIZE_MASK]);
	}
#if 0
      printf("%c %12d\n", c, freq[c & A_SIZE_MASK]);
#endif
    }

  seq_generator->total_residues = 0;
  seq_generator->no_letters = A_SIZE;
  for (i=0; i < A_SIZE; i++)
    {
      seq_generator->cum_freq[i] = 
	seq_generator->total_residues + freq[i];
      seq_generator->total_residues += freq[i];
    } 

  fclose(stream);
  return seq_generator;
}

/* Destructor */
void SEQ_GENERATOR_destroy(SEQ_GENERATOR_t *seq_generator)
{
  free(seq_generator);
}

/* Generate random sequence */
void SEQ_GENERATOR_rand_seq(SEQ_GENERATOR_t *seq_generator, 
			    BIOSEQ *seq, ULINT len, char *seq_heap)
{
  ULINT i, j;
  ULINT r;
  seq->len = len;
  

  for (i=0; i < len; i++)
    {
      r = lrand48() % seq_generator->total_residues;
      for (j=1; j < seq_generator->no_letters; j++)
	{
	  if (r < seq_generator->cum_freq[j])
	    {
	      seq_heap[i] = j + 64;
	      break;
	    }
	}
    }
  seq_heap[len] = '\0';
  seq->start = seq_heap;
}
