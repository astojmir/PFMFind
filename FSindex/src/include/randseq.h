/********************************************************************/    
/*                                                                  */
/*                     RANDOM SEQUENCES module                      */ 
/*                                                                  */
/********************************************************************/    

#ifndef _RANDSEQ_H
#define _RANDSEQ_H

#include "fastadb.h"
#include "partition.h"


typedef struct
{
  ULINT cum_freq[A_SIZE];
  ULINT total_residues;
  ULINT no_letters;
} SEQ_GENERATOR_t;

/* Constructor */
SEQ_GENERATOR_t *SEQ_GENERATOR_create(const char *filename,
				      FS_PARTITION_t *ptable);

/* Destructor */
void SEQ_GENERATOR_destroy(SEQ_GENERATOR_t *seq_generator);

/* Generate random sequence */
void SEQ_GENERATOR_rand_seq(SEQ_GENERATOR_t *seq_generator, 
			    BIOSEQ *seq, ULINT len, char *seq_heap);


#endif /* #ifndef _RANDSEQ_H */
