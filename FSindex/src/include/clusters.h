/********************************************************************/    
/*                                                                  */
/*                     CLUSTERS module                              */ 
/*                                                                  */
/********************************************************************/    

#ifndef _CLUSTERS_H
#define _CLUSTERS_H

#include "misclib.h"
#include <stdio.h>
#include <stdlib.h>
#include "smatrix.h"
#include "fastadb.h"
#include "partition.h"
#include "quadtree.h"

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


/* The seq_cluster structure is used for storing clusters in a linked
   list. We copy them into array before using quadtree */ 

struct seq_cluster
{
  ULINT ccount;     /* Sequences/ cluster counter */
  ULINT cfirst;     /* Points to the first sequence in each cluster
		      */ 
  ULINT clast;      /* Points to the last sequence in each cluster */
  struct seq_cluster *next;
  struct seq_cluster *previous;
};

/* Note that we have a two-dimensional linked list. Then we copy into
   an array of linked lists. */


typedef struct
{
  SEQUENCE_DB *s_db; /* Sequence database */
  ULINT frag_len;    /* Fragment Length */
  ULINT no_seqs;     /* Number of sequences loaded */  
  ULINT max_no_seqs; /* Number of allocated spaces for sequences */
  ULINT no_clusters; /* Current number of clusters */
  ULINT max_no_clusters; /* Allocated space for clusters */
  ULINT no_clusters0;/* Number of clusters with identical sequences
			collected */
  SEQ_index_t *seqs; /* Array of fragments */
  ULINT *next_seq;   /* Array of pointers to next fragment - for
			linked list */
  ULINT *ccount;     /* Sequences/ cluster counter */
  ULINT *cfirst;     /* Points to the first sequence in each cluster
		      */ 
  ULINT *clast;      /* Points to the last sequence in each cluster
		      */
  struct seq_cluster *rfirst;
  struct seq_cluster *rseqs;
  QUAD_TREE *qtree;  /* Quad tree to dynamically merge clusters */
  SCORE_MATRIX_t *D; /* Symmetric distance matrix on sequences */
  FS_PARTITION_t *ptable; /* Partition table */
  unsigned char *dist;
  ULINT max_frag_len; /* Maximum size of the arrays below */
  int *pattern_temp; /* Used for generating patterns */
  int *offset_temp; /* Used for generating patterns */
} SEQ_CLUSTERS;

typedef struct
{
  ULINT no_seq;
  ULINT no_clusters;
  ULINT *cluster_size;
  char *cluster_pattern;
  SEQ_index_t *seqs;
} CLUSTER_BIN;

/* TO DO:

   - nice bin printing routine
   - search routine

*/
/* Main Constructor */
SEQ_CLUSTERS *SEQ_CLUSTERS_create(ULINT no_seqs, ULINT frag_len,
				 SEQ_index_t *seqs,
				 SCORE_MATRIX_t *D, 
				 SEQUENCE_DB *s_db,
				 FS_PARTITION_t *ptable);

/* Reuse Constructor */
void SEQ_CLUSTERS_reset(SEQ_CLUSTERS *sclusters, ULINT no_seqs,
			ULINT frag_len, SEQ_index_t *seqs, 
			SCORE_MATRIX_t *D, SEQUENCE_DB *s_db,
			FS_PARTITION_t *ptable);

  
/* Destructor */
void SEQ_CLUSTERS_destroy(SEQ_CLUSTERS *sclusters);


/* Complete linkage threshold clustering: cluster x,y if d(x,y) < T */  
void SEQ_CLUSTERS_CLT_cluster(SEQ_CLUSTERS *sclusters, 
			      CLUSTER_BIN *cbin, 	       
			      unsigned int T, ULINT K);

/* Complete linkage pattern clustering: 
   Clusters A, B, suppose d(A,B) is minimal.
   P(C) = regular grammar patter of sequence set C.
   Let C = A \cup B. Merge A and B if:
     diam(P(C)) < T */  

void SEQ_CLUSTERS_CLP_cluster(void);





void CLUSTER_BIN_print(CLUSTER_BIN *cbin, int frag_len, 
		       SEQUENCE_DB *s_db, FS_PARTITION_t *ptable, 
		       FILE *stream);


#endif /* #ifndef _CLUSTERS_H */
