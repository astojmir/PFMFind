/********************************************************************/    
/*                                                                  */
/*                     CLUSTERS module                              */ 
/*                                                                  */
/********************************************************************/    

#ifndef _CLUSTERS_H
#define _CLUSTERS_H

#include "misclib.h"
#include "smatrix.h"
#include "fastadb.h"
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

typedef struct
{
  SEQUENCE_DB *s_db; /* Sequence database */
  ULINT frag_len;    /* Fragment Length */
  ULINT no_seqs;     /* Number of sequences loaded */  
  ULINT max_no_seqs; /* Number of allocated spaces for sequences */
  ULINT no_clusters; /* Current number of clusters */
  SEQ_index_t *seqs; /* Array of fragments */
  ULINT *next_seq;   /* Array of pointers to next fragment - for
			linked list */
  ULINT *ccount;     /* Sequences/ cluster counter */
  ULINT *cfirst;     /* Points to the first sequence in each cluster
		      */ 
  ULINT *clast;      /* Points to the last sequence in each cluster
		      */ 
  QUAD_TREE *qtree;  /* Quad tree to dynamically merge clusters */
  SCORE_MATRIX_t *D; /* Symmetric distance matrix on sequences */
  unsigned char *dist;
} SEQ_CLUSTERS;

/* Main Constructor */
SEQ_CLUSTERS *SEQ_CLUSTERS_create(ULINT no_seqs, SEQ_index_t *seqs,
				  SCORE_MATRIX_t *D);

/* Reuse Constructor */
void SEQ_CLUSTERS_reset(SEQ_CLUSTERS *sclusters, ULINT no_seqs,
			SEQ_index_t *seqs, SCORE_MATRIX_t *D); 
  
/* Destructor */
void SEQ_CLUSTERS_destroy(SEQ_CLUSTERS *sclusters);


/* Complete linkage threshold clustering: cluster x,y if d(x,y) < T */  
void SEQ_CLUSTERS_CLT_cluster(unsigned char T, ULINT K);

/* Complete linkage pattern clustering: 
   Clusters A, B, suppose d(A,B) is minimal.
   P(C) = regular grammar patter of sequence set C.
   Let C = A \cup B. Merge A and B if:
     diam(P(C)) < T */  

void SEQ_CLUSTERS_CLP_cluster();


#endif /* #ifndef _CLUSTERS_H */
