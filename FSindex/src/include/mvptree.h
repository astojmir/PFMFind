#ifndef _MVPTREE_H
#define _MVPTREE_H


#include <stdio.h>
#include <math.h>
#include "misclib.h"
#include "bioseq.h"
#include "smatrix.h"
#include "hit_list.h"
#include "fastadb.h"

/* DIST_TYPE must be signed because of underflow!! */
#define DIST_TYPE signed short

/* M is the number of partitions with respect to first vantage
   point in a node, N is the number of partitions with respect
   to the second vantage point (for each of the first M partitions).
   In the paper, M and N are taken to be equal and referred as the
   parameter m. K is the maximum fanout for a leaf node (referred
   to as "k" in the paper. */

//#define MVPT_M 3	
//#define MVPT_N 3	
#define MVPT_M 2	
#define MVPT_N 2	
//#define K 53
#define MVPT_K 16
#define INTERNAL 1
#define LEAF 0


struct LEAF_STRUCT {
  DIST_TYPE dist1[MVPT_K];	/* distances from the first vp */
  DIST_TYPE dist2[MVPT_K];	/* distances from the second vp */
  int datapoint[MVPT_K];		
};
struct INTERNAL_STRUCT {
  DIST_TYPE cutoff1[MVPT_M+1];	/* cutoff points for the first vp */
  DIST_TYPE cutoff2[MVPT_M][MVPT_N+1];	/* cutoff points for the second vp */
  int child[MVPT_M][MVPT_N];              /* pointers to children */
};
union nodedata {
  struct LEAF_STRUCT leaf;
  struct INTERNAL_STRUCT internal;
};
struct node_struct {
  int node_type;
  int size;
  int vp1;	                /* the oid of the first vantage point */
  int vp2;	                /* the oid of the second vantage point */
  DIST_TYPE  vp1_to_vp2;	/* the distance from the first vp to the
			           second vp. It is only used if the leaf
			           node keeps only 2 points (vp1 and vp2) */
  union nodedata CD;
};
typedef struct node_struct NODE;


#define NODES_INCR 4096
#define MAX_HITS 100

typedef void dpt_func (void *, ULINT, BIOSEQ *);
typedef int mtr_func (void *, BIOSEQ *, BIOSEQ *);

typedef struct
{
  void *db;             /* Pointer to the dataset */
  void *matrix;         /* Pointer to the score matrix */

  dpt_func *data_pt;    /* Gets the i-th point */
  mtr_func *metric;     /* Distance between two sequences */


  int num_ref;
  DIST_TYPE **refdist;
  DIST_TYPE *SPdist;    /* the array SPdist keeps the
			   distances from the query point to 
			   the first NUM_REF vantage points
			   along any search path. It is
			   referred as the array "PATH[]" in
			   the paper. */ 

  int nodes_avail;
  int no_pts;
  int nodes_used;
  NODE *nodes;

  ULINT num_calc_con;   /* Number of distance computations for
			   construction */
  ULINT seqs_visited;
  ULINT seqs_hit;
  ULINT nodes_visited;
  ULINT nodes_hit;

  int hits_avail;
  int hits_used;
  int *hits;
  DIST_TYPE *hit_dists;

} MVP_TREE;  


MVP_TREE *MVP_TREE_init(void *db, void *matrix,  dpt_func *data_pt, 
			mtr_func *metric, int num_ref, 
			DIST_TYPE **refdist);

void MVP_TREE_create(MVP_TREE *MVPT, int *oid, int oid_size);

void MVP_TREE_destroy(MVP_TREE *MVPT);

int MVP_TREE_write(MVP_TREE *MVPT, FILE *fp);
int MVP_TREE_read(MVP_TREE *MVPT, FILE *fp);




/* Printing functions */
char *MVP_TREE_sprint_stats(MVP_TREE *MVPT, int options);
void MVP_TREE_fprint_stats(MVP_TREE *MVPT, FILE *fp, int options);

/* Search functions */

void MVP_TREE_rng_srch(MVP_TREE *MVPT, BIOSEQ *query, int range);

#if 0
HIT_LIST_t *MVP_TREE_kNN_srch(MVP_TREE *MVPT, BIOSEQ *query, int kNN, 
			      HIT_LIST_t *HL);
#endif


/* Fragment metric indexes */

MVP_TREE *MVP_TREE_init2(const char *db_name, ULINT frag_len, 
			 const char *matrix_name, int num_ref); 
void MVP_TREE_del2(MVP_TREE *MVPT2);
void MVP_TREE_save2(MVP_TREE *MVPT2, const char *filename);
MVP_TREE *MVP_TREE_load2(const char *filename);

HIT_LIST_t *MVP_TREE_rng_srch2(MVP_TREE *MVPT, BIOSEQ *query, int range,
			       HIT_LIST_t *HL);



#endif /* #ifndef _MVPTREE_H */
 
