/********************************************************************/    
/*                                                                  */
/*                     QUAD_TREE module                             */ 
/*                                                                  */
/*                                                                  */
/*  Based on the Quad Tree algorithm of David Eppstein              */
/*  (http://www.ics.uci.edu/~eppstein/), modified for the needs of  */
/*  our indexing scheme.                                            */
/*                                                                  */
/*  The half-matrix of (symmetric) distances is passed to the       */
/*  constructor (d(i,j) = half_array(i,j), j < i). Distances must   */
/*  be able to be expressed as 8 bit integers (0 - 255).            */
/*                                                                  */
/*  No distances are re-evaluated. When the clusters i,j are merged */
/*  we have for each other active k, d(i,k) = max{d(i,k), d(j,k)}.  */
/*                                                                  */
/********************************************************************/    

#ifndef _QUADTREE_H
#define _QUADTREE_H

#include "misclib.h"
#include <stdio.h>

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
  unsigned char *min_dist;
  unsigned char *cell;
} QUAD_TREE_LEVEL;

typedef struct
{
  ULINT no_levels;
  ULINT max_levels;
  QUAD_TREE_LEVEL *level;
  ULINT no_clusters;
  ULINT max_used_clusters;
  ULINT max_clusters;
  unsigned char *active;
} QUAD_TREE;

/* Main Constructor */
QUAD_TREE *QUAD_TREE_create(ULINT no_pts, unsigned char *dist);

/* Reuse Constructor */
void QUAD_TREE_reset(QUAD_TREE *qtree, ULINT no_pts, 
		     unsigned char *dist); 

/* Destructor */
void QUAD_TREE_destroy(QUAD_TREE *qtree);

/* Delete cluster */
void QUAD_TREE_delete_pt(QUAD_TREE *qtree, ULINT pt);

/* Merge two clusters */
void QUAD_TREE_merge_pts(QUAD_TREE *qtree, ULINT pt1, ULINT pt2);

/* Find pair with minimal distance */
int QUAD_TREE_find_minimal_pair(QUAD_TREE *qtree, ULINT *pt1, 
				ULINT *pt2); 

/* Half-array function */
ULINT half_array(ULINT i, ULINT j);

void print_levels(QUAD_TREE *qtree, FILE *stream);

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS                    ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

#define half_array(i, j) \
       /* Assume j < i */ \
       (((i)*((i)-1))/2 + (j))

#endif /* #ifndef _QUADTREE_H */
