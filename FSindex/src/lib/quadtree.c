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
 

#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "quadtree.h"

/* Static functions */ 
static inline
unsigned char min_dist_cell_update(QUAD_TREE *qtree, 
				   unsigned char *mcell, ULINT i,
				   ULINT j, ULINT l)
{
  unsigned char mdist;
  /* Find the minimal cell */
  mdist = qtree->level[l+1].min_dist[half_array(2*i-1, 2*j)];
  *mcell = 0;

      if (qtree->level[l+1].min_dist[half_array(2*i, 2*j)] < mdist)
	{
	  mdist = qtree->level[l+1].min_dist[half_array(2*i, 2*j)];
	  *mcell = 1;
	}
      if (qtree->level[l+1].min_dist[half_array(2*i, 2*j+1)] < mdist) 
	{
	  mdist = qtree->level[l+1].min_dist[half_array(2*i, 2*j+1)];
	  *mcell = 3;
	}
      if (i > j+1 &&
	  qtree->level[l+1].min_dist[half_array(2*i-1, 2*j+1)] < mdist) 
	{
	  mdist = qtree->level[l+1].min_dist[half_array(2*i-1, 2*j+1)];
	  *mcell = 2;
	}
      return mdist;
}



static inline
void QUAD_TREE_update_row(QUAD_TREE *qtree, ULINT row, ULINT lvl)
{
  ULINT k;
  ULINT k0;
  unsigned char mcell;

  k0 = half_array(row, 0);
  for(k=0; k < row; k++, k0++)
    {
      qtree->level[lvl].min_dist[k0] = 
	min_dist_cell_update(qtree, &mcell, row, k, lvl);
      qtree->level[lvl].cell[k0] = mcell;
    }
}

static inline
void QUAD_TREE_update_col(QUAD_TREE *qtree, ULINT col, ULINT lvl)
{
  ULINT k;
  ULINT k0;
  unsigned char mcell;
  int no_rows = (1 << lvl);

  for(k=col+1; k <= no_rows; k++)
    {
      k0 = half_array(k, col);
      qtree->level[lvl].min_dist[k0] = 
	min_dist_cell_update(qtree, &mcell, k, col, lvl);
      qtree->level[lvl].cell[k0] = mcell;
    }
}



/* Main Constructor */
QUAD_TREE *QUAD_TREE_create(ULINT no_pts, unsigned char *dist)
{
  QUAD_TREE *qtree = callocec(1, sizeof(QUAD_TREE));
  QUAD_TREE_reset(qtree, no_pts, dist);
  return qtree;
}

/* Reuse Constructor */
void QUAD_TREE_reset(QUAD_TREE *qtree, ULINT no_pts, 
		     unsigned char *dist)
{
  int i = 0; 
  int j = no_pts - 1;          /* Level counter */
  ULINT no_levels1;
  ULINT max_clusters1;
  ULINT dist_size = j*(j+1)/2; /* Size of distance half matrix */  
  ULINT w;                     /* Width of half matrices */

  qtree->no_clusters = no_pts;
  qtree->max_used_clusters = no_pts;

  /* Find smallest i such that 2^i >= no_pts - 1;
     Then no_levels (as usual including the 0th level) = i+1. */
  while (j > (1 << i))
    i++;
  no_levels1 = i+1;

  /* max_clusters is the maximum number of points that can be loaded
     into quad tree. It is equal to  2^no_levels1 + 1 because each
     level is a half-matrix without the main diagonal */
  max_clusters1 = (1 << no_levels1) + 1;

  /* Enlarge quadtree if necessary */
  if (qtree->max_levels < no_levels1)
    {
      qtree->level = reallocec(qtree->level, 
		       no_levels1 * sizeof(QUAD_TREE_LEVEL)); 
      for (j = qtree->max_levels; j < no_levels1; j++)
	{
	  w = 1 << j;
	  qtree->level[j].min_dist = mallocec(w*(w+1)/2);
	  qtree->level[j].cell = mallocec(w*(w+1)/2);
	}
      qtree->active = reallocec(qtree->active, max_clusters1);
      qtree->max_levels = no_levels1;
      qtree->max_clusters = max_clusters1;
    }
  qtree->no_levels = no_levels1;

  /* Set active flags */
  memset(qtree->active, 1, no_pts);
  memset(qtree->active+no_pts, 0, qtree->max_clusters - no_pts);
  
  /* Set top level (copy distance half-matrix passed as argument) */
  j = no_levels1 - 1;
  w = 1 << j;
  memcpy(qtree->level[j].min_dist, dist, dist_size);
  /* Note that since min_dist is unsigned, -1 == UCHAR_MAX */
  memset(qtree->level[j].min_dist+dist_size, -1, 
	 w*(w+1)/2 - dist_size); 

  /* Update all other levels */
  for (j = no_levels1 - 2; j >= 0; j--)
    {
      w = 1 << j;
      for (i=1; i <= w; i++)
	QUAD_TREE_update_row(qtree, i, j);
    }
}

/* Destructor */
void QUAD_TREE_destroy(QUAD_TREE *qtree)
{
  ULINT i;
  free(qtree->active);
  for (i = 0; i < qtree->max_levels; i++)
    {
      free(qtree->level[i].min_dist);
      free(qtree->level[i].cell);
    }
  free(qtree->level);
  free(qtree);
}

/* Delete cluster */
void QUAD_TREE_delete_pt(QUAD_TREE *qtree, ULINT pt)
{
  ULINT i;
  ULINT j;
  ULINT l = qtree->no_levels - 1;

  if (qtree->active[pt] == 0)
    return;

  /* Top level */
  qtree->active[pt] = 0;

  /* Set the pt row to the maximum distance */
  memset(qtree->level[l].min_dist + pt*(pt-1)/2, -1, pt);

  /* Set the pt col to the maximum distance */
  for(i=pt+1; i < qtree->max_used_clusters; i++)
    qtree->level[l].min_dist[half_array(i,pt)] = -1;
    
  qtree->no_clusters--;

  /* Update all other levels */
  l--;
  i = pt;
  j = pt;
  do
    {
      i = (i+1)/2;
      j = j/2;
      QUAD_TREE_update_row(qtree, i, l);
      QUAD_TREE_update_col(qtree, j, l);
      l--;
    } while(l);
}

/* Merge two clusters */
void QUAD_TREE_merge_pts(QUAD_TREE *qtree, ULINT pt1, ULINT pt2)
{
  ULINT i;
  ULINT j;
  ULINT k;
  ULINT k0;
  ULINT k1;
  int l = qtree->no_levels - 1;
  ULINT pt1_i = pt1;
  ULINT pt1_j = pt1;
  ULINT pt2_i = pt2;
  ULINT pt2_j = pt2;
  

  /* Pick the point with larger index as the one pt1 and pt2 are
     getting merged into. Traversing by column is more efficient. */
  if (pt1 < pt2)
    {
      j = pt1;
      i = pt2;
    }
  else
    {
      i = pt1;
      j = pt2;
    }


  /* Update merged row - traverse by column */
  k0 = half_array(i,0);
  for(k=0; k < i; k++, k0++)
    {
      if (j < k)
	k1 = half_array(k,j);
      else if (k < j)
	k1 = half_array(j,k);
      else /* j == k */
	continue;

      if (qtree->level[l].min_dist[k0] < qtree->level[l].min_dist[k1])  
	qtree->level[l].min_dist[k0] = qtree->level[l].min_dist[k1]; 	
    }

  /* Update merged col - traverse by row */

  for(k=i+1; k < qtree->max_used_clusters; k++)
    {
      k0 = half_array(k,i);
      if (j < k)
	k1 = half_array(k,j);
      else if (k < j)
	k1 = half_array(j,k);
      else /* j == k */
	continue;
      
      if (qtree->level[l].min_dist[k0] < qtree->level[l].min_dist[k1]) 
	qtree->level[l].min_dist[k0] = qtree->level[l].min_dist[k1]; 
    } 

  /* Delete the other point */
  qtree->active[j] = 0;

  /* Set the jth row to the maximum distance */
  memset(qtree->level[l].min_dist + j*(j-1)/2, -1, j);

  /* Set the pt col to the maximum distance */
  for(k=j+1; k < qtree->max_used_clusters; k++)
    qtree->level[l].min_dist[half_array(k,j)] = -1;
  
  qtree->no_clusters--;

  /* Update all other levels */
  for (l = qtree->no_levels - 2; l >= 0; l--)
    {
      pt1_i = (pt1_i+1)/2;
      pt1_j = pt1_j/2;
      pt2_i = (pt2_i+1)/2;
      pt2_j = pt2_j/2;

      QUAD_TREE_update_row(qtree, pt1_i, l);
      if (pt2_i != pt1_i)
	QUAD_TREE_update_row(qtree, pt2_i, l);
      QUAD_TREE_update_col(qtree, pt1_j, l);
      if (pt2_j != pt1_j)
	QUAD_TREE_update_col(qtree, pt2_j, l);
   }
}

/* Find pair with minimal distance */
int QUAD_TREE_find_minimal_pair(QUAD_TREE *qtree, ULINT *i, 
				ULINT *j)
{
  ULINT k=0;
  int cell;  
  *i = 1;
  *j = 0;

  for (;k < qtree->no_levels - 1; k++)
    {      
      cell = qtree->level[k].cell[half_array(*i, *j)];
      *i = (cell & 1) + 2*(*i) - 1;
      *j = (cell >> 1) + 2*(*j);
#if DEBUG > 10
      printf("level = %d, cell = %d, (%d,%d)\n", k,
	     cell, *i, *j);      
#endif      
    }
  
  return (int) qtree->level[0].min_dist[0];
} 


void print_levels(QUAD_TREE *qtree, FILE *stream)
{
  int i;
  int j;
  int l;

  for(l = 0; l < qtree->no_levels; l++)
    {
      fprintf(stream, "LEVEL %d\n", l);  
      fprintf(stream, "Minimum distance matrix\n"); 
      for (i=1; i <= (1 << l); i++)
	{
	  fprintf(stream, "%2d ", i);
	  for (j=0; j < i; j++)
	    fprintf(stream, "%3d ", qtree->level[l].min_dist[half_array(i, j)]);
	  fprintf(stream, "\n"); 
	}  
      fprintf(stream, "\n   "); 
      for (j=0; j < (1 << l); j++)
	fprintf(stream, "%3d ", j);
      fprintf(stream, "\n"); 
      
#if 0      
      fprintf(stream, "Minimum distance cells\n"); 
      for (i=1; i <= (1 << l); i++)
	{
	  fprintf(stream, "%2d ", i);
	  for (j=0; j < i; j++)
	    fprintf(stream, "%3d ", qtree->level[l].cell[half_array(i, j)]);
	  fprintf(stream, "\n"); 
	}  
      fprintf(stream, "\n   "); 
      for (j=0; j < (1 << l); j++)
	fprintf(stream, "%3d ", j);
      fprintf(stream, "\n");
#endif
    }   
}
