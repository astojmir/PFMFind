/********************************************************************/    
/*                                                                  */
/*                     CLUSTERS module                              */ 
/*                                                                  */
/********************************************************************/    

#include "misclib.h"
#include <string.h>
#include "quadtree.h"
#include "smatrix.h"
#include "fastadb.h"
#include "clusters.h"



/* Static functions */
static inline
void eval_pairwise_distances(SEQ_CLUSTERS *sclusters)
{
  ULINT i;
  ULINT j;
  ULINT k = 0;
  BIOSEQ query;
  BIOSEQ subject;

  for (i=1; i < sclusters->no_seqs; i++)
    {
      fastadb_get_Ffrag(sclusters->s_db, sclusters->frag_len, 
			&query, sclusters->seqs[i]);      
      for (j=0; j < i; j++)
	{
	  fastadb_get_Ffrag(sclusters->s_db, sclusters->frag_len, 
			    &subject, sclusters->seqs[j]);
	  sclusters->dist[k] = SCORE_MATRIX_evaluate(sclusters->D, 
				 &query, &subject);
	  k++;
	}
    }
} 

static inline
void merge_clusters(SEQ_CLUSTERS *sclusters, ULINT cl1, ULINT cl2)
{
  ULINT i;
  ULINT j;

  if (cl1 < cl2)
    {
      j = cl1;
      i = cl2;
    }
  else
    {
      i = cl1;
      j = cl2;
    }
  if (sclusters->ccount[i] == 0)
    sclusters->cfirst[i] = sclusters->cfirst[j]; 
  else
    sclusters->next_seq[sclusters->clast[i]] = sclusters->cfirst[j];
  sclusters->clast[i] = sclusters->clast[j];
  sclusters->ccount[i] += sclusters->ccount[j];
  sclusters->no_clusters --;
  sclusters->ccount[j] = 0;
  sclusters->cfirst[j] = -1;
  sclusters->clast[j] = -1;
}


static inline
ULINT collect_small_clusters(SEQ_CLUSTERS *sclusters, ULINT K)
{
  /* The largest index of a small, (ccount < K) cluster. 
   Has to be the largest index because merge_clusters() merges into
   the cluster with larger index. Note also that there will be one
   such cluster (possibly empty) for K > 1. */
  ULINT first_small = scluster->no_seqs-1;  			    
  ULINT i;

  while (sclusters.ccount[first_small] >= K && first_small >= 0)
    first_small--;

  for (i = first_small; i >=0; i--)
    {
      if (scluster->ccount[i] == 0)
	continue;
      if (scluster->ccount[i] < K)
	merge_clusters(sclusters, first_small, i);	
    }      
  return first_small;
}

static inline
void get_patterns(SEQ_CLUSTERS *sclusters)
{

}


static inline
void *pack_clusters(SEQ_CLUSTERS *sclusters)
{


}

/* Main Constructor */
SEQ_CLUSTERS *SEQ_CLUSTERS_create(ULINT no_seqs, SEQ_index_t *seqs,
				  SCORE_MATRIX_t *D, 
				  SEQUENCE_DB *s_db,  ULINT frag_len)  
{
  SEQ_CLUSTERS *sclusters = callocec(sizeof(SEQ_CLUSTERS));
  SEQ_CLUSTERS_reset(sclusters, no_seqs, seqs);
  return sclusters;
}

/* Reuse Constructor */
void SEQ_CLUSTERS_reset(SEQ_CLUSTERS *sclusters, ULINT no_seqs,
			SEQ_index_t *seqs, SCORE_MATRIX_t *D,
			SEQUENCE_DB *s_db, ULINT frag_len)
{
  int create_qtree = 0;
  sclusters->no_seqs = no_seqs;
  sclusters->no_clusters = no_seqs;
  sclusters->seqs = seqs;
  sclusters->D = D;
  sclusters->s_db = s_db;
  sclusters->frag_len = frag_len;


  if (sclusters->max_no_seqs < no_seqs)
    {
      if (sclusters->max_no_seqs == 0)
	create_qtree = 1; 
      sclusters->max_no_seqs = no_seqs;
      sclusters->next_seq = reallocec(sclusters->next_seq, no_seqs *
				      sizeof(ULINT)); 
      sclusters->ccount = reallocec(sclusters->ccount, no_seqs *
				      sizeof(ULINT)); 
      sclusters->cfirst = reallocec(sclusters->cfirst, no_seqs *
				      sizeof(ULINT)); 
      sclusters->clast = reallocec(sclusters->last, no_seqs *
				   sizeof(ULINT)); 
      sclusters->dist = reallocec(sclusters->dist, 
				  half_array(no_seqs, 0)); 
    }

  memset(sclusters->next_seq, -1, no_seqs * sizeof(ULINT));
  for (i = 0; i < no_seqs; i++)
    sclusters->ccount = 1;
  for (i = 0; i < no_seqs; i++)
    {
      sclusters->cfirst = i;    
      sclusters->clast = i;    
    }

  eval_pairwise_distances(sclusters);

  if (create_qtree)
    sclusters->qtree = QUAD_TREE_create(no_seqs, sclusters->dist);
  else
    QUAD_TREE_reset(sclusters->qtree, no_seqs, sclusters->dist);
} 

/* Destructor */
void SEQ_CLUSTERS_destroy(SEQ_CLUSTERS *sclusters)
{
  free(sclusters->next_seq);
  free(sclusters->ccount);
  free(sclusters->cfirst);
  free(sclusters->dist);
  QUAD_TREE_destroy(sclusters->qtree);
  free(sclusters);
}


void SEQ_CLUSTERS_CLT_cluster(SEQ_CLUSTERS *sclusters, 
			      unsigned char T, ULINT K)
{
  ULINT cl1;
  ULINT cl2;
  int d;
  int T0 = T;
  ULINT unclassified;

  while (sclusters->no_clusters > 1)
    {
      d = QUAD_TREE_find_minimal_pair(sclusters->qtree, &cl1, &cl2); 
      if (d > T0)
	break;
      else
	{
	  QUAD_TREE_merge_pts(sclusters->qtree, cl1, cl2);
	  merge_clusters(sclusters, cl1, cl2);
	}      
    }

  unclassified = collect_small_clusters(sclusters, K);
  get_patterns(sclusters);
  pack_clusters(sclusters);

}
