/********************************************************************/    
/*                                                                  */
/*                     CLUSTERS module                              */ 
/*                                                                  */
/********************************************************************/    

#include "misclib.h"
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "quadtree.h"
#include "smatrix.h"
#include "fastadb.h"
#include "partition.h"
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
  ULINT first_small = sclusters->no_seqs-1;  			    
  int i;

  while (sclusters->ccount[first_small] >= K && first_small >= 0)
    first_small--;

  for (i = first_small-1; i >=0; i--)
    {
      if (sclusters->ccount[i] == 0)
	continue;
      if (sclusters->ccount[i] < K)
	merge_clusters(sclusters, first_small, i);	
    }      
  return first_small;
}


static inline
void get_patterns(SEQ_CLUSTERS *sclusters, CLUSTER_BIN *cbin, 
		  ULINT unclassified)
{
  ULINT i = 0; /* Cluster counter */
  ULINT j = 0; /* Non-empty cluster counter */ 
  ULINT k = 0; /* Sequence position counter */
  ULINT s = 0; /* Sequence counter */

  int p;       /* Letter position within partition */
  ULINT q;     /* Index of the current sequence */
  BIOSEQ query;
  int *offset = sclusters->offset_temp;
  int *pattern = sclusters->pattern_temp;
  FS_PARTITION_t *ptable = sclusters->ptable;

  /* Allocate entries of cbin */
  cbin->no_seq = sclusters->no_seqs;
  cbin->no_clusters = sclusters->no_clusters - 1;
  cbin->cluster_size = mallocec(cbin->no_clusters * sizeof(ULINT));
  cbin->cluster_pattern = mallocec(cbin->no_clusters *
				   sclusters->frag_len);    
  cbin->seqs = mallocec(cbin->no_seq * sizeof(SEQ_index_t));

  
  /* Process clusters with patterns */
  for (i=0; i < sclusters->no_seqs; i++)
    {
      if (sclusters->ccount[i] == 0 || unclassified == i)
	continue;
      
      /* Set pattern to 0 */
      memset(pattern, 0, sclusters->frag_len * sizeof(int));

      cbin->cluster_size[j] = sclusters->ccount[i];

      q = sclusters->cfirst[i]; 
      while (1)
	{
	  cbin->seqs[s++] = sclusters->seqs[q];

	  fastadb_get_Ffrag(sclusters->s_db, sclusters->frag_len, 
			    &query, sclusters->seqs[q]);

	  for (k=0; k < sclusters->frag_len; k++)
	    {
	      p = FS_PARTITION_get_pos(ptable, query.start[k]);
	      pattern[k] = pattern[k] | (1 << p);
	    }

	  if (q == sclusters->clast[i])
	    break;
	  q = sclusters->next_seq[q];
	}
      
      for (k=0; k < sclusters->frag_len; k++)
	{
	  pattern[k] += offset[k];
	  cbin->cluster_pattern[j*sclusters->frag_len + k] = 
	    (char) pattern[k];
	}
 
     j++;
    }

  /* Now put unclassified sequences */

  i = unclassified;
  assert(s == cbin->no_seq - sclusters->ccount[i]);

  q = sclusters->cfirst[i]; 
  while (1)
    {
      cbin->seqs[s++] = sclusters->seqs[q];
      
      if (q == sclusters->clast[i])
	break;
      q = sclusters->next_seq[q];
    }

  assert(s == cbin->no_seq);
}

/* Main Constructor */
SEQ_CLUSTERS *SEQ_CLUSTERS_create(ULINT no_seqs, ULINT frag_len,
				  SEQ_index_t *seqs,
				  SCORE_MATRIX_t *D, 
				  SEQUENCE_DB *s_db,
				  FS_PARTITION_t *ptable)  
{
  SEQ_CLUSTERS *sclusters = callocec(1, sizeof(SEQ_CLUSTERS));
  SEQ_CLUSTERS_reset(sclusters, no_seqs, frag_len, seqs, D, s_db, ptable);  
  return sclusters;
}

/* Reuse Constructor */
void SEQ_CLUSTERS_reset(SEQ_CLUSTERS *sclusters, ULINT no_seqs,
			ULINT frag_len, SEQ_index_t *seqs, 
			SCORE_MATRIX_t *D, SEQUENCE_DB *s_db,
			FS_PARTITION_t *ptable)
{
  int create_qtree = 0;
  int pttn;
  BIOSEQ query;
  ULINT i;
  ULINT k;

  sclusters->no_seqs = no_seqs;
  sclusters->no_clusters = no_seqs;
  sclusters->seqs = seqs;
  sclusters->D = D;
  sclusters->s_db = s_db;
  sclusters->frag_len = frag_len;
  sclusters->ptable = ptable;

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
      sclusters->clast = reallocec(sclusters->clast, no_seqs *
				   sizeof(ULINT)); 
      sclusters->dist = reallocec(sclusters->dist, 
				  half_array(no_seqs, 0)); 
    }

  memset(sclusters->next_seq, -1, no_seqs * sizeof(ULINT));
  for (i = 0; i < no_seqs; i++)
    sclusters->ccount[i] = 1;
  for (i = 0; i < no_seqs; i++)
    {
      sclusters->cfirst[i] = i;    
      sclusters->clast[i] = i;    
    }

  eval_pairwise_distances(sclusters);

  if (create_qtree)
    sclusters->qtree = QUAD_TREE_create(no_seqs, sclusters->dist);
  else
    QUAD_TREE_reset(sclusters->qtree, no_seqs, sclusters->dist);

  /* Allocate pattern_temp and offset_temp */
  if (sclusters->max_frag_len < sclusters->frag_len)
    {
      sclusters->max_frag_len = sclusters->frag_len;
      sclusters->pattern_temp = mallocec(sclusters->max_frag_len *
					 sizeof(int));
      sclusters->offset_temp = mallocec(sclusters->max_frag_len *
					sizeof(int));
    }      
  /* Now get the partition offsets */
  fastadb_get_Ffrag(sclusters->s_db, sclusters->frag_len, 
		    &query, sclusters->seqs[0]);
  for (k=0; k < sclusters->frag_len; k++)
    {
      pttn = FS_PARTITION_get_pttn(ptable, query.start[k]);
      sclusters->offset_temp[k] = 
	FS_PARTITION_get_poffset(ptable, pttn);
    }  

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
			      CLUSTER_BIN *cbin, 	       
			      unsigned int T, ULINT K)
{
  ULINT cl1;
  ULINT cl2;
  int d;
  int T0 = T;
  ULINT unclassified;

#if DEBUG > 10
  int l;
  printf("INITIALISED QUAD TREE\n");
  printf("no_clusters = %d\n", sclusters->no_clusters);
  print_levels(sclusters->qtree, stdout);
#endif	  

  while (sclusters->no_clusters > 1)
    {
      d = QUAD_TREE_find_minimal_pair(sclusters->qtree, &cl1, &cl2); 
      if (d > T0)
	break;
      else
	{
	  QUAD_TREE_merge_pts(sclusters->qtree, cl1, cl2);
	  merge_clusters(sclusters, cl1, cl2);
#if DEBUG > 10
	  printf("MERGE CLUSTERS (%d,%d)\n", cl1, cl2);
	  printf("no_clusters = %d\n", sclusters->no_clusters);
	  print_levels(sclusters->qtree, stdout);
#endif	  
	}      
    }

  unclassified = collect_small_clusters(sclusters, K);
  get_patterns(sclusters, cbin, unclassified);
}


void CLUSTER_BIN_print(CLUSTER_BIN *cbin, int frag_len, 
		       SEQUENCE_DB *s_db, FS_PARTITION_t *ptable, 
		       FILE *stream)
{
  ULINT i;
  ULINT j;
  ULINT k;
  char *pattern;
  char s[C_SIZE+1];
  ULINT seq=0;
  BIOSEQ query;

  fprintf(stream, "Total sequences: %ld\n", cbin->no_seq);
  fprintf(stream, "Total clusters: %ld\n", cbin->no_clusters);
  fprintf(stream, "\n");

  for (i=0; i < cbin->no_clusters; i++)
    {
      pattern = cbin->cluster_pattern + i*frag_len;
      for (j=0; j < frag_len; j++)
	{
	  FS_PARTITION_pletter_2_string(ptable, *(pattern+j), &s[0]);
	  if (strlen(s) > 1)    
	    fprintf(stream,"[");
	  fprintf(stream,"%s",s);
	  if (strlen(s) > 1)    
	    fprintf(stream,"]");
	  if (j < frag_len-1)
	    fprintf(stream,"-");
	}
      fprintf(stream, "\n");
      fprintf(stream, "\n");

      for (k=0; k < cbin->cluster_size[i]; k++)
	{
      	  fastadb_get_Ffrag(s_db, frag_len, &query,
			    cbin->seqs[seq++]);
	  fprintf(stream, "%*.*s\n", frag_len, frag_len, query.start);
	}

      fprintf(stream, "\n");

    }

  fprintf(stream, "Unclassified sequences\n");
  while (seq < cbin->no_seq)
    {
      fastadb_get_Ffrag(s_db, frag_len, &query,
			cbin->seqs[seq++]);
      fprintf(stream, "%*.*s\n", frag_len, frag_len, query.start);
    }
  fprintf(stream, "-------------------------------------------------\n");
}
