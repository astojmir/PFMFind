#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include "fastadb.h"
#include "blastlib.h"
#include "FSindex.h"
#include "clusters.h"


int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;

int main(int argc, char **argv)
{
  /* Index data */
  const char *filename;
  FS_HASH_TABLE_t *HT;
  ULINT no_bins;
  ULINT no_seqs; 
  ULINT frag_len;
  SEQ_index_t *seqs;
  SEQUENCE_DB *s_db;
  
  /* Matrices */
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  char *matrix_full;
  char *matrix_base;
  char *matrix_dir;

  /* Cluster data */
  SEQ_CLUSTERS *sclusters;
  CLUSTER_BIN *cbin = callocec(1, sizeof(CLUSTER_BIN));

  /* Clustering params */
  unsigned int T;
  ULINT K;

  /* Stats */
  SEED_HIST *clusters_per_bin = seedhist_init(0, 10);
  SEED_HIST *cluster_size = seedhist_init(0, 500);
  SEED_HIST *unclassified_size = seedhist_init(0, 5000);
  SEED_HIST *overall_size = seedhist_init(0, 5000);
  
  
  

  ULINT i;
  int j;
  ULINT s;
  ULINT one_percent_bins;
  const char *outfile;
  const char *outfile1;
  FILE *stream;
  FILE *stream1;

  /* Input params;
     - index_file
     - matrix
     - T
     - K
     - final stats output file
     - interesting bins output file
  */

  if (argc < 7)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix_file T K outfile"
	      " outfile1 \n", argv[0]);  
      exit(EXIT_FAILURE);
    }

  filename = argv[1];
  matrix_full = argv[2];
  T = atoi(argv[3]);
  K = atoi(argv[4]);
  outfile = argv[5];
  outfile1 = argv[6];
   
  if((stream = fopen(outfile, "w")) == NULL)
    {
      fprintf(stderr, 
	      "Could not open output file %s!\n", outfile);
      exit(EXIT_FAILURE);
    }
  if((stream1 = fopen(outfile1, "w")) == NULL)
    {
      fprintf(stderr, 
	      "Could not open output file %s!\n", outfile1);
      exit(EXIT_FAILURE);
    }


  fprintf(stdout, "Loading index... \n");
  FS_INDEX_load(filename);
  ptable = FS_INDEX_get_ptable();
  s_db = FS_INDEX_get_database();

  fprintf(stdout, "Loading matrices... \n");
  split_base_dir(matrix_full, &matrix_base, &matrix_dir);  
  S = SCORE_MATRIX_create(matrix_full, ptable); 
  D = SCORE_MATRIX_S_2_Dmax(S);

  fprintf(stdout, "Clustering... \n");
  i=0;

  HT = FS_INDEX_get_hash_table();
  frag_len = FS_INDEX_get_frag_len();
  no_bins = FS_HASH_TABLE_get_no_bins(HT);
  one_percent_bins = no_bins / 100;
  no_seqs = FS_HASH_TABLE_get_no_seqs(HT, i); 
  seqs = FS_HASH_TABLE_get_all_seqs(HT, i);

  sclusters = SEQ_CLUSTERS_create(no_seqs, frag_len, seqs,
				  D, s_db, ptable);
  while (1)
    {
      /* Print bar */
      printbar(stdout, i+1, one_percent_bins, 50);  

      if (no_seqs > 0)
	{
	  SEQ_CLUSTERS_CLT_cluster(sclusters, cbin, T, K); 	       
	  seedhist_add_acount(clusters_per_bin, cbin->no_clusters);
	  s=0;
	  for (j=0; j < cbin->no_clusters; j++)
	    {
	      s += cbin->cluster_size[j];
	      seedhist_add_acount(cluster_size,
				  cbin->cluster_size[j]);
	      seedhist_add_acount(overall_size,
				  cbin->cluster_size[j]);
	    }
	  seedhist_add_acount(unclassified_size, cbin->no_seq - s);
	  seedhist_add_acount(overall_size, cbin->no_seq - s);

	  /* If the bin is interesting, print it */
	  if (cbin->no_seq > 100)
	    CLUSTER_BIN_print(cbin, frag_len, s_db, ptable, stream1);

	  free(cbin->cluster_size);
	  free(cbin->cluster_pattern);
	  free(cbin->seqs); 
	}
      else
	seedhist_add_acount(overall_size, 0);

      if (++i >= no_bins)
	break;
      no_seqs = FS_HASH_TABLE_get_no_seqs(HT, i); 
      if (no_seqs > 0)
	{
	  seqs = FS_HASH_TABLE_get_all_seqs(HT, i);
	  SEQ_CLUSTERS_reset(sclusters, no_seqs, frag_len, seqs, 
			     D, s_db, ptable);
	}      
    }

  /* Print FSindex stats */
  FS_INDEX_print_stats(stream, 0, 0); 

  /* Print params */
  fprintf(stream, "\n");
  fprintf(stream, "Distance matrix: %s\n", matrix_base);
  fprintf(stream, "Maximum diameter: %d\n", T);
  fprintf(stream, "Minumum cluster size: %ld\n", K);
  fprintf(stream, "\n");
   
  /* Print histograms */
  fprintf(stream, "****** DISTRIBUTION OF CLUSTERS PER BIN ******\n");
  seedhist_print(clusters_per_bin, stream);
  fprintf(stream, "****** DISTRIBUTION OF CLUSTER SIZES ******\n");
  seedhist_print(cluster_size, stream);
  fprintf(stream, "****** DISTRIBUTION OF UNCLASSIFIED SIZES ******\n");
  seedhist_print(unclassified_size, stream);
  fprintf(stream, "****** DISTRIBUTION OF OVERALL SIZES ******\n");
  seedhist_print(overall_size, stream);


 return EXIT_SUCCESS;
}
