#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "fastadb.h"
#include "blastlib.h"
#include "FSindex.h"
#include "clusters.h"

#define TOO_BIG (1 << 14) 

int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;

static 
void print_help(const char *progname)
{
  fprintf(stderr, "Usage: %s options\n"
	  "\n"
	  "Mandatory:\n"
	  "-F  f   Use FSindex index file f\n"
	  "-M  f   Use score matrix file f\n"
	  "-T  n   Maximum diameter n\n"
	  "-K  n   Minimum seqs per cluster n\n"
	  "-o  f   Write results to file f\n"
	  "Optional:\n"
	  "-d  f   Write debug info to file f\n"
	  "-c  f   Write clusters (seqs) to file f\n"
	  "-b  f   Don't print progress bar\n"
	  "\n"
	  ,progname);
}





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
  
  /* Getopt variables */
  int c;                              /* Option character          */ 
  extern char *optarg;                /* Option argument           */
  extern int optind;
  extern int optopt;
  extern int opterr;  
  int errflg = 0;                     /* Option error flag         */
  
  ULINT i;
  int j;
  ULINT s;
  ULINT one_percent_bins;
  const char *outfile = NULL;
  const char *outfile1 = NULL;
  const char *debug_file = NULL;
  FILE *stream;
  FILE *stream1;
  FILE *stream2;


  int index_flag = 0;
  int matrix_flag = 0;
  int diameter_flag = 0;
  int minsize_flag = 0;
  int outfile_flag = 0;
  int debug_flag = 0;
  int clusters_flag = 0;
  int nobar_flag = 0;


  /* Input params;
     - index_file
     - matrix
     - T
     - K
     - final stats output file
     - interesting bins output file
  */

  while ((c = getopt(argc, argv, "F:M:T:K:o:d:c:b")) != EOF)
    switch (c) 
      {
      case 'F':
	filename = optarg;
	index_flag++;
	break;
      case 'M':
	matrix_full = optarg;
	matrix_flag++;
	break;
      case 'T':
	T = atoi(optarg);
	diameter_flag = 1;
	break;
      case 'K':
	K = atoi(optarg);
	minsize_flag = 1;
	break;
      case 'o':
	outfile = optarg;
	if((stream = fopen(outfile, "w")) == NULL)
	  {
	    fprintf(stderr, 
		    "Could not open output file %s!\n", outfile);
	    exit(EXIT_FAILURE);
	  }
	outfile_flag = 1;
	break;
      case 'c':
	outfile1 = optarg;
	if((stream1 = fopen(outfile1, "w")) == NULL)
	  {
	    fprintf(stderr, 
		    "Could not open output file %s!\n", outfile1);
	    exit(EXIT_FAILURE);
	  }
	clusters_flag = 1;
	break;
      case 'd':
	debug_file = optarg;
	if((stream2 = fopen(debug_file, "w")) == NULL)
	  {
	    fprintf(stderr, 
		    "Could not open output file %s!\n", debug_file);
	    exit(EXIT_FAILURE);
	  }
	debug_flag = 1;
	break;
      case 'b':
	nobar_flag = 1;
	break;
      case '?':
      default:
	errflg++;
	break;
      }

  if (!(index_flag || matrix_flag || diameter_flag || minsize_flag ||
	outfile_flag))
    errflg++;

  if (errflg)
    {
      fprintf(stderr,"Error: Insufficient or invalid arguments \n");
      print_help(argv[0]);
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
      if (!nobar_flag)
	printbar(stdout, i+1, one_percent_bins, 50);  
      /* Print debug info */
      if (debug_flag)
	{
	  rewind(stream2);
	  fprintf(stream2, "Last bin:%d, Size: %ld\n", i, no_seqs);
	}
      if (no_seqs > 0 && no_seqs < TOO_BIG)
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

	  if (clusters_flag && (cbin->no_seq > 100))
	    CLUSTER_BIN_print(cbin, frag_len, s_db, ptable, stream1);

	  free(cbin->cluster_size);
	  free(cbin->cluster_pattern);
	  free(cbin->seqs); 
	}
      else if (no_seqs == 0)
	seedhist_add_acount(overall_size, 0);
      else if (no_seqs >= TOO_BIG)
	{
	  seedhist_add_ucount(overall_size);
	  seedhist_add_acount(overall_size, no_seqs);
	}
      
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
