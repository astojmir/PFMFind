#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "hit_list.h"
#include "partition.h"
#include "smatrix.h"
#include "FSindex.h"

#ifdef TRY_CLUSTERS
#include "clusters.h"
#endif


/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               LIBRARY FUNCTIONS                              ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    



/********************************************************************/    
/*                                                                  */
/*                     FS_HASH_TABLE module                         */ 
/*                                                                  */
/********************************************************************/    


/* Main constructor */
FS_HASH_TABLE_t *FS_HASH_TABLE_create(ULINT no_bins, ULINT def_size) 
{
  FS_HASH_TABLE_t *FS_HT = mallocec(sizeof(FS_HASH_TABLE_t));
  ULINT i;

  FS_HT->no_bins = no_bins;
  FS_HT->no_seqs = 0;
  FS_HT->max_seqs_per_bin = 5 * def_size;
  FS_HT->largest_bin = 0;
  FS_HT->heap = NULL;
  
  FS_HT->freq_seqs_per_bin = callocec(FS_HT->max_seqs_per_bin,
				      sizeof(ULINT));
  FS_HT->freq_seqs_per_bin[0] = no_bins;
  FS_HT->bin_size = mallocec(no_bins * sizeof(ULINT));
  FS_HT->max_bin_size = mallocec(no_bins * sizeof(ULINT)); 
  FS_HT->bin = mallocec(no_bins * sizeof(SEQ_index_t *));
  for (i=0; i < no_bins; i++)
    {
      FS_HT->bin_size[i] = 0;      
      FS_HT->max_bin_size[i] = def_size;
      FS_HT->bin[i] = mallocec(def_size * sizeof(SEQ_index_t));
    }
 return FS_HT;  
}

/* Destructor */
void FS_HASH_TABLE_destroy(FS_HASH_TABLE_t *FS_HT)
{
  ULINT i;
  if (FS_HT->heap == NULL)
    for (i=0; i < FS_HT->no_bins; i++)
      free(FS_HT->bin[i]);
  else
    free(FS_HT->heap);
  free(FS_HT->bin);
  free(FS_HT->bin_size);
  free(FS_HT->freq_seqs_per_bin);
  if (FS_HT->max_bin_size != NULL)
    free(FS_HT->max_bin_size);
  free(FS_HT);
}

/* Insertion */
void FS_HASH_TABLE_insert_seq(FS_HASH_TABLE_t *FS_HT, ULINT i,
			      SEQ_index_t seq_i)
{
  ULINT old_seqs_per_bin;
  FS_HT->freq_seqs_per_bin[FS_HT->bin_size[seq_i]]--;

  FS_HT->bin[seq_i][FS_HT->bin_size[seq_i]] = i;
  FS_HT->bin_size[seq_i]++;  
  if (FS_HT->bin_size[seq_i] >=  FS_HT->max_bin_size[seq_i])
    {
      FS_HT->max_bin_size[seq_i] *= 2;
      FS_HT->bin[seq_i] = reallocec(FS_HT->bin[seq_i], 
			FS_HT->max_bin_size[seq_i] *
				    sizeof(SEQ_index_t));   
    }
  FS_HT->no_seqs++;

  if (FS_HT->bin_size[seq_i] >= FS_HT->max_seqs_per_bin)
    {
      old_seqs_per_bin = FS_HT->max_seqs_per_bin;
      FS_HT->max_seqs_per_bin = FS_HT->bin_size[seq_i] + 1;       
      FS_HT->freq_seqs_per_bin = reallocec(FS_HT->freq_seqs_per_bin,
				   FS_HT->max_seqs_per_bin *
					   sizeof(ULINT)); 
      memset(FS_HT->freq_seqs_per_bin + old_seqs_per_bin, 0, 
	     (FS_HT->max_seqs_per_bin - old_seqs_per_bin) 
	     * sizeof(ULINT));
      FS_HT->largest_bin = seq_i;
    }
  FS_HT->freq_seqs_per_bin[FS_HT->bin_size[seq_i]]++;
}

/* Trimming */
void FS_HASH_TABLE_resize(FS_HASH_TABLE_t *FS_HT)
{
  ULINT i;
  for (i=0; i < FS_HT->no_bins; i++)
    {
      if (FS_HT->bin_size[i] > 0)
	{
	  FS_HT->bin[i] = reallocec(FS_HT->bin[i], 
	  		FS_HT->bin_size[i] * sizeof(SEQ_index_t)); 
	  FS_HT->max_bin_size[i] = FS_HT->bin_size[i];
	}      
    }
}


/* Reading, writing to file */
int FS_HASH_TABLE_write(FS_HASH_TABLE_t *FS_HT, FILE *stream)
{
  ULINT i;
#if 0
  ULINT one_percent = FS_HT->no_bins/100;
#endif
  fwrite(&(FS_HT->no_bins), sizeof(ULINT), 1, stream);
  fwrite(&(FS_HT->no_seqs), sizeof(ULINT), 1, stream);
  fwrite(&(FS_HT->max_seqs_per_bin), sizeof(ULINT), 1, stream);
  fwrite(&(FS_HT->largest_bin), sizeof(SEQ_index_t), 1, stream);

  fwrite(FS_HT->freq_seqs_per_bin, sizeof(ULINT),
	 FS_HT->max_seqs_per_bin, stream); 
  fwrite(FS_HT->bin_size, sizeof(ULINT), FS_HT->no_bins, stream);  

  for (i=0; i < FS_HT->no_bins; i++)
    {
      fwrite(FS_HT->bin[i], sizeof(SEQ_index_t), FS_HT->bin_size[i],
	     stream);
#if 0
      printbar(stdout, i+1, one_percent, 50);  
#endif
    }   
  return 1;
}

FS_HASH_TABLE_t *FS_HASH_TABLE_read(FILE *stream)
{
  FS_HASH_TABLE_t *FS_HT = mallocec(sizeof(FS_HASH_TABLE_t));
  ULINT i;
  SEQ_index_t *current;

  fread(&(FS_HT->no_bins), sizeof(ULINT), 1, stream);
  fread(&(FS_HT->no_seqs), sizeof(ULINT), 1, stream);
  fread(&(FS_HT->max_seqs_per_bin), sizeof(ULINT), 1, stream);
  fread(&(FS_HT->largest_bin), sizeof(SEQ_index_t), 1, stream);

  FS_HT->freq_seqs_per_bin = mallocec(FS_HT->max_seqs_per_bin
				      * sizeof(ULINT));
  FS_HT->bin_size = mallocec(FS_HT->no_bins * sizeof(ULINT));
  FS_HT->max_bin_size = NULL;
  FS_HT->bin = mallocec(FS_HT->no_bins * sizeof(SEQ_index_t *));

  fread(FS_HT->freq_seqs_per_bin, sizeof(ULINT),
	 FS_HT->max_seqs_per_bin, stream); 
  fread(FS_HT->bin_size, sizeof(ULINT), FS_HT->no_bins, stream);  

  FS_HT->heap = mallocec(FS_HT->no_seqs * sizeof(SEQ_index_t)); 
  current =  FS_HT->heap;
  for (i=0; i < FS_HT->no_bins; i++)
    {
#if 0
      if (FS_HT->bin_size[i] > 0)
	{
	  FS_HT->bin[i] = mallocec(FS_HT->bin_size[i] *
				   sizeof(SEQ_index_t));
	} 
#endif
      FS_HT->bin[i] = current;
      fread(FS_HT->bin[i], sizeof(SEQ_index_t), FS_HT->bin_size[i],
	    stream);
      current += FS_HT->bin_size[i];
    }   
  return FS_HT;  
}

/* Statistics */
void FS_HASH_TABLE_print_stats(FS_HASH_TABLE_t *FS_HT, FILE *stream,
			       FS_PARTITION_t *FS_partition, 
			       ULINT frag_len)
{
  ULINT i;
  
  fprintf(stream, "Total number of fragments in index: %ld\n", 
	  FS_HT->no_seqs);
  fprintf(stream, "Total number of index entries: %ld\n", 
	  FS_HT->no_bins);
  fprintf(stream, "Average size of index entry: %.2f\n", 
	  (float) FS_HT->no_seqs / (float) FS_HT->no_bins);
  fprintf(stream, "Maximum size of index entry: %ld\n", 
	  FS_HT->max_seqs_per_bin-1);
  fprintf(stream, "Largest index entry: %s\n",
	  FS_seq_print(FS_HT->largest_bin, FS_partition, frag_len));

  fprintf(stream, "\n* Distribution of numbers of fragments per"
	  " index entry *\n");
  fprintf(stream, "%17s %10s\n","Index entry size", "Frequency");
  for (i = 0; i < FS_HT->max_seqs_per_bin; i++)
    {
      if (FS_HT->freq_seqs_per_bin[i] > 0)
	fprintf(stream, "%17ld %10ld\n", i, 
		FS_HT->freq_seqs_per_bin[i]);
    }
  fprintf(stream, "\n");
}

/********************************************************************/ 
/*                                                                  */
/*                     FS_INDEX module                              */ 
/*                                                                  */
/********************************************************************/ 

/********************************************************************/ 
/*      Static variables and functions                              */
/********************************************************************/ 

/* FSindex variables - essentially constant throughout all searches
   with the same index, Can be a problem if running searches with
   other indexes in parallel with shared memory, but this is very
   unlikely. */ 

/* These were members of FS_INDEX_t structures before but it seems
   better this way as we will probably not do searhces with
   different databases at the same time. 
*/  

static char *db_name;
static char *matrix;
static char *alphabet; 
static char *index_name;
static char separator;
static int frag_len; 
static ULINT db_no_frags;
static SEQUENCE_DB *s_db;
static FS_PARTITION_t *ptable;
static FS_HASH_TABLE_t *HT;
static int K; /* Number of letters in FS alphabet */
static ULINT *FS_KK; /* FS_KK[i] = FS_K ^ i */
static int *K_modtable; /* K_modtable[i] = i % K for i < 2*K */

/* Scoring matrices and thresholds - search instance specific */
static SCORE_MATRIX_t *S;
static SCORE_MATRIX_t *D;
static int S_cutoff;
static int D_cutoff;
static FS_INDEX_process_func *pfunc; 

static POS_MATRIX *PS;
static POS_MATRIX *PD;
static FS_INDEX_profile_process_func *ppfunc; 
ULINT Pfrom;
ULINT Pto;

/* Query specific variables - can be shared for the same query */ 
static HIT_LIST_t *hit_list;
static BIOSEQ *query;
static int *FS_Q; /* partition numbers */
static ULINT *FS_DC;
static ULINT *FS_REM;
static int *qind; /* sequence converted to indices */



/* Main constructor */
void FS_INDEX_create(const char *database, ULINT flen,
		     const char *abet, 
		     const char sepchar, int skip)
 
{
  time_t time0 = time(NULL);
  time_t time1;
  double dt;
  ULINT one_percent_fragments;
  char *db_dir;
  char *db_base;

  ULINT i, j;
  ULINT no_FS;
  ULINT no_frags;
  BIOSEQ *frag = mallocec(sizeof(BIOSEQ));
  FS_SEQ_t FS_seq;
  ULINT added_frags = 0;

  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  ULINT bin_size;

  /* Initialise FS_index */
  split_base_dir(database, &db_base, &db_dir);
  db_name = strdup(db_base);
  alphabet = strdup(abet);
  separator = sepchar;
  frag_len = flen;

  FS_KK = mallocec(frag_len * sizeof(ULINT));
  FS_Q = mallocec(frag_len * sizeof(int));
  FS_DC = mallocec(frag_len * sizeof(ULINT));
  FS_REM = mallocec(frag_len * sizeof(ULINT));
  qind = mallocec(frag_len * sizeof(int));
  
  /* Load database */
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = NONE;
  fastadb_argv[0].access_type = MEMORY;
  fastadb_argv[1].retrieve_deflines = NO;

  s_db = fastadb_open(database, fastadb_argt, fastadb_argv); 
 
  /* Make fragment database */
  no_frags = fastadb_count_Ffrags(s_db, frag_len);
  one_percent_fragments = (ULINT) (((double) no_frags/skip) / 100);
  fastadb_init_Ffrags(s_db, frag_len);

  /* Create partition table */
  ptable = FS_PARTITION_create(alphabet, separator); 

  /* Error if the sequence is too long for the partitioning*/
  K = FS_PARTITION_get_no_partitions(ptable);
  K_modtable = mallocec(2*K*sizeof(int));
  for (i = 0; i < K; i++)
    {
      K_modtable[i] = i;
      K_modtable[i+K] = i;
    }
  if ((K-1) >> ((sizeof(ULINT)*8)/frag_len))
    {
      fprintf(stderr, 
	      "FS_INDEX_create(): Fragment length too large for"
	      " the given\n partitions and 32 bit words!\n"
	      "Partitions: %s\n", alphabet);
      /* This is a fatal error so just exit. */
      exit(EXIT_FAILURE);
    }

  /* Number of bins to allocate in hash table */
  no_FS = 1;
  for (i=0; i < frag_len; i++)
    {
      FS_KK[i] = no_FS;
      no_FS *= K;
    }

  /* Make initial bin size slightly bigger than average */
  bin_size = ((5 * no_frags) / (4 * no_FS)) + 1;

  /* Initialise hash table */
  HT = FS_HASH_TABLE_create(no_FS, bin_size);

  /* For each fragment */ 
  j = 0;
  while (fastadb_get_next_Ffrag(s_db, frag_len, frag, &i, skip))  
    {
      j++;
      /* Calculate its FS_seq */
      if (BIOSEQ_2_FS_SEQ(frag, ptable, &FS_seq))
	{
	  /* Add to appropriate FS_bin */
	  added_frags++;
	  FS_HASH_TABLE_insert_seq(HT, i, FS_seq);
	}

      /* Print progress bar */
      if (FS_INDEX_PRINT_BAR > 0)
	printbar(stdout, j+1, one_percent_fragments, 50);  
    }  

  /* Resize bins */
  FS_HASH_TABLE_resize(HT); 

  /* Take time */
  time1 = time(NULL);
  dt = difftime(time1, time0)/60.0;
  db_no_frags = no_frags;

  /* Print statistics */
  if (FS_INDEX_VERBOSE > 0)
    FS_INDEX_print_stats(stdout, j, dt); 

  /* Clean up */
  free(frag);
  free(db_dir);
  free(db_base);

  /* We have an index ready for search */
  index_name = strdup("New Index");
  return; 
}


/* Destructor */
void FS_INDEX_destroy(void)
{
  /* Free all bins */
  if (HT != NULL)
    FS_HASH_TABLE_destroy(HT); 

  /* Free partition table */
  FS_PARTITION_destroy(ptable);  

  /* Free sequence database */
  fastadb_close(s_db);

  /* Free names */
  free(db_name);
  free(alphabet);
}

/* Load, Save */
int FS_INDEX_save(const char *filename)
{ 
  FILE *stream = fopen(filename, "wb");
  int len;
 
  /* Open file */
  if(stream == NULL)
    {
      fprintf(stderr, 
	      "FS_INDEX_save(): Could not open the file %s\n",
	      filename);
      /* This is a for now a fatal error as well */
      exit(EXIT_FAILURE);
    }

  len = strlen(db_name) + 1;
  fwrite(&len, sizeof(int), 1, stream);
  fwrite(db_name, sizeof(char), len, stream);

  len = strlen(alphabet) + 1;
  fwrite(&len, sizeof(int), 1, stream);
  fwrite(alphabet, sizeof(char), len, stream);

  fwrite(&(separator), sizeof(char), 1, stream);
  fwrite(&(frag_len), sizeof(int), 1, stream);
  fwrite(&( db_no_frags), sizeof(ULINT), 1, stream);

  FS_PARTITION_write(ptable, stream); 
  FS_HASH_TABLE_write(HT, stream);

  fclose(stream);
  return 1;
}

int FS_INDEX_load(const char *filename)
{
  FILE *stream = fopen(filename, "rb");
  int len;
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  char *basename;
  char *dirname;
  char *full_dbname;
  int i;
 
  /* Open file */
  if(stream == NULL)
    {
      fprintf(stderr, 
	      "FS_INDEX_load(): Could not open index file %s!\n", 
	      filename);
      exit(EXIT_FAILURE);
    }
  index_name = strdup(filename);

  fread(&len, sizeof(int), 1, stream);
  db_name = mallocec(len);
  fread(db_name, sizeof(char), len, stream);

  split_base_dir(filename, &basename, &dirname);
  cat_base_dir(&full_dbname, db_name, dirname);

  fread(&len, sizeof(int), 1, stream);
  alphabet = mallocec(len);
  fread(alphabet, sizeof(char), len, stream);

  fread(&separator, sizeof(char), 1, stream);
  fread(&frag_len, sizeof(int), 1, stream);
  fread(&db_no_frags, sizeof(ULINT), 1, stream);

  /* Load database */
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = NONE;
  fastadb_argv[0].access_type = MEMORY;
  fastadb_argv[1].retrieve_deflines = YES;

  s_db = fastadb_open(full_dbname, fastadb_argt, fastadb_argv); 
  if (s_db == NULL)
    {
      fprintf(stderr, 
	      "FS_INDEX_load(): Could not open fasta database %s!\n", 
	      full_dbname);
      exit(EXIT_FAILURE);
    }
      
  ptable = FS_PARTITION_read(stream); 
  HT = FS_HASH_TABLE_read(stream);

  /* Allocate arrays and calculate FS_KK */
  FS_KK = mallocec(frag_len * sizeof(ULINT));
  FS_Q = mallocec(frag_len * sizeof(int));
  FS_DC = mallocec(frag_len * sizeof(ULINT));
  FS_REM = mallocec(frag_len * sizeof(ULINT));
  qind = mallocec(frag_len * sizeof(int));

  K = FS_PARTITION_get_no_partitions(ptable);
  K_modtable = mallocec(2*K*sizeof(int));
  for (i = 0; i < K; i++)
    {
      K_modtable[i] = i;
      K_modtable[i+K] = i;
    }
  FS_KK[0] = 1;
  for (i=1; i < frag_len; i++)
    FS_KK[i] = K * FS_KK[i-1];
    
  fclose(stream);
  free(basename);
  free(dirname);
  return 1;
}

/* Print */
void FS_INDEX_print_stats(FILE *stream, ULINT count, double dtime) 
{
  fprintf(stream, "\n*** FS_INDEX Statistics ***\n");
  fprintf(stream, "Database: %s\n", db_name);
  fprintf(stream, "Database size: %ld\n", s_db->length); 
  fprintf(stream, "Number of sequences: %ld\n\n", s_db->no_seq);  
  fprintf(stream, "Partitions: %s\n", alphabet);
  FS_PARTITION_print(ptable, stream);
  fprintf(stream, "Fragment Length: %d\n", frag_len);
  fprintf(stream, "Total Fragments in database: %ld\n", db_no_frags); 
  if (count > 0)
    fprintf(stream, "Total Counted Fragments in database: %ld\n", 
	    count);
  if (dtime > 0)
    fprintf(stream, "Creation Time: %.2f mins\n\n", dtime);
  FS_HASH_TABLE_print_stats(HT, stream, ptable, frag_len); 
}

/********************************************************************/ 
/*                 Access functions                                 */
/********************************************************************/ 

SEQUENCE_DB *FS_INDEX_get_database(void)
{
  return s_db;
}

FS_PARTITION_t *FS_INDEX_get_ptable(void)
{
  return ptable;
}

int FS_INDEX_get_frag_len(void)
{
  return frag_len;
}

FS_HASH_TABLE_t *FS_INDEX_get_hash_table(void)
{
  return HT;
}


/* Search */


/********************************************************************/ 
/*      Recursive tree traversal function                           */
/********************************************************************/ 
static int R;
static int R_offset;
static int dist1;
static int FS_neighbour1;

static
void check_bins(FS_SEQ_t FS_neighbour, int i, int dist)
{
  int j;
  int k;

  for (k=i+1; k < frag_len; k++)
    {
      if (dist + D->pMclosest[qind[k]] <= D_cutoff)
	{      
	  for (j=1; j < K; j++)
	    {
	      R = K_modtable[FS_Q[k] + j];
	      R_offset = 
		FS_PARTITION_get_poffset(ptable, R);
	      FS_neighbour1 = FS_neighbour + R * FS_KK[k];
	      dist1 = dist + D->pM[qind[k]][R_offset]; 
	      HIT_LIST_count_FS_seq_visited(hit_list, 1);
	      if (dist1 <= D_cutoff)
		{
		  HIT_LIST_count_FS_seq_hit(hit_list, 1);
		  pfunc(FS_neighbour1 + FS_REM[k]);
		  check_bins(FS_neighbour1, k, dist1);
		}
	    }
	}
      FS_neighbour += FS_DC[k];
    }
  return;
}


/********************************************************************/ 
/*      Processing functions                                        */
/********************************************************************/ 

void FS_INDEX_QD_process_bin(FS_SEQ_t FS_query)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(HT, FS_query, j);

      fastadb_get_Ffrag(s_db, frag_len, &subject, frag_offset);
      if (SCORE_MATRIX_evaluate_max(D, query, &subject, 
				    D_cutoff, &D_value))
	{
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  HIT_LIST_insert_seq_hit(hit_list, &subject, 
				  (float) D_value); 
	}
      HIT_LIST_count_seq_visited(hit_list, 1);
    }
}

void FS_INDEX_S_process_bin(FS_SEQ_t FS_query)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int S_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(HT, FS_query, j); 

      fastadb_get_Ffrag(s_db, frag_len, &subject, frag_offset);
      if (SCORE_MATRIX_evaluate_min(S, query, &subject, 
				    S_cutoff, &S_value))
	{
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  HIT_LIST_insert_seq_hit(hit_list, &subject, 
				  (float) S_value); 
	}
      HIT_LIST_count_seq_visited(hit_list, 1);
    }
}


static inline
int FS_INDEX_is_empty_process_bin(FS_SEQ_t FS_query)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(HT, FS_query, j);

      fastadb_get_Ffrag(s_db, frag_len, &subject, frag_offset);

      if (SCORE_MATRIX_evaluate_max(D, query, &subject, 
				    D_cutoff, &D_value))
	return 1;
    }
  return 0;
}

void FS_INDEX_count_only_process_bin(FS_SEQ_t FS_query)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(HT, FS_query);
  HIT_LIST_count_seq_visited(hit_list, n);
}



/********************************************************************/ 
/*                 Conversion functions                             */
/********************************************************************/ 

int FS_INDEX_identity_convert(SCORE_MATRIX_t *S0, BIOSEQ *query0,
			      int cutoff)
{
  return cutoff;
}

int FS_INDEX_S2QD_convert(SCORE_MATRIX_t *S0, BIOSEQ *query0,
			  int cutoff)
{
  return SCORE_MATRIX_evaluate(S0, query0, query0) - cutoff;
}




/********************************************************************/ 
/*                 Search functions                                 */
/********************************************************************/ 


int FS_INDEX_search(HIT_LIST_t *hit_list0, BIOSEQ *query0, 
		    SCORE_MATRIX_t *S0, SCORE_MATRIX_t *D0,
		    int cutoff, FS_INDEX_process_func *pfunc0, 
		    FS_INDEX_range_convert_func *cfunc)
{
  int i;
  ULINT powKtmp = 1;
  FS_SEQ_t FS_query;

  /* Copy arguments into global variables */
  hit_list = hit_list0;
  query = query0;
  S = S0;
  D = D0;
  pfunc = pfunc0;
  S_cutoff = cutoff;
  D_cutoff = cfunc(S, query, cutoff);
  HIT_LIST_set_converted_range(hit_list, D_cutoff);
  HIT_LIST_set_index_data(hit_list, index_name, alphabet, HT->no_bins,
			  HT->no_seqs);



  /* Check length of query */
  if (query->len != frag_len)
    {
      fprintf(stderr, 
	      "FS_INDEX_search(): Query fragment length (%ld)"
	      " different from\n index fragment length (%d)."
	      "Query: %*s\n", query->len, frag_len, (int) query->len,
	      query->start);
      exit(EXIT_FAILURE);
    }

  /* Calculate FS_seq from seq */
  if (!BIOSEQ_2_FS_SEQ(query, ptable, &FS_query))
    {
      fprintf(stderr, 
	      "FS_INDEX_search(): Could not convert query"
	      " to FS fragment\n Query: %*s\n", (int) query->len, 
	      query->start);
      exit(EXIT_FAILURE);
    }
  
  /* Initialise variables */
  for (i = 0; i < frag_len; i++)
    {
      qind[i] = query->start[i] & A_SIZE_MASK;
      FS_Q[i] = FS_PARTITION_get_pttn(ptable, query->start[i]); 
      powKtmp *= K;      
      FS_DC[i] = FS_query % powKtmp - FS_query % FS_KK[i];
    }
  FS_REM[frag_len-1] = 0;  
  for (i = frag_len-2; i >= 0; i--)
    FS_REM[i] = FS_REM[i+1] + FS_DC[i+1];

  /* Process the bin associated with the query sequence */
  pfunc(FS_query); 

  /* Process neighbouring bins - recursive */
  check_bins(0, -1, 0); 

  HIT_LIST_stop_timer(hit_list);
  return 1;
}


int FS_INDEX_kNN_search(HIT_LIST_t *hit_list0, BIOSEQ *query0,
			SCORE_MATRIX_t *S0, SCORE_MATRIX_t *D0,
			ULINT k, int V0,
			FS_INDEX_process_func *pfunc0, 
			FS_INDEX_range_convert_func *cfunc)
{
  /* For now this works only for distance searches */
  int high = cfunc(S0, query0, V0);
  int l = high/2;
  ULINT size = 0;
  SEQ_HIT_t *hit_k;
  SEQ_HIT_t *hit_largest;
  float d_k;
  float d_largest;

  if (high != V0)
    {
      fprintf(stderr, "Can only run distance kNN searches!\n");
      exit(EXIT_FAILURE);
    }

  do
    {
      HIT_LIST_reset(hit_list0, query0, s_db,
		     matrix, high);
      FS_INDEX_search(hit_list0, query0, S0, D0, high,
		      pfunc0, cfunc);
      size = HIT_LIST_get_seqs_hits(hit_list0);
      if (size >= k)
	{
	  HIT_LIST_sort_incr(hit_list0);
	  hit_k = HIT_LIST_get_hit(hit_list0, k-1);
	  hit_largest = HIT_LIST_get_hit(hit_list0, size-1);
	  d_k = hit_k->value;
	  d_largest = hit_largest->value;
	  if (d_k < d_largest)
	    {
	      HIT_LIST_reset(hit_list0, query0, s_db,
			     matrix, (int) d_k);
	      FS_INDEX_search(hit_list0, query0, S0, D0, (int) d_k,
			      pfunc, cfunc);
	    }
	  return 1;
	}
      else
	  high += l;
    }
  while (1);
}

#if 0
/* Is there an epsilon neighbour */

int FS_INDEX_has_neighbour(FS_INDEX_t *FS_index, BIOSEQ *query, 
			   ULINT D_cutoff)
{
  FS_SEQ_t FS_query;
  FS_SEQ_t FS_neighbour;
  int dist;
  ULINT i;
  int FS_range;

  /* Calculate FS_seq from seq */
  if (!FS_INDEX_verify_query(FS_index, query, &FS_query))
    exit(EXIT_FAILURE);

 /* Initialise variables */
  FS_K = ptable->no_partitions;
  FS_KK[0] = 1;
  for (i = 1; i <= frag_len; i++)
    FS_KK[i] = FS_KK[i-1] * FS_K;

  /* Calculate the maximum number of changes */
  FS_range = FS_INDEX_get_no_changes(FS_index, query, D_cutoff); 

  /* Process the bin associated with the query sequence */
  if (FS_INDEX_is_empty_process_bin(FS_index, query, FS_query,
				    D_cutoff))
    return 1;

  /* Process neighbouring bins */
  i = 1;
  do 
    {
      /* Initialize combination generator*/
      kSubsetLexInit(i);
      do 
	{
	  
	  if (!SCORE_MATRIX_verify_pos(D, query, 
				       TT, i, D_cutoff))
	    continue;
	  FS_INDEX_iter_init(FS_index, FS_query);

	  while (FS_INDEX_iter_next(FS_index, query, FS_query, 
				    &dist, &FS_neighbour))
	    {
	      if (dist > D_cutoff)
		continue;
	      if (FS_INDEX_is_empty_process_bin(FS_index, 
			       query, FS_neighbour, D_cutoff))
		return 1;
	    }
	} 
      while (kSubsetLexSuccessor(i, frag_len));
      i++;
    } 
  while (i <= FS_range);

  return 0;
}
#endif


/********************************************************************/
/********************************************************************/

/********************************************************************/    
/*                                                                  */
/*                     PROFILE BASED SEARCHES                       */ 
/*                                                                  */
/********************************************************************/    

/* Profile specific functions:

   FS_INDEX_profile_S_process_bin();
   FS_INDEX_profile_D_process_bin();
   check_bin_profile();
 */

/********************************************************************/ 
/*      Processing functions                                        */
/********************************************************************/ 

void FS_INDEX_profile_S_process_bin(FS_SEQ_t FS_query)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int S_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(HT, FS_query, j);

      fastadb_get_Ffrag(s_db, frag_len, &subject, frag_offset);

      if (POS_MATRIX_evaluate_min(PS, &subject, Pfrom, Pto, 
				  S_cutoff, &S_value))
	{
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  HIT_LIST_insert_seq_hit(hit_list, &subject, 
				  (float) S_value); 
	}
      HIT_LIST_count_seq_visited(hit_list, 1);
    }
}

/* This is no longer active but left just in case */
void FS_INDEX_profile_D_process_bin(FS_SEQ_t FS_query)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(HT, FS_query, j);

      fastadb_get_Ffrag(s_db, frag_len, &subject, frag_offset);

      if (POS_MATRIX_evaluate_max(PD, &subject, Pfrom, Pto, 
				  D_cutoff, &D_value))
	{
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  HIT_LIST_insert_seq_hit(hit_list, &subject, 
				  (float) D_value); 
	}
      HIT_LIST_count_seq_visited(hit_list, 1);
    }
}

/********************************************************************/ 
/*      Recursive tree traversal function - profiles                */
/********************************************************************/ 

static
void check_bins_profile(FS_SEQ_t FS_neighbour, int i, int dist)
{
  int j;
  int k;

  for (k=i+1; k < frag_len; k++)
    {
      if (dist + PS->pMclosest[k] <= D_cutoff)
	{      
	  for (j=1; j < K; j++)
	    {
	      R = K_modtable[FS_Q[k] + j];
	      R_offset = 
		FS_PARTITION_get_poffset(ptable, R);
	      FS_neighbour1 = FS_neighbour + R * FS_KK[k];
	      dist1 = dist + PS->pM[PM_pM(k, R_offset)]; 
	      HIT_LIST_count_FS_seq_visited(hit_list, 1);
	      if (dist1 <= D_cutoff)
		{
		  HIT_LIST_count_FS_seq_hit(hit_list, 1);
		  ppfunc(FS_neighbour1 + FS_REM[k]);
		  check_bins_profile(FS_neighbour1, k, dist1);
		}
	    }
	}
      FS_neighbour += FS_DC[k];
    }
  return;
}

/********************************************************************/ 
/*                 Search functions  profiles                       */
/********************************************************************/ 

int FS_INDEX_profile_search(HIT_LIST_t *hit_list0, POS_MATRIX *PS0, 
			    ULINT Pfrom0, ULINT Pto0, int cutoff, 
			    FS_INDEX_profile_process_func *ppfunc0, 
			    FS_INDEX_range_convert_func *cfunc)
{
  int i;
  ULINT powKtmp = 1;
  FS_SEQ_t FS_query;

  /* Copy arguments into global variables */
  hit_list = hit_list0;
  query = PS0->query;
  PS = PS0;
  Pfrom = Pfrom0;
  Pto = Pto0;
  ppfunc = ppfunc0;
  S_cutoff = cutoff;
  D_cutoff = POS_MATRIX_evaluate(PS, query, Pfrom, Pto) - cutoff;
  HIT_LIST_set_converted_range(hit_list, D_cutoff);
  HIT_LIST_set_index_data(hit_list, index_name, alphabet, 
			  HT->no_bins, HT->no_seqs);

  /* Check length of query */
  if (query->len != frag_len)
    {
      fprintf(stderr, 
	      "FS_INDEX_search(): Query fragment length (%ld)"
	      " different from\n index fragment length (%d)."
	      "Query: %*s\n", query->len, frag_len, (int) query->len,
	      query->start);
      exit(EXIT_FAILURE);
    }

  /* Calculate FS_seq from seq */
  if (!BIOSEQ_2_FS_SEQ(query, ptable, &FS_query))
    {
      fprintf(stderr, 
	      "FS_INDEX_search(): Could not convert query"
	      " to FS fragment\n Query: %*s\n", (int) query->len, 
	      query->start);
      exit(EXIT_FAILURE);
    }
  
  /* Initialise variables */
  for (i = 0; i < frag_len; i++)
    {
      FS_Q[i] = FS_PARTITION_get_pttn(ptable, query->start[i]); 
      powKtmp *= K;      
      FS_DC[i] = FS_query % powKtmp - FS_query % FS_KK[i];
    }
  FS_REM[frag_len-1] = 0;  
  for (i = frag_len-2; i >= 0; i--)
    FS_REM[i] = FS_REM[i+1] + FS_DC[i+1];

  /* Process the bin associated with the query sequence */
  ppfunc(FS_query); 

  /* Process neighbouring bins - recursive */
  check_bins_profile(0, -1, 0); 

  HIT_LIST_stop_timer(hit_list);
  return 1;
}


/*********************************************************************/
/*********************************************************************/
