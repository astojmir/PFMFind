#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "hit_list.h"
#include "partition.h"
#include "smatrix.h"
#include "FSindex.h"

#ifdef TRY_CLUSTERS
#include "clusters.h"
#endif

#define MERGE_DUPLICATES

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
  FS_HASH_TABLE_t *HT = mallocec(sizeof(FS_HASH_TABLE_t));
 
  HT->no_bins = no_bins;
  HT->no_seqs = 0;
  HT->no_useqs = 0;

  HT->shist_len = 5 * def_size;
  HT->shist = callocec(HT->shist_len, sizeof(ULINT)); 
  HT->uhist_len = 5 * def_size;
  HT->uhist = callocec(HT->uhist_len, sizeof(ULINT)); 
  HT->bin_size = callocec(no_bins, sizeof(ULINT));
  HT->max_bin_size = callocec(no_bins, sizeof(ULINT)); 
  HT->bin = mallocec(no_bins * sizeof(SEQ_index_t *));
  HT->binL = 0;
  HT->heap = NULL;

  HT->u_size = callocec(no_bins, sizeof(ULINT));
  HT->u = mallocec(no_bins * sizeof(int));
  HT->u_heap = NULL;
  HT->uL = 0;

  return HT;  
}

/* Destructor */
void FS_HASH_TABLE_destroy(FS_HASH_TABLE_t *HT)
{
  ULINT i;

  if (HT->heap == NULL)
    for (i=0; i < HT->no_bins; i++)
      free(HT->bin[i]);
  else
    free(HT->heap);

  /* TO DO: Free unique sequneces */

  free(HT->bin);
  free(HT->bin_size);
  free(HT->shist);
  if (HT->max_bin_size != NULL)
    free(HT->max_bin_size);
  free(HT);
}

/* Insertion */
void FS_HASH_TABLE_count_seq(FS_HASH_TABLE_t *HT, ULINT i,
			      SEQ_index_t seq_i)
{
  /* Count the sequence */
  HT->max_bin_size[seq_i]++;
  HT->no_seqs++;
}

void FS_HASH_TABLE_allocate_bins(FS_HASH_TABLE_t *HT)
{
  ULINT i;
  for (i=0; i < HT->no_bins; i++)
    HT->bin[i] = mallocec(HT->max_bin_size[i] 
			     * sizeof(SEQ_index_t));
}

void FS_HASH_TABLE_insert_seq(FS_HASH_TABLE_t *HT, ULINT i,
			      SEQ_index_t seq_i)
{
  /* Insert the sequence */
  HT->bin[seq_i][HT->bin_size[seq_i]] = i;
  HT->bin_size[seq_i]++;
}

void FS_HASH_TABLE_add_bin_stats(FS_HASH_TABLE_t *HT,
				 SEQ_index_t seq_i)
{
  ULINT s0; /* Old length of frequency vector */

  if (HT->bin_size[seq_i] > HT->bin_size[HT->binL])
    HT->binL = seq_i;
    
  if (HT->bin_size[seq_i] >= HT->shist_len)
    {
      s0 = HT->shist_len;
      HT->shist_len = HT->bin_size[seq_i]+1;
      HT->shist = reallocec(HT->shist, HT->shist_len * sizeof(ULINT)); 
      memset(HT->shist + s0, 0, (HT->shist_len - s0) * sizeof(ULINT));
    }
  HT->shist[HT->bin_size[seq_i]]++;

#ifdef MERGE_DUPLICATES
  if (HT->u_size[seq_i] > HT->u_size[HT->uL])
    HT->uL = seq_i;

  if (HT->u_size[seq_i] >= HT->uhist_len)
    {
      s0 = HT->uhist_len;
      HT->uhist_len = HT->u_size[seq_i]+1;
      HT->uhist = reallocec(HT->uhist, HT->uhist_len * sizeof(ULINT));  
      memset(HT->uhist + s0, 0, (HT->uhist_len - s0) * sizeof(ULINT));
    }
  HT->uhist[HT->u_size[seq_i]]++;
#endif
}

/* Trimming */
void FS_HASH_TABLE_resize(FS_HASH_TABLE_t *HT)
{
  ULINT i;
  for (i=0; i < HT->no_bins; i++)
    {
      if (HT->bin_size[i] > 0)
	{
	  HT->bin[i] = reallocec(HT->bin[i], 
	  		HT->bin_size[i] * sizeof(SEQ_index_t)); 
	  HT->max_bin_size[i] = HT->bin_size[i];
	}      
    }
}


/* Reading, writing to file */
int FS_HASH_TABLE_write(FS_HASH_TABLE_t *HT, FILE *stream)
{
  ULINT i;

  fwrite(&(HT->no_bins), sizeof(ULINT), 1, stream);

  fwrite(&(HT->no_seqs), sizeof(ULINT), 1, stream);
  fwrite(&(HT->shist_len), sizeof(ULINT), 1, stream);
  fwrite(&(HT->binL), sizeof(SEQ_index_t), 1, stream);
  fwrite(HT->shist, sizeof(ULINT), HT->shist_len, stream); 
  fwrite(HT->bin_size, sizeof(ULINT), HT->no_bins, stream);  

  for (i=0; i < HT->no_bins; i++)
    fwrite(HT->bin[i], sizeof(SEQ_index_t), HT->bin_size[i], stream);

  fwrite(&(HT->no_useqs), sizeof(ULINT), 1, stream);
  fwrite(&(HT->uhist_len), sizeof(ULINT), 1, stream);
  fwrite(&(HT->uL), sizeof(SEQ_index_t), 1, stream);
  fwrite(HT->uhist, sizeof(ULINT), HT->uhist_len, stream); 
  fwrite(HT->u_size, sizeof(ULINT), HT->no_bins, stream);  
  
  for (i=0; i < HT->no_bins; i++)
    fwrite(HT->u[i], sizeof(int), HT->u_size[i], stream);

  return 1;
}

FS_HASH_TABLE_t *FS_HASH_TABLE_read(FILE *stream)
{
  FS_HASH_TABLE_t *HT = mallocec(sizeof(FS_HASH_TABLE_t));
  ULINT i;
  SEQ_index_t *current;
  int *t2;

  fread(&(HT->no_bins), sizeof(ULINT), 1, stream);
  HT->bin_size = mallocec(HT->no_bins * sizeof(ULINT));
  HT->u_size = mallocec(HT->no_bins * sizeof(ULINT));
  HT->bin = mallocec(HT->no_bins * sizeof(SEQ_index_t *));
  HT->u = mallocec(HT->no_bins * sizeof(int *));

  fread(&(HT->no_seqs), sizeof(ULINT), 1, stream);
  HT->heap = mallocec(HT->no_seqs * sizeof(SEQ_index_t)); 

  fread(&(HT->shist_len), sizeof(ULINT), 1, stream);
  HT->shist = mallocec(HT->shist_len * sizeof(ULINT));

  fread(&(HT->binL), sizeof(SEQ_index_t), 1, stream);

  HT->max_bin_size = NULL;

  fread(HT->shist, sizeof(ULINT), HT->shist_len, stream); 
  fread(HT->bin_size, sizeof(ULINT), HT->no_bins, stream);  

  current =  HT->heap;
  for (i=0; i < HT->no_bins; i++)
    {
      HT->bin[i] = current;
      fread(HT->bin[i], sizeof(SEQ_index_t), HT->bin_size[i],
	    stream);
      current += HT->bin_size[i];
    }

  fread(&(HT->no_useqs), sizeof(ULINT), 1, stream);
  HT->u_heap = mallocec(HT->no_useqs * sizeof(int)); 

  fread(&(HT->uhist_len), sizeof(ULINT), 1, stream);
  HT->uhist = mallocec(HT->uhist_len * sizeof(ULINT)); 

  fread(&(HT->uL), sizeof(SEQ_index_t), 1, stream);

  fread(HT->uhist, sizeof(ULINT), HT->uhist_len, stream); 
  fread(HT->u_size, sizeof(ULINT), HT->no_bins, stream);  
  
  t2 = HT->u_heap;
  for (i=0; i < HT->no_bins; i++)
    {
      HT->u[i] = t2;
      fread(HT->u[i], sizeof(int), HT->u_size[i], stream);
      t2 += HT->u_size[i];
    }
   
  return HT;  
}

/* Statistics */
void FS_HASH_TABLE_print_stats(FS_HASH_TABLE_t *HT, FILE *stream,
			       FS_PARTITION_t *FS_partition, 
			       ULINT frag_len)
{
  ULINT i;
  ULINT CF;
  ULINT CS;
  
  fprintf(stream, "Total number of fragments in index: %ld\n", 
	  HT->no_seqs);
  fprintf(stream, "Total number of index entries: %ld\n", 
	  HT->no_bins);
  fprintf(stream, "Average size of index entry: %.2f\n", 
	  (float) HT->no_seqs / (float) HT->no_bins);
  fprintf(stream, "Largest index entry: %s (%ld)\n",
	  FS_seq_print(HT->binL, FS_partition, frag_len),
	  HT->bin_size[HT->binL]);
  fprintf(stream, "Total number of distinct fragments in index: %ld\n", 
	  HT->no_useqs);
  fprintf(stream, "Average 'distinct' size of index entry: %.2f\n", 
	  (float) HT->no_useqs / (float) HT->no_bins);
  fprintf(stream, "Largest 'distinct' index entry: %s (%ld)\n",
	  FS_seq_print(HT->uL, FS_partition, frag_len),
	  HT->u_size[HT->uL]);

  fprintf(stream, "\n* Distribution of numbers of fragments per"
	  " index entry *\n");
  fprintf(stream, "%17s %10s\n","Index entry size", "Frequency");
  CF = 0;
  CS = 0;
  for (i = 0; i < HT->shist_len; i++)
    {
      CF += HT->shist[i];
      CS += HT->shist[i] * i;    
      if (HT->shist[i] > 0)
	fprintf(stream, "%17ld %10ld %10ld %10ld\n", i, HT->shist[i], 
		CF, CS);
    }
  fprintf(stream, "\n");

  fprintf(stream, "\n* Distribution of numbers of distinct fragments per"
	  " index entry *\n");
  fprintf(stream, "%17s %10s\n","Index entry size", "Frequency");
  CF = 0;
  CS = 0;
  for (i = 0; i < HT->uhist_len; i++)
    {
      CF += HT->uhist[i];
      CS += HT->uhist[i] * i;
      if (HT->uhist[i] > 0)
	fprintf(stream, "%17ld %10ld %10ld %10ld\n", i, HT->uhist[i], 
		CF, CS);
    }
  fprintf(stream, "\n");
}

ULINT FS_HASH_TABLE_get_total_seqs(FS_HASH_TABLE_t *HT)
{
  return HT->no_seqs;
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
/* static char *matrix; */
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
static int kNN;
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


#ifdef MERGE_DUPLICATES
/* Fragment comparison routine */
static 
int comp_frags(const void *fg1, const void *fg2)
{
  SEQ_index_t *f1 = (SEQ_index_t *) fg1;
  SEQ_index_t *f2 = (SEQ_index_t *) fg2;
  BIOSEQ s1;
  BIOSEQ s2;

  fastadb_get_Ffrag(s_db, frag_len, &s1, *f1);
  fastadb_get_Ffrag(s_db, frag_len, &s2, *f2);
  return memcmp(s1.start, s2.start, frag_len);
}
#endif


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

#ifdef MERGE_DUPLICATES
  int k;
#endif


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

  /* Make initial bin size half as large as mean size */
  bin_size = (int)((float) no_frags / no_FS) + 1;

  /* Initialise hash table */
  HT = FS_HASH_TABLE_create(no_FS, bin_size);

  /* We do it in two passes - first we count the fragments, then put
     them */
  if (FS_INDEX_VERBOSE || FS_INDEX_PRINT_BAR)
    printf("Counting fragments ...\n");
  /* Counting */ 
  j = 0;
  while (fastadb_get_next_Ffrag(s_db, frag_len, frag, &i, skip))  
    {
      j++;
      /* Calculate its FS_seq */
      if (BIOSEQ_2_FS_SEQ(frag, ptable, &FS_seq))
	FS_HASH_TABLE_count_seq(HT, i, FS_seq);
	
      /* Print progress bar */
      if (FS_INDEX_PRINT_BAR > 0)
	printbar(stdout, j+1, one_percent_fragments, 50);  
    }  
  /* Allocating */
  if (FS_INDEX_VERBOSE || FS_INDEX_PRINT_BAR)
    printf("Allocating bins ...\n");
  FS_HASH_TABLE_allocate_bins(HT);

  /* Commiting */
  if (FS_INDEX_VERBOSE || FS_INDEX_PRINT_BAR)
    printf("Collecting fragments ...\n");
  j = 0;
  fastadb_init_Ffrags(s_db, frag_len);
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
  if (FS_INDEX_VERBOSE || FS_INDEX_PRINT_BAR)
    printf("Collecting duplicate fragments ...\n");

  for (j=0; j < HT->no_bins; j++)
    {

#ifdef MERGE_DUPLICATES
      if (HT->bin_size[j] != 0) {

      /* Sort fragments */
      qsort(HT->bin[j], HT->bin_size[j], sizeof(SEQ_index_t),
	    comp_frags);
      /* Count unique fragments */
      HT->u_size[j] = 1;
      for (i=1; i < HT->bin_size[j]; i++)
	if (comp_frags(HT->bin[j]+i, HT->bin[j]+i-1) != 0)
	  HT->u_size[j]++;

      HT->no_useqs += HT->u_size[j];
      /* Allocate unique fragments */
      HT->u[j] = callocec(HT->u_size[j], sizeof(int));

      /* Fill unique array */
 
      HT->u[j][0] = 0;
      k = 0;
      for (i=1; i < HT->bin_size[j]; i++)
	if (comp_frags(HT->bin[j]+i, HT->bin[j]+i-1) != 0)
	  {
	    k++;
	    HT->u[j][k] = i; 
	  }
      }
      if (FS_INDEX_PRINT_BAR > 0)
	printbar(stdout, j+1, HT->no_bins/100, 50);  
      
#endif
      /* Add statistics */
      FS_HASH_TABLE_add_bin_stats(HT, j);     
    }

#if 0
  /* Resize bins */
  FS_HASH_TABLE_resize(HT); 
#endif

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

const char *FS_INDEX_get_db_name(void)
{
  return db_name;
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

const char *FS_index_get_alphabet(void)
{
  return alphabet;
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
  ULINT n = FS_HASH_TABLE_get_no_useqs(HT, FS_query);
  ULINT N = FS_HASH_TABLE_get_no_seqs(HT, FS_query);
  int i;
  int j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (i = 0; i < n; i++)
    {
      j = HT->u[FS_query][i];

      frag_offset = FS_HASH_TABLE_retrieve_seq(HT, FS_query, j);
      fastadb_get_Ffrag(s_db, frag_len, &subject, frag_offset);

      if (SCORE_MATRIX_evaluate_max(D, query, &subject, 
				    D_cutoff, &D_value))
	{
	  HIT_LIST_count_useq_hit(hit_list, 1);
	  HIT_LIST_insert_seq_hit(hit_list, &subject, (float) D_value); 
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  j++;
	  while ((j < N && i == n-1) || j < HT->u[FS_query][i+1])
	    {
	      HIT_LIST_count_seq_hit(hit_list, 1);
	      frag_offset = FS_HASH_TABLE_retrieve_seq(HT, FS_query, j);
	      fastadb_get_Ffrag(s_db, frag_len, &subject, frag_offset);
	      HIT_LIST_insert_seq_hit(hit_list, &subject, (float) D_value); 
	      j++;
	    }

	}
      HIT_LIST_count_seq_visited(hit_list,
	i == n-1 ? N - HT->u[FS_query][i] :
                   HT->u[FS_query][i+1]- HT->u[FS_query][i]);
      HIT_LIST_count_useq_visited(hit_list, 1);
    }
}

static
void FS_INDEX_kNN_QD_process_bin(FS_SEQ_t FS_query)
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
      if (SCORE_MATRIX_evaluate_max(D, query, &subject, D_cutoff, &D_value) ||
	  HIT_LIST_get_seqs_hits(hit_list) < kNN)
	{
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  D_cutoff = 
	    HIT_LIST_insert_seq_hit_queue(hit_list, &subject, 
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
  HIT_LIST_count_FS_seq_visited(hit_list, 1);
  HIT_LIST_count_FS_seq_hit(hit_list, 1);

  /* Process neighbouring bins - recursive */
  check_bins(0, -1, 0); 

  HIT_LIST_stop_timer(hit_list);
  return 1;
}


int FS_INDEX_kNN_search(HIT_LIST_t *hit_list0, BIOSEQ *query0,
			SCORE_MATRIX_t *S0, SCORE_MATRIX_t *D0,
			ULINT k)
{
  int i;
  ULINT powKtmp = 1;
  FS_SEQ_t FS_query;

  /* Copy arguments into global variables */
  kNN = k;
  hit_list = hit_list0;
  HIT_LIST_set_kNN(hit_list, kNN); 
  query = query0;
  S = S0;
  D = D0;
  pfunc = FS_INDEX_kNN_QD_process_bin;
  S_cutoff = 0;
  D_cutoff = INT_MAX;
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
  HIT_LIST_count_FS_seq_visited(hit_list, 1);
  HIT_LIST_count_FS_seq_hit(hit_list, 1);

  /* Process neighbouring bins - recursive */
  check_bins(0, -1, 0); 

  HIT_LIST_sort_kNN(hit_list);
  HIT_LIST_set_converted_range(hit_list, D_cutoff);
  HIT_LIST_stop_timer(hit_list);

  return 1;
}


/********************************************************************/
/********************************************************************/

/********************************************************************/    
/*                                                                  */
/*                     FS_INDEX_has_neighbour                       */ 
/*                                                                  */
/********************************************************************/    

/* Is there an epsilon neighbour? */
/* We use a copy of check_bin for efficiency */


static int empty_nbhd = 1;

static
void check_bins_empty(FS_SEQ_t FS_neighbour, int i, int dist)
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
	      if (dist1 <= D_cutoff)
		{
		  if (FS_INDEX_is_empty_process_bin(FS_neighbour1 
						    + FS_REM[k]))
		    {
		      empty_nbhd = 0;
		      return;
		    }

		  check_bins_empty(FS_neighbour1, k, dist1);
		  if (empty_nbhd == 0)
		    return;
		}
	    }
	}
      FS_neighbour += FS_DC[k];
    }
  return;
}


int FS_INDEX_has_neighbour(BIOSEQ *query0, SCORE_MATRIX_t *D0, 
			   int cutoff)
{
  int i;
  ULINT powKtmp = 1;
  FS_SEQ_t FS_query;

  /* Copy arguments into global variables */
  query = query0;
  D = D0;

  D_cutoff = cutoff;
  empty_nbhd = 1;

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

  if (FS_INDEX_is_empty_process_bin(FS_query))
    return 1;

  /* Process neighbouring bins - recursive */
  check_bins_empty(0, -1, 0); 


  if (empty_nbhd)
    return 0;
  else
    return 1;
}



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
