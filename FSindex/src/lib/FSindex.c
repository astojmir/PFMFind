#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "hit_list.h"
#include "partition.h"
#include "smatrix.h"


/* If the extern inline functions are not to be inlined, they must be
   here */ 
#ifdef DEBUG
#define FS_INDEX_INLINE
#endif

#include "FSindex.h"

#ifdef DEBUG
#undef FS_INDEX_INLINE
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

/* Main constructor */
FS_INDEX_t *FS_INDEX_create(const char *db_name, ULINT frag_len,
			    const char *alphabet, 
			    const char separator)
 
{
  time_t time0 = time(NULL);
  time_t time1;
  double dt;
  ULINT one_percent_fragments;
  char *db_dir;
  char *db_base;

  ULINT i, j;
  ULINT K;
  ULINT no_FS;
  ULINT no_frags;
  BIOSEQ *frag = mallocec(sizeof(BIOSEQ));
  FS_SEQ_t FS_seq;

  FS_INDEX_t *FS_index = callocec(1, sizeof(FS_INDEX_t)); 
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  ULINT bin_size;

  /* Initialise FS_index */
  split_base_dir(db_name, &db_base, &db_dir);
  FS_index->db_name = strdup(db_base);
  FS_index->alphabet = strdup(alphabet);
  FS_index->separator = separator;
  FS_index->frag_len = frag_len;

  /* Load database */
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = NONE;
  fastadb_argv[0].access_type = MEMORY;
  fastadb_argv[1].retrieve_deflines = NO;

  FS_index->s_db = fastadb_open(db_name, fastadb_argt, fastadb_argv); 
 
  /* Make fragment database */
  no_frags = fastadb_count_Ffrags(FS_index->s_db, frag_len);
  one_percent_fragments = no_frags / 100;
  fastadb_init_Ffrags(FS_index->s_db, frag_len);

  /* Create partition table */
  FS_index->ptable = FS_PARTITION_create(FS_index->alphabet,
					 separator); 

  /* Error if the sequence is too long for the partitioning*/
  K = FS_PARTITION_get_no_partitions(FS_index->ptable);
  if ((K-1) >> ((sizeof(ULINT)*8)/frag_len))
    {
      /* Clean up */
      FS_INDEX_destroy(FS_index);
      return NULL;
    }

  /* Number of bins to allocate in hash table */
  no_FS = 1;
  for (i=0; i < frag_len; i++)
    no_FS *= K;

  /* Make initial bin size slightly bigger than average */
  bin_size = ((5 * no_frags) / (4 * no_FS)) + 1;

  /* Initialise hash table */
  FS_index->HT = FS_HASH_TABLE_create(no_FS, bin_size);

  /* For each fragment */ 
  j = 0;
  while (fastadb_get_next_Ffrag(FS_index->s_db, frag_len, frag, &i))  
    {
      j++;
      /* Calculate its FS_seq */
      if (BIOSEQ_2_FS_SEQ(frag, FS_index->ptable, &FS_seq))
	{
	  /* Add to appropriate FS_bin */
	  FS_HASH_TABLE_insert_seq(FS_index->HT, i, FS_seq); 
	}
      /* Print progress bar */
      if (FS_INDEX_PRINT_BAR > 0)
	printbar(stdout, j+1, one_percent_fragments, 50);  
    }  
  /* Resize bins */
  FS_HASH_TABLE_resize(FS_index->HT); 

  /* Take time */
  time1 = time(NULL);
  dt = difftime(time1, time0)/60.0;
  FS_index->db_no_frags = no_frags;

  /* Print statistics */
  if (FS_INDEX_VERBOSE > 0)
    FS_INDEX_print_stats(FS_index, stdout, j, dt); 

  /* Clean up */
  free(frag);
  free(db_dir);
  free(db_base);

  /* Return an index ready for search */
  return FS_index; 
}


/* Destructor */
void FS_INDEX_destroy(FS_INDEX_t *FS_index)
{
  /* Free all bins */
  if (FS_index->HT != NULL)
    FS_HASH_TABLE_destroy(FS_index->HT); 

  /* Free partition table */
  FS_PARTITION_destroy(FS_index->ptable);  

  /* Free sequence database */
  fastadb_close(FS_index->s_db);

  /* Free names */
  free(FS_index->db_name);
  free(FS_index->alphabet);

  /* Free FS_index */
  free(FS_index);
}

/* Load, Save */
int FS_INDEX_save(FS_INDEX_t *FS_index, const char *filename)
{ 
  FILE *stream = fopen(filename, "wb");
  UINT_t len;
 
  /* Open file */
  if(stream == NULL)
    return 0;

  len = strlen(FS_index->db_name) + 1;
  fwrite(&len, sizeof(UINT_t), 1, stream);
  fwrite(FS_index->db_name, sizeof(char), len, stream);

  len = strlen(FS_index->alphabet) + 1;
  fwrite(&len, sizeof(UINT_t), 1, stream);
  fwrite(FS_index->alphabet, sizeof(char), len, stream);

  fwrite(&(FS_index->separator), sizeof(char), 1, stream);
  fwrite(&(FS_index->frag_len), sizeof(ULINT), 1, stream);
  fwrite(&(FS_index-> db_no_frags), sizeof(ULINT), 1, stream);

  FS_PARTITION_write(FS_index->ptable, stream); 
  FS_HASH_TABLE_write(FS_index->HT, stream);

  fclose(stream);
  return 1;
}

FS_INDEX_t *FS_INDEX_load(const char *filename)
{
  FILE *stream = fopen(filename, "rb");
  UINT_t len;
  FS_INDEX_t *FS_index = callocec(1, sizeof(FS_INDEX_t)); 
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  char *basename;
  char *dirname;
  char *full_dbname;
 
  /* Open file */
  if(stream == NULL)
    {
      fprintf(stderr, "Could not open index file %s!\n",
	      filename);
      exit(EXIT_FAILURE);
    }
  fread(&len, sizeof(UINT_t), 1, stream);
  FS_index->db_name = mallocec(len);
  fread(FS_index->db_name, sizeof(char), len, stream);

  split_base_dir(filename, &basename, &dirname);
  cat_base_dir(&full_dbname, FS_index->db_name, dirname);

  fread(&len, sizeof(UINT_t), 1, stream);
  FS_index->alphabet = mallocec(len);
  fread(FS_index->alphabet, sizeof(char), len, stream);

  fread(&(FS_index->separator), sizeof(char), 1, stream);
  fread(&(FS_index->frag_len), sizeof(ULINT), 1, stream);
  fread(&(FS_index-> db_no_frags), sizeof(ULINT), 1, stream);

  /* Load database */
  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = NONE;
  fastadb_argv[0].access_type = MEMORY;
  fastadb_argv[1].retrieve_deflines = NO;

  FS_index->s_db = fastadb_open(full_dbname, fastadb_argt, 
				fastadb_argv); 
  if (FS_index->s_db == NULL)
    {
      fprintf(stderr, "Could not open database file %s!\n",
	      full_dbname);
      exit(EXIT_FAILURE);
    }
      
  FS_index->ptable = FS_PARTITION_read(stream); 
  FS_index->HT = FS_HASH_TABLE_read(stream);

  fclose(stream);
  free(basename);
  free(dirname);
  return FS_index;
}

/* Print */
void FS_INDEX_print_stats(FS_INDEX_t *FS_index, FILE *stream, 
			  ULINT count, double dtime)
{
  fprintf(stream, "\n*** FS_INDEX Statistics ***\n");
  fprintf(stream, "Database: %s\n", FS_index->db_name);
  fprintf(stream, "Database size: %ld\n", FS_index->s_db->length); 
  fprintf(stream, "Number of sequences: %ld\n",
	  FS_index->s_db->no_seq);  
  fprintf(stream, "Partitions: %s\n", FS_index->alphabet);
  fprintf(stream, "Fragment Length: %ld\n", FS_index->frag_len);
  fprintf(stream, "Total Fragments in database: %ld\n",
	  FS_index->db_no_frags); 
  if (count > 0)
    fprintf(stream, "Total Counted Fragments in database: %ld\n", 
	    count);
  if (dtime > 0)
    fprintf(stream, "Creation Time: %.2lf mins\n\n", dtime);
  FS_HASH_TABLE_print_stats(FS_index->HT, stdout, FS_index->ptable,
			    FS_index->frag_len); 
}


/* Search */



/* Static variables and functions declarations */

static inline
int FS_INDEX_verify_query(FS_INDEX_t *FS_index, BIOSEQ *query, 
			  FS_SEQ_t *FS_seq)
{
  if (query->len != FS_index->frag_len)
    return 0;

  if (!BIOSEQ_2_FS_SEQ(query, FS_index->ptable, FS_seq))
    return 0;

  return 1;
}


static 
int FS_INDEX_get_no_changes(FS_INDEX_t *FS_index, BIOSEQ *query, 
			    int D_cutoff)
{
  UINT_t D[32];
  int i;
  int S = 0;

  for (i = 0; i < FS_index->frag_len; i++)
    {
      D[i] = 
	FS_index->D->pMclosest[query->start[i] & A_SIZE_MASK];
    }  
  qsort(D, FS_index->frag_len, sizeof(UINT_t), compare_int);
  for (i = 0; i < FS_index->frag_len; i++)
    {
      S += D[i];
      if (S > D_cutoff)
	break;
    }
  return i;
}


/* Combinations -lexicographic order */

static USINT TT[32]; /* Vector of positions to be changed */ 
static USINT UU[32]; /* Copy of TT */

static int FS_K; /* Number of letters in FS alphabet */
static int FS_r; /* Rank of the current set of changed residues */
static int FS_l; /* Number of changed residues */
static USINT FS_Q[32]; /* Residues in the FS_query sequence in
			  the changed positions given by TT.
			  i.e. FS_Q[TT[i]] is the i-th residue
			  to be changed. */
static FS_SEQ_t FS_S0; /* Query sequence where the positions to be
			  changed are set to 0 */
static ULINT FS_KTT[32]; /* FS_KTT[i] = FS_K ^ TT[i] */
static ULINT FS_KK[32]; /* FS_KK[i] = FS_K ^ i */
static ULINT FS_Kl; /*  (FS_K-1) ^ FS_l */


static inline
void kSubsetLexInit(int k)
{
  ULINT i;

  FS_Kl = 1;
  FS_l = k;

  for (i = 0; i < k; i++)
    {
      TT[i] = i;
      FS_Kl *= (FS_K-1);
    }
  memcpy(UU, TT, k*sizeof(USINT));
}

static inline
int kSubsetLexSuccessor(int k, int n)
{
  int j;
  int i = k - 1;
  
  while ((i >= 0) && TT[i] == n - k + i)
    i--;
  if (i == -1)
    return 0;
  
  for (j=i; j < k; j++)
    UU[j] = TT[i] + 1 + j - i;

  memcpy(TT, UU, k*sizeof(USINT));
  return 1;
}

static inline
void FS_INDEX_iter_init(FS_INDEX_t *FS_index, FS_SEQ_t FS_query)
{
  ULINT j;

  FS_r = 0;
  FS_S0 = FS_query;
  j = 0;
  for (j = 0; j < FS_l; j++)
    {
      FS_KTT[j] = FS_KK[TT[j]];
      FS_Q[j] = (FS_query / FS_KTT[j]) % FS_K;
      FS_S0 -= FS_Q[j] * FS_KTT[j];
    }
}

static inline
int FS_INDEX_iter_next(FS_INDEX_t *FS_index, BIOSEQ *query, 
		       FS_SEQ_t FS_query, int *dist,  
		       FS_SEQ_t *FS_neighbour) 
{
  ULINT i;
  int r = FS_r;
  int R;

  if (FS_r >= FS_Kl) return 0;
  *dist = 0;
  *FS_neighbour = FS_S0;
  for (i = 0; i < FS_l; i++)
    {
      R = ((r % (FS_K-1)) + FS_Q[i] + 1) % FS_K;
      *dist += FS_index->D->
	pM[query->start[TT[i]] & A_SIZE_MASK][R];
      *FS_neighbour += R * FS_KTT[i];
      r /= (FS_K-1);
    }
  FS_r++;
  return 1;
}

inline
void FS_INDEX_QD_process_bin(FS_INDEX_t *FS_index, 
			     HIT_LIST_t *hit_list, 
			     BIOSEQ *query,
			     FS_SEQ_t FS_query,
			     int D_cutoff)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(FS_index->HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(FS_index->HT, 
					       FS_query, j);

      fastadb_get_Ffrag(FS_index->s_db, FS_index->frag_len, 
			&subject, frag_offset);
      if (SCORE_MATRIX_evaluate_max(FS_index->D, query, &subject, 
				    D_cutoff, &D_value))
	{
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  HIT_LIST_insert_seq_hit(hit_list, &subject, 
				  (float) D_value); 
	}
      HIT_LIST_count_seq_visited(hit_list, 1);
    }
}

void FS_INDEX_S_process_bin(FS_INDEX_t *FS_index, 
			    HIT_LIST_t *hit_list, 
			    BIOSEQ *query,
			    FS_SEQ_t FS_query,
			    int S_cutoff)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(FS_index->HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int S_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(FS_index->HT, 
					       FS_query, j);

      fastadb_get_Ffrag(FS_index->s_db, FS_index->frag_len, 
			&subject, frag_offset);
      if (SCORE_MATRIX_evaluate_min(FS_index->S, query, &subject, 
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
int FS_INDEX_is_empty_process_bin(FS_INDEX_t *FS_index, 
				   BIOSEQ *query,
				   FS_SEQ_t FS_query,
				   int D_cutoff)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(FS_index->HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(FS_index->HT, 
					       FS_query, j);

      fastadb_get_Ffrag(FS_index->s_db, FS_index->frag_len, 
			&subject, frag_offset);

      if (SCORE_MATRIX_evaluate_max(FS_index->D, query, &subject, 
				    D_cutoff, &D_value))
	return 1;
    }
  return 0;
}

void FS_INDEX_count_only_process_bin(FS_INDEX_t *FS_index, 
				     HIT_LIST_t *hit_list, 
				     BIOSEQ *query,
				     FS_SEQ_t FS_query,
				     int D_cutoff)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(FS_index->HT, FS_query);
  HIT_LIST_count_seq_visited(hit_list, n);
}



/* Conversion functions */
int FS_INDEX_identity_convert(FS_INDEX_t *FS_index,
			      BIOSEQ *query,
			      int cutoff)
{
  return cutoff;
}

int FS_INDEX_S2QD_convert(FS_INDEX_t *FS_index,
			  BIOSEQ *query,
			  int cutoff)
{
  return SCORE_MATRIX_evaluate(FS_index->S, query, query) -
    cutoff;
}




/* Searches */

HIT_LIST_t *FS_INDEX_QD_search(FS_INDEX_t *FS_index, 
			    BIOSEQ *query, ULINT D_cutoff)
{
  FS_SEQ_t FS_query;
  FS_SEQ_t FS_neighbour;
  int dist;
  HIT_LIST_t *hit_list;

  ULINT i;
  int FS_range;
#if DEBUG > 8
  ULINT j;
  int debug_flag =0;
#endif

  /* Calculate FS_seq from seq */
  if (!FS_INDEX_verify_query(FS_index, query, &FS_query))
    return NULL;

#if DEBUG > 8
  fprintf(stdout, "Query sequence: %.10s\n", query->start);
  fprintf(stdout, "Query sequence: %s\n",
	  FS_seq_print(FS_query, FS_index->ptable, 10));
#endif


 /* Initialise HIT_LIST and other variables */
  hit_list = HIT_LIST_create(query, FS_index->s_db, 
			     FS_index->matrix, D_cutoff);
  FS_K = FS_index->ptable->no_partitions;
  FS_KK[0] = 1;
  for (i = 1; i <= FS_index->frag_len; i++)
    FS_KK[i] = FS_KK[i-1] * FS_K;

  /* Calculate the maximum number of changes */
  FS_range = FS_INDEX_get_no_changes(FS_index, query, D_cutoff); 

#if DEBUG > 8
  fprintf(stdout, "Maximum number of changes: %d\n", FS_range);
#endif
  
  /* Process the bin associated with the query sequence */
  FS_INDEX_QD_process_bin(FS_index, hit_list, query, FS_query, D_cutoff); 

  /* Process neighbouring bins */
  i = 1;
  do {
    /* Initialize combination generator*/
  kSubsetLexInit(i);
    do {

#if DEBUG > 8
    fprintf(stdout, "TT    : ");
    for (j=0; j < i; j++)
      fprintf(stdout, "%d ", TT[j]);
    fprintf(stdout, "\n");
#endif

    if (!SCORE_MATRIX_verify_pos(FS_index->D, query, 
				 TT, i, D_cutoff))
      continue;
    FS_INDEX_iter_init(FS_index, FS_query);

#if DEBUG > 8
    fprintf(stdout, "FS_Q  : ");
    for (j=0; j < i; j++)
      fprintf(stdout, "%d ", FS_Q[j]);
      fprintf(stdout, "\n");
#endif

#if DEBUG > 8
      fprintf(stdout, "S0: %s\n",
	      FS_seq_print(FS_S0, FS_index->ptable, 10));
#endif
      while (FS_INDEX_iter_next(FS_index, query, FS_query, 
				&dist, &FS_neighbour))
	{
	  HIT_LIST_count_FS_seq_visited(hit_list, 1);
	  if (dist > D_cutoff)
	    continue;
	  HIT_LIST_count_FS_seq_hit(hit_list, 1);
	  FS_INDEX_QD_process_bin(FS_index, hit_list, query,
				  FS_neighbour, D_cutoff);
	}
    } while (kSubsetLexSuccessor(i, FS_index->frag_len));
    i++;
  } while (i <= FS_range);


  HIT_LIST_stop_timer(hit_list);
  return hit_list;
}

int FS_INDEX_search(FS_INDEX_t *FS_index, HIT_LIST_t *hit_list,
		    BIOSEQ *query, ULINT cutoff,
		    FS_INDEX_process_func *pfunc, 
		    FS_INDEX_range_convert_func *cfunc)
{
  FS_SEQ_t FS_query;
  FS_SEQ_t FS_neighbour;
  int dist;

  ULINT i;
  int FS_range;
  int S_cutoff = cutoff;
  int D_cutoff = cfunc(FS_index, query, cutoff);

  HIT_LIST_set_converted_range(hit_list, D_cutoff);
  /* Calculate FS_seq from seq */
  if (!FS_INDEX_verify_query(FS_index, query, &FS_query))
    return 0;

 /* Initialise variables */
  FS_K = FS_index->ptable->no_partitions;
  FS_KK[0] = 1;
  for (i = 1; i <= FS_index->frag_len; i++)
    FS_KK[i] = FS_KK[i-1] * FS_K;

  /* Calculate the maximum number of changes */
  FS_range = FS_INDEX_get_no_changes(FS_index, query, D_cutoff); 

  /* Process the bin associated with the query sequence */
  pfunc(FS_index, hit_list, query, FS_query, S_cutoff); 

  /* Process neighbouring bins */
  i = 1;
  do 
    {
      /* Initialize combination generator*/
      kSubsetLexInit(i);
      do 
	{
	  
	  if (!SCORE_MATRIX_verify_pos(FS_index->D, query, 
				       TT, i, D_cutoff))
	    continue;
	  FS_INDEX_iter_init(FS_index, FS_query);

	  while (FS_INDEX_iter_next(FS_index, query, FS_query, 
				    &dist, &FS_neighbour))
	    {
	      HIT_LIST_count_FS_seq_visited(hit_list, 1);
	      if (dist > D_cutoff)
		continue;
	      HIT_LIST_count_FS_seq_hit(hit_list, 1);
	      pfunc(FS_index, hit_list, query, FS_neighbour, 
		    S_cutoff);
	    }
	} 
      while (kSubsetLexSuccessor(i, FS_index->frag_len));
      i++;
    } 
  while (i <= FS_range);


  HIT_LIST_stop_timer(hit_list);
  return 1;
}


int FS_INDEX_kNN_search(FS_INDEX_t *FS_index, HIT_LIST_t *hit_list, 
			BIOSEQ *query, ULINT k, int V0,
			FS_INDEX_process_func *pfunc, 
			FS_INDEX_range_convert_func *cfunc)
{
  /* For now this works only for distance searches */
  int high = cfunc(FS_index, query, V0);
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
      HIT_LIST_reset(hit_list, query, FS_index->s_db,
		     FS_index->matrix, high);
      FS_INDEX_search(FS_index, hit_list, query, high,
		      pfunc, cfunc);
      size = HIT_LIST_get_seqs_hits(hit_list);
      if (size >= k)
	{
	  HIT_LIST_sort_incr(hit_list);
	  hit_k = HIT_LIST_get_hit(hit_list, k-1);
	  hit_largest = HIT_LIST_get_hit(hit_list, size-1);
	  d_k = hit_k->value;
	  d_largest = hit_largest->value;
	  if (d_k < d_largest)
	    {
	      HIT_LIST_reset(hit_list, query, FS_index->s_db,
			     FS_index->matrix, (int) d_k);
	      FS_INDEX_search(FS_index, hit_list, query, (int) d_k,
			      pfunc, cfunc);
	    }
	  return 1;
	}
      else
	  high += l;
    }
  while (1);
}


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
  FS_K = FS_index->ptable->no_partitions;
  FS_KK[0] = 1;
  for (i = 1; i <= FS_index->frag_len; i++)
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
	  
	  if (!SCORE_MATRIX_verify_pos(FS_index->D, query, 
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
      while (kSubsetLexSuccessor(i, FS_index->frag_len));
      i++;
    } 
  while (i <= FS_range);

  return 0;
}

HIT_LIST_t *SSCAN_QD_search(SEQUENCE_DB *s_db, const char *matrix, 
			    BIOSEQ *query, ULINT D_cutoff)
{
  ULINT one_percent_fragments;

  ULINT i, j;
  BIOSEQ *frag = mallocec(sizeof(BIOSEQ));


  ULINT frag_len = query->len;
  FS_PARTITION_t *ptable;
  SCORE_MATRIX_t *S;
  SCORE_MATRIX_t *D;
  int D_value;
  HIT_LIST_t *hit_list;
  ULINT no_frags;
  FS_SEQ_t FS_seq;
  ULINT len;


  no_frags = fastadb_count_Ffrags(s_db, frag_len);
  one_percent_fragments = no_frags / 100;
  fastadb_init_Ffrags(s_db, frag_len);

  /* Create partition table */
  ptable = FS_PARTITION_create("STAN#ILVM#KR#EDQ#WFYH#GPC", '#'); 

  /* Load matrices, create quasi-metric */
  S = SCORE_MATRIX_create(matrix, ptable); 
  D = SCORE_MATRIX_S_2_Dquasi(S); 

  hit_list = HIT_LIST_create(query, s_db, strdup(matrix), D_cutoff);
  frag_len = query->len;



  /* For each fragment */ 
  j = 0;
#if 0
  while (fastadb_get_next_Ffrag(s_db, frag_len, frag, &i))  
    {
      if (!BIOSEQ_2_FS_SEQ(frag, ptable, &FS_seq))
	continue;
      if (SCORE_MATRIX_evaluate_max(D, query, frag, 
				    D_cutoff, &D_value))
	{
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  HIT_LIST_insert_seq_hit(hit_list, frag, 
				  (float) D_value); 
	}
      HIT_LIST_count_seq_visited(hit_list, 1);
      /* Print progress bar */
      if (FS_INDEX_PRINT_BAR > 0)
	printbar(stdout, j+1, one_percent_fragments, 50);  
      j++;
   }
#endif
  i = 0;
  while (i < s_db->seq_data_len)  
    {
      len = strlen(s_db->seq_data + i);
      if (len < frag_len)
	{
	  i += len + 1;
	  continue;
	}
      
      frag->len = frag_len;
      frag->start = s_db->seq_data + i;
      i++;
      if (!BIOSEQ_2_FS_SEQ(frag, ptable, &FS_seq))
	continue;
      if (SCORE_MATRIX_evaluate_max(D, query, frag, 
				    D_cutoff, &D_value))
	{
	  HIT_LIST_count_seq_hit(hit_list, 1);
	  HIT_LIST_insert_seq_hit(hit_list, frag, 
				  (float) D_value); 
	}
      HIT_LIST_count_seq_visited(hit_list, 1);
      /* Print progress bar */
      if (FS_INDEX_PRINT_BAR > 0)
	printbar(stdout, j+1, one_percent_fragments, 50);  
      j++;
   }

  /* Take time */
  HIT_LIST_stop_timer(hit_list);

  /* Clean up */
  free(frag);

  /* Free matrices */
  SCORE_MATRIX_destroy(D);
  SCORE_MATRIX_destroy(S);

  /* Free partition table */
  FS_PARTITION_destroy(ptable);  

  return hit_list; 
}

int SSCAN_has_neighbour(SEQUENCE_DB *s_db,  SCORE_MATRIX_t *D, 
			FS_PARTITION_t *ptable,
			BIOSEQ *query, ULINT D_cutoff)
{
  ULINT i;
  BIOSEQ frag;


  ULINT frag_len = query->len;
  int D_value;
  ULINT no_frags;
  FS_SEQ_t FS_seq;
  ULINT len;


  no_frags = fastadb_count_Ffrags(s_db, frag_len);
  fastadb_init_Ffrags(s_db, frag_len);
  frag_len = query->len;



  /* For each fragment */ 
  i = 0;
  while (i < s_db->seq_data_len)  
    {
      len = strlen(s_db->seq_data + i);
      if (len < frag_len)
	{
	  i += len + 1;
	  continue;
	}
      
      frag.len = frag_len;
      frag.start = s_db->seq_data + i;
      i++;
      if (!BIOSEQ_2_FS_SEQ(&frag, ptable, &FS_seq))
	continue;
      if (SCORE_MATRIX_evaluate_max(D, query, &frag, 
				    D_cutoff, &D_value))
	{
#if 1
	  fprintf(stderr, "%d\n%.*s\n%.*s\n\n", D_value,
		  (int) query->len, query->start, 
		  (int) frag_len, frag.start);
#endif
	  return 1;
	}
   }
  return 0; 
}









/********************************************************************/
/********************************************************************/

/********************************************************************/    
/*                                                                  */
/*                     PROFILE BASED SEARCHES                       */ 
/*                                                                  */
/********************************************************************/    

/* FSINDEX general (static) functions and variables used here:

   FS_INDEX_verify_query();
   FS_INDEX_iter_init();
   kSubsetLexSuccessor();
   All static variables associated with FS generator.
*/

/* Profile specific functions:

   POS_MATRIX_get_no_changes() - may be moved to pmatrix module;
   FS_INDEX_profile_S_process_bin();
   FS_INDEX_profile_D_process_bin();
   FS_INDEX_profile_iter_next();
 */

static 
int POS_MATRIX_get_no_changes(POS_MATRIX *PD, int D_cutoff, 
			      ULINT Pfrom, ULINT Pto)
{
  int *D = mallocec((Pto - Pfrom + 1) * sizeof(int));
  ULINT i;
  int S = 0;

  for (i = 0; i <= Pto - Pfrom; i++)
    D[i] = PD->pMclosest[i+Pfrom];
  
  qsort(D, Pto - Pfrom + 1, sizeof(int), compare_int);

  for (i = 0; i <= Pto - Pfrom; i++)
    {
      S += D[i];
      if (S > D_cutoff)
	break;
    }
  free(D);
  return i;
}

void FS_INDEX_profile_S_process_bin(FS_INDEX_t *FS_index, 
				    HIT_LIST_t *hit_list,
				    POS_MATRIX *PS,
				    ULINT Pfrom, ULINT Pto,
				    FS_SEQ_t FS_query, 
				    int S_cutoff)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(FS_index->HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int S_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(FS_index->HT, 
					       FS_query, j);

      fastadb_get_Ffrag(FS_index->s_db, FS_index->frag_len, 
			&subject, frag_offset);

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


void FS_INDEX_profile_D_process_bin(FS_INDEX_t *FS_index, 
				    HIT_LIST_t *hit_list,
				    POS_MATRIX *PD,
				    ULINT Pfrom, ULINT Pto,
				    FS_SEQ_t FS_query,
				    int D_cutoff)
{
  ULINT n = FS_HASH_TABLE_get_no_seqs(FS_index->HT, FS_query);
  ULINT j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (j = 0; j < n; j++)
    {
      frag_offset = FS_HASH_TABLE_retrieve_seq(FS_index->HT, 
					       FS_query, j);

      fastadb_get_Ffrag(FS_index->s_db, FS_index->frag_len, 
			&subject, frag_offset);

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

static inline
int FS_INDEX_profile_iter_next(FS_INDEX_t *FS_index, POS_MATRIX *PD, 
			       ULINT Pfrom, ULINT Pto, 
			       FS_SEQ_t FS_query, int *dist,  
			       FS_SEQ_t *FS_neighbour) 
{
  ULINT i;
  int r = FS_r;
  int R;

  if (FS_r >= FS_Kl) return 0;
  *dist = 0;
  *FS_neighbour = FS_S0;
  for (i = 0; i < FS_l; i++)
    {
      R = ((r % (FS_K-1)) + FS_Q[i] + 1) % FS_K;
      *dist += PD->pM[PM_pM(Pfrom + TT[i], R)];
      *FS_neighbour += R * FS_KTT[i];
      r /= (FS_K-1);
    }
  FS_r++;
  return 1;
}






int FS_INDEX_profile_search(FS_INDEX_t *FS_index, 
			    HIT_LIST_t *hit_list, 
			    POS_MATRIX *PS, POS_MATRIX *PD, 
			    ULINT Pfrom, ULINT Pto, BIOSEQ *query,
			    ULINT cutoff, 
			    FS_INDEX_profile_process_func *pfunc, 
			    POS_MATRIX_range_convert_func *cfunc)
{
  FS_SEQ_t FS_query;
  FS_SEQ_t FS_neighbour;
  int dist;

  ULINT i;
  int FS_range;
  int S_cutoff = cutoff;
  int D_cutoff = cfunc(PS, query, cutoff, Pfrom, Pto);

  HIT_LIST_set_converted_range(hit_list, D_cutoff);
  /* Calculate FS_seq from seq */
  if (!FS_INDEX_verify_query(FS_index, query, &FS_query))
    return 0;

 /* Initialise variables */
  FS_K = FS_index->ptable->no_partitions;
  FS_KK[0] = 1;
  for (i = 1; i <= FS_index->frag_len; i++)
    FS_KK[i] = FS_KK[i-1] * FS_K;

  /* Calculate the maximum number of changes */
  FS_range =  POS_MATRIX_get_no_changes(PD, D_cutoff, Pfrom, Pto);

  /* Process the bin associated with the query sequence */
  pfunc(FS_index, hit_list, PS, Pfrom, Pto, FS_query, S_cutoff);

  /* Process neighbouring bins */
  i = 1;
  do 
    {
      /* Initialize combination generator*/
      kSubsetLexInit(i);
      do 
	{
	  if (!POS_MATRIX_verify_pos(PD, Pfrom, Pto, TT, i, D_cutoff))
	    continue;

	  FS_INDEX_iter_init(FS_index, FS_query);

	  while (FS_INDEX_profile_iter_next(FS_index, PD, Pfrom, Pto,  
		   FS_query, &dist, &FS_neighbour)) 
	    {
	      HIT_LIST_count_FS_seq_visited(hit_list, 1);
	      if (dist > D_cutoff)
		continue;
	      HIT_LIST_count_FS_seq_hit(hit_list, 1);
	      pfunc(FS_index, hit_list, PS, Pfrom, Pto, FS_neighbour,
		    S_cutoff); 
	    }
	} 
      while (kSubsetLexSuccessor(i, FS_index->frag_len));
      i++;
    } 
  while (i <= FS_range);


  HIT_LIST_stop_timer(hit_list);
  return 1;
}

/*********************************************************************/
/*********************************************************************/
