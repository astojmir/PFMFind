#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <assert.h>
#include "hit_list.h"
#include "keyword.h"
#include "partition.h"
#include "smatrix.h"
#include "FSindex.h"

#ifdef TRY_CLUSTERS
#include "clusters.h"
#endif

#define MERGE_DUPLICATES

/********************************************************************/ 
/*                                                                  */
/*                     FSINDX object                                */ 
/*                                                                  */
/********************************************************************/ 

#ifdef MERGE_DUPLICATES
static SEQUENCE_DB *s_db1;
static int m1; 

/* Fragment comparison routine */
static 
int comp_frags(const void *fg1, const void *fg2)
{
  SEQ_index_t *f1 = (SEQ_index_t *) fg1;
  SEQ_index_t *f2 = (SEQ_index_t *) fg2;
  BIOSEQ s1;
  BIOSEQ s2;

  fastadb_get_Ffrag(s_db1, m1, &s1, *f1);
  fastadb_get_Ffrag(s_db1, m1, &s2, *f2);
  return memcmp(s1.start, s2.start, m1);
}
#endif


/* Main constructor */
FSINDX *FS_INDEX_create(const char *database, ULINT flen,
			const char *abet, 
			const char sepchar, int skip)
 
{
  EXCEPTION *except;
  FSINDX *FSI;

  time_t time0 = time(NULL);
  time_t time1;
  double dt;
  ULINT one_percent_fragments;
  char *db_dir;
  char *db_base;

  ULINT i, j;
  ULINT no_frags;
  BIOSEQ *frag;
  FS_SEQ_t FS_seq;
  ULINT added_frags = 0;

  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  ULINT bs;

#ifdef MERGE_DUPLICATES
  int k;
#endif

  FSI = callocec(1, sizeof(FSINDX));

  Try {
    frag = mallocec(sizeof(BIOSEQ));

    /* Initialise FS_index */
    split_base_dir(database, &db_base, &db_dir);
    FSI->db_name = strdup(db_base);
    FSI->alphabet = strdup(abet);
    FSI->separator = sepchar;
    FSI->m = flen;
    FSI->KK = mallocec(FSI->m * sizeof(ULINT));
  
    /* Load database */
    fastadb_argt[0] = ACCESS_TYPE;
    fastadb_argt[1] = RETREIVE_DEFLINES;
    fastadb_argt[2] = NONE;
    fastadb_argv[0].access_type = MEMORY;
    fastadb_argv[1].retrieve_deflines = YES;
    FSI->s_db = fastadb_open(database, fastadb_argt, fastadb_argv); 
 
    /* Make fragment database */
    no_frags = fastadb_count_Ffrags(FSI->s_db, FSI->m);
    one_percent_fragments = (ULINT) (((double) no_frags/skip) / 100);
    fastadb_init_Ffrags(FSI->s_db, FSI->m);

    /* This is only necessary for Z-scores */
    fastadb_init_frags(FSI->s_db, FSI->m, FSI->m);

    /* Create partition table */
    FSI->ptable = FS_PARTITION_create(FSI->alphabet, FSI->separator); 

    /* Error if the sequence is too long for the partitioning*/
    FSI->K = FS_PARTITION_get_no_partitions(FSI->ptable);
    FSI->K_modtable = mallocec(2*FSI->K*sizeof(int));
    for (i = 0; i < FSI->K; i++)
      {
	FSI->K_modtable[i] = i;
	FSI->K_modtable[i + FSI->K] = i;
      }
    if ((FSI->K-1) >> ((sizeof(ULINT)*8)/FSI->m))
      Throw FSexcept(BAD_ARGS, "FS_INDEX_create(): Fragment length"
		     " too long for 32 bit word size!");

    /* Number of bins to allocate in hash table */
    FSI->no_bins = 1;
    for (i=0; i < FSI->m; i++)
      {
	FSI->KK[i] = FSI->no_bins;
	FSI->no_bins *= FSI->K;
      }

    /* Make initial bin size half as large as mean size */
    bs = (int)((float) no_frags / FSI->no_bins) + 1;

    /* Initialise hash table */
 
    FSI->no_seqs = 0;
    FSI->no_useqs = 0;
    FSI->shist_len = 5 * bs;
    FSI->shist = callocec(FSI->shist_len, sizeof(ULINT)); 
    FSI->uhist_len = 5 * bs;
    FSI->uhist = callocec(FSI->uhist_len, sizeof(ULINT)); 
    FSI->bin_size = callocec(FSI->no_bins, sizeof(ULINT));
    FSI->max_bin_size = callocec(FSI->no_bins, sizeof(ULINT)); 
    FSI->bin = mallocec(FSI->no_bins * sizeof(SEQ_index_t *));
    FSI->binL = 0;
    FSI->heap = NULL;

    FSI->u_size = callocec(FSI->no_bins, sizeof(ULINT));
    FSI->u = mallocec(FSI->no_bins * sizeof(int));
    FSI->u_heap = NULL;
    FSI->uL = 0;

    /* Two passes - first count fragments, then put them */

    printf("Counting fragments ...\n");

    /* Counting */ 
    j = 0;
    while (fastadb_get_next_Ffrag(FSI->s_db, FSI->m, frag, &i, skip))  
      {
	j++;

	/* Calculate its FS_seq */
	if (BIOSEQ_2_FS_SEQ(frag, FSI->ptable, &FS_seq))
	  {
	    /* Count the sequence */
	    FSI->max_bin_size[FS_seq]++;
	    FSI->no_seqs++;
	  }
	/* Print progress bar */
	printbar(stdout, j+1, one_percent_fragments, 50);  
      }  

    /* Allocating */
    printf("Allocating bins ...\n");

    for (j=0; j < FSI->no_bins; j++)
      FSI->bin[j] = mallocec(FSI->max_bin_size[j] 
			     * sizeof(SEQ_index_t));
    
    /* Commiting */
    printf("Collecting fragments ...\n");

    j = 0;
    fastadb_init_Ffrags(FSI->s_db, FSI->m);
    while (fastadb_get_next_Ffrag(FSI->s_db, FSI->m, frag, &i, skip))  
      {
	j++;
	/* Calculate its FS_seq */
	if (BIOSEQ_2_FS_SEQ(frag, FSI->ptable, &FS_seq))
	  {
	    /* Add to appropriate FS_bin */
	    added_frags++;
	    FSI->bin[FS_seq][FSI->bin_size[FS_seq]] = i;
	    FSI->bin_size[FS_seq]++;
	  }
	
	/* Print progress bar */
	printbar(stdout, j+1, one_percent_fragments, 50);  
      }  

#ifdef MERGE_DUPLICATES
    printf("Collecting duplicate fragments ...\n");

    for (j=0; j < FSI->no_bins; j++)
      {
	
	m1 = FSI->m;
	s_db1 = FSI->s_db;
	if (FSI->bin_size[j] != 0) {
	  
	  /* Sort fragments */
	  qsort(FSI->bin[j], FSI->bin_size[j], sizeof(SEQ_index_t),
		comp_frags);
	  /* Count unique fragments */
	  FSI->u_size[j] = 1;
	  for (i=1; i < FSI->bin_size[j]; i++)
	    if (comp_frags(FSI->bin[j]+i, FSI->bin[j]+i-1) != 0)
	      FSI->u_size[j]++;

	  FSI->no_useqs += FSI->u_size[j];
	  /* Allocate unique fragments */
	  FSI->u[j] = callocec(FSI->u_size[j], sizeof(int));
	  
	  /* Fill unique array */
	  
	  FSI->u[j][0] = 0;
	  k = 0;
	  for (i=1; i < FSI->bin_size[j]; i++)
	    if (comp_frags(FSI->bin[j]+i, FSI->bin[j]+i-1) != 0)
	      {
		k++;
		FSI->u[j][k] = i; 
	      }
	}
	printbar(stdout, j+1, FSI->no_bins/100, 50);  
      
#endif
	
	/* Add statistics */
	{
	  ULINT s0; /* Old length of frequency vector */
	  
	  if (FSI->bin_size[j] > FSI->bin_size[FSI->binL])
	    FSI->binL = j;
	  
	  if (FSI->bin_size[j] >= FSI->shist_len)
	    {
	      s0 = FSI->shist_len;
	      FSI->shist_len = FSI->bin_size[j]+1;
	      FSI->shist = reallocec(FSI->shist, 
				     FSI->shist_len * sizeof(ULINT)); 
	      memset(FSI->shist + s0, 0, 
		     (FSI->shist_len - s0) * sizeof(ULINT));
	    }
	  FSI->shist[FSI->bin_size[j]]++;
	  
#ifdef MERGE_DUPLICATES
	  if (FSI->u_size[j] > FSI->u_size[FSI->uL])
	    FSI->uL = j;
	  
	  if (FSI->u_size[j] >= FSI->uhist_len)
	    {
	      s0 = FSI->uhist_len;
	      FSI->uhist_len = FSI->u_size[j]+1;
	      FSI->uhist = reallocec(FSI->uhist, FSI->uhist_len * sizeof(ULINT));  
	      memset(FSI->uhist + s0, 0, (FSI->uhist_len - s0) * sizeof(ULINT));
	    }
	  FSI->uhist[FSI->u_size[j]]++;
#endif
      }
    }


    /* Take time */
    time1 = time(NULL);
    dt = difftime(time1, time0)/60.0;
    FSI->db_no_frags = no_frags;
    FSI->dtime = dt;

    /* Clean up */
    free(frag);
    free(db_dir);
    free(db_base);
  
    /* We have an index ready for search */
    FSI->index_name = strdup("New Index");
  }
  Catch (except) {
    FS_INDEX_destroy(FSI);
    Throw except;
  }
  return FSI; 
}


/* Destructor */
void FS_INDEX_destroy(FSINDX *FSI)
{
  /* Free all bins */
  ULINT i;

  if (FSI->heap == NULL)
    for (i=0; i < FSI->no_bins; i++)
      free(FSI->bin[i]);
  else
    free(FSI->heap);

  /* TO DO: Free unique sequneces */

  free(FSI->bin);
  free(FSI->bin_size);
  free(FSI->shist);
  if (FSI->max_bin_size != NULL)
    free(FSI->max_bin_size);


  /* Free partition table */
  if (FSI->ptable != NULL)
    FS_PARTITION_destroy(FSI->ptable);  

  /* Free sequence database */
  if (FSI->s_db != NULL)
    fastadb_close(FSI->s_db);

  /* Free names */
  free(FSI->db_name);
  free(FSI->alphabet);

  free(FSI);
}

/* Load, Save */
int FS_INDEX_save(FSINDX *FSI, const char *filename)
{ 
  FILE *stream = fopen(filename, "wb");
  int len;
  ULINT i;
 
  /* Open file */
  if(stream == NULL)
    Throw FSexcept(FOPEN_ERR, 
		   "FS_INDEX_save(): Could not open the file %s.",
		   filename);

  len = strlen(FSI->db_name) + 1;
  fwrite(&len, sizeof(int), 1, stream);
  fwrite(FSI->db_name, sizeof(char), len, stream);

  len = strlen(FSI->alphabet) + 1;
  fwrite(&len, sizeof(int), 1, stream);
  fwrite(FSI->alphabet, sizeof(char), len, stream);

  fwrite(&(FSI->separator), sizeof(char), 1, stream);
  fwrite(&(FSI->m), sizeof(int), 1, stream);
  fwrite(&(FSI->db_no_frags), sizeof(ULINT), 1, stream);

  FS_PARTITION_write(FSI->ptable, stream); 

  fwrite(&(FSI->no_bins), sizeof(ULINT), 1, stream);
  fwrite(&(FSI->no_seqs), sizeof(ULINT), 1, stream);
  fwrite(&(FSI->shist_len), sizeof(ULINT), 1, stream);
  fwrite(&(FSI->binL), sizeof(SEQ_index_t), 1, stream);
  fwrite(FSI->shist, sizeof(ULINT), FSI->shist_len, stream); 
  fwrite(FSI->bin_size, sizeof(ULINT), FSI->no_bins, stream);  

  for (i=0; i < FSI->no_bins; i++)
    fwrite(FSI->bin[i], sizeof(SEQ_index_t), FSI->bin_size[i], stream);

  fwrite(&(FSI->no_useqs), sizeof(ULINT), 1, stream);
  fwrite(&(FSI->uhist_len), sizeof(ULINT), 1, stream);
  fwrite(&(FSI->uL), sizeof(SEQ_index_t), 1, stream);
  fwrite(FSI->uhist, sizeof(ULINT), FSI->uhist_len, stream); 
  fwrite(FSI->u_size, sizeof(ULINT), FSI->no_bins, stream);  
  
  for (i=0; i < FSI->no_bins; i++)
    fwrite(FSI->u[i], sizeof(int), FSI->u_size[i], stream);
  fwrite(&(FSI->dtime), sizeof(double), 1, stream);

  fclose(stream);
  return 1;
}

FSINDX *FS_INDEX_load(const char *filename)
{
  EXCEPTION *except;
  FSINDX *FSI;
  FILE *stream = fopen(filename, "rb");
  int len;
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  char *basename;
  char *dirname;
  char *full_dbname;
  int i;
  SEQ_index_t *current;
  int *t2;
 


  FSI = callocec(1, sizeof(FSINDX));
  Try {

    /* Open file */
    if(stream == NULL)
      Throw FSexcept(FOPEN_ERR, 
		     "FS_INDEX_load(): Could not open the file %s.",
		     filename);
    FSI->index_name = strdup(filename);

    fread(&len, sizeof(int), 1, stream);
    FSI->db_name = mallocec(len);
    fread(FSI->db_name, sizeof(char), len, stream);

    split_base_dir(filename, &basename, &dirname);
    cat_base_dir(&full_dbname, FSI->db_name, dirname);
    
    fread(&len, sizeof(int), 1, stream);
    FSI->alphabet = mallocec(len);
    fread(FSI->alphabet, sizeof(char), len, stream);

    fread(&FSI->separator, sizeof(char), 1, stream);
    fread(&FSI->m, sizeof(int), 1, stream);
    fread(&FSI->db_no_frags, sizeof(ULINT), 1, stream);

    /* Load database */
    fastadb_argt[0] = ACCESS_TYPE;
    fastadb_argt[1] = RETREIVE_DEFLINES;
    fastadb_argt[2] = NONE;
    fastadb_argv[0].access_type = MEMORY;
    fastadb_argv[1].retrieve_deflines = YES;
    
    FSI->s_db = fastadb_open(full_dbname, fastadb_argt, fastadb_argv); 
      
    /* This is only necessary for Z-scores */
    fastadb_init_frags(FSI->s_db, FSI->m, FSI->m);

    FSI->ptable = FS_PARTITION_read(stream); 
    
    fread(&(FSI->no_bins), sizeof(ULINT), 1, stream);
    FSI->bin_size = mallocec(FSI->no_bins * sizeof(ULINT));
    FSI->u_size = mallocec(FSI->no_bins * sizeof(ULINT));
    FSI->bin = mallocec(FSI->no_bins * sizeof(SEQ_index_t *));
    FSI->u = mallocec(FSI->no_bins * sizeof(int *));

    fread(&(FSI->no_seqs), sizeof(ULINT), 1, stream);
    FSI->heap = mallocec(FSI->no_seqs * sizeof(SEQ_index_t)); 

    fread(&(FSI->shist_len), sizeof(ULINT), 1, stream);
    FSI->shist = mallocec(FSI->shist_len * sizeof(ULINT));

    fread(&(FSI->binL), sizeof(SEQ_index_t), 1, stream);

    FSI->max_bin_size = NULL;

    fread(FSI->shist, sizeof(ULINT), FSI->shist_len, stream); 
    fread(FSI->bin_size, sizeof(ULINT), FSI->no_bins, stream);  

    current =  FSI->heap;
    for (i=0; i < FSI->no_bins; i++)
      {
	FSI->bin[i] = current;
	fread(FSI->bin[i], sizeof(SEQ_index_t), FSI->bin_size[i],
	      stream);
	current += FSI->bin_size[i];
      }

    fread(&(FSI->no_useqs), sizeof(ULINT), 1, stream);
    FSI->u_heap = mallocec(FSI->no_useqs * sizeof(int)); 

    fread(&(FSI->uhist_len), sizeof(ULINT), 1, stream);
    FSI->uhist = mallocec(FSI->uhist_len * sizeof(ULINT)); 

    fread(&(FSI->uL), sizeof(SEQ_index_t), 1, stream);

    fread(FSI->uhist, sizeof(ULINT), FSI->uhist_len, stream); 
    fread(FSI->u_size, sizeof(ULINT), FSI->no_bins, stream);  
  
    t2 = FSI->u_heap;
    for (i=0; i < FSI->no_bins; i++)
      {
	FSI->u[i] = t2;
	fread(FSI->u[i], sizeof(int), FSI->u_size[i], stream);
	t2 += FSI->u_size[i];
      }
    fread(&(FSI->dtime), sizeof(double), 1, stream);
    
    /* Allocate arrays and calculate FS_KK */
    FSI->KK = mallocec(FSI->m * sizeof(ULINT));
    FSI->K = FS_PARTITION_get_no_partitions(FSI->ptable);
    FSI->K_modtable = mallocec(2*FSI->K*sizeof(int));
    for (i = 0; i < FSI->K; i++)
      {
	FSI->K_modtable[i] = i;
	FSI->K_modtable[i + FSI->K] = i;
      }
    FSI->KK[0] = 1;
    for (i=1; i < FSI->m; i++)
      FSI->KK[i] = FSI->K * FSI->KK[i-1];
    
    fclose(stream);
    free(basename);
    free(dirname);
  }
  Catch(except) {
    FS_INDEX_destroy(FSI);
    Throw except;
  }
  return FSI;
}

/* Check query */

static
FS_SEQ_t FSSRCH_check_query(FSINDX *FSI, BIOSEQ *query)
{
  FS_SEQ_t FS_query;

  /* Check length of query */
  if (query->len != FSI->m)
    Throw FSexcept(BAD_ARGS, "FSSRCH_check_query():"
		   " Could not convert query to FS fragment\n"
		   " Query: %*s\n", (int) query->len, 
		   query->start);

  /* Calculate FS_seq from seq */
  if (!BIOSEQ_2_FS_SEQ(query, FSI->ptable, &FS_query))
    Throw FSexcept(BAD_ARGS, "FSSRCH_check_query():"
		   " Could not convert query"
		   " to FS fragment\n Query: %*s\n", 
		   (int) query->len, query->start);

  /* Also check maximum size of bit-field vs (K-1)m */
  if (sizeof(long long)*8 < ((FSI->K - 1) * FSI->m))
    Throw FSexcept(BAD_ARGS, "FSSRCH_check_query():"
		   " Too small bit field\n");

  return FS_query;
}


/* Print */
void FS_INDEX_print_stats(FSINDX *FSI, FILE *stream, int options) 
{
  ULINT i;
  ULINT CF;
  ULINT CS;
  
  fprintf(stream, "\n*** FS_INDEX Statistics ***\n");
  fprintf(stream, "Database: %s\n", FSI->db_name);
  fprintf(stream, "Database size: %ld\n", FSI->s_db->length); 
  fprintf(stream, "Number of sequences: %ld\n\n", FSI->s_db->no_seq);  
  fprintf(stream, "Partitions: %s\n", FSI->alphabet);
  fprintf(stream, "Fragment Length: %d\n", FSI->m);
  fprintf(stream, "Total Fragments in database: %ld\n", 
	  FSI->db_no_frags); 
  fprintf(stream, "Creation Time: %.2f mins\n\n", FSI->dtime);

 
  fprintf(stream, "Total number of fragments in index: %ld\n", 
	  FSI->no_seqs);
  fprintf(stream, "Total number of index entries: %ld\n", 
	  FSI->no_bins);
  fprintf(stream, "Average size of index entry: %.2f\n", 
	  (float) FSI->no_seqs / (float) FSI->no_bins);
  fprintf(stream, "Largest index entry: %s\n",
	  FS_seq_print(FSI->binL, FSI->ptable, FSI->m));
  fprintf(stream, "Size of largest index entry: %ld\n",
       	  FSI->bin_size[FSI->binL]);
  fprintf(stream, "Total number of distinct fragments in index: %ld\n", 
	  FSI->no_useqs);
  fprintf(stream, "Average 'distinct' size of index entry: %.2f\n", 
	  (float) FSI->no_useqs / (float) FSI->no_bins);
  fprintf(stream, "Largest 'distinct' index entry: %s\n",
	  FS_seq_print(FSI->uL, FSI->ptable, FSI->m));
  fprintf(stream, "Size of largest 'distinct' index entry: %ld\n",
	  FSI->u_size[FSI->uL]);

  if (options & 1) {
    fprintf(stream, "\n* Distribution of numbers of fragments per"
	    " index entry *\n");
    fprintf(stream, "%17s %10s\n","Index entry size", "Frequency");
    CF = 0;
    CS = 0;
    for (i = 0; i < FSI->shist_len; i++)
      {
	CF += FSI->shist[i];
	CS += FSI->shist[i] * i;    
	if (FSI->shist[i] > 0)
	  fprintf(stream, "%17ld %10ld %10ld %10ld\n", i, FSI->shist[i], 
		  CF, CS);
      }
    fprintf(stream, "\n");
  }

  if (options & 2) {
    fprintf(stream, "\n* Distribution of numbers of distinct"
	    " fragments per index entry *\n");
    fprintf(stream, "%17s %10s\n","Index entry size", "Frequency");
    CF = 0;
    CS = 0;
    for (i = 0; i < FSI->uhist_len; i++)
      {
	CF += FSI->uhist[i];
	CS += FSI->uhist[i] * i;
	if (FSI->uhist[i] > 0)
	  fprintf(stream, "%17ld %10ld %10ld %10ld\n", i, FSI->uhist[i], 
		  CF, CS);
      }
    fprintf(stream, "\n");
  }
}

void FS_INDEX_print_bin(FSINDX *FSI, BIOSEQ *query, FILE *stream,
			int options)
{
  int j;
  ULINT id;
  ULINT from;
  char *defline;
  ULINT frag_offset;
  BIOSEQ subject;

  FS_SEQ_t FS_query = FSSRCH_check_query(FSI, query);
  char *s = FS_seq_print(FS_query, FSI->ptable, FSI->m);

  fprintf(stream, "Bin: %s (%ld)\n", s, FS_query);
  fprintf(stream, "Size: %ld\n", FSI->bin_size[FS_query]); 
  fprintf(stream, "'Unique' size: %ld\n",
	  FSI->u_size[FS_query]);
  for (j=0; j < FSI->bin_size[FS_query]; j++)
    {
      frag_offset = FSI->bin[FS_query][j];
      fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, frag_offset);
      fastadb_find_Ffrag_seq(FSI->s_db, &subject, &id, &from);
      defline = FSI->s_db->seq[id].id.defline;
      if (options)
	fprintf(stream, "%8ld %5ld ", id, from+1);
      fprintf(stream, "%*.*s ", FSI->m, FSI->m, subject.start);
      if (options)
	fprintf(stream, "%5ld %.*s", from+FSI->m , 
		58-FSI->m, defline);
      fprintf(stream, "\n");
    }

  free(s);

}


/********************************************************************/ 
/*      Access functions                                            */
/********************************************************************/ 

ULINT FS_INDEX_get_bin_size(FSINDX *FSI, ULINT bin)
{
  if (bin >= FSI->no_bins)
    Throw FSexcept(INDEX_OUT_OF_RANGE,
		   "FS_INDEX_get_bin_size(): Index out of range.");
  else
    return FSI->bin_size[bin];
}

ULINT FS_INDEX_get_unique_bin_size(FSINDX *FSI, ULINT bin)
{
  if (bin >= FSI->no_bins)
    Throw FSexcept(INDEX_OUT_OF_RANGE,
		   "FS_INDEX_get_bin_size(): Index out of range.");
  else
    return FSI->u_size[bin];
}

BIOSEQ *FS_INDEX_get_seq(FSINDX *FSI, ULINT bin, ULINT seq)
{
  BIOSEQ *subject = mallocec(sizeof(BIOSEQ));
  ULINT frag_offset;
  ULINT id;
  ULINT from;

  if (bin >= FSI->no_bins || seq >= FSI->bin_size[bin])
    Throw FSexcept(INDEX_OUT_OF_RANGE,
		   "FS_INDEX_get_bin_size(): Index out of range.");

  frag_offset = FSI->bin[bin][seq];
  fastadb_get_Ffrag(FSI->s_db, FSI->m, subject, frag_offset);
  fastadb_find_Ffrag_seq(FSI->s_db, subject, &id, &from); 
  subject->id.defline = FSI->s_db->seq[id].id.defline;

  return subject;
}

/********************************************************************/ 
/*                                                                  */
/*                     FSSRCH object                                */ 
/*                                                                  */
/********************************************************************/ 

/* We only accept distance matrices. Conversion to be done outside.  */


#define BITSET(a, b) ((a) |= (1LL << (b)))
#define BITTEST(a, b) ((a) & (1LL << (b)))


typedef struct
{
  int pos;                 /* Position within fragment              */
  int letter;              /* Letter changed to                     */
  int dist;                /* Distance of the tranformation         */
  FS_SEQ_t binval;         /* Bin value                             */  
} FSTRANSF;


typedef struct
{
  void *M;                 /* Scoring matrix, normal or profile     */
  BIOSEQ *query;           /* Query sequence                        */
  HIT_LIST_t *HL;          /* Search results                        */
  FSINDX *FSI;             /* Index                                 */
  int *eps;                /* Cutoff value (pointer)                */  
  int *kNN;                /* kNN value (pointer)                   */  
} PFUNC_ARGS;

typedef void FSINDX_scan_func(PFUNC_ARGS *args, FS_SEQ_t bin);

typedef struct
{
  int m;                   /* Fragment length                       */
  FS_SEQ_t qbin;           /* Full query bin                        */
  FS_SEQ_t *qbval;         /* Query bin values by position          */  
  FSINDX_scan_func *pfunc; /* Bin scan function pter                */
  PFUNC_ARGS *args;        /* Arguments to the pfunc                */
  int eps;                 /* Cutoff value                          */
  int kNN;                 /* Number of nearest neighbours          */
  int Tn;                  /* Size of T array                       */
  FSTRANSF *T;             /* Transformations ordered by priority   */
  long long *UT;           /* Allowed transformations by position   */
  int p;                   /* Number of permutations for Z-score    */
  FPARAMS *fp;             /* Filtering parameters                  */
} FSSRCH; 






/********************************************************************/ 
/*      Auxillary functions                                         */
/********************************************************************/ 



/*  FSTRANSF comparison routines */
static 
int FSTRANSF_comp_incr(const void *tr1, const void *tr2)
{
  FSTRANSF *t1 = (FSTRANSF *) tr1;
  FSTRANSF *t2 = (FSTRANSF *) tr2;
  return (t1->dist - t2->dist);
}

static 
int FSTRANSF_comp_decr(const void *tr1, const void *tr2)
{
  FSTRANSF *t1 = (FSTRANSF *) tr1;
  FSTRANSF *t2 = (FSTRANSF *) tr2;
  return (t2->dist - t1->dist);
}

/********************************************************************/ 
/*      Processing functions                                        */
/********************************************************************/ 
static
void process_bin(PFUNC_ARGS *args, FS_SEQ_t bin0)
{
  FSINDX *FSI = args->FSI;  
  ULINT n = FSI->u_size[bin0];
  ULINT N = FSI->bin_size[bin0];
  
  int i;
  int j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (i = 0; i < n; i++)
    {
      j = FSI->u[bin0][i];
      frag_offset = FSI->bin[bin0][j];
      fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, frag_offset);

      if (SCORE_MATRIX_evaluate_max(args->M, args->query, &subject, 
				    *(args->eps), &D_value))
	{
	  HIT_LIST_count_useq_hit(args->HL, 1);
	  HIT_LIST_insert_seq_hit(args->HL, &subject, (float) D_value); 
	  HIT_LIST_count_seq_hit(args->HL, 1);
	  j++;
	  while ((j < N && i == n-1) || j < FSI->u[bin0][i+1])
	    {
	      HIT_LIST_count_seq_hit(args->HL, 1);
	      frag_offset = FSI->bin[bin0][j];
	      fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, 
				frag_offset);
	      HIT_LIST_insert_seq_hit(args->HL, &subject, 
				      (float) D_value); 
	      j++;
	    }

	}
      HIT_LIST_count_seq_visited(args->HL,
	i == n-1 ? N - FSI->u[bin0][i] :
                   FSI->u[bin0][i+1]- FSI->u[bin0][i]);
      HIT_LIST_count_useq_visited(args->HL, 1);
    }
}

static
void kNN_process_bin(PFUNC_ARGS *args, FS_SEQ_t bin0)
{
  FSINDX *FSI = args->FSI;  
  ULINT n = FSI->u_size[bin0];
  ULINT N = FSI->bin_size[bin0];
  
  int i;
  int j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (i = 0; i < n; i++)
    {
      j = FSI->u[bin0][i];
      frag_offset = FSI->bin[bin0][j];
      fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, frag_offset);

      if (SCORE_MATRIX_evaluate_max(args->M, args->query, &subject, 
				    *(args->eps), &D_value) ||
	  HIT_LIST_get_seqs_hits(args->HL) < *(args->kNN))
	{
	  HIT_LIST_count_useq_hit(args->HL, 1);
	  *(args->eps) = 
	    HIT_LIST_insert_seq_hit_queue(args->HL, &subject, 
					  (float) D_value); 	  
	  HIT_LIST_count_seq_hit(args->HL, 1);
	  j++;
	  while ((j < N && i == n-1) || j < FSI->u[bin0][i+1])
	    {
	      HIT_LIST_count_seq_hit(args->HL, 1);
	      frag_offset = FSI->bin[bin0][j];
	      fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, 
				frag_offset);
	      *(args->eps) = 
		HIT_LIST_insert_seq_hit_queue(args->HL, &subject, 
					      (float) D_value); 	  
	      j++;
	    }

	}
      HIT_LIST_count_seq_visited(args->HL,
	i == n-1 ? N - FSI->u[bin0][i] :
                   FSI->u[bin0][i+1]- FSI->u[bin0][i]);
      HIT_LIST_count_useq_visited(args->HL, 1);
    }
}

/*      if (POS_MATRIX_evaluate_min(PS, &subject, Pfrom, Pto, 
	S_cutoff, &S_value)) */



static
void profile_process_bin(PFUNC_ARGS *args, FS_SEQ_t bin0)
{
  FSINDX *FSI = args->FSI;  
  ULINT n = FSI->u_size[bin0];
  ULINT N = FSI->bin_size[bin0];
  
  int i;
  int j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (i = 0; i < n; i++)
    {
      j = FSI->u[bin0][i];
      frag_offset = FSI->bin[bin0][j];
      fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, frag_offset);

      if (POS_MATRIX_evaluate_max(args->M, &subject, 0, FSI->m-1,
				  *(args->eps), &D_value)) 
	{
	  HIT_LIST_count_useq_hit(args->HL, 1);
	  HIT_LIST_insert_seq_hit(args->HL, &subject, (float) D_value); 
	  HIT_LIST_count_seq_hit(args->HL, 1);
	  j++;
	  while ((j < N && i == n-1) || j < FSI->u[bin0][i+1])
	    {
	      HIT_LIST_count_seq_hit(args->HL, 1);
	      frag_offset = FSI->bin[bin0][j];
	      fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, 
				frag_offset);
	      HIT_LIST_insert_seq_hit(args->HL, &subject, 
				      (float) D_value); 
	      j++;
	    }

	}
      HIT_LIST_count_seq_visited(args->HL,
	i == n-1 ? N - FSI->u[bin0][i] :
                   FSI->u[bin0][i+1]- FSI->u[bin0][i]);
      HIT_LIST_count_useq_visited(args->HL, 1);
    }
}



static
void profile_kNN_process_bin(PFUNC_ARGS *args, FS_SEQ_t bin0)
{
  FSINDX *FSI = args->FSI;  
  ULINT n = FSI->u_size[bin0];
  ULINT N = FSI->bin_size[bin0];
  
  int i;
  int j;
  ULINT frag_offset;
  BIOSEQ subject;
  int D_value;

  for (i = 0; i < n; i++)
    {
      j = FSI->u[bin0][i];
      frag_offset = FSI->bin[bin0][j];
      fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, frag_offset);

      if (POS_MATRIX_evaluate_max(args->M, &subject, 0, FSI->m-1,
				  *(args->eps), &D_value) ||
	  HIT_LIST_get_seqs_hits(args->HL) < *(args->kNN))
	{
	  HIT_LIST_count_useq_hit(args->HL, 1);
	  *(args->eps) = 
	    HIT_LIST_insert_seq_hit_queue(args->HL, &subject, 
					  (float) D_value); 	  
	  HIT_LIST_count_seq_hit(args->HL, 1);
	  j++;
	  while ((j < N && i == n-1) || j < FSI->u[bin0][i+1])
	    {
	      HIT_LIST_count_seq_hit(args->HL, 1);
	      frag_offset = FSI->bin[bin0][j];
	      fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, 
				frag_offset);
	      *(args->eps) = 
		HIT_LIST_insert_seq_hit_queue(args->HL, &subject, 
					      (float) D_value); 	  
	      j++;
	    }

	}
      HIT_LIST_count_seq_visited(args->HL,
	i == n-1 ? N - FSI->u[bin0][i] :
                   FSI->u[bin0][i+1]- FSI->u[bin0][i]);
      HIT_LIST_count_useq_visited(args->HL, 1);
    }
}





/********************************************************************/ 
/*         FSSRCH creation / cleanup functions                      */
/********************************************************************/ 
static
void FSSRCH_clean(FSSRCH *FSS)
{
  free(FSS->qbval);
  free(FSS->T);
  free(FSS->UT);
  free(FSS->args);
  free(FSS);
}

static
FSSRCH *FSSRCH_init(FSINDX *FSI, BIOSEQ *query, SCORE_MATRIX_t *D,
		    int d0, int kNN, HIT_LIST_t *HL)
{
  EXCEPTION *except;
  FSSRCH *FSS = callocec(1, sizeof(FSSRCH));
  char *q = query->start;

  FS_SEQ_t FS_query = 0;
  int i;
  int j;
  int k;
  int R;

  /* Check correctness of query */
  Try {
    FSS->qbin = FSSRCH_check_query(FSI, query);

    FSS->args = mallocec(sizeof(PFUNC_ARGS));
    FSS->args->M = (void *) D;
    FSS->args->query = query;

    if (HL == NULL)  
      FSS->args->HL = HIT_LIST_create(query, FSI->s_db, 
			 SCORE_MATRIX_filename(D), D, 
			 smatrix_eval, d0);
    else
      HIT_LIST_reset(HL, query, FSI->s_db, SCORE_MATRIX_filename(D),
		     D, smatrix_eval, d0);
		  
    FSS->args->FSI = FSI;
    FSS->args->eps = &FSS->eps;
    FSS->args->kNN = &FSS->kNN;
 
    FSS->m = query->len;
    FSS->qbval = mallocec(FSS->m * sizeof(FS_SEQ_t));
    for (i=0; i < FSS->m; i++)
      {
	j = FS_PARTITION_get_pttn(FSI->ptable, q[i]);
	FSS->qbval[i] = j * FSI->KK[i];
	FS_query += FSS->qbval[i];
      }

    if (kNN <= 0)
      FSS->pfunc = process_bin;
    else
      {
	FSS->pfunc = kNN_process_bin;
	HIT_LIST_set_kNN(FSS->args->HL, kNN); 
      }
    FSS->eps = d0;
    FSS->kNN = kNN;
    
    /* All transformations */
    FSS->Tn = (FSI->K - 1) * FSI->m;
    FSS->T = mallocec(FSS->Tn * sizeof(FSTRANSF));
    FSS->UT = callocec(FSI->m, sizeof(long long));
    
    for (i=0, k=0; i < FSI->m; i++)
      for (j=0; j < FSI->K; j++)
	if (j != FS_PARTITION_get_pttn(FSI->ptable, q[i]))
	  {
	    FSS->T[k].pos = i;
	    FSS->T[k].letter = j; 
	    R = FS_PARTITION_get_poffset(FSI->ptable, j);
	    FSS->T[k].dist = D->pM[q[i] & A_SIZE_MASK][R]; 
	    FSS->T[k].binval = j * FSI->KK[i];
	    k++;
	  }
    /* Sort */
    if (kNN <= 0)
      qsort(FSS->T, FSS->Tn, sizeof(FSTRANSF), FSTRANSF_comp_decr);
    else
      qsort(FSS->T, FSS->Tn, sizeof(FSTRANSF), FSTRANSF_comp_incr);
    /* Flag transformations by position */
    
    for (k=0; k < FSS->Tn; k++)
      {
	BITSET(FSS->UT[FSS->T[k].pos], k);
      }  

    /* Set some info for hit list */
    HIT_LIST_set_converted_range(FSS->args->HL, FSS->eps);
    HIT_LIST_set_index_data(FSS->args->HL, FSI->index_name, 
			    FSI->alphabet, FSI->no_bins,
			    FSI->no_seqs);
  }
  Catch (except) {
    FSSRCH_clean(FSS);
    Throw except;
  }

  return FSS;
}

static
FSSRCH *FSSRCH_profile_init(FSINDX *FSI, POS_MATRIX *PD,
		    int d0, int kNN, HIT_LIST_t *HL)
{
  EXCEPTION *except;
  FSSRCH *FSS = callocec(1, sizeof(FSSRCH));
  char *q = PD->query->start;

  FS_SEQ_t FS_query = 0;
  int i;
  int j;
  int k;
  int R;

  /* Check correctness of query */
  Try {
    FSS->qbin = FSSRCH_check_query(FSI, PD->query);

    FSS->args = mallocec(sizeof(PFUNC_ARGS));
    FSS->args->M = PD;
    FSS->args->query = PD->query;
    if (HL == NULL)  
      FSS->args->HL = 
	HIT_LIST_create(PD->query, FSI->s_db, 
			POS_MATRIX_filename(PD), PD, 
			profile_eval, d0);
    else
      HIT_LIST_reset(HL, PD->query, FSI->s_db, 
		     POS_MATRIX_filename(PD), 
		     PD, profile_eval, d0);

    FSS->args->FSI = FSI;
    FSS->args->eps = &FSS->eps;
    FSS->args->kNN = &FSS->kNN;
    
    FSS->m = PD->query->len;
    FSS->qbval = mallocec(FSS->m * sizeof(FS_SEQ_t));
    for (i=0; i < FSS->m; i++)
      {
	j = FS_PARTITION_get_pttn(FSI->ptable, q[i]);
	FSS->qbval[i] = j * FSI->KK[i];
	FS_query += FSS->qbval[i];
      }
    
    if (kNN <= 0)
      FSS->pfunc = profile_process_bin;
    else
      {
	FSS->pfunc = profile_kNN_process_bin;
	HIT_LIST_set_kNN(FSS->args->HL, kNN); 
      }
    FSS->eps = d0;
    FSS->kNN = kNN;

    /* All transformations */
    FSS->Tn = (FSI->K - 1) * FSI->m;
    FSS->T = mallocec(FSS->Tn * sizeof(FSTRANSF));
    FSS->UT = callocec(FSI->m, sizeof(long long));

    for (i=0, k=0; i < FSI->m; i++)
      for (j=0; j < FSI->K; j++)
	if (j != FS_PARTITION_get_pttn(FSI->ptable, q[i]))
	{
	  FSS->T[k].pos = i;
	  FSS->T[k].letter = j; 
	  R = FS_PARTITION_get_poffset(FSI->ptable, j);
	  FSS->T[k].dist = PD->pM[PM_pM(i, R)]; 
	  FSS->T[k].binval = j * FSI->KK[i];
	  k++;
	}
  /* Sort */
    if (kNN <= 0)
      qsort(FSS->T, FSS->Tn, sizeof(FSTRANSF), FSTRANSF_comp_incr);
    else
      qsort(FSS->T, FSS->Tn, sizeof(FSTRANSF), FSTRANSF_comp_decr);
    
    /* Flag transformations by position */
    for (k=0; k < FSS->Tn; k++)
    BITSET(FSS->UT[FSS->T[k].pos], k);
    

    /* Set some info for hit list */
    HIT_LIST_set_converted_range(FSS->args->HL, FSS->eps);
    HIT_LIST_set_index_data(FSS->args->HL, FSI->index_name, 
			    FSI->alphabet, FSI->no_bins,
			    FSI->no_seqs);
  }
  Catch (except) {
    FSSRCH_clean(FSS);
    Throw except;
  }

  return FSS;
}


/********************************************************************/ 
/*      Recursive tree traversal function                           */
/********************************************************************/ 

static
void check_bins(FSSRCH *FSS, int dist, FS_SEQ_t cbin, long long AT,
		int k0)
{
  int k;
  int dist1;
  FS_SEQ_t bin1;
  long long AT1;

#if 0
  for (k=k0; k < FSS->Tn; k++)
    if (BITTEST(AT, k))
      printf("%d ", FSS->T[k].dist);
  putchar('\n');
#endif

  for (k=k0; k < FSS->Tn; k++)
    if (BITTEST(AT, k))
    {
      HIT_LIST_count_FS_seq_visited(FSS->args->HL, 1);
      if ((dist1 = dist + FSS->T[k].dist) <= FSS->eps)
	{
	  bin1 = cbin + FSS->T[k].binval - FSS->qbval[FSS->T[k].pos]; 
	  HIT_LIST_count_FS_seq_hit(FSS->args->HL, 1);
	  FSS->pfunc(FSS->args, bin1); 
	  AT1 = AT & ~FSS->UT[FSS->T[k].pos];
	  if (AT1 > 0)
	    check_bins(FSS, dist1, bin1, AT1,k+1);	  
	}
    }



  return;
} 

/********************************************************************/ 
/*      Search functions                                            */
/********************************************************************/ 

static
HIT_LIST_t *FSSRCH_search(FSSRCH *FSS)
{
  long long AT = 0; /* Allowed transformations */
  int k;

  /* Set all relevant AT bits to 1 */
  for (k=0; k < FSS->Tn; k++)
    BITSET(AT, k);

  /* Process the bin associated with the query sequence */
  HIT_LIST_count_FS_seq_visited(FSS->args->HL, 1);
  HIT_LIST_count_FS_seq_hit(FSS->args->HL, 1);
  FSS->pfunc(FSS->args, FSS->qbin); 

  /* Process neighbouring bins - recursive */
  check_bins(FSS, 0, FSS->qbin, AT, 0); 

  /* Set off timer */
  HIT_LIST_stop_timer(FSS->args->HL);

  /* Post-processing */

  /* Pull out kNN sequences */
  if (FSS->kNN > 0)
    HIT_LIST_sort_kNN(FSS->args->HL);

  /* Convert scores */

  /* Cleanup */
  FSSRCH_clean(FSS);
  return FSS->args->HL;
}


HIT_LIST_t *FSINDX_rng_srch(FSINDX *FSI, BIOSEQ *query, SCORE_MATRIX_t *D,
		     int d0, HIT_LIST_t *HL)  
{
  return FSSRCH_search(FSSRCH_init(FSI, query, D, d0, -1, HL));
}

HIT_LIST_t *FSINDX_kNN_srch(FSINDX *FSI, BIOSEQ *query, SCORE_MATRIX_t *D,
		     int kNN, HIT_LIST_t *HL)
{
  return FSSRCH_search(FSSRCH_init(FSI, query, D, INT_MAX, kNN, HL));
}

HIT_LIST_t *FSINDX_prof_rng_srch(FSINDX *FSI, BIOSEQ *query, 
				 SCORE_MATRIX_t *D, int s0, 
				 double lambda, double A,
				 const char *freq_filename,
				 int iters, int s1,
				 HIT_LIST_t *HL)  
{
  int i;
  int d0;
  POS_MATRIX *PS = SCORE_2_POS_MATRIX(D, query);
 
  POS_MATRIX_init(PS, lambda, POS_MATRIX_simple_pseudo_counts,
		  freq_filename);
  POS_MATRIX_simple_pseudo_counts_init(PS, A);
	  
   d0 = PS->qS - s0;
   HL = FSSRCH_search(FSSRCH_profile_init(FSI, PS, d0, -1, HL));
   POS_MATRIX_convert(PS, HL);
   HIT_LIST_sort_decr(HL);
   HIT_LIST_print(HL, stdout, 0); 

  i = iters;
  while (i--)
    {
     
      /* Filter here */

      HIT_LIST_get_hit_seqs(HL, &PS->seq, s0, 
			    &PS->no_seqs, &PS->max_no_seqs);
      if (PS->no_seqs == 0)
	break;
      POS_MATRIX_Henikoff_weights(PS);
      POS_MATRIX_update(PS);
      
      d0 = PS->qS - s1;
      HL = FSSRCH_search(FSSRCH_profile_init(FSI, PS, d0, -1, HL));
      POS_MATRIX_convert(PS, HL);
      HIT_LIST_sort_decr(HL);
      HIT_LIST_print(HL, stdout, 0); 
    }
  return HL;
}

HIT_LIST_t *FSINDX_prof_kNN_srch(FSINDX *FSI, BIOSEQ *query, 
				 SCORE_MATRIX_t *D, int kNN, 
				 double A,
				 const char *freq_filename,
				 int iters, 
				 HIT_LIST_t *HL)
{
  int i;

  POS_MATRIX *PS = SCORE_2_POS_MATRIX(D, query);
 
  POS_MATRIX_init(PS, 1.0, POS_MATRIX_simple_pseudo_counts,
		  freq_filename);
  POS_MATRIX_simple_pseudo_counts_init(PS, A);
	  
  HL = FSSRCH_search(FSSRCH_profile_init(FSI, PS, INT_MAX, kNN, HL));

  POS_MATRIX_convert(PS, HL);
  HIT_LIST_print(HL, stdout, 0); 

  i = iters;
  while (i--)
    {
      /* Filter here */

      HIT_LIST_get_hit_seqs(HL, &PS->seq, -100000, 
			    &PS->no_seqs, &PS->max_no_seqs);
      if (PS->no_seqs == 0)
	break;
      POS_MATRIX_Henikoff_weights(PS);
      POS_MATRIX_update(PS);

      HL = FSSRCH_search(FSSRCH_profile_init(FSI, PS, INT_MAX, kNN, HL));

      POS_MATRIX_convert(PS, HL);
      HIT_LIST_print(HL, stdout, 0); 
    }
  return HL;
}
