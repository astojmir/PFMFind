#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "sarray.h"
#include "FSindex.h"
#include "misclib.h"

#define MERGE_DUPLICATES
/********************************************************************/ 
/*                                                                  */
/*                     FSSRCH object                                */ 
/*                                                                  */
/********************************************************************/ 

#define MAX_HITS 100

#define PQDATUM HIT 
#define PQPRIO(p) (p.dist) 

#include "pqueue.h"

/********************************************************************/ 
/*                                                                  */
/*                     FSINDX object                                */ 
/*                                                                  */
/********************************************************************/ 

/* Allocate search specific arrays */
static
void alloc_FSSRCH(FSSRCH *FSS, int m, int *nbpt)
{
  int i;

  FSS->nbpt_closest = mallocec(m * sizeof(int));
  FSS->Mlen = m * 3 / 2;
  FSS->M = callocec(FSS->Mlen, sizeof(int *));
  FSS->cdist = callocec(FSS->Mlen+1, sizeof(int));
  FSS->cdist++;
  FSS->nbpt_diff = mallocec(m * sizeof(long *));
  FSS->nbpt_dist = mallocec(m * sizeof(int *));

  for (i=0; i < FSS->Mlen; i++) 
    FSS->M[i] = mallocec(A_SIZE * sizeof(int));
  for (i=0; i < m; i++) {
    FSS->nbpt_diff[i] = mallocec(nbpt[i] * sizeof(long));
    FSS->nbpt_dist[i] = mallocec(nbpt[i] * sizeof(int));
  }
  FSS->p_queue = pqinit(NULL, MAX_HITS);
  FSS->hits_size = MAX_HITS;
  FSS->hits = mallocec(MAX_HITS * sizeof(HIT));
}

static char *_base_;
static char *_end_;
static int _len_; 
static int **_rank_;

/* Fragment comparison routine */
static 
int comp_frags(const void *fg1, const void *fg2)
{
  int i;
  int **r = _rank_;
  const uint32_t *f1 = fg1;
  const uint32_t *f2 = fg2;
  char *s1 = _base_ + *f1;
  char *s2 = _base_ + *f2;
  for (i=0; i < _len_; i++, r++, s1++, s2++) {
    if ((*r)[*s1] < (*r)[*s2])
      return -1;
    else if ((*r)[*s1] > (*r)[*s2])
      return 1;
  }
  return 0;
}

static inline 
int get_unique_size(FSINDX *FSI, uint32_t len, FS_SEQ_t j) 
{
  ULINT k;
  ULINT c;
  if (FSI->bins[j+1] - FSI->bins[j]) {
    c = 1;
    for (k=FSI->bins[j]+1; k < FSI->bins[j+1]; k++)
      if (FSI->lcpb[k] < len) c++;
  }
  else c = 0;
  return c;
}

static 
uint32_t FS_INDEX_get_distinct_hist(FSINDX *FSI, uint32_t **h,
				 uint32_t *size, uint32_t len,
				 uint32_t *uL)
{
  ULINT j;
  ULINT c;
  ULINT M = 0;
  ULINT N = 0;
  ULINT L = 0;

  *h = callocec(*size, sizeof(uint32_t)); 
  for (j=0; j < FSI->no_bins; j++) {
    c = get_unique_size(FSI, len, j); 
    N += c;
    if (M < c) {
      M = c;
      L = j;
    }
    (*h)[c]++;
  }
  *size = M+1;
  *h = reallocec(*h, M*sizeof(uint32_t));
  if(uL) *uL = L;
  return N;
}

static
unsigned char *compute_lcpb(int *a, const char *base, uint32_t n)
{
  unsigned char *lcpb = mallocec(n);
  int i;
  lcpb[0] = 0;
  for (i=1; i < n; i++) {
    int h = 0;
    const char *s = base + a[i-1];
    const char *t = base + a[i];
    while (*s++ == *t++) h++;
    lcpb[i] = min(h, 255);
  }
  return lcpb;
}

/* Main constructor */
FSINDX *FS_INDEX_create(const char *database, uint32_t len, const char **sepn, int use_sa,
			int print_flag)
 {
  EXCEPTION *except;
  FSINDX *FSI = NULL;

  FS_TABLE *ptable;

  time_t time0 = time(NULL);
  time_t time1;
  double dt;
  char *db_dir = NULL;
  char *db_base = NULL;

  ULINT i, j;
  char *s;
  FS_SEQ_t FS_seq;

  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  uint32_t bs = 0;

  uint32_t no_bins = 0;
  int **tmp_bin = NULL;
  uint32_t *tmp_size = NULL;
  int *tmp_lcpx = NULL;

  Try {
    /* Create partition table */
    ptable = FS_TABLE_init(sepn, len);
    if (ptable == NULL)
      Throw FSexcept(BAD_ARGS, "FS_INDEX_create():"
		     " Could not create partition table\n");

    /* Initialise FS_index */
    FSI = callocec(1, sizeof(FSINDX));
    FSI->ptable = ptable;
    split_base_dir(database, &db_base, &db_dir);
    FSI->db_name = strdup(db_base);
    FSI->index_name = strdup("New Index");
    FSI->m = len;
    FSI->no_bins = no_bins = ptable->bins;
    FSI->use_sa = use_sa;

    /* Load database */
    fastadb_argt[0] = ACCESS_TYPE;
    fastadb_argt[1] = RETREIVE_DEFLINES;
    fastadb_argt[2] = NONE;
    fastadb_argv[0].access_type = MEMORY;
    fastadb_argv[1].retrieve_deflines = YES;
    FSI->s_db = fastadb_open(database, fastadb_argt, fastadb_argv); 
 
    /* Global variables */
    _base_ = fastadb_data_pter(FSI->s_db, 0);  
    _end_ = fastadb_end_heap(FSI->s_db);
    _len_ = len; 
    _rank_ = ptable->rank;
    
    tmp_size = callocec(FSI->no_bins+1, sizeof(uint32_t));
    tmp_size++;
    tmp_size[-1] = 0;

    if (use_sa) {
      *_end_ = '_';
      FSI->sa_size = _end_ - _base_ + 1;
      if (print_flag) {
	printf("***** Suffix Array *****\n");
	printf("Converting sequences ...\n");
      }
      FSI->sa = FS_TABLE_order_convert_heap(FSI->ptable, _base_,
				       FSI->sa_size);
      if (print_flag) {
	printf("Creating suffix array ...\n");
      }
      if (sarray(FSI->sa, FSI->sa_size) == -1)
	Throw FSexcept(RUNTIME_ERR, "FS_INDEX_create():"
		       " Could not create suffix array.\n");
      if (print_flag) {
	printf("Creating lcp array ...\n");
      }
      tmp_lcpx = lcp(FSI->sa, _base_, FSI->sa_size); 
      if (!tmp_lcpx)
	Throw FSexcept(RUNTIME_ERR, "FS_INDEX_create():"
		       " Could not create lcp array.\n");
      FSI->lcpx = mallocec(FSI->sa_size);
      for (i=0; i < FSI->sa_size; i++)
	FSI->lcpx[i] = min(tmp_lcpx[i], 255);
      free(tmp_lcpx);
      *_end_='\0';
    }

    if (print_flag) {
      printf("***** FSindex bins *****\n");
      printf("Counting fragments ...\n");
    }
    for (s=_base_; s <= _end_; s++) {
      FS_seq = FS_SEQ(ptable, s, len); 
      if (FS_seq != -1) {
	tmp_size[FS_seq]++;
	FSI->no_seqs++;
      }
    }  
    FSI->db_no_frags = FSI->no_seqs;
    if (print_flag) {
      printf("Allocating bins ...\n");
    }    
    tmp_bin = callocec(no_bins, sizeof(uint32_t *));    
    for (j=0; j < FSI->no_bins; j++)
      if (tmp_size[j]) {
	tmp_bin[j] = callocec(tmp_size[j], sizeof(uint32_t));
	tmp_size[j] = 0;
      }

    if (print_flag) {
      printf("Collecting fragments ...\n");
    }
    for (s=_base_; s <= _end_; s++) {
      FS_seq = FS_SEQ(ptable, s, len);
      if (FS_seq != -1) {
	tmp_bin[FS_seq][tmp_size[FS_seq]]
	  = s - _base_;
	tmp_size[FS_seq]++;
      }
    }

    if (print_flag) {
      printf("Sorting each bin ...\n");
    }
    FSI->oa = mallocec(FSI->no_seqs * sizeof(int));
    for (j=0; j < no_bins; j++) {
      if (tmp_size[j]) {
	qsort(tmp_bin[j], tmp_size[j], sizeof(int), comp_frags);
	if (bs < tmp_size[j]) {
	  FSI->binL = j;
	  bs = tmp_size[j];
	}
      }
      memcpy(FSI->oa+tmp_size[j-1], tmp_bin[j], tmp_size[j]*sizeof(int));
      tmp_size[j] += tmp_size[j-1];
      if (tmp_bin[j]) {
	free(tmp_bin[j]);
	tmp_bin[j] = NULL;
      }
    }
    free(tmp_bin);
    tmp_bin = NULL;
    tmp_size--;
    FSI->bins = tmp_size;

    if (print_flag) {
      printf("Computing lcp array ...\n");
    }

    FSI->lcpb = compute_lcpb(FSI->oa, _base_, FSI->no_seqs);

    /* Statistics */
    /* Allocate histograms */
    FSI->shist_len = bs+1;
    FSI->shist = callocec(FSI->shist_len, sizeof(uint32_t)); 

    /* Get bin size stats */
    if (print_flag) {
      printf("Collecting bin size stats ...\n");
    }
    for (j=0; j < no_bins; j++) 
      FSI->shist[tmp_size[j+1]-tmp_size[j]]++;

    /* Get histogram of sizes of bins in terms
       of distinct fragments */
    FSI->uhist_len = FSI->shist_len;
    FSI->no_useqs = FS_INDEX_get_distinct_hist(FSI, &FSI->uhist,
					       &FSI->uhist_len,
					       FSI->m, &FSI->uL);
    /* Take time */
    time1 = time(NULL);
    dt = difftime(time1, time0)/60.0;
    FSI->dtime = dt;

    /* Allocate search specific variables */
    for (i=0; i < THREADS; i++)
      alloc_FSSRCH(FSI->FSS+i, FSI->m, FSI->ptable->nbpt);
  }
  Catch (except) {
    if (tmp_bin) { 
      for (j=0; j < no_bins; j++)
	free(tmp_bin[j]);
      free(tmp_bin);
    }
    free(tmp_lcpx);
    if (tmp_size) {
      tmp_size--;
      free(tmp_size);
    }
    FS_INDEX_destroy(FSI);
    Throw except;
  }
  free(db_dir);
  free(db_base);
  return FSI; 
}

/* Destructor */
void FS_INDEX_destroy(FSINDX *FSI)
{
  /* printf("FS_INDEX_destroy run.\n"); */
  /* Free all bins */
  ULINT i, j;
  if (FSI == NULL) return;

  /* Free FSS structures */
  for (i=0; i < THREADS; i++) {
    if (FSI->FSS[i].M != NULL) {
      for (j=0; j < FSI->FSS[i].Mlen; j++)
	free(FSI->FSS[i].M[j]);
      for (j=0; j < FSI->m; j++) {
	free(FSI->FSS[i].nbpt_diff[j]);
	free(FSI->FSS[i].nbpt_dist[j]);
      }
      free(FSI->FSS[i].M);
    }
    if (FSI->FSS[i].cdist) 
      FSI->FSS[i].cdist--;
    free(FSI->FSS[i].cdist);
    free(FSI->FSS[i].nbpt_diff);
    free(FSI->FSS[i].nbpt_dist);
    free(FSI->FSS[i].nbpt_closest);
    free(FSI->FSS[i].hits);
    if (FSI->FSS[i].p_queue) {
      free(FSI->FSS[i].p_queue->d);
      free(FSI->FSS[i].p_queue);
    }
  }

  free(FSI->db_name);
  free(FSI->index_name);
  if (FSI->s_db != NULL)
    fastadb_close(FSI->s_db);
  FS_TABLE_del(FSI->ptable);

  free(FSI->shist);
  free(FSI->uhist);

  if (FSI->sa) {
    free(FSI->sa);
    free(FSI->lcpx);
  }
  free(FSI->bins);
  free(FSI->oa);
  free(FSI->lcpb);

  free(FSI);
}

/* Load, Save */
int FS_INDEX_save(FSINDX *FSI, const char *filename)
{ 
  FILE *fp = fopen(filename, "wb");

  if (fp == NULL)
    Throw FSexcept(FOPEN_ERR, 
		   "FS_INDEX_save(): Could not open the file %s.",
		   filename);
  
  fwrite_string(FSI->db_name, fp);
  fwrite_string(FSI->index_name, fp);
  FS_TABLE_write(FSI->ptable, fp); 

  fwrite(&(FSI->m), sizeof(int), 1, fp);
  fwrite(&(FSI->dtime), sizeof(double), 1, fp);
  fwrite(&(FSI->db_no_frags), sizeof(uint32_t), 1, fp);
  fwrite(&(FSI->no_bins), sizeof(uint32_t), 1, fp);

  fwrite(&(FSI->no_seqs), sizeof(uint32_t), 1, fp);
  fwrite(&(FSI->use_sa), sizeof(int), 1, fp);

  fwrite(&(FSI->binL), sizeof(SEQ_index_t), 1, fp);
  fwrite(&(FSI->uL), sizeof(SEQ_index_t), 1, fp);

  fwrite(&(FSI->shist_len), sizeof(uint32_t), 1, fp);
  fwrite(FSI->shist, sizeof(uint32_t), FSI->shist_len, fp); 
  fwrite(&(FSI->uhist_len), sizeof(uint32_t), 1, fp);
  fwrite(FSI->uhist, sizeof(uint32_t), FSI->uhist_len, fp); 

  fwrite(FSI->bins, sizeof(uint32_t), FSI->no_bins+1, fp);  
  fwrite(FSI->oa, sizeof(int), FSI->no_seqs, fp);
  fwrite(FSI->lcpb, sizeof(char), FSI->no_seqs, fp);

  if (FSI->use_sa) {
    fwrite(&(FSI->sa_size), sizeof(int), 1, fp);
    fwrite(FSI->sa, sizeof(int), FSI->sa_size, fp); 
    fwrite(FSI->lcpx, sizeof(char), FSI->sa_size, fp); 
  }

  fclose(fp);
  return 1;
}

FSINDX *FS_INDEX_load(const char *filename)
{
  EXCEPTION *except;
  FSINDX *FSI = callocec(1, sizeof(FSINDX));

  FILE *fp = fopen(filename, "rb");
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  char *basename;
  char *dirname;
  char *full_dbname;
  int i;

  fastadb_argt[0] = ACCESS_TYPE;
  fastadb_argt[1] = RETREIVE_DEFLINES;
  fastadb_argt[2] = NONE;
  fastadb_argv[0].access_type = MEMORY;
  fastadb_argv[1].retrieve_deflines = YES;

  Try {
    if(fp == NULL)
      Throw FSexcept(FOPEN_ERR, 
		     "FS_INDEX_load(): Could not open the file %s.",
		     filename);
    fread_string(&(FSI->db_name), fp);
    fread_string(&(FSI->index_name), fp);
    FSI->ptable = FS_TABLE_read(fp);
    if (FSI->ptable == NULL)
      Throw FSexcept(BAD_ARGS, 
		     "FS_INDEX_load(): Could not init ptable.");

    split_base_dir(filename, &basename, &dirname);
    cat_base_dir(&full_dbname, FSI->db_name, dirname);
    FSI->s_db = fastadb_open(full_dbname, fastadb_argt, fastadb_argv); 

    fread(&FSI->m, sizeof(int), 1, fp);
    fread(&(FSI->dtime), sizeof(double), 1, fp);
    fread(&FSI->db_no_frags, sizeof(uint32_t), 1, fp);
    fread(&(FSI->no_bins), sizeof(uint32_t), 1, fp);

    fread(&(FSI->no_seqs), sizeof(uint32_t), 1, fp);
    fread(&(FSI->use_sa), sizeof(int), 1, fp);

    fread(&(FSI->binL), sizeof(SEQ_index_t), 1, fp);
    fread(&(FSI->uL), sizeof(SEQ_index_t), 1, fp);

    fread(&(FSI->shist_len), sizeof(uint32_t), 1, fp);
    FSI->shist = mallocec(FSI->shist_len * sizeof(uint32_t));
    fread(FSI->shist, sizeof(uint32_t), FSI->shist_len, fp); 
    fread(&(FSI->uhist_len), sizeof(uint32_t), 1, fp);
    FSI->uhist = mallocec(FSI->uhist_len * sizeof(uint32_t)); 
    fread(FSI->uhist, sizeof(uint32_t), FSI->uhist_len, fp); 

    FSI->bins = mallocec((FSI->no_bins+1) * sizeof(uint32_t));
    FSI->oa = mallocec(FSI->no_seqs * sizeof(uint32_t));
    FSI->lcpb = mallocec(FSI->no_seqs * sizeof(char));
    fread(FSI->bins, sizeof(uint32_t), FSI->no_bins+1, fp);
    fread(FSI->oa, sizeof(int), FSI->no_seqs, fp);
    fread(FSI->lcpb, sizeof(char), FSI->no_seqs, fp);

    if (FSI->use_sa) {
      fread(&(FSI->sa_size), sizeof(int), 1, fp);
      FSI->sa = mallocec(FSI->sa_size * sizeof(int)); 
      FSI->lcpx = mallocec(FSI->sa_size); 
      fread(FSI->sa, sizeof(int), FSI->sa_size, fp); 
      fread(FSI->lcpx, sizeof(char), FSI->sa_size, fp); 
    }

    fclose(fp);
    free(basename);
    free(dirname);
    free(full_dbname);

    /* Allocate search specific variables */
    for (i=0; i < THREADS; i++)
      alloc_FSSRCH(FSI->FSS+i, FSI->m, FSI->ptable->nbpt);
  }
  Catch(except) {
    FS_INDEX_destroy(FSI);
    Throw except;
  }
  return FSI;
}


/* Print */
char *FS_INDEX_sprint_stats(FSINDX *FSI, int options) 
{
  char *buf = NULL;
  int size = 0;
  int len = 0;
  ULINT i;
  ULINT CF;
  ULINT CS;
  char *s;

  absprintf(&buf, &size, &len, "\n*** FS_INDEX Statistics ***\n");
  absprintf(&buf, &size, &len, "Database: %s\n", FSI->db_name);
  absprintf(&buf, &size, &len, "Database size: %ld\n", FSI->s_db->length); 
  absprintf(&buf, &size, &len, "Number of sequences: %ld\n\n", FSI->s_db->no_seq);  
  absprintf(&buf, &size, &len, "Partitions:\n");
  s = FS_TABLE_sprint(FSI->ptable);
  absprintf(&buf, &size, &len, "%s", s);
  free(s);
  absprintf(&buf, &size, &len, "Fragment Length: %d\n", FSI->m);
  absprintf(&buf, &size, &len, "Total Fragments in database: %ld\n", 
	  FSI->db_no_frags); 
  absprintf(&buf, &size, &len, "Creation Time: %.2f mins\n\n", FSI->dtime);

 
  absprintf(&buf, &size, &len, "Total number of fragments in index: %ld\n", 
	  FSI->no_seqs);
  absprintf(&buf, &size, &len, "Total number of index entries: %ld\n", 
	  FSI->no_bins);
  absprintf(&buf, &size, &len, "Average size of index entry: %.2f\n", 
	  (float) FSI->no_seqs / (float) FSI->no_bins);
  s = FS_SEQ_print(FSI->ptable, FSI->binL);
  absprintf(&buf, &size, &len, "Largest index entry: %s\n", s);
  free(s);
  absprintf(&buf, &size, &len, "Size of largest index entry: %ld\n",
       	  FSI->bins[FSI->binL+1] - FSI->bins[FSI->binL]);

  absprintf(&buf, &size, &len, "Total number of distinct fragments in index: %ld\n", 
	    FSI->no_useqs);
  absprintf(&buf, &size, &len, "Average 'distinct' size of index entry: %.2f\n", 
	    (float) FSI->no_useqs / (float) FSI->no_bins);
  s = FS_SEQ_print(FSI->ptable, FSI->uL);
  absprintf(&buf, &size, &len, "Largest 'distinct' index entry: %s\n", s);
  free(s);
  absprintf(&buf, &size, &len, "Size of largest 'distinct' index entry: %ld\n",
	    get_unique_size(FSI, FSI->m, FSI->uL));

  if (options & 1) {
    absprintf(&buf, &size, &len, "\n* Distribution of numbers of fragments per"
	    " index entry *\n");
    absprintf(&buf, &size, &len, "%17s %10s\n","Index entry size", "Frequency");
    CF = 0;
    CS = 0;
    for (i = 0; i < FSI->shist_len; i++) {
      CF += FSI->shist[i];
      CS += FSI->shist[i] * i;    
      if (FSI->shist[i] > 0)
	absprintf(&buf, &size, &len, "%17ld %10ld %10ld %10ld\n", i, FSI->shist[i], 
		CF, CS);
    }
    absprintf(&buf, &size, &len, "\n");
  }

  if (options & 2) {
    absprintf(&buf, &size, &len, "\n* Distribution of numbers of distinct"
	    " fragments per index entry *\n");
    absprintf(&buf, &size, &len, "%17s %10s\n","Index entry size", "Frequency");
    CF = 0;
    CS = 0;
    for (i = 0; i < FSI->uhist_len; i++) {
      CF += FSI->uhist[i];
      CS += FSI->uhist[i] * i;
      if (FSI->uhist[i] > 0)
	absprintf(&buf, &size, &len, "%17ld %10ld %10ld %10ld\n", i, FSI->uhist[i], 
		CF, CS);
    }
    absprintf(&buf, &size, &len, "\n");
  }
  return buf;
}

void FS_INDEX_fprint_stats(FSINDX *FSI, FILE *fp, int options) 
{
  char *buf = FS_INDEX_sprint_stats(FSI, options);
  fprintf(fp, "%s", buf);
  free(buf);
}

char *FS_INDEX_sprint_bin(FSINDX *FSI, FS_SEQ_t bin, int options)
{
  char *buf = NULL;
  int size = 0;
  int len = 0;
  int j;
  ULINT id;
  ULINT from;
  const char *defline;
  BIOSEQ subject;

  char *s;

  if (bin >= FSI->no_bins)
    return NULL;
  s = FS_SEQ_print(FSI->ptable, bin);

  absprintf(&buf, &size, &len, "Bin: %s (%ld)\n", s, bin);
  absprintf(&buf, &size, &len, "Size: %ld\n", 
	    FSI->bins[bin+1] - FSI->bins[bin]); 
  absprintf(&buf, &size, &len, "'Unique' size: %ld\n",
	    get_unique_size(FSI, FSI->m, bin)); 

  for (j=FSI->bins[bin]; j < FSI->bins[bin+1]; j++) {
    subject.start = fastadb_data_pter(FSI->s_db, FSI->oa[j]);
    subject.len = min(FSI->m, strlen(subject.start));
    fastadb_find_Ffrag_seq(FSI->s_db, &subject, &id, &from);
    defline = FSI->s_db->seq[id].id.defline;
    if (options)
      absprintf(&buf, &size, &len, "%8ld %5ld ", id, from+1);
    absprintf(&buf, &size, &len, "%*.*s ", (int) subject.len, (int) subject.len, 
	    subject.start);
    if (options)
      absprintf(&buf, &size, &len, "%5ld %.*s", from+subject.len, 
	      58-(int)subject.len, defline);
    absprintf(&buf, &size, &len, "\n");
  }
  free(s);
  return buf;
}

void FS_INDEX_fprint_bin(FSINDX *FSI, FS_SEQ_t bin, FILE *fp,
			int options)
{
  char *buf = FS_INDEX_sprint_bin(FSI, bin, options);
  fprintf(fp, "%s", buf);
  free(buf);
}


/********************************************************************/ 
/*      Insertion functions                                         */
/********************************************************************/ 

static 
void FSSRCH_insert(FSSRCH *FSS, uint32_t offset, int dist) 
{
  if (FSS->hits_len >= FSS->hits_size) {
    FSS->hits_size += MAX_HITS;
    FSS->hits = reallocec(FSS->hits, FSS->hits_size * sizeof(HIT));
  }
  FSS->hits[FSS->hits_len].offset = offset;
  FSS->hits[FSS->hits_len++].dist = dist;
  FSS->no_hits++;
}

static
void FSSRCH_insert_queue(FSSRCH *FSS, uint32_t offset, int dist) 
{
  struct pqueue *p_queue = FSS->p_queue;
  HIT top;
  HIT temp;
  HIT new_hit;

  new_hit.offset = offset;
  new_hit.dist = dist;

  if (FSS->no_hits < FSS->kNN) {
    pqinsert(p_queue, new_hit);
    if (++FSS->no_hits == FSS->kNN) {
      pqpeek(p_queue, &top);
      FSS->eps = top.dist;
    }
    else
      FSS->eps = INT_MAX;
  }
  else if (dist < FSS->eps) {
    pqremove(p_queue, &temp);
    pqinsert(p_queue, new_hit);
    pqpeek(p_queue, &top);
    FSS->eps = top.dist;

    if (temp.dist == top.dist) {
      if (FSS->hits_len >= FSS->hits_size) {
	FSS->hits_size += FSS->kNN;
	FSS->hits = reallocec(FSS->hits, FSS->hits_size*sizeof(HIT));
      }
      FSS->hits[FSS->hits_len++] = temp;
      FSS->no_hits++; 
    }
    else {
      FSS->hits_len = 0;
      FSS->no_hits = FSS->kNN; 
    }
  }
  else { /* if (dist == FSS->eps) */
    if (FSS->hits_len >= FSS->hits_size) {
      FSS->hits_size += FSS->kNN;
      FSS->hits = reallocec(FSS->hits, FSS->hits_size*sizeof(HIT));
    }

    FSS->hits[FSS->hits_len++] = new_hit;
    FSS->no_hits++;
  }
}

/********************************************************************/ 
/*      Recursive tree traversal functions                          */
/********************************************************************/ 
static
void check_bins(FSSRCH *FSS, FS_SEQ_t cbin, int dist, int i)
{
  int j;
  int k;
  int dist1;
  FS_SEQ_t bin1;

  //  for (k=i; k < FSS->len; k++) {
  for (k=FSS->len-1; k >= i; k--) {
    if (dist + FSS->nbpt_closest[k] > FSS->eps) continue;      
    for (j=0; j < FSS->nbpt[k]; j++) {
      HIT_LIST_count_FS_seq_visited(FSS->HL, 1);
      bin1 = cbin + FSS->nbpt_diff[k][j];
      dist1 = dist + FSS->nbpt_dist[k][j];
      if (dist1 > FSS->eps) continue;
      HIT_LIST_count_FS_seq_hit(FSS->HL, 1);
      FSS->pfunc(FSS, bin1);
      check_bins(FSS, bin1, dist1, k+1);
    }
  }
}

/********************************************************************/ 
/*      Processing functions                                        */
/********************************************************************/ 

static inline
int eval_dist(int **M, char *s, int qlen)
{
  int d=0;
  int i;
  for (i=0; i < qlen; i++, M++, s++)
    d += (*M)[*s & A_SIZE_MASK];
  return d;

}


/* Scan all without taking account of suffix array structure */
static
void process_bin_all(FSSRCH *FSS, FS_SEQ_t bin)
{
  ULINT n;
  char *base;
  int **M;
  int qlen;
  int dist;
  int i;
  int *a;

  n = FSS->FSI->bins[bin+FSS->bincr] - FSS->FSI->bins[bin];
  if (n == 0) return;
  a = FSS->FSI->oa + FSS->FSI->bins[bin]; 
  M = FSS->M;
  qlen = FSS->qlen;
  base = FSS->base;
  for (i=0; i < n; i++, a++) {
    dist = eval_dist(M, base + *a, qlen);
    if (dist > FSS->eps) continue;
    FSS->ifunc(FSS, *a, dist);
    HIT_LIST_count_seq_hit(FSS->HL, 1);
    HIT_LIST_count_useq_hit(FSS->HL, qlen);
  }
  HIT_LIST_count_useq_visited(FSS->HL, n * qlen);
  HIT_LIST_count_seq_visited(FSS->HL, n);
}


static
void process_bin_dups(FSSRCH *FSS, FS_SEQ_t bin)
{
  ULINT n;
  char *base;
  int **M;
  int qlen;
  int dist;
  int i;

  int *a;
  uint8_t *lcpx;


  n = FSS->FSI->bins[bin+FSS->bincr] - FSS->FSI->bins[bin];
  if (n == 0) return;
  a = FSS->FSI->oa + FSS->FSI->bins[bin]; 
  M = FSS->M;
  qlen = FSS->qlen;
  base = FSS->base;
  /* Do first sequence in any case */
  dist = eval_dist(M, base + *a, qlen);
  HIT_LIST_count_useq_visited(FSS->HL, qlen);
  if (dist <= FSS->eps) {
    FSS->ifunc(FSS, *a, dist);
    HIT_LIST_count_useq_hit(FSS->HL, qlen);
    HIT_LIST_count_seq_hit(FSS->HL, 1);
  }
  a++;
  lcpx = FSS->FSI->lcpb + FSS->FSI->bins[bin] + 1;
  for (i=1; i < n; i++, a++, lcpx++) {
    if (qlen > *lcpx) {
      dist = eval_dist(M, base + *a, qlen);
      HIT_LIST_count_useq_visited(FSS->HL, qlen);
    }
    if (dist > FSS->eps) continue;
    FSS->ifunc(FSS, *a, dist);
    HIT_LIST_count_useq_hit(FSS->HL, qlen);
    HIT_LIST_count_seq_hit(FSS->HL, 1);
  }
  HIT_LIST_count_seq_visited(FSS->HL, n);
}

static inline
void traverse_sarray(FSSRCH *FSS, uint32_t offset, uint32_t n)

{
  int qlen = FSS->qlen; 
  int **MD = FSS->M; 
  char *base = FSS->base;
  int *sa = FSS->sa + offset; 
  uint8_t *lcpx = FSS->lcpx + offset + 1;
  int *cdist = FSS->cdist;

  int i, j;
  int k = qlen+1;
  char *s;

  if (n==0) return;
  /* First sequence in any case */

  s = base + *sa;
  for (j=0; j < qlen; j++, s++) 
    cdist[j] = cdist[j-1] + MD[j][*s & A_SIZE_MASK];
  HIT_LIST_count_useq_visited(FSS->HL, qlen);

  if (cdist[j-1] <= FSS->eps) {
    FSS->ifunc(FSS, *sa, cdist[j-1]);
    HIT_LIST_count_useq_hit(FSS->HL, qlen);
    HIT_LIST_count_seq_hit(FSS->HL, 1);
  }
  sa++;

  for (i=1; i < n; i++, sa++, lcpx++) {
    j=min(*lcpx, qlen);
    if (j >= k) continue;
    else k = qlen+1;
    s = base + *sa + j;
    if (*(lcpx+1) < qlen) {
      for (; j < *(lcpx+1); j++, s++) {
	HIT_LIST_count_useq_visited(FSS->HL, 1);
	cdist[j] = cdist[j-1] + MD[j][*s & A_SIZE_MASK];
      }
      if (cdist[j-1] > FSS->eps) {
	k = j;
	continue;
      }
    }
    for (; j < qlen; j++, s++) { 
      HIT_LIST_count_useq_visited(FSS->HL, 1);
      cdist[j] = cdist[j-1] + MD[j][*s & A_SIZE_MASK];
    }
    if (cdist[j-1] <= FSS->eps) {
      FSS->ifunc(FSS, *sa, cdist[j-1]);
      HIT_LIST_count_useq_hit(FSS->HL, qlen);
      HIT_LIST_count_seq_hit(FSS->HL, 1);
    }
  }
  HIT_LIST_count_seq_visited(FSS->HL, n);
}

static
void process_bin_sarray(FSSRCH *FSS, FS_SEQ_t bin)
{
  uint32_t n = FSS->FSI->bins[bin+FSS->bincr] - FSS->FSI->bins[bin];
  if (n) traverse_sarray(FSS, FSS->FSI->bins[bin], n);
}

/********************************************************************/ 
/*      Main Search functions                                       */
/********************************************************************/ 

static 
void FS_search(FSSRCH *FSS, void *qseq0)
{
  int i, j, p;
  char *qseq = qseq0;
  FSINDX *FSI = FSS->FSI;

  FSS->qbin = FS_SEQ(FSI->ptable, qseq, FSS->qlen);

  FSS->sa = FSS->FSI->oa;
  FSS->lcpx = FSS->FSI->lcpb;

  if (FSS->qlen >= FSI->m) {
    FSS->len = FSI->m;
    FSS->bincr = 1;
  }
  else {
    FSS->len = FSS->qlen;
    FSS->bincr = FSI->ptable->hash[FSS->qlen-1][1];
  }

  for (i=0; i < FSS->len; i++) {
    int p0 = FSI->ptable->pttn[i][qseq[i] & A_SIZE_MASK];
    int h0 = FSI->ptable->hash[i][p0];
    for (p=0; p < FSS->nbpt[i]; p++) {
      FSS->nbpt_dist[i][p] = INT_MAX;
      FSS->nbpt_closest[i] = INT_MAX;
      if (p >= p0)
	FSS->nbpt_diff[i][p] = FSI->ptable->hash[i][p+1] - h0;
      else 
	FSS->nbpt_diff[i][p] = FSI->ptable->hash[i][p] - h0;
    }
    for (j=0; j < A_SIZE; j++) {
      p = FSI->ptable->pttn[i][j];
      if(p == -1 || p == p0) continue;
      else if (p > p0) p--;
      if (FSS->nbpt_dist[i][p] > FSS->M[i][j]) {
	FSS->nbpt_dist[i][p] = FSS->M[i][j];
	if (FSS->nbpt_closest[i] > FSS->nbpt_dist[i][p]) {
	  FSS->nbpt_closest[i] = FSS->nbpt_dist[i][p];
	}
      }
    }
  }
  /* Process the bin associated with the query sequence */
  HIT_LIST_count_FS_seq_visited(FSS->HL, 1);
  HIT_LIST_count_FS_seq_hit(FSS->HL, 1);
  FSS->pfunc(FSS, FSS->qbin); 
  /* Process neighbouring bins - recursive */
  check_bins(FSS, FSS->qbin, 0, 0); 
}


static
void sarray_search(FSSRCH *FSS, void *ptr)
{
  ULINT n = FSS->FSI->bins[FSS->FSI->no_bins] 
    - FSS->FSI->bins[1];
  FSS->sa = FSS->FSI->sa;
  FSS->lcpx = FSS->FSI->lcpx;
  if (FSS->FSI->use_sa)
    traverse_sarray(FSS, 1, n);
}

static
void seq_scan_search(FSSRCH *FSS, void *ptr)
{
  int j=0;
  int dist;

  int **M = FSS->M;
  int qlen = FSS->qlen;
  char *base = FSS->base;

  char *sseq;
  ULINT slen;
  MIXED_ITER *f_iter;

  /* Database iterator */
  f_iter = MIXED_ITER_init(FSS->FSI->s_db, qlen, qlen, 1);

  /* Sequential scan */
  while ((sseq = MIXED_ITER_next(f_iter, &slen))) {
    j++;
    dist = eval_dist(M, sseq, qlen);
    if (dist > FSS->eps) continue;
    FSS->ifunc(FSS, sseq-base, dist);
    HIT_LIST_count_seq_hit(FSS->HL, 1);
  }
  HIT_LIST_count_seq_visited(FSS->HL, j);
  MIXED_ITER_del(f_iter);
}

static
HIT_LIST_t *run_search(FSINDX *FSI, int thread, SCORE_MATRIX *M, 
		       BIOSEQ *query, int range, int kNN, 
		       int conv_type, HIT_LIST_t *HL, 
		       SFUNC_TYPE stype, PFUNC_TYPE ptype)
{
  int i, j;
  FSSRCH *FSS = (FSSRCH *)(FSI->FSS) + thread;
  HIT top;
  BIOSEQ subject;
  int sim;
  char *qseq;
  int qlen;
  BIOSEQ q1;
  void (*search_func) (FSSRCH *, void *);

  /* Check: index alphabet subset of matrix alphabet */
  if (!check_word_alphabet(FSI->ptable->alphabet, strlen(FSI->ptable->alphabet),
			   M->alphabet, M->alen))
    return NULL;

  /* Assign query */
  if (M->Mtype == SCORE) {
    qseq = query->start;
    qlen = query->len;
  }
  else {
    q1.start = qseq = M->qseq;
    q1.len = qlen = M->len;
    q1.id.defline = "Positional Matrix";
    query = &q1;
  }

  /* Check query alphabet */
  if (!check_word_alphabet(qseq, qlen, FSI->ptable->alphabet, 
			   strlen(FSI->ptable->alphabet)))
      return NULL;

  /* Init hit list */
  if (HL == NULL) {
    FSS->HL = HIT_LIST_create(query, FSI->s_db, M, range, 
			      kNN, conv_type);
  }
  else {
    HIT_LIST_reset(HL, query, FSI->s_db, M, range, kNN, conv_type);
    FSS->HL = HL;
  }


  HIT_LIST_index_data(FSS->HL, FSI->index_name, FSI->ptable->sepn, 
		      FSI->m, FSI->no_bins, FSI->no_seqs,
		      FSI->no_useqs);

  /* ************ Setup search **************************/
  /* Search-specific params */
  FSS->FSI = FSI;
  FSS->nbpt = FSI->ptable->nbpt;
  FSS->base = fastadb_data_pter(FSI->s_db, 0);
  FSS->qlen = qlen;

  if (FSS->qlen > FSS->Mlen) {
    FSS->M = reallocec(FSS->M, FSS->qlen * sizeof(int *));
    FSS->cdist--;
    FSS->cdist = reallocec(FSS->cdist, (FSS->qlen+1) * sizeof(int));
    FSS->cdist++;
    for (i=FSS->Mlen; i < FSS->qlen; i++) 
      FSS->M[i] = mallocec(A_SIZE * sizeof(int));
    FSS->Mlen = FSS->qlen;
  }

  FSS->kNN = kNN;

  switch (stype) 
    {
    case FS_BINS:
      search_func = FS_search;
      switch (ptype)
	{
	case SARRAY:
	  FSS->pfunc = process_bin_sarray;
	  break;
	case DUPS_ONLY:
	  FSS->pfunc = process_bin_dups;
	  break;
	case FULL_SCAN:
	default:
	  FSS->pfunc = process_bin_all;
	  break;
	}
      break;
    case SUFFIX_ARRAY:
      search_func = sarray_search;
      FSS->pfunc = NULL;
      break;
    case SEQ_SCAN:
    default:
      search_func = seq_scan_search;
      FSS->pfunc = NULL;
      break;
    }

  if (FSS->kNN == -1) {
    M->set_conv_type(M, conv_type);
    FSS->eps = M->range_conv(M, qseq, qlen, range);
    FSS->ifunc = FSSRCH_insert;
  }
  else {
    FSS->eps = INT_MAX;
    FSS->ifunc = FSSRCH_insert_queue;
  }

  /* Score matrix */
  /* Remark - if qseq[i] or j are not in the matrix alphabet
     M->item_conv(M, qseq[i], i, j) will return VERY_LARGE_DISTANCE,
     so that all such fragments will be rejected. Of course,
     this may fail for pathological queries. In this case, the only
     important value of j is '\0' - hence strings which go over the
     sequence boundary will be rejected.
  */
  for (i=0; i < qlen; i++)
    for (j=0; j < A_SIZE; j++)
      FSS->M[i][j] = M->item_conv(M, qseq[i], i, j);

  /* ************ Search proper *************************/
  search_func(FSS, qseq);

  /* Set off timer */
  HIT_LIST_stop_timer(FSS->HL);

  /********** Insert hits into hit list *****************/
  /* 1. hits in the FSS->hits array */
  for (i=0; i < FSS->hits_len; i++) {
    subject.start = FSS->base + FSS->hits[i].offset;
    subject.len = qlen;
    if (M->Stype == SIMILARITY) 
       sim = M->eval_score(M, qseq, subject.start, qlen);
    else
      sim = 0;
    HIT_LIST_insert_hit(FSS->HL, &subject, FSS->hits[i].dist, sim);
  }
  FSS->hits_len = 0;
  /* 2. hits in priority queue */
  if (FSS->kNN != -1) {
    while (pqremove(FSS->p_queue, &top) != NULL) {
      subject.start = FSS->base + top.offset;
      subject.len = qlen;
      if (M->Stype == SIMILARITY) 
	sim = M->eval_score(M, qseq, subject.start, qlen);
      else
	sim = 0;
      HIT_LIST_insert_hit(FSS->HL, &subject, top.dist, sim);
    }
  }
  FSS->no_hits = 0;
  return FSS->HL;
}

HIT_LIST_t *FSINDX_rng_srch(FSINDX *FSI, BIOSEQ *query, SCORE_MATRIX *M,
			    int range, int conv_type, HIT_LIST_t *HL,  
			    SFUNC_TYPE stype, PFUNC_TYPE ptype)
{
  HL = run_search(FSI, 0, M, query, range, -1, conv_type, HL, stype, ptype);
  if (HL == NULL) 
    Throw FSexcept(RUNTIME_ERR,
		   "FSINDX_rng_srch(): Query not consistent with index.");
  return HL;
}

HIT_LIST_t *FSINDX_kNN_srch(FSINDX *FSI, BIOSEQ *query, SCORE_MATRIX *M,
			    int kNN, HIT_LIST_t *HL, SFUNC_TYPE stype, 
			    PFUNC_TYPE ptype)
{
  HL = run_search(FSI, 0, M, query, 0, kNN, 0, HL, stype, ptype);
  if (HL == NULL) 
    Throw FSexcept(RUNTIME_ERR,
		   "FS_INDX_kNN_srch(): Query not consistent with index.");
  return HL;
}

/********************************************************************/ 
/*      Threaded Search                                             */
/********************************************************************/ 
int next_avail_srch;

#if THREADS > 1
#include <pthread.h>
static pthread_t callThd[THREADS];
static pthread_mutex_t mutexnext;

struct THREAD_ARGS_s
{
  FSINDX *FSI;
  SRCH_ARGS *args;
  int n;  
  HIT_LIST_t *HL;
  int id;
};


static
void *run_thread(void *targs0)
{
  struct THREAD_ARGS_s *targs = targs0;
  int cs; /* Current search */
  BIOSEQ *query;
  SCORE_MATRIX *M;
  int range;
  int kNN;
  int conv_type;
  HIT_LIST_t *HL;
  SFUNC_TYPE stype; 
  PFUNC_TYPE ptype;  

  while (1) {
    /* Lock mutex before accessing next search */
    pthread_mutex_lock(&mutexnext);
    cs = next_avail_srch++;
    /*printf("Thread #%d, item %d out of %d\n", targs->id, cs, targs->n); */
    pthread_mutex_unlock(&mutexnext);

    if (cs >= targs->n)
      break;

    /* This was done for easier debugging */
    query = &(targs->args[cs].query);
    M = targs->args[cs].M;
    range = targs->args[cs].range;
    kNN = targs->args[cs].kNN;
    conv_type = M->conv_type;
    HL = targs->HL + cs;
    stype = targs->args[cs].stype; 
    ptype = targs->args[cs].ptype;  

    run_search(targs->FSI, targs->id, M, query, range, kNN, 
	       conv_type, HL, stype, ptype);
  }
  return (void*) 0;
}

HIT_LIST_t *FSINDX_threaded_srch(FSINDX *FSI, SRCH_ARGS *args, int n)
{
  int i;
  struct THREAD_ARGS_s targs[THREADS];
  HIT_LIST_t *HL = callocec(n, sizeof(HIT_LIST_t));
  pthread_attr_t attr;
  int status;

  pthread_mutex_init(&mutexnext, NULL);  
  next_avail_srch = 0;

  /* Create threads */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  for (i=0; i < THREADS; i++) {
    targs[i].FSI = FSI;
    targs[i].args = args;
    targs[i].HL = HL;
    targs[i].n = n;
    targs[i].id = i;

    pthread_create(&callThd[i], &attr, run_thread, (void *)(&targs[i]));
  }

  /* Wait on remaining threads */
  for (i=0; i < THREADS; i++) 
    pthread_join(callThd[i], (void **)&status);

  /* Destroy mutex and attribute */
  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&mutexnext);
  return HL;
}
#else
HIT_LIST_t *FSINDX_threaded_srch(FSINDX *FSI, SRCH_ARGS *args, int n)
{
  return NULL;
}
#endif /* #if THREADS > 1 */

/********************************************************************/ 
/*      Access functions with error checking                        */
/********************************************************************/ 

uint32_t FS_INDEX_get_bin_size(FSINDX *FSI, FS_SEQ_t bin)
{
  if (bin >= FSI->no_bins)
    Throw FSexcept(INDEX_OUT_OF_RANGE,
		   "FS_INDEX_get_bin_size(): Index out of range.");
  else
    return FSI->bins[bin+1] - FSI->bins[bin];
}

uint32_t FS_INDEX_get_unique_bin_size(FSINDX *FSI, FS_SEQ_t bin)
{
  if (bin >= FSI->no_bins)
    Throw FSexcept(INDEX_OUT_OF_RANGE,
		   "FS_INDEX_get_bin_size(): Index out of range.");
  else
    return get_unique_size(FSI, FSI->m, bin);
}
#if 0 
BIOSEQ *FS_INDEX_get_seq(FSINDX *FSI, uint32_t bin, uint32_t seq)
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




static 
void pull_duplicates(FSINDX *FSI, HIT_LIST_t *HL)
{
  ULINT n;
  ULINT N;
  ULINT frag_offset;
  BIOSEQ subject;
  ULINT i;
  ULINT j;
  ULINT k;
  ULINT M = HL->actual_seqs_hits;
  SEQ_HIT_t *h;
  ULINT bin0;
  float val;

  for (k=0; k < M; k++)
    {
      h = HIT_LIST_get_hit(HL, k);
      bin0 = h->bin;
      i = h->pos;
      val = h->value;
      n = FSI->u_size[bin0];
      N = FSI->bin_size[bin0];

      j = FSI->u[bin0][i] + 1;
      while ((j < N && i == n-1) || j < FSI->u[bin0][i+1])
	{
	  HIT_LIST_count_seq_hit(HL, 1);
	  frag_offset = FSI->bin[bin0][j];
	  fastadb_get_Ffrag(FSI->s_db, FSI->m, &subject, 
			    frag_offset);
	  
	  HIT_LIST_insert_seq_hit(HL, &subject, val); 
	  j++;
	}
    }
}


#endif
