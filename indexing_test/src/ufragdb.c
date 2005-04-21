#include "ufragdb.h"
#include "avl.h"

int cmp_ufrags(const void *S1, const void *S2)
{
  ULINT *T1 = (ULINT *)S1;
  ULINT *T2 = (ULINT *)S2;

  if (*T1 > *T2)
    return 1;
  else if (*T1 < *T2)
    return -1;
  else
    return 0;
}


UFRAG_DB *ufragdb_init(void)
{
  UFRAG_DB *udb = callocec(1, sizeof(UFRAG_DB));
  return udb;
}

void ufragdb_del(UFRAG_DB *udb)
{
  if (udb == NULL) return;
  free(udb->db_name);
  if (udb->sdb != NULL)
    fastadb_close(udb->sdb);
  if (udb->fptr != NULL)
    fclose(udb->fptr);
  free(udb->pts);
  free(udb->dup_offset);
  free(udb->dup_heap);
  free(udb);
  udb = NULL;
}

typedef struct
{
  ULINT id;
  ULINT pt;
} DUP_FRAG;

static
int compare_seqstr(const void *pa, const void *pb, void *param)
{
  ULINT len = *((ULINT *)param);
  const char *a = pa;
  const char *b = pb;

  return memcmp(a, b, len);
}

static
int cmp_dupfrag(const void *pa, const void *pb)
{
  const DUP_FRAG *a = pa;
  const DUP_FRAG *b = pb;

  if (a->pt > b->pt)
    return 1;
  else if (a->pt < b->pt)
    return -1;
  else
    return 0;
}

#define DFRG_INCR 65536

void ufragdb_create(UFRAG_DB *udb, const char *db_name, 
		    ULINT frag_len, int skip)
{
  EXCEPTION *except;

  char *seqstr;
  char **tree_item;
  struct avl_table *AVLtree = NULL;

  ULINT i, j, k, l;
  ULINT no_frags;
  ULINT one_percent_fragments;
  int valid = 1;

  DUP_FRAG *dfrg = NULL;
  ULINT max_dfrg = DFRG_INCR;

  char *db_dir = NULL;
  char *db_base = NULL;

  char *_base_;
  char *_end_;

  /* We want to start in clean state in case something fails */
  if (udb->frag_len > 0)
    return;
  
  Try {
    path_split(db_name, &db_dir, &db_base);

    udb->db_name_len = strlen(db_base);
    udb->db_name = strdup(db_base);
    udb->frag_len = frag_len;
    udb->skip = skip;
 
    udb->sdb = fastadb_open(db_name); 

    no_frags = fastadb_count_Ffrags(udb->sdb, frag_len);
    one_percent_fragments = (ULINT) (((double) no_frags/skip) / 100);

    /* Allocate to maximum possible number so we don't have to
       grow the array later */
    udb->pts = mallocec(no_frags * sizeof(ULINT));

    dfrg = mallocec(max_dfrg * sizeof(DUP_FRAG));

    AVLtree = avl_create(compare_seqstr, &frag_len, NULL);

    _base_ = fastadb_data_pter(udb->sdb, 0);  
    _end_ = fastadb_end_heap(udb->sdb);

    j=0;
    fprintf(stdout, "Sorting fragments.\n");
    for (seqstr=_base_; seqstr <= _end_; seqstr++) {
      j++;

      /* Check if the fragment is valid */
      valid = 1;
      for (l=0; l < frag_len; l++, seqstr++) {
	if (*seqstr == '\0' || *seqstr == 'B' || 
	    *seqstr == 'Z' || *seqstr == 'X' )
	  valid = 0;
	  break;
      } 

      if (valid) {
	i = seqstr - _base_;
	tree_item = (char **)avl_probe(AVLtree, seqstr);	  

	if ((*tree_item) == seqstr) {
	  udb->pts[udb->pts_size++] = i;
	}
	else { 
	  dfrg[udb->dup_heap_size].id = i;
	  dfrg[udb->dup_heap_size++].pt = 
	    (*tree_item) - udb->sdb->seq_data;
	  if (udb->dup_heap_size >= max_dfrg) {
	    max_dfrg += DFRG_INCR;
	    dfrg = reallocec(dfrg, max_dfrg * sizeof(DUP_FRAG));
	  }
	}
      }
      // Progress bar
      printbar(stdout, j+1, one_percent_fragments, 50);  
    }  

    avl_destroy (AVLtree, NULL);
    AVLtree = NULL;

    /* This algorithm works because the udb->pts array is already 
       sorted. We sort dfrg first by pt (i.e. unique point the duplicate
       is equal to). Then, we traverse it, swapping elements of udb->pts
       as we go. Thus we get what we want, udb->pts whose first
       udb->dpts_size elements represent fragments with duplicates, all
       sorted by id.
    */

    qsort(dfrg, udb->dup_heap_size, sizeof(DUP_FRAG), cmp_dupfrag);

    i=-1;
    j=0;
    
    udb->dup_heap = mallocec(udb->dup_heap_size * sizeof(ULINT));

    fprintf(stdout, "Copying duplicates.\n");
    for (k=0; k < udb->dup_heap_size; k++) {
      udb->dup_heap[k] = dfrg[k].id;
      if (k==0 || dfrg[k].pt != dfrg[k-1].pt) {
	i++;
	while (udb->pts[j] != dfrg[k].pt)
	  j++;
	udb->pts[j++] = udb->pts[i];	
	udb->pts[i] = dfrg[k].pt;
      }
    }

    udb->dpts_size = i+1;
    udb->upts_size = udb->pts_size - udb->dpts_size;

    udb->dup_offset = mallocec(udb->dpts_size * sizeof(ULINT));
    udb->dup_offset[0]=0;
    for (i=0, k=1; k < udb->dup_heap_size; k++)
      if (dfrg[k].pt != dfrg[k-1].pt) {
	udb->dup_offset[++i] = k;
      }

    free(dfrg);
    free(db_dir);
    free(db_base);
  }
  Catch(except) {
    free(dfrg);
    free(db_dir);
    free(db_base);
    
    free(udb->db_name);
    if (udb->sdb != NULL)
      fastadb_close(udb->sdb);
    if (udb->fptr != NULL)
      fclose(udb->fptr);
    free(udb->pts);
    free(udb->dup_offset);
    free(udb->dup_heap);
    memset(udb, 0, sizeof(UFRAG_DB));

    Throw except;
  }

  return;
}

void ufragdb_read(UFRAG_DB *udb, FILE *fp, char *db_dir)
{
  char *db_full;

  fprintf(stdout, "Loading fragments.\n");
  if (udb->fptr != NULL)
    fclose(udb->fptr);
  udb->fptr = fp;

  fread(&udb->db_name_len, sizeof(int), 1, udb->fptr);
  fread(&udb->frag_len, sizeof(ULINT), 1, udb->fptr);
  fread(&udb->skip, sizeof(ULINT), 1, udb->fptr);
  fread(&udb->pts_size, sizeof(ULINT), 1, udb->fptr);
  fread(&udb->upts_size, sizeof(ULINT), 1, udb->fptr);
  fread(&udb->dpts_size, sizeof(ULINT), 1, udb->fptr);
  fread(&udb->dup_heap_size, sizeof(ULINT), 1, udb->fptr);
  
  udb->db_name = mallocec(udb->db_name_len+1);
  udb->db_name[udb->db_name_len] = '\0';

  fread(udb->db_name, 1, udb->db_name_len, udb->fptr);

  if (db_dir != NULL) {
    db_full = mallocec(strlen(db_dir)+strlen(udb->db_name)+2);
    strcpy(db_full, db_dir);
    strcat(db_full, "/");
    strcat(db_full, udb->db_name);
  } 
  else
    db_full = udb->db_name; 

  udb->sdb = fastadb_open(db_full); 

  udb->pts = mallocec(udb->pts_size * sizeof(ULINT));
  udb->dup_offset = mallocec(udb->dpts_size * sizeof(ULINT));
  udb->dup_heap = mallocec(udb->dup_heap_size * sizeof(ULINT));

  fread(udb->pts, sizeof(ULINT), udb->pts_size, udb->fptr);
  fread(udb->dup_offset, sizeof(ULINT), udb->dpts_size, udb->fptr);
  fread(udb->dup_heap, sizeof(ULINT), udb->dup_heap_size, udb->fptr);

  udb->fptr = NULL;
  
  if (db_dir != NULL)
    free(db_full);
}


void ufragdb_write(UFRAG_DB *udb, FILE *fp)
{
  if (udb->fptr != NULL)
    fclose(udb->fptr);
  udb->fptr = fp;

  fwrite(&udb->db_name_len, sizeof(int), 1, udb->fptr);
  fwrite(&udb->frag_len, sizeof(ULINT), 1, udb->fptr);
  fwrite(&udb->skip, sizeof(ULINT), 1, udb->fptr);
  fwrite(&udb->pts_size, sizeof(ULINT), 1, udb->fptr);
  fwrite(&udb->upts_size, sizeof(ULINT), 1, udb->fptr);
  fwrite(&udb->dpts_size, sizeof(ULINT), 1, udb->fptr);
  fwrite(&udb->dup_heap_size, sizeof(ULINT), 1, udb->fptr);
  
  fwrite(udb->db_name, 1, udb->db_name_len, udb->fptr);

  fwrite(udb->pts, sizeof(ULINT), udb->pts_size, udb->fptr);
  fwrite(udb->dup_offset, sizeof(ULINT), udb->dpts_size, udb->fptr);
  fwrite(udb->dup_heap, sizeof(ULINT), udb->dup_heap_size, udb->fptr);

  udb->fptr = NULL;
}



void ufragdb_load(UFRAG_DB *udb, const char *udb_name,
		  int options)
{
  FILE *fp;
  char *db_dir = NULL;
  char *db_base = NULL;

  fp = fopen(udb_name, "rb");
  if(fp == NULL)
    Throw FSexcept(FOPEN_ERR, 
		   "ufragdb_load(): Could not open the file %s.",
		   udb_name);

  path_split(udb_name, &db_dir, &db_base);

  ufragdb_read(udb, fp, db_dir);

  fclose(fp);
  free(db_dir);
  free(db_base);
}

void ufragdb_save(UFRAG_DB *udb, const char *udb_name)
{
  FILE *fp;

  fp = fopen(udb_name, "wb");
  if(fp == NULL)
    Throw FSexcept(FOPEN_ERR, 
		   "ufragdb_save(): Could not open the file %s.",
		   udb_name);
  ufragdb_write(udb, fp);
  fclose(fp);
}

ULINT ufragdb_get_nopts(UFRAG_DB *udb)
{
  return udb->pts_size;
}

int ufragdb_get_ufrag(UFRAG_DB *udb, ULINT i)
{
  if (i >= udb->pts_size) 
    return -1;
  return udb->pts[i];
}

void ufragdb_init_ufrags(UFRAG_DB *udb, ULINT i0)
{
  if (i0 >= udb->pts_size)
    udb->cur_pt = udb->pts_size;
  else
    udb->cur_pt = i0;
}

int ufragdb_get_next_ufrag(UFRAG_DB *udb)
{
  if (udb->cur_pt >= udb->pts_size) 
    return -1;
  return udb->pts[udb->cur_pt++];
}

int ufragdb_get_dfrags(UFRAG_DB *udb, ULINT frag_offset, 
		       ULINT **dfrags, ULINT *dsize)
{
  ULINT *item;
  ULINT i;

  /* find if there is a duplicate with such offset */
  item = bsearch(&frag_offset, udb->pts, udb->dpts_size,
		 sizeof(ULINT), cmp_ufrags);

  if (item == NULL) {
    *dsize = 0;
    *dfrags = NULL;
    return 0;
  }
  i = item - udb->pts;  
  *dsize = i==0 ? 0 : udb->dup_offset[i] - udb->dup_offset[i-1];
  *dfrags = udb->dup_heap + udb->dup_offset[i];
  return 1;
}

int ufragdb_get_dfrags2(UFRAG_DB *udb, ULINT i, 
			ULINT **dfrags, ULINT *dsize)
{
  if (i >= udb->dpts_size) {
    *dsize = 0;
    *dfrags = NULL;
    return 0;
  }
  else if (i == udb->dpts_size-1) 
    *dsize = udb->dpts_size - udb->dup_offset[i];
  else 
    *dsize = udb->dup_offset[i+1] - udb->dup_offset[i];

  *dfrags = udb->dup_heap + udb->dup_offset[i];
  return 1;
}

void ufragdb_print(UFRAG_DB *udb, FILE *stream)
{
  fprintf(stream, "\n*** UFRAG_DB Statistics ***\n");
  fprintf(stream, "Database: %s\n", udb->db_name);
  fprintf(stream, "Fragment length: %ld\n", udb->frag_len);
  fprintf(stream, "Points: %ld\n", udb->pts_size);
  fprintf(stream, "Unique points: %ld\n", udb->upts_size);
  fprintf(stream, "Points with duplicates: %ld\n", udb->dpts_size);
  fprintf(stream, "Total duplicates: %ld\n", udb->dup_heap_size);
}
