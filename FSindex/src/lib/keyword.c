#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "misclib.h"
#include "keyword.h"

#define HEAP_CHUNK 4096
KW_INDEX *KW_INDEX_load_txt(const char *basename)
{
  /*extern int getline(char **_lineptr, size_t *_n, FILE *_stream);*/
  /* four files:
     1. basename.scl: maps sequence to its cluster
       
     2. basename.cfd: maps cluster to its frequency and string.

     3. basename.skw: maps sequence to its keywords
       
     4. basename.kfd: maps keyword to its frequency and string.

  */
  KW_INDEX *KWI = callocec(1, sizeof(KW_INDEX));

  char *linebuf = NULL;
  int bufsize = 0;
  char *bptr;
  int len;

  int i;
  int s;
  int t;

  ULINT M;
  

  /* Text files to open */
  char *scl_name = mallocec(strlen(basename)+5);
  char *cfd_name = mallocec(strlen(basename)+5);
  char *skw_name = mallocec(strlen(basename)+5);
  char *kfd_name = mallocec(strlen(basename)+5);
  FILE *scl;
  FILE *cfd;
  FILE *skw;
  FILE *kfd;


  /* Concatenate names */
  strcpy(scl_name, basename);  
  strcat(scl_name, ".scl");
  strcpy(cfd_name, basename);  
  strcat(cfd_name, ".cfd");
  strcpy(skw_name, basename);  
  strcat(skw_name, ".skw");
  strcpy(kfd_name, basename);  
  strcat(kfd_name, ".kfd");


  /* Open all four files, fail if missing */
  if(((scl = fopen(scl_name, "r")) == NULL) ||
     ((cfd = fopen(cfd_name, "r")) == NULL) ||
     ((skw = fopen(skw_name, "r")) == NULL) ||
     ((kfd = fopen(kfd_name, "r")) == NULL))
    {
      fprintf(stderr, "Could not open one of keyword files!\n");
      exit(EXIT_FAILURE);
    }


  /* IMPORTANT: we count sequences from 0 in our C library routines while in
     keyword files the count goes form 1. 
  */

  /* basename.scl */

  getline(&linebuf, &bufsize, scl);
  sscanf(linebuf, "%d", &KWI->ns);
  KWI->seq2cl = callocec(KWI->ns, sizeof(int));
  KWI->seq2kw0 = callocec(KWI->ns, sizeof(int));
  KWI->seq2kw1 = callocec(KWI->ns, sizeof(int));

  KWI->nsc = 0;
  while (getline(&linebuf, &bufsize, scl) != EOF)
    {
      sscanf(linebuf, "%d\t%d", &s, &t);
      s--;
      KWI->seq2cl[s] = t;
      KWI->nsc++;
    }
  fclose(scl);

  /* basename.skw */

  getline(&linebuf, &bufsize, skw);
  KWI->nsk = 0;
  KWI->skw_hs = 0;
  M = HEAP_CHUNK;
  KWI->skw_h = mallocec(M * sizeof(int));
  while (getline(&linebuf, &bufsize, skw) != EOF)
    {
      sscanf(linebuf, "%d\t%d", &s, &t);
      s--;
      KWI->seq2kw1[s] = t;
      KWI->seq2kw0[s] = KWI->skw_hs;
      KWI->skw_hs += t;
      if (KWI->skw_hs >= M)
	{
	  M += HEAP_CHUNK;
	  KWI->skw_h = reallocec(KWI->skw_h, M * sizeof(int));
	}
      bptr = linebuf;
      bptr = strstr(bptr, "\t") + 1;
      for (i=0; i < t; i++)
	{
	  bptr = strstr(bptr, "\t") + 1;
	  sscanf(bptr, "%d", KWI->skw_h+KWI->seq2kw0[s]+i);
	}
      KWI->nsk++;
    }
  fclose(skw);

  /* basename.cfd */

  getline(&linebuf, &bufsize, cfd);
  sscanf(linebuf, "%d", &KWI->nc);
  KWI->cld_hs = 0;
  M = HEAP_CHUNK;
  KWI->cld_h = mallocec(M);
  KWI->clf = mallocec((KWI->nc+1) * sizeof(int));
  KWI->cld = mallocec((KWI->nc+1) * sizeof(ULINT));
  
  s = 1;
  while (getline(&linebuf, &bufsize, cfd) != EOF)
    {
      sscanf(linebuf, "%d\t", &t);
      KWI->clf[s] = t;
      bptr = strstr(linebuf, "\t") + 1;
      KWI->cld[s] = KWI->cld_hs;
     
      /* Remember there is '\n' at the end that we don't want */
      len = strlen(bptr);
      KWI->cld_hs += len;
      if (KWI->cld_hs >= M)
	{
	  M += HEAP_CHUNK;
	  KWI->cld_h = reallocec(KWI->cld_h, M);
	}
      memcpy(KWI->cld_h+KWI->cld[s], bptr, len-1);
      KWI->cld_h[KWI->cld[s] + len - 1] = '\0';
      s++;
    }
  fclose(cfd);

  /* basename.kfd */

  getline(&linebuf, &bufsize, kfd);
  sscanf(linebuf, "%d", &KWI->nk);
  KWI->kwd_hs = 0;
  M = HEAP_CHUNK;
  KWI->kwd_h = mallocec(M);
  KWI->kwf = mallocec((KWI->nk+1) * sizeof(int));
  KWI->kwd = mallocec((KWI->nk+1) * sizeof(ULINT));
  
  s = 1;
  while (getline(&linebuf, &bufsize, kfd) != EOF)
    {
      sscanf(linebuf, "%d\t", &t);
      KWI->kwf[s] = t;
      bptr = strstr(linebuf, "\t") + 1;
      KWI->kwd[s] = KWI->kwd_hs;
     
      /* Remember there is '\n' at the end that we don't want */
      len = strlen(bptr);
      KWI->kwd_hs += len;
      if (KWI->kwd_hs >= M)
	{
	  M += HEAP_CHUNK;
	  KWI->kwd_h = reallocec(KWI->kwd_h, M);
	}
      memcpy(KWI->kwd_h+KWI->kwd[s], bptr, len-1);
      KWI->kwd_h[KWI->kwd[s] + len - 1] = '\0';
      s++;
    }
  fclose(kfd);
  return KWI;
}


void KW_INDEX_save_bin(KW_INDEX *KWI, const char *filename)
{
  FILE *stream = fopen(filename, "wb");

  if(stream == NULL)
    {
      fprintf(stderr, 
	      "KW_INDEX_save_bin(): Could not open the file %s\n",
	      filename);
      exit(EXIT_FAILURE);
    }
  
  fwrite(&KWI->ns, sizeof(int), 1, stream);
  fwrite(&KWI->nsc, sizeof(int), 1, stream);
  fwrite(&KWI->nsk, sizeof(int), 1, stream);
  fwrite(&KWI->nc, sizeof(int), 1, stream);
  fwrite(&KWI->nk, sizeof(int), 1, stream);
  fwrite(&KWI->cld_hs, sizeof(int), 1, stream);
  fwrite(&KWI->kwd_hs, sizeof(int), 1, stream);
  fwrite(&KWI->skw_hs, sizeof(int), 1, stream);


  fwrite(KWI->clf, sizeof(int), KWI->nc, stream);
  fwrite(KWI->kwf, sizeof(int), KWI->nk, stream);
  fwrite(KWI->cld, sizeof(ULINT), KWI->nc, stream);
  fwrite(KWI->kwd, sizeof(ULINT), KWI->nk, stream);

  fwrite(KWI->seq2cl, sizeof(int), KWI->ns, stream);
  fwrite(KWI->seq2kw0, sizeof(int), KWI->ns, stream);
  fwrite(KWI->seq2kw1, sizeof(int), KWI->ns, stream);

  fwrite(KWI->cld_h, sizeof(char), KWI->cld_hs, stream);
  fwrite(KWI->kwd_h, sizeof(char), KWI->kwd_hs, stream);
  fwrite(KWI->skw_h, sizeof(int), KWI->skw_hs, stream);

  fclose(stream);
}

KW_INDEX *KW_INDEX_load_bin(const char *filename)
{
  KW_INDEX *KWI = mallocec(sizeof(KW_INDEX));
  FILE *stream = fopen(filename, "r");

  if(stream == NULL)
    {
      fprintf(stderr, "KW_INDEX_load_bin():"
	      " Could not open the file %s\n", filename);
      exit(EXIT_FAILURE);
    }
  
  fread(&KWI->ns, sizeof(int), 1, stream);
  fread(&KWI->nsc, sizeof(int), 1, stream);
  fread(&KWI->nsk, sizeof(int), 1, stream);
  fread(&KWI->nc, sizeof(int), 1, stream);
  fread(&KWI->nk, sizeof(int), 1, stream);
  fread(&KWI->cld_hs, sizeof(int), 1, stream);
  fread(&KWI->kwd_hs, sizeof(int), 1, stream);
  fread(&KWI->skw_hs, sizeof(int), 1, stream);

  KWI->clf = mallocec(KWI->nc * sizeof(int));
  KWI->kwf = mallocec(KWI->nk * sizeof(int));
  KWI->cld = mallocec(KWI->nc * sizeof(ULINT));
  KWI->kwd = mallocec(KWI->nk * sizeof(ULINT));
  KWI->seq2cl = mallocec(KWI->ns * sizeof(int));
  KWI->seq2kw0 = mallocec(KWI->ns * sizeof(int));
  KWI->seq2kw1 = mallocec(KWI->ns * sizeof(int));
  KWI->cld_h = mallocec(KWI->cld_hs * sizeof(char));
  KWI->kwd_h = mallocec(KWI->kwd_hs * sizeof(char));
  KWI->skw_h = mallocec(KWI->skw_hs * sizeof(int));

  fread(KWI->clf, sizeof(int), KWI->nc, stream);
  fread(KWI->kwf, sizeof(int), KWI->nk, stream);
  fread(KWI->cld, sizeof(ULINT), KWI->nc, stream);
  fread(KWI->kwd, sizeof(ULINT), KWI->nk, stream);

  fread(KWI->seq2cl, sizeof(int), KWI->ns, stream);
  fread(KWI->seq2kw0, sizeof(int), KWI->ns, stream);
  fread(KWI->seq2kw1, sizeof(int), KWI->ns, stream);

  fread(KWI->cld_h, sizeof(char), KWI->cld_hs, stream);
  fread(KWI->kwd_h, sizeof(char), KWI->kwd_hs, stream);
  fread(KWI->skw_h, sizeof(int), KWI->skw_hs, stream);

  fclose(stream);
  return KWI;
}


void KW_INDEX_clear(KW_INDEX *KWI)
{
  free(KWI->clf);
  free(KWI->kwf);
  free(KWI->cld);
  free(KWI->kwd);
  free(KWI->seq2cl);
  free(KWI->seq2kw0);
  free(KWI->seq2kw1);
  free(KWI->cld_h);
  free(KWI->kwd_h);
  free(KWI->skw_h);
  
  free(KWI);
}


#if 0
/********************************************************************/ 
/*      Filtering                                                   */
/********************************************************************/ 

static
void FSSRCH_filter(FSSRCH *FSS)
{
  int offset1;
  int offset2;
  int no_hits =  FSS->args->HL->actual_seqs_hits;
  FPARAMS *fp = FSS->fp;
  int i;

  if (fp == NULL) return;
  

  /* Z-scores */
  srand48(time(NULL)); 
  hit_Zscores(FSS);

  HIT_LIST_sort_evalue(FSS->args->HL, 0);

  i=0;
  while (FSS->args->HL->hits[i].evalue < fp->ZE0 && i < no_hits)
    i++;
  offset1 = i;
  while (FSS->args->HL->hits[i].evalue <= fp->ZE1 && i < no_hits)
    i++;
  FSS->args->HL->actual_seqs_hits = i;
  no_hits =  FSS->args->HL->actual_seqs_hits;
  

  /* TO DO:
     Deallocate hits so deleted - this is a memory leak */

  /* Conservation ratios */
  HIT_LIST_process_cl(FSS->args->HL, offset1);
  /* Must re-sort hit list */
  HIT_LIST_sort_evalue(FSS->args->HL, 0);
  HIT_LIST_sort_cratio(FSS->args->HL, offset1);
  i=offset1;
  while (FSS->args->HL->hits[i].cratio > fp->CR0 && i < no_hits)
    i++;
  offset2 = i;
  while (FSS->args->HL->hits[i].cratio >= fp->CR1 && i < no_hits)
    i++;
  FSS->args->HL->actual_seqs_hits = i;
  no_hits =  FSS->args->HL->actual_seqs_hits;


  /* Keywords */
  HIT_LIST_process_kw(FSS->args->HL, offset2);
  /* Must re-sort hit list */
  HIT_LIST_sort_evalue(FSS->args->HL, 0);
  HIT_LIST_sort_cratio(FSS->args->HL, offset1);
  HIT_LIST_sort_kwscore(FSS->args->HL, offset2);
  i=offset2;

  while (FSS->args->HL->hits[i].kw_score >= fp->KW0 && i < no_hits)
    i++;
  FSS->args->HL->actual_seqs_hits = i;
}


/********************************************************************/    
/*                                                                  */
/*                   WRD_CNTR module                                */ 
/*                                                                  */
/********************************************************************/    

#define WC_STEP 500

void WRD_CNTR_add(WRD_CNTR **wc, int wc_size, int *wc_max)
{
  while (wc_size >= *wc_max)
    {
      *wc_max += WC_STEP;
      *wc = reallocec(*wc, *wc_max * sizeof(WRD_CNTR));
    } 
}

void WRD_CNTR_clear(WRD_CNTR *wc)
{
  free(wc);
}

static
int WRD_CNTR_cmp_score(const void *S1, const void *S2)
{
  const WRD_CNTR *T1 = S1;
  const WRD_CNTR *T2 = S2;

  /* Sorting in decreasing order */
  if (T2->score > T1->score)
    return 1;
  else if (T1->score > T2->score)
    return -1;
  else
    return 0;
}

void WRD_CNTR_sort_score(WRD_CNTR *wc, int wc_size)
{
  qsort(wc, wc_size, sizeof(WRD_CNTR), WRD_CNTR_cmp_score);
}



/* Keywords */

/* Use AVL tree to store clusters / keywords */

static 
int compare_items(const void *pa, const void *pb, void *param)
{
  const WRD_CNTR *a = pa;
  const WRD_CNTR *b = pb;
  
  return a->item - b->item;
}

void HIT_LIST_process_cl(HIT_LIST_t *HL, int offset)
{
  int i;
  KW_INDEX *KWI = HL->KWI;
  ULINT seqid;
  ULINT seqid0;
  struct avl_table *tree;
  WRD_CNTR **tree_item;
  WRD_CNTR *wc_tmp;
  WRD_CNTR wc_tmp1;

  
  if (KWI == NULL)
    return;
  HIT_LIST_sort_by_sequence(HL);
  seqid0 = KWI->ns;
  /*******************************************/
  /*          ORTHOLOGOUS CLUSTERS           */
  /*******************************************/

  tree = avl_create(compare_items, NULL, NULL);
  WRD_CNTR_add(&HL->oc, HL->oc_size, &HL->oc_max);

  /* Load tree */
  for (i=0; i < HL->actual_seqs_hits; i++)
    {
      seqid = HL->hits[i].sequence_id;
      if (seqid != seqid0)
	{
	  seqid0 = seqid;
	  if (KWI->seq2cl[seqid] != 0)
	    {
	      wc_tmp = HL->oc + HL->oc_size;
	      wc_tmp->item = KWI->seq2cl[seqid];
	      wc_tmp->count = 0;
	      wc_tmp->score = 0.0;
	      tree_item = (WRD_CNTR **) avl_probe (tree, wc_tmp);	  
	      (*tree_item)->count++;
	      if ((*tree_item) == wc_tmp)
		{
		  HL->oc_size++;
		  WRD_CNTR_add(&HL->oc, HL->oc_size, &HL->oc_max);
		}
	
	    }
	}
    }

  /* Calculate conservation ratio */
  for (i=0; i < HL->oc_size; i++)
    HL->oc[i].score = 
      (double) HL->oc[i].count / KWI->clf[HL->oc[i].item];

  /* Assign ratios to hits */
  for (i=0; i < HL->actual_seqs_hits; i++)
    {
      seqid = HL->hits[i].sequence_id;
      if (KWI->seq2cl[seqid] != 0)
	{
	  HL->hits[i].oc_cluster = KWI->seq2cl[seqid];
	  wc_tmp1.item = KWI->seq2cl[seqid];
	  tree_item = (WRD_CNTR **) avl_probe (tree, &wc_tmp1);	  
	  HL->hits[i].cratio = (*tree_item)->score;
	}
      else
	{
	  HL->hits[i].oc_cluster = 0;
	  HL->hits[i].cratio = 0.0;
	}
    }

  /* Clean up tree */
  avl_destroy (tree, NULL);
  WRD_CNTR_sort_score(HL->oc, HL->oc_size);
}

void HIT_LIST_process_kw(HIT_LIST_t *HL, int offset)
{
  int i;
  int j;
  KW_INDEX *KWI = HL->KWI;
  ULINT seqid;
  ULINT seqid0;
  struct avl_table *tree;
  WRD_CNTR **tree_item;
  WRD_CNTR *wc_tmp;
  WRD_CNTR wc_tmp1;
  int kw_hits = 0;
  double s;

  if (KWI == NULL)
    return;
  HIT_LIST_sort_by_sequence(HL);
  seqid0 = KWI->ns;

  /*******************************************/
  /*               KEYWORDS                  */
  /*******************************************/

  tree = avl_create(compare_items, NULL, NULL);
  WRD_CNTR_add(&HL->kw, HL->kw_size, &HL->kw_max);

  /* Load tree */
  for (i=0; i < HL->actual_seqs_hits; i++)
    {
      seqid = HL->hits[i].sequence_id;
      if (seqid != seqid0)
	{
	  seqid0 = seqid;
	  if (KWI->seq2kw1[seqid] != 0)
	    {
	      kw_hits++;
	      for (j=0; j < KWI->seq2kw1[seqid]; j++)
		{
		  wc_tmp = HL->kw + HL->kw_size;
		  wc_tmp->item = KWI->skw_h[KWI->seq2kw0[seqid]+j];
		  wc_tmp->count = 0;
		  wc_tmp->score = 0.0;
		  tree_item = (WRD_CNTR **) avl_probe (tree, wc_tmp);	  
		  (*tree_item)->count++;
		  if ((*tree_item) == wc_tmp)
		    {
		      HL->kw_size++;	
		      WRD_CNTR_add(&HL->kw, HL->kw_size, &HL->kw_max);
		    }
		}
	    }
	}
    }

  /* Calculate keyword scores */
  s = log((double) KWI->nsk) - log((double) kw_hits);

  for (i=0; i < HL->kw_size; i++)
    {
      HL->kw[i].score = s + log((double) HL->kw[i].count)
	- log((double) KWI->kwf[HL->kw[i].item]);
      HL->kw[i].score /= log(2.0);
    }

  /* Assign scores to hits */
  for (i=0; i < HL->actual_seqs_hits; i++)
    {
      seqid = HL->hits[i].sequence_id;
      if (KWI->seq2kw1[seqid] != 0)
	{
	  HL->hits[i].kw_score = 0.0;
	  for (j=0; j < KWI->seq2kw1[seqid]; j++)
	    {
	      wc_tmp1.item = KWI->skw_h[KWI->seq2kw0[seqid]+j];
	      tree_item = (WRD_CNTR **) avl_probe (tree, &wc_tmp1);	  
	      HL->hits[i].kw_score += (*tree_item)->score;
	    }
	}
      else
	{
	  HL->hits[i].kw_score = 0.0;
	}
    }

  /* Clean up tree */
  avl_destroy (tree, NULL);

  /* Sort HL->kw and HL->oc */
  /* Perhaps this can be done directly using avl-tree but this seems
     easier right now. */

  WRD_CNTR_sort_score(HL->kw, HL->kw_size);
}













#endif



