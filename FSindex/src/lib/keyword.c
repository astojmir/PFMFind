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


