/********************************************************************/    
/*                                                                  */
/*                    FS_PARTITION module                           */ 
/*                                                                  */
/********************************************************************/    
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "partition.h"
#include "misclib.h"

/* Main constructor */
FS_TABLE *FS_TABLE_init(const char **sepn, int len)
{
  int i, j, l, m, p;
  const char *c;
  char *abet;
  ULINT KK = 1;
  FS_TABLE *ptable; 

  /* Make sure that alphabet is in the right format and
     alphabets at all positions are the same. We do it
     here rather than later so that we can abort without
     allocating anything. The code following the check 
     assumes everything is as it should be. */

  /* Check format */
  if (len == 0) return NULL;
  for (i=len-1; i >= 0; i--) {
    l = strlen(sepn[i]);
    if (l == 0 || sepn[i][0] == '#' ||  sepn[i][l-1] == '#')
      return NULL;
    for (j=0; j < l-1; j++) {
      if (sepn[i][j] == '#' && sepn[i][j+1] == '#')
	return NULL;
    }
  }

  /* Get alphabet at the first position */
  abet = callocec(l, 1);
  for (c = sepn[0], l=0; *c; c++)
    if (*c != '#') abet[l++] = *c;
  /* Sort */
  qsort(abet, l, 1, compare_char);

  /* Check that each letter appears only once */
  for (j=1; j < l; j++) 
    if (abet[j] == abet[j-1]) {
      free(abet);
      return NULL;
    }

  /* Check if other positions have the same alphabet as the first. 
     It is sufficient to check that 
     - the lengths of alphabets match, and
     - each letter in abet matches one and only one letter 
     in sepn[i]. */

  for (i=1; i < len; i++) {
    for (j=0; j < l; j++) {
      for (c=sepn[i], p=0; *c; c++) 
	if (*c == abet[j]) p++;
      if (p != 1) {
	free(abet);
	return NULL;
      }
    }
    for(m=0, c = sepn[i]; *c; c++)
      if (*c == '#') continue;
      else m++;
    if (m != l) {
      free(abet);
      return NULL;
    }
  }

  ptable  = callocec(1, sizeof(FS_TABLE));
  ptable->len = len;
  ptable->sepn = callocec(len, sizeof(char *));
  ptable->pttn = callocec(len, sizeof(int *));
  ptable->hash = callocec(len, sizeof(int *));
  ptable->K = mallocec(len * sizeof(int));
  ptable->nbpt = mallocec(len * sizeof(int));
  ptable->rank = callocec(len, sizeof(int *));

  ptable->alphabet = abet;
  /* Get the unsorted alphabet orders */
  /* 0 is reserved for the end of sequence data ('_');
     1 to A_SIZE - l is reserved for any unreckognised characters +
     other '\0' */
  for (i=0; i < len; i++) {
    ptable->rank[i] = callocec(A_SIZE, sizeof(int));
    for (c = sepn[i], j= A_SIZE - l; *c; c++)
      if (*c != '#') 
	ptable->rank[i][*c & A_SIZE_MASK] = j++;
    for (j=0, m=1; j < A_SIZE; j++) {
      if (j == ('_' & A_SIZE_MASK)) continue;
      if (!ptable->rank[i][j]) ptable->rank[i][j] = m++;
    }
  }

  /* Check if order is the same at each position */
  ptable->same_order = 1;
  for (i=1; i < len; i++) {
    for (j=0; j < A_SIZE; j++) 
      if (ptable->rank[i][j] != ptable->rank[0][j]) {
	ptable->same_order = 0;
	break;
      }
    if (!ptable->same_order) break;
  }


  for (i=len-1; i >= 0; i--) {
    c = sepn[i];
    ptable->sepn[i] = strdup(c);
    ptable->pttn[i] = mallocec(A_SIZE * sizeof(int));
    memset(ptable->pttn[i], -1, A_SIZE * sizeof(int));
    p = 0;
    while (*c) {
      if (*c == '#') p++;
      else ptable->pttn[i][*c & A_SIZE_MASK] = p;
      c++;
    }
    ptable->K[i] = ++p;
    ptable->nbpt[i] = p-1;
    ptable->hash[i] = mallocec(p * sizeof(int));
    for (j=0; j < p; j++) {
      ptable->hash[i][j] = j*KK;
    }
    if (KK > ULONG_MAX / p) {
      FS_TABLE_del(ptable);
      return NULL;
    }
    KK *= p;
  }   
  ptable->bins = KK;
  return ptable;
}

/* Destructor */

void FS_TABLE_del(FS_TABLE *ptable)
{
  int i;
  if (ptable == NULL) return;
  free(ptable->alphabet);
  free(ptable->K);
  free(ptable->nbpt);
  for (i=0; i < ptable->len; i++) {
    free(ptable->sepn[i]);
    free(ptable->pttn[i]);
    free(ptable->hash[i]);
    free(ptable->rank[i]);
  }
  free(ptable->sepn);
  free(ptable->pttn);
  free(ptable->hash);
  free(ptable->rank);
  free(ptable);
}

void FS_TABLE_write(FS_TABLE *ptable, FILE *fp)
{
  int i;
  fwrite(&(ptable->len), sizeof(int), 1, fp);
  for (i=0; i < ptable->len; i++) 
    fwrite_string(ptable->sepn[i], fp);
}

FS_TABLE *FS_TABLE_read(FILE *fp)
{
  char **sepn;
  int i, len;
  FS_TABLE *ptable;

  fread(&len, sizeof(int), 1, fp);
  sepn = mallocec(len * sizeof(char *));
  for (i=0; i < len; i++) 
    fread_string(sepn + i, fp);
  
  ptable = FS_TABLE_init((const char **)sepn, len);

  for (i=0; i < len; i++)
    free(sepn[i]);
  free(sepn);
  return ptable;
}

char *FS_TABLE_sprint(FS_TABLE *ptable)
{
  char *buf = NULL;
  int size = 0;
  int len = 0;
  int i;
  for (i=0; i < ptable->len; i++)
    absprintf(&buf, &size, &len, "%4d %s\n", i, ptable->sepn[i]);
  return buf;
}

void FS_TABLE_fprint(FS_TABLE *ptable, FILE *fp)
{
  char *buf = FS_TABLE_sprint(ptable);
  fprintf(fp, "%s", buf);
  free(buf);
}


#define SEQ_BUF_SIZE 32
/* Pretty printing */
char *FS_SEQ_print(FS_TABLE *ptable, FS_SEQ_t FSseq)
{
  int i, l;
  int n = 0;
  int N = SEQ_BUF_SIZE;
  int q;
  char *c;
  char *s = mallocec(N);
  int *p = mallocec(ptable->len*sizeof(int));

  for (i=ptable->len-1; i >= 0; i--) {
    p[i] = FSseq % ptable->K[i];
    FSseq /= ptable->K[i];
  }
  for (i=0; i < ptable->len; i++) {
    for (q=0, c=ptable->sepn[i]; *c; c++) 
      if (*c == '#') q++;
      else if (p[i] == q) break;
    for (l=0; c[l] != '#'; l++);
    while (n+l+3 >= N) {
      N += SEQ_BUF_SIZE;
      s = reallocec(s, N);
    }
    s[n++] = '[';
    while (*c && *c != '#' )
      s[n++] = *c++;
    s[n++] = ']';
  }
  s[n] = '\0';
  free(p);
  return s;
}

int *FS_TABLE_order_convert_heap(FS_TABLE *ptable, const char *s, int len)
{
  int i, k;
  int *r;
  int amap[A_SIZE];
  int freq[A_SIZE];
  int ainv[A_SIZE];
  const char *t = s;

  r= mallocec(len*sizeof(int));
  memset(freq, 0, A_SIZE*sizeof(int));

  /* Must scan the whole dataset to determine what letters are used
     so we don't have any 'holes' in order. This is a requirement
     of sarray function */
  for (i=0; i < A_SIZE; i++) 
    ainv[ptable->rank[0][i]] = i;
  
  for (i=0; i < len; i++, t++)
    freq[*t & A_SIZE_MASK]++;

  for (k=0, i=0; i < A_SIZE; i++) {
    if (freq[ainv[i]]) amap[ainv[i]] = k++;
    else amap[ainv[i]] = -1;
  }

  for (i=0; i < len; i++, s++)
    r[i] = amap[*s & A_SIZE_MASK];
  return r;
} 
