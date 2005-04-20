#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "misclib.h"

/********************************************************************/ 
/*                                                                  */
/*                     Exception handling                           */ 
/*                                                                  */
/********************************************************************/ 

EXCEPTION *FSexcept(int code, const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);

  FS_EXCEPTION->code = code;
  vsnprintf(FS_EXCEPTION->msg, MSG_SIZE, fmt, ap);
  va_end(ap);
  return FS_EXCEPTION;
}

  
/********************************************************************/ 
/*                                                                  */
/*      Error checking wrappers to malloc, calloc and realloc       */ 
/*                                                                  */
/********************************************************************/ 

void *callocec(long int number, long int size)
{
  void *memp;

  if (number < 0)
    Throw FSexcept(NEG_MEM_REQ, "Negative memory request.");
  memp=calloc(number, size);
  if (memp==NULL)
    Throw FSexcept(NO_MEM, "No Memory!");
  return(memp);
}

void *mallocec(long int size)
{
  void *memp;

  if (size < 0)
    Throw FSexcept(NEG_MEM_REQ, "Negative memory request.");
  memp=malloc(size);
  if (memp==NULL)
    Throw FSexcept(NO_MEM, "No Memory!");
  return(memp);
}

void *reallocec(void *pt, long int size)
{
  void *memp;

  if (size < 0)
    Throw FSexcept(NEG_MEM_REQ, "Negative memory request.");
  memp=realloc(pt, size);
  if (memp==NULL)
    Throw FSexcept(NO_MEM, "No Memory!");
  return(memp);
}

/********************************************************************/ 
/*                                                                  */
/*                    Unix path processing                          */ 
/*                                                                  */
/********************************************************************/ 

void path_split(const char *p, char **h, char **t)
{
  /* Equivalent to Python os.path.split except that we don't strip 
     multiple '/' 
  */

  int L = strlen(p);
  const char *s = p + L - 1;

  while (*s != '/' && s >= p) s--;

  if (s >= p) {
    /* Using calloc thus guaranteeing null-termination */
    *t = callocec(1, p + L - s);
    memcpy(*t, s+1, p + L - s - 1);

    if (s == p) { /* path starts with '/' */
      *h = callocec(1, 2);
      **h = '/';
    }
    else {
      *h = callocec(1, s - p + 1);
      memcpy(*h, p, s - p);
    }
  }
  else if (s < p) {
    *t = strdup(p);
    *h = callocec(1, 1);
  }
}

char *path_join(const char *h, const char *t)
{
  /* Equivalent to Python os.path.join */

  char *p;
  int m = strlen(h);
  int n = strlen(t);

  if (*t == '/')
    p = strdup(t);
  else if (m == 0 || *(h+m-1) == '/') {
    p = callocec(1, m+n+1);
    memcpy(p, h, m);
    memcpy(p+m, t, n);
  }
  else {
    p = callocec(1, m+n+2);
    memcpy(p, h, m);
    *(p+m) = '/';
    memcpy(p+m+1, t, n);
  }
  return p;
}

/********************************************************************/ 
/*                                                                  */
/*                    String I/O                                    */ 
/*                                                                  */
/********************************************************************/ 

void fwrite_string(char *s, FILE *fp) 
{
  int len = strlen(s) + 1;    
  fwrite(&len, sizeof(int), 1, fp);           
  fwrite(s, sizeof(char), len, fp);           
} 

void fread_string(char **s, FILE *fp)
{
  int len;                                        
  fread(&len, sizeof(int), 1, fp);            
  *s = mallocec(len);                         
  fread(*s, sizeof(char), len, fp);           
}


/********************************************************************/ 
/*                                                                  */
/*                Comparison functions for qsort                    */ 
/*                                                                  */
/********************************************************************/ 

int compare_int(const void *M1, const void *M2)
{
  const int *i1 = M1;
  const int *i2 = M2;
  if (*i1 > *i2)
    return 1;
  else if (*i1 < *i2)
    return -1;
  else
    return 0;
}

int compare_char(const void *M1, const void *M2)
{
  const char *i1 = M1;
  const char *i2 = M2;
  if (*i1 > *i2)
    return 1;
  else if (*i1 < *i2)
    return -1;
  else
    return 0;
}

int compare_dbl(const void *M1, const void *M2)
{
  const double *i1 = M1;
  const double *i2 = M2;
  if (*i1 > *i2)
    return 1;
  else if (*i1 < *i2)
    return -1;
  else
    return 0;
}

/********************************************************************/ 
/*                                                                  */
/*                Check validity of strings                         */ 
/*                                                                  */
/********************************************************************/ 

int check_word_alphabet(const char *word, int word_len,
			const char *alphabet, int alphabet_len)
{
  /* Checks if all letters from word belong to alphabet.
     alphabet must be sorted using compare_char */
  int i;

  for (i=0; i < word_len; i++, word++) 
    if (bsearch(word, alphabet, alphabet_len, 1, 
		compare_char) == NULL)
      return 0;

  return 1;
}

/********************************************************************/ 
/*                                                                  */
/*                Print into buffer that can grow                   */ 
/*                                                                  */
/********************************************************************/ 

#define ABSPRINTF_BUF 256
void absprintf(char **buf, int *size, int *len, const char *fmt, ...)
{
  va_list ap;
  int r; 
  int n;

  if (*buf == NULL) {
    *buf = mallocec(ABSPRINTF_BUF);
    *size = ABSPRINTF_BUF;
    *len = 0;
  }
  r = *size - *len;
  va_start(ap, fmt);
  n = vsnprintf(*buf+*len, r, fmt, ap);
  if (n >= r) {
    *size += max(ABSPRINTF_BUF, n+1-r);
    *buf = reallocec(*buf, *size);
    vsnprintf(*buf+*len, n+1, fmt, ap);
  }
  *len += n;
  va_end(ap);
} 
