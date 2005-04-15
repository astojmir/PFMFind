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


int split_base_dir(const char *full_name, char **basename, 
		   char **dirname)
{
  int len = strlen(full_name);
  
  int i;

  for (i = len-1; i > 0 && full_name[i] != '/'; i--);

  if (i == 0)
    {
      *basename = mallocec(len);
      *dirname = mallocec(2);
      if (full_name[i] == '/')
	{
	  (*dirname)[0] = '/';
	  (*dirname)[1] = '\0';
	  strncpy(*basename, full_name+1, len);
	}
      else 
	{
	  (*dirname)[0] = '\0';
	  strncpy(*basename, full_name, len);
	}
    }
  else
    {
      *basename = mallocec(len-i+1);
      *dirname = mallocec(i+1);
      strncpy(*basename, full_name+i+1, len-i+1);
      strncpy(*dirname, full_name, i);
      (*dirname)[i] = '\0';
    }      
  return 1;      
}

int cat_base_dir(char **full_name, const char *basename, 
		 const char *dirname)
{
  int base_len = strlen(basename);
  int dir_len = strlen(dirname);

  if (dir_len == 1 && dirname[0] == '/')
    {
      *full_name = callocec(base_len + 2, 1);
      *full_name[0] = '/';
      memcpy(*full_name+1, basename, base_len);
    }
  else if (dir_len > 0)
    {
      *full_name = callocec(dir_len + base_len + 2, 1);
      memcpy(*full_name, dirname, dir_len);
      (*full_name)[dir_len] = '/';
      memcpy(*full_name + dir_len + 1, basename, base_len);
      (*full_name)[dir_len + base_len + 1] = '\0';
    }
  else
    {
      *full_name = callocec(base_len + 1, 1);
      memcpy(*full_name, basename, base_len);
    }

  return 1;
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
