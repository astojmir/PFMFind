#include "misclib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TOTAL_MEM
#ifdef TOTAL_MEM
static unsigned long int total_mem = 0; 
#endif

void *callocec(long int number, long int size)
{
  void *memp;

  if (number < 0)
    { 
      fprintf(stderr, "Negative memory request \n");
      exit(-1);
    }
  memp=calloc(number, size);
  if (memp==NULL)
    { 
      fprintf(stderr,"Memory Full! %ld requested.\n", number*size);
#ifdef TOTAL_MEM
      fprintf(stderr,"%ld bytes allocated.\n", total_mem);
#endif
      exit(-1);
    }
#ifdef TOTAL_MEM
  total_mem += (number * size);
#endif
  return(memp);
}

void *mallocec(long int size)
{
  void *memp;

  if (size < 0)
    { 
      fprintf(stderr, "Negative memory request \n");
      exit(-1);
    }  
  memp=malloc(size);
  if (memp==NULL)
    { 
      fprintf(stderr,"Memory Full! %ld requested.\n", size);
#ifdef TOTAL_MEM
      fprintf(stderr,"%ld bytes allocated.\n", total_mem);
#endif
      exit(-1);
    }
#ifdef TOTAL_MEM
  total_mem += size;
#endif
  return(memp);
}

void *reallocec(void *pt, long int size)
{
  void *memp;
  if (size < 0)
    { 
      fprintf(stderr, "Negative memory request \n");
      exit(-1);
    }  
  memp=realloc(pt, size);
  if (memp==NULL)
    { 
      fprintf(stderr,"realloc error! %ld requested.\n", size);
      exit(-1);
    }
  return(memp);
}

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
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS                    ***/ 
/***                     (library copies)                         ***/
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

int max(int a, int b)
{
  return a > b ? a : b;
}

int min(int a, int b)
{
  return a < b ? a : b;
}

void printbar(FILE *outstream, ULINT cntr, ULINT dispc, 
	      USINT dpr)  
{
  char dot = '.';

  if((cntr % dispc) == 0)
    {
      /* print '*' as every 10th dot */ 
      if ((cntr % (10*dispc)) == 0) dot = '*'; 
      fputc(dot, outstream);
      if ((cntr % (dpr*dispc)) == 0) 
        fprintf(outstream, " %ld\n", cntr);
      fflush(outstream);
    }
}         

int compare_int(const void *M1, const void *M2)
{
  int *B1;
  int *B2;

  B1 = (int *)M1;
  B2 = (int *)M2;

  return *B1 - *B2;
}

int compare_dbl(const void *M1, const void *M2)
{
  double *B1;
  double *B2;

  B1 = (double *)M1;
  B2 = (double *)M2;

  if (*B1 > *B2)
    return 1;
  else if (*B1 < *B2)
    return -1;
  else
    return 0;
}
