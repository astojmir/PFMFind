#include "misclib.h"
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void discrete_meanvar(int n, SSINT *values, double *probs, 
		      double *mean, double *var)
{
  int i = 0;
  double T = 0.0;
  *mean = 0.0;
  *var = 0.0;

  for (i=0; i < n; i++)
    {
      T = values[i] * probs[i];
      *mean += T;
      T *= values[i];
      *var += T;
    }
  *var = *var - (*mean * *mean); 
}


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


/* Getline implementation */

/* GNU mailutils - a suite of utilities for electronic mail
   Copyright (C) 1999, 2000, 2001, 2002 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Library Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Library General Public License for more details.

   You should have received a copy of the GNU Library General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */

/* First implementation by Alain Magloire */

ssize_t
getline (char **lineptr, size_t *n, FILE *stream)
{
  return getdelim (lineptr, n, '\n', stream);
}

#ifndef HAVE_GETDELIM

/* Default value for line length.  */
static const int line_size = 128;

ssize_t
getdelim (char **lineptr, size_t *n, int delim, FILE *stream)
{
  int indx = 0;
  int c;

  /* Sanity checks.  */
  if (lineptr == NULL || n == NULL || stream == NULL)
    return -1;

  /* Allocate the line the first time.  */
  if (*lineptr == NULL)
    {
      *lineptr = malloc (line_size);
      if (*lineptr == NULL)
	return -1;
      *n = line_size;
    }

  while ((c = getc (stream)) != EOF)
    {
      /* Check if more memory is needed.  */
      if (indx >= *n)
	{
	  *lineptr = realloc (*lineptr, *n + line_size);
	  if (*lineptr == NULL)
	    return -1;
	  *n += line_size;
	}

      /* Push the result in the line.  */
      (*lineptr)[indx++] = c;

      /* Bail out.  */
      if (c == delim)
	break;
    }

  /* Make room for the null character.  */
  if (indx >= *n)
    {
      *lineptr = realloc (*lineptr, *n + line_size);
      if (*lineptr == NULL)
       return -1;
      *n += line_size;
    }

  /* Null terminate the buffer.  */
  (*lineptr)[indx++] = 0;

  /* The last line may not have the delimiter, we have to
   * return what we got and the error will be seen on the
   * next iteration.  */
  return (c == EOF && (indx - 1) == 0) ? -1 : indx - 1;
}

#endif /* HAVE_GETDELIM */

