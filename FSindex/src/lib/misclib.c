#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "misclib.h"

EXCEPTION *FSexcept(int code, const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);

  FS_EXCEPTION->code = code;
  vsnprintf(FS_EXCEPTION->msg, MSG_SIZE, fmt, ap);
  va_end(ap);
  return FS_EXCEPTION;
}
  
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

/* ************* seedhist module *************** */


SEED_HIST *seedhist_init(SLINT min_bin, SLINT max_bin)
{
  SEED_HIST *hist;
  if (min_bin > max_bin) return NULL;
  hist = callocec(1, sizeof(SEED_HIST));
  hist->min_bin = min_bin;
  hist->no_bins = 1 + max_bin - min_bin;
  hist->freq = callocec(hist->no_bins, sizeof(ULINT));
  return hist;
}

void seedhist_init1(SEED_HIST *hist, SLINT min_bin, SLINT max_bin)
{
  if (min_bin > max_bin) 
    {
#if DEBUG > 0
      printf("Invalid seedhist_init1 range.\n");
#endif
      return;
    }
  hist->min_bin = min_bin;
  hist->no_bins = 1 + max_bin - min_bin;
  hist->freq = (ULINT *) callocec(hist->no_bins, sizeof(ULINT));
}


int seedhist_clear(SEED_HIST *hist)
{
  free(hist->freq);
  free(hist);
  return 1;
}

int seedhist_clear1(SEED_HIST *hist)
{
  free(hist->freq);
  return 1;
}

void seedhist_add_acount(SEED_HIST *hist, SLINT bin)
{
  SLINT old_no_bins;
  SLINT dsize;
  ULINT *old_freq;

  if (bin < hist->min_bin)
    {
#if DEBUG > 0
      printf("Smaller bin!\n");
#endif
      old_no_bins = hist->no_bins;
      dsize = hist->min_bin - bin;
      hist->min_bin = bin;
      hist->no_bins += dsize;
      old_freq = hist->freq;
      hist->freq = callocec(hist->no_bins, sizeof(ULINT));
      memcpy(hist->freq+dsize, old_freq, old_no_bins*sizeof(ULINT)); 
      free(old_freq);
    }
  else if (bin > (hist->no_bins+hist->min_bin-1))
    {
#if DEBUG > 0
      printf("Larger bin!\n");
#endif
      dsize = bin - hist->no_bins - hist->min_bin +1;
      old_no_bins = hist->no_bins;
      old_freq = hist->freq;
      hist->no_bins += dsize;
      hist->freq = reallocec(hist->freq, hist->no_bins*sizeof(ULINT));
      memset(hist->freq+old_no_bins, 0, dsize*sizeof(ULINT));
    }
  hist->freq[bin - hist->min_bin]++;
  hist->attainable_pts++;
  return ;
}


int seedhist_merge(SEED_HIST *hist1, SEED_HIST *hist2)
{
  ULINT i;
  if ((hist1->no_bins != hist2->no_bins) || 
      (hist1->min_bin != hist2->min_bin))
    return 0;

  hist1->attainable_pts += hist2->attainable_pts;
  hist1->unattainable_pts += hist2->unattainable_pts;
  for (i=0; i < hist1->no_bins; i++)
    hist1->freq[i] += hist2->freq[i];
  return 1;
}

int seedhist_read(SEED_HIST *hist, FILE *stream)
{
  fread(&hist->no_bins, sizeof(SLINT), 1, stream);
  fread(&hist->min_bin, sizeof(SLINT), 1, stream);
  fread(&hist->attainable_pts, sizeof(ULINT), 1, stream);
  fread(&hist->unattainable_pts, sizeof(ULINT), 1, stream);
  hist->freq = callocec(hist->no_bins, sizeof(ULINT));
  fread(hist->freq, sizeof(ULINT), hist->no_bins, stream);
  return 1;
}

int seedhist_write(SEED_HIST *hist, FILE *stream)
{
  fwrite(&hist->no_bins, sizeof(SLINT), 1, stream);
  fwrite(&hist->min_bin, sizeof(SLINT), 1, stream);
  fwrite(&hist->attainable_pts, sizeof(ULINT), 1, stream);
  fwrite(&hist->unattainable_pts, sizeof(ULINT), 1, stream);
  fwrite(hist->freq, sizeof(ULINT), hist->no_bins, stream);
  return 1;
}

int seedhist_print(SEED_HIST *hist, FILE *stream)
{
  ULINT i;
  ULINT min0 = 0;
  ULINT max0 = hist->no_bins -1;
  ULINT counted_seq = 0;

  while(hist->freq[min0] == 0 && min0 < hist->no_bins)
    min0++;
  while(hist->freq[max0] == 0 && max0 > 0)
    max0--;
#if DEBUG > 5
  printf("min0 = %ld\n", min0);
  printf("max0 = %ld\n", max0);
  fflush(stdout);
#endif
  for (i=min0; i <= max0; i++)
    counted_seq += hist->freq[i];

  fprintf(stream, "Points: %ld %ld\n", hist->attainable_pts,
	  hist->unattainable_pts);
#if DEBUG > 5
  fprintf(stream, "Counted points: %ld\n", counted_seq);
#endif
  fprintf(stream, "Values: %ld\n", max0+1-min0);
  fprintf(stream, "%12s %12s\n", "Value", "Frequency");

  for (i=min0; i <= max0; i++)
    fprintf(stream, "%12ld %12ld\n", hist->min_bin+i, hist->freq[i]);
  fprintf(stream, "\n");
  return 1;
}




/* ************* End of seedhist module *********** */


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

int
getline (char **lineptr, size_t *n, FILE *stream)
{
  return getdelim (lineptr, n, '\n', stream);
}

#ifndef HAVE_GETDELIM

/* Default value for line length.  */
static const int line_size = 128;

int
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

