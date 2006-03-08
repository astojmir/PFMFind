/*
 * Copyright (C) 2004-2006 Victoria University of Wellington
 *
 * This file is part of the PFMFind module.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2,
 * or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */



#ifndef _MISCLIB_H
#define _MISCLIB_H

#include <stdio.h>
#include <stdarg.h>
#include "cexcept.h"

#ifdef USE_MPATROL
#include <mpatrol.h>
#endif

#ifdef USE_DMALLOC
#include <dmalloc.h>
#endif

#ifdef GCC_INLINE
#define MY_INLINE extern inline
#else
#define MY_INLINE static inline
#endif

/* Some typedefs */

typedef unsigned long int ULINT;
typedef signed long int SLINT;
typedef unsigned short int USINT;
typedef signed short int SSINT;


/* Exception handling - using cexept */

#define MSG_SIZE 512
typedef struct 
{
  int code;
  char msg[MSG_SIZE+1];
} EXCEPTION;

extern EXCEPTION FSexcept_array[1];
#define FS_EXCEPTION (FSexcept_array)
EXCEPTION *FSexcept(int code, const char *fmt, ...); 

define_exception_type(EXCEPTION *);
extern struct exception_context the_exception_context[1];

#define NO_ERR 0
#define NO_MEM 1
#define NEG_MEM_REQ 2
#define FOPEN_ERR 3
#define FCLOSE_ERR 4
#define BAD_ARGS 5
#define EOF_REACHED 6
#define GETLINE_ERR 7
#define INDEX_OUT_OF_RANGE 8
#define INCONSITENT_DATA 9
#define IO_ERR 10
#define RUNTIME_ERR 11


/* Printing into buffers */

void absprintf(char **buf, int *size, int *len, const char *fmt, ...);


/* Memory managment functions - error checking */

void *callocec(long int number, long int size);
void *mallocec(long int size);
void *reallocec(void *pt, long int size);


/* Auxillaries */
int max(int a, int b);
int min(int a, int b);

/* I/O */

void fwrite_string(char *s, FILE *fp);
void fread_string(char **s, FILE *fp);

void swrite_int32(void *ptr, size_t nmemb, unsigned char *s);
void sread_int32(void *ptr, size_t nmemb, unsigned char *s);
void fwrite_int32(void *ptr, size_t nmemb, FILE *stream);
void fwrite_int32(void *ptr, size_t nmemb, FILE *stream);

/* Progress bar */

void printbar(FILE *outstream, ULINT cntr, ULINT dispc, 
	      USINT dpr);


/* Comparison functions */

int compare_int(const void *M1, const void *M2);
int compare_char(const void *M1, const void *M2);
int compare_dbl(const void *M1, const void *M2);

int check_word_alphabet(const char *word, int word_len,
			const char *alphabet, int alphabet_len);

/* Directory names - almost equivalent to corresponding Python
   os.path functions */

void path_split(const char *path, char **head, char **tail);
char *path_join(const char *head, const char *tail);



/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS (Macros)           ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    


#define max(a, b)  ((a) > (b) ? (a) : (b))
#define min(a, b)  ((a) < (b) ? (a) : (b))

#define printbar(outstream, cntr, dispc, dpr) do { \
  char dot = '.';                                  \
  if(((cntr) % (dispc)) == 0)                      \
    {                                              \
      if (((cntr) % (10*(dispc))) == 0) dot = '*'; \
      fputc(dot, (outstream));                     \
      if (((cntr) % ((dpr)*(dispc))) == 0)         \
        {                                          \
          fprintf((outstream), " %ld\n", (cntr));  \
        }                                          \
      fflush((outstream));                         \
    }                                              \
  } while(0)        

#endif /* _MISCLIB_H */
