
#ifndef _MISCLIB_H
#define _MISCLIB_H

#include <stdio.h>
#include <stdarg.h>
#include "cexcept.h"

#ifdef USE_MPATROL
#include <mpatrol.h>
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

typedef unsigned int UINT_t;
typedef signed int SINT_t;

/* Exception handling */

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








/* Some numerical typedefs and functions from gamerf package and
   otherwise */

#ifndef DCOMPLEX
struct dcomplex_ {
    double re;
    double im;
};
#define DCOMPLEX struct dcomplex_
#define DREAL(x) (x).re
#define DIMAG(x) (x).im
#define DCMPLX(x,y,z) (z).re = x, (z).im = y
#endif

DCOMPLEX cdgamma(DCOMPLEX x);
double dcbrt(double x);
double derf(double x);
double derfc(double x);
double dgamma(double x);
double dierfc(double y); 
double dlgamma(double x);

double dibeta(double p, double a, double b);

/* Statistical functions - my own */
void discrete_meanvar(int n, SSINT *values, double *probs, 
		      double *mean, double *var); 

/* Memory managment functions - error checking */

char *newmem(long int number, long int size);
void *callocec(long int number, long int size);
void *mallocec(long int size);
void *reallocec(void *pt, long int size);

/* File functions - error checking */




/* Auxillaries */
int max(int a, int b);
int min(int a, int b);

void printbar(FILE *outstream, ULINT cntr, ULINT dispc, 
	      USINT dpr);

int compare_int(const void *M1, const void *M2);
int compare_dbl(const void *M1, const void *M2);


/* Directory names */
int split_base_dir(const char *full_name, char **basename, 
		   char **dirname);
int cat_base_dir(char **full_name, const char *basename, 
		 const char *dirname);
 
/* XML tag prinitng */
void xml_open_tag(FILE *stream, const char *tag, int indent,
		  int eol);


typedef struct
  {
    SLINT no_bins;
    SLINT min_bin;
    ULINT *freq;
    ULINT attainable_pts;
    ULINT unattainable_pts;
} SEED_HIST;

SEED_HIST *seedhist_init(SLINT min_bin, SLINT max_bin);
void seedhist_init1(SEED_HIST *hist, SLINT min_bin, SLINT max_bin);
int seedhist_clear(SEED_HIST *hist);
int seedhist_clear1(SEED_HIST *hist);
void seedhist_add_acount(SEED_HIST *hist, SLINT bin);
void seedhist_add_ucount(SEED_HIST *hist);
int seedhist_merge(SEED_HIST *hist1, SEED_HIST *hist2);
int seedhist_read(SEED_HIST *hist, FILE *stream);
int seedhist_write(SEED_HIST *hist, FILE *stream);
int seedhist_print(SEED_HIST *hist, FILE *stream);

        
/*#define seedhist_add_acount(hist, bin) \
        ((hist)->freq[(bin)]++), \
        ((hist)->attainable_pts++)*/

#define seedhist_add_ucount(hist) \
        ((hist)->unattainable_pts++)




/* Getline implementation  from mailutils */
extern int getline(char **_lineptr, size_t *_n, FILE *_stream);

extern int getdelim(char **_lineptr, size_t *_n, int _delimiter, FILE *_stream);

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS                    ***/ 
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

#define xml_open_tag(fp, tag, ind, eol) do {       \
  fprintf((fp),"%*.*s<%s>",(ind),(ind),"",(tag));  \
  if ((eol)) fprintf((fp), "\n");                  \
  } while(0)                                       \

#define xml_close_tag(fp, tag, ind) do {           \
  fprintf((fp),"%*.*s</%s>\n",(ind),(ind),"",(tag));\
  } while(0)                                       \



#if 0
#define e_malloc(ptr, size) do {                  \
  if ((size) < 0)                                 \
    return NEG_MEM_REQ;                           \
  if (((ptr) = malloc(size)) == NULL)             \
    return NO_MEM;                                \
  } while(0)

#define e_calloc(ptr, num, size) do {             \
  if ((size)*(num) < 0)                           \
    return NEG_MEM_REQ;                           \
  if ((ptr = calloc((num), (size))) == NULL)      \
    return NO_MEM;                                \
  } while(0)

#define e_realloc(ptr, size) do {                 \
  if ((size) < 0)                                 \
    return NEG_MEM_REQ;                           \
  if (((ptr) = realloc((ptr), (size))) == NULL)   \
    return NO_MEM;                                \
  } while(0)

#define e_free(ptr) do {                          \
  free((ptr));                                    \
  ptr = NULL;                                     \
  } while(0)

#define e_fopen(stream, path, mode) do {          \
  if (((stream) = fopen((path), (mode))) == NULL) \
    return FOPEN_ERR;                             \
  } while(0)

#define e_fclose(stream) do {                     \
  if (fclose((stream)) == EOF)                    \
    return FCLOSE_ERR;                            \
  } while(0)

#define check_err(status) do {                    \
  if (status)                                     \
    return status;                                \
  } while(0)
#endif


#ifdef  MY_INLINE

#if 0
MY_INLINE
int max(int a, int b)
{
  return a > b ? a : b;
}

MY_INLINE
int min(int a, int b)
{
  return a < b ? a : b;
}

MY_INLINE
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
#endif

#endif /* MISCLIB_INLINE */   


#endif
