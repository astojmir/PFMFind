
#ifndef _MISCLIB_H
#define _MISCLIB_H

#include <stdio.h>

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

/* Memory managment functions */

char *newmem(long int number, long int size);
void *callocec(long int number, long int size);
void *mallocec(long int size);
void *reallocec(void *pt, long int size);


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
