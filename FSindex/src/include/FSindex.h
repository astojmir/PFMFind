#ifndef _FS_INDEX_H
#define _FS_INDEX_H

#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fastadb.h"
#include "partition.h"
#include "smatrix.h"
#include "pmatrix.h"
#include "hit_list.h"
#include "keyword.h"
#ifdef USE_MPATROL
#include <mpatrol.h>
#endif


/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               PROTOTYPES                                     ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

/********************************************************************/    
/*                                                                  */
/*                     UFRAG                                        */ 
/*                                                                  */
/********************************************************************/    

typedef union
{
  SEQ_index_t *p; /* Pointer to the array of fragments */
  SEQ_index_t s;  /* Only fragment */
} FRAG_HNDL; 

typedef struct
{
  int n;
  FRAG_HNDL frag; 
} UFRAG;

/********************************************************************/    
/*                                                                  */
/*                     Filtering parameters                         */ 
/*                                                                  */
/********************************************************************/    

typedef struct
{
  double ZE0;
  double ZE1;
  double CR0;
  double CR1;
  double KW0;
} FPARAMS;



/********************************************************************/    
/*                                                                  */
/*                     FS_INDEX module                              */ 
/*                                                                  */
/********************************************************************/    

#if 0
extern int FS_INDEX_VERBOSE;
extern int FS_INDEX_PRINT_BAR;
#endif

/* FSindex variables - essentially constant throughout all searches
   with the same index .*/ 

typedef struct
{
  char *db_name;           /* Name of the database                  */
  char *index_name;        /* Name of the index                     */
  SEQUENCE_DB *s_db;       /* Pointer to the fasta database         */
  double dtime;            /* Creation time                         */
  char *alphabet;          /* Amino acid alphabet - partitioned     */
  char separator;          /* Partitioning character                */
  int K;                   /* Number of letters in reduced alphabet */
  FS_PARTITION_t *ptable;  /* Partition table                       */
  ULINT *KK;               /* KK[i] = K ^ i, i = 0, 1 ... frag_len  */
  int *K_modtable;         /* K_modtable[i] = i % K for i < 2*K     */

  int m;                   /* Length of indexed fragments           */
  ULINT db_no_frags;       /* Number of fragments in full database  */
  ULINT no_bins;           /* Number of bins = K^frag_len           */
  ULINT no_seqs;           /* Number of indexed fragments           */
  ULINT no_useqs;          /* Number of indexed unique fragments    */

  ULINT shist_len;         /* Length of shist                       */
  ULINT *shist;            /* Histogram of bin sizes                */
  ULINT uhist_len;         /* Length of uhist                       */
  ULINT *uhist;            /* Histogram of bin sizes - unique seqs  */

  ULINT *bin_size;


  ULINT *max_bin_size;
  SEQ_index_t **bin;
  SEQ_index_t binL;        /* Largest bin                           */
  SEQ_index_t *heap;
  ULINT *u_size;
  int **u;
  int *u_heap;
  SEQ_index_t uL;          /* Largest bin - unique seqs             */

} FSINDX;



FSINDX *FS_INDEX_create(const char *database, ULINT flen,
			const char *abet, const char sepchar, 
			int skip);
void FS_INDEX_destroy(FSINDX *FSI);

int FS_INDEX_save(FSINDX *FSI, const char *filename);
FSINDX *FS_INDEX_load(const char *filename);

void FS_INDEX_print_stats(FSINDX *FSI, FILE *stream, int options); 

void FS_INDEX_print_bin(FSINDX *FSI, BIOSEQ *query, FILE *stream,
			int options);

/* Search functions */

HIT_LIST_t *FSINDX_rng_srch(FSINDX *FSI, BIOSEQ *query, SCORE_MATRIX_t *D,
			    int d0, HIT_LIST_t *HL);

HIT_LIST_t *FSINDX_kNN_srch(FSINDX *FSI, BIOSEQ *query, SCORE_MATRIX_t *D,
			    int kNN, HIT_LIST_t *HL);

HIT_LIST_t *FSINDX_prof_rng_srch(FSINDX *FSI, BIOSEQ *query, 
				 SCORE_MATRIX_t *D, int s0, 
				 double lambda, double A,
				 const char *freq_filename,
				 int iters, int s1,
				 HIT_LIST_t *HL);

HIT_LIST_t *FSINDX_prof_kNN_srch(FSINDX *FSI, BIOSEQ *query, 
				 SCORE_MATRIX_t *D, int kNN, 
				 double A,
				 const char *freq_filename,
				 int iters, 
				 HIT_LIST_t *HL);


/* Experiments */
HIT_LIST_t *SSCAN_QD_search(SEQUENCE_DB *s_db, const char *matrix, 
			    BIOSEQ *query, ULINT D_cutoff);
#if 0
int FS_INDEX_has_neighbour(BIOSEQ *query, SCORE_MATRIX_t *D0, 
			   int cutoff);
#endif
int SSCAN_has_neighbour(SEQUENCE_DB *s_db,  SCORE_MATRIX_t *D, 
			FS_PARTITION_t *ptable,
			BIOSEQ *query, ULINT D_cutoff);




/* Access functions */
ULINT FS_INDEX_get_bin_size(FSINDX *FSI, ULINT bin);
ULINT FS_INDEX_get_unique_bin_size(FSINDX *FSI, ULINT bin);
BIOSEQ *FS_INDEX_get_seq(FSINDX *FSI, ULINT bin, ULINT seq);

#define FS_INDEX_get_ptable(FSI) \
        ((FSI)->ptable)

#define FS_INDEX_get_database(FSI) \
        ((FSI)->s_db)

#define FS_INDEX_get_frag_len(FSI) \
        ((FSI)->m)

#define FS_INDEX_get_no_bins(FSI) \
        ((FSI)->no_bins)

#endif   

