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

extern int FS_INDEX_VERBOSE;
extern int FS_INDEX_PRINT_BAR;

/* FSindex variables - essentially constant throughout all searches
   with the same index .*/ 

typedef struct
{
  char *db_name;           /* Name of the database                  */
  char *index_name;        /* Name of the index                     */
  SEQUENCE_DB *s_db;       /* Pointer to the fasta database         */

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

void FS_INDEX_print_stats(FSINDX *FSI, FILE *stream, ULINT count, 
			  double dtime); 


/* Search functions */

HIT_LIST_t *FSINDX_rng_srch(FSINDX *FSI, BIOSEQ *query, SCORE_MATRIX_t *D,
			    int d0, HIT_LIST_t *HL, KW_INDEX *KWI, 
			    FPARAMS *fp);

HIT_LIST_t *FSINDX_kNN_srch(FSINDX *FSI, BIOSEQ *query, SCORE_MATRIX_t *D,
			    int kNN, HIT_LIST_t *HL, KW_INDEX *KWI, 
			    FPARAMS *fp);

HIT_LIST_t *FSINDX_prof_rng_srch(FSINDX *FSI, BIOSEQ *query, 
				 SCORE_MATRIX_t *D, int s0, 
				 double lambda, double A,
				 const char *freq_filename,
				 int iters, int s1,
				 HIT_LIST_t *HL,
				 KW_INDEX *KWI, FPARAMS *fp);

HIT_LIST_t *FSINDX_prof_kNN_srch(FSINDX *FSI, BIOSEQ *query, 
				 SCORE_MATRIX_t *D, int kNN, 
				 double A,
				 const char *freq_filename,
				 int iters, 
				 HIT_LIST_t *HL, 
				 KW_INDEX *KWI,
				 FPARAMS *fp);

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
#define FS_INDEX_get_ptable(FSI) \
        ((FSI)->ptable)

#define FS_INDEX_get_database(FSI) \
        ((FSI)->s_db)

#define FS_INDEX_get_frag_len(FSI) \
        ((FSI)->m)

#define FS_INDEX_get_no_bins(FSI) \
        ((FSI)->no_bins)

#define FS_INDEX_get_no_seqs(FSI, i) \
        ((FSI)->bin_size[(i)])

#define FS_INDEX_get_no_useqs(FSI, i) \
        ((FSI)->u_size[(i)])

#define FS_INDEX_retrieve_seq(FSI, i, j) \
        ((FSI)->bin[(i)][(j)])

#define FS_INDEX_get_all_seqs(FSI, i) \
        ((FSI)->bin[(i)])


#if 0













/* Main constructor */
FS_HASH_TABLE_t *FS_HASH_TABLE_create(ULINT no_bins, ULINT def_size); 

/* Destructor */
void FS_HASH_TABLE_destroy(FS_HASH_TABLE_t *HT);

/* Insertion */
void FS_HASH_TABLE_count_seq(FS_HASH_TABLE_t *HT, ULINT i,
			     SEQ_index_t seq_i);
void FS_HASH_TABLE_allocate_bins(FS_HASH_TABLE_t *HT);
void FS_HASH_TABLE_insert_seq(FS_HASH_TABLE_t *HT, ULINT i,
			      SEQ_index_t seq_i);
void FS_HASH_TABLE_add_bin_stats(FS_HASH_TABLE_t *HT,
				 SEQ_index_t seq_i);
void FS_HASH_TABLE_resize(FS_HASH_TABLE_t *HT);

/* Element access */
ULINT FS_HASH_TABLE_get_no_bins(FS_HASH_TABLE_t *HT);
ULINT FS_HASH_TABLE_get_total_seqs(FS_HASH_TABLE_t *HT);
ULINT FS_HASH_TABLE_get_no_seqs(FS_HASH_TABLE_t *HT, ULINT i);
SEQ_index_t FS_HASH_TABLE_retrieve_seq(FS_HASH_TABLE_t *HT, ULINT i, 
				       ULINT j);
SEQ_index_t *FS_HASH_TABLE_get_all_seqs(FS_HASH_TABLE_t *HT, 
					ULINT i);


/* Reading, writing to file */
int FS_HASH_TABLE_write(FS_HASH_TABLE_t *HT, FILE *stream);
FS_HASH_TABLE_t *FS_HASH_TABLE_read(FILE *stream);

/* Statistics */
void FS_HASH_TABLE_print_stats(FS_HASH_TABLE_t *HT, FILE *stream,
			       FS_PARTITION_t *FS_partition, 
			       ULINT frag_len);





typedef void FS_INDEX_process_func(FS_SEQ_t FS_neighbour);

/* These functions should convert the cutoff value to distances for
   bin search */

typedef int FS_INDEX_range_convert_func(SCORE_MATRIX_t *S0, 
					BIOSEQ *query0,
					int cutoff);	

/* Main constructor */
void FS_INDEX_create(const char *database, ULINT flen,
		     const char *abet, 
		     const char sepchar, int skip);

/* Destructor */
void FS_INDEX_destroy(void);

/* Load, Save */
int FS_INDEX_save(const char *filename);
int FS_INDEX_load(const char *filename);

/* Print */
void FS_INDEX_print_stats(FILE *stream, ULINT count, double dtime); 


/* Get / Set members */
SEQUENCE_DB *FS_INDEX_get_database(void);
const char *FS_INDEX_get_db_name(void);
FS_PARTITION_t *FS_INDEX_get_ptable(void);
int FS_INDEX_get_frag_len(void);
FS_HASH_TABLE_t *FS_INDEX_get_hash_table(void);
const char *FS_index_get_alphabet(void);

/* Main general range search */
int FS_INDEX_search(HIT_LIST_t *hit_list0, BIOSEQ *query0, 
		    SCORE_MATRIX_t *S0, SCORE_MATRIX_t *D0,
		    int cutoff, FS_INDEX_process_func *pfunc0, 
		    FS_INDEX_range_convert_func *cfunc);

/* Other Searches */
int FS_INDEX_kNN_search(HIT_LIST_t *hit_list0, BIOSEQ *query0,
			SCORE_MATRIX_t *S0, SCORE_MATRIX_t *D0,
			ULINT k);


/* Process bin functions */

FS_INDEX_process_func FS_INDEX_QD_process_bin;
FS_INDEX_process_func FS_INDEX_S_process_bin;
FS_INDEX_process_func FS_INDEX_count_only_process_bin;

/* Conversion functions */
FS_INDEX_range_convert_func FS_INDEX_identity_convert;
FS_INDEX_range_convert_func FS_INDEX_S2QD_convert;




/********************************************************************/    
/*                                                                  */
/*                     PROFILE BASED SEARCHES                       */ 
/*                                                                  */
/********************************************************************/    
typedef void FS_INDEX_profile_process_func(FS_SEQ_t FSneighbour); 

typedef int FS_INDEX_profile_range_convert_func(POS_MATRIX *PS0, 
						BIOSEQ *query0,
						int cutoff);	

int FS_INDEX_profile_search(HIT_LIST_t *hit_list0, POS_MATRIX *PS0, 
			    ULINT Pfrom0, ULINT Pto0, int cutoff, 
			    FS_INDEX_profile_process_func *ppfunc0, 
			    FS_INDEX_range_convert_func *cfunc);

/* Process bin functions */
FS_INDEX_profile_process_func FS_INDEX_profile_S_process_bin;
FS_INDEX_profile_process_func FS_INDEX_profile_D_process_bin;

/* Conversion functions */


/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTIONS AND MACROS                    ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    


/********************************************************************/    
/*                                                                  */
/*                     FS_HASH_TABLE module                         */ 
/*                                                                  */
/********************************************************************/    

/* Element access */

ULINT FS_HASH_TABLE_get_total_seqs(FS_HASH_TABLE_t *FSI)
{
  return FSI->no_seqs;
}

/********************************************************************/ 
/*                 Access functions                                 */
/********************************************************************/ 


const char *FS_INDEX_get_db_name(void)
{
  return db_name;
}

FS_PARTITION_t *FS_INDEX_get_ptable(void)
{
  return ptable;
}



const char *FS_index_get_alphabet(void)
{
  return alphabet;
}
#endif /* #if 0 */ 
#endif   

