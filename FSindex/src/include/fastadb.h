#ifndef _FASTADB_H
#define _FASTADB_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bioseq.h"

#include "misclib.h"


#define MAX_STORAGE_ROWS 2
#define MAX_SEQUENCES 10000

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               fastadb                                        ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    


typedef struct
{
  const char *db_name;
  const char *db_path;
  ULINT length; 
  ULINT no_seq; 
  ULINT no_frags;
  ULINT min_len; 
  ULINT max_len; 
} SEQDB_STATS;


typedef enum {MEMORY, RNDM, SEQUENTIAL} access_type_t;

typedef struct
{
  /* General */
  const char *db_name;  /* Name of the database */
  FILE *dbfile;         /* Stream to read seqs from */
  int real_file;        /* Non zero if dbfile is not stdin */
  ULINT length;         /* Total Length of Database in residues*/
  char *linebuf;        /* Getline reading buffer */
  int bufsize;          /* Size of the buffer     */

  /* Sequences */
  ULINT no_seq;         /* Number of sequences in the database */
  ULINT *start_pts;     /* Starts of sequences in terms of residues */
  ULINT max_seqs;       /* Maximum number of sequences that can be put
			   in */
  ULINT current_seq;    /* Current sequence */
  BIOSEQ *seq;          /* Array of sequences                  */
  long *offset;         /* Array of offsets of sequences in 
			   opened file (dbfile) */ 

  /* Fragments */

  ULINT no_frags;       /* Number of fragments         */
  ULINT min_len;        /* Minimum length of fragments */
  ULINT max_len;        /* Maximum length of fragments */
  ULINT *L_start_pts;   /* Start of fragments by length */
  ULINT **f_start_pts;  /* Start of fragments by sequence for each
			   length */  
  ULINT NL;             /* Number of different fragment lengths */
  ULINT current_frag;   /* Current fragment */

  /* Ffragments */
  SLINT Ff_c;            /* Offset of current residue in heap */
  ULINT Ff_l;            /* Length of fragment */
  ULINT Ff_q;            /* Current sequence */
  SLINT Ff_r;            /* Distance of current residue from the end
			    of the current sequence */

  /* Storage Heap */
  char *seq_data;       /* Data heap */
  ULINT seq_data_len;   /* Length of data heap */
  ULINT seq_data_max_len;/* Maximum length of data heap */

  char *deflines;       /* Defline heap */
  ULINT deflines_len;   /* Length of defline heap */
  ULINT deflines_max_len;/* Maximum length of defline heap */ 

  /* Options */
  access_type_t access_type; /* How to access sequences */
  yn_bool_t keep_oldseqs;    /* Keep_previous sequence Y/N */
  yn_bool_t retrieve_deflines; /* Retrieve comment lines Y/N */

} SEQUENCE_DB;

/* These two types enable passing options without accessing
   SEQUENCE_DB directly */

typedef enum {NONE, ACCESS_TYPE, KEEP_OLDSEQS, RETREIVE_DEFLINES,
	      MAX_ROW_LENGTH} fastadb_arg;

typedef union
{
  ULINT max_row_length; 
  access_type_t access_type; 
  yn_bool_t keep_oldseqs;
  yn_bool_t retrieve_deflines;  
} fastadb_argv_t;


/* ************** fastadb functions ******************* */

/* General Purpose Functions */
SEQUENCE_DB *fastadb_open(const char *db_name, fastadb_arg *argt, 
			  fastadb_argv_t *argv);
int fastadb_close(SEQUENCE_DB *s_db);
int fastadb_get_length(SEQUENCE_DB *s_db, ULINT *length);

/* Sequence functions */
int fastadb_get_noseqs(SEQUENCE_DB *s_db, ULINT *no_seqs);
int fastadb_get_seq(SEQUENCE_DB *s_db, ULINT seq_no, BIOSEQ **seq);
int fastadb_get_next_seq(SEQUENCE_DB *s_db, BIOSEQ **seq);

/* Offset <--> char pter conversion functions (as macros) */
char *fastadb_data_pter(SEQUENCE_DB *s_db, ULINT offset);
ULINT fastadb_data_offset(SEQUENCE_DB *s_db, char *s);

/* Fragment functions */
int fastadb_init_frags(SEQUENCE_DB *s_db, ULINT min_len, 
		       ULINT max_len);
int fastadb_clear_frags(SEQUENCE_DB *s_db);

int fastadb_get_nofrags(SEQUENCE_DB *s_db, ULINT *no_frags, 
			ULINT min_len, ULINT max_len);
int fastadb_get_frag(SEQUENCE_DB *s_db, BIOSEQ *frag, ULINT frag_no,
		     ULINT min_len, ULINT max_len);

/* Ffragment functions */
MY_INLINE
void fastadb_init_Ffrags(SEQUENCE_DB *s_db, ULINT len);
MY_INLINE
ULINT fastadb_count_Ffrags(SEQUENCE_DB *s_db, ULINT len);
MY_INLINE
int fastadb_get_next_Ffrag(SEQUENCE_DB *s_db, ULINT len, 
			   BIOSEQ *frag, ULINT *offset,
			   int skip);
MY_INLINE
int fastadb_get_Ffrag(SEQUENCE_DB *s_db, ULINT len, 
		      BIOSEQ *frag, ULINT offset);
MY_INLINE
int fastadb_find_Ffrag_seq(SEQUENCE_DB *s_db, BIOSEQ *frag, 
			   ULINT *seq_id, ULINT *rel_offset);
MY_INLINE
int fastadb_get_Ffrag_seq(SEQUENCE_DB *s_db, BIOSEQ *frag, 
			  ULINT seq_id, ULINT from, ULINT to);

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               MIXED_ITER                                     ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

/* Fragment iterator for FS_INDEX */

typedef struct
{
  int l;
  int m;
  int skip;
  BIOSEQ *first_seq;
  BIOSEQ *last_seq;
  BIOSEQ *seq;
  char *c;
  int c_len;
} MIXED_ITER;


MY_INLINE
MIXED_ITER *MIXED_ITER_init(SEQUENCE_DB *s_db, int min_len,
			    int len, int skip);
MY_INLINE
void MIXED_ITER_del(MIXED_ITER *iterator);
MY_INLINE
char *MIXED_ITER_next(MIXED_ITER *iterator, ULINT *len);
MY_INLINE
void MIXED_ITER_reset(MIXED_ITER *iterator);

#ifdef UFRAGDB
/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               ufragdb                                        ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

typedef struct 
{
  char *db_name;         /* Fasta database name                     */
  int db_name_len;       /* Fasta db name length                    */
  SEQUENCE_DB *sdb;      /* Fasta database                          */
  char *abet;            /* Valid alphabet                          */  
  int abet_len;          /* Alphabet string length                  */
  ULINT frag_len;        /* Fragment length                         */
  int skip;              /* How many fragments to skip (1=none)     */
  FILE *fptr;            /* File handle                             */

  ULINT *pts;            /* Offsets of unique sequences             */
  ULINT pts_size;        /* Full Number of 'unique' points          */
  ULINT upts_size;       /* Number of points without duplicates     */
  ULINT dpts_size;       /* Number of points with duplicates        */
  ULINT *dup_offset;     /* Offsets of duplicates in heap           */
  ULINT *dup_heap;       /* Heap containing all the duplicates      */ 
  ULINT dup_heap_size;			    
  ULINT cur_pt;          /* Current point                           */

  ULINT pts0;            /* Buffered read - starting index          */
  ULINT pts1;            /* Buffered read - one past end index      */
} UFRAG_DB;


UFRAG_DB *ufragdb_init(void);

void ufragdb_del(UFRAG_DB *udb);

void ufragdb_create(UFRAG_DB *udb, const char *db_name, 
		    ULINT frag_len, const char *abet, int skip);

void ufragdb_load(UFRAG_DB *udb, const char *udb_name,
		  int options);

void ufragdb_save(UFRAG_DB *udb, const char *udb_name);

ULINT ufragdb_get_nopts(UFRAG_DB *udb);

int ufragdb_get_ufrag(UFRAG_DB *udb, ULINT i);

void ufragdb_init_ufrags(UFRAG_DB *udb, ULINT i0);

int ufragdb_get_next_ufrag(UFRAG_DB *udb);

int ufragdb_get_dfrags(UFRAG_DB *udb, ULINT frag_offset, 
		       ULINT **dfrags, ULINT *dsize);

void ufragdb_print(UFRAG_DB *udb, FILE *stream);

int cmp_ufrags(const void *S1, const void *S2);
#endif /* #ifdef UFRAGDB */

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS                    ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

#define fastadb_data_pter(db, i) ((db)->seq_data + (i))
#define fastadb_data_offset(db, s) ((s) - (db)->seq_data)
#define fastadb_end_heap(db) ((db)->seq[(db)->no_seq-1].start        \
        + (db)->seq[(db)->no_seq-1].len)


/********************************************************************/    
/*                                                                  */
/*                     Ffragment functions                          */ 
/*                                                                  */
/********************************************************************/    

MY_INLINE
void fastadb_init_Ffrags(SEQUENCE_DB *s_db, ULINT len)
{
  s_db->Ff_q = 0;
  s_db->Ff_c = 0;
  s_db->Ff_r = s_db->seq[0].len - len + 1;
}


MY_INLINE
ULINT fastadb_count_Ffrags(SEQUENCE_DB *s_db, ULINT len)
{
  ULINT i;
  ULINT S = 0;
  for (i=0; i < s_db->no_seq; i++)
    if (s_db->seq[i].len >= len) 
      {
	S += s_db->seq[i].len - len + 1;
      }
  return S;
}

MY_INLINE
int fastadb_get_next_Ffrag(SEQUENCE_DB *s_db, ULINT len, 
			   BIOSEQ *frag, ULINT *offset,
			   int skip)
{
  if (s_db->Ff_c >= s_db->Ff_r)
    {
      do 
	{
	  s_db->Ff_q++;
	  if (s_db->Ff_q == s_db->no_seq)
	    {
#if 0
	      Throw FSexcept(INDEX_OUT_OF_RANGE, 
			     "fastadb_get_next_Ffrag():"
			     " No more fragments. "
			     "Ff_q=%ld, no_seq=%ld",
			     s_db->Ff_q, s_db->no_seq);
#endif
	      return 0;
	    }
	}  
      while (s_db->seq[s_db->Ff_q].len < len);
      s_db->Ff_c = s_db->seq[s_db->Ff_q].start - s_db->seq_data;
      s_db->Ff_r = s_db->Ff_c + s_db->seq[s_db->Ff_q].len - len + 1;
    }

  frag->start = s_db->seq_data + s_db->Ff_c;
  frag->len = len;
  *offset = s_db->Ff_c;
  s_db->Ff_c+= skip;
  return 1;
}

MY_INLINE
int fastadb_get_Ffrag(SEQUENCE_DB *s_db, ULINT len, 
		      BIOSEQ *frag, ULINT offset)
{
  frag->start = s_db->seq_data + offset;
  frag->len = len;
  /* Should really check length ! */
  return 1;
 }

MY_INLINE
int fastadb_find_Ffrag_seq(SEQUENCE_DB *s_db, BIOSEQ *frag, 
			   ULINT *seq_id, ULINT *rel_offset) 
{
  ULINT a = 0;
  ULINT b = s_db->no_seq - 1;
  ULINT k = b/2; 

  if (frag->start > s_db->seq[s_db->no_seq - 1].start)
    {
      *seq_id = s_db->no_seq - 1;
      *rel_offset = 
	frag->start - s_db->seq[s_db->no_seq - 1].start;
      return 1;
    }

  while (a<=b)
    {
      k=(a+b)/2;
      if (s_db->seq[k].start <= frag->start)
        {
          if (s_db->seq[k+1].start > frag->start)
	    {
	      *seq_id = k;
	      *rel_offset = 
		frag->start - s_db->seq[k].start;
	      return 1;
	    }
	  else
	    a=k+1;
	}
      else
	b=k;
    }
  *seq_id = k;
  *rel_offset = 
    frag->start - s_db->seq[k].start;
  return 1;
}

MY_INLINE
int fastadb_get_Ffrag_seq(SEQUENCE_DB *s_db, BIOSEQ *frag, 
			  ULINT seq_id, ULINT from, ULINT to)
{
  frag->len = to - from + 1;
  frag->start = s_db->seq[seq_id].start + from;
  return 1;
}

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               MIXED_ITER functions                           ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

MY_INLINE
MIXED_ITER *MIXED_ITER_init(SEQUENCE_DB *s_db, int min_len,
			    int len, int skip)
{
  MIXED_ITER *I = mallocec(sizeof(MIXED_ITER));
  I->l = len;
  I->m = min_len > len ? len : min_len;
  I->skip = skip;
  I->first_seq = s_db->seq;
  I->last_seq = s_db->seq + s_db->no_seq - 1;
  I->seq = I->first_seq;
  I->c = I->first_seq->start;
  I->c_len = I->first_seq->len;
  return I;
}

MY_INLINE
void MIXED_ITER_del(MIXED_ITER *I)
{
  free(I);
}

MY_INLINE
char *MIXED_ITER_next(MIXED_ITER *I, ULINT *len)
{
  char *s;
  if (I->c_len < I->m) {
    do {
      I->seq++;
      if (I->seq > I->last_seq)
	return NULL;
    }  
    while (I->seq->len < I->m);
    I->c = I->seq->start;
    I->c_len = I->seq->len;
  }
  s = I->c;
  *len = I->c_len > I->l ? I->l : I->c_len;
  I->c += I->skip;
  I->c_len -= I->skip;
  return s;
}

MY_INLINE
void MIXED_ITER_reset(MIXED_ITER *I)
{
  I->seq = I->first_seq;
  I->c = I->first_seq->start;
  I->c_len = I->first_seq->len;
}

#endif /* #ifndef _FASTADB_H */ 

