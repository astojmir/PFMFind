#ifndef _FASTADB_H
#define _FASTADB_H

#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef USE_MPATROL
#include <mpatrol.h>
#endif



#define MAX_STORAGE_ROWS 2
#define MAX_SEQUENCES 10000

/* TO DO:

 - fragment functions
 - function descriptions here 
 - void fastadb_convert(SEQUENCE_DB *s_db, char ctable[256])
 - what can we do for non-char conversions ? - void * as start in
                                               BIOSEQ + type field?
 - converting BIOSEQ itself to something is quite easy - maybe we
 should have another object which would contain such object: now we
 miss C++ and inheritance !! 
*/


typedef union 
{
  ULINT id_num;
  char *defline;
} seq_id_t;


typedef struct
  {
    seq_id_t id;
    ULINT len;
    char *start;
  } BIOSEQ; 

typedef enum {MEMORY, RANDOM, SEQUENTIAL} access_type_t;
typedef enum {NO, YES} yn_bool_t; 



typedef struct
{
  /* General */
  const char *db_name;  /* Name of the database */
  FILE *dbfile;         /* Stream to read seqs from */
  int real_file;        /* Non zero if dbfile is not stdin */
  ULINT length;         /* Total Length of Database in residues*/

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


/* ************* bioseq functions ********************* */

int cmp_bioseq(const void *S1, const void *S2);
int bioseq_parse(BIOSEQ *seq, char *filename, yn_bool_t defline);
int bioseq_get_frag(BIOSEQ *seq, BIOSEQ *frag, ULINT start, ULINT
		    length, ULINT id_no);
void string2seq(BIOSEQ *seq, char *string);
void bioseq_random(BIOSEQ *seq, ULINT seq_len, char *alphabet,
		   ULINT a_len);
void bioseq_seq2string(BIOSEQ *seq, char *string, ULINT from, 
		       ULINT to);

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

/* Fragment functions */
int fastadb_init_frags(SEQUENCE_DB *s_db, ULINT min_len, 
		       ULINT max_len);
int fastadb_clear_frags(SEQUENCE_DB *s_db);

int fastadb_get_nofrags(SEQUENCE_DB *s_db, ULINT *no_frags, 
			ULINT min_len, ULINT max_len);
int fastadb_get_frag(SEQUENCE_DB *s_db, BIOSEQ *frag, ULINT frag_no,
		     ULINT min_len, ULINT max_len);

/* Ffragment functions */

void fastadb_init_Ffrags(SEQUENCE_DB *s_db, ULINT len);

ULINT fastadb_count_Ffrags(SEQUENCE_DB *s_db, ULINT len);

int fastadb_get_next_Ffrag(SEQUENCE_DB *s_db, ULINT len, 
			   BIOSEQ *frag, ULINT *offset,
			   int skip);

int fastadb_get_Ffrag(SEQUENCE_DB *s_db, ULINT len, 
		      BIOSEQ *frag, ULINT offset);

int fastadb_find_Ffrag_seq(SEQUENCE_DB *s_db, BIOSEQ *frag, 
			   ULINT *seq_id, ULINT *rel_offset);

int fastadb_get_Ffrag_seq(SEQUENCE_DB *s_db, BIOSEQ *frag, 
			  ULINT seq_id, ULINT from, ULINT to);



#ifdef  MY_INLINE
/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               INLINE FUNCTION DEFINITIONS                    ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

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
	    return 0;
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
#endif /* #ifdef  FASTA_DB_INLINE */



#endif /* #ifndef _FASTADB_H */ 

