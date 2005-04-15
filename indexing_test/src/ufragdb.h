#ifndef _UFRAGDB_H
#define _UFRAGDB_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "bioseq.h"
#include "misclib.h"
#include "fastadb.h"

typedef struct 
{
  char *db_name;         /* Fasta database name                     */
  int db_name_len;       /* Fasta db name length                    */
  SEQUENCE_DB *sdb;      /* Fasta database                          */
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
		    ULINT frag_len, int skip);

void ufragdb_load(UFRAG_DB *udb, const char *udb_name,
		  int options);

void ufragdb_save(UFRAG_DB *udb, const char *udb_name);

void ufragdb_read(UFRAG_DB *udb, FILE *fp, char *db_dir);

void ufragdb_write(UFRAG_DB *udb, FILE *fp);



ULINT ufragdb_get_nopts(UFRAG_DB *udb);

int ufragdb_get_ufrag(UFRAG_DB *udb, ULINT i);

void ufragdb_init_ufrags(UFRAG_DB *udb, ULINT i0);

int ufragdb_get_next_ufrag(UFRAG_DB *udb);

int ufragdb_get_dfrags(UFRAG_DB *udb, ULINT frag_offset, 
		       ULINT **dfrags, ULINT *dsize);

int ufragdb_get_dfrags2(UFRAG_DB *udb, ULINT i, 
			ULINT **dfrags, ULINT *dsize);

void ufragdb_print(UFRAG_DB *udb, FILE *stream);

int cmp_ufrags(const void *S1, const void *S2);




#endif /* _UFRAGDB_H */
