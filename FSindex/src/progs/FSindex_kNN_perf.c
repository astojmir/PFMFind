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


#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "FSindex.h"
#include "hit_list.h"
#include "misclib.h"

#define MAX_RUNS 500
#define MAX_kNN 10
#define BUF_SIZE 1023

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int main(int argc, char **argv)
{
  BIOSEQ query;
  const char *index_file;
  const char *query_file;
  char *matrix_full;
  char *input_file = NULL;

  time_t time0;
  time_t time1;
  double dt;

  int n = 0;        /* Number of runs */
  int N = MAX_RUNS;

  char **queries;
  int *r;
  int *k;
  int *bins0;
  int *bins1;
  int *seqs0;
  int *seqs1;
  ULINT *resd0;
  ULINT *resd1;
  
  int M = MAX_kNN;
  int m = 0;
  int *kNN; 
  char **out_files;
  int l;
  int *sfunc_type;
  int *pfunc_type;

  int i, j;

  FILE *fp;
  FILE *outstream = stdout;
  FILE *instream = stdin;

  SCORE_MATRIX *S;
  SCORE_MATRIX *D;
  HIT_LIST_t *HL = NULL;
  FSINDX *I;
  char buffer[BUF_SIZE+1];

  int no_args = 3;

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file query_file matrix_file [input_file]\n", argv[0]);  
      fprintf(stderr,"Input file contains on each line the number"
	      " of nearest neighbours to search\nand the path of output file"
	      "separated by one (or more) blanks\n");
      exit(EXIT_FAILURE);
    }
  index_file = argv[1];
  query_file = argv[2];
  matrix_full = argv[3];

  if (argc == 4 + 1) input_file = argv[4];

  Try {
    /* Load kNN */
    if (input_file != NULL) { 
      if((instream = fopen(input_file, "r")) == NULL)
	Throw FSexcept(FOPEN_ERR, "Could not open the file %s.",
		       input_file);
    }
    else instream = stdin; 

    kNN = callocec(M, sizeof(int)); 
    out_files = callocec(M, sizeof(char *));
    sfunc_type = callocec(M, sizeof(int));
    pfunc_type = callocec(M, sizeof(int));
	
    while (fgets(buffer, BUF_SIZE, instream)) {
      if (buffer[0] == '#') continue;
      if (m >= M) {
	M += MAX_kNN;
	kNN = reallocec(kNN, M * sizeof(int));
	out_files = reallocec(out_files, M * sizeof(char *));
	sfunc_type = reallocec(sfunc_type, M * sizeof(int));
	pfunc_type = reallocec(pfunc_type, M * sizeof(int));
      }
      l = strlen(buffer);
      buffer[l-1] = '\0';
      out_files[m] = mallocec(l+1);
      sscanf(buffer, "%d %d %d %s", kNN+m, 
	     sfunc_type+m, pfunc_type+m, 
	     out_files[m]);
      m++;
    }
    if (input_file != NULL) 
      fclose(instream);


    /* Load queries */
    if((fp = fopen(query_file, "r")) == NULL)
      Throw FSexcept(FOPEN_ERR, "Could not open the file %s.",
		     query_file);

    queries = mallocec(N * sizeof(char *));
    while (fgets(buffer, BUF_SIZE, fp)) {
      if (buffer[0] == '#') continue;
      if (n >= N) {
	N += MAX_RUNS;
	queries = reallocec(queries, N * sizeof(char *));
      }
      buffer[strlen(buffer)-1] = '\0';
      queries[n++] = strdup(buffer);
    }
    fclose(fp);

    query.id.defline = "Query sequence";

    r = mallocec(n * sizeof(int));
    k = mallocec(n * sizeof(int));
    bins0 = mallocec(n * sizeof(int));
    bins1 = mallocec(n * sizeof(int));
    seqs0 = mallocec(n * sizeof(int));
    seqs1 = mallocec(n * sizeof(int));
    resd0 = mallocec(n * sizeof(ULINT));
    resd1 = mallocec(n * sizeof(ULINT));

    /* Load matrix */
    S = SCORE_MATRIX_from_file(matrix_full);
    S->set_conv_type(S, QUASI);
    D = S->matrix_conv(S, NULL, 0);

    /* Load index */
    I = FS_INDEX_load(index_file);

    for (i=0; i < m; i++) {
      printf("NN=%d\n", kNN[i]);
      memset(r, 0, n * sizeof(int));
      memset(k, 0, n * sizeof(int));
      memset(bins0, 0, n * sizeof(int));
      memset(bins1, 0, n * sizeof(int));
      memset(seqs0, 0, n * sizeof(int));
      memset(seqs1, 0, n * sizeof(int));
      memset(resd0, 0, n * sizeof(ULINT));
      memset(resd1, 0, n * sizeof(ULINT));

      time0 = time(NULL);
      for (j=0; j < n; j++) {
	query.start = queries[j];
	query.len = strlen(queries[j]);
	HL = FSINDX_kNN_srch(I, &query, D, kNN[i], HL, sfunc_type[i], pfunc_type[i]);
	HIT_LIST_sort_incr(HL);
	r[j] = HL->hits[HL->actual_seqs_hits-1].dist;
	k[j] = HL->actual_seqs_hits;
	bins0[j] = HL->FS_seqs_visited;
	bins1[j] = HL->FS_seqs_hits;
	seqs0[j] = HL->seqs_visited;
	seqs1[j] = HL->seqs_hits;
	resd0[j] = HL->useqs_visited;
	resd1[j] = HL->useqs_hits;
      }
      time1 = time(NULL);
      dt = difftime(time1, time0);

      if (out_files[i] != NULL) { 
	if((outstream = fopen(out_files[i], "w")) == NULL)
	  Throw FSexcept(FOPEN_ERR, "Could not open the file %s.",
			 out_files[i]);
      }
      else outstream = stdout; 

      fprintf(outstream, "# time = %.2f usecs\n", dt);
      for (j=0; j < n; j++) {
	fprintf(outstream, "%s %5d %5d %9d %9d %9d %9d %12ld %12ld\n",
		queries[j], r[j], k[j], bins0[j], bins1[j],
		seqs0[j], seqs1[j], resd0[j], resd1[j]);
      }

      if (out_files[i] != NULL)
	fclose(outstream);
    }

    HIT_LIST_destroy(HL);
    FS_INDEX_destroy(I);
    for (j=0; j < m; j++)
      free(out_files[j]);
    free(out_files);
    free(kNN);
    free(sfunc_type);
    free(pfunc_type);
    for (j=0; j < n; j++) 
      free(queries[j]);
    free(queries);
    free(r);
    free(k);
    free(bins0);
    free(bins1);
    free(seqs0);
    free(seqs1);
    free(resd0);
    free(resd1);
    SCORE_MATRIX_del(S);
    SCORE_MATRIX_del(D);
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
