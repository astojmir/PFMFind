#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "FSindex.h"
#include "hit_list.h"
#include "misclib.h"

#define MAX_RUNS 500
#define MAX_LINES 10
#define BUF_SIZE 1023

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int main(int argc, char **argv)
{
  BIOSEQ query;
  const char *index_file;
  char *matrix_full;
  char *input_file = NULL;

  time_t time0;
  time_t time1;
  double dt;

  int n = 0;        /* Number of runs */
  int N = MAX_RUNS;
  char *cp;

  char **queries = NULL;
  int *r = NULL;
  int *k = NULL;
  int *bins0 = NULL;
  int *bins1 = NULL;
  int *seqs0 = NULL;
  int *seqs1 = NULL;
  ULINT *resd0 = NULL;
  ULINT *resd1 = NULL;
  
  int M = MAX_LINES;
  int m = 0;
  int *rng; 
  char **out_files;
  char **query_files;
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

  int no_args = 2;

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file matrix_file [input_file]\n", argv[0]);  
      fprintf(stderr,"Input file contains on each line the query file path"
	      "\nand the path of output file"
	      "separated by one (or more) blanks\n");
      exit(EXIT_FAILURE);
    }
  index_file = argv[1];
  matrix_full = argv[2];

  if (argc == no_args + 2) input_file = argv[3];

  Try {
    query.id.defline = "Query sequence";
    /* Load input data:
       Each line of input_file is of the form: query_file sfunc_type pfunc_type output_file
       where query_file contains at each line at least the fields sequence and radius.
       These are meant to come from running the kNN routines.
    */
    if (input_file != NULL) { 
      if((instream = fopen(input_file, "r")) == NULL)
	Throw FSexcept(FOPEN_ERR, "Could not open the file %s.",
		       input_file);
    }
    else instream = stdin; 

    query_files = callocec(M, sizeof(char *));
    out_files = callocec(M, sizeof(char *));
    sfunc_type = callocec(M, sizeof(int));
    pfunc_type = callocec(M, sizeof(int));
	
    while (fgets(buffer, BUF_SIZE, instream)) {
      if (buffer[0] == '#') continue;
      if (m >= M) {
	M += MAX_LINES;
	query_files = reallocec(query_files, M * sizeof(char *));
	out_files = reallocec(out_files, M * sizeof(char *));
	sfunc_type = reallocec(sfunc_type, M * sizeof(int));
	pfunc_type = reallocec(pfunc_type, M * sizeof(int));
      }
      l = strlen(buffer);
      buffer[l-1] = '\0';
      /* Too much but safe */
      query_files[m] = mallocec(l+1);
      out_files[m] = mallocec(l+1);
      sscanf(buffer, "%s %d %d %s", query_files[m], 
	     sfunc_type+m, pfunc_type+m, out_files[m]);
      m++;
    }
    if (input_file != NULL) 
      fclose(instream);

    /* Load matrix */
    S = SCORE_MATRIX_from_file(matrix_full);
    S->set_conv_type(S, QUASI);
    D = S->matrix_conv(S, NULL, 0);

    /* Load index */
    I = FS_INDEX_load(index_file);

    /* Allocate data arrays */
    queries = mallocec(N * sizeof(char *));
    rng = mallocec(N * sizeof(int)); 

    r = mallocec(N * sizeof(int));
    k = mallocec(N * sizeof(int));
    bins0 = mallocec(N * sizeof(int));
    bins1 = mallocec(N * sizeof(int));
    seqs0 = mallocec(N * sizeof(int));
    seqs1 = mallocec(N * sizeof(int));
    resd0 = mallocec(N * sizeof(ULINT));
    resd1 = mallocec(N * sizeof(ULINT));

    for (i=0; i < m; i++) {
      printf("%s\n", query_files[i]);

      /* Load queries and ranges from each query_file */
      n=0;
      if((fp = fopen(query_files[i], "r")) == NULL)
	Throw FSexcept(FOPEN_ERR, "Could not open the file %s.",
		       query_files[i]);

      while (fgets(buffer, BUF_SIZE, fp)) {
	if (buffer[0] == '#') continue;
	if (n >= N) {
	  N += MAX_RUNS;
	  queries = reallocec(queries, N * sizeof(char *));
	  rng = reallocec(rng, N * sizeof(int)); 

	  r = reallocec(r, N * sizeof(int));
	  k = reallocec(k, N * sizeof(int));
	  bins0 = reallocec(bins0, N * sizeof(int));
	  bins1 = reallocec(bins1, N * sizeof(int));
	  seqs0 = reallocec(seqs0, N * sizeof(int));
	  seqs1 = reallocec(seqs1, N * sizeof(int));
	  resd0 = reallocec(resd0, N * sizeof(ULINT));
	  resd1 = reallocec(resd1, N * sizeof(ULINT));
	}
	buffer[strlen(buffer)-1] = '\0'; /* Removes '\n' */

	/* We assume fields are separated by ' ' */
	cp = strchr(buffer, ' ');
	*cp = '\0';
	queries[n] = strdup(buffer);
	sscanf(cp+1, "%d", rng+n);
	n++;
      }
      fclose(fp);

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
	HL = FSINDX_rng_srch(I, &query, D, rng[j], QUASI, HL, sfunc_type[i], pfunc_type[i]);
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

      for (j=0; j < n; j++) 
	free(queries[j]);
    }

    HIT_LIST_destroy(HL);
    FS_INDEX_destroy(I);
    for (i=0; i < m; i++) {
      free(query_files[i]);
      free(out_files[i]);
    }
    free(query_files);
    free(out_files);
    free(sfunc_type);
    free(pfunc_type);
    free(queries);
    free(rng);
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
