#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "FSindex.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 0;
int POS_MATRIX_VERBOSE;
FILE *POS_MATRIX_STREAM;

struct rndix
{
  unsigned long i;
  unsigned long rnd;
};

static
int cmp_i(const void *S1, const void *S2)
{
  const struct rndix *T1 = S1;
  const struct rndix *T2 = S2;

  if (T1->i > T2->i)
    return 1;
  if (T1->i < T2->i)
    return -1;
  else
    return 0;
}

static
int cmp_rnd(const void *S1, const void *S2)
{
  const struct rndix *T1 = S1;
  const struct rndix *T2 = S2;

  if (T1->rnd > T2->rnd)
    return 1;
  if (T1->rnd < T2->rnd)
    return -1;
  else
    return 0;
}
      
static 
void print_help(const char *progname)
{
  fprintf(stderr, "Usage: %s database samples\n"
	  ,progname);
}


int main(int argc, char **argv)
{
  const char *db_name = NULL;
  int n;
  ULINT N;
  BIOSEQ *seq;
  int no_args = 2;
  fastadb_arg fastadb_argt[3];
  fastadb_argv_t fastadb_argv[3];
  SEQUENCE_DB *s_db;
  int i;
  struct rndix *p;

  db_name = argv[1];
  n = atoi(argv[2]);
  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      print_help(argv[0]);
      exit(EXIT_FAILURE);
    }

  Try {
    fastadb_argt[0] = ACCESS_TYPE;
    fastadb_argt[1] = RETREIVE_DEFLINES;
    fastadb_argt[2] = NONE;
    fastadb_argv[0].access_type = MEMORY;
    fastadb_argv[1].retrieve_deflines = YES;
    s_db = fastadb_open(db_name, fastadb_argt, fastadb_argv); 
    N = s_db->no_seq;
    n = n > N ? N : n;
    p = mallocec(N*sizeof(struct rndix));
    srand48(time(NULL)); 
    for (i=0; i < N; i++) {
      p[i].i=i;
      p[i].rnd=lrand48();
    }
    qsort(p, N, sizeof(struct rndix), cmp_rnd);
    qsort(p, n, sizeof(struct rndix), cmp_i);
    for (i=0; i < n; i++) {
      seq = s_db->seq + p[i].i;
      bioseq_parse(seq, NULL, YES);
    }
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
