#include "misclib.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <unistd.h>
#include "partition.h"

#define BUF_SIZE 256 

int FS_PARTITION_VERBOSE = 0;

static 
void print_help(const char *progname)
{
  fprintf(stderr, "Usage: %s options\n"
	  "\n"
	  "Mandatory:\n"
	  "-p  s   Alphabet partitions\n"
	  "-f  f   Amino acid counts file f\n"
	  "-n  n   Number of points n\n"
	  "-m  n   Fragment length n\n"
	  "-M  n   Maximum size to compute probabilities\n"
	  "Optional:\n"
	  "-o  f   Write results to file f\n"
	  "-b      Don't print progress bar\n"
	  "\n"
	  ,progname);
}

/* Global variables */
static ULINT n; 
static int m;
static int K;
static char *alphabet;
static double Prob[A_SIZE];
static double *Hist;
static double *logfact;
static int max_hist;
static ULINT bins = 0;
static ULINT one_percent_bins;
static FILE *outstream = stdout;
static int nobar_flag = 0;
static double LOG_DBL_MIN;

static
void init_hist(void)
{
  int k;
  double sum = 0.0;
  Hist = callocec(max_hist, sizeof(double));
  logfact = callocec(max_hist, sizeof(double));
  for (k=1; k < max_hist; k++)
    {
      sum += log(k);
      logfact[k] = sum;
    }
}

static
void init_cluster_probs(const char *filename, FS_PARTITION_t *ptable) 
{
  FILE *stream;
  char buffer[BUF_SIZE];
  ULINT i;
  ULINT nn; 
  char c;
  ULINT freq[A_SIZE];
  ULINT cf;
  ULINT total = 0;
  int cluster;
  
  memset(freq, 0, A_SIZE * sizeof(ULINT));
  memset(Prob, 0, A_SIZE * sizeof(double));
  if ((stream = fopen(filename, "r")) == NULL)
    {
      fprintf(stderr, "Unable to open %s\n ", filename);
      exit(EXIT_FAILURE);
    }  

  /* Skip first lines */
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);
  sscanf(buffer+15, "%ld", &nn); 
  fgets(buffer, BUF_SIZE, stream);
  fgets(buffer, BUF_SIZE, stream);

  /* Now get frequencies */
  for (i=0; i < nn; i++)
    {
      fgets(buffer, BUF_SIZE, stream);
      sscanf(buffer+15, "%c", &c);
      cluster = ptable->partition_table[c & A_SIZE_MASK];
      if (cluster != -1)
	{
	  sscanf(buffer+23, "%ld", &cf);
	  freq[cluster] += cf;
	  total += cf;
	}
    }

  /* Divide to get probabilities */
  for (i=0; i < K; i++)
    Prob[i] = ((double) freq[i]) / ((double) total);
  fclose(stream);
}


static inline
void eval_probs(double P)
{
  int k = 0;
  double nP = n*P;
  double lognP = log(nP);
  double logPk;
  double Pk;

  if (!nobar_flag)
    printbar(outstream, bins+1, one_percent_bins, 50);

  for (k=0; k < max_hist; k++)
    {
      logPk = k*lognP - logfact[k] - nP;
      Pk = logPk < LOG_DBL_MIN ? 0.0 : exp(logPk);
      Hist[k] += P*Pk;
    }  
  bins++;
}

static
void compute_probs(int depth, double P)
{
  int j;
  if (depth == m)
    {
      eval_probs(P);
      return;
    }
  for (j=0; j < K; j++)
    {
      compute_probs(depth+1, P*Prob[j]);
    }
  return;
}

static
void print_probs(void)
{
  int i;
  double lambda = (double) n / pow (K,m);
  double Pav;
  double logPav;
  double log_lambda = log(lambda);
  double total = 0.0;
  fprintf(outstream, "*** SIMULATED BIN SIZE PROBABILITIES ***\n");
  fprintf(outstream, "Number of fragments: %ld\n", n);
  fprintf(outstream, "Fragment length: %d\n", m);
  fprintf(outstream, "Number of clusters: %d\n", K);
  fprintf(outstream, "Number of bins: %d\n", (int) pow(K,m));
  fprintf(outstream, "Average bin size: %.2f\n", lambda);
  fprintf(outstream, "Clusters: %s\n\n", alphabet);
  fprintf(outstream, "Cluster probabilities\n");
  for (i=0; i < K; i++)
    fprintf(outstream, "  %d    ", i);
  fprintf(outstream, "\n");
  for (i=0; i < K; i++)
    fprintf(outstream, "%6.4f ", Prob[i]);
  fprintf(outstream, "\n\n");
  fprintf(outstream, "%5s %10s %10s %10s %10s\n", "k", "P(X=k)", 
	  "log P(X=k)", "Pav", "log Pav");
  for (i=0; i < max_hist; i++)
    {
      logPav = i*log_lambda - logfact[i] - lambda;
      if (logPav < LOG_DBL_MIN)
	Pav = 0;
      else
	Pav = exp(logPav);
      fprintf(outstream, "%5d %-10.4e %-10.4e %-10.4e %-10.4e\n", 
	      i, Hist[i], log(Hist[i]), Pav, logPav);
      total += Hist[i];
    }
  fprintf(outstream, "\n");
  fprintf(outstream, "Total probability: %.5f\n", total);
}


int main(int argc, char **argv)
{
  FS_PARTITION_t *ptable;
  char separator = '#';

  /* Getopt variables */
  int c;                              /* Option character          */ 
  extern char *optarg;                /* Option argument           */
  extern int optind;
  extern int optopt;
  extern int opterr;  
  int errflg = 0;                     /* Option error flag         */
  
  const char *count_file = NULL;
  const char *outfile = NULL;


  int partition_flag = 0;
  int freq_flag = 0;
  int n_flag = 0;
  int m_flag = 0;
  int M_flag = 0;
  int outfile_flag = 0;

  /* Input params;
     - partitions
     - amino acid probabilities
     - n (number of database sequences)
     - m (fragment length)
     - max_hist (maximum k)
     - output file
  */

  while ((c = getopt(argc, argv, "p:f:n:m:o:M:b")) != EOF)
    switch (c) 
      {
      case 'p':
	alphabet = optarg;
	partition_flag++;
      case 'f':
	count_file = optarg;
	freq_flag++;
	break;
      case 'n':
	n = atoi(optarg);
	n_flag = 1;
	break;
      case 'm':
	m = atoi(optarg);
	m_flag = 1;
	break;
      case 'M':
	max_hist = atoi(optarg);
	M_flag = 1;
	break;
      case 'o':
	outfile = optarg;
	if((outstream = fopen(outfile, "w")) == NULL)
	  {
	    fprintf(stderr, 
		    "Could not open output file %s!\n", outfile);
	    exit(EXIT_FAILURE);
	  }
	outfile_flag = 1;
	break;
      case 'b':
	nobar_flag = 1;
	break;
      case '?':
      default:
	errflg++;
	break;
      }

  if (!(partition_flag || freq_flag || n_flag || m_flag || M_flag))
    errflg++;

  if (errflg)
    {
      fprintf(stderr,"Error: Insufficient or invalid arguments \n");
      print_help(argv[0]);
      exit(EXIT_FAILURE);
    }

  LOG_DBL_MIN = -2 * log(n);

  /* Initialise partitions */
  ptable = FS_PARTITION_create(alphabet, separator); 
  K = FS_PARTITION_get_no_partitions(ptable);

  /* Prepare progress bar */
  one_percent_bins = pow(K, m)/100; 

  /* Load amino acid frequencies - get alphabet cluster probs */ 
  init_cluster_probs(count_file, ptable);

  /* Initialise histogram, compute factoriels */
  init_hist();

  /* Compute histogram */
  compute_probs(0, 1.0);

  /* Print histogram (normal, log) */
  print_probs();

 return EXIT_SUCCESS;
}
