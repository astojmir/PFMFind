#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "misclib.h"
#include "blastlib.h"
#include "fastadb.h"


/*#define STRDUP*/

# define BUFFER_ALLOC 640
/* ************* smatrix module ****************** */

SMATRIX smatrix_load(const char *filename)
{
  /* Here we leave security totally behind, to pick it up later if
     necessary. We assume:
     - 23 amino acid alphabet (standard + B Z X - we ignore *)
     - all letters are capitals so the first usable index is 65 
     
     This is known to be true for SwissProt and BLOSUM62.
     
  */


  char linebuf[512]; /* line buffer - hopefully 512 bytes should be
		       enough for now */
  int alphabet[32];
  int header_read = 0;
  ULINT i, j;
  SMATRIX S = callocec(32, sizeof(int *));
  FILE *stream = fopen(filename, "r");
  if (stream == NULL)
    return NULL;
  for (i = 0; i < 32; i++)
    {
      S[i] = callocec(32, sizeof(int));
      S[i] -= 65;
    }
  S -= 65;

  j = 0;
  while (fgets(linebuf,512,stream) != NULL)
    {
      if (linebuf[0] == '#') continue;
      if (linebuf[0] == ' ' && !header_read)
	{
	  for (i=0; i < 23; i++)
	    alphabet[i] = linebuf[3+3*i]; 
	  header_read = 1;
	}
      else
	{
	  for (i=0; i < 23; i++)
	    sscanf(linebuf+2+3*i, "%d", S[alphabet[j]]+alphabet[i]);
	  j++;
	  if (j >= 23) break;
	}
    }
  fclose(stream);
  return S;
}

void smatrix_remove(SMATRIX S)
{
  int i;
  for (i = 0; i < 32; i++)
    free(S[i+65]+65);
  free(S+65);
}

SLINT smatrix_max(SMATRIX S)
{
  ULINT i,j;
  SLINT max_S = 0;
  for (i=0; i < 32; i++)
    for (j=i; j < 32; j++)
      if (S[i+65][j+65] > max_S) 
	max_S = S[i+65][j+65];
  return max_S;

}

SMATRIX smatrix2dist(SMATRIX S)
{
  ULINT i, j;
  SMATRIX D = callocec(32, sizeof(int *));
  for (i = 0; i < 32; i++)
    {
      D[i] = callocec(32, sizeof(int));
      D[i] -= 65;
    }
  D -= 65;
  for (i = 65; i < 97; i++)
    for (j = 65; j < 97; j++)
      D[i][j] = S[i][i] + S[j][j] - 2*S[i][j];
  return D;
}

inline
int Hamming_dist(BIOSEQ *seq1, BIOSEQ *seq2, SMATRIX D)
{
  /* Make sure you can compare two sequences - select the minimum
     length */
  ULINT len = seq1->len > seq2->len ? seq2->len : seq1->len;
  ULINT i;
  int d = 0;
  for (i = 0; i < len; i++)
    d += D[seq1->start[i]][seq2->start[i]];
  return d;
}



/* ************* End of smatrix module *********** */




/* ************* seedhist module *************** */


SEED_HIST *seedhist_init(SLINT min_bin, SLINT max_bin)
{
  SEED_HIST *hist;
  if (min_bin > max_bin) return NULL;
  hist = callocec(1, sizeof(SEED_HIST));
  hist->min_bin = min_bin;
  hist->no_bins = 1 + max_bin - min_bin;
  hist->freq = callocec(hist->no_bins, sizeof(ULINT));
  return hist;
}

void seedhist_init1(SEED_HIST *hist, SLINT min_bin, SLINT max_bin)
{
  if (min_bin > max_bin) 
    {
#if DEBUG > 0
      printf("Invalid seedhist_init1 range.\n");
#endif
      return;
    }
  hist->min_bin = min_bin;
  hist->no_bins = 1 + max_bin - min_bin;
  hist->freq = (ULINT *) callocec(hist->no_bins, sizeof(ULINT));
}


int seedhist_clear(SEED_HIST *hist)
{
  free(hist->freq);
  free(hist);
  return 1;
}

int seedhist_clear1(SEED_HIST *hist)
{
  free(hist->freq);
  return 1;
}

inline 
void seedhist_add_acount(SEED_HIST *hist, SLINT bin)
{
  SLINT old_no_bins;
  SLINT dsize;
  ULINT *old_freq;

  if (bin < hist->min_bin)
    {
#if DEBUG > 0
      printf("Smaller bin!\n");
#endif
      old_no_bins = hist->no_bins;
      dsize = hist->min_bin - bin;
      hist->min_bin = bin;
      hist->no_bins += dsize;
      old_freq = hist->freq;
      hist->freq = callocec(hist->no_bins, sizeof(ULINT));
      memcpy(hist->freq+dsize, old_freq, old_no_bins*sizeof(ULINT)); 
      free(old_freq);
    }
  else if (bin > (hist->no_bins+hist->min_bin-1))
    {
#if DEBUG > 0
      printf("Larger bin!\n");
#endif
      dsize = bin - hist->no_bins - hist->min_bin +1;
      old_no_bins = hist->no_bins;
      old_freq = hist->freq;
      hist->no_bins += dsize;
      hist->freq = reallocec(hist->freq, hist->no_bins*sizeof(ULINT));
      memset(hist->freq+old_no_bins, 0, dsize*sizeof(ULINT));
    }
  hist->freq[bin - hist->min_bin]++;
  hist->attainable_pts++;
  return ;
}


int seedhist_merge(SEED_HIST *hist1, SEED_HIST *hist2)
{
  ULINT i;
  if ((hist1->no_bins != hist2->no_bins) || 
      (hist1->min_bin != hist2->min_bin))
    return 0;

  hist1->attainable_pts += hist2->attainable_pts;
  hist1->unattainable_pts += hist2->unattainable_pts;
  for (i=0; i < hist1->no_bins; i++)
    hist1->freq[i] += hist2->freq[i];
  return 1;
}

int seedhist_read(SEED_HIST *hist, FILE *stream)
{
  fread(&hist->no_bins, sizeof(SLINT), 1, stream);
  fread(&hist->min_bin, sizeof(SLINT), 1, stream);
  fread(&hist->attainable_pts, sizeof(ULINT), 1, stream);
  fread(&hist->unattainable_pts, sizeof(ULINT), 1, stream);
  hist->freq = callocec(hist->no_bins, sizeof(ULINT));
  fread(hist->freq, sizeof(ULINT), hist->no_bins, stream);
  return 1;
}

int seedhist_write(SEED_HIST *hist, FILE *stream)
{
  fwrite(&hist->no_bins, sizeof(SLINT), 1, stream);
  fwrite(&hist->min_bin, sizeof(SLINT), 1, stream);
  fwrite(&hist->attainable_pts, sizeof(ULINT), 1, stream);
  fwrite(&hist->unattainable_pts, sizeof(ULINT), 1, stream);
  fwrite(hist->freq, sizeof(ULINT), hist->no_bins, stream);
  return 1;
}

int seedhist_print(SEED_HIST *hist, FILE *stream)
{
  ULINT i;
  ULINT min0 = 0;
  ULINT max0 = hist->no_bins -1;
  ULINT counted_seq = 0;

  while(hist->freq[min0] == 0 && min0 < hist->no_bins)
    min0++;
  while(hist->freq[max0] == 0 && max0 > 0)
    max0--;
#if DEBUG > 5
  printf("min0 = %ld\n", min0);
  printf("max0 = %ld\n", max0);
  fflush(stdout);
#endif
  for (i=min0; i <= max0; i++)
    counted_seq += hist->freq[i];

  fprintf(stream, "Points: %ld %ld\n", hist->attainable_pts,
	  hist->unattainable_pts);
#if DEBUG > 5
  fprintf(stream, "Counted points: %ld\n", counted_seq);
#endif
  fprintf(stream, "Values: %ld\n", max0+1-min0);
  fprintf(stream, "%12s %12s\n", "Value", "Frequency");

  for (i=min0; i <= max0; i++)
    fprintf(stream, "%12ld %12ld\n", hist->min_bin+i, hist->freq[i]);
  fprintf(stream, "\n");
  return 1;
}




/* ************* End of seedhist module *********** */



/* ************* blast_params module ****************** */

BPARAMS *default_blast_params(void)
{
  BPARAMS *BP = mallocec(sizeof(BPARAMS));

  BP->program_name = strdup("blastp");
  BP->expectation_value = 10.0;
  BP->filter_query_sequence = strdup("T");
  BP->cost_to_open_gap = 11;
  BP->cost_to_extend_gap = 1;
  BP->X_dropoff_value = 0;
  BP->nt_mismatch_penalty = -3;
  BP->nt_match_reward = 1;
  BP->no_one_line_descriptions = 500;
  BP->no_allignments = 50000;
  BP->threshold_for_extending_hits = 0;
  BP->gapped_alignment = strdup("T");
  BP->query_genetic_code = 1;
  BP->db_genetic_code = 1;
  BP->believe_defline = strdup("F");
  BP->matrix = strdup("BLOSUM62");
  BP->word_size = 0;
  BP->effective_db_length = 0;
  BP->no_best_hits_to_keep = 0;
  BP->no_passes = 0;
  BP->effective_sspace_length = 0;
  BP->query_strands_to_search = 3;
  BP->dropoff_for_extensions = 0.0;
  BP->X_dropoff_for_final_gapped_alignment = 0;

  return BP;
}

void read_blast_params(BPARAMS *BP, 
		       const char *BPfilename)
{
  char linebuf[1024];
  char dumpstr[40];
  FILE *BPfile = fopen(BPfilename, "r");
  if (BPfile == NULL)
    {
      fprintf(stderr, "Unable to open %s\n ", BPfilename);
      exit(EXIT_FAILURE);
    }

  /* Allocate character strings to maximum lengths */
  BP->program_name = mallocec(20);  
  BP->filter_query_sequence = mallocec(2); /* "T"/"F" */
  BP->gapped_alignment = mallocec(2);
  BP->believe_defline = mallocec(2);
  BP->matrix = mallocec(64);

  fgets(linebuf,1024, BPfile);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %19s", dumpstr, 
	 BP->program_name);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %f", dumpstr, 
	 &BP->expectation_value);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %1s", dumpstr, 
	 BP->filter_query_sequence);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->cost_to_open_gap);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->cost_to_extend_gap);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->X_dropoff_value);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->nt_mismatch_penalty);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->nt_match_reward);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->no_one_line_descriptions);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->no_allignments);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->threshold_for_extending_hits);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %1s", dumpstr, 
	 BP->gapped_alignment);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->query_genetic_code);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->db_genetic_code);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %1s", dumpstr, 
	 BP->believe_defline);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %63s", dumpstr, 
	 BP->matrix);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->word_size);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %f", dumpstr, 
	 &BP->effective_db_length);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->no_best_hits_to_keep);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->no_passes);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %f", dumpstr, 
	 &BP->effective_sspace_length);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->query_strands_to_search);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %f", dumpstr, 
	 &BP->dropoff_for_extensions);
  fgets(linebuf,1024, BPfile);
  sscanf(linebuf, "%36s %d", dumpstr, 
	 &BP->X_dropoff_for_final_gapped_alignment);
  
  fclose(BPfile);
}

void write_blast_params(BPARAMS *BP, 
		       const char *BPfilename)
{
  FILE *BPfile = fopen(BPfilename, "w");
  if (BPfile == NULL)
    {
      fprintf(stderr, "Unable to open %s\n ", BPfilename);
      exit(EXIT_FAILURE);
    }
  fprintf(BPfile, "***** Blast Options *****\n");
  fprintf(BPfile, "%-36s %s\n", "Program Name:", 
	  BP->program_name);
  fprintf(BPfile, "%-36s %.4f\n", "Expectation Value:",
	  BP->expectation_value);
  fprintf(BPfile, "%-36s %s\n", "Filter query sequence:",
	  BP->filter_query_sequence);
  fprintf(BPfile, "%-36s %d\n", "Cost to open gap:",
	  BP->cost_to_open_gap);
  fprintf(BPfile, "%-36s %d\n", "Cost to extend gap:",
	  BP->cost_to_extend_gap);
  fprintf(BPfile, "%-36s %d\n", "X dropoff value:", 
	  BP->X_dropoff_value);
  fprintf(BPfile, "%-36s %d\n", "Nucleotide mismatch penalty:", 
	  BP->nt_mismatch_penalty);
  fprintf(BPfile, "%-36s %d\n", "Nucleotide match reward:", 
	  BP->nt_match_reward);
  fprintf(BPfile, "%-36s %d\n", "Number of one line descriptions:", 
	  BP->no_one_line_descriptions);
  fprintf(BPfile, "%-36s %d\n", "Number of alignments:", 
	  BP->no_allignments);
  fprintf(BPfile, "%-36s %d\n", "Threshold for extending hits:", 
	  BP->threshold_for_extending_hits);
  fprintf(BPfile, "%-36s %s\n", "Gapped alignement:", 
	  BP->gapped_alignment);
  fprintf(BPfile, "%-36s %d\n", "Query genetic code:", 
	  BP->query_genetic_code);
  fprintf(BPfile, "%-36s %d\n", "Database genetic code:", 
	  BP->db_genetic_code);
  fprintf(BPfile, "%-36s %s\n", "Believe query defline:", 
	  BP->believe_defline);
  fprintf(BPfile, "%-36s %s\n", "Scoring matrix:", 
	  BP->matrix);
  fprintf(BPfile, "%-36s %d\n", "Word size:", 
	  BP->word_size);
  fprintf(BPfile, "%-36s %.4f\n", "Effective db length:", 
	  BP->effective_db_length);
  fprintf(BPfile, "%-36s %d\n", "Number of best hits to keep:", 
	  BP->no_best_hits_to_keep);
  fprintf(BPfile, "%-36s %d\n", "Number of passes:", 
	  BP->no_passes);
  fprintf(BPfile, "%-36s %.4f\n", "Effective search space length:", 
	  BP->effective_sspace_length);
  fprintf(BPfile, "%-36s %d\n", "Query strands to search:", 
	  BP->query_strands_to_search);
  fprintf(BPfile, "%-36s %.4f\n", "Dropoff for extensions:", 
	  BP->dropoff_for_extensions);
  fprintf(BPfile, "%-36s %d\n", "Dropoff for final gapped"
	  " alignment:", 
	  BP->X_dropoff_for_final_gapped_alignment);
  
  fclose(BPfile);
    
}

void clean_blast_params(BPARAMS *BP)
{
  free(BP);
}

/* ************* End of blast_params module *********** */




/* ************* blastlib_params module *************** */


BLPARAMS *blastlib_params_init(char *fasta_db, ULINT max_seqs,
			       float max_time, USINT paralel_instance,
			       USINT paralel_total, 
			       USINT progress_report_freq,
			       blast_xml_output_func *blast_xml_parse, 
			       void *parser_params)
{
  BLPARAMS *BLP = mallocec(sizeof(BLPARAMS));
  memset(BLP, 0, sizeof(BLPARAMS));

  BLP->fasta_db = fasta_db;
  BLP->max_seqs = max_seqs;
  BLP->max_time = max_time;
  BLP->paralel_instance = paralel_instance;
  BLP->paralel_total = paralel_total;
  BLP->progress_report_freq = progress_report_freq;
  BLP->blast_xml_parse = blast_xml_parse;
  BLP->parser_params = parser_params;  
  return BLP;
}

int blastlib_params_clear(BLPARAMS *BLP)
{
  free(BLP);
  return 1;
}

static
void print_progress(char *prog_tmpfile, BLPARAMS *BLP, int finished)
{
  FILE *progress_file;

  progress_file = fopen(prog_tmpfile, "w");
  if (progress_file == NULL)
    return;

  fprintf(progress_file, "Progress Report\n");
  fprintf(progress_file, "Database: %s\n", BLP->fasta_db);

  if (BLP->paralel_total > 0)
    {
      fprintf(progress_file, "Paralel processing\n");
      fprintf(progress_file, "Instance: %d out of %d \n",
	      BLP->paralel_instance, BLP->paralel_total);
      fprintf(progress_file, "Sequences: from %ld to %ld\n",
	      BLP->first_seq, BLP->last_seq);
    }
  else
    {
      fprintf(progress_file, "Random sampling\n");
      fprintf(progress_file, "Maximum running time: %.2f mins\n",
	      BLP->max_time);
      fprintf(progress_file, "Maximum sequences to sample: %ld\n",
	      BLP->max_seqs);
    }
  fprintf(progress_file, "Running time: %.2f mins\n",
	  BLP->current_time);
  fprintf(progress_file, "Last sequence no. processed: %ld\n",
	  BLP->current_seq);
  fprintf(progress_file, "Total number of BLAST runs so far: %ld\n",
	  BLP->current_run+1);
  fprintf(progress_file, "Average rate: %.2f runs per minute\n",
	  (float) (BLP->current_run+1)/BLP->current_time);

  if (finished)
      fprintf(progress_file, "***** FINISHED *****\n");

  fclose(progress_file);
}

void blastlib_run_blast_expt(BLPARAMS *BLP, BPARAMS *BP, 
			     SEQUENCE_DB* s_db)
{
  char seq_tmpfile[64]; 
  char xml_tmpfile[64];
  char prog_tmpfile[64];
  pid_t current_pid = getpid();
  time_t t0 = time(NULL);
  time_t t1;
  ULINT no_seqs; 
  BIOSEQ *seq;


  sprintf(seq_tmpfile,".blastlib.%ld.seq", (ULINT) current_pid);  
  sprintf(xml_tmpfile,".blastlib.%ld.xml", (ULINT) current_pid);  
  sprintf(prog_tmpfile,"blastlib.%ld.progress", (ULINT) current_pid);
  BLP->seq_tmpfile = seq_tmpfile;
  BLP->xml_tmpfile = xml_tmpfile;

  BLP->current_run = 0;
  fastadb_get_noseqs(s_db, &no_seqs);

  if (BLP->paralel_total > 0)
    {
      BLP->first_seq = 
	no_seqs * BLP->paralel_instance / BLP->paralel_total;
      BLP->last_seq = 
	(no_seqs * (BLP->paralel_instance+1) / BLP->paralel_total)-1;   
 
      BLP->current_seq = BLP->first_seq;
      while (BLP->current_seq <= BLP->last_seq)
	{
	  fastadb_get_seq(s_db, BLP->current_seq, &seq);
	  bioseq_parse(seq, BLP->seq_tmpfile,
		       s_db->retrieve_deflines); 	     
	  run_blast(BP, BLP);
	  BLP->current_seq++;
	  BLP->current_run++;
	  t1 = time(NULL);
	  BLP->current_time = difftime(t1,t0)/60;
	  if ((BLP->current_run+1)%BLP->progress_report_freq == 0)
	    print_progress(prog_tmpfile, BLP, 0);
#ifdef DEBUG
	  printf("Run #%ld, Point #%ld\n", BLP->current_run,
		 BLP->current_seq); 
#endif
	}
    }
  else
    {
      srand48(time(NULL));
      
      while (BLP->current_time < BLP->max_time && 
	     BLP->current_run < BLP->max_seqs)
	{
	  BLP->current_seq = lrand48()%no_seqs;
	  fastadb_get_seq(s_db, BLP->current_seq, &seq);
	  bioseq_parse(seq, BLP->seq_tmpfile,
		       s_db->retrieve_deflines); 	     
	  run_blast(BP, BLP);
	  BLP->current_run++;
	  t1 = time(NULL);
	  BLP->current_time = difftime(t1,t0)/60;
	  if ((BLP->current_run+1)%BLP->progress_report_freq == 0)
	    print_progress(prog_tmpfile, BLP, 0);
#ifdef DEBUG
	  printf("Run #%ld, Point #%ld\n", BLP->current_run,
		 BLP->current_seq); 
#endif
	}
    }
  /* Print final progress report */
  print_progress(prog_tmpfile, BLP, 1);
	
  /* Remove temporary files */
  remove(BLP->seq_tmpfile);
  remove(BLP->xml_tmpfile);
}

void *run_blast(struct blast_params *BP, struct blastlib_params *BLP)
{
  char  blast_command[1024];

  /* Parse blast command */
  /* We assume default matrices - BLOSUM62 for proteins
     -3,1 for DNA */
  sprintf(&blast_command[0], "./blastall -p %s -d %s -i %s -o %s "
	  "-e %f -m 7 -F %s -G %d -E %d -X %d -v %d -b %d -f %d "
	  "-g %s -J %s -W %d -z %f -K %d -P %d -Y %f -y %f -Z %d",
	  BP->program_name,
	  BLP->fasta_db,
	  BLP->seq_tmpfile,
	  BLP->xml_tmpfile,
	  BP->expectation_value,
	  BP->filter_query_sequence,
	  BP->cost_to_open_gap, 
	  BP->cost_to_extend_gap,
	  BP->X_dropoff_value,
	  BP->no_one_line_descriptions,
	  BP->no_allignments,
	  BP->threshold_for_extending_hits, 
	  BP->gapped_alignment,
	  BP->believe_defline,
	  BP->word_size, 
	  BP->effective_db_length,
	  BP->no_best_hits_to_keep,
	  BP->no_passes,
	  BP->effective_sspace_length,
	  BP->dropoff_for_extensions,
	  BP->X_dropoff_for_final_gapped_alignment);

  /* Run blast */
  system(blast_command);
  /* Process output */
  return BLP->blast_xml_parse(BLP->xml_tmpfile, BLP->parser_params);
}


static
int get_next_tag_item(FILE *stream, unsigned char **tag, int *taglen, 
                      unsigned char **value, int *valuelen, int *level)
{
  int c;
  int i;
  enum {LEADING_BLANK, OPEN_TAG, IS_VALUE, VALUE, CLOSE_TAG,
	TRAILING_BLANK}
      part;
 

  /* Ensure that at there is space for at least one character in
     arrays */
  if ((*taglen) < 1)
    {
      (*taglen) = BUFFER_ALLOC;
      *tag = mallocec(*taglen);
    }
  if ((*valuelen) < 1)
    {
      (*valuelen) = BUFFER_ALLOC;
      *value = mallocec(*valuelen);
    }
#ifndef OLD_GET_TAG_LEN

  part = LEADING_BLANK;
  while ((c = fgetc(stream)) != EOF)
    {
      if (part == LEADING_BLANK)
	{
	  if (c == '<')
	    {
	      part = OPEN_TAG;
	      (*level)++;
	      i = 0;
	    }
	}
      else if (part == OPEN_TAG)
	{
	  if (c == '>')
	    {
	      part = IS_VALUE;
	      (*tag)[i]='\0';
	      if ((*tag)[0] == '/')
		{
		  (*level)--;
		  (*level)--;
		}
	      i = 0;
	      continue;
	    }
	  if (i >= (*taglen)-1)
	    {
	      (*taglen) += BUFFER_ALLOC;
	      *tag = reallocec(*tag, (*taglen));
	    }	    
	  (*tag)[i++]=c;
	}
      else if (part == IS_VALUE)
	{
	  if (c == '\n')
	    {
	      (*value)[i++]='\0';      
	      break;
	    }
	  else
	    {
	      (*value)[i++]=c;
	      part = VALUE;
	    }
	}
      else if (part == VALUE)
	{
	  if (c == '<')
	    {
	      (*value)[i]='\0';
	      (*level)--;
	      part = CLOSE_TAG;
	      continue;
	    }
	  if (i >= (*valuelen)-1)
	    {
	      (*valuelen) += BUFFER_ALLOC;
	      *value = reallocec(*value, (*valuelen));
	    }		
	  (*value)[i++]=c;
	}
      else if (part == CLOSE_TAG)
	{	  
	  if (c == '\n')
	    break;
	}
    }

    if (c == EOF)
      {
	free(*tag);
	*tag = NULL;
	free(*value);
	*value = NULL;
	*taglen=0;
	*valuelen=0;
	return 0;
      }
    else
      return 1;
#else


  /* Skip leading blanks */

  while ((c=fgetc(stream)) != '<')
    if (c==EOF)
      {
	free(*tag);
	*tag = NULL;
	free(*value);
	*value = NULL;
	*taglen=0;
	*valuelen=0;
	return 0;
      }   
  (*level)++;
  /* Read opening tag */
  i = 0;
  while (1)
    {
      while (i < (*taglen)-1 && (c=fgetc(stream)) != '>')
	{
	  (*tag)[i++]=c;
	}
      if (i >= (*taglen)-1)
	{
	  (*taglen) += BUFFER_ALLOC;
	  *tag = reallocec(*tag, (*taglen));
	}
      else 
	break;
    }
  (*tag)[i]='\0';
  if ((*tag)[0] == '/')
    {
      (*level)--;
      (*level)--;
    }
  /* Read value */
  i=0;
  if ((c=fgetc(stream))=='\n')
    {
      (*value)[i++]='\0';      
      return 1;
    }
  else 
    (*value)[i++]=c;
  
  while (1)
    {
      while (i < (*valuelen)-1 && (c=fgetc(stream)) != '<')
	{
	  (*value)[i++]=c;
	}
      if (i >= (*valuelen)-1)
	{
	  (*valuelen) += BUFFER_ALLOC;
	  *value = reallocec(*value, (*valuelen));
	}
      else 
	break;
    }
  (*value)[i]='\0';
  /* Skip closing tag */
  (*level)--;
  while ((c=fgetc(stream)) != '\n');
  return 1;
#endif
}

/* ************* End of blastlib_params module ********* */





/* ************* seed_stats module ********************* */

SEED_STATS *seedstats_init(char *fasta_db, BPARAMS *BP, 
			    float *e, ULINT ne, ULINT *m, ULINT nm)
{
  ULINT i, j;
  SLINT max_S , min_S, max_D;
  

  SEED_STATS* Stats = callocec(1, sizeof(SEED_STATS));

  Stats->ne = ne;
  Stats->e = e;
  Stats->nm = nm;
  Stats->m = m;

  Stats->fasta_db = fasta_db;
  Stats->matrix = BP->matrix;
  Stats->SS = smatrix_load(BP->matrix);
  Stats->alpha = BP->cost_to_open_gap;
  Stats->beta = BP->cost_to_extend_gap; 

#if DEBUG > 2 
  max_S = 20;
  min_S = 10;
  max_D = 1;
#else
  max_S = smatrix_max(Stats->SS);
  min_S = -Stats->alpha;
  max_D = max_S + 2*Stats->alpha;
#endif

  Stats->Ig = mallocec(sizeof(SEED_HIST *) * ne);
  Stats->Mg = mallocec(sizeof(SEED_HIST *) * ne);
  Stats->Sg = mallocec(sizeof(SEED_HIST *) * ne);
  Stats->Dg = mallocec(sizeof(SEED_HIST *) * ne);
  Stats->Iu = mallocec(sizeof(SEED_HIST *) * ne);
  Stats->Mu = mallocec(sizeof(SEED_HIST *) * ne);
  Stats->Su = mallocec(sizeof(SEED_HIST *) * ne);
  Stats->Du = mallocec(sizeof(SEED_HIST *) * ne);

  for (i=0; i < ne; i++)
    {
      Stats->Ig[i] = callocec(sizeof(SEED_HIST), nm);
      Stats->Mg[i] = callocec(sizeof(SEED_HIST), nm);
      Stats->Sg[i] = callocec(sizeof(SEED_HIST), nm);
      Stats->Dg[i] = callocec(sizeof(SEED_HIST), nm);
      Stats->Iu[i] = callocec(sizeof(SEED_HIST), nm);
      Stats->Mu[i] = callocec(sizeof(SEED_HIST), nm);
      Stats->Su[i] = callocec(sizeof(SEED_HIST), nm);
      Stats->Du[i] = callocec(sizeof(SEED_HIST), nm);      

      for (j=0; j < nm; j++)
	{
	  seedhist_init1(Stats->Ig[i]+j, 0, m[j]);
	  seedhist_init1(Stats->Mg[i]+j, 0, m[j]);
	  seedhist_init1(Stats->Sg[i]+j, m[j]*min_S, m[j]*max_S);
	  seedhist_init1(Stats->Dg[i]+j, 0, m[j]*max_D);
	  seedhist_init1(Stats->Iu[i]+j, 0, m[j]);
	  seedhist_init1(Stats->Mu[i]+j, 0, m[j]);
	  seedhist_init1(Stats->Su[i]+j, m[j]*min_S, m[j]*max_S);
	  seedhist_init1(Stats->Du[i]+j, 0, m[j]*max_D);
	}

    }
  return Stats;
}

int seedstats_clear(SEED_STATS *Stats)
{
  ULINT i, j;

  for (i=0; i < Stats->ne; i++)
    {
      for (j=0; j < Stats->nm; j++)
	{
	  seedhist_clear1(Stats->Ig[i]+j);
	  seedhist_clear1(Stats->Mg[i]+j);
	  seedhist_clear1(Stats->Sg[i]+j);
	  seedhist_clear1(Stats->Dg[i]+j);
	  seedhist_clear1(Stats->Iu[i]+j);
	  seedhist_clear1(Stats->Mu[i]+j);
	  seedhist_clear1(Stats->Su[i]+j);
	  seedhist_clear1(Stats->Du[i]+j);
	}
      free(Stats->Ig[i]);
      free(Stats->Mg[i]);
      free(Stats->Sg[i]);
      free(Stats->Dg[i]);
      free(Stats->Iu[i]);
      free(Stats->Mu[i]);
      free(Stats->Su[i]);
      free(Stats->Du[i]);
    }
  free(Stats->Ig);
  free(Stats->Mg);
  free(Stats->Sg);
  free(Stats->Dg);
  free(Stats->Iu);
  free(Stats->Mu);
  free(Stats->Su);
  free(Stats->Du);

  if (Stats->SS+65 != NULL)
    smatrix_remove(Stats->SS);

  free(Stats);
  return 1;
}


void seedstats_write(SEED_STATS *Stats, char *filename)
{
  ULINT i, j, n;
  FILE *stream = fopen(filename, "wb");

  if (stream == NULL)
    {
      fprintf(stderr, "Unable to open %s\n ", filename);
      exit(EXIT_FAILURE);
    }

  n = strlen(Stats->fasta_db);
  fwrite(&n, sizeof(ULINT), 1, stream);
  fwrite(Stats->fasta_db, sizeof(char), n, stream);
  n = strlen(Stats->matrix);
  fwrite(&n, sizeof(ULINT), 1, stream);
  fwrite(Stats->matrix, sizeof(char), n, stream);
  fwrite(&Stats->alpha, sizeof(ULINT), 1, stream);
  fwrite(&Stats->beta, sizeof(ULINT), 1, stream);

  fwrite(&(Stats->ne), sizeof(ULINT), 1, stream);
  fwrite(Stats->e, sizeof(ULINT), Stats->ne, stream);
  fwrite(&(Stats->nm), sizeof(ULINT), 1, stream);
  fwrite(Stats->m, sizeof(ULINT), Stats->nm, stream);

  for (i=0; i < Stats->ne; i++)
    for (j=0; j < Stats->nm; j++)
      {
	seedhist_write(Stats->Ig[i]+j, stream);
	seedhist_write(Stats->Mg[i]+j, stream);
	seedhist_write(Stats->Sg[i]+j, stream);
	seedhist_write(Stats->Dg[i]+j, stream);
	seedhist_write(Stats->Iu[i]+j, stream);
	seedhist_write(Stats->Mu[i]+j, stream);
	seedhist_write(Stats->Su[i]+j, stream);
	seedhist_write(Stats->Du[i]+j, stream);
      }
  fclose(stream);
}

SEED_STATS *seedstats_read(char *filename)
{
  ULINT i, j, n;
  FILE *stream = fopen(filename, "rb");
  SEED_STATS* Stats;

  if (stream == NULL)
    {
      fprintf(stderr, "Unable to open %s\n ", filename);
      exit(EXIT_FAILURE);
    }

  Stats = callocec(1, sizeof(SEED_STATS));
  Stats->SS = NULL;

  fread(&n, sizeof(ULINT), 1, stream);
  Stats->fasta_db = mallocec(n+1);
  Stats->fasta_db[n] = '\0';
  fread(Stats->fasta_db, sizeof(char), n, stream);

  fread(&n, sizeof(ULINT), 1, stream);
  Stats->matrix = mallocec(n+1);
  Stats->matrix[n] = '\0';
  fread(Stats->matrix, sizeof(char), n, stream);

  fread(&Stats->alpha, sizeof(ULINT), 1, stream);
  fread(&Stats->beta, sizeof(ULINT), 1, stream);

  fread(&(Stats->ne), sizeof(ULINT), 1, stream);
  Stats->e = mallocec(Stats->ne*sizeof(float));
  fread(Stats->e, sizeof(ULINT), Stats->ne, stream);

  fread(&(Stats->nm), sizeof(ULINT), 1, stream);
  Stats->m = mallocec(Stats->nm*sizeof(ULINT));
  fread(Stats->m, sizeof(ULINT), Stats->nm, stream);

  Stats->Ig = mallocec(sizeof(SEED_HIST *) * Stats->ne);
  Stats->Mg = mallocec(sizeof(SEED_HIST *) * Stats->ne);
  Stats->Sg = mallocec(sizeof(SEED_HIST *) * Stats->ne);
  Stats->Dg = mallocec(sizeof(SEED_HIST *) * Stats->ne);
  Stats->Iu = mallocec(sizeof(SEED_HIST *) * Stats->ne);
  Stats->Mu = mallocec(sizeof(SEED_HIST *) * Stats->ne);
  Stats->Su = mallocec(sizeof(SEED_HIST *) * Stats->ne);
  Stats->Du = mallocec(sizeof(SEED_HIST *) * Stats->ne);

  for (i=0; i < Stats->ne; i++)
    {
      Stats->Ig[i] = callocec(sizeof(SEED_HIST), Stats->nm);
      Stats->Mg[i] = callocec(sizeof(SEED_HIST), Stats->nm);
      Stats->Sg[i] = callocec(sizeof(SEED_HIST), Stats->nm);
      Stats->Dg[i] = callocec(sizeof(SEED_HIST), Stats->nm);
      Stats->Iu[i] = callocec(sizeof(SEED_HIST), Stats->nm);
      Stats->Mu[i] = callocec(sizeof(SEED_HIST), Stats->nm);
      Stats->Su[i] = callocec(sizeof(SEED_HIST), Stats->nm);
      Stats->Du[i] = callocec(sizeof(SEED_HIST), Stats->nm);      

      for (j=0; j < Stats->nm; j++)
	{
	  seedhist_read(Stats->Ig[i]+j, stream);
	  seedhist_read(Stats->Mg[i]+j, stream);
	  seedhist_read(Stats->Sg[i]+j, stream);
	  seedhist_read(Stats->Dg[i]+j, stream);
	  seedhist_read(Stats->Iu[i]+j, stream);
	  seedhist_read(Stats->Mu[i]+j, stream);
	  seedhist_read(Stats->Su[i]+j, stream);
	  seedhist_read(Stats->Du[i]+j, stream);
	}

    }
  fclose(stream);
  return Stats;
}


void seedstats_print(SEED_STATS *Stats, FILE *stream)
{
  ULINT i, j;

  fprintf(stream, "******* BLAST seed statistics *******\n");
  fprintf(stream, "Database: %s\n", Stats->fasta_db);
  fprintf(stream, "Matrix: %s\n", Stats->matrix);
  fprintf(stream, "Open: %ld\n", Stats->alpha);
  fprintf(stream, "Extend: %ld\n", Stats->beta);
  fprintf(stream, "Size: %ld %ld\n\n", Stats->ne, Stats->nm);

  for (i=0; i <  Stats->ne; i++)
    {
      for (j=0; j < Stats->nm; j++)
      {
	fprintf(stream, "Expect: %8.2e\n",Stats->e[i]);
	fprintf(stream, "Length: %ld\n\n",Stats->m[j]);
       
	fprintf(stream, "\"Hamming Distance - Gapped\"\n");
	seedhist_print(Stats->Ig[i]+j ,stream);
	fprintf(stream, "\"Hamming Distance - Ungapped\"\n");
	seedhist_print(Stats->Iu[i]+j ,stream);

	fprintf(stream, "\"Modified Hamming Distance - Gapped\"\n");
	seedhist_print(Stats->Mg[i]+j ,stream);
	fprintf(stream, "\"Modified Hamming Distance - Ungapped\"\n");
	seedhist_print(Stats->Mu[i]+j ,stream);

	fprintf(stream, "\"%s Similarity - Gapped\"\n", Stats->matrix);
	seedhist_print(Stats->Sg[i]+j ,stream);
	fprintf(stream, "\"%s Similarity - Ungapped\"\n", Stats->matrix);
	seedhist_print(Stats->Su[i]+j ,stream);
	
	fprintf(stream, "\"%s Distance - Gapped\"\n", Stats->matrix);
	seedhist_print(Stats->Dg[i]+j ,stream);
	fprintf(stream, "\"%s Distance - Ungapped\"\n", Stats->matrix);
	seedhist_print(Stats->Du[i]+j ,stream);
      }
    }
}

void seedstats_collect(int no_filenames, char **in_filenames, 
		       char *out_filename)
{
  ULINT i, j, k;
  ULINT Err_flag = 0;
  SEED_STATS *Stats;
  SEED_STATS *Stats0; 

  /* Perhaps (quite) a bit inefficient because we allocate and free
     seed_stats all the time. */

  if (no_filenames < 1) return;

  Stats = seedstats_read(in_filenames[0]);

  for (k=1; k < no_filenames; k++)
    {
      Stats0 = seedstats_read(in_filenames[k]);

      /* We strictly enforce consistency by checking database name
	 and matrix */

      if ((strcmp(Stats->fasta_db, Stats0->fasta_db) != 0) ||
	  (strcmp(Stats->matrix, Stats0->matrix) != 0))
	{
	  fprintf(stderr, "Error: "
		  "inconsistent database or matrix name in %s.\n",
		  in_filenames[k]);
	  exit(EXIT_FAILURE);
	}

      if ((Stats0->nm != Stats->nm) || (Stats0->ne != Stats->ne) ||
	  (Stats0->alpha != Stats->alpha) ||
	  (Stats0->beta != Stats->beta))
	{
	  Err_flag = 1;
	}
      else
	{
	  for (i=0; i < Stats->nm; i++)
	    if (Stats->m[i] != Stats0->m[i]) 
	      {
		Err_flag = 1;
		break;
	      }
	  for (i=0; i < Stats->ne; i++)
	    if (Stats->e[i] != Stats0->e[i]) 
	      {
		Err_flag = 1;
		break;
	      }
	}
      if (Err_flag)
	{
	  fprintf(stderr, "Error: "
		  "incompatible results in %s.\n",
		  in_filenames[k]);
	  exit(EXIT_FAILURE);	  
	}

      for (i=0; i < Stats->ne; i++)
	for (j=0; j < Stats->nm; j++)
	  {
	    seedhist_merge(Stats->Ig[i]+j, Stats0->Ig[i]+j);
	    seedhist_merge(Stats->Mg[i]+j, Stats0->Mg[i]+j);
	    seedhist_merge(Stats->Sg[i]+j, Stats0->Sg[i]+j);
	    seedhist_merge(Stats->Dg[i]+j, Stats0->Dg[i]+j);
	    seedhist_merge(Stats->Iu[i]+j, Stats0->Iu[i]+j);
	    seedhist_merge(Stats->Mu[i]+j, Stats0->Mu[i]+j);
	    seedhist_merge(Stats->Su[i]+j, Stats0->Su[i]+j);
	    seedhist_merge(Stats->Du[i]+j, Stats0->Du[i]+j);
	  }
      seedstats_clear(Stats0);
    }

  if (out_filename == NULL)
    seedstats_print(Stats, stdout);
  else
    seedstats_write(Stats, out_filename);

  return;
}


void *get_seed_stat(char *out_filename, void *params)
{
  ULINT i, j;
  SLINT r;

  SEED_STATS *Stats = params;

  /* Binary search variables */
  ULINT k, low, mid;

  float evalue;

#define ALIGN_SIZE 20
  ULINT qseq_len = ALIGN_SIZE;
  unsigned char *qseq = mallocec(qseq_len);
  ULINT hseq_len = ALIGN_SIZE;
  unsigned char *hseq = mallocec(hseq_len);
  ULINT line_len = 0;

  unsigned char *qseq0;
  unsigned char *hseq0;


  /* xml field variables */
  unsigned char *tag;
  int taglen = 0;
  unsigned char *value;
  int valuelen = 0;
  int level = 0;

  /* Cumulative scores */
  ULINT aloc_size = ALIGN_SIZE;
  SLINT *cum_I = callocec(aloc_size, sizeof(SLINT));
  SLINT *cum_M = callocec(aloc_size, sizeof(SLINT));
  SLINT *cum_S = callocec(aloc_size, sizeof(SLINT));
  SLINT *cum_D = callocec(aloc_size, sizeof(SLINT));
  SLINT *cum_G = callocec(aloc_size, sizeof(SLINT));
  SLINT *cum_qL = callocec(aloc_size, sizeof(SLINT));


  /* Current scores */
  SLINT curr_I;
  SLINT curr_M;
  SLINT curr_S;
  SLINT curr_D;

  SLINT D0;

  SLINT gapS;
  ULINT last_gap;

  /* Allocate cumulative score vectors  */
  SLINT *min_Ig = mallocec(Stats->nm*sizeof(SLINT));
  SLINT *min_Mg = mallocec(Stats->nm*sizeof(SLINT));
  SLINT *max_Sg = mallocec(Stats->nm*sizeof(SLINT));
  SLINT *min_Dg = mallocec(Stats->nm*sizeof(SLINT));

  SLINT *min_Iu = mallocec(Stats->nm*sizeof(SLINT));
  SLINT *min_Mu = mallocec(Stats->nm*sizeof(SLINT));
  SLINT *max_Su = mallocec(Stats->nm*sizeof(SLINT));
  SLINT *min_Du = mallocec(Stats->nm*sizeof(SLINT));

  FILE *out_stream;

#if DEBUG > 5
  printf("******************************\n"); 
#endif
  /* Open temporary file */
  out_stream = fopen(out_filename, "r");
  if (out_stream == NULL)
    {
      fprintf(stderr, "Unable to open %s\n ", out_filename);
      exit(EXIT_FAILURE);
    }

  /* Run over all Hsps in file */

  while (get_next_tag_item(out_stream, &tag, &taglen, &value, &valuelen,
			   &level)) 
    {
      /* We're interested only in Hsp_evalue, Hsp_qseq, Hsp_hseq */ 

      /* Get evalue */
      if (strcmp("Hsp_evalue", tag) == 0)
	{
	  evalue = atof(value);
	  continue;
	}

     /* Get qseq */
     if (strcmp("Hsp_qseq", tag) == 0)
	{
	  if (qseq_len < valuelen)
	    {
#if DEBUG > 5
	      printf("Enlarging qseq. qseq_len = %ld,"
		     " value_len = %d\n", qseq_len, valuelen); 
#endif
	      qseq_len = valuelen;
	      qseq = reallocec(qseq, qseq_len);
	    }
	  memcpy(qseq, value, qseq_len);
	  line_len = strlen(qseq);
	  continue;
	}
     /* Get hseq */
     if (strcmp("Hsp_hseq", tag) == 0)
       {
	 if (hseq_len < valuelen)
	   {
#if DEBUG > 5
	      printf("Enlarging hseq. hseq_len = %ld,"
		     " value_len = %d\n", hseq_len, valuelen); 
#endif
	     hseq_len = valuelen;
	     hseq = reallocec(hseq, hseq_len);

	     aloc_size = valuelen+1;
	     cum_I = reallocec(cum_I, aloc_size*sizeof(SLINT));
	     cum_M = reallocec(cum_M, aloc_size*sizeof(SLINT));
	     cum_S = reallocec(cum_S, aloc_size*sizeof(SLINT));
	     cum_D = reallocec(cum_D, aloc_size*sizeof(SLINT));
	     cum_G = reallocec(cum_G, aloc_size*sizeof(SLINT));
	     cum_qL = reallocec(cum_qL, aloc_size*sizeof(SLINT));
	     cum_I[0] = 0;
	     cum_M[0] = 0;
	     cum_S[0] = 0;
	     cum_D[0] = 0;
	     cum_G[0] = 0;
	     cum_qL[0] = 0;
#if DEBUG > 5
	      printf("Enlarging lines."
		     " aloc_size = %d\n", aloc_size); 
#endif
	   }	

	 memcpy(hseq, value, hseq_len);

	 /* Process data */

	 /* Find k such that Stats->e[k-1] < evalue <= Stats->e[k] 
	    by binary search. Assume Stats->e[-1] = 0.0 */
	 k = Stats->ne-1;
	 low = 0;

	 if (evalue > Stats->e[k]) 
	   continue; /* Skip if not significant enough */
	 else if (evalue <= Stats->e[low])
	   k = 0;
	 else /* Binary search */
	   while (low < k-1)
	     {
	       mid = (k+low)/2;
	       if (evalue > Stats->e[mid])  
		 low = mid;
	       else
		 k = mid; 	  
	     }
	  
	  /* Set current and minimum distance vectors */
	 for (i = 0; i < Stats->nm; i++)
	   {
	     min_Ig[i] = LONG_MAX;
	     min_Mg[i] = LONG_MAX;
	     max_Sg[i] = LONG_MIN;
	     min_Dg[i] = LONG_MAX;
	     min_Iu[i] = LONG_MAX;
	     min_Mu[i] = LONG_MAX;
	     max_Su[i] = LONG_MIN;
	     min_Du[i] = LONG_MAX;	
	   }

	 /* Now calculate the minimum distance or similarity score 
	    (over the alignment) of fragments of lengths Stats->m */
	 qseq0 = qseq-1;
	 hseq0 = hseq-1;
	 for (i=1; i <= line_len; i++)
	   {
	     cum_I[i] = cum_I[i-1];
	     cum_M[i] = cum_M[i-1];
	     cum_S[i] = cum_S[i-1];
	     cum_D[i] = cum_D[i-1];
	     cum_G[i] = cum_G[i-1];
	     cum_qL[i] = cum_qL[i-1];
	     
	     /* For now, the modified distance is the same as Hamming
		distance */  
	     if (qseq0[i] != hseq0[i])
	       {
		 cum_I[i] += 1;
	       }		  
	     
	     if (qseq0[i] != '-')
	       {
		 cum_qL[i] += 1;
		 if (hseq0[i] != '-') /* No gap */
		   {
		     last_gap = 0;
		     cum_S[i] += Stats->SS[qseq0[i]][hseq0[i]];
		     D0 =  Stats->SS[qseq0[i]][qseq0[i]]+
		       Stats->SS[hseq0[i]][hseq0[i]]-
		       2*Stats->SS[qseq0[i]][hseq0[i]];
		     cum_D[i] += D0;
		     if (D0 > 7)
		       cum_M[i] += 1;
		   }
		 else
		   {
		     cum_M[i] += 1;
		     cum_G[i] += 1;
		     if (last_gap == 0) 
		       gapS = -Stats->alpha;
		     else
		       gapS = -Stats->beta;
		     last_gap = 1;
		     cum_S[i] += gapS;
		     cum_D[i] += Stats->SS[qseq0[i]][qseq0[i]] - 
		       2*gapS;
		   }
	       }
	     else
	       {
		 cum_M[i] += 1;
		 cum_G[i] += 1;
		 if (last_gap == 0) 
		   gapS = -Stats->alpha;
		 else
		   gapS = -Stats->beta;
		 last_gap = 1;
		 cum_S[i] += gapS;
		 cum_D[i] += Stats->SS[hseq0[i]][hseq0[i]] - 2*gapS; 
		 /* We are at a gap in qseq0, and so
		    we need to skip the thing altogether
		     because the gapped fragment would just be a
		     duplicate of the previous one (the ungapped
		     would skip anyway). */
		 continue;
	       }
	     
	     /* Determine actual optimal scores. */
	     for(j=0; j < Stats->nm; j++)
	       {
		 r = i-Stats->m[j];
		 if (r >= 0)
		   {
		     curr_I = cum_I[i] - cum_I[r];
		     curr_M = cum_M[i] - cum_M[r];
		     curr_S = cum_S[i] - cum_S[r];
		     curr_D = cum_D[i] - cum_D[r];
		     
		     /* Ungapped scores */
		     if (cum_G[i] - cum_G[r] == 0)
		       {
			  if (curr_I < min_Iu[j])
			    min_Iu[j] = curr_I;
			  if (curr_M < min_Mu[j])
			    min_Mu[j] = curr_M; 
			  if (curr_S > max_Su[j])
			    max_Su[j] = curr_S;
			  if (curr_D < min_Du[j])
			    min_Du[j] = curr_D;
			}
		     else /* We have a gap ! */ 
		       /* For gapped scores, we want the length of the
			  querry fragment to be m, whatever the length
			  of alignemnt might be. If there are gaps
			  in the query sequence, we will need to extend
			  the sequence to the left if possible. We use
			 cum_qL[i], the true (i.e without gaps) length
			 of qseq0 up to i, to find r, the index of qseq0
			 such that the true length of the qseq0 fragment
			 starting at r and ending at i is in fact
			 m. */
		       

		       {
			 while ((cum_qL[i] - cum_qL[r] < Stats->m[j])
				&& r > 0)
			   r--;
			 /* If r is unchanged this is unnecessary */
			 curr_I = cum_I[i] - cum_I[r];
			 curr_M = cum_M[i] - cum_M[r];
			 curr_S = cum_S[i] - cum_S[r];
			 curr_D = cum_D[i] - cum_D[r];
		       }
		     
		     /* Ungapped scores */
		     if (cum_qL[i] - cum_qL[r] == Stats->m[j])
		       /* Otherwise, r == -1 and 
			  cum_qL[i] - cum_qL[r] < Stats->m[j] so we
			  skip */
		       
		       {
			 /* Correction for affine gaps */
			 if ((cum_G[r+1] - cum_G[r] > 0) && 
			     (cum_S[r+1] - cum_S[r] != -Stats->alpha))
			   {
			     curr_S += Stats->beta - Stats->alpha;
			     /*curr_D += 2*(Stats->alpha - Stats->beta);*/
			   }
			 
			 if (curr_I < min_Ig[j])
			   min_Ig[j] = curr_I;
			 if (curr_M < min_Mg[j])
			   min_Mg[j] = curr_M; 
			 if (curr_S > max_Sg[j])
			   max_Sg[j] = curr_S;
			 if (curr_D < min_Dg[j])
			   min_Dg[j] = curr_D;   
		       }
		   }
	       }
	   }
	 
	  /* Update histograms */
	 for (i=k; i < Stats->ne; i++) 
	   for (j=0; j < Stats->nm; j++)
	     {
#if DEBUG > 6
	       printf("%s\n", qseq);
	       printf("%s\n", hseq);
#endif 
#if DEBUG > 5
	       printf("min_Iu[%d] = %d\n",j, min_Iu[j]); 
	       printf("min_Ig[%d] = %d\n\n",j, min_Ig[j]);
	       printf("min_Mu[%d] = %d\n",j, min_Mu[j]); 
	       printf("min_Mg[%d] = %d\n\n",j, min_Mg[j]);
	       printf("max_Su[%d] = %d\n",j, max_Su[j]); 
	       printf("max_Sg[%d] = %d\n\n",j, max_Sg[j]);
	       printf("min_Du[%d] = %d\n",j, min_Du[j]); 
	       printf("min_Dg[%d] = %d\n\n",j, min_Dg[j]);
#endif
 
#ifndef SD_ONLY
	       if (min_Ig[j] != LONG_MAX)
		 seedhist_add_acount(Stats->Ig[i]+j, min_Ig[j]);
	       else
		 seedhist_add_ucount(Stats->Ig[i]+j);
	       if (min_Iu[j] != LONG_MAX)
		 seedhist_add_acount(Stats->Iu[i]+j, min_Iu[j]);
	       else
		 seedhist_add_ucount(Stats->Iu[i]+j);
	       if (min_Mg[j] != LONG_MAX)
		 seedhist_add_acount(Stats->Mg[i]+j, min_Mg[j]);
	       else
		 seedhist_add_ucount(Stats->Mg[i]+j);
	       if (min_Mu[j] != LONG_MAX)
		 seedhist_add_acount(Stats->Mu[i]+j, min_Mu[j]);
	       else
		 seedhist_add_ucount(Stats->Mu[i]+j);
#endif
	       if (max_Sg[j] != LONG_MIN)
		 seedhist_add_acount(Stats->Sg[i]+j, max_Sg[j]);
	       else
		  seedhist_add_ucount(Stats->Sg[i]+j);
	       if (max_Su[j] != LONG_MIN)
		 seedhist_add_acount(Stats->Su[i]+j, max_Su[j]);
	       else
		 seedhist_add_ucount(Stats->Su[i]+j);
	       if (min_Dg[j] != LONG_MAX)
		 seedhist_add_acount(Stats->Dg[i]+j, min_Dg[j]);
	       else
		 seedhist_add_ucount(Stats->Dg[i]+j);
		if (min_Du[j] != LONG_MAX)
		  seedhist_add_acount(Stats->Du[i]+j, min_Du[j]);
		else
		  seedhist_add_ucount(Stats->Du[i]+j);
	      }
       }
    }
  
  fclose(out_stream);

  /* Free memory and exit */
  free(min_Ig);
  free(min_Mg);
  free(max_Sg);
  free(min_Dg);
  free(min_Iu);
  free(min_Mu);
  free(max_Su);
  free(min_Du);

  free(cum_I);
  free(cum_M);
  free(cum_S);
  free(cum_D);
  free(cum_G);
  free(cum_qL);

  free(qseq);
  free(hseq);

  return NULL;
}
