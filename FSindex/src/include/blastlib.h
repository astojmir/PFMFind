#ifndef BLASTLIB_H 
#define BLASTLIB_H 1

#include <stdio.h>
#include "misclib.h"
#include "fastadb.h"
#ifdef USE_MPATROL
#include <mpatrol.h>
#endif

typedef void *blast_xml_output_func(char *out_filename, void *params);


typedef int ** SMATRIX;

SMATRIX smatrix_load(const char *filename);
void smatrix_remove(SMATRIX S);
SLINT smatrix_max(SMATRIX S);
SMATRIX smatrix2dist(SMATRIX S);
inline
int Hamming_dist(BIOSEQ *seq1, BIOSEQ *seq2, SMATRIX D);

typedef struct
  {
    SLINT no_bins;
    SLINT min_bin;
    ULINT *freq;
    ULINT attainable_pts;
    ULINT unattainable_pts;
} SEED_HIST;

SEED_HIST *seedhist_init(SLINT min_bin, SLINT max_bin);
void seedhist_init1(SEED_HIST *hist, SLINT min_bin, SLINT max_bin);
int seedhist_clear(SEED_HIST *hist);
int seedhist_clear1(SEED_HIST *hist);
inline
void seedhist_add_acount(SEED_HIST *hist, SLINT bin);
void seedhist_add_ucount(SEED_HIST *hist);
int seedhist_merge(SEED_HIST *hist1, SEED_HIST *hist2);
int seedhist_read(SEED_HIST *hist, FILE *stream);
int seedhist_write(SEED_HIST *hist, FILE *stream);
int seedhist_print(SEED_HIST *hist, FILE *stream);

        
/*#define seedhist_add_acount(hist, bin) \
        ((hist)->freq[(bin)]++), \
        ((hist)->attainable_pts++)*/

#define seedhist_add_ucount(hist) \
        ((hist)->unattainable_pts++)

typedef struct blast_params
  {
    char *program_name;
    float expectation_value;
    char *filter_query_sequence;
    int cost_to_open_gap;
    int cost_to_extend_gap;
    int X_dropoff_value;
    int nt_mismatch_penalty;
    int nt_match_reward;
    int no_one_line_descriptions;
    int no_allignments;
    int threshold_for_extending_hits;
    char *gapped_alignment;
    int query_genetic_code;
    int db_genetic_code;
    char *believe_defline;
    char *matrix;
    int word_size;
    float effective_db_length;
    int no_best_hits_to_keep;
    int no_passes;
    float effective_sspace_length;
    int query_strands_to_search;
    float dropoff_for_extensions;
    int X_dropoff_for_final_gapped_alignment;
  } BPARAMS;

BPARAMS *default_blast_params(void);
void read_blast_params(BPARAMS *BP, const char *BPfilename); 
void write_blast_params(BPARAMS *BP, const char *BPfilename);
void clean_blast_params(BPARAMS *BP);


typedef struct blastlib_params
  {
    char *fasta_db;
    ULINT first_seq; /* in thousands - sequence to start and to stop */ 
    ULINT last_seq;
    ULINT current_seq;
    ULINT current_run;
    ULINT max_seqs;
    float max_time; /* in minutes */
    float current_time;
    USINT paralel_instance;
    USINT paralel_total;
    USINT progress_report_freq;
    char *seq_tmpfile;
    char *xml_tmpfile;
    blast_xml_output_func *blast_xml_parse;
    void *parser_params;
  } BLPARAMS;

BLPARAMS *blastlib_params_init(char *fasta_db, ULINT max_seqs,
			       float max_time, USINT paralel_instance,
			       USINT paralel_total, 
			       USINT progress_report_freq,
			       blast_xml_output_func *blast_xml_parse, 
			       void *parser_params);

int blastlib_params_clear(BLPARAMS *BLP);

void blastlib_run_blast_expt(BLPARAMS *BLP, BPARAMS *BP, 
			     SEQUENCE_DB* s_db);

void *run_blast(BPARAMS *BP, BLPARAMS *BLP);






typedef struct 
  {
    char *fasta_db;
    char *matrix;
    SMATRIX SS;
    ULINT alpha;
    ULINT beta;

    ULINT nm;
    ULINT *m;
    ULINT ne;
    float *e;

    SEED_HIST **Ig;
    SEED_HIST **Mg;
    SEED_HIST **Sg;
    SEED_HIST **Dg;
    SEED_HIST **Iu;
    SEED_HIST **Mu;
    SEED_HIST **Su;
    SEED_HIST **Du;

  } SEED_STATS;

SEED_STATS *seedstats_init(char *fasta_db, BPARAMS *BP, 
			   float *e, ULINT ne, ULINT *m, ULINT nm);

int seedstats_clear(SEED_STATS *Stats);

void seedstats_write(SEED_STATS *Stats, char *filename);
SEED_STATS *seedstats_read(char *filename);
void seedstats_print(SEED_STATS *Stats, FILE *stream);
void seedstats_collect(int no_filenames, char **in_filenames, 
		       char *out_filename);
void *get_seed_stat(char *out_filename, void *params);







#endif
      
      

