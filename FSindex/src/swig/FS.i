%module FS

%{
#include "misclib.h"
#include "bioseq.h"
#include "fastadb.h"
#include "partition.h"
#include "smatrix.h"
#include "pmatrix.h"
#include "randseq.h"
#include "hit_list.h"
#include "FSindex.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;

struct exception_context the_exception_context[1];

 typedef FS_PARTITION_t ptable;
 typedef SEQUENCE_DB db;
 typedef SCORE_MATRIX_t smatrix;
 typedef SEQ_HIT_t hit;
 typedef HIT_LIST_t hit_list;
 typedef FSINDX index;


  typedef struct {
    BIOSEQ *seq;
    int free_flag;
    int defline_flag;
  } seqn;

  int seqn_len_get(seqn *self) {
    return self->seq->len;
  }
  const char *seqn_seq_get(seqn *self) {
    return self->seq->start;
  }
  const char *seqn_defline_get(seqn *self) {
    if (!self->defline_flag || self->seq->id.defline == NULL)
      return "";
    else 
      return self->seq->id.defline;
  }

  int ptable_no_pttn_get(ptable *self) {
    return FS_PARTITION_get_no_partitions(self);
  }

  ptable *ptable_read(FILE *stream) {
    return FS_PARTITION_read(stream);
  }

  typedef struct {
    SEQ_GENERATOR_t *sgen;
    char *heap;
    ULINT max_len;
    seqn *Fseq;
  } fgen;

  index *index_load(const char *filename) {
   return FS_INDEX_load(filename);
  }

  hit_list *sscan_qd_srch(db *s_db, const char *matrix, 
			  seqn *query, ULINT D_cutoff) {
    return SSCAN_QD_search(s_db, matrix, query->seq, 
			  D_cutoff);
  }
 int sscan_has_nbr(db *s_db, smatrix *D, ptable *pt, 
		   seqn *query, ULINT D_cutoff) {
   return SSCAN_has_neighbour(s_db, D, pt, query->seq, D_cutoff);
 }



%}

%typemap(in) FILE *stream {
  $1 = PyFile_AsFile($input);
}


%exception {
  Try {
    $action
  }
  Catch(except) {
    switch (except->code)
      {
      case NO_ERR:
	break;
      case NO_MEM:
	PyErr_SetString(PyExc_MemoryError, except->msg);
	break;
      case NEG_MEM_REQ:
	PyErr_SetString(PyExc_RuntimeError, except->msg);
	break;
      case FOPEN_ERR:
	PyErr_SetString(PyExc_IOError, except->msg);
	break;
      case FCLOSE_ERR:
	PyErr_SetString(PyExc_IOError, except->msg);
	break;
      case BAD_ARGS:
	PyErr_SetString(PyExc_RuntimeError, except->msg);
	break;
      case EOF_REACHED:
	PyErr_SetString(PyExc_IOError, except->msg);
	break;
      case GETLINE_ERR:
	PyErr_SetString(PyExc_IOError, except->msg);
	break;
      case INDEX_OUT_OF_RANGE:
	PyErr_SetString(PyExc_IndexError, except->msg);
	break;
      default:
	PyErr_SetString(PyExc_RuntimeError, except->msg);
	break;
      }
    return NULL;
  }
}

%nodefault;

typedef unsigned long ULINT;
typedef signed short SSINT;


/********************************************************************/ 
/*         seqn                                                    */
/********************************************************************/ 

typedef struct {} 
seqn;


%extend seqn {
  seqn(char *sqn, char *defline) {
    seqn *Fseq;
    
    Fseq = mallocec(sizeof(seqn));
    Fseq->seq = mallocec(sizeof(BIOSEQ));
    Fseq->seq->start = sqn;
    Fseq->seq->len = strlen(sqn);
    Fseq->seq->id.defline = defline;
    Fseq->free_flag = 1;
    Fseq->defline_flag = 1;
    return Fseq;
  }
  ~seqn() {
    if (self->free_flag)
      free(self->seq);
    free(self);
  }
  const char *__str__() {
    return self->seq->start;
  }
  %immutable;
  int len;
  const char *seq;
  const char *defline;
}

/********************************************************************/ 
/*         db                                                       */
/********************************************************************/ 

typedef struct
{
%immutable;
  char *db_name; 
  ULINT no_seq;  
} db;


%extend db {
  db(const char *db_name) {
    fastadb_arg fastadb_argt[3];
    fastadb_argv_t fastadb_argv[3];

    fastadb_argt[0] = ACCESS_TYPE;
    fastadb_argt[1] = RETREIVE_DEFLINES;
    fastadb_argt[2] = NONE;
    fastadb_argv[0].access_type = MEMORY;
    fastadb_argv[1].retrieve_deflines = YES;

    return fastadb_open(db_name, fastadb_argt, fastadb_argv);
  }
  ~db() {
    fastadb_close(self);
  }
  %newobject __str__;
  char *__str__() {
    char *s;
    ULINT len;
    len = fastadb_count_Ffrags(self, 1);
    s = mallocec(strlen(self->db_name)+7+20+9+10+12);
    sprintf(s,"File: %s\nLength: %d\nSequences: %d",
	    self->db_name, len, self->no_seq);
    return s;
  }
  seqn *get_seq(ULINT seq_no) {
    seqn *Fseq;
    Fseq = mallocec(sizeof(seqn));
    Fseq->free_flag = 0;
    Fseq->defline_flag = 1;

    fastadb_get_seq(self, seq_no, &(Fseq->seq));
    return Fseq;
  }
  /* Fragment functions */
  void init_frags(ULINT min_len, ULINT max_len) {
    fastadb_init_frags(self, min_len, max_len);
  }
  void clear_frags() {
     fastadb_clear_frags(self);
  }
  ULINT get_nofrags(ULINT min_len, ULINT max_len) {
    ULINT n = 0;

    fastadb_get_nofrags(self, &n, min_len, max_len);
    return n;
  }
  seqn *get_frag(ULINT n, ULINT min_len, ULINT max_len) {
    seqn *Fseq;
    BIOSEQ *ts;

    Fseq = mallocec(sizeof(seqn));
    Fseq->free_flag = 1;
    Fseq->defline_flag = 0;

    fastadb_get_frag(self, ts, n, min_len, max_len);
    Fseq->seq = bioseq_copy(ts);
    return Fseq;
  }

  /* Ffragment functions */
  void init_Ffrags(ULINT len) {
    fastadb_init_Ffrags(self, len);
  }
  ULINT count_Ffrags(ULINT len) {
    ULINT c;
    c = fastadb_count_Ffrags(self, len);
    return c;
  }
  seqn *get_Ffrag(ULINT len, ULINT n) {
    seqn *Fseq;
    BIOSEQ *ts;
    Fseq = mallocec(sizeof(seqn));
    Fseq->free_flag = 1;
    Fseq->defline_flag = 0;

    fastadb_get_Ffrag(self, len, ts, n);
    Fseq->seq = bioseq_copy(ts);
    return Fseq;
  }
  /* TO DO:
     fastadb_find_Ffrag_seq()
     fastadb_get_Ffrag_seq()
  */
}

/********************************************************************/ 
/*         ptable                                                   */
/********************************************************************/ 

typedef struct {}
ptable;

%extend ptable {
  ptable(const char *alphabet) {
    return FS_PARTITION_create(alphabet, '#');
  }
  ~ptable() {
    FS_PARTITION_destroy(self);
  }
  void print_table(FILE *stream=stdout) {
    FS_PARTITION_print(self, stream);
  }
  char * __str__() {
    return self->alphabet;
  }
  %immutable;
  int no_pttn;

  int get_pttn(char letter) {
    return FS_PARTITION_get_pttn(self, letter);
  }
  int get_posn(char letter) {
    return FS_PARTITION_get_pos(self, letter);
  }
  int get_pttn_size(int pttn) {
    return FS_PARTITION_get_size(self, pttn);
  }
  int get_poffset(int pttn) {
    return FS_PARTITION_get_poffset(self, pttn);
  }
  char get_letter(int pttn, int posn) {
    return (char) FS_PARTITION_get_letter(self, pttn, posn);
  }
  int check_seqn(seqn *Fseq) {
    return FS_PARTITION_check_seq(Fseq->seq, self);
  }
  void write(FILE *stream=stdout) {
    FS_PARTITION_write(self, stream);
  }
  %newobject seqn2reduced;
  seqn *seqn2reduced(seqn *Fseq) {
    char c;
    int i;
    int j=0;
    char *newseq = mallocec(Fseq->seq->len+1);
    seqn *Rseq;

    for (j=0; j < Fseq->seq->len; j++)
      {
	c = Fseq->seq->start[j];
	i = self->partition_table[c & A_SIZE_MASK];
	if (i == -1) 
	  {
	    free(newseq);
	    Throw FSexcept(BAD_ARGS, "Sequence contains"
			   " letters not in ptable.\n");
	  }
	sprintf(newseq+j, "%1.1d", i);
      }
    newseq[j]='\0';

    Rseq = mallocec(sizeof(seqn));
    Rseq->seq = mallocec(sizeof(BIOSEQ));
    Rseq->seq->start = newseq;
    Rseq->seq->len = Fseq->seq->len;
    Rseq->seq->id.defline = Fseq->seq->id.defline;
    Rseq->free_flag = 1;
    Rseq->defline_flag = 1;
    return Rseq;
  }
}
  extern ptable *ptable_read(FILE *stream=stdin);

/********************************************************************/ 
/*         randseq                                                  */
/********************************************************************/ 

typedef struct {
%immutable;
  ULINT max_len;
} fgen;

%extend fgen {
  fgen(const char *filename, ptable *ptable, ULINT max_len=30) {
    fgen *FG = mallocec(sizeof(fgen));
    FG->sgen = SEQ_GENERATOR_create(filename, ptable);
    FG->max_len = max_len;
    FG->heap = mallocec(max_len+1);
    FG->Fseq = mallocec(sizeof(seqn));
    FG->Fseq->free_flag = 0;
    FG->Fseq->defline_flag = 1;
    FG->Fseq->seq = mallocec(sizeof(BIOSEQ));
    FG->Fseq->seq->id.defline = "Random Sequence";
    return FG;
  }
  ~fgen() {
    SEQ_GENERATOR_destroy(self->sgen);
    free(self->heap);
    free(self->Fseq->seq);
    free(self->Fseq);
    free(self);
  }
  seqn *rand_seq(ULINT len) {
    if (len >= self->max_len) {
      self->max_len = len;
      self->heap = reallocec(self->heap, self->max_len+1);
    }
    SEQ_GENERATOR_rand_seq(self->sgen, self->Fseq->seq, len, 
			   self->heap);
    return self->Fseq;
  }
  double freq(char c) {
    int i = c & A_SIZE_MASK;
    if (i==0)
      return ((double) self->sgen->cum_freq[i])/ self->sgen->total_residues;  
    else
      return ((double) (self->sgen->cum_freq[i] - self->sgen->cum_freq[i-1])) /
	self->sgen->total_residues;
  }


}

/********************************************************************/ 
/*         smatrix                                                  */
/********************************************************************/ 

typedef struct {
%immutable;
  const char *filename;
  int similarity_flag;
}
smatrix;

%extend smatrix {
  smatrix(const char *filename, ptable *ptable) {
    return SCORE_MATRIX_create(filename, ptable);
  }
  ~smatrix() {
    SCORE_MATRIX_destroy(self);
  }
  void set_M(char row, char col, SSINT val) {
    SCORE_MATRIX_set_M(self, row, col, val);
  }
  void set_SS(char row, SSINT val) {
    SCORE_MATRIX_set_SS(self, row, val);
  }
  int get_M(char row, char col) {
    return SCORE_MATRIX_get_M(self, row, col);
  }
  int get_pM(char row, int group) {
    return SCORE_MATRIX_get_pM(self, row, group);
  }
  int get_pMc(char row) {
    return SCORE_MATRIX_get_pMc(self, row);
  }
  int get_SS(char row) {
    return SCORE_MATRIX_get_SS(self, row);
  }
  int score(seqn *query, seqn *subject) {
    return SCORE_MATRIX_evaluate(self, query->seq, subject->seq);
  }
  void print_matrix(FILE *stream=stdout, 
		    const char *title="SCORE MATRIX") {
    SCORE_MATRIX_print(self, stream, title); 
  }
  smatrix *S2Dmax() {
    return SCORE_MATRIX_S_2_Dmax(self);
  }
  smatrix *S2Davg() {
    return SCORE_MATRIX_S_2_Davg(self);
  }
  smatrix *S2Dquasi() {
    return SCORE_MATRIX_S_2_Dquasi(self);
  }
}

/********************************************************************/ 
/*         hit_list                                                 */
/********************************************************************/ 

typedef struct 
{
%immutable;
  ULINT sequence_id;
  ULINT sequence_from;
  int rejected;
  float value;
  double pvalue;
  double evalue;
  double zvalue;
  int oc_cluster;
  double cratio;
  double kw_score;
} hit;

%extend hit {
  seqn *get_subject() {
    seqn *Fseq;
    
    Fseq = mallocec(sizeof(seqn));
    Fseq->seq = bioseq_copy(&self->subject);
    Fseq->free_flag = 1;
    Fseq->defline_flag = 1;
    return Fseq;
  }
}

typedef struct
{
%immutable;
  ULINT frag_len;
  db *s_db;
  const char *matrix;
  int range;
  int converted_range;
  ULINT kNN;
  const char *index_name;
  const char *alphabet;
  ULINT FS_seqs_total;
  ULINT FS_seqs_visited;
  ULINT FS_seqs_hits;
  ULINT index_seqs_total;
  ULINT seqs_visited;
  ULINT seqs_hits;
  ULINT useqs_visited;
  ULINT useqs_hits;
  double start_time;
  double end_time;
  double search_time;

  ULINT max_hits;
  ULINT actual_seqs_hits;
  ULINT accepted;

  double shape;
  double rate;
  double Zmin; 
} hit_list;

%extend hit_list {
  /* No constructor - can get it only as results of searches */
  ~hit_list() {
    HIT_LIST_destroy(self);
  }
  void print_list(FILE *stream=stdout) {
    HIT_LIST_print(self, stream, NULL);
  } 
  void print_xml(FILE *stream=stdout) {
    HIT_LIST_xml(self, stream, 0);
  }

  hit *get_hit(ULINT i) {
    return HIT_LIST_get_hit(self, i);
  }
  %newobject get_query;
  seqn *get_query() {
    seqn *Fseq;
    
    Fseq = mallocec(sizeof(seqn));
    Fseq->seq = bioseq_copy(&self->query);
    Fseq->free_flag = 1;
    Fseq->defline_flag = 1;
    return Fseq;
  }
  void Z_scores() {
    HIT_LIST_Zscores(self);
  }
  void sort_decr() {
    HIT_LIST_sort_decr(self);
  }
  void sort_incr() {
    HIT_LIST_sort_incr(self);
  }
  void sort_by_seq() {
    HIT_LIST_sort_by_sequence(self);
  }
  void sort_by_oc() {
    HIT_LIST_sort_oc(self); 
  }
  void sort_by_evalue() {
    HIT_LIST_sort_evalue(self, 0);    
  }
  void sort_by_cratio() {
    HIT_LIST_sort_cratio(self, 0);  
  }
  void sort_by_kwscore() {
    HIT_LIST_sort_kwscore(self, 0);
  }    
}

/********************************************************************/ 
/*         pmatrix                                                  */
/********************************************************************/ 




/********************************************************************/ 
/*         index                                                    */
/********************************************************************/ 

typedef struct
{
  %immutable;
  char *db_name;           /* Name of the database                  */
  char *index_name;        /* Name of the index                     */
  db *s_db;                /* Pointer to the fasta database         */
  ptable *ptable;          /* Partition table                       */
  int m;                   /* Length of indexed fragments           */
  ULINT db_no_frags;       /* Number of fragments in full database  */
  ULINT no_bins;           /* Number of bins = K^frag_len           */
  ULINT no_seqs;           /* Number of indexed fragments           */
  ULINT no_useqs;          /* Number of indexed unique fragments    */
} index;

extern index *index_load(const char *filename);


%extend index {
  index(const char *database, ULINT flen, const char *abet, 
	 int skip=1) {
    return FS_INDEX_create(database, flen, abet, '#', skip);
  }
  ~index() {
    FS_INDEX_destroy(self);
  }
  void save(const char *filename) {
    FS_INDEX_save(self, filename);
  }
  void print_stats(int options=3, FILE *stream=stdout) {
    FS_INDEX_print_stats(self, stream, options);
  }  
  void print_bin(seqn *Fseq, int options=1, FILE *stream=stdout) {
    FS_INDEX_print_bin(self, Fseq->seq, stream, options);
  }
  ULINT get_bin_size(ULINT bin) {
    return FS_INDEX_get_bin_size(self, bin); 
  }
  ULINT get_unique_bin_size(ULINT bin) {
    return FS_INDEX_get_unique_bin_size(self, bin);
  }
  %newobject get_seq;
  seqn *get_seq(ULINT bin, ULINT i) {
    seqn *Fseq;
    
    Fseq = mallocec(sizeof(seqn));
    Fseq->seq = bioseq_copy(FS_INDEX_get_seq(self, bin, i));
    Fseq->free_flag = 1;
    Fseq->defline_flag = 1;
    return Fseq;
  }
  hit_list *rng_srch(seqn *query, smatrix *D, int d0) {
    return FSINDX_rng_srch(self, query->seq, D, d0, NULL);
  }
  hit_list *kNN_srch(seqn *query, smatrix *D, int k) {
    return FSINDX_kNN_srch(self, query->seq, D, k, NULL);
  }


  /* TO DO:
     - profile searches
  */
}

extern hit_list *sscan_qd_srch(db *s_db, const char *matrix, 
			       seqn *query, ULINT D_cutoff);

extern int sscan_has_nbr(db *s_db, smatrix *D, ptable *pt, 
			 seqn *query, ULINT D_cutoff);


