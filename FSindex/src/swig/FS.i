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


%module FS


%{
#include "misclib.h"
#include "bioseq.h"
#include "fastadb.h"
#include "smatrix.h"
#include "FSindex.h"
#include "hit_list.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;

struct exception_context the_exception_context[1];

 typedef SCORE_MATRIX Smatrix;
 typedef SEQUENCE_DB db;
 typedef FSINDX Index;
 typedef SEQ_HIT_t hit;
 typedef HIT_LIST_t hit_list;
 typedef HIT_LIST_t * hit_list_array;

extern int next_avail_srch;
%}

%include cstring.i
%typemap(newfree) char * {
  free($1);
}
%typemap(in) FILE *stream {
  $1 = PyFile_AsFile($input);
}
%apply (char *STRING, int LENGTH) { (char *s1, int l1) };
%apply (char *STRING, int LENGTH) { (char *s2, int l2) };
%apply (char *STRING, int LENGTH) { (char *qseq, int qlen) };

/* Converting matrix as dictionary into M + alphabet */
%typemap(in) (int **M, char *alphabet) {
  int i, j;
  int a_set[A_SIZE];
  int a_len = 0;
  char *a = NULL;
  char *b = NULL;
  int val;
  PyObject *key, *value;
  PyObject *Pa, *Pb;
  Py_ssize_t pos = 0;

  if (!PyDict_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a dictionary");
    return NULL;
  }
  memset(a_set, 0 , A_SIZE * sizeof(int));
  $1 = callocec(A_SIZE, sizeof(int *));
  for (i=0; i < A_SIZE; i++) {
    $1[i] = callocec(A_SIZE, sizeof(int));
  }
  while (PyDict_Next($input, &pos, &key, &value)) {
    if (PyTuple_Check(key) && PyTuple_GET_SIZE(key) == 2) {
      Pa = PyTuple_GetItem(key, 0);
      Pb = PyTuple_GetItem(key, 1);
      if (PyString_Check(Pa) && PyString_Check(Pb)) {
	if (PyString_Size(Pa) == 1 && PyString_Size(Pb) == 1) {
	  a = PyString_AsString(Pa);
	  b = PyString_AsString(Pb);
	}
      }
      else continue;
      if (PyInt_Check(value)) {
	val = (int) PyInt_AsLong(value);
      }
      else if (PyFloat_Check(value)) {
	val = (int) PyFloat_AsDouble(value);
      }
      else continue;
      $1[*a & A_SIZE_MASK][*b & A_SIZE_MASK] = val;
      $1[*b & A_SIZE_MASK][*a & A_SIZE_MASK] = val;
      if (++(a_set[*a & A_SIZE_MASK]) == 1) {
	a_len++;
      }
      if (++(a_set[*b & A_SIZE_MASK]) == 1) {
	a_len++;
      }
    }
  }
  if (!a_len) {
    PyErr_SetString(PyExc_ValueError, "No valid entry given");
    for (i=0; i < A_SIZE; i++) {
      free($1[i]);
    }
    free($1);
    return NULL;
  }
  $2 = callocec(a_len+1,1);
  for (j=0, i=0; i < A_SIZE; i++)
    if (a_set[i]) {
      $2[j++] = 64+i;
    }
}

/* Converting pssm as list of dictionaries into M + alphabet + len */
%typemap(in) (int **M, int len, char *alphabet) {
  int i, j;
  int a_set[A_SIZE];
  int a_len = 0;
  char *a;
  int val;
  PyObject *key, *value;
  PyObject *row;
  PyObject *dict;
  Py_ssize_t pos = 0;

  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  memset(a_set, 0 , A_SIZE * sizeof(int));
  $2 = PyList_GET_SIZE($input);
  $1 = callocec($2, sizeof(int *));

  for (i=0; i < $2; i++) {
    $1[i] = callocec(A_SIZE, sizeof(int));
    row = PyList_GET_ITEM($input, i);
    if (!(PyTuple_Check(row) && PyTuple_GET_SIZE(row) == 2))
      continue;
    dict = PyTuple_GetItem(row, 1);
    if (!PyDict_Check(dict)) continue;
    pos = 0;
    while (PyDict_Next(dict, &pos, &key, &value)) {
      if (PyString_Check(key) && PyString_GET_SIZE(key)==1)
	a = PyString_AsString(key);
      else continue;
      if (PyInt_Check(value))
	val = (int) PyInt_AsLong(value);
      else if (PyFloat_Check(value))
	val = (int) PyFloat_AsDouble(value);
      else continue;
      $1[i][*a & A_SIZE_MASK] = val;
      if (++(a_set[*a & A_SIZE_MASK]) == 1)
	a_len++;
    }
  }
  if (!a_len) {
    PyErr_SetString(PyExc_ValueError, "No valid entry given");
    for (i=0; i < $2; i++)
      free($1[i]);
    free($1);
    return NULL;
  }
  $3 = callocec(a_len+1,1);
  for (j=0, i=0; i < A_SIZE; i++)
    if (a_set[i]) $3[j++] = 64+i;
}

/* Converting list of strings into char ** */
%typemap(in) (ULINT len, const char **sepn) {
  int i;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $1 = PyList_Size($input);
  $2 = (char **) malloc($1*sizeof(char *));
  for (i = 0; i < $1; i++) {
    PyObject *s = PyList_GetItem($input,i);
    if (!PyString_Check(s)) {
        free($2);
        PyErr_SetString(PyExc_ValueError, "List items must be strings");
        return NULL;
    }
    $2[i] = PyString_AsString(s);
  }
}

%typemap(freearg) (int len, const char **pttn) {
   if ($2) free($2);
}

%typemap(in,numinputs=0) (ULINT *len, char ***sepn) (ULINT temp_len, char **temp_sepn) {
  $1 = &temp_len;
  $2 = &temp_sepn;
}

%typemap(argout) (ULINT *len, char ***sepn) {
 int i;

 $result = PyList_New(*$1);
 for (i=0; i < *$1; i++) {
   PyObject *s = PyString_FromString((*$2)[i]);
   PyList_SetItem($result, i, s);
 }
}
/* Returning a string + size - without allocation */
%typemap(in,numinputs=0) (char **s, int *l) (char *stemp, int ltemp) {
  $1 = &stemp;
  $2 = &ltemp;
}
%typemap(argout) (char **s, int *l) {
  $result = PyString_FromStringAndSize(*$1, *$2);
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
      case FOPEN_ERR:
      case FCLOSE_ERR:
      case EOF_REACHED:
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
/*         db                                                       */
/********************************************************************/

typedef struct
{
%immutable;
  char *db_name;
  ULINT length;
  ULINT no_seq;
} db;


%extend db {
  db(const char *db_name) {
    return fastadb_open(db_name);
  }
  ~db() {
    fastadb_close(self);
  }
  %pythoncode %{
  def __str__(self):
    return 'FASTA SEQUENCE DATABASE\nFile: %s\nLength: %d\nSequences: %d\n' % (self.db_name, self.length , self.no_seq)
  %}
 const char *get_seq(ULINT i) {
   if (i < self->no_seq)
     return self->seq[i].start;
   else
     return NULL;
 }
 const char *get_def(ULINT i) {
   if (i < self->no_seq)
     return self->seq[i].id.defline;
   else
     return NULL;
 }
 /* No error checking here - this function is meant only
    for getting the sequences of hits */
 %newobject get_frag;
 PyObject *get_frag(ULINT i, ULINT from, ULINT to) {
   int len = to - from;
   return PyString_FromStringAndSize(self->seq[i].start + from, len);
 }
}


/********************************************************************/
/*         smatrix                                                  */
/********************************************************************/

typedef enum {SCORE, POSITIONAL} MATRIX_TYPE;
typedef enum {SIMILARITY, DISTANCE} SCORE_TYPE;

#define POS 1
#define QUASI 2
#define MAX 4
#define AVG 6

typedef struct {
%immutable;
  MATRIX_TYPE Mtype;
  SCORE_TYPE Stype;
  int len;
  char *qseq;
  char *alphabet;
%mutable;
}
Smatrix;

%newobject Smatrix_load;
extern Smatrix *Smatrix_load(const char *filename);
%newobject Smatrix_pssm;
extern Smatrix *Smatrix_pssm(int **M, int len, char *alphabet);
%{
  Smatrix *Smatrix_load(const char *filename) {
    return SCORE_MATRIX_from_file(filename);
  }
  Smatrix *Smatrix_pssm(int **M, int len, char *alphabet) {
    return SCORE_MATRIX_init(M, len, POSITIONAL, SIMILARITY, alphabet);
  }
%}

%extend Smatrix {
  Smatrix(int **M, char *alphabet, SCORE_TYPE Stype=SIMILARITY) {
    return SCORE_MATRIX_init(M, 0, SCORE, Stype, alphabet);
  }
  ~Smatrix() {
    SCORE_MATRIX_del(self);
  }
  %newobject copy;
  Smatrix *copy() {
    return SCORE_MATRIX_copy(self);
  }
  %newobject __str__;
  char *__str__() {
    return SCORE_MATRIX_sprint(self);
  }
  void fprint(FILE *fp=stdout) {
    SCORE_MATRIX_fprint(self, fp);
  }
  int score_get_item(char c1, char c2) {
    return self->M[c1 & A_SIZE_MASK][c2 & A_SIZE_MASK];
  }
  int profile_get_item(int i1, char c2) {
    return self->M[i1][c2 & A_SIZE_MASK];
  }
  void score_set_item(char c1, char c2, int val) {
    self->M[c1 & A_SIZE_MASK][c2 & A_SIZE_MASK] = val;
  }
  void profile_set_item(int i1, char c2, int val) {
    self->M[i1][c2 & A_SIZE_MASK] = val;
  }
  int eval_score(char *s1, int l1, char *s2, int l2) {
    l2 = l1 < l2 ? l1 : l2;
    return self->eval_score(self, s1, s2, l2);
  }
  int range_conv(char *s1, int l1, int r) {
    return self->range_conv(self, s1, l1, r);
  }
  int item_conv(char c, int i, int j) {
    return self->item_conv(self, c, i, j);
  }
  %newobject matrix_conv;
  Smatrix *matrix_conv(char *s1, int l1) {
    return self->matrix_conv(self, s1, l1);
  }
  int conv_type;
}
%{
int Smatrix_conv_type_get(Smatrix *S) {
  return S->conv_type;
}
void Smatrix_conv_type_set(Smatrix *S, int ctype) {
  S->set_conv_type(S, ctype);
}
%}


/********************************************************************/
/*         Index                                                    */
/********************************************************************/

/* Conversion of hit lists into Python dictionaries */
%{
PyObject *Dict_FromHitList(hit_list *HL) {

    PyObject *dict = PyDict_New();
    PyObject *obj;
    PyObject *list;
    PyObject *hdict;
    SEQ_HIT_t *hit;
    int i;

    HIT_LIST_sort_incr(HL);

    /* Query data */
    obj = PyString_FromStringAndSize(HL->query.start,
				     HL->query.len);
    PyDict_SetItemString(dict, "query_seq", obj);
    Py_DECREF(obj);

    obj = PyString_FromString(HL->query.id.defline);
    PyDict_SetItemString(dict, "query_def", obj);
    Py_DECREF(obj);

    /* Matrix is deliberately missing - the caller of
       search function should have it anyway. */

    PyDict_SetItemString(dict, "conv_type",
			 PyInt_FromLong(HL->conv_type));
    PyDict_SetItemString(dict, "sim_range",
			 PyInt_FromLong(HL->sim_range));
    PyDict_SetItemString(dict, "dist_range",
			 PyInt_FromLong(HL->dist_range));
    PyDict_SetItemString(dict, "kNN",
			 PyInt_FromLong(HL->kNN));

    /* Index performance data */
    PyDict_SetItemString(dict, "bins_visited",
			 PyInt_FromLong(HL->FS_seqs_visited));
    PyDict_SetItemString(dict, "bins_hit",
			 PyInt_FromLong(HL->FS_seqs_hits));
    PyDict_SetItemString(dict, "frags_visited",
			 PyInt_FromLong(HL->seqs_visited));
    PyDict_SetItemString(dict, "frags_hit",
			 PyInt_FromLong(HL->seqs_hits));
    PyDict_SetItemString(dict, "unique_frags_visited",
			 PyInt_FromLong(HL->useqs_visited));
    PyDict_SetItemString(dict, "unique_frags_hit",
			 PyInt_FromLong(HL->useqs_hits));
    PyDict_SetItemString(dict, "search_time",
			 PyFloat_FromDouble(HL->search_time));


    list = PyList_New(HL->actual_seqs_hits);
    for (i=0, hit=HL->hits; i < HL->actual_seqs_hits; i++, hit++) {
      hdict = PyDict_New();
      /* Sequence and defline are not here - saving memory */
#if 0
      obj = PyString_FromStringAndSize(hit->subject.start,
				       hit->subject.len);
      PyDict_SetItemString(hdict, "seq", obj);
#endif
      PyDict_SetItemString(hdict, "seq_id",
			   PyInt_FromLong(hit->sequence_id));
      PyDict_SetItemString(hdict, "seq_from",
			   PyInt_FromLong(hit->sequence_from));
      PyDict_SetItemString(hdict, "seq_to",
			   PyInt_FromLong(hit->sequence_to));
      PyDict_SetItemString(hdict, "dist",
			   PyInt_FromLong(hit->dist));
      PyDict_SetItemString(hdict, "sim",
			   PyInt_FromLong(hit->sim));

      PyList_SET_ITEM(list, i, hdict);
    }
    PyDict_SetItemString(dict, "hits", list);
    Py_DECREF(list);

    return dict;
}
%}

%typemap(out) hit_list * {
  $result = Dict_FromHitList($1);
  /* Free the hit list */
  HIT_LIST_destroy($1);
}


/* Search arguments as list */
%typemap(in, numinputs=1) (SRCH_ARGS *args, int n, hit_list_array *res) (hit_list_array tmp) {
  int i;
  int process_error = 0;
  int list_size;
  PyObject *obj;
  PyObject *list;

  $3 = &tmp;
  if (!PyList_Check($input)) {
    PyErr_SetString(PyExc_ValueError, "Expecting a list");
    return NULL;
  }
  $2 = PyList_GET_SIZE($input);
  $1 = mallocec($2 * sizeof(SRCH_ARGS));
  for (i=0; i < $2; i++) {
    list = PyList_GetItem($input, i);
    if (!PyList_Check(list) || (list_size = PyList_GET_SIZE(list)) < 5) {
      process_error++;
      break;
    }
    /* Default values for optional arguments */
    $1[i].ptype = FS_BINS;
    $1[i].stype = SARRAY;

    /* First item - query sequence */
    obj = PyList_GetItem(list, 0);
    if (!PyString_Check(obj)) {
      process_error++;
      break;
    }
    $1[i].query.start = PyString_AsString(obj);
    $1[i].query.len = PyString_GET_SIZE(obj);

    /* Second item - matrix */
    obj = PyList_GetItem(list, 1);
    if ((SWIG_ConvertPtr(obj, (void **) &($1[i].M),
			 $descriptor(Smatrix *),
			 SWIG_POINTER_EXCEPTION | 0 )) == -1) {
      process_error++;
      break;
    }

    /* Third item - range */
    obj = PyList_GetItem(list, 2);
    if (PyInt_Check(obj)) {
      $1[i].range = PyInt_AS_LONG(obj);
    }
    else if (PyLong_Check(obj)) {
      $1[i].range = PyLong_AsLong(obj);
    }
    else {
      process_error++;
      break;
    }

    /* Fourth item - kNN */
    obj = PyList_GetItem(list, 3);
    if (PyInt_Check(obj)) {
      $1[i].kNN = PyInt_AS_LONG(obj);
    }
    else if (PyLong_Check(obj)) {
      $1[i].kNN = PyLong_AsLong(obj);
    }
    else {
      process_error++;
      break;
    }

    /* Fifth item - query defline */
    obj = PyList_GetItem(list, 4);
    if (!PyString_Check(obj)) {
      process_error++;
      break;
    }
    $1[i].query.id.defline = PyString_AsString(obj);

    /* Optional arguments */

    /* SFUNC */
    if (list_size < 6) continue;
    obj = PyList_GetItem(list, 5);
    if (PyInt_Check(obj)) {
      $1[i].ptype = PyInt_AS_LONG(obj);
    }

    /* PFUNC */
    if (list_size < 7) continue;
    obj = PyList_GetItem(list, 6);
    if (PyInt_Check(obj)) {
      $1[i].stype = PyInt_AS_LONG(obj);
    }
  }
  /* Error - free array */
  if (process_error) {
    free($1);
    PyErr_SetString(PyExc_ValueError, "Wrong Arguments");
    return NULL;
  }
}

%typemap(argout) (SRCH_ARGS *args, int n, hit_list_array *res) {
  int i;

  $result = PyList_New($2);
  for (i=0; i < $2; i++) {
    /* printf("conversion #%d, hits %d\n", i, (*$3)[i].actual_seqs_hits); */
    PyList_SET_ITEM($result, i, Dict_FromHitList(*$3+i));
  }
}

%typemap(freearg) (SRCH_ARGS *args, int n, hit_list_array *res) {
  int i;
  if ($1) free($1);
  if (*$3) {
    for (i=0; i < $2; i++) {
      HIT_LIST_cleanup(*$3+i);
    }
    free(*$3);
  }
}
/* ***** End of hit list conversion code ***** */

extern int next_avail_srch;

typedef struct
{
  %immutable;
  db *s_db;                /* Pointer to the fasta database         */
} Index;

/* Only database attribute is available - everything else can be
   obtained as Python dictionary by calling get_data() method */

typedef enum {FS_BINS, SUFFIX_ARRAY, SEQ_SCAN} SFUNC_TYPE;
typedef enum {SARRAY, DUPS_ONLY, FULL_SCAN} PFUNC_TYPE;

%extend Index {
  /* Call with [] (empty list) as second argument to load from file */
 Index(const char *database, ULINT len,
       const char **sepn, int use_sa = 1,
       int print_flag=0) {
   if (len)
     return FS_INDEX_create(database, len, sepn, use_sa, print_flag);
   else
     return FS_INDEX_load(database);
  }
  ~Index() {
    FS_INDEX_destroy(self);
  }
  void save(const char *filename) {
    FS_INDEX_save(self, filename);
  }
  %newobject __str__;
  char *__str__() {
    return FS_INDEX_sprint_stats(self, 0);
  }
  ULINT seq2bin(char *s1, int l1) {
    return FS_SEQ(self->ptable, s1, l1);
  }
  %newobject print_bin;
  char *print_bin(ULINT bin, int options=1) {
    return FS_INDEX_sprint_bin(self, bin, options);
  }
  %newobject print_stats;
  char *print_stats(int options=3) {
    return FS_INDEX_sprint_stats(self, options);
  }
  ULINT get_bin_size(ULINT bin) {
    return FS_INDEX_get_bin_size(self, bin);
  }
  ULINT get_unique_bin_size(ULINT bin) {
    return FS_INDEX_get_unique_bin_size(self, bin);
  }

  /* Index data as Python dictionary */
  %newobject get_data;
  PyObject *get_data() {
    PyObject *dict = PyDict_New();
    PyObject *obj;
    SEQUENCE_DB *sdb = self->s_db;
    FS_TABLE *ptable = self->ptable;
    int i;

    /* FASTA Sequence Database data */
    PyDict_SetItemString(dict, "db_name",
			 PyString_FromString(sdb->db_name));
    PyDict_SetItemString(dict, "db_length",
			 PyInt_FromLong(sdb->length));
    PyDict_SetItemString(dict, "db_no_seq",
			 PyInt_FromLong(sdb->no_seq));

    /* Index data */
    PyDict_SetItemString(dict, "index_name",
			 PyString_FromString(self->index_name));
    obj = PyList_New(ptable->len);
    for (i=0; i < ptable->len; i++) {
      PyList_SetItem(obj, i, PyString_FromString(ptable->sepn[i]));
    }
    PyDict_SetItemString(dict, "ptable", obj);
    Py_DECREF(obj);

    PyDict_SetItemString(dict, "alphabet",
			 PyString_FromString(ptable->alphabet));
    PyDict_SetItemString(dict, "bins",
			 PyInt_FromLong(self->no_bins));
    PyDict_SetItemString(dict, "fragments",
			 PyInt_FromLong(self->no_seqs));
    PyDict_SetItemString(dict, "unique_fragments",
			 PyInt_FromLong(self->no_useqs));
    PyDict_SetItemString(dict, "indexed_fragment_length",
			 PyInt_FromLong(self->m));
    PyDict_SetItemString(dict, "largest_bin",
			 PyInt_FromLong(self->binL));
    PyDict_SetItemString(dict, "largest_bin_unique",
			 PyInt_FromLong(self->uL));

    return dict;
  }

  hit_list *rng_srch(char *qseq, int qlen, Smatrix *M, int range,
		     int conv_type, SFUNC_TYPE stype = FS_BINS,
		     PFUNC_TYPE ptype = SARRAY, char *qdef ="") {
    BIOSEQ query;
    query.start = qseq;
    query.len = qlen;
    query.id.defline = qdef;
    return FSINDX_rng_srch(self, &query, M, range, conv_type, NULL,
			   stype, ptype);
  }
  hit_list *kNN_srch(char *qseq, int qlen, Smatrix *M, int kNN,
		     SFUNC_TYPE stype = FS_BINS,
		     PFUNC_TYPE ptype = SARRAY, char *qdef ="") {
    BIOSEQ query;
    query.start = qseq;
    query.len = qlen;
    query.id.defline = qdef;
    return FSINDX_kNN_srch(self, &query, M, kNN, NULL, stype, ptype);
  }
  void threaded_search(SRCH_ARGS *args, int n, hit_list_array *res) {
    *res = FSINDX_threaded_srch(self, args, n);
  }
}
