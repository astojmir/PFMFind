#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "fastadb.h"
#include "partition.h"
#include "misclib.h"

/* Allocate memory in 8 Mb blocks */
#ifndef ALLOC_BLOCK_SIZE
#define ALLOC_BLOCK_SIZE (1 << 23)
#endif

/****************************************************************
 ****************************************************************/

/* static inline ULINT find_residue_seq(ULINT residue,
           const ULINT *boundaries, const ULINT length)

   This little routine assumes that boundaries is a pointer to an
   array of ULINTs which is sorted in ascending order and whose
   elements represent the starting points of sequences put one after
   the other. The goal is to find the index of the sequence the
   residue belongs to.  

   The function can return unexpected results if not used as
   intended (i.e. nothing is checked). 

   Parameters: 

   ULINT residue - number of sequence to find index of

   const ULINT* boundaries - the boundaries of sequences.

   const ULINT length - length of boundaries array.

*/

static inline
ULINT find_residue_seq(ULINT residue, const ULINT *boundaries,
		       const ULINT length)
{
  ULINT a = 0;
  ULINT b = length-1;
  ULINT k = b/2; 

  while (a<=b)
    {
      k=(a+b)/2;
      if (boundaries[k] <= residue)
        {
          if (boundaries[k+1] > residue)
	    return k;
	  else
	    a=k+1;
	}
      else
	b=k;
    }
  return k;
}


static inline
int fastadb_put_line(char *buffer, char **heap, ULINT *len, 
		     ULINT *max_len)
{
  int l;
  l = strlen(buffer);

  if (*len + l + 2 >= *max_len)
    {
      *max_len += ALLOC_BLOCK_SIZE;
      *heap = reallocec(*heap, *max_len); 
    }
  if (buffer[l-1] == '\n')
    l--;
  
  if (l > 0)
    memcpy(*heap+*len, buffer, l);

  *len += l;
  return l;  
}

static inline
void fastadb_update_deflines(SEQUENCE_DB *s_db, char **oldheap)
{
  ULINT i;
  ULINT from = 0;
  ULINT to = s_db->current_seq + 1;

  if (*oldheap == s_db->deflines)
    return;
  if (s_db->keep_oldseqs == NO)
    from = s_db->current_seq;
      
  for (i=from; i < to; i++)
    s_db->seq[i].id.defline += (s_db->deflines - *oldheap);
  *oldheap = s_db->deflines; 
}

static inline
void fastadb_update_seq_data(SEQUENCE_DB *s_db, char **oldheap)
{
  ULINT i;
  ULINT from = 0;
  ULINT to = s_db->current_seq + 1;

  if (*oldheap == s_db->seq_data)
    return;
  if (s_db->keep_oldseqs == NO)
    from = s_db->current_seq;
      
  for (i=from; i < to; i++)
    s_db->seq[i].start += (s_db->seq_data - *oldheap);
  *oldheap = s_db->seq_data;
}

static inline
void fastadb_put_char(char c, char **heap, ULINT *len)
{
  (*heap)[(*len)++] = c;
}


/* static int fastadb_check_array_size(SEQUENCE_DB *s_db)

   Checks if there is allocated space for a new sequence and if there
   isn't it allocates some. This is based on the value in
   s_db->no_seq rather than s_db->current_seq.

   Returns 1 if the memory is succesfully allocated (otherwise the
   program is exited).

   Parameters: 

   SEQUENCE_DB *s_db - Sequence database.

*/

static inline
int fastadb_check_array_size(SEQUENCE_DB *s_db)
{
  if (s_db->no_seq >= s_db->max_seqs)
    {
      s_db->max_seqs += MAX_SEQUENCES;
      s_db->start_pts = (ULINT *) reallocec(s_db->start_pts,
				  (s_db->max_seqs+1) * sizeof(ULINT)); 
      s_db->seq = (BIOSEQ *) reallocec(s_db->seq, 
			    s_db->max_seqs * sizeof(BIOSEQ));
      s_db->offset = (long * ) reallocec(s_db->offset, 
			       s_db->max_seqs * sizeof(long));

    }
  return 1;
}

/* static int fastadb_reduce_array_size(SEQUENCE_DB *s_db)

   Trims the allocated space for sequences.

   Returns 1 if the memory is succesfully reallocated 
   (otherwise the program is exited).

   Parameters: 

   SEQUENCE_DB *s_db - Sequence database.

*/

static inline
int fastadb_reduce_array_size(SEQUENCE_DB *s_db)
{
  s_db->start_pts = reallocec(s_db->start_pts, (s_db->no_seq + 1)
			      * sizeof(ULINT)); 
  s_db->seq = reallocec(s_db->seq, s_db->no_seq * sizeof(BIOSEQ));
  s_db->offset = reallocec(s_db->offset, s_db->no_seq * sizeof(long));
  return 1;
}


/* static int fastadb_load_next_seq(SEQUENCE_DB *s_db)

   Loads the next sequence to memory. The pointers to this sequence
   are set in s_db->seq[s_db->current_seq] and then s_db->current_seq
   is incremented by one. The previous sequence is either kept or
   deleted depending on s_db->keep_oldseqs (TO DO). The sequence 
   definiton line (the one starting with '>') can be loaded depending
   on s_db->retrieve_deflines.

   Returns 0 if there is no next sequence (EOF reached) or 1 if the
   sequence was succesfully loaded.

   Parameters: 

   SEQUENCE_DB *s_db - Sequence database.

*/
#define BUF_SIZE 512
static inline
int fastadb_load_next_seq(SEQUENCE_DB *s_db)
{
  int c;
  int l;
  static char buffer[BUF_SIZE+1];
  char *old_deflines;
  char *old_seq_data;

  /* Just in case the line is bigger than buffer, this would give
     BUF_SIZE as strlen */
  buffer[BUF_SIZE] = '\0';

  if (feof(s_db->dbfile))
    return 0;

  if (s_db->keep_oldseqs == NO)
    {
      s_db->deflines_len = 0;
      s_db->seq_data_len = 0;      
      s_db->length = 0;
    }

  /* The heap rows as well as the s_db->seq must be allocated before
     entering here */
  fastadb_check_array_size(s_db);

  /* Process the comment line */
  if (s_db->retrieve_deflines)
    {
      /* We know that the first character is '>' because this is a new
	 sequence. We also know the new sequence exists (we exit at
	 the beginning if EOF */
      old_deflines = s_db->deflines;
      s_db->seq[s_db->current_seq].id.defline 
	= s_db->deflines + s_db->deflines_len;

      while (1)
	{
	  fgets(buffer, BUF_SIZE, s_db->dbfile);
	  l = fastadb_put_line(buffer+1, &(s_db->deflines), 
			   &(s_db->deflines_len),
			       &(s_db->deflines_max_len));
	  fastadb_update_deflines(s_db, &old_deflines);
	  if (buffer[l+1] == '\n')	  
	    break;
	}
      fastadb_put_char('\0', &(s_db->deflines),
		       &(s_db->deflines_len)); 
    }
  else
    {
      while (1)
	{
	  fgets(buffer, BUF_SIZE, s_db->dbfile);
	  if (buffer[strlen(buffer)-1] == '\n')
	    break;
	}
    }

  /* Now read and store the sequence itself */
  old_seq_data = s_db->seq_data;
  s_db->seq[s_db->current_seq].start = 
    s_db->seq_data + s_db->seq_data_len; 
  s_db->seq[s_db->current_seq].len = 0;
  while (1)
    {
      c = fgetc(s_db->dbfile);
      if (c == EOF)
	break;

      if (c == '>')
	{ 
	  ungetc(c, s_db->dbfile);
	  break;
	}

      if (c != '\n')
	{
	  buffer[0] = c;
	  fgets(buffer+1, BUF_SIZE, s_db->dbfile);
	  l = fastadb_put_line(buffer, &(s_db->seq_data), 
			   &(s_db->seq_data_len),
			   &(s_db->seq_data_max_len)); 
	  fastadb_update_seq_data(s_db, &old_seq_data);
	  s_db->seq[s_db->current_seq].len += l;
	}
    }
  fastadb_put_char('\0', &(s_db->seq_data),
		   &(s_db->seq_data_len)); 
  s_db->current_seq++;
  return 1;
}

/* static int fastadb_count_next_seq(SEQUENCE_DB *s_db)

   Counts the length of the next sequence without loading to memory. 
   s_db->current_seq is incremented by one. 

   Returns 0 if there is no next sequence (EOF reached) or 1 if the
   sequence was succesfully loaded.

   Parameters: 

   SEQUENCE_DB *s_db - Sequence database.

*/

static inline
int fastadb_count_next_seq(SEQUENCE_DB *s_db)
{
  int c;
  int c_old = 0;

  if (feof(s_db->dbfile))
    return 0;

  /* s_db->seq must be allocated before
     entering here */
  fastadb_check_array_size(s_db);

  /* Process the comment line */
  while (1)
    if ((c = fgetc(s_db->dbfile))=='\n') 
      break;

  /* Now read and count the sequence itself */
  s_db->seq[s_db->current_seq].len = 0;
  while (1)
    {
      c_old = c;
      c = fgetc(s_db->dbfile);
      if (c==EOF)
	break;
      if (c =='>' && c_old == '\n')
	{ 
	  ungetc(c, s_db->dbfile);
	  break;
	}
      if (c != '\n')
	s_db->seq[s_db->current_seq].len++;
    }
  s_db->current_seq++;
  return 1;
}

/* SEQUENCE_DB *fastadb_open(const char *db_name, 
                             fastadb_arg *argt, 
			     fastadb_argv_t *argv)

   Returns a pointer to an instance of SEQUENCE_DB to be used
   with other fastadb_ routines. It must be closed when finished
   with fastadb_close. If there is a problem, a NULL pointer is
   returned. 

   Parameters: 

   const char *db_name - The name (relative to the working directory)
                         of the FASTA format database to load.
   
   fastadb_arg *argt - A 'NONE'-terminated list of fastadb_arg
                       arguments. The first entry can be NONE in which
		       case there are no arguments and the defaults
		       are used. See below.

   fastadb_argv_t *argv - A list of values of the arguments whose type
                          is specified in argt. For each argt term
			  except the last one there must be an argv
			  term of the type specified by the
			  corresponding type in argt. 

   The valid arguments, their types and their legal values are below.
   If specific values are not listed then all legal values for the
   type specified are allowed.

   The type fastadb_argv_t is a union having the fields:
     ULINT max_row_length - Maximum length of a heap row 
     access_type_t access_type - Access type 
     yn_bool_t keep_oldseqs - Keep old sequences ?
     yn_bool_t retrieve_deflines - Retrieve deflines  

   Keep old sequences option is not implemented at the moment: all
   sequneces are stored as they are read and the space is only
   released by fastadb_close()

   The values in argt[i] determine the type store in argv[i].
   They can be:

   NONE - the last term in the list. The first term can be NONE.

   ACCESS_TYPE - Associated with access_type field. The valid values
                 are:
		 * MEMORY - store all sequences and comment lines in
		            memory;
		 * RNDM - retrieve each sequence from the file as
		            needed.
		 * SEQUENTIAL - retrieve sequences one by one. Use
		                single pass only. 
			    
     Default: MEMORY.

   KEEP_OLDSEQS - Associated with keep_oldseqs field. Only relevant if
                  access type is SEQUENTIAL. The values are:
		 * NO - Overwrite the previous sequence in memory as
	 	   the next one is loaded.
		 * YES - Keep all sequences in memory until the
		         database is closed.
     Default: YES. 

   RETRIEVE_DEFLINES - Retrieve the comment lines from the FASTA file.
                       Associated with retrieve_deflines field.
		      * NO - Do not load comment lines. Use 
		        id.id_num field of BIOSEQ to give the number
			of the sequence in the database.
		      * YES - Load the comment line and return a
			pointer to it (as a '\0' terminated string) in
			the id.defline field of BIOSEQ. 
     Default: NO.

   MAX_ROW_LENGTH - Length of a segment of memory first
                    allocated for storing data. One such segment (row)
		    is allocated for residues and one for
		    deflines.  This is the default length of 
		    seq_data segment. The length of defline segment is
		    1/4 of the seq_data segment. They are grown if
		    needed. Associated with max_row_length field. 

     Default: 2^26-1 = 67108863 bytes = 64 MB.

*/
  
SEQUENCE_DB *fastadb_open(const char *db_name, fastadb_arg *argt, 
			  fastadb_argv_t *argv)
{
  EXCEPTION *except;
  FILE *stream;
  SEQUENCE_DB *s_db;  
  int real_file = 1;
  char *old_deflines;
  char *old_seq_data;

  if (db_name == NULL)
    {
      stream = stdin;
      db_name = "stdin";
      real_file = 0;
    }
  else if ((stream = fopen(db_name, "r")) == NULL)
    Throw FSexcept(FOPEN_ERR, 
		   "fastadb_open(): Could not open fasta database %s",
		   db_name);

  /* Allocate and set defaults */
  s_db = callocec(1, sizeof(SEQUENCE_DB));

  Try {

    s_db->db_name = strdup(db_name);
    s_db->length = 0;
    s_db->no_seq = 0;

    s_db->max_seqs = MAX_SEQUENCES;
    s_db->start_pts = mallocec((s_db->max_seqs+1) * sizeof(ULINT));
    s_db->seq = mallocec(s_db->max_seqs * sizeof(BIOSEQ));
    s_db->offset = mallocec(s_db->max_seqs * sizeof(long));

    s_db->seq_data = NULL;
    s_db->seq_data_len = 0;
    s_db->seq_data_max_len = ALLOC_BLOCK_SIZE;

    s_db->deflines = NULL;
    s_db->deflines_len = 0;
    s_db->deflines_max_len = ALLOC_BLOCK_SIZE;

    s_db->dbfile = stream;
    s_db->real_file = real_file;
    s_db->current_seq=0;
    s_db->access_type = MEMORY;
    s_db->keep_oldseqs = YES;
    s_db->retrieve_deflines = NO;

    s_db->no_frags = 0;
    s_db->min_len = 0;
    s_db->max_len = 0;
    s_db->L_start_pts = NULL;
    s_db->f_start_pts = NULL;
    s_db->NL = 0;
    s_db->current_frag = 0;

    /* Check parameter list and change params */
    while (*argt != NONE)
      {
	switch (*argt)
	  {
	  case (ACCESS_TYPE):
	    s_db->access_type = argv->access_type;
	    break;
	  case (KEEP_OLDSEQS):
	    s_db->keep_oldseqs = argv->keep_oldseqs;
	    break;
	  case (RETREIVE_DEFLINES):
	    s_db->retrieve_deflines = argv->retrieve_deflines;
	    break;
	  case (MAX_ROW_LENGTH):
	    s_db->seq_data_max_len = argv->max_row_length;
	    break;
	  default:
	    Throw FSexcept(BAD_ARGS, 
			   "fastadb_open(): Unreckognised option!");
	  }
	argt++;
	argv++;
      }
    if (s_db->access_type == MEMORY)
      s_db->keep_oldseqs = YES;
    else if (s_db->access_type == RNDM)
      s_db->keep_oldseqs = NO;
    
    /* Allocate the storage heap */
    
    s_db->seq_data = mallocec(s_db->seq_data_max_len);

    if (s_db->retrieve_deflines)
      {
	s_db->deflines = mallocec(s_db->deflines_max_len);
      }

    /* Either load or just pass sequences in the database */
    if (s_db->access_type == MEMORY)
      {
	while(1)
	  {
	    if(!fastadb_load_next_seq(s_db))
	      break;
	    s_db->start_pts[s_db->no_seq] = s_db->length;
	    s_db->length += s_db->seq[s_db->current_seq-1].len;
	    s_db->no_seq++;
	  }
	if (s_db->real_file)
	  fclose(s_db->dbfile);
	s_db->dbfile = NULL;
	s_db->start_pts[s_db->no_seq] = s_db->length;
      }
    else if (s_db->access_type == RNDM)
      {
	while(1)
	  {
	    s_db->offset[s_db->current_seq]=ftell(s_db->dbfile);
	    if(!fastadb_count_next_seq(s_db))
	      break;
	    s_db->start_pts[s_db->no_seq] = s_db->length;
	    s_db->length += s_db->seq[s_db->current_seq-1].len;
	    s_db->no_seq++;
	  }
	rewind(s_db->dbfile);
	s_db->start_pts[s_db->no_seq] = s_db->length;
	s_db->current_seq = 0;
      }
    
    if (s_db->access_type != SEQUENTIAL)
      {
	/* Trim memory segments - we add 30 chars more to seq_data
	   so that we can have fast fragment comparison near the end
	*/
	fastadb_reduce_array_size(s_db);
	old_deflines = s_db->deflines;
	old_seq_data = s_db->seq_data;
	s_db->seq_data = reallocec(s_db->seq_data, s_db->seq_data_len+30);
	s_db->seq_data_max_len = s_db->seq_data_len+30;
	s_db->deflines = reallocec(s_db->deflines, s_db->deflines_len);
	s_db->deflines_max_len = s_db->deflines_len;
	fastadb_update_deflines(s_db, &old_deflines);
	fastadb_update_seq_data(s_db, &old_seq_data);
      }
  }
  Catch (except) {
    fastadb_close(s_db);
    Throw except;
  }
  return s_db;
}


/* int fastadb_close(SEQUENCE_DB *s_db)

   Cleans up the database, frees all memory and close the file.
   Returns zero if something failed and 1 on success.

   Parameters: 

   SEQUENCE_DB *s_db - Database to be closed.

*/
int fastadb_close(SEQUENCE_DB *s_db)
{
  if (s_db == NULL) return 1;
  free((char *)s_db->db_name);
  /* Free heap storage */
  free(s_db->seq_data);
  free(s_db->deflines);
  s_db->seq_data_max_len = 0;
  s_db->seq_data_len = 0;
  s_db->deflines_max_len = 0;
  s_db->deflines_len = 0;

  /* Free sequence data */
  free(s_db->seq);
  s_db->seq = NULL;
  free(s_db->offset);
  s_db->offset = NULL;
  free(s_db->start_pts);
  s_db->start_pts = NULL;

  /* Close the file if needed */
  if (s_db->dbfile != NULL && s_db->real_file)
    fclose(s_db->dbfile);

  /* Free s_db */
  free(s_db);

  return 1;
}

/* int fastadb_get_length(SEQUENCE_DB *s_db, ULINT *length)

   Gives the number of letters stored in the database.
   Returns 1 if successful, otherwise 0. 

   Parameters: 

   SEQUENCE_DB *s_db - Sequence database.

   ULINT *no_seqs    - Pointer to the space to store the length of
                       database.

*/
int fastadb_get_length(SEQUENCE_DB *s_db, ULINT *length)
{
  *length = s_db->length;
  if (s_db->length > 0)
    return 1;
  else
    return 0;
}


/* int fastadb_get_noseqs(SEQUENCE_DB *s_db, ULINT *no_seqs)

   Gives the number of sequences stored in the database.
   Returns 1 if successful, otherwise 0. 

   Parameters: 

   SEQUENCE_DB *s_db - Sequence database.

   ULINT *no_seqs    - Pointer to the space to store the number of
                       sequences.

 
*/
int fastadb_get_noseqs(SEQUENCE_DB *s_db, ULINT *no_seqs)
{
  *no_seqs = s_db->no_seq;
  if ( s_db->no_seq > 0)
    return 1;
  else
    return 0;
}


/* int fastadb_get_seq(SEQUENCE_DB *s_db, ULINT seq_no, BIOSEQ **seq)

   Gets a sequence from the database. Returns 0 if the sequence does
   not exist and 1 if succesfully obtained.

   Parameters: 

   SEQUENCE_DB *s_db - Sequence database.

   ULINT seq_no      - The sequence number to retrieve. It ranges from
                       0 to the number of sequences retrieved by
		       fastadb_get_noseqs().  
  
   BIOSEQ *seq       - Pointer to the space where the information
                       about the sequence is stored. This space must be
		       allocated before entering this function.
                       NOTE: The field seq->start points to the actual
		       sequence and should not be changed. If another
		       instance is necessary then it must be copied.

*/
int fastadb_get_seq(SEQUENCE_DB *s_db, ULINT seq_no, BIOSEQ **seq)
{
  if (!(seq_no < s_db->no_seq))
    Throw FSexcept(INDEX_OUT_OF_RANGE, "Index out of range.");

  if (s_db->access_type == SEQUENTIAL)
    {
      if (s_db->current_seq != seq_no)
	Throw FSexcept(INCONSITENT_DATA, 
		       "Requested fragment is not current.");
      else return fastadb_get_next_seq(s_db, seq);
    }

  s_db->current_seq = seq_no;
  if (s_db->access_type == MEMORY)
    {
      s_db->current_seq++;
    }
  else if (s_db->access_type == RNDM)
    {
      fseek(s_db->dbfile, s_db->offset[seq_no], SEEK_SET); 
      fastadb_load_next_seq(s_db);
    }

  *seq = s_db->seq+seq_no;
  return 1;
}
 
/* int fastadb_get_next_seq(SEQUENCE_DB *s_db, BIOSEQ *seq)

   Gets the next sequence from the database. The counter can be reset by
   using fastadb_get_seq(). Returns 0 if the sequence does not exist
   and 1 if succesfully obtained. 

   Parameters: 

   SEQUENCE_DB *s_db - Sequence database.
  
   BIOSEQ *seq       - Pointer to the space where the information
                       about the sequence is stored. This space must be
		       allocated before entering this function.
                       NOTE: The field seq->start points to the actual
		       sequence and should not be changed. If another
		       instance is necessary then it must be copied.

*/

int fastadb_get_next_seq(SEQUENCE_DB *s_db, BIOSEQ **seq)
{
  ULINT i = s_db->current_seq;

  if (s_db->access_type == SEQUENTIAL)
    {
      if(!fastadb_load_next_seq(s_db))
	Throw FSexcept(EOF_REACHED, "EOF reached.");

      s_db->start_pts[s_db->no_seq] = s_db->length;
      s_db->length += s_db->seq[s_db->current_seq-1].len;
      s_db->no_seq++;
    }
  else
    {
      if (!(s_db->current_seq < s_db->no_seq))
	Throw FSexcept(INDEX_OUT_OF_RANGE, "Index out of range.");

      if (s_db->access_type == MEMORY)
	s_db->current_seq++;
      else if (s_db->access_type == RNDM)
	fastadb_load_next_seq(s_db);
    }
  
  *seq = s_db->seq+i;
  return 1;
}

