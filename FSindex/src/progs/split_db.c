/********************************************************************

  count_letters.c

  counts letters in molecular sequence database

 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "misclib.h"



int main(int argc, char **argv)
{
  const char *fasta_db;
  char *out_filename;
  int def_flag = 0;
  FILE *stream;
  FILE *outstream;
  int c;
  int max_seqs;
  int i=0;
  int j=1;


  if (argc > 2)
    {
      fasta_db = argv[1];
      max_seqs = atoi(argv[2]);
    }
  else
    {
      fprintf(stderr, "Use %s fasta_db max_seqs\n", argv[0]);
      exit(EXIT_FAILURE);
    }

  if ((stream = fopen(fasta_db, "r")) == NULL)
    {
      fprintf(stderr, "Unable to open %s\n ", fasta_db);
      exit(EXIT_FAILURE);
    }
  
  out_filename = mallocec(strlen(fasta_db)+20);
  sprintf(out_filename, "%s_p%05d", fasta_db, j);
  if ((outstream = fopen(out_filename, "w")) == NULL)
    {
      fprintf(stderr, "Unable to open %s\n ", out_filename);
      exit(EXIT_FAILURE);
    }  
 
  while ((c=fgetc(stream)) != EOF)
    {
      if (c=='>' && def_flag == 0)
	{
	  i++;
	  def_flag = 1;
	  if (i > max_seqs)
	    {
	      ungetc(c, stream);
	      fclose(outstream);	      
	      j++;
	      i=1;
	      sprintf(out_filename, "%s_p%05d", fasta_db, j);
	      if ((outstream = fopen(out_filename, "w")) == NULL)
		{
		  fprintf(stderr, "Unable to open %s\n ", out_filename);
		  exit(EXIT_FAILURE);
		}
	      continue;
	    }
	}
      else if (c == '\n')
	{
	  def_flag = 0;
	}
      fputc(c, outstream);      
    }
  exit(EXIT_SUCCESS);

}
