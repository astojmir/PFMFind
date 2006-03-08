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
  char *fasta_db;

  ULINT letter_count[256];
  ULINT i=0;
  ULINT sum=0;
  ULINT seqs=0;
  int def_flag = 0;
  ULINT alphabet=0;
  FILE *stream;
  int c;


  if (argc > 1)
    fasta_db = argv[1];
  else 
    exit(EXIT_FAILURE);

  if ((stream = fopen(fasta_db, "r")) == NULL)
    {
      fprintf(stderr, "Unable to open %s\n ", fasta_db);
      exit(EXIT_FAILURE);
    }  
  memset(&letter_count[0], 0, 256*sizeof(ULINT));

  while ((c=fgetc(stream)) != EOF)
    {
      if (c=='>' && def_flag == 0)
	{
	  def_flag = 1;
	  seqs++;
	  continue;
	}
      else if (c == '\n')
	{
	  def_flag = 0;
	  continue;
	}
      if (!def_flag)
	{
	  letter_count[c]++;
	  sum++;
	}
    }

  for (i=0, alphabet=0; i < 256; i++)
    if (letter_count[i] > 0)
      alphabet++;

  printf("******* count_letters Results *******\n");
  printf("Database: %s\n", fasta_db);
  printf("Total sequences: %ld\n", seqs);
  printf("Total letters: %ld\n", sum);
  printf("Alphabet size: %ld\n", alphabet);

  printf("\n%10s %10s %10s %10s\n", "ASCII code", "Letter",
	 "Frequency", "Percentage");

  for (i=0; i < 256; i++)
    {
      if (letter_count[i] > 0)
	{
	  if (i <= 32) 
	    c = 32;
	  else 
	    c = i;
	  printf("%5ld%5s %5c%5s %10ld %10.4f\n", i, "", c, "",
		 letter_count[i], 
		 (((float)letter_count[i])/((float) sum))*100); 
	}
    }
  return EXIT_SUCCESS;
}
