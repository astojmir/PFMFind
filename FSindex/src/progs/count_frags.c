/*
 * Copyright (C) 2005-2006 Victoria University of Wellington
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


#include <stdlib.h>
#include <string.h>
#include "FSindex.h"
#include "misclib.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int main(int argc, char **argv)
{
  const char *index_file;
  int top_len = 20;
  FILE *outstream = stdout;
  FSINDX *I;
  int unique=0;
  int total=0;
  char *s;
  int j,k,m;


  int no_args = 2;

  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file max_length\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  index_file = argv[1];
  top_len = atoi(argv[2]);

  Try {
    I = FS_INDEX_load(index_file);
    if (!I->use_sa) {
      Throw FSexcept(FOPEN_ERR, "Index must have suffix array built.");
    }
    fprintf(outstream, "%8.8s %12.12s %12.12s\n", "Length", "Unique", "Total");
    for (m=1; m<=top_len; m++) {
      unique = 0;
      total = 0;
      for (k=0; k < I->sa_size; k++) {
	/* Check if valid */
	s = fastadb_data_pter(I->s_db, I->sa[k]);
	for (j=0; *s && j < m; j++, s++) {
	  if (*s == 'B' || *s == 'Z' || *s == 'X')
	    break;
	} 

	if (j == m) {
	  total++;
	  /* Check if unique */
	  if (I->lcpx[k] < m) unique++;
	}
      }
      fprintf(outstream, "%8d %12d %12d\n", m, unique, total);
    }

  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
