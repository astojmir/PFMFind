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


#include "misclib.h"
#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include "FSindex.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int FS_PARTITION_VERBOSE = 0;
int SCORE_MATRIX_VERBOSE = 0;
int FS_INDEX_VERBOSE = 0;
int FS_INDEX_PRINT_BAR = 1;
int POS_MATRIX_VERBOSE = 0;
FILE *POS_MATRIX_STREAM;

int main(int argc, char **argv)
{
  const char *filename;
  FSINDX *FSI;

  if (argc < 2)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s index_file \n",
	      argv[0]);  
      exit(EXIT_FAILURE);
    }
  filename = argv[1];
  Try {
    FSI = FS_INDEX_load(filename);
    FS_INDEX_print_stats(FSI, stdout, 2); 
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }
 return EXIT_SUCCESS;
}
