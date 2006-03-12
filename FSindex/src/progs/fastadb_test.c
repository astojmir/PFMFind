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


#include "fastadb.h"
#include <stdlib.h>

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int main(int argc, char **argv)
{
  char *db_name;
  SEQUENCE_DB *sdb;
  int i;
  FILE *fp;

  if (argc < 2)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s database [binfile]\n", argv[0]); 
      exit(EXIT_FAILURE);
    }
  db_name = argv[1];

  Try {
    sdb = fastadb_open(db_name);

    for (i=0; i < sdb->no_seq; i++) {
      printf("%12d %12ld %12d\n", fastadb_data_offset(sdb, sdb->seq[i].start),
	     sdb->seq[i].len, (sdb->seq[i].id.defline - sdb->deflines));

    }
    if (argc >= 3) {
      fp = fopen(argv[2], "wb");
      if (fp == NULL)
	Throw FSexcept(FOPEN_ERR, "Could not open the file %s.", argv[2]);

      fwrite(&sdb->deflines_len, sizeof(ULINT), 1, fp);
      fwrite(&sdb->seq_data_len, sizeof(ULINT), 1, fp);
      fwrite(sdb->deflines, 1, sdb->deflines_len, fp);
      fwrite(sdb->seq_data, 1, sdb->seq_data_len, fp);
      fclose(fp);
    }
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }

  return 0;
}
