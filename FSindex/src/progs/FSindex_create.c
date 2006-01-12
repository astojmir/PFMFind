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


#include "fastadb.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "FSindex.h"
#include "misclib.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

#define BUF_SIZE 1024
static 
FSINDX *FS_INDEX_create_from_file(FILE *fp)
{
  FSINDX *FSI = NULL;
  char buf[BUF_SIZE];
  char *database = NULL;
  char *outfile = NULL;
  ULINT len = 0;
  int use_sa = 1;
  const char **sepn = NULL;
  int row = 0;
  int j = 0;

  while(!feof(fp)) {
    fgets(buf, BUF_SIZE, fp);
    if (row == 0) {
      database = strdup(buf);
      database[strlen(database)-1]='\0';
    }
    else if (row == 1) {
      outfile = strdup(buf);
      outfile[strlen(outfile)-1]='\0';
    }
    else if (row == 2) {
      sscanf(buf, "%ld %d", &len, &use_sa);
      if (len == 0) break;
      sepn = callocec(len, sizeof(char *));
    }
    else {
      if (row > len+2) break;
      buf[strlen(buf)-1] = '\0';
      sepn[row-3] = strdup(buf);
      j++;
    }
    row++;
  }

  Try {
    if (database && outfile && len && j == len) {
      printf("Creating index...\n");
      FSI = FS_INDEX_create(database, len, sepn, use_sa, 1);
      printf("Saving index...\n");
      FS_INDEX_save(FSI, outfile);
      for (j=0; j < len; j++)
	free((char *)sepn[j]);
      free(sepn);
      free(database);
      free(outfile);
    }
    else {
      Throw FSexcept(BAD_ARGS, "Insufficient arguments supplied.");
    }
  }
  Catch(except) {
    for (j=0; j < len; j++)
      free((char *)sepn[j]);
    free(sepn);
    free(database);
    free(outfile);
    Throw except;
  }
  return FSI;
}

int main(int argc, char **argv)
{
  FSINDX *FSI = NULL;
  FILE *fp = stdin;

  if (argc > 1 && (fp=fopen(argv[1],"r"))==NULL) {
    fprintf(stderr, "Could not open %s\n", argv[1]);
    return EXIT_FAILURE;
  } 
  Try {
    FSI = FS_INDEX_create_from_file(fp);
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }
  if (argc > 1) fclose(fp);
  FS_INDEX_destroy(FSI);
#ifdef USE_DMALLOC
  dmalloc_shutdown();
#endif
  return EXIT_SUCCESS;
}
