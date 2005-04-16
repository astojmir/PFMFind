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

  if (argc < 2)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s database \n", argv[0]); 
      exit(EXIT_FAILURE);
    }
  db_name = argv[1];

  Try {
    sdb = fastadb_open(db_name);

    for (i=0; i < sdb->no_seq; i++) {
      printf("%12d %12ld %12d\n", fastadb_data_offset(sdb, sdb->seq[i].start),
	     sdb->seq[i].len, (sdb->seq[i].id.defline - sdb->deflines));
    }
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }

  return 0;
}
