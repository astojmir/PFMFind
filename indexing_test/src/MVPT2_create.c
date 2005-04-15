#include <stdlib.h>
#include <string.h>
#include "mvptree.h"


EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int main(int argc, char **argv)
{
  char *db_name;
  ULINT frag_len;
  char *matrix_name;
  int num_ref;
  char *save_name;
  MVP_TREE *MVPT2;

  int no_args = 5;
  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s db_name matrix_name frag_len num_ref save_name\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  db_name = argv[1];
  matrix_name = argv[2];
  frag_len = atoi(argv[3]);
  num_ref = atoi(argv[4]);
  save_name = argv[5];

  Try {
    MVPT2 = MVP_TREE_init2(db_name, frag_len, matrix_name, num_ref);
    MVP_TREE_save2(MVPT2, save_name);    
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

