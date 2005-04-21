#include "SAindex.h"
#include "misclib.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int main(int argc, char **argv)
{
  const char *db_name;
  const char *save_name;
  SAINDX *SAI = NULL;

  int no_args = 2;
  if (argc < no_args+1)
    {
      fprintf(stderr,"Insufficient arguments \n");
      fprintf(stderr,"Usage: %s db_name save_name\n", argv[0]);  
      exit(EXIT_FAILURE);
    }
  db_name = argv[1];
  save_name = argv[2];

  Try {
    SAI = SA_INDEX_create(db_name, 1);
    SA_INDEX_save(SAI, save_name);    
  }
  Catch(except) {
    fprintf(stderr, "Error %d: %s\n", except->code, except->msg);  
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
