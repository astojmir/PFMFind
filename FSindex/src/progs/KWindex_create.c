#include "misclib.h"
#include "keyword.h"
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>


static 
void print_help(const char *progname)
{
  fprintf(stderr, "Usage: %s basename\n",progname);
}


int main(int argc, char **argv)
{
  char *bin_name;
  const char *basename;
  KW_INDEX *KWI;

  if (argc < 2) 
    {
      print_help(argv[0]);
      return EXIT_FAILURE;
    }
  basename = argv[1];

  bin_name = mallocec(strlen(basename)+5);
  strcpy(bin_name, basename);  
  strcat(bin_name, ".kbn");
  
  KWI = KW_INDEX_load_txt(basename);
  KW_INDEX_save_bin(KWI, bin_name);

  return EXIT_SUCCESS;
}
