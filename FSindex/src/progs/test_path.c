
#include "misclib.h"



int main(int argc, char **argv)
{
  char *filename;
  char *basename;
  char *dirname;

  if (argc > 1)
    {
      filename = argv[1];
      split_base_dir(filename, &basename, &dirname);
      printf("Filename: %s\n", filename);
      printf("Basename: %s\n", basename);
      printf("Dirname: %s\n", dirname);
      cat_base_dir(&filename, basename, dirname);
      printf("Concatenated: %s\n", filename);
    }
  exit(1);

}
