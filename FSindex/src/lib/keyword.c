#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "misclib.h"
#include "keyword.h"
#include "hit_list.h"



/********************************************************************/    
/*                                                                  */
/*             OC (orthologous clusters) module                     */ 
/*                                                                  */
/********************************************************************/    

OC *load_oc(const char *filename)
{
  /* two files:
     1. maps sequence to cluster
     2. maps cluster to its frequency
  
     Both are text files, numbers separated by '\n',
     the first number is the size of the vector.
     
     Filenames obtained by tagging .scl (1.) or .cfq (2.)
     to filename. 

  */
#define BUF_SIZE 512
  char buffer[BUF_SIZE+1];
  int i;

  char *scl_name = mallocec(strlen(filename)+5);
  char *cfq_name = mallocec(strlen(filename)+5);
  FILE *scl;
  FILE *cfq;
  OC *oc;

  strcpy(scl_name, filename);  
  strcat(scl_name, ".scl");
  strcpy(cfq_name, filename);  
  strcat(cfq_name, ".cfq");


  if((scl = fopen(scl_name, "r")) == NULL)
    {
      fprintf(stderr, 
	      "Could not open cluster file %s!\n", scl_name);
      exit(EXIT_FAILURE);
    }

  fgets(buffer, BUF_SIZE, scl);
  sscanf(buffer, "%d", &oc->no_seq);
  oc->oclstr = mallocec(oc->no_seq * sizeof(int));
  for (i=0; i < oc->no_seq; i++)
    {
      fgets(buffer, BUF_SIZE, scl);
      sscanf(buffer, "%d", oc->oclstr+i);
    }
  fclose(scl);

  if((cfq = fopen(cfq_name, "r")) == NULL)
    {
      fprintf(stderr, 
	      "Could not open cluster file %s!\n", cfq_name);
      exit(EXIT_FAILURE);
    }

  fgets(buffer, BUF_SIZE, cfq);
  sscanf(buffer, "%d", &oc->no_clstrs);
  oc->clfreq = mallocec(oc->no_clstrs * sizeof(int));
  for (i=0; i < oc->no_clstrs; i++)
    {
      fgets(buffer, BUF_SIZE, cfq);
      sscanf(buffer, "%d", oc->clfreq+i);
    }
  fclose(cfq);

  free(scl_name);
  free(cfq_name);

  return oc;
}


void clean_oc(OC *oc)
{
  free(oc->oclstr);
  free(oc->clfreq);
}

void get_cons_ratio(OC *oc, HIT_LIST_t *HL)
{
  int i;                   
  int j;                   
  int c;                   /* count                                 */
  int s;                   /* start of cluster                      */
  double r;                /* conservation ratio                    */

  /* Assign clusters to all sequencies */
  for (i=0; i < HL->actual_seqs_hits; i++)
    HL->hits[i].oc_cluster = oc->oclstr[HL->hits[i].sequence_id];
  HIT_LIST_sort_oc(HL);

  s=0;
  c=0;
  for (i=1; i < HL->actual_seqs_hits; i++)
    {
      if (HL->hits[i].oc_cluster == -1)
	continue;
      if (HL->hits[i].sequence_id != HL->hits[i-1].sequence_id)
	{
	  if (HL->hits[i].oc_cluster == HL->hits[i-1].oc_cluster)
	    c++;
	  else
	    {
	      r = (double) c / oc->clfreq[HL->hits[i-1].oc_cluster];
	      for (j=s; j < i; j++)
		HL->hits[j].cratio = r;
	      s=i;
	      c=0;
	    }
	}
    }
  r = (double) c / oc->clfreq[HL->hits[i-1].oc_cluster];
  for (j=s; j < i; j++)
    HL->hits[j].cratio = r;
}

