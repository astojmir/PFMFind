#ifndef _KEYWORD_H
#define _KEYWORD_H

#ifdef USE_MPATROL
#include <mpatrol.h>
#endif

/********************************************************************/    
/********************************************************************/    
/***                                                              ***/
/***               PROTOTYPES                                     ***/ 
/***                                                              ***/
/********************************************************************/    
/********************************************************************/    

typedef struct 
{
  int ns;                  /* Number of sequences                   */
  int nsc;                 /* Number of sequences with clusters     */
  int nsk;                 /* Number of sequences with keywords     */
  int nc;                  /* Number of clusters                    */
  int nk;                  /* Number of keywords                    */
  int *clf;                /* Cluster frequency                     */
  int *kwf;                /* Keyword frequency                     */
  ULINT *cld;              /* Cluster descriptor (string) offset    */
  ULINT *kwd;              /* Keyword descriptor (string) offset    */
  int *seq2cl;             /* Sequence to cluster map               */
  int *seq2kw0;            /* Sequence to keywords offset           */
  int *seq2kw1;            /* Sequence to keywords length           */
  ULINT cld_hs;            /* Cluster descriptor heap size          */
  ULINT kwd_hs;            /* Keyword descriptor heap size          */
  ULINT skw_hs;            /* Sequence to keywords list heap size   */
  char *kwd_h;             /* Keyword descriptor heap               */
  char *cld_h;             /* Cluster descriptor heap               */
  int *skw_h;              /* Sequence to keywords list heap        */
} KW_INDEX;


KW_INDEX *KW_INDEX_load_txt(const char *basename);
void KW_INDEX_save_bin(KW_INDEX *KWI, const char *filename);
KW_INDEX *KW_INDEX_load_bin(const char *filename);
void KW_INDEX_clear(KW_INDEX *KWI);


#endif /* #ifndef  _KEYWORD_H */
