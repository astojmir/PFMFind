#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "smatrix.h"
#include "misclib.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int main(int argc, char **argv)
{

  /* Not too thorough but should cover all */
  SCORE_MATRIX *SS;
  SCORE_MATRIX *SD;
  SCORE_MATRIX *PS;
  SCORE_MATRIX *PD;
  const char *s1 = "AGNMKRDE";
  const char *s2 = "PQSTVMFY";
  int l = 8;
 
  printf("s1: %s\n", s1);
  printf("s2: %s\n", s2);
  printf("BLOSUM62\n");
  SS = SCORE_MATRIX_from_file("../matrix/BLOSUM62");
  SCORE_MATRIX_fprint(SS, stdout);
  SCORE_MATRIX_del(SS);
  SS = SCORE_MATRIX_from_file("../matrix/BLOSUM62");
  printf("eval_score: %d\n", SS->eval_score(SS,s1,s2,l));

  printf("DEFAULT CONVERSION TYPE: %d\n", SS->conv_type);
  printf("range_conv (r=57): %d\n", SS->range_conv(SS,s1,l,57));
  printf("matrix_conv\n");
  SD = SS->matrix_conv(SS, NULL, 0);
  SCORE_MATRIX_fprint(SD, stdout);
  printf("eval_score: %d\n", SD->eval_score(SD,s1,s2,l));
  SCORE_MATRIX_del(SD);

  SS->set_conv_type(SS, POS);
  printf("CONVERSION TYPE POS: %d\n", SS->conv_type);
  printf("range_conv (r=57): %d\n", SS->range_conv(SS,s1,l,57));
  printf("matrix_conv\n");
  PS = SS->matrix_conv(SS, s1, l);
  SCORE_MATRIX_fprint(PS, stdout);
  printf("eval_score: %d\n", PS->eval_score(PS,s1,s2,l));
  SCORE_MATRIX_del(PS);

  SS->set_conv_type(SS, AVG);
  printf("CONVERSION TYPE AVG: %d\n", SS->conv_type);
  printf("range_conv (r=57): %d\n", SS->range_conv(SS,s1,l,57));
  printf("matrix_conv\n");
  SD = SS->matrix_conv(SS, NULL, 0);
  SCORE_MATRIX_fprint(SD, stdout);
  printf("eval_score: %d\n", SD->eval_score(SD,s1,s2,l));
 
  SD->set_conv_type(SS, POS);
  printf("CONVERSION TYPE AVG->POS: %d\n", SD->conv_type);
  printf("range_conv (r=57): %d\n", SD->range_conv(SD,s1,l,57));
  printf("matrix_conv\n");
  PD = SD->matrix_conv(SD, s1, l);
  SCORE_MATRIX_fprint(PD, stdout);
  printf("eval_score: %d\n", PD->eval_score(PD,s1,s2,l));

  SCORE_MATRIX_del(SS);
  SCORE_MATRIX_del(SD);
  SCORE_MATRIX_del(PD);

  return EXIT_SUCCESS;
}
