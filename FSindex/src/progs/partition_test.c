#include "misclib.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include "partition.h"

EXCEPTION FSexcept_array[1];
EXCEPTION *except;
struct exception_context the_exception_context[1];

int main(int argc, char **argv)
{
  FS_TABLE *ptable;
  const char **sepn;
  int len = 5;
  FILE *fp;
  char *cFSseq;
  const char *seq;
  ULINT seq_test;
  FS_SEQ_t FSseq;

  sepn = mallocec(len*sizeof(char *));

  /* 1. Normal input - should have no problems */
  sepn[0] = "S#T#A#N#I#L#V#M#K#R#D#E#Q#W#F#Y#H#G#P#C";
  sepn[1] = "S#T#A#N#I#L#V#M#K#R#D#E#Q#W#F#Y#H#G#P#C";
  sepn[2] = "STA#N#ILV#M#KR#DEQ#W#FY#H#G#P#C";
  sepn[3] = "STAN#ILVM#KR#DEQ#WFYH#GPC";
  sepn[4] = "STAN#ILVM#KRDEQ#WFYH#GPC";
  if ((ptable = FS_TABLE_init(sepn, len)) == NULL) {
    printf("TEST 1: Normal input - FAILED!\n");
    return EXIT_FAILURE;
  }
  FS_TABLE_del(ptable);
  if ((ptable = FS_TABLE_init(sepn, len)) == NULL) {
    printf("TEST 1: Normal input - FAILED!\n");
    return EXIT_FAILURE;
  }
  FS_TABLE_fprint(ptable, stdout);
  fp = fopen("temp_test.dat", "w");
  FS_TABLE_write(ptable, fp);
  FS_TABLE_del(ptable);
  fclose(fp);
  fp = fopen("temp_test.dat", "r");
  ptable = FS_TABLE_read(fp);
  fclose(fp);
  printf("TEST 1: Normal input - PASSED!\n");
  FS_TABLE_del(ptable);
  
  /* 2. Bad format - '#' at start */
  sepn[3] = "#STAN#ILVM#KR#DEQ#WFYH#GPC";
  if ((ptable = FS_TABLE_init(sepn, len)) != NULL) {
    printf("TEST 2: Bad format - '#' at start - FAILED!\n");
    return EXIT_FAILURE;
  }
  printf("TEST 2: Bad format - '#' at start - PASSED!\n");
  /* 3. Bad format - '#' at end */
  sepn[3] = "STAN#ILVM#KR#DEQ#WFYH#GPC#";
  if ((ptable = FS_TABLE_init(sepn, len)) != NULL) {
    printf("TEST 3: Bad format - '#' at end - FAILED!\n");
    return EXIT_FAILURE;
  }
  printf("TEST 3: Bad format - '#' at end - PASSED!\n");
  /* 4. Bad format - repeating letter at posn 0 */
  sepn[0] = "S#T#A#N#I#L#VVVV#M#K#R#D#E#Q#W#F#Y#H#G#P#C";
  sepn[3] = "STAN#ILVM#KR#DEQ#WFYH#GPC";
  if ((ptable = FS_TABLE_init(sepn, len)) != NULL) {
    printf("TEST 4: Bad format - repeating letter at posn 0 - FAILED!\n");
    return EXIT_FAILURE;
  }
  printf("TEST 4: Bad format - repeating letter at posn 0 - PASSED!\n");
  /* 5. Bad format - mismatching alphabets */
  sepn[0] = "S#T#A#N#I#L#V#M#K#R#D#E#Q#W#F#Y#H#G#P#C";
  sepn[3] = "STAN#ILVMX";
  if ((ptable = FS_TABLE_init(sepn, len)) != NULL) {
    printf("TEST 5: Bad format - mismatching alphabets - FAILED!\n");
    return EXIT_FAILURE;
  }
  printf("TEST 5: Bad format - mismatching alphabets - PASSED!\n");
  /* 6. FS sequence - example 1 */
  sepn[0] = "S#T#A#N#I#L#V#M#K#R#D#E#Q#W#F#Y#H#G#P#C";
  sepn[1] = "S#T#A#N#I#L#V#M#K#R#D#E#Q#W#F#Y#H#G#P#C";
  sepn[2] = "STA#N#ILV#M#KR#DEQ#W#FY#H#G#P#C";
  sepn[3] = "STAN#ILVM#KR#DEQ#WFYH#GPC";
  sepn[4] = "STAN#ILVM#KRDEQ#WFYH#GPC";
  ptable = FS_TABLE_init(sepn, len);
  printf("TEST 6: sequence conversion\n");
  seq = "AMQYG";
  seq_test = 136542;
  if (!FS_SEQ(ptable, seq, 5, &FSseq)) {
    printf("TEST 6: sequence conversion - FAILED!\n");
    return EXIT_FAILURE;
  }
  cFSseq = FS_SEQ_print(ptable, FSseq);
  printf("seq: %s\n", seq);
  printf("seq_test: %ld\n", seq_test);
  printf("FSseq: %ld\n", FSseq);
  printf("%s\n", cFSseq);
  free(cFSseq);
  seq = "AMQYGPE";
  if (!FS_SEQ(ptable, seq, 7, &FSseq)) {
    printf("TEST 6: sequence conversion - FAILED!\n");
    return EXIT_FAILURE;
  }
  cFSseq = FS_SEQ_print(ptable, FSseq);
  printf("seq: %s\n", seq);
  printf("seq_test: %ld\n", seq_test);
  printf("FSseq: %ld\n", FSseq);
  printf("%s\n", cFSseq);
  free(cFSseq);
  seq = "AMQ";
  seq_test = 2142;
  if (!FS_SEQ(ptable, seq, 3, &FSseq)) {
    printf("TEST 6: sequence conversion - FAILED!\n");
    return EXIT_FAILURE;
  }
  cFSseq = FS_SEQ_print(ptable, FSseq);
  printf("seq: %s\n", seq);
  printf("seq_test: %ld\n", seq_test);
  printf("FSseq: %ld\n", FSseq);
  printf("%s\n", cFSseq);
  free(cFSseq);
  seq = "AMQXO";
  printf("seq: %s - bad alphabet\n", seq);
  if (FS_SEQ(ptable, seq, 5, &FSseq)) {
    printf("TEST 6: sequence conversion - FAILED!\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
