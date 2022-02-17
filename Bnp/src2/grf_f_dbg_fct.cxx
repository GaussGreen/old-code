/************************************************

        functions for printing Grfn Structures

**********************************************/

#include "grf_h_all.h"

void print_sprdsht(long numeventdates, Date *eventdates, GrfnCell **sprdsht,
                   FILE *f) {
  int i = 0;

  if (!sprdsht || !f)
    return;

  for (i = 0; i < numeventdates; i++) {
    switch (sprdsht[i][0].type) {
    case GRFNSCELL:
      fprintf(f, "row: %2d date %d cmd: %s\n", i, eventdates[i],
              sprdsht[i][0].sval);
      break;
    default:
      fprintf(f, "row: %2d date %d cmd:%.4lf\n", i, eventdates[i],
              sprdsht[i][0].dval);
      break;
    }
  }
}

/* sets *answer to be equal to d        , prints out ref rate name */
Err grf_f_testhistfct(Date d, char *ref_rate_name, double *answer) {
  if (answer) {
    *answer = 0.0;
    return NULL;
  } else
    return serror("grf_f_testhistfct called with NULL *answer");
}
