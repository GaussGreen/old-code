#ifndef SRTGRFNMAINQTO_H
#define SRTGRFNMAINQTO_H

#include "grf_h_all.h"
#include "srt_h_all.h"

char *SrtGrfnQuantoPDE(char *und3dfx, int numeventdates, long *eventdates,
                       long tableauRows, long tableauCols,
                       char ***tableauStrings, int **tableauMask, long auxWidth,
                       long *auxLen, double **aux, long num_stp, long num_stpx,
                       int *num_prod, double **prod_val);

char *SrtGrfnQuantoPDEWithCorrelTS(char *und3dfx, int numeventdates,
                                   long *eventdates, long tableauRows,
                                   long tableauCols, char ***tableauStrings,
                                   int **tableauMask, long auxWidth,
                                   long *auxLen, double **aux, long num_stp,
                                   long num_stpx, int *num_prod,
                                   double **prod_val);

char *SrtGrfnQuantoPDEWithCorrelTS2(char *und3dfx, int numeventdates,
                                    long *eventdates, long tableauRows,
                                    long tableauCols, char ***tableauStrings,
                                    int **tableauMask, long auxWidth,
                                    long *auxLen, double **aux, long num_stp,
                                    long num_stpx, int *num_prod,
                                    double **prod_val);

#endif
