#ifndef _TARGET_NOTE_H_
#define _TARGET_NOTE_H_
// ------------------------------------------------------------------------------------------------------------------
// //
//
// TargetNote.h
//
//
#include "TargetNoteProdStruct.h"
#include "utTypes.h"

char *TARN_TimeSteps_2F(TARN_Struct *tarn, TARN_AUX *aux, TARN_MC_AUX *aux_mc);

// ------------------------------------------------------------------------------------------------------------------
// //
//
// Different coupon PV calculators
//
char *fillTargetNoteFunding(int i, int flag, long *lvEventDates,
                            long *lvFundStartDates, long *lvFundEndDates,
                            double *dvFundMargin, double *dvFundSpread,
                            SrtBasisCode basisFloat, TARN_MC_AUX *aux_mc,
                            TARN_EVENT *targetNotePtr);

char *fillTargetNoteCoupon(int i, TARN_Struct *tarn, TARN_AUX *aux,
                           TARN_MC_AUX *mc_aux, TARN_EVENT *targetNotePtr);

// ------------------------------------------------------------------------------------------------------------------
// //
//
// LGM2F pricer
char *TargetNoteMC_2F(char *szUnd, TARN_Struct *tarn, TARN_AUX *aux,
                      double **prod_val);

char *TargetNoteMC_KO_2F(char *szUnd, TARN_Struct *tarn, TARN_AUX *aux,
                         double **prod_val);

#endif // #ifndef _TARGET_NOTE_H_