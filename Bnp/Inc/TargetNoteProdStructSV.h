#ifndef _TARGETNOTE_PRODSTRUCT_SV_
#define _TARGETNOTE_PRODSTRUCT_SV_

#include "TargetNoteProdStruct.h"
#include "TargetNoteSV.h"

/* -----------------------------------------------------------------------------------------------------------------
 */
/* The SV auxiliary MC structure */

typedef struct {
  /* LGM SV parameters */
  double *sigma;
  double *alpha;
  double *rho;
  double *lameps;
  double *lvleps;
  /* 2F parameters */
  double *dLGMalpha;
  double *dLGMgamma;
  double *dLGMrho;
  double *dRho2;
  /* derived quantities */
  double *logdff_star;
  double *gam1_star;
  double *gam2_star;
  double *gam1_2_star;
  double *gam2_2_star;
  double *gam12_star;
  /* underlying model */
  LGMSV_model model;
  LGMSVParam Params;
} TARN_MC_AUX_SV;

char *init_TARN_MC_AUX_SV(TARN_MC_AUX_SV **mc_aux_sv);
void free_TARN_MC_AUX_SV(TARN_MC_AUX_SV *ptrTN);

/* The 2FSV auxiliary MC structure */
/* -----------------------------------------------------------------------------------------------------------------
 */

char *alloc_TARN_MC_AUX_SV(TARN_Struct *tarn, TARN_AUX *aux,
                           TARN_MC_AUX *mc_aux);

char *fillTargetNoteFundingSV(int i, TARN_Struct *tarn, TARN_AUX *aux,
                              TARN_MC_AUX *mc_aux, TARN_EVENT *event);

char *fillTargetNoteCouponSV(int iEvt, TARN_AUX *aux, TARN_MC_AUX *mc_aux,
                             TARN_EVENT *event);

char *set_TARN_TimeSteps_SV(TARN_Struct *tarn, TARN_AUX *aux,
                            TARN_MC_AUX *mc_aux);

#endif /* _TARGETNOTE_PRODSTRUCT_SV_ */