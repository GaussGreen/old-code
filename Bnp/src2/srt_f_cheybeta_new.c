/*--------------------------------------------------------------
        FILE: srt_f_cheybeta_new.c
        PURPOSE: Cheyette beta model interface (new)
        AUTHOR: Dimitri Mayevski
        DATE: 20/09/2002
  --------------------------------------------------------------*/

#include "srt_h_all.h"
#include "srt_h_cheybeta_new.h"

double CheyBeta_Vol(double x, double phi, double sig, double beta, double gamT,
                    double fr, double xcut_u, double xcut_d) {
  double vol, xx;

  xx = x - phi * gamT;
  if (xx > xcut_u)
    xx = xcut_u;
  if (xx < xcut_d)
    xx = xcut_d;

  vol = sig * (1.0 + beta * xx / fr);
  if (vol < 0.0001)
    vol = 0.0001;
  if (vol > 0.03)
    vol = 0.03;

  return vol;
}

Err CheyBeta_Init(SCheyBeta *pmdl, char *und_name) {
  SrtUndPtr und;
  SrtIrDesc *ird;
  TermStruct *ts;
  SrtLst *ls;
  int i;
  IrmTermStructVal *tsv;

  und = lookup_und(und_name);
  if (!und)
    return serror("Underlying not found");

  ird = (SrtIrDesc *)(und->spec_desc);
  ts = (TermStruct *)(ird->ts);

  pmdl->today = get_today_from_underlying(und);
  pmdl->ycname = ird->yc_name;

  /* Calculate nsig */
  for (ls = ts->head, pmdl->nsig = 0; ls; ls = ls->next) {
    tsv = (IrmTermStructVal *)(ls->element->val.pval);
    if ((tsv->val_origin == SIGMA_DATE) || (tsv->val_origin == BOTH_DATE))
      pmdl->nsig++;
  }
  if (pmdl->nsig == 0)
    return serror("No sigma dates found");

  /* Allocate and get sigmas and betas */
  pmdl->sig = (double *)calloc(pmdl->nsig, sizeof(double));
  pmdl->beta = (double *)calloc(pmdl->nsig, sizeof(double));
  pmdl->sigtms = (double *)calloc(pmdl->nsig, sizeof(double));
  if (!pmdl->sig || !pmdl->beta || !pmdl->sigtms)
    return serror("Memory failure");

  for (ls = ts->head, i = 0; ls; ls = ls->next) {
    tsv = (IrmTermStructVal *)(ls->element->val.pval);
    if ((tsv->val_origin == SIGMA_DATE) || (tsv->val_origin == BOTH_DATE)) {
      pmdl->sigtms[i] = tsv->time;
      pmdl->sig[i] = tsv->sig;
      pmdl->beta[i] = tsv->beta;
      i++;
    }
  }

  /* Get lambda */
  for (ls = ts->head; ls; ls = ls->next) {
    tsv = (IrmTermStructVal *)(ls->element->val.pval);
    if ((tsv->val_origin == TAU_DATE) || (tsv->val_origin == BOTH_DATE))
      break;
  }
  if (!ls)
    return serror("No tau dates found");
  pmdl->lambda = 1.0 / tsv->tau;

  return NULL;
}

Err CheyBeta_Free(SCheyBeta *pmdl) {
  free(pmdl->sig);
  free(pmdl->beta);
  free(pmdl->sigtms);
  memset(pmdl, 0, sizeof(SCheyBeta));
  return NULL;
}
