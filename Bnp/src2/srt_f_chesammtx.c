/* ================================================================================

   FILENAME:      srt_f_chesammtx.c

   PURPOSE:       different limits for the Phi varaible in Cheyette tree

   ================================================================================
 */

#include "math.h"
#include "srt_h_all.h"

/*-----------------------------------------------------------------------
  FUNCNAME			:srt_f_chephilim
  AUTHOR			:E.Auld
  DESCRIPTION		:determine max and min phi for a given sample
                                         tminf and trinf should correspond to
stp does both "trapezoidal" and "rectangular" sampling in phi MODIFIES
:sam->phi is set to the phi centered on *phimin and *phimax are also set if they
are not NULL
---------------------------------------------------------------------------*/

/* spacing in phi */
#define PHITRAPEZOID 0
#define PHIRECTANGLE 1
#define PHIMETHOD PHIRECTANGLE

void srt_f_chephilim(SrtStpPtr top, SrtStpPtr stp, SrtIRMTmInf *tminf,
                     SrtSample *sam, double *phimin, double *phimax)

{
  double centerphi;

  switch (PHIMETHOD) {
  case PHITRAPEZOID: {
    srt_f_chelinrphi(top, stp, sam);
    centerphi = samptr_get(sam, 0, PHI);
    break;
  }
  case PHIRECTANGLE:
  default: {
    centerphi = sam_get(tminf->fwd_sam, 0, PHI);
    break;
  }
  }

  if (phimax) {
    *phimax = centerphi * DMAX(2.0, stp->time);
  }
  if (phimin) {
    *phimin = centerphi * DMIN(0.5, 1.0 / stp->time);
  }
  return;
}

/****************************************************************************
  FUNCNAME        :srt_f_chesammatrix
  AUTHOR          :E.Auld
  DESCRIPTION     :populate a matrix of samples in the cheyette model;
                   indices of arrays are determined inside stp.
                   Does various different kinds of spacing in r and
                   phi depending on internal flags (see below).
  MODIFIES        :cur_r and cur_phi
******************************************************************************/

/* creates matrix of state variables in cheyette tree */
void srt_f_chesammatrix(SrtStpPtr top,    /* ptr to top of lst */
                        SrtStpPtr stp,    /* ptr to current step in lst */
                        double *cur_r,    /* space for rates to be stored */
                        double **cur_phi, /* space for phis to be stored */
                        SrtGrfnParam *grfnparam /* info about model */
)

{
  int i, j, max_ri, min_ri;
  double logu;
  double logr;
  double loguphi;
  double logphi;
  double phimin;
  double phimax;
  SrtSample sam;
  int prflg = 0;
  SrtCheTreInf *trinf;
  SrtIRMTmInf *tminf;

  /* get trinf and tminf */
  trinf = stp->trinf;
  tminf = stp->tminf[0];

  /* if only one node */
  /* ie. if step 1 */
  if (trinf->max_r_index == trinf->min_r_index) {
    cur_r[trinf->min_r_index] = sam_get(tminf->fwd_sam, 0, SHORT_RATE);
    cur_phi[trinf->min_r_index][trinf->min_phi_index] = 0.0;
    return;
  }

  /* otherwise */

  /* modified AS */
  max_ri = IMAX(trinf->max_r_index,
                ((SrtCheTreInf *)(stp->prev->trinf))->max_r_index + 1);
  min_ri = IMIN(trinf->min_r_index,
                ((SrtCheTreInf *)(stp->prev->trinf))->min_r_index - 1);

  /** generate rates **/
  for (i = min_ri, logu = log(trinf->u),
      logr = log(sam_get(tminf->fwd_sam, 0, F_0_t)) + min_ri * logu;
       i <= max_ri; i++, logr += logu) {
    cur_r[i] = exp(logr);
  }

  /** generate phis **/
  for (i = min_ri; i <= max_ri; i++) {
    sam_get(sam, 0, SHORT_RATE) = cur_r[i];
    srt_f_chephilim(top, stp, tminf, &sam, &phimin, &phimax);

    if (trinf->max_phi_index == trinf->min_phi_index) {
      loguphi = 0;
    } else {
      loguphi =
          log(phimax / phimin) / (trinf->max_phi_index - trinf->min_phi_index);
    }

    for (j = trinf->min_phi_index, logphi = j * loguphi;
         j <= trinf->max_phi_index; j++, logphi += loguphi) {
      cur_phi[i][j] = phimin * exp(logphi);
    }
  }
  return;
}
