/* =================================================================================

   FILENAME: srt_f_cheydynamics.c

   PURPOSE:  Give all the Equations required to descretise the model:
                                - drift of r
                                - expectation of r at t+1 knowing t
                                - local volatility of r
                                - varaince of r at t+1 knowing t (and same for
   phi)
   =================================================================================
 */

#include "srt_h_all.h"
#include "srt_h_cheydynamics.h"

/* --------------------------------------------------------------------------
  FUNCNAME        :srt_f_chedrfatsam
  AUTHOR          :E.Auld
  DESCRIPTION     :in Chey  , computes expectation of sv at stp
                                   and stores in drift_sam
  MODIFIES        :drift_sam
 ---------------------------------------------------------------------------- */

void srt_f_chedrfatsam(SrtStpPtr stp, SrtSample *cur_sam, SrtSample *drift_sam,
                       int index)

{
  double mu;
  double dphi;
  SrtIRMTmInf *tminf;
  SrtIRMTmInf *nxttminf;

  if (!stp->next) {
    return;
  }

  tminf = stp->tminf[index];
  nxttminf = stp->next->tminf[index];

  /* The STATEVAR is X = r - f(0  ,t) */
  samptr_get(cur_sam, 0, STATEVAR) =
      samptr_get(cur_sam, 0, SHORT_RATE) - sam_get(tminf->fwd_sam, 0, F_0_t);
  mu = -tminf->ev.onef.lambda * samptr_get(cur_sam, 0, STATEVAR) +
       samptr_get(cur_sam, 0, PHI);
  mu *= stp->delta_t;

  dphi = samptr_get(cur_sam, 0, SHORT_RATE) * tminf->ev.onef.sig;
  dphi =
      dphi * dphi - 2.0 * tminf->ev.onef.lambda * samptr_get(cur_sam, 0, PHI);
  dphi *= stp->delta_t;

  samptr_get(drift_sam, 0, STATEVAR) = samptr_get(cur_sam, 0, STATEVAR) + mu;
  samptr_get(drift_sam, 0, SHORT_RATE) =
      sam_get(nxttminf->fwd_sam, 0, F_0_t) + samptr_get(drift_sam, 0, STATEVAR);

  samptr_get(drift_sam, 0, PHI) = samptr_get(cur_sam, 0, PHI) + dphi;
  return;
}

/* -------------------------------------------------------------------------------
 */

/* -------------------------------------------------------------------------------

   FUNCNAME        :srt_f_chevaratsam

   DESCRIPTION     :calculates variance at stp  , cur_sam  , that is
sig^2*r^2*dt (stp should correspond to cur_sam)

   MODIFIES        :var_at_sam

   CALL            :

  AMENDMENTS      :
  Reference       :
  Author          :
  Date            :
  Description     :
  -------------------------------------------------------------------------------
*/

void srt_f_chevaratsam(SrtStpPtr stp, SrtSample *cur_sam, double *var_at_sam,
                       int index) {
  double var;
  SrtIRMTmInf *tminf;

  tminf = (SrtIRMTmInf *)stp->tminf[index];

  /* see comment in srt_f_chedrfatsam concerning the index in the following
     macro call */

  var = samptr_get(cur_sam, 0, SHORT_RATE) * tminf->ev.onef.sig;
  var *= var * stp->delta_t;
  *var_at_sam = var;
  return;
}

/* -------------------------------------------------------------------------------
 */
