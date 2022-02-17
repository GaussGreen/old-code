/* ==========================================================================
   FILE_NAME:	FXLGMSVGrfn.c

   PURPOSE:		Monte Carlo FX Black-Scholes / IR Dom LGMSV2F / IR for
   LGMSV2F

   DATE:		01/10/03
   ========================================================================== */

#include "FXLGMSVGrfn.h"
#include "FXLGMSVMC.h"
#include "FXLGMSVUnd.h"
#include "MCEBOptimisation.h"
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

Err payoff_fxlgmsv_mc(long path_index, double evt_date, double evt_time,
                      void *func_parm,

                      // Domestic
                      double domft1, double domft2, double dompsi1,
                      double dompsi2, double dompsi12, double domv,

                      // Foreign
                      double forft1, double forft2, double forpsi1,
                      double forpsi2, double forpsi12, double forv,

                      // FX
                      double fx_spot,

                      /* Vector of results to be updated */
                      int nprod,
                      /* Result	*/
                      double *prod_val, int *stop_path) {
  GRFNPARMMCFXLGMSV total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  static double temp, tmpx1, tmpx2, df;
  static int k;

  static Err err = NULL;

  /*	Get the event		*/
  total = (GRFNPARMMCFXLGMSV)func_parm;
  global = total->global;
  local = total->local;

  memset(prod_val, 0, nprod * sizeof(double));

  /*	Evaluation */
  /* Calc market data */
  local->smp.und[total->fx_idx].sv[SPOT] = fx_spot;

  if (total->do_dom) {
    for (k = 0; k < total->num_dom_df; k++) {
      local->evt->df[total->dom_idx][k] =
          exp(total->dom_dff[k] - total->dom_gam1[k] * domft1 -
              total->dom_gam2[k] * domft2 - total->dom_gam1_2[k] * dompsi1 -
              total->dom_gam2_2[k] * dompsi2 - total->dom_gam12[k] * dompsi12);
    }

    // Previous code:
    // local->smp.und[0].sv[PHI1] = dompsi1;
    // local->smp.und[0].sv[PHI2] = dompsi2;
    // local->smp.und[0].sv[CROSSPHI] = dompsi12;
    // local->smp.und[0].sv[SIGMA] = domv;
  }
  if (total->do_fx) {
    for (k = 0; k < total->num_fx_df; k++) {
      local->evt->df[total->fx_idx][k] =
          exp(total->fx_dff[k] - total->fx_gam1[k] * domft1 -
              total->fx_gam2[k] * domft2 - total->fx_gam1_2[k] * dompsi1 -
              total->fx_gam2_2[k] * dompsi2 - total->fx_gam12[k] * dompsi12);
    }
  }
  if (total->do_for) {
    for (k = 0; k < total->num_for_df; k++) {
      local->evt->df[total->for_idx][k] =
          exp(total->for_dff[k] - total->for_gam1[k] * forft1 -
              total->for_gam2[k] * forft2 - total->for_gam1_2[k] * forpsi1 -
              total->for_gam2_2[k] * forpsi2 - total->for_gam12[k] * forpsi12);
    }

    // Put the state variables on the foreign stochvol
    local->smp.und[0].sv[PHI1] = dompsi1;
    local->smp.und[0].sv[PHI2] = dompsi2;
    local->smp.und[0].sv[CROSSPHI] = dompsi12;
    local->smp.und[0].sv[SIGMA] = domv;
  }

  err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, prod_val, &temp);

  return err;
}

Err payoff_qtolgmsv1f_mc(long path_index, double evt_date, double evt_time,
                         void *func_parm,

                         // Domestic
                         double domft, double dompsi,

                         // Foreign
                         double forft, double forpsi, double forv,

                         /* Vector of results to be updated */
                         int nprod,
                         /* Result	*/
                         double *prod_val, int *stop_path) {
  GRFNPARMMCFXLGMSV total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  static double temp, tmpx1, tmpx2, df;
  static int k;

  static Err err = NULL;

  /*	Get the event		*/
  total = (GRFNPARMMCFXLGMSV)func_parm;
  global = total->global;
  local = total->local;

  memset(prod_val, 0, nprod * sizeof(double));

  /*	Evaluation */
  /* Calc market data */
  local->smp.und[total->fx_idx].sv[SPOT] = -999999;

  if (total->do_dom) {
    for (k = 0; k < total->num_dom_df; k++) {
      local->evt->df[total->dom_idx][k] =
          exp(total->dom_dff[k] - total->dom_gam1[k] * domft -
              total->dom_gam1_2[k] * dompsi);
    }

    // Previous code:
    // local->smp.und[0].sv[PHI1] = dompsi;
  }
  if (total->do_fx) {
    for (k = 0; k < total->num_fx_df; k++) {
      local->evt->df[total->fx_idx][k] =
          exp(total->fx_dff[k] - total->fx_gam1[k] * domft -
              total->fx_gam1_2[k] * dompsi);
    }
  }
  if (total->do_for) {
    for (k = 0; k < total->num_for_df; k++) {
      local->evt->df[total->for_idx][k] =
          exp(total->for_dff[k] - total->for_gam1[k] * forft -
              total->for_gam1_2[k] * forpsi);
    }

    // Put the state variables on the foreign stochvol
    local->smp.und[0].sv[PHI1] = forpsi;
    local->smp.und[0].sv[SIGMA] = forv;
  }

  err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, prod_val, &temp);

  return err;
}

Err payoff_qtolgmsv2f_mc(long path_index, double evt_date, double evt_time,
                         void *func_parm,

                         // Domestic
                         double domft, double dompsi,

                         // Foreign
                         double forft1, double forft2, double forpsi1,
                         double forpsi2, double forpsi12, double forv,

                         /* Vector of results to be updated */
                         int nprod,
                         /* Result	*/
                         double *prod_val, int *stop_path) {
  GRFNPARMMCFXLGMSV total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;

  static double temp, tmpx1, tmpx2, df;
  static int k;

  static Err err = NULL;

  // TEST J.M.L.
  /*
  smessage (" forft1 = %f\n"  , (float)forft1);
  smessage (" forft2 = %f\n"  , (float)forft2);
  smessage (" forpsi1 = %f\n"  , (float)forpsi1);
  smessage (" forpsi2 = %f\n"  , (float)forpsi2);
  smessage (" forpsi12 = %f\n"  , (float)forpsi12);
  smessage (" forv = %f\n"  , (float)forv);
  */

  /*	Get the event		*/
  total = (GRFNPARMMCFXLGMSV)func_parm;

  // TEST J.M.L.
  /*
  if (total == NULL)
          smessage (" total = NULL\n");
  else
          smessage (" total = %X\n"  , total);
  */

  global = total->global;
  local = total->local;

  memset(prod_val, 0, nprod * sizeof(double));

  /*	Evaluation */
  /* Calc market data */
  local->smp.und[total->fx_idx].sv[SPOT] = -999999;

  // TEST J.M.L.
  /*
  if (total != NULL)
  {
          smessage (" total->do_dom = %f\n"  , (float)total->do_dom);
          smessage (" total->do_fx = %f\n"  , (float)total->do_fx);
          smessage (" total->do_for = %f\n"  , (float)total->do_for);
  }
  */

  if (total->do_dom) {
    for (k = 0; k < total->num_dom_df; k++) {
      local->evt->df[total->dom_idx][k] =
          exp(total->dom_dff[k] - total->dom_gam1[k] * domft -
              total->dom_gam1_2[k] * dompsi);
    }

    // Previous code:
    // local->smp.und[0].sv[PHI1] = dompsi;
  }
  if (total->do_fx) {
    for (k = 0; k < total->num_fx_df; k++) {
      local->evt->df[total->fx_idx][k] =
          exp(total->fx_dff[k] - total->fx_gam1[k] * domft -
              total->fx_gam1_2[k] * dompsi);
    }
  }
  if (total->do_for) {
    for (k = 0; k < total->num_for_df; k++) {
      local->evt->df[total->for_idx][k] =
          exp(total->for_dff[k] - total->for_gam1[k] * forft1 -
              total->for_gam2[k] * forft2 - total->for_gam1_2[k] * forpsi1 -
              total->for_gam2_2[k] * forpsi2 - total->for_gam12[k] * forpsi12);

      // TEST J.M.L.
      /*
      smessage (" total->for_dff[%d]=%f\n"  , k  , (float)total->for_dff[k]);
      smessage (" total->for_gam1[%d]=%f\n"  , k  , (float)total->for_gam1[k]);
      smessage (" total->for_gam2[%d]=%f\n"  , k  , (float)total->for_gam2[k]);
      smessage (" total->for_gam1_2[%d]=%f\n"  , k  ,
      (float)total->for_gam1_2[k]); smessage (" total->for_gam2_2[%d]=%f\n"  , k
      , (float)total->for_gam2_2[k]); smessage (" total->for_gam12[%d]=%f\n"  ,
      k  , (float)total->for_gam12[k]); smessage (" DF[%d]=%f\n"  , k  ,
      (float)local->evt->df[total->for_idx][k]);
      */
    }

    // Put the state variables on the foreign stochvol
    local->smp.und[0].sv[PHI1] = forpsi1;
    local->smp.und[0].sv[PHI2] = forpsi2;
    local->smp.und[0].sv[CROSSPHI] = forpsi12;
    local->smp.und[0].sv[SIGMA] = forv;
  }

  err = FIRSTEvalEvent(global, local, nprod, 2, NULL, NULL, prod_val, &temp);

  return err;
}
