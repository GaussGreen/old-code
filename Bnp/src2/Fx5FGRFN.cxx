
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

Err grfn_payoff_4_5dfx_mc(
    /* Event */
    double evt_date, double evt_time, void *func_parm,
    /* Market data */
    double spot_fx, double R1D, double R2D, double R1F, double R2F, double Z,
    /* Results */
    int num_col, double *res, int *stop_path) {
  GRFNPARMMC5F total;
  GRFNCOMMSTRUCT global;
  FIRSTMktAtT *local;
  int l;
  double temp;
  Err err;

  if (func_parm == NULL) {
    memset(res, 0, num_col * sizeof(double));
    *stop_path = 0;
    return NULL;
  } else {
    err = NULL;
  }

  memset(res, 0, num_col * sizeof(double));

  /* Get the event */
  total = (GRFNPARMMC5F)func_parm;
  global = total->global;
  local = total->local;

  /* Calc market data */
  local->smp.und[total->fx_idx].sv[SPOT] = Z;

  if (total->do_dom) {
    for (l = 0; l < total->num_dom_df; l++) {
      local->evt->df[total->dom_idx][l] =
          total->dom_dff[l] *
          exp(total->dom_gam12[l] - total->dom_gam1[l] * R1D +
              -total->dom_gam2[l] * R2D);
    }
  }
  if (total->do_fx) {
    for (l = 0; l < total->num_fx_df; l++) {
      local->evt->df[total->fx_idx][l] =
          total->fx_dff[l] * exp(total->fx_gam12[l] - total->fx_gam1[l] * R1D +
                                 -total->fx_gam2[l] * R2D);
    }
  }
  if (total->do_for) {
    for (l = 0; l < total->num_for_df; l++) {
      local->evt->df[total->for_idx][l] =
          total->for_dff[l] *
          exp(total->for_gam12[l] - total->for_gam1[l] * R1F +
              -total->for_gam2[l] * R2F);
    }
  }

  err = FIRSTEvalEvent(global, local, num_col, 2, NULL, NULL, res, &temp);

  *stop_path = 0;
  return err;
}