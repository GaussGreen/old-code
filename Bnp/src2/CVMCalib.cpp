//---------------------------------------------------------------------
//-------------	Continuous Volatility Model	(CVM)	-------------------
//---------------------------------------------------------------------
#include "CVMCalib.h"
#include "CVM.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

static cvm_str *cvm_static = NULL;
static int calDateIndex_static = 0;
static long calStart_static = 0;
static double bump_vol_static = 0.0001;
static double bump_rho_static = 0.01;

Err cvm_init_static_struct(cvm_str *cvm) {
  Err err = NULL;

  cvm_static = (cvm_str *)calloc(1, sizeof(cvm_str));
  err = cvm_copy_struct(cvm_static, cvm);
  return err;
}

Err cvm_free_static_struct() {
  Err err = NULL;

  err = cvm_free_struct(cvm_static);
  free(cvm_static);

  return err;
}

Err cvm_vol_for_Levenberg_Marquardt(double end, double *pillarsVols,
                                    double *vol, double *gradient,
                                    int nPillars) {
  Err err = NULL;
  int i;
  double *pillarsVolsLoc = NULL;
  double voltemp;
  int out_range;

  out_range = 1;
  if (cvm_static->allow_negativeVols == 0) {
    for (i = 0; i < nPillars; ++i) {
      if (pillarsVols[i + 1] < 0) {
        out_range = 0;
        i = nPillars;
      }
    }
  }

  if (out_range == 0) {
    *vol = -1000000;
    for (i = 0; i < nPillars; i++) {
      gradient[i + 1] = 1000000;
    }
  } else {

    pillarsVolsLoc = dvector(0, nPillars - 1);
    for (i = 0; i < nPillars; ++i) {
      pillarsVolsLoc[i] = pillarsVols[i + 1];
    }

    err = cvm_set_pillars_vols(calDateIndex_static, nPillars, pillarsVolsLoc,
                               cvm_static);
    if (err)
      goto FREE_RETURN;
    err = cvm_vol(calStart_static, end, vol, cvm_static);
    if (err)
      goto FREE_RETURN;
    *vol = *vol * 10000.0;

    for (i = 0; i < nPillars; ++i) {
      pillarsVolsLoc[i] = pillarsVolsLoc[i] + bump_vol_static;
      err = cvm_set_pillars_vols(calDateIndex_static, nPillars, pillarsVolsLoc,
                                 cvm_static);
      if (err)
        goto FREE_RETURN;
      err = cvm_vol(calStart_static, end, &voltemp, cvm_static);
      if (err)
        goto FREE_RETURN;
      voltemp = voltemp * 10000.0;
      pillarsVolsLoc[i] = pillarsVolsLoc[i] - bump_vol_static;
      gradient[i + 1] = (voltemp - *vol) / bump_vol_static;
    }
    err = cvm_set_pillars_vols(calDateIndex_static, nPillars, pillarsVolsLoc,
                               cvm_static);
  }

FREE_RETURN:
  if (pillarsVolsLoc) {
    free_dvector(pillarsVolsLoc, 0, nPillars - 1);
    pillarsVolsLoc = NULL;
  }
  return err;
}
/*
Err bump_correl_matrix(int nPillars  , double **correlMat_source  , double
**correlMat_target  , double *bumpCorrel)
{
        Err err=NULL;
        int i  ,k;

        correlMat_target[0][0] = 1.0;
        for(i=1;i<nPillars;++i)
        {
                correlMat_target[i][i] = 1.0;
                for(k=0;k<nPillars-i;++k)
                {
                        correlMat_target[k][i+k] = DMAX(-0.9999  , DMIN(0.9999
, correlMat_source[k][i+k] + bumpCorrel[i-1])); correlMat_target[i+k][k] =
correlMat_target[k][i+k];
                }
        }

        PositiveMatrix(correlMat_target  , nPillars);

        return err;
}
*/

Err bump_correl_matrix(int nPillars, double **correlMat_source,
                       double **correlMat_target, double *bumpCorrel) {
  Err err = NULL;
  int i, k;

  for (i = 0; i < nPillars; ++i) {
    for (k = 0; k < nPillars; ++k) {
      correlMat_target[k][i] =
          DMAX(-0.999, DMIN(0.999, correlMat_source[k][i] + bumpCorrel[0]));
    }
  }

  PositiveMatrix(correlMat_target, nPillars);

  return err;
}

Err cvm_vol_for_Levenberg_Marquardt_for_correl(double end, double *bumpCorrel,
                                               double *vol, double *gradient,
                                               int nPillarsMinusOne) {
  Err err = NULL;
  int i;
  double *bump = NULL;
  double **correlLoc = NULL;
  double voltemp;

  bump = dvector(0, nPillarsMinusOne - 1);
  for (i = 0; i < nPillarsMinusOne; ++i) {
    bump[i] = bumpCorrel[i + 1];
  }

  correlLoc = dmatrix(0, nPillarsMinusOne, 0, nPillarsMinusOne);
  err = bump_correl_matrix(cvm_static->nbOfPillars,
                           cvm_static->pillars_correls[calDateIndex_static],
                           correlLoc, bump);
  PositiveMatrix(correlLoc, cvm_static->nbOfPillars);

  err = cvm_set_pillars_correls(calDateIndex_static, cvm_static->nbOfPillars,
                                correlLoc, cvm_static);
  if (err)
    goto FREE_RETURN;
  err = cvm_vol(calStart_static, end, vol, cvm_static);
  if (err)
    goto FREE_RETURN;
  *vol = *vol * 10000.0;

  for (i = 0; i < nPillarsMinusOne; ++i) {
    bump[i] = bump[i] + bump_rho_static;
    err =
        bump_correl_matrix(cvm_static->nbOfPillars, correlLoc, correlLoc, bump);
    err = cvm_set_pillars_correls(calDateIndex_static, cvm_static->nbOfPillars,
                                  correlLoc, cvm_static);
    if (err)
      goto FREE_RETURN;
    err = cvm_vol(calStart_static, end, &voltemp, cvm_static);
    if (err)
      goto FREE_RETURN;
    voltemp = voltemp * 10000.0;
    bump[i] = bump[i] - bump_rho_static;
    err =
        bump_correl_matrix(cvm_static->nbOfPillars, correlLoc, correlLoc, bump);
    gradient[i + 1] = (voltemp - *vol) / bump_rho_static;
  }

  err = bump_correl_matrix(cvm_static->nbOfPillars, correlLoc, correlLoc, bump);
  err = cvm_set_pillars_correls(calDateIndex_static, cvm_static->nbOfPillars,
                                correlLoc, cvm_static);

FREE_RETURN:
  if (correlLoc) {
    free_dmatrix(correlLoc, 0, cvm_static->nbOfPillars - 1, 0,
                 cvm_static->nbOfPillars - 1);
    correlLoc = NULL;
  }

  if (bump) {
    free_dvector(bump, 0, nPillarsMinusOne - 1);
    bump = NULL;
  }

  return err;
}

Err cvm_calib_pillarsVols(cvm_str *cvm) {
  int i, calDateIndex;
  Err err = NULL;
  int nbr_iter = 25;
  double fitting_error;
  long start, end;
  double *ends = NULL;
  double *market_vols = NULL;
  double *model_vols = NULL;
  double *weights = NULL;
  long *use_param = NULL;
  double *pillarsVols = NULL;

  err = cvm_init_static_struct(cvm);
  if (err)
    goto FREE_RETURN;

  ends = dvector(1, cvm->nbOfUnderlyingTenors);
  market_vols = dvector(1, cvm->nbOfUnderlyingTenors);
  model_vols = dvector(1, cvm->nbOfUnderlyingTenors);
  weights = dvector(1, cvm->nbOfUnderlyingTenors);
  use_param = lngvector(1, cvm->nbOfUnderlyingTenors);
  pillarsVols = dvector(1, cvm->nbOfPillars);

  for (calDateIndex = 0; calDateIndex < cvm->nbOfOptionTenors; ++calDateIndex) {
    for (i = 0; i < cvm->nbOfPillars; ++i) {
      if (calDateIndex > 0) {
        pillarsVols[i + 1] = cvm->pillars_vols_norm[calDateIndex - 1][i];
      } else {
        pillarsVols[i + 1] = cvm->market_vols[calDateIndex][0];
      }
    }

    calDateIndex_static = calDateIndex;
    err = add_tenor(cvm->spot, cvm->option_tenors[calDateIndex],
                    MODIFIED_SUCCEEDING, &start);
    if (err)
      goto FREE_RETURN;
    calStart_static = start;

    for (i = 0; i < cvm->nbOfUnderlyingTenors; ++i) {
      err = add_tenor(start, cvm->underlying_tenors[i], NO_BUSDAY_CONVENTION,
                      &end);
      if (err)
        goto FREE_RETURN;
      ends[i + 1] = end;
      market_vols[i + 1] = 10000.0 * cvm->market_vols[calDateIndex][i];
      model_vols[i + 1] = 10000.0 * cvm->model_vols[calDateIndex][i];
      weights[i + 1] = 1.0;
      use_param[i + 1] = 1;
    }

    err = levenberg_marquardt_select(
        ends, market_vols, weights, cvm->nbOfUnderlyingTenors, pillarsVols,
        use_param, cvm->nbOfPillars, nbr_iter, cvm_vol_for_Levenberg_Marquardt,
        &fitting_error);
    if (err)
      goto FREE_RETURN;
  }

  err = cvm_copy_struct(cvm, cvm_static);

  err = cvm_free_static_struct();

FREE_RETURN:

  if (err) {
    cvm_free_static_struct();
  }

  if (ends)
    free_dvector(ends, 1, cvm->nbOfUnderlyingTenors);
  if (market_vols)
    free_dvector(market_vols, 1, cvm->nbOfUnderlyingTenors);
  if (model_vols)
    free_dvector(model_vols, 1, cvm->nbOfUnderlyingTenors);
  if (pillarsVols)
    free_dvector(pillarsVols, 1, cvm->nbOfPillars);
  if (weights)
    free_dvector(weights, 1, cvm->nbOfUnderlyingTenors);
  if (use_param)
    free_lngvector(use_param, 1, cvm->nbOfUnderlyingTenors);

  return err;
}

Err cvm_calib_pillarsVolsAndCorrels(cvm_str *cvm) {
  int i, calDateIndex, calibcorrel;
  Err err = NULL;
  int nbr_iter = 25;
  double fitting_error;
  long start, end;
  double *ends = NULL;
  double *market_vols = NULL;
  double *model_vols = NULL;
  double *weights = NULL;
  long *use_param = NULL;
  double *pillarsVols = NULL;
  double *correlBumps = NULL;

  err = cvm_init_static_struct(cvm);
  if (err)
    goto FREE_RETURN;

  ends = dvector(1, cvm->nbOfUnderlyingTenors);
  market_vols = dvector(1, cvm->nbOfUnderlyingTenors);
  model_vols = dvector(1, cvm->nbOfUnderlyingTenors);
  weights = dvector(1, cvm->nbOfUnderlyingTenors);
  use_param = lngvector(1, cvm->nbOfPillars);
  pillarsVols = dvector(1, cvm->nbOfPillars);
  correlBumps = dvector(1, cvm->nbOfPillars - 1);

  calibcorrel = 0;

  for (calDateIndex = 0; calDateIndex < cvm->nbOfOptionTenors; ++calDateIndex) {
    for (i = 0; i < cvm->nbOfPillars; ++i) {
      if (calDateIndex > 0) {
        pillarsVols[i + 1] = cvm_static->pillars_vols_norm[calDateIndex - 1][i];
      } else {
        pillarsVols[i + 1] = cvm->market_vols[calDateIndex][i];
      }
    }

    calDateIndex_static = calDateIndex;
    err = add_tenor(cvm->spot, cvm->option_tenors[calDateIndex],
                    MODIFIED_SUCCEEDING, &start);
    if (err)
      goto FREE_RETURN;
    calStart_static = start;

    for (i = 0; i < cvm->nbOfUnderlyingTenors; ++i) {
      err = add_tenor(start, cvm->underlying_tenors[i], NO_BUSDAY_CONVENTION,
                      &end);
      if (err)
        goto FREE_RETURN;
      ends[i + 1] = end;
      market_vols[i + 1] = 10000.0 * cvm->market_vols[calDateIndex][i];
      model_vols[i + 1] = 10000.0 * cvm->model_vols[calDateIndex][i];
      weights[i + 1] = 1.0;
    }

    for (i = 0; i < cvm->nbOfPillars; ++i) {
      use_param[i + 1] = 1;
      correlBumps[i + 1] = 0;
    }
    use_param[cvm->nbOfPillars] = 1;

    err = levenberg_marquardt_select(
        ends, market_vols, weights, cvm->nbOfUnderlyingTenors, pillarsVols,
        use_param, cvm->nbOfPillars, nbr_iter, cvm_vol_for_Levenberg_Marquardt,
        &fitting_error);

    if ((fabs(fitting_error) > 1e-10) && (calibcorrel == 0)) {
      err = levenberg_marquardt_select(
          ends, market_vols, weights, cvm->nbOfUnderlyingTenors, correlBumps,
          use_param, cvm->nbOfPillars - 1, nbr_iter,
          cvm_vol_for_Levenberg_Marquardt_for_correl, &fitting_error);
      calibcorrel = 1;
      calDateIndex = calDateIndex - 1;
    } else {
      calibcorrel = 0;
    }

    if (err)
      goto FREE_RETURN;
  }

  err = cvm_copy_struct(cvm, cvm_static);

  err = cvm_free_static_struct();

FREE_RETURN:

  if (err) {
    cvm_free_static_struct();
  }

  if (ends)
    free_dvector(ends, 1, cvm->nbOfUnderlyingTenors);
  if (market_vols)
    free_dvector(market_vols, 1, cvm->nbOfUnderlyingTenors);
  if (model_vols)
    free_dvector(model_vols, 1, cvm->nbOfUnderlyingTenors);
  if (pillarsVols)
    free_dvector(pillarsVols, 1, cvm->nbOfPillars);
  if (weights)
    free_dvector(weights, 1, cvm->nbOfUnderlyingTenors);
  if (use_param)
    free_lngvector(use_param, 1, cvm->nbOfPillars);

  return err;
}

Err cvm_calib_und(int nbOptionTenors, char **option_tenors,
                  int nbUnderlyingTenors, char **underlying_tenors,
                  double **market_vols, cvm_str *cvm) {
  Err err = NULL;
  int i, j;

  cvm->nbOfOptionTenors = nbOptionTenors;
  cvm->nbOfUnderlyingTenors = nbUnderlyingTenors;

  cvm->option_tenors =
      (char256 *)calloc(cvm->nbOfOptionTenors, sizeof(char256));
  for (i = 0; i < cvm->nbOfOptionTenors; ++i) {
    strcpy(cvm->option_tenors[i], option_tenors[i]);
  }

  cvm->underlying_tenors =
      (char256 *)calloc(cvm->nbOfUnderlyingTenors, sizeof(char256));
  for (i = 0; i < cvm->nbOfUnderlyingTenors; ++i) {
    strcpy(cvm->underlying_tenors[i], underlying_tenors[i]);
  }

  cvm->market_vols =
      dmatrix(0, cvm->nbOfOptionTenors - 1, 0, cvm->nbOfUnderlyingTenors - 1);
  cvm->model_vols =
      dmatrix(0, cvm->nbOfOptionTenors - 1, 0, cvm->nbOfUnderlyingTenors - 1);
  for (i = 0; i < cvm->nbOfOptionTenors; ++i) {
    for (j = 0; j < cvm->nbOfUnderlyingTenors; ++j) {
      cvm->market_vols[i][j] = market_vols[i][j];
    }
  }

  err = cvm_calib_pillarsVols(cvm);
  //	err = cvm_calib_pillarsVolsAndCorrels(cvm);

  return err;
}

char *
SrtCalibCVMUnd(char *undName,
               // Market ID
               char *ycname, char *vcname,
               // Swaption Details
               SrtCompounding srtFreq, SrtBasisCode srtBasis, char *refRate,

               // Interpollation Type
               double shift, int converging_sliding, int interpollation_mode,

               // Pillars
               int nbOfPillars, char **pillars_tenors,

               // Pillars Correlations
               double ***pillars_correls,

               // Sensitivities
               int sensitype, double **lambdas,

               int nbOfSensiTenors, char **sensi_tenors, double **sensitivities,

               // Numerical Parameter
               int nQuadLegendre,

               // Calibration Parameters
               int nbOptionTenors, char **option_tenors, int nbUnderlyingTenors,
               char **underlying_tenors, double **market_vols,
               int allow_negativeVols) {
  Err err = NULL;
  int bCleanUpUndFlag = 1;
  SrtUndPtr pUndDesc;
  SrtUndListPtr und_list;
  SrtCurvePtr pYieldCurve;
  char *ccy;
  cvm_str *cvm = NULL;
  double *vol_dates = NULL;
  double **pillars_vols = NULL;
  int i, j, nbOfDates;
  long lToday, spotDate, end;

  cvm = calloc(1, sizeof(cvm_str));

  // Get the yield curve and the spot date from the yield curve name
  pYieldCurve = lookup_curve(ycname);
  if (!pYieldCurve) {
    err = "Cannot find Yield Curve";
    return err;
  }
  ccy = get_curve_ccy(pYieldCurve);
  lToday = (Date)get_today_from_curve(pYieldCurve);

  // Create the new CVM underlying
  und_list = get_underlying_list();
  pUndDesc = (SrtUndPtr)calloc(1, sizeof(SrtUndDesc));
  strcpy(pUndDesc->underl_name, undName);
  strupper(pUndDesc->underl_name);
  strip_white_space(pUndDesc->underl_name);
  strcpy(pUndDesc->underl_lbl, "CVM_UND");
  pUndDesc->underl_ccy = ccy;
  pUndDesc->underl_type = CVM_UND;

  nbOfDates = nbOptionTenors;
  vol_dates = dvector(0, nbOptionTenors - 1);
  pillars_vols = dmatrix(0, nbOptionTenors - 1, 0, nbOfPillars - 1);
  spotDate = add_unit(lToday, 2, SRT_BDAY, MODIFIED_SUCCEEDING);
  for (i = 0; i < nbOptionTenors; ++i) {
    err = add_tenor(spotDate, option_tenors[i], MODIFIED_SUCCEEDING, &end);
    end = add_unit(end, -2, SRT_BDAY, MODIFIED_SUCCEEDING);
    vol_dates[i] = end;
    for (j = 0; j < nbOfPillars; ++j) {
      pillars_vols[i][j] = market_vols[0][0];
    }
  }

  // Init CVM Structure
  err = cvm_fill_struct(
      // Market ID
      ycname, vcname,
      // Swaption Details
      srtFreq, srtBasis, refRate,

      // TS Dates
      shift, converging_sliding, interpollation_mode, nbOfDates, vol_dates,

      // Pillars
      nbOfPillars, pillars_tenors,

      // Pillars Correlations and Vols
      pillars_correls, pillars_vols,

      // Sensitivities
      sensitype, lambdas, nbOfSensiTenors, sensi_tenors, sensitivities,

      // Numerical Parameter
      nQuadLegendre,

      cvm);
  if (err)
    return err;

  cvm->allow_negativeVols = allow_negativeVols;

  err = cvm_calib_und(nbOptionTenors, option_tenors, nbUnderlyingTenors,
                      underlying_tenors, market_vols, cvm);
  if (err)
    return err;

  pUndDesc->spec_desc = cvm;

  // Put the underlying into the depot
  err = srt_f_lstins(und_list, pUndDesc->underl_name, 0.0, OBJ_PTR_UND,
                     (void *)pUndDesc, &cvm_free_und_struct,
                     &(pUndDesc->underl_ticker));

  if (vol_dates)
    free_dvector(vol_dates, 0, nbOptionTenors - 1);
  if (pillars_vols)
    free_dmatrix(pillars_vols, 0, nbOptionTenors - 1, 0, nbOfPillars - 1);

  return err;
}
