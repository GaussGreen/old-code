//---------------------------------------------------------------------
//-------------	Continuous Volatility Model	(CVM)	-------------------
//---------------------------------------------------------------------

#include "CVM.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

Err cvm_fill_struct(
    // Market ID
    char *ycid, char *vcid,
    // Swaption Details
    SrtCompounding srtFreq, SrtBasisCode srtBasis, char *refRate,

    // TS Dates
    double shift, int converging_sliding, int interpollation_mode,
    int nbOfDates, double *vol_dates,

    // Pillars
    int nbOfPillars, char **pillars_tenors,

    // Pillars Correlations and Vols
    double ***pillars_correls, double **pillars_vols_norm,

    // Sensitivities
    int sensitype,

    // Parametric Sensitivities
    double **lambdas,

    // Custom Sensitivities
    int nbOfSensiTenors, char **sensi_tenors, double **sensitivities,

    // Numerical Parameter
    int nQuadLegendre,

    cvm_str *cvm) {
  Err err = NULL;
  int i, j, k;
  Date lToday, end;
  SrtCurvePtr pCurve;

  cvm->shift = shift;

  cvm->sensi_type = sensitype;

  cvm->vol_dates = NULL;
  cvm->vol_times = NULL;

  cvm->pillars_tenors = NULL;
  cvm->pillars_maturities = NULL;
  cvm->pillars_times = NULL;

  cvm->pillars_vols = NULL;
  cvm->pillars_vols_norm = NULL;

  cvm->pillars_correls = NULL;
  cvm->pillars_angles = NULL;
  cvm->pillars_chol = NULL;

  cvm->lambdas = NULL;
  cvm->sensitivities = NULL;
  cvm->sensitivities2 = NULL;
  cvm->sensi_tenors = NULL;
  cvm->sensi_times = NULL;

  cvm->option_tenors = NULL;
  cvm->underlying_tenors = NULL;
  cvm->market_vols = NULL;
  cvm->model_vols = NULL;

  cvm->converging_sliding = converging_sliding;
  cvm->interpollation_mode = interpollation_mode;

  pCurve = lookup_curve(ycid);
  if (!pCurve) {
    err = "Cannot find Yield Curve";
    goto FREE_RETURN;
  }
  lToday = (Date)get_today_from_curve(pCurve);
  cvm->today = lToday;

  strcpy(cvm->ycid, ycid);
  strcpy(cvm->vcid, vcid);

  cvm->srtFreq = srtFreq;
  cvm->srtBasis = srtBasis;
  strcpy(cvm->refRate, refRate);
  err = srt_f_get_spot_lag_from_refrate(refRate, &cvm->spotLag);
  if (err)
    goto FREE_RETURN;
  cvm->spot = add_unit(lToday, cvm->spotLag, SRT_BDAY, MODIFIED_SUCCEEDING);

  cvm->nbOfDates = nbOfDates;
  cvm->vol_dates = dvector(0, nbOfDates - 1);
  cvm->vol_times = dvector(0, nbOfDates - 1);
  for (i = 0; i < nbOfDates; ++i) {
    cvm->vol_dates[i] = vol_dates[i];
    cvm->vol_times[i] = (vol_dates[i] - lToday) / 365.0;
  }

  cvm->nbOfPillars = nbOfPillars;
  cvm->pillars_tenors = (char256 *)calloc(nbOfPillars, sizeof(char256));
  cvm->pillars_maturities = dvector(0, nbOfPillars - 1);
  cvm->pillars_times = dvector(0, nbOfPillars - 1);
  for (i = 0; i < nbOfPillars; ++i) {
    strcpy(cvm->pillars_tenors[i], pillars_tenors[i]);
    err = add_tenor(cvm->spot, pillars_tenors[i], NO_BUSDAY_CONVENTION, &end);
    cvm->pillars_maturities[i] = end;
    if (err)
      goto FREE_RETURN;
    cvm->pillars_times[i] = (end - lToday) / 365.0;
  }

  cvm->pillars_vols = (double ***)calloc(nbOfDates, sizeof(double **));
  cvm->pillars_vols_norm = (double **)calloc(nbOfDates, sizeof(double *));
  cvm->pillars_correls = (double ***)calloc(nbOfDates, sizeof(double **));
  cvm->pillars_angles = (double ***)calloc(nbOfDates, sizeof(double **));
  cvm->pillars_chol = (double ***)calloc(nbOfDates, sizeof(double **));
  for (i = 0; i < nbOfDates; ++i) {
    cvm->pillars_vols_norm[i] = (double *)calloc(nbOfPillars, sizeof(double));
    cvm->pillars_vols[i] = (double **)calloc(nbOfPillars, sizeof(double *));
    cvm->pillars_correls[i] = (double **)calloc(nbOfPillars, sizeof(double *));
    cvm->pillars_angles[i] = (double **)calloc(nbOfPillars, sizeof(double *));
    cvm->pillars_chol[i] = (double **)calloc(nbOfPillars, sizeof(double *));
    for (j = 0; j < nbOfPillars; ++j) {
      cvm->pillars_vols_norm[i][j] = pillars_vols_norm[i][j];
      cvm->pillars_vols[i][j] = (double *)calloc(nbOfPillars, sizeof(double));
      cvm->pillars_correls[i][j] =
          (double *)calloc(nbOfPillars, sizeof(double));
      cvm->pillars_angles[i][j] = (double *)calloc(nbOfPillars, sizeof(double));
      cvm->pillars_chol[i][j] = (double *)calloc(nbOfPillars, sizeof(double));
      for (k = 0; k < nbOfPillars; ++k) {
        cvm->pillars_correls[i][j][k] = pillars_correls[i][j][k];
      }
    }

    PositiveMatrix(cvm->pillars_correls[i], nbOfPillars);

    CholDec(&CholAlg, cvm->pillars_correls[i], cvm->pillars_chol[i],
            nbOfPillars);

    err = cvm_set_pillars_vols(i, nbOfPillars, pillars_vols_norm[i], cvm);
    if (err)
      goto FREE_RETURN;
  }

  if ((sensitype == 0) || (sensitype == 2) ||
      ((sensitype == 1) && ((nbOfSensiTenors == 0) || (sensi_tenors == NULL) ||
                            (sensitivities == NULL)))) {
    cvm->nbOfSensiTenors = nbOfPillars;

    cvm->sensitivities = dmatrix(0, nbOfPillars - 1, 0, nbOfPillars - 1);
    cvm->sensitivities2 = dmatrix(0, nbOfPillars - 1, 0, nbOfPillars - 1);

    for (i = 0; i < nbOfPillars; ++i) {
      for (j = 0; j < nbOfPillars; ++j) {
        cvm->sensitivities[i][j] = 0.0;
      }
      cvm->sensitivities[i][i] = 1.0;
    }

    for (i = 0; i < nbOfPillars; ++i) {
      spline(cvm->pillars_times - 1, cvm->sensitivities[i] - 1, nbOfPillars,
             1.0e30, 1.0e30, cvm->sensitivities2[i] - 1);
    }

    cvm->lambdas = dmatrix(0, nbOfDates - 1, 0, nbOfPillars - 1);
    for (i = 0; i < nbOfDates; ++i) {
      for (j = 0; j < nbOfPillars; ++j) {
        cvm->lambdas[i][j] = lambdas[i][j];
      }
    }
  } else if (sensitype == 1) {
    cvm->nbOfSensiTenors = nbOfSensiTenors;

    cvm->sensitivities = dmatrix(0, nbOfPillars - 1, 0, nbOfSensiTenors - 1);
    cvm->sensitivities2 = dmatrix(0, nbOfPillars - 1, 0, nbOfSensiTenors - 1);
    cvm->sensi_tenors = (char256 *)calloc(nbOfSensiTenors, sizeof(char256));
    cvm->sensi_times = dvector(0, nbOfSensiTenors - 1);

    for (i = 0; i < nbOfSensiTenors; ++i) {
      for (j = 0; j < nbOfPillars; ++j) {
        cvm->sensitivities[j][i] = sensitivities[j][i];
      }
      strcpy(cvm->sensi_tenors[i], sensi_tenors[i]);
      err = add_tenor(cvm->spot, sensi_tenors[i], NO_BUSDAY_CONVENTION, &end);
      if (err)
        goto FREE_RETURN;
      cvm->sensi_times[i] = (end - lToday) / 365.0;
    }

    for (i = 0; i < nbOfPillars; ++i) {
      spline(cvm->sensi_times - 1, cvm->sensitivities[i] - 1, nbOfSensiTenors,
             1.0e30, 1.0e30, cvm->sensitivities2[i] - 1);
    }
  }

  cvm->nQuadLegendre = nQuadLegendre;

FREE_RETURN:

  if (err) {
    cvm_free_struct(cvm);
  }

  return err;
}

/*
Err cvm_fill_struct(
                                        //Market ID
                                        char			*ycid        ,
                                        char			*vcid        ,
                                        //Swaption Details
                                        SrtCompounding	srtFreq        ,
                                        SrtBasisCode	srtBasis        ,
                                        char			*refRate ,

                                        //TS Dates
                                        double			shift        ,
                                        int
converging_sliding        , int interpollation_mode        , int nbOfDates ,
double *vol_dates        ,

                                        //Pillars
                                        int nbOfPillars        , char
**pillars_tenors        ,

                                        //Pillars Correlations and Vols
                                        double ***pillars_correls        ,
                                        double **pillars_vols_norm        ,

                                        //Sensitivities
                                        int		sensitype        ,
                                        double	**lambdas        ,

                                        //Numerical Parameter
                                        int		nQuadLegendre        ,

                                        cvm_str *cvm)
{
        Err err = NULL;
        int i        ,j        ,k;
        Date lToday        , end;
        SrtCurvePtr pCurve;

        cvm->shift = shift;

        cvm->sensi_type = sensitype;

        cvm->vol_dates = NULL;
        cvm->vol_times = NULL;

        cvm->pillars_tenors = NULL;
        cvm->pillars_maturities = NULL;
        cvm->pillars_times = NULL;

        cvm->pillars_vols = NULL;
        cvm->pillars_vols_norm = NULL;

        cvm->pillars_correls = NULL;
        cvm->pillars_angles = NULL;
        cvm->pillars_chol = NULL;

        cvm->lambdas = NULL;
        cvm->sensitivities = NULL;
        cvm->sensitivities2 = NULL;

        cvm->option_tenors = NULL;
        cvm->underlying_tenors = NULL;
        cvm->market_vols = NULL;
        cvm->model_vols = NULL;

        cvm->converging_sliding = converging_sliding;
        cvm->interpollation_mode = interpollation_mode;

        pCurve = lookup_curve (ycid);
        if(!pCurve)
        {
                err = "Cannot find Yield Curve";
                goto FREE_RETURN;
        }
        lToday = (Date) get_today_from_curve (pCurve);
        cvm->today = lToday;

        strcpy(cvm->ycid        , ycid);
        strcpy(cvm->vcid        , vcid);

        cvm->srtFreq = srtFreq;
        cvm->srtBasis = srtBasis;
        strcpy(cvm->refRate        , refRate);
        err =  srt_f_get_spot_lag_from_refrate(refRate        , &cvm->spotLag);
        if (err) goto FREE_RETURN;
        cvm->spot = add_unit(lToday        , cvm->spotLag        , SRT_BDAY ,
MODIFIED_SUCCEEDING);

        cvm->nbOfDates = nbOfDates;
        cvm->vol_dates = dvector(0        , nbOfDates-1);
        cvm->vol_times = dvector(0        , nbOfDates-1);
        for(i=0;i<nbOfDates;++i)
        {
                cvm->vol_dates[i] = vol_dates[i];
                cvm->vol_times[i] = (vol_dates[i] - lToday) / 365.0;
        }

        cvm->nbOfPillars = nbOfPillars;
        cvm->pillars_tenors = (char256*) calloc(nbOfPillars        ,
sizeof(char256)); cvm->pillars_maturities = dvector(0        , nbOfPillars-1);
        cvm->pillars_times = dvector(0        , nbOfPillars-1);
        for(i=0;i<nbOfPillars;++i)
        {
                strcpy(cvm->pillars_tenors[i]        , pillars_tenors[i]);
                err = add_tenor(cvm->spot        , pillars_tenors[i]        ,
NO_BUSDAY_CONVENTION        , &end); cvm->pillars_maturities[i] = end; if(err)
goto FREE_RETURN; cvm->pillars_times[i] = (end - lToday) / 365.0;
        }

        cvm->pillars_vols = (double***) calloc(nbOfDates        ,
sizeof(double**)); cvm->pillars_vols_norm = (double**) calloc(nbOfDates        ,
sizeof(double*)); cvm->pillars_correls = (double***) calloc(nbOfDates        ,
sizeof(double**)); cvm->pillars_angles = (double***) calloc(nbOfDates        ,
sizeof(double**)); cvm->pillars_chol = (double***) calloc(nbOfDates        ,
sizeof(double**)); for(i=0;i<nbOfDates;++i)
        {
                cvm->pillars_vols_norm[i] = (double*) calloc(nbOfPillars ,
sizeof(double)); cvm->pillars_vols[i] = (double**) calloc(nbOfPillars        ,
sizeof(double*)); cvm->pillars_correls[i] = (double**) calloc(nbOfPillars ,
sizeof(double*)); cvm->pillars_angles[i] = (double**) calloc(nbOfPillars ,
sizeof(double*)); cvm->pillars_chol[i] = (double**) calloc(nbOfPillars        ,
sizeof(double*)); for(j=0;j<nbOfPillars;++j)
                {
                        cvm->pillars_vols_norm[i][j] = pillars_vols_norm[i][j];
                        cvm->pillars_vols[i][j] = (double*) calloc(nbOfPillars
      , sizeof(double)); cvm->pillars_correls[i][j] = (double*)
calloc(nbOfPillars        , sizeof(double)); cvm->pillars_angles[i][j] =
(double*) calloc(nbOfPillars        , sizeof(double)); cvm->pillars_chol[i][j] =
(double*) calloc(nbOfPillars        , sizeof(double));
for(k=0;k<nbOfPillars;++k)
                        {
                                cvm->pillars_correls[i][j][k] =
pillars_correls[i][j][k];
                        }
                }

                CholDec	(&CholAlg        ,
                                 cvm->pillars_correls[i]        ,
                                 cvm->pillars_chol[i]        ,
                                 nbOfPillars);
        }

        cvm->lambdas = dmatrix(0        , nbOfDates-1        , 0        ,
nbOfPillars-1); cvm->sensitivities = dmatrix(0        , nbOfPillars-1        , 0
, nbOfPillars-1); cvm->sensitivities2 = dmatrix(0        , nbOfPillars-1 , 0 ,
nbOfPillars-1);

        for(i=0;i<nbOfPillars;++i)
        {
                for(j=0;j<nbOfPillars;++j)
                {
                        cvm->sensitivities[i][j]=0.0;
                }
                cvm->sensitivities[i][i]=1.0;
        }

        for(i=0;i<nbOfPillars;++i)
        {
                spline(cvm->pillars_times-1        , cvm->sensitivities[i]-1 ,
nbOfPillars        , 1.0e30        , 1.0e30        , cvm->sensitivities2[i]-1);
        }

        cvm->nQuadLegendre = nQuadLegendre;

        for(i=0;i<nbOfDates;++i)
        {
                for(j=0;j<nbOfPillars;++j)
                {
                        cvm->lambdas[i][j] = lambdas[i][j];
                }
                err = cvm_set_pillars_vols(i        , nbOfPillars        ,
pillars_vols_norm[i]        , cvm); if(err) goto FREE_RETURN;
        }

FREE_RETURN:

        if(err)
        {
                cvm_free_struct(cvm);
        }

        return err;
}
*/

Err cvm_copy_struct(cvm_str *cvm_target, cvm_str *cvm_source) {
  Err err = NULL;
  char **pillars_tenors = NULL;
  char **sensi_tenors = NULL;
  int i;

  pillars_tenors = svector_size(0, cvm_source->nbOfPillars - 1, 256);
  sensi_tenors = svector_size(0, cvm_source->nbOfSensiTenors - 1, 256);

  cvm_target->allow_negativeVols = cvm_source->allow_negativeVols;

  for (i = 0; i < cvm_source->nbOfPillars; ++i) {
    strcpy(pillars_tenors[i], cvm_source->pillars_tenors[i]);
  }

  if (cvm_source->sensi_type == 1) {
    for (i = 0; i < cvm_source->nbOfSensiTenors; ++i) {
      strcpy(sensi_tenors[i], cvm_source->sensi_tenors[i]);
    }
  }

  err = cvm_fill_struct(
      cvm_source->ycid, cvm_source->vcid, cvm_source->srtFreq,
      cvm_source->srtBasis, cvm_source->refRate, cvm_source->shift,
      cvm_source->converging_sliding, cvm_source->interpollation_mode,
      cvm_source->nbOfDates, cvm_source->vol_dates, cvm_source->nbOfPillars,
      pillars_tenors, cvm_source->pillars_correls,
      cvm_source->pillars_vols_norm, cvm_source->sensi_type,
      cvm_source->lambdas, cvm_source->nbOfSensiTenors, sensi_tenors,
      cvm_source->sensitivities, cvm_source->nQuadLegendre, cvm_target);

  if (pillars_tenors)
    free_svector_size(pillars_tenors, 0, cvm_source->nbOfPillars - 1, 256);
  if (sensi_tenors)
    free_svector_size(sensi_tenors, 0, cvm_source->nbOfSensiTenors - 1, 256);

  return err;
}

Err cvm_set_pillars_vols(int dateIndex, int nPillars, double *pillarsVols,
                         cvm_str *cvm) {
  Err err = NULL;
  int i, k;

  for (i = 0; i < nPillars; ++i) {
    cvm->pillars_vols_norm[dateIndex][i] = pillarsVols[i];
    for (k = 0; k < nPillars; ++k) {
      cvm->pillars_vols[dateIndex][k][i] =
          cvm->pillars_chol[dateIndex][k][i] * pillarsVols[k];
    }
  }

  return err;
}

Err cvm_bump_pillars_correls(int dateIndex, int nPillarsMinusOne,
                             double *correls, cvm_str *cvm) {
  Err err = NULL;
  int i, k;

  for (i = 0; i < nPillarsMinusOne; ++i) {
    for (k = 1; k < nPillarsMinusOne - i; ++k) {
      cvm->pillars_correls[dateIndex][k][i + k] =
          DMAX(-0.9999, DMIN(0.9999, cvm->pillars_correls[dateIndex][k][i + k] +
                                         correls[i]));
    }
  }

  PositiveMatrix(cvm->pillars_correls[dateIndex], cvm->nbOfPillars);

  CholDec(&CholAlg, cvm->pillars_correls[dateIndex],
          cvm->pillars_chol[dateIndex], cvm->nbOfPillars);

  for (i = 0; i < cvm->nbOfPillars; ++i) {
    for (k = 0; k < cvm->nbOfPillars; ++k) {
      cvm->pillars_vols[dateIndex][k][i] = cvm->pillars_chol[dateIndex][k][i] *
                                           cvm->pillars_vols_norm[dateIndex][k];
    }
  }

  return err;
}

Err cvm_set_pillars_correls(int dateIndex, int nPillars, double **correls,
                            cvm_str *cvm) {
  Err err = NULL;
  int i, k;

  for (i = 0; i < nPillars; ++i) {
    for (k = 0; k < nPillars; ++k) {
      cvm->pillars_correls[dateIndex][i][k] = correls[i][k];
    }
  }

  PositiveMatrix(cvm->pillars_correls[dateIndex], cvm->nbOfPillars);

  CholDec(&CholAlg, cvm->pillars_correls[dateIndex],
          cvm->pillars_chol[dateIndex], cvm->nbOfPillars);

  for (i = 0; i < cvm->nbOfPillars; ++i) {
    for (k = 0; k < cvm->nbOfPillars; ++k) {
      cvm->pillars_vols[dateIndex][k][i] = cvm->pillars_chol[dateIndex][k][i] *
                                           cvm->pillars_vols_norm[dateIndex][k];
    }
  }

  return err;
}

Err cvm_compute_fwdvol_at_voldate(int dateIndex, double maturity, int nPillars,
                                  double *fwdvol, cvm_str *cvm) {
  Err err = NULL;
  int i, k;
  double *sensi = NULL;
  double date, time;

  date = cvm->vol_dates[dateIndex];
  time = cvm->vol_times[dateIndex];
  if (time > maturity) {
    if (dateIndex > 0) {
      err = cvm_compute_fwdvol_at_voldate(dateIndex - 1, maturity, nPillars,
                                          fwdvol, cvm);
      return err;
    } else {
      err = "This case should not occur !!!";
      return err;
    }
  }

  sensi = dvector(0, cvm->nbOfPillars - 1);
  err = cvm_compute_sensitivities(time, maturity, cvm->nbOfPillars, sensi, cvm);
  if (err)
    goto FREE_RETURN;

  for (i = 0; i < nPillars; ++i) {
    fwdvol[i] = 0.0;
  }

  for (i = 0; i < nPillars; ++i) {
    for (k = 0; k < nPillars; ++k) {
      fwdvol[i] += cvm->pillars_vols[dateIndex][k][i] * sensi[k];
    }
  }

FREE_RETURN:
  if (sensi)
    free_dvector(sensi, 0, cvm->nbOfPillars - 1);

  return err;
}

Err cvm_compute_pillarsvols(double time, double maturity, int nPillars,
                            double **pillarsvols, cvm_str *cvm) {
  Err err = NULL;
  int i, j, k;

  i = 0;
  if (time <= cvm->vol_times[0]) {
    for (j = 0; j < nPillars; ++j) {
      for (k = 0; k < nPillars; ++k) {
        pillarsvols[j][k] = cvm->pillars_vols[0][j][k];
      }
    }
  } else {
    while ((i < cvm->nbOfDates) && (time > cvm->vol_times[i])) {
      i = i + 1;
    }

    if (i == cvm->nbOfDates) {
      for (j = 0; j < nPillars; ++j) {
        for (k = 0; k < nPillars; ++k) {
          pillarsvols[j][k] = cvm->pillars_vols[cvm->nbOfDates - 1][j][k];
        }
      }
    } else {
      for (j = 0; j < nPillars; ++j) {
        for (k = 0; k < nPillars; ++k) {
          if (cvm->interpollation_mode == 0) {
            pillarsvols[j][k] = cvm->pillars_vols[i][j][k];
          } else {
            pillarsvols[j][k] =
                cvm->pillars_vols[i - 1][j][k] +
                (cvm->pillars_vols[i][j][k] - cvm->pillars_vols[i - 1][j][k]) *
                    (time - cvm->vol_times[i - 1]) /
                    (cvm->vol_times[i] - cvm->vol_times[i - 1]);
          }
        }
      }
    }
  }

  return err;
}

Err cvm_compute_fwdvol_sliding(double time, double maturity, int nPillars,
                               double *fwdvol, cvm_str *cvm) {
  Err err = NULL;
  double **vol = NULL;
  double *sensi = NULL;
  int i, j;

  vol = dmatrix(0, cvm->nbOfPillars - 1, 0, cvm->nbOfPillars - 1);
  sensi = dvector(0, cvm->nbOfPillars - 1);

  err = cvm_compute_pillarsvols(time, maturity, cvm->nbOfPillars, vol, cvm);
  if (err)
    goto FREE_RETURN;
  err = cvm_compute_sensitivities(time, maturity, cvm->nbOfPillars, sensi, cvm);
  if (err)
    goto FREE_RETURN;

  for (i = 0; i < cvm->nbOfPillars; ++i) {
    fwdvol[i] = 0.0;
  }
  for (i = 0; i < cvm->nbOfPillars; ++i) {
    for (j = 0; j < cvm->nbOfPillars; ++j) {
      fwdvol[i] += sensi[j] * vol[j][i];
    }
  }

FREE_RETURN:

  if (vol)
    free_dmatrix(vol, 0, cvm->nbOfPillars - 1, 0, cvm->nbOfPillars - 1);
  if (sensi)
    free_dvector(sensi, 0, cvm->nbOfPillars - 1);

  return err;
}

Err cvm_compute_fwdvol_converging(double time, double maturity, int nPillars,
                                  double *fwdvol, cvm_str *cvm) {
  Err err = NULL;
  double *vol1 = NULL;
  double *vol2 = NULL;
  int i, k;

  if (maturity < cvm->vol_times[0]) {
    err = cvm_compute_fwdvol_converging(time + cvm->vol_times[0] - maturity,
                                        cvm->vol_times[0], cvm->nbOfPillars,
                                        fwdvol, cvm);
    return err;
  }

  vol1 = dvector(0, cvm->nbOfPillars - 1);
  vol2 = dvector(0, cvm->nbOfPillars - 1);

  i = 0;
  if (time <= cvm->vol_times[0]) {
    err = cvm_compute_fwdvol_at_voldate(0, maturity, nPillars, vol1, cvm);
    if (err)
      goto FREE_RETURN;
    for (k = 0; k < nPillars; ++k) {
      fwdvol[k] = vol1[k];
    }
  } else {
    while ((i < cvm->nbOfDates) && (time > cvm->vol_times[i])) {
      i = i + 1;
    }

    if (i == cvm->nbOfDates) {
      err = cvm_compute_fwdvol_at_voldate(cvm->nbOfDates - 1, maturity,
                                          nPillars, vol2, cvm);
      if (err)
        goto FREE_RETURN;
      for (k = 0; k < nPillars; ++k) {
        fwdvol[k] = vol2[k];
      }
    } else {
      err = cvm_compute_fwdvol_at_voldate(i - 1, maturity, nPillars, vol1, cvm);
      if (err)
        goto FREE_RETURN;
      err = cvm_compute_fwdvol_at_voldate(i, maturity, nPillars, vol2, cvm);
      if (err)
        goto FREE_RETURN;
      for (k = 0; k < nPillars; ++k) {
        if (cvm->interpollation_mode == 0) {
          fwdvol[k] = vol2[k]; // vol1[k] + (vol2[k] - vol1[k]) * (time -
                               // cvm->vol_times[i-1]) / (cvm->vol_times[i] -
                               // cvm->vol_times[i-1]);
        } else {
          fwdvol[k] = vol1[k] + (vol2[k] - vol1[k]) *
                                    (time - cvm->vol_times[i - 1]) /
                                    (cvm->vol_times[i] - cvm->vol_times[i - 1]);
        }
      }
    }
  }

FREE_RETURN:

  if (vol1)
    free_dvector(vol1, 0, cvm->nbOfPillars - 1);
  if (vol2)
    free_dvector(vol2, 0, cvm->nbOfPillars - 1);

  return err;
}

Err cvm_compute_fwdvol(double time, double maturity, int nPillars,
                       double *fwdvol, cvm_str *cvm) {
  Err err = NULL;

  if (cvm->converging_sliding == 0) {
    err = cvm_compute_fwdvol_converging(time, maturity, nPillars, fwdvol, cvm);
  } else {
    err = cvm_compute_fwdvol_sliding(time, maturity, nPillars, fwdvol, cvm);
  }

  return err;
}

Err cvm_vol(double start, double end, double *vol, cvm_str *cvm) {
  Err err = NULL;
  int i, j, k, l, dateIndex;
  SwapDP Swap;
  Date lToday;

  long exercise;
  double exercise_time;

  long iNFixPayDates, iNFixDates;
  long *lFixPayDates = NULL, *lFixStartDates = NULL, *lFixEndDates = NULL;
  double *dFixCoverages = NULL;

  long iNFloatPayDates, iNFloatDates;
  long *lFloatFixingDates = NULL, *lFloatPayDates = NULL,
       *lFloatStartDates = NULL, *lFloatEndDates = NULL;
  double *dFloatCoverages = NULL;
  double *dFloatSpreads = NULL;

  double level, fwdswap, frai;
  double *weights = NULL;
  double *volvect = NULL;
  double *volloc = NULL;
  double volscal;
  double spreadleg;

  double *x = NULL;
  double *w = NULL;

  lToday = cvm->today;

  err = swp_f_setSwapDP((long)(start), (long)(end), cvm->srtFreq, cvm->srtBasis,
                        &Swap);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap, lToday, &lFixPayDates, &iNFixPayDates, &lFixStartDates,
      &lFixEndDates, &dFixCoverages, &iNFixDates);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap, lToday, cvm->refRate, &lFloatPayDates, &iNFloatPayDates,
      &lFloatFixingDates, &lFloatStartDates, &lFloatEndDates, &dFloatCoverages,
      &dFloatSpreads, &iNFloatDates);
  if (err)
    goto FREE_RETURN;

  level = 0.0;
  for (i = 0; i < iNFixDates; ++i) {
    level +=
        dFixCoverages[i] * swp_f_df(lToday, lFixPayDates[i + 1], cvm->ycid);
  }

  weights = dvector(0, iNFloatDates - 1);
  if (cvm->shift != -1) {
    spreadleg = 0;
  }
  for (i = 0; i < iNFloatDates; ++i) {
    if (cvm->shift == -1) {
      weights[i] = dFloatCoverages[i] *
                   swp_f_df(lToday, lFloatPayDates[i + 1], cvm->ycid) / level;
    } else {
      frai = (swp_f_df(lToday, lFloatPayDates[i], cvm->ycid) /
                  swp_f_df(lToday, lFloatPayDates[i + 1], cvm->ycid) -
              1) /
                 dFloatCoverages[i] +
             dFloatSpreads[i];

      spreadleg += dFloatCoverages[i] * dFloatSpreads[i] *
                   swp_f_df(lToday, lFloatPayDates[i + 1], cvm->ycid);

      weights[i] = (cvm->shift + frai) * dFloatCoverages[i] *
                   swp_f_df(lToday, lFloatPayDates[i + 1], cvm->ycid) / level;
    }
  }

  if (cvm->shift != -1) {
    fwdswap = (swp_f_df(lToday, lFloatPayDates[0], cvm->ycid) -
               swp_f_df(lToday, lFloatPayDates[iNFloatDates], cvm->ycid) +
               spreadleg) /
              level;
  }

  volvect = dvector(0, cvm->nbOfPillars - 1);
  volloc = dvector(0, cvm->nbOfPillars - 1);

  exercise =
      add_unit((long)(start), -cvm->spotLag, SRT_BDAY, MODIFIED_SUCCEEDING);
  exercise_time = (exercise - lToday) / 365.0;

  x = dvector(0, cvm->nQuadLegendre - 1);
  w = dvector(0, cvm->nQuadLegendre - 1);
  dateIndex = 0;
  volscal = 0.0;
  if (exercise <= cvm->vol_dates[0]) {
    GaussLeg(0, exercise_time, x - 1, w - 1, cvm->nQuadLegendre);
    for (k = 0; k < cvm->nQuadLegendre; ++k) {
      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volvect[j] = 0.0;
      }

      for (i = 0; i < iNFloatDates; ++i) {
        err = cvm_compute_fwdvol(x[k], (lFloatStartDates[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc, cvm);
        if (err)
          goto FREE_RETURN;
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          if (cvm->shift == -1) {
            volvect[j] += weights[i] * volloc[j];
          } else {
            volvect[j] += weights[i] * volloc[j] / (cvm->shift + fwdswap);
          }
        }
      }

      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volscal += w[k] * volvect[j] * volvect[j];
      }
    }
  } else {
    while ((dateIndex < cvm->nbOfDates) &&
           (exercise > cvm->vol_dates[dateIndex])) {
      dateIndex = dateIndex + 1;
    }

    GaussLeg(0, cvm->vol_times[0], x - 1, w - 1, cvm->nQuadLegendre);
    for (k = 0; k < cvm->nQuadLegendre; ++k) {
      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volvect[j] = 0.0;
      }
      for (i = 0; i < iNFloatDates; ++i) {
        err = cvm_compute_fwdvol(x[k], (lFloatStartDates[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc, cvm);
        if (err)
          goto FREE_RETURN;
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          if (cvm->shift == -1) {
            volvect[j] += weights[i] * volloc[j];
          } else {
            volvect[j] += weights[i] * volloc[j] / (cvm->shift + fwdswap);
          }
        }
      }

      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volscal += w[k] * volvect[j] * volvect[j];
      }
    }

    for (l = 1; l < dateIndex; ++l) {
      GaussLeg(cvm->vol_times[l - 1], cvm->vol_times[l], x - 1, w - 1,
               cvm->nQuadLegendre);
      for (k = 0; k < cvm->nQuadLegendre; ++k) {
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect[j] = 0.0;
        }
        for (i = 0; i < iNFloatDates; ++i) {
          err = cvm_compute_fwdvol(x[k], (lFloatStartDates[i] - lToday) / 365.0,
                                   cvm->nbOfPillars, volloc, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            if (cvm->shift == -1) {
              volvect[j] += weights[i] * volloc[j];
            } else {
              volvect[j] += weights[i] * volloc[j] / (cvm->shift + fwdswap);
            }
          }
        }

        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volscal += w[k] * volvect[j] * volvect[j];
        }
      }
    }

    GaussLeg(cvm->vol_times[dateIndex - 1], exercise_time, x - 1, w - 1,
             cvm->nQuadLegendre);
    for (k = 0; k < cvm->nQuadLegendre; ++k) {
      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volvect[j] = 0.0;
      }
      for (i = 0; i < iNFloatDates; ++i) {
        err = cvm_compute_fwdvol(x[k], (lFloatStartDates[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc, cvm);
        if (err)
          goto FREE_RETURN;
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          if (cvm->shift == -1) {
            volvect[j] += weights[i] * volloc[j];
          } else {
            volvect[j] += weights[i] * volloc[j] / (cvm->shift + fwdswap);
          }
        }
      }

      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volscal += w[k] * volvect[j] * volvect[j];
      }
    }
  }

  *vol = sqrt(volscal / exercise_time);

FREE_RETURN:

  if (lFixPayDates)
    free(lFixPayDates);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixCoverages)
    free(dFixCoverages);

  if (lFloatPayDates)
    free(lFloatPayDates);
  if (lFloatFixingDates)
    free(lFloatFixingDates);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dFloatCoverages)
    free(dFloatCoverages);
  if (dFloatSpreads)
    free(dFloatSpreads);

  if (weights)
    free_dvector(weights, 0, iNFloatDates - 1);

  if (volvect)
    free_dvector(volvect, 0, cvm->nbOfPillars - 1);
  if (volloc)
    free_dvector(volloc, 0, cvm->nbOfPillars - 1);

  if (x)
    free_dvector(x, 0, cvm->nQuadLegendre - 1);
  if (w)
    free_dvector(w, 0, cvm->nQuadLegendre - 1);

  return err;
}

Err cvm_partial_forward_vol(double date, double expiry, double start,
                            double end, double *vol, cvm_str *cvm) {
  Err err = NULL;
  int i, j, k, l, dateIndex, dateIndex2;
  SwapDP Swap;
  Date lToday;

  long exercise;
  double exercise_time;

  double time, expiry_time;

  long iNFixPayDates, iNFixDates;
  long *lFixPayDates = NULL, *lFixStartDates = NULL, *lFixEndDates = NULL;
  double *dFixCoverages = NULL;

  long iNFloatPayDates, iNFloatDates;
  long *lFloatFixingDates = NULL, *lFloatPayDates = NULL,
       *lFloatStartDates = NULL, *lFloatEndDates = NULL;
  double *dFloatCoverages = NULL;
  double *dFloatSpreads = NULL;

  double level;
  double *weights = NULL;
  double *volvect = NULL;
  double *volloc = NULL;
  double volscal;

  double *x = NULL;
  double *w = NULL;

  lToday = cvm->today;

  time = (date - lToday) / 365.0;
  expiry_time = (expiry - lToday) / 365.0;

  err = swp_f_setSwapDP((long)(start), (long)(end), cvm->srtFreq, cvm->srtBasis,
                        &Swap);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap, lToday, &lFixPayDates, &iNFixPayDates, &lFixStartDates,
      &lFixEndDates, &dFixCoverages, &iNFixDates);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap, lToday, cvm->refRate, &lFloatPayDates, &iNFloatPayDates,
      &lFloatFixingDates, &lFloatStartDates, &lFloatEndDates, &dFloatCoverages,
      &dFloatSpreads, &iNFloatDates);
  if (err)
    goto FREE_RETURN;

  level = 0.0;
  for (i = 0; i < iNFixDates; ++i) {
    level +=
        dFixCoverages[i] * swp_f_df(lToday, lFixPayDates[i + 1], cvm->ycid);
  }

  weights = dvector(0, iNFloatDates - 1);
  for (i = 0; i < iNFloatDates; ++i) {
    weights[i] = dFloatCoverages[i] *
                 swp_f_df(lToday, lFloatPayDates[i + 1], cvm->ycid) / level;
  }

  volvect = dvector(0, cvm->nbOfPillars - 1);
  volloc = dvector(0, cvm->nbOfPillars - 1);

  exercise =
      add_unit((long)(start), -cvm->spotLag, SRT_BDAY, MODIFIED_SUCCEEDING);
  exercise_time = (exercise - lToday) / 365.0;

  exercise = (long)(DMIN(exercise, expiry));
  exercise_time = DMIN(exercise_time, expiry_time);

  x = dvector(0, cvm->nQuadLegendre - 1);
  w = dvector(0, cvm->nQuadLegendre - 1);
  dateIndex = 0;
  dateIndex2 = 0;
  volscal = 0.0;
  if (exercise <= cvm->vol_dates[0]) {
    GaussLeg(time, exercise_time, x - 1, w - 1, cvm->nQuadLegendre);
    for (k = 0; k < cvm->nQuadLegendre; ++k) {
      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volvect[j] = 0.0;
      }

      for (i = 0; i < iNFloatDates; ++i) {
        err = cvm_compute_fwdvol(x[k], (lFloatStartDates[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc, cvm);
        if (err)
          goto FREE_RETURN;
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect[j] += weights[i] * volloc[j];
        }
      }

      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volscal += w[k] * volvect[j] * volvect[j];
      }
    }
  } else {
    while ((dateIndex < cvm->nbOfDates) &&
           (exercise > cvm->vol_dates[dateIndex])) {
      dateIndex = dateIndex + 1;
    }

    while ((dateIndex2 < cvm->nbOfDates) &&
           (date > cvm->vol_dates[dateIndex2])) {
      dateIndex2 = dateIndex2 + 1;
    }

    if ((dateIndex2 < cvm->nbOfDates) && (dateIndex2 < dateIndex)) {
      GaussLeg(time, cvm->vol_times[dateIndex2], x - 1, w - 1,
               cvm->nQuadLegendre);
      for (k = 0; k < cvm->nQuadLegendre; ++k) {
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect[j] = 0.0;
        }
        for (i = 0; i < iNFloatDates; ++i) {
          err = cvm_compute_fwdvol(x[k], (lFloatStartDates[i] - lToday) / 365.0,
                                   cvm->nbOfPillars, volloc, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect[j] += weights[i] * volloc[j];
          }
        }

        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volscal += w[k] * volvect[j] * volvect[j];
        }
      }

      for (l = dateIndex2 + 1; l < dateIndex; ++l) {
        GaussLeg(cvm->vol_times[l - 1], cvm->vol_times[l], x - 1, w - 1,
                 cvm->nQuadLegendre);
        for (k = 0; k < cvm->nQuadLegendre; ++k) {
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect[j] = 0.0;
          }
          for (i = 0; i < iNFloatDates; ++i) {
            err =
                cvm_compute_fwdvol(x[k], (lFloatStartDates[i] - lToday) / 365.0,
                                   cvm->nbOfPillars, volloc, cvm);
            if (err)
              goto FREE_RETURN;
            for (j = 0; j < cvm->nbOfPillars; ++j) {
              volvect[j] += weights[i] * volloc[j];
            }
          }

          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volscal += w[k] * volvect[j] * volvect[j];
          }
        }
      }

      GaussLeg(cvm->vol_times[dateIndex - 1], exercise_time, x - 1, w - 1,
               cvm->nQuadLegendre);
      for (k = 0; k < cvm->nQuadLegendre; ++k) {
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect[j] = 0.0;
        }
        for (i = 0; i < iNFloatDates; ++i) {
          err = cvm_compute_fwdvol(x[k], (lFloatStartDates[i] - lToday) / 365.0,
                                   cvm->nbOfPillars, volloc, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect[j] += weights[i] * volloc[j];
          }
        }

        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volscal += w[k] * volvect[j] * volvect[j];
        }
      }
    } else {
      GaussLeg(time, exercise_time, x - 1, w - 1, cvm->nQuadLegendre);
      for (k = 0; k < cvm->nQuadLegendre; ++k) {
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect[j] = 0.0;
        }
        for (i = 0; i < iNFloatDates; ++i) {
          err = cvm_compute_fwdvol(x[k], (lFloatStartDates[i] - lToday) / 365.0,
                                   cvm->nbOfPillars, volloc, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect[j] += weights[i] * volloc[j];
          }
        }

        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volscal += w[k] * volvect[j] * volvect[j];
        }
      }
    }
  }

  *vol = sqrt(volscal / (exercise_time - time));

FREE_RETURN:

  if (lFixPayDates)
    free(lFixPayDates);
  if (lFixStartDates)
    free(lFixStartDates);
  if (lFixEndDates)
    free(lFixEndDates);
  if (dFixCoverages)
    free(dFixCoverages);

  if (lFloatPayDates)
    free(lFloatPayDates);
  if (lFloatFixingDates)
    free(lFloatFixingDates);
  if (lFloatStartDates)
    free(lFloatStartDates);
  if (lFloatEndDates)
    free(lFloatEndDates);
  if (dFloatCoverages)
    free(dFloatCoverages);
  if (dFloatSpreads)
    free(dFloatSpreads);

  if (weights)
    free_dvector(weights, 0, iNFloatDates - 1);

  if (volvect)
    free_dvector(volvect, 0, cvm->nbOfPillars - 1);
  if (volloc)
    free_dvector(volloc, 0, cvm->nbOfPillars - 1);

  if (x)
    free_dvector(x, 0, cvm->nQuadLegendre - 1);
  if (w)
    free_dvector(w, 0, cvm->nQuadLegendre - 1);

  return err;
}

Err cvm_partial_forward_cov(double date, double expiry, double start1,
                            double end1, double start2, double end2,
                            double *vol1, double *vol2, double *cov,
                            cvm_str *cvm) {
  Err err = NULL;
  int i, j, k, l, dateIndex, dateIndex2;
  SwapDP Swap1, Swap2;
  Date lToday;

  long exercise;
  double exercise_time;

  long exercise1;
  double exercise_time1;

  long exercise2;
  double exercise_time2;

  double time, expiry_time;

  long iNFixPayDates1, iNFixDates1;
  long *lFixPayDates1 = NULL, *lFixStartDates1 = NULL, *lFixEndDates1 = NULL;
  double *dFixCoverages1 = NULL;

  long iNFloatPayDates1, iNFloatDates1;
  long *lFloatFixingDates1 = NULL, *lFloatPayDates1 = NULL,
       *lFloatStartDates1 = NULL, *lFloatEndDates1 = NULL;
  double *dFloatCoverages1 = NULL;
  double *dFloatSpreads1 = NULL;

  long iNFixPayDates2, iNFixDates2;
  long *lFixPayDates2 = NULL, *lFixStartDates2 = NULL, *lFixEndDates2 = NULL;
  double *dFixCoverages2 = NULL;

  long iNFloatPayDates2, iNFloatDates2;
  long *lFloatFixingDates2 = NULL, *lFloatPayDates2 = NULL,
       *lFloatStartDates2 = NULL, *lFloatEndDates2 = NULL;
  double *dFloatCoverages2 = NULL;
  double *dFloatSpreads2 = NULL;

  double level1;
  double *weights1 = NULL;
  double *volvect1 = NULL;
  double *volloc1 = NULL;
  double volscal1;

  double level2;
  double *weights2 = NULL;
  double *volvect2 = NULL;
  double *volloc2 = NULL;
  double volscal2;

  double covscal;

  double *x = NULL;
  double *w = NULL;

  lToday = cvm->today;

  time = (date - lToday) / 365.0;
  expiry_time = (expiry - lToday) / 365.0;

  // Underlying 1
  err = swp_f_setSwapDP((long)(start1), (long)(end1), cvm->srtFreq,
                        cvm->srtBasis, &Swap1);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap1, lToday, &lFixPayDates1, &iNFixPayDates1, &lFixStartDates1,
      &lFixEndDates1, &dFixCoverages1, &iNFixDates1);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap1, lToday, cvm->refRate, &lFloatPayDates1, &iNFloatPayDates1,
      &lFloatFixingDates1, &lFloatStartDates1, &lFloatEndDates1,
      &dFloatCoverages1, &dFloatSpreads1, &iNFloatDates1);
  if (err)
    goto FREE_RETURN;

  level1 = 0.0;
  for (i = 0; i < iNFixDates1; ++i) {
    level1 +=
        dFixCoverages1[i] * swp_f_df(lToday, lFixPayDates1[i + 1], cvm->ycid);
  }

  weights1 = dvector(0, iNFloatDates1 - 1);
  for (i = 0; i < iNFloatDates1; ++i) {
    weights1[i] = dFloatCoverages1[i] *
                  swp_f_df(lToday, lFloatPayDates1[i + 1], cvm->ycid) / level1;
  }

  volvect1 = dvector(0, cvm->nbOfPillars - 1);
  volloc1 = dvector(0, cvm->nbOfPillars - 1);

  exercise1 =
      add_unit((long)(start1), -cvm->spotLag, SRT_BDAY, MODIFIED_SUCCEEDING);
  exercise_time1 = (exercise1 - lToday) / 365.0;

  // Underlying 2
  err = swp_f_setSwapDP((long)(start2), (long)(end2), cvm->srtFreq,
                        cvm->srtBasis, &Swap2);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FixedLegDatesAndCoverages(
      &Swap2, lToday, &lFixPayDates2, &iNFixPayDates2, &lFixStartDates2,
      &lFixEndDates2, &dFixCoverages2, &iNFixDates2);
  if (err)
    goto FREE_RETURN;

  err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
      &Swap2, lToday, cvm->refRate, &lFloatPayDates2, &iNFloatPayDates2,
      &lFloatFixingDates2, &lFloatStartDates2, &lFloatEndDates2,
      &dFloatCoverages2, &dFloatSpreads2, &iNFloatDates2);
  if (err)
    goto FREE_RETURN;

  level2 = 0.0;
  for (i = 0; i < iNFixDates2; ++i) {
    level2 +=
        dFixCoverages2[i] * swp_f_df(lToday, lFixPayDates2[i + 1], cvm->ycid);
  }

  weights2 = dvector(0, iNFloatDates2 - 1);
  for (i = 0; i < iNFloatDates2; ++i) {
    weights2[i] = dFloatCoverages2[i] *
                  swp_f_df(lToday, lFloatPayDates2[i + 1], cvm->ycid) / level2;
  }

  volvect2 = dvector(0, cvm->nbOfPillars - 1);
  volloc2 = dvector(0, cvm->nbOfPillars - 1);

  exercise2 =
      add_unit((long)(start2), -cvm->spotLag, SRT_BDAY, MODIFIED_SUCCEEDING);
  exercise_time2 = (exercise2 - lToday) / 365.0;

  exercise = (long)(DMIN(exercise1, DMIN(exercise2, expiry)));
  exercise_time = DMIN(exercise_time1, DMIN(exercise_time2, expiry_time));

  x = dvector(0, cvm->nQuadLegendre - 1);
  w = dvector(0, cvm->nQuadLegendre - 1);
  dateIndex = 0;
  dateIndex2 = 0;
  covscal = 0.0;
  volscal1 = 0.0;
  volscal2 = 0.0;
  if (exercise <= cvm->vol_dates[0]) {
    GaussLeg(time, exercise_time, x - 1, w - 1, cvm->nQuadLegendre);
    for (k = 0; k < cvm->nQuadLegendre; ++k) {
      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volvect1[j] = 0.0;
        volvect2[j] = 0.0;
      }

      for (i = 0; i < iNFloatDates1; ++i) {
        err = cvm_compute_fwdvol(x[k], (lFloatStartDates1[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc1, cvm);
        if (err)
          goto FREE_RETURN;
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect1[j] += weights1[i] * volloc1[j];
        }
      }

      for (i = 0; i < iNFloatDates2; ++i) {
        err = cvm_compute_fwdvol(x[k], (lFloatStartDates2[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc2, cvm);
        if (err)
          goto FREE_RETURN;
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect2[j] += weights2[i] * volloc2[j];
        }
      }

      for (j = 0; j < cvm->nbOfPillars; ++j) {
        volscal1 += w[k] * volvect1[j] * volvect1[j];
        volscal2 += w[k] * volvect2[j] * volvect2[j];
        covscal += w[k] * volvect1[j] * volvect2[j];
      }
    }
  } else {
    while ((dateIndex < cvm->nbOfDates) &&
           (exercise > cvm->vol_dates[dateIndex])) {
      dateIndex = dateIndex + 1;
    }

    while ((dateIndex2 < cvm->nbOfDates) &&
           (date > cvm->vol_dates[dateIndex2])) {
      dateIndex2 = dateIndex2 + 1;
    }

    if ((dateIndex2 < cvm->nbOfDates) && (dateIndex2 < dateIndex)) {
      GaussLeg(time, cvm->vol_times[dateIndex2], x - 1, w - 1,
               cvm->nQuadLegendre);
      for (k = 0; k < cvm->nQuadLegendre; ++k) {
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect1[j] = 0.0;
          volvect2[j] = 0.0;
        }
        for (i = 0; i < iNFloatDates1; ++i) {
          err =
              cvm_compute_fwdvol(x[k], (lFloatStartDates1[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc1, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect1[j] += weights1[i] * volloc1[j];
          }
        }
        for (i = 0; i < iNFloatDates2; ++i) {
          err =
              cvm_compute_fwdvol(x[k], (lFloatStartDates2[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc2, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect2[j] += weights2[i] * volloc2[j];
          }
        }

        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volscal1 += w[k] * volvect1[j] * volvect1[j];
          volscal2 += w[k] * volvect2[j] * volvect2[j];
          covscal += w[k] * volvect1[j] * volvect2[j];
        }
      }

      for (l = dateIndex2 + 1; l < dateIndex; ++l) {
        GaussLeg(cvm->vol_times[l - 1], cvm->vol_times[l], x - 1, w - 1,
                 cvm->nQuadLegendre);
        for (k = 0; k < cvm->nQuadLegendre; ++k) {
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect1[j] = 0.0;
            volvect2[j] = 0.0;
          }
          for (i = 0; i < iNFloatDates1; ++i) {
            err = cvm_compute_fwdvol(x[k],
                                     (lFloatStartDates1[i] - lToday) / 365.0,
                                     cvm->nbOfPillars, volloc1, cvm);
            if (err)
              goto FREE_RETURN;
            for (j = 0; j < cvm->nbOfPillars; ++j) {
              volvect1[j] += weights1[i] * volloc1[j];
            }
          }
          for (i = 0; i < iNFloatDates2; ++i) {
            err = cvm_compute_fwdvol(x[k],
                                     (lFloatStartDates2[i] - lToday) / 365.0,
                                     cvm->nbOfPillars, volloc2, cvm);
            if (err)
              goto FREE_RETURN;
            for (j = 0; j < cvm->nbOfPillars; ++j) {
              volvect2[j] += weights2[i] * volloc2[j];
            }
          }

          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volscal1 += w[k] * volvect1[j] * volvect1[j];
            volscal2 += w[k] * volvect2[j] * volvect2[j];
            covscal += w[k] * volvect1[j] * volvect2[j];
          }
        }
      }

      GaussLeg(cvm->vol_times[dateIndex - 1], exercise_time, x - 1, w - 1,
               cvm->nQuadLegendre);
      for (k = 0; k < cvm->nQuadLegendre; ++k) {
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect1[j] = 0.0;
          volvect2[j] = 0.0;
        }
        for (i = 0; i < iNFloatDates1; ++i) {
          err =
              cvm_compute_fwdvol(x[k], (lFloatStartDates1[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc1, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect1[j] += weights1[i] * volloc1[j];
          }
        }
        for (i = 0; i < iNFloatDates2; ++i) {
          err =
              cvm_compute_fwdvol(x[k], (lFloatStartDates2[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc2, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect2[j] += weights2[i] * volloc2[j];
          }
        }

        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volscal1 += w[k] * volvect1[j] * volvect1[j];
          volscal2 += w[k] * volvect2[j] * volvect2[j];
          covscal += w[k] * volvect1[j] * volvect2[j];
        }
      }
    } else {
      GaussLeg(time, exercise_time, x - 1, w - 1, cvm->nQuadLegendre);
      for (k = 0; k < cvm->nQuadLegendre; ++k) {
        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volvect1[j] = 0.0;
          volvect2[j] = 0.0;
        }
        for (i = 0; i < iNFloatDates1; ++i) {
          err =
              cvm_compute_fwdvol(x[k], (lFloatStartDates1[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc1, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect1[j] += weights1[i] * volloc1[j];
          }
        }
        for (i = 0; i < iNFloatDates2; ++i) {
          err =
              cvm_compute_fwdvol(x[k], (lFloatStartDates2[i] - lToday) / 365.0,
                                 cvm->nbOfPillars, volloc2, cvm);
          if (err)
            goto FREE_RETURN;
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            volvect2[j] += weights2[i] * volloc2[j];
          }
        }

        for (j = 0; j < cvm->nbOfPillars; ++j) {
          volscal1 += w[k] * volvect1[j] * volvect1[j];
          volscal2 += w[k] * volvect2[j] * volvect2[j];
          covscal += w[k] * volvect1[j] * volvect2[j];
        }
      }
    }
  }

  *vol1 = sqrt(volscal1 / (exercise_time - time));
  *vol2 = sqrt(volscal2 / (exercise_time - time));
  *cov = covscal / (exercise_time - time);

FREE_RETURN:

  if (lFixPayDates1)
    free(lFixPayDates1);
  if (lFixStartDates1)
    free(lFixStartDates1);
  if (lFixEndDates1)
    free(lFixEndDates1);
  if (dFixCoverages1)
    free(dFixCoverages1);

  if (lFloatPayDates1)
    free(lFloatPayDates1);
  if (lFloatFixingDates1)
    free(lFloatFixingDates1);
  if (lFloatStartDates1)
    free(lFloatStartDates1);
  if (lFloatEndDates1)
    free(lFloatEndDates1);
  if (dFloatCoverages1)
    free(dFloatCoverages1);
  if (dFloatSpreads1)
    free(dFloatSpreads1);

  if (weights1)
    free_dvector(weights1, 0, iNFloatDates1 - 1);

  if (volvect1)
    free_dvector(volvect1, 0, cvm->nbOfPillars - 1);
  if (volloc1)
    free_dvector(volloc1, 0, cvm->nbOfPillars - 1);

  if (lFixPayDates2)
    free(lFixPayDates2);
  if (lFixStartDates2)
    free(lFixStartDates2);
  if (lFixEndDates2)
    free(lFixEndDates2);
  if (dFixCoverages2)
    free(dFixCoverages2);

  if (lFloatPayDates2)
    free(lFloatPayDates2);
  if (lFloatFixingDates2)
    free(lFloatFixingDates2);
  if (lFloatStartDates2)
    free(lFloatStartDates2);
  if (lFloatEndDates2)
    free(lFloatEndDates2);
  if (dFloatCoverages2)
    free(dFloatCoverages2);
  if (dFloatSpreads2)
    free(dFloatSpreads2);

  if (weights2)
    free_dvector(weights2, 0, iNFloatDates2 - 1);

  if (volvect2)
    free_dvector(volvect2, 0, cvm->nbOfPillars - 1);
  if (volloc2)
    free_dvector(volloc2, 0, cvm->nbOfPillars - 1);

  if (x)
    free_dvector(x, 0, cvm->nQuadLegendre - 1);
  if (w)
    free_dvector(w, 0, cvm->nQuadLegendre - 1);

  return err;
}

Err cvm_compute_sensitivities(double time, double maturity, int nPillars,
                              double *sensi, cvm_str *cvm) {
  Err err = NULL;
  int i, dateIndex, index;
  long start, end;
  double sum;
  double *lambda = NULL;

  if ((maturity - time) <= cvm->pillars_times[0]) {
    sensi[0] = 1.0;
    for (i = 1; i < cvm->nbOfPillars; ++i) {
      sensi[i] = 0.0;
    }
    goto FREE_RETURN;
  } else if (maturity - time >= cvm->pillars_times[cvm->nbOfPillars - 1]) {
    sensi[cvm->nbOfPillars - 1] = 1.0;
    for (i = 0; i < cvm->nbOfPillars - 1; ++i) {
      sensi[i] = 0.0;
    }
    goto FREE_RETURN;
  }

  for (i = 0; i < cvm->nbOfPillars; ++i) {
    start = (long)(cvm->today + 365 * time);
    err = add_tenor(start, cvm->pillars_tenors[i], NO_BUSDAY_CONVENTION, &end);
    cvm->pillars_times[i] = (end - start) / 365.0;
    cvm->pillars_maturities[i] = cvm->today + 365.0 * cvm->pillars_times[i];
  }

  if (cvm->sensi_type == 0) // Parametric Sensitivities
  {
    lambda = dvector(0, cvm->nbOfPillars - 1);

    if (time <= cvm->vol_times[0]) {
      for (i = 0; i < cvm->nbOfPillars - 1; ++i) {
        lambda[i] = cvm->lambdas[0][i];
      }
    } else {
      dateIndex = 0;
      while ((dateIndex < cvm->nbOfDates) &&
             (time > cvm->vol_times[dateIndex])) {
        dateIndex = dateIndex + 1;
      }

      if (dateIndex == cvm->nbOfDates) {
        for (i = 0; i < cvm->nbOfPillars - 1; ++i) {
          lambda[i] = cvm->lambdas[dateIndex - 1][i];
        }
      } else {
        for (i = 0; i < cvm->nbOfPillars - 1; ++i) {
          lambda[i] = cvm->lambdas[dateIndex][i];
        }
      }
    }

    sum = 0.0;
    for (i = 0; i < cvm->nbOfPillars - 1; ++i) {
      err = splint(cvm->pillars_times - 1, cvm->sensitivities[i] - 1,
                   cvm->sensitivities2[i] - 1, cvm->nbOfPillars,
                   maturity - time, &sensi[i]);
      if (err)
        goto FREE_RETURN;
      sensi[i] = sensi[i] *
                 exp(-lambda[i] * (maturity - time - cvm->pillars_times[i]));
      sum += sensi[i];
    }
    sensi[i] = 1.0 - sum;
  } else if (cvm->sensi_type == 1) // Custom Sensitivities
  {
    sum = 0.0;
    for (i = 0; i < cvm->nbOfPillars - 1; ++i) {
      err = splint(cvm->sensi_times - 1, cvm->sensitivities[i] - 1,
                   cvm->sensitivities2[i] - 1, cvm->nbOfSensiTenors,
                   maturity - time, &sensi[i]);
      if (err)
        goto FREE_RETURN;
      sum += sensi[i];
    }
    sensi[i] = 1.0 - sum;
  } else // Linear Sensitivities
  {
    index = 0;
    while ((index < cvm->nbOfPillars) &&
           ((maturity - time) > cvm->pillars_times[index])) {
      index = index + 1;
    }

    if (index == cvm->nbOfPillars) {
      for (i = 0; i < cvm->nbOfPillars - 1; ++i) {
        sensi[i] = 0.0;
      }
      sensi[cvm->nbOfPillars - 1] = 1.0;
    } else {
      for (i = 0; i < cvm->nbOfPillars; ++i) {
        sensi[i] = 0;
      }
      sensi[index] =
          cvm->sensitivities[index][index - 1] +
          (cvm->sensitivities[index][index] -
           cvm->sensitivities[index][index - 1]) /
              (cvm->pillars_times[index] - cvm->pillars_times[index - 1]) *
              (maturity - time - cvm->pillars_times[index - 1]);
      sensi[index - 1] = 1 - sensi[index];
    }
  }

FREE_RETURN:
  if (lambda)
    free_dvector(lambda, 0, cvm->nbOfPillars - 1);

  return err;
}

Err cvm_free_und_struct(SrtUndPtr pUndDesc) {
  Err err = NULL;
  err = cvm_free_struct(pUndDesc->spec_desc);
  free(pUndDesc->spec_desc);
  free(pUndDesc);
  pUndDesc = NULL;
  return err;
}

Err cvm_free_struct(cvm_str *cvm) {
  Err err = NULL;
  int i, j;

  if (cvm) {
    if (cvm->vol_dates) {
      free_dvector(cvm->vol_dates, 0, cvm->nbOfDates - 1);
      cvm->vol_dates = NULL;
    }

    if (cvm->vol_times) {
      free_dvector(cvm->vol_times, 0, cvm->nbOfDates - 1);
      cvm->vol_times = NULL;
    }

    if (cvm->pillars_tenors) {
      free(cvm->pillars_tenors);
      cvm->pillars_tenors = NULL;
    }

    if (cvm->pillars_maturities) {
      free_dvector(cvm->pillars_maturities, 0, cvm->nbOfPillars - 1);
      cvm->pillars_maturities = NULL;
    }

    if (cvm->pillars_times) {
      free_dvector(cvm->pillars_times, 0, cvm->nbOfPillars - 1);
      cvm->pillars_times = NULL;
    }

    if (cvm->pillars_vols_norm) {
      for (i = 0; i < cvm->nbOfDates; ++i) {
        if (cvm->pillars_vols_norm[i]) {
          free(cvm->pillars_vols_norm[i]);
          cvm->pillars_vols_norm[i] = NULL;
        }
      }
      free(cvm->pillars_vols_norm);
      cvm->pillars_vols_norm = NULL;
    }

    if (cvm->pillars_vols) {
      for (i = 0; i < cvm->nbOfDates; ++i) {
        if (cvm->pillars_vols[i]) {
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            if (cvm->pillars_vols[i][j]) {
              free(cvm->pillars_vols[i][j]);
              cvm->pillars_vols[i][j] = NULL;
            }
          }
          free(cvm->pillars_vols[i]);
          cvm->pillars_vols[i] = NULL;
        }
      }
      free(cvm->pillars_vols);
      cvm->pillars_vols = NULL;
    }

    if (cvm->pillars_correls) {
      for (i = 0; i < cvm->nbOfDates; ++i) {
        if (cvm->pillars_correls[i]) {
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            if (cvm->pillars_correls[i][j]) {
              free(cvm->pillars_correls[i][j]);
              cvm->pillars_correls[i][j] = NULL;
            }
          }
          free(cvm->pillars_correls[i]);
          cvm->pillars_correls[i] = NULL;
        }
      }
      free(cvm->pillars_correls);
      cvm->pillars_correls = NULL;
    }

    if (cvm->pillars_angles) {
      for (i = 0; i < cvm->nbOfDates; ++i) {
        if (cvm->pillars_angles[i]) {
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            if (cvm->pillars_angles[i][j]) {
              free(cvm->pillars_angles[i][j]);
              cvm->pillars_angles[i][j] = NULL;
            }
          }
          free(cvm->pillars_angles[i]);
          cvm->pillars_angles[i] = NULL;
        }
      }
      free(cvm->pillars_angles);
      cvm->pillars_angles = NULL;
    }

    if (cvm->pillars_chol) {
      for (i = 0; i < cvm->nbOfDates; ++i) {
        if (cvm->pillars_chol[i]) {
          for (j = 0; j < cvm->nbOfPillars; ++j) {
            if (cvm->pillars_chol[i][j]) {
              free(cvm->pillars_chol[i][j]);
              cvm->pillars_chol[i][j] = NULL;
            }
          }
          free(cvm->pillars_chol[i]);
          cvm->pillars_chol[i] = NULL;
        }
      }
      free(cvm->pillars_chol);
      cvm->pillars_chol = NULL;
    }

    if (cvm->lambdas) {
      free_dmatrix(cvm->lambdas, 0, cvm->nbOfDates - 1, 0,
                   cvm->nbOfPillars - 1);
      cvm->lambdas = NULL;
    }

    if (cvm->sensitivities) {
      free_dmatrix(cvm->sensitivities, 0, cvm->nbOfPillars - 1, 0,
                   cvm->nbOfSensiTenors - 1);
      cvm->sensitivities = NULL;
    }

    if (cvm->sensitivities2) {
      free_dmatrix(cvm->sensitivities2, 0, cvm->nbOfPillars - 1, 0,
                   cvm->nbOfSensiTenors - 1);
      cvm->sensitivities2 = NULL;
    }

    if (cvm->sensi_tenors) {
      free(cvm->sensi_tenors);
      cvm->sensi_tenors = NULL;
    }

    if (cvm->sensi_times) {
      free_dvector(cvm->sensi_times, 0, cvm->nbOfSensiTenors - 1);
    }

    if (cvm->option_tenors) {
      free(cvm->option_tenors);
      cvm->option_tenors = NULL;
    }

    if (cvm->underlying_tenors) {
      free(cvm->underlying_tenors);
      cvm->underlying_tenors = NULL;
    }

    if (cvm->market_vols) {
      free_dmatrix(cvm->market_vols, 0, cvm->nbOfOptionTenors - 1, 0,
                   cvm->nbOfUnderlyingTenors - 1);
      cvm->market_vols = NULL;
    }

    if (cvm->model_vols) {
      free_dmatrix(cvm->model_vols, 0, cvm->nbOfOptionTenors - 1, 0,
                   cvm->nbOfUnderlyingTenors - 1);
      cvm->model_vols = NULL;
    }

    cvm = NULL;
  }

  return err;
}

char *SrtInitCVMUnd(char *undName,
                    // Market ID
                    char *ycname, char *vcname,
                    // Swaption Details
                    SrtCompounding srtFreq, SrtBasisCode srtBasis,
                    char *refRate,

                    // TS Dates
                    double shift, int converging_sliding,
                    int interpollation_mode, int nbOfDates, double *vol_dates,

                    // Pillars
                    int nbOfPillars, char **pillars_tenors,

                    // Pillars Correlations and Vols
                    double ***pillars_correls, double **pillars_vols,

                    // Sensitivities
                    int sensitype, double **lambdas,

                    int nbOfSensiTenors, char **sensi_tenors,
                    double **sensitivities,

                    // Numerical Parameter
                    int nQuadLegendre) {
  Err err = NULL;
  int bCleanUpUndFlag = 1;
  SrtUndPtr pUndDesc;
  SrtUndListPtr und_list;
  SrtCurvePtr pYieldCurve;
  char *ccy;
  cvm_str *cvm = NULL;

  cvm = calloc(1, sizeof(cvm_str));

  // Get the yield curve and the spot date from the yield curve name
  pYieldCurve = lookup_curve(ycname);
  ccy = get_curve_ccy(pYieldCurve);

  // Create the new CVM underlying
  und_list = get_underlying_list();
  pUndDesc = (SrtUndPtr)calloc(1, sizeof(SrtUndDesc));
  strcpy(pUndDesc->underl_name, undName);
  strupper(pUndDesc->underl_name);
  strip_white_space(pUndDesc->underl_name);
  strcpy(pUndDesc->underl_lbl, "CVM_UND");
  pUndDesc->underl_ccy = ccy;
  pUndDesc->underl_type = CVM_UND;

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
      sensitype, lambdas,

      nbOfSensiTenors, sensi_tenors, sensitivities,

      // Numerical Parameter
      nQuadLegendre,

      cvm);

  pUndDesc->spec_desc = cvm;

  // Put the underlying into the depot
  err = srt_f_lstins(und_list, pUndDesc->underl_name, 0.0, OBJ_PTR_UND,
                     (void *)pUndDesc, &cvm_free_und_struct,
                     &(pUndDesc->underl_ticker));

  return err;
}

Err cvm_get_struct_from_und(char *und, cvm_str **cvm) {
  Err err = NULL;
  SrtUndPtr pUndDesc;

  // Get the underlying through its name and check it exists
  // Check on the underlying type
  pUndDesc = lookup_und(und);
  if (!pUndDesc) {
    err = "Undefined underlying";
    return err;
  }
  if (!(ISUNDTYPE(pUndDesc, CVM_UND))) {
    err = "Not a CVM Underlying";
  }

  // Extract the information from the underlying
  *cvm = ((cvm_str *)((SrtUndDesc *)(pUndDesc->spec_desc)));

  return err;
}

Err sortspline(int nbasis, double *xbasis, double *ybasis, int n, double *x,
               double *y) {
  Err err = NULL;
  int i;
  double *ybasis2 = NULL;

  ybasis2 = (double *)calloc(nbasis, sizeof(double));

  spline(xbasis - 1, ybasis - 1, nbasis, 0, 0, ybasis2 - 1);

  for (i = 0; i < n; ++i) {
    err = splint(xbasis - 1, ybasis - 1, ybasis2 - 1, nbasis, x[i], &y[i]);
  }

  return err;
}
