/* ==========================================================================
   FILE_NAME:	McEBOptimsation.cxx

   PURPOSE:		Monte Carlo Exercise boundary optimisation for Callable

   DATE:		01/09/03
   ========================================================================== */

#include "MCEBOptimisation.h"
#include "math.h"
#include "srt_h_all.h"

void mceb_set_default_params(MCEBPARAMS sParams) {
  /* GRFN Read Params */
  sParams->iColBound = 0;
  sParams->iColPay = 0;

  /* MCEB Settings */
  sParams->lNbDates = 0;
  sParams->iIsKO = 1;
  sParams->iCallCurrent = 1;
  sParams->iMultiIndex = 0;
  sParams->iNbIndex = 1;

  sParams->iDoSmoothing = 1;
  sParams->dCallSpread = 0.0;

  /* Main Outputs */
  sParams->dBarrier = NULL;
  sParams->dCoefLin = NULL;

  /* Fees for Call */
  sParams->iHasFees = 0;
  sParams->iColFees = 0;
  sParams->iPayAllFeeCol = 0;

  /* Fwd IV calculation / adjustment */
  sParams->iCalcIV = 0;
  sParams->iAdjustIV = 0;
  sParams->iHasNumeraire = 0;
  sParams->iColNumeraire = 0;
  sParams->iAddMultFee = 0;
  sParams->dFee = NULL;
  sParams->dMarketFwdIV = NULL;
  sParams->dModelFwdIV = NULL;

  /* All the extra informations */
  sParams->iCalcOneTime = 0;
  sParams->iCalcOneTimePartial = 0;
  sParams->iCalcExeProba = 0;
  sParams->dOneTimeCall = NULL;
  sParams->dOneTimePartial = NULL;
  sParams->dExeProba = NULL;

  /* Extra Parameters */
  sParams->iDoInfos = 0;
  sParams->iKnockInCol = 0;
  sParams->iFindBestOptim = 0;
  sParams->iRemoveLastOnLast = 0;
  sParams->iAddNonOptimisedForKI = 0;
  sParams->iRescaleLinCoefs = 1;

  /* File Name where to save all infos */
  sParams->sFilePath[0] = '\0';
}

void mceb_copy_nondynamic_params(MCEBPARAMS sParamsCpy, MCEBPARAMS sParamsSrc) {
  /* GRFN Read Params */
  sParamsCpy->iColBound = sParamsSrc->iColBound;
  sParamsCpy->iColPay = sParamsSrc->iColPay;

  /* MCEB Settings */
  sParamsCpy->lNbDates = sParamsSrc->lNbDates;
  sParamsCpy->iIsKO = sParamsSrc->iIsKO;
  sParamsCpy->iCallCurrent = sParamsSrc->iCallCurrent;
  sParamsCpy->iMultiIndex = sParamsSrc->iMultiIndex;
  sParamsCpy->iNbIndex = sParamsSrc->iNbIndex;

  sParamsCpy->iDoSmoothing = sParamsSrc->iDoSmoothing;
  sParamsCpy->dCallSpread = sParamsSrc->dCallSpread;

  /* Main Outputs */

  /* Fwd IV calculation / adju		stment */
  sParamsCpy->iCalcIV = sParamsSrc->iCalcIV;
  sParamsCpy->iAdjustIV = sParamsSrc->iAdjustIV;
  sParamsCpy->iHasNumeraire = sParamsSrc->iHasNumeraire;
  sParamsCpy->iColNumeraire = sParamsSrc->iColNumeraire;
  sParamsCpy->iAddMultFee = sParamsSrc->iAddMultFee;

  sParamsCpy->iHasFees = sParamsSrc->iHasFees;
  sParamsCpy->iColFees = sParamsSrc->iColFees;
  sParamsCpy->iPayAllFeeCol = sParamsSrc->iPayAllFeeCol;

  /* All the extra informations */
  sParamsCpy->iCalcOneTime = sParamsSrc->iCalcOneTime;
  sParamsCpy->iCalcOneTimePartial = sParamsSrc->iCalcOneTimePartial;
  sParamsCpy->iCalcExeProba = sParamsSrc->iCalcExeProba;

  /* Extra Parameters */
  sParamsCpy->iDoInfos = sParamsSrc->iDoInfos;
  sParamsCpy->iKnockInCol = sParamsSrc->iKnockInCol;
  sParamsCpy->iFindBestOptim = sParamsSrc->iFindBestOptim;
  sParamsCpy->iRemoveLastOnLast = sParamsSrc->iRemoveLastOnLast;
  sParamsCpy->iAddNonOptimisedForKI = sParamsSrc->iAddNonOptimisedForKI;
  sParamsCpy->iRescaleLinCoefs = sParamsSrc->iRescaleLinCoefs;

  /* File Name where to save all infos */
  memcpy(sParamsCpy->sFilePath, sParamsSrc->sFilePath,
         sizeof(sParamsSrc->sFilePath) * sizeof(char));
}

Err mceb_allocate_params(MCEBPARAMS sParams, long lNbDates) {
  if (sParams) {
    sParams->lNbDates = lNbDates;

    sParams->dBarrier = calloc(lNbDates, sizeof(double));
    sParams->dCoefLin = dmatrix(0, lNbDates - 1, 0, sParams->iNbIndex);

    if (!sParams->dBarrier || !sParams->dCoefLin) {
      return "Memory allocation faillure in mceb_allocate_params";
    }

    if (sParams->iCalcIV || sParams->iAdjustIV) {
      sParams->dModelFwdIV = calloc(lNbDates, sizeof(double));
      sParams->dFee = calloc(lNbDates, sizeof(double));
      sParams->dMarketFwdIV = calloc(lNbDates, sizeof(double));

      if (!sParams->dModelFwdIV || !sParams->dFee || !sParams->dMarketFwdIV) {
        return "Memory allocation faillure in mceb_allocate_params";
      }
    }

    if (sParams->iCalcExeProba) {
      sParams->dExeProba = calloc(lNbDates, sizeof(double));

      if (!sParams->dExeProba) {
        return "Memory allocation faillure in mceb_allocate_params";
      }
    }

    if (sParams->iCalcOneTime) {
      sParams->dOneTimeCall = calloc(lNbDates, sizeof(double));

      if (!sParams->dOneTimeCall) {
        return "Memory allocation faillure in mceb_allocate_params";
      }
    }

    if (sParams->iCalcOneTimePartial) {
      sParams->dOneTimePartial = calloc(lNbDates, sizeof(double));

      if (!sParams->dOneTimePartial) {
        return "Memory allocation faillure in mceb_allocate_params";
      }
    }
  }

  return NULL;
}

void mceb_shift_extrainfos(MCEBPARAMS sParams) {
  int i;

  for (i = 0; i < sParams->lNbDates - 1; i++) {
    if (sParams->iAdjustIV || sParams->iCalcIV) {
      sParams->dModelFwdIV[i] = sParams->dModelFwdIV[i + 1];
    }

    if (sParams->iAdjustIV) {
      sParams->dFee[i] = sParams->dFee[i + 1];
    }

    if (sParams->iCalcOneTime) {
      sParams->dOneTimeCall[i] = sParams->dOneTimeCall[i + 1];
    }

    if (sParams->iCalcOneTimePartial) {
      sParams->dOneTimePartial[i] = sParams->dOneTimePartial[i + 1];
    }

    if (sParams->iCalcExeProba) {
      sParams->dExeProba[i] = sParams->dExeProba[i + 1];
    }
  }
}

void mceb_free_params(MCEBPARAMS sParams) {
  if (sParams) {
    if (sParams->dBarrier)
      free(sParams->dBarrier);
    if (sParams->dCoefLin)
      free_dmatrix(sParams->dCoefLin, 0, sParams->lNbDates - 1, 0,
                   sParams->iNbIndex);
    if (sParams->dFee)
      free(sParams->dFee);
    if (sParams->dMarketFwdIV)
      free(sParams->dMarketFwdIV);
    if (sParams->dModelFwdIV)
      free(sParams->dModelFwdIV);
    if (sParams->dOneTimeCall)
      free(sParams->dOneTimeCall);
    if (sParams->dOneTimePartial)
      free(sParams->dOneTimePartial);
    if (sParams->dExeProba)
      free(sParams->dExeProba);
  }
}

Err mceb_allocate_savevalues_for_GRFN(long iNbPaths, long iNbEvent,
                                      MCEBPARAMS sParams,
                                      double ****dSaveValues) {
  *dSaveValues =
      f3tensor(0, iNbEvent - 1, 0,
               sParams->iNbIndex + sParams->iHasNumeraire + sParams->iHasFees,
               0, iNbPaths - 1);

  if (!(*dSaveValues)) {
    return "Memory allocation faillure in mceb_allocate_savevalues_for_GRFN";
  } else {
    return NULL;
  }
}

void mceb_free_savevalues_for_GRFN(double ***dSaveValues, long iNbPaths,
                                   long iNbEvent, MCEBPARAMS sParams)

{
  if (dSaveValues) {
    free_f3tensor(dSaveValues, 0, iNbEvent - 1, 0,
                  sParams->iNbIndex + sParams->iHasNumeraire +
                      sParams->iHasFees,
                  0, iNbPaths - 1);
  }
}

void mceb_fill_savevalues_from_GRFN(double **dSaveValues, double *dResEvt,
                                    long lSimulIdx, double dDF,
                                    MCEBPARAMS sParams) {
  int k;

  for (k = 0; k < sParams->iNbIndex; k++) {
    dSaveValues[k][lSimulIdx] = dResEvt[sParams->iColBound + k];
  }

  if (!sParams->iKnockInCol) {
    dSaveValues[sParams->iNbIndex][lSimulIdx] = dResEvt[sParams->iColPay] / dDF;
  } else {
    dSaveValues[sParams->iNbIndex][lSimulIdx] = dResEvt[sParams->iColPay];
  }

  if (sParams->iHasNumeraire) {
    dSaveValues[sParams->iNbIndex + 1][lSimulIdx] =
        dResEvt[sParams->iColNumeraire] / dDF;
  }

  if (sParams->iHasFees) {
    dSaveValues[sParams->iNbIndex + sParams->iHasNumeraire + 1][lSimulIdx] =
        dResEvt[sParams->iColFees] / dDF;
  }
}

Err mceb_adjust_fwdiv(double ***save_values, long start_date_idx,
                      long end_date_idx, long nb_paths, int *optimise,
                      MCEBPARAMS params) {
  long i, j;
  double value, numeraire, fee;
  double *save_valuesjpay, *savevaluesjnum;
  Err err = NULL;

  if (params->iCalcIV || params->iAdjustIV) {
    for (j = start_date_idx; j <= end_date_idx; j++) {
      save_valuesjpay = save_values[j][params->iNbIndex];

      /* Calculate IV */
      value = 0.0;

      for (i = 0; i < nb_paths; i++) {
        value += save_valuesjpay[i];
      }

      value /= nb_paths;

      params->dModelFwdIV[j] = value;

      if (params->iAdjustIV) {
        /* Calculate Numeraire */
        if (params->iHasNumeraire && params->iAddMultFee == 0) {
          savevaluesjnum = save_values[j][params->iNbIndex + 1];

          numeraire = 0.0;

          for (i = 0; i < nb_paths; i++) {
            numeraire += savevaluesjnum[i];
          }

          numeraire /= nb_paths;
        } else {
          numeraire = 1.0;
        }

        if (params->iAddMultFee == 0) {
          if (fabs(numeraire) > 1.0E-10) {
            fee = (params->dMarketFwdIV[j] - value) / numeraire;
          } else if (fabs(params->dMarketFwdIV[j] - value) < 1.0E-10) {
            fee = 0.0;
          } else {
            err = "Cannot adjust Fwd IV in MCEB because of null numeraire";
            return err;
          }

          if (params->iHasNumeraire) {
            for (i = 0; i < nb_paths; i++) {
              save_valuesjpay[i] += fee * savevaluesjnum[i];
            }
          } else {
            for (i = 0; i < nb_paths; i++) {
              save_valuesjpay[i] += fee;
            }
          }
        } else {
          fee = params->dMarketFwdIV[j] / value;

          for (i = 0; i < nb_paths; i++) {
            save_valuesjpay[i] *= fee;
          }
        }

        params->dFee[j] = fee;
      }
    }
  }

  return err;
}

Err find_and_optimise_boundary(double ***save_values, long nb_dates,
                               long nb_paths, int *optimise, MCEBPARAMS params,
                               double *value, double *error) {
  double *next_payoff = NULL, **cur_payoff = NULL, **x = NULL, *y = NULL,
         **mat = NULL, *total_payoff = NULL, *b, **save_valuesj;
  int *call_col = NULL;

  double value_max, payoff, sum_next;
  double res, res2;

  long i, j, k, index;
  int add_one_current;
  int iNbIndex;
  long first_optim;
  double inv_numpaths;

  double *dStdCallSpread = NULL;
  int iDoSmoothing;
  double dNotional, dValue, dValue2, dIndex;

  FILE *stream = NULL;

  Err err = NULL;

  /* Adjust Fwd IV */
  err = mceb_adjust_fwdiv(save_values, 0, nb_dates - 1, nb_paths, optimise,
                          params);

  if (err)
    goto FREE_RETURN;

  inv_numpaths = 1.0 / nb_paths;

  if (params->iCallCurrent) {
    add_one_current = 0;
  } else {
    add_one_current = 1;
  }

  /* Allocate new variables */
  next_payoff = dvector(0, nb_paths - 1);
  cur_payoff = dmatrix(0, 3 - params->iIsKO, 0, nb_paths - 1);
  mat = dmatrix(0, params->iNbIndex, 0, nb_paths - 1);

  if (params->iMultiIndex) {
    x = dmatrix(0, params->iNbIndex, 0, nb_paths - 1);
    y = dvector(0, nb_paths - 1);
  }

  if (!next_payoff || !cur_payoff || !mat ||
      (params->iMultiIndex && (!x || !y))) {
    err = "Memory allocation failure in optimise_boundary";
    goto FREE_RETURN;
  }

  if (params->iCalcOneTime) {
    total_payoff = calloc(nb_paths, sizeof(double));

    if (!total_payoff) {
      err = "Memory allocation failure in optimise_boundary";
      goto FREE_RETURN;
    }
  }

  if (params->iCalcExeProba) {
    call_col = calloc(nb_paths, sizeof(int));

    if (!call_col) {
      err = "Memory allocation failure in optimise_boundary";
      goto FREE_RETURN;
    }

    for (i = 0; i < nb_paths; i++) {
      call_col[i] = -1;
    }
  }

  iDoSmoothing = 0;

  if (params->iDoSmoothing && fabs(params->dCallSpread) > 1.0E-08) {
    iDoSmoothing = 1;

    dStdCallSpread = calloc(nb_dates, sizeof(double));

    if (!dStdCallSpread) {
      err = "Memory allocation failure in optimise_boundary";
      goto FREE_RETURN;
    }
  }

  if (strlen(params->sFilePath) > 0) {
    stream = fopen(params->sFilePath, "w+");
    if (!stream)
      goto FREE_RETURN;
  }

  first_optim = 1;

  if (params->iIsKO) {
    /* Knock Out cas */

    /*	Perform the optimisation	*/
    for (j = nb_dates - 1; j >= 0; j--) {
      if (first_optim && params->iNbIndex > 1 &&
          params->iRemoveLastOnLast == 1) {
        /* Remove last index for last date */
        iNbIndex = params->iNbIndex - 1;
      } else {
        iNbIndex = params->iNbIndex;
      }

      if (optimise[j]) {
        /* First compute the partial one time call if required */
        if (!first_optim && params->iCalcOneTimePartial) {
          /* Regress the payoff */
          b = params->dCoefLin[j];

          if (params->iMultiIndex) {
            save_valuesj = save_values[j];

            for (i = 0; i < nb_paths; i++) {
              cur_payoff[1][i] = total_payoff[i];

              if (params->iHasFees) {
                /* Substract Paid Fees */
                cur_payoff[1][i] += save_valuesj[params->iNbIndex +
                                                 params->iHasNumeraire + 1][i];
              }

              y[i] = -cur_payoff[1][i];
              x[0][i] = 1.0;
            }

            for (k = 0; k < iNbIndex; k++) {
              memcpy(x[k + 1], save_valuesj[k], nb_paths * sizeof(double));
            }

            err = find_linear_coef(iNbIndex + 1, nb_paths, y, x, b, mat);

            while (err && iNbIndex > 1) {
              iNbIndex--;
              err = find_linear_coef(iNbIndex + 1, nb_paths, y, x, b, mat);
            }

            if (err)
              goto FREE_RETURN;

            /* Rescale the coeff */
            if (params->iRescaleLinCoefs) {
              res = 0.0;
              for (k = 0; k < iNbIndex; k++) {
                res += b[k + 1];
              }

              if (fabs(res) < 1.0E-08) {
                res = 1.0;
              } else {
                res = fabs(1.0 / res);
              }

              for (k = 0; k < iNbIndex; k++) {
                b[k + 1] *= res;
              }
            }

            /* Update the payoff */

            /* add the current PV's */
            for (i = 0; i < nb_paths; i++) {
              cur_payoff[0][i] = 0.0;

              for (k = 0; k < iNbIndex; k++) {
                cur_payoff[0][i] += b[k + 1] * save_valuesj[k][i];
              }
            }
          } else {
            for (i = 0; i < nb_paths; i++) {
              cur_payoff[0][i] = save_values[j][0][i];
              cur_payoff[1][i] = total_payoff[i];
            }

            b[1] = 1.0;
          }

          /* Sort this vector by increasing index */
          num_f_sort_matrix(nb_paths, 2, 0, cur_payoff);

          /* Find the maximum */
          value_max = payoff = 0.0;
          index = -1;

          for (i = 0; i < nb_paths; i++) {
            payoff += cur_payoff[1][i];

            if (payoff > value_max) {
              value_max = payoff;
              index = i;
            }
          }

          params->dOneTimePartial[j] = (value_max - payoff) / nb_paths;
        }

        /* Not anymore first optim */
        first_optim = 0;

        /* Update the coupons for One Time Callable */
        if (params->iCalcOneTime || params->iCalcOneTimePartial) {
          k = j;
          while (k < nb_dates - add_one_current && (k == j || !optimise[k])) {
            for (i = 0; i < nb_paths; i++) {
              total_payoff[i] +=
                  save_values[k + add_one_current][params->iNbIndex][i];
            }

            k++;
          }
        }

        if (params->iCalcOneTime) {
          /* Regress the payoff */
          b = params->dCoefLin[j];

          if (params->iMultiIndex) {
            save_valuesj = save_values[j];

            for (i = 0; i < nb_paths; i++) {
              cur_payoff[1][i] = total_payoff[i];

              if (params->iHasFees) {
                /* Substract Paid Fees */
                cur_payoff[1][i] += save_valuesj[params->iNbIndex +
                                                 params->iHasNumeraire + 1][i];
              }

              y[i] = -cur_payoff[1][i];
              x[0][i] = 1.0;
            }

            for (k = 0; k < iNbIndex; k++) {
              memcpy(x[k + 1], save_valuesj[k], nb_paths * sizeof(double));
            }

            err = find_linear_coef(iNbIndex + 1, nb_paths, y, x, b, mat);

            while (err && iNbIndex > 1) {
              iNbIndex--;
              err = find_linear_coef(iNbIndex + 1, nb_paths, y, x, b, mat);
            }

            if (err)
              goto FREE_RETURN;

            /* Rescale the coeff */
            if (params->iRescaleLinCoefs) {
              res = 0.0;
              for (k = 0; k < iNbIndex; k++) {
                res += b[k + 1];
              }

              if (fabs(res) < 1.0E-08) {
                res = 1.0;
              } else {
                res = fabs(1.0 / res);
              }

              for (k = 0; k < iNbIndex; k++) {
                b[k + 1] *= res;
              }
            }

            /* Update the payoff */

            /* add the current PV's */
            for (i = 0; i < nb_paths; i++) {
              cur_payoff[0][i] = 0.0;

              for (k = 0; k < iNbIndex; k++) {
                cur_payoff[0][i] += b[k + 1] * save_valuesj[k][i];
              }
            }
          } else {
            for (i = 0; i < nb_paths; i++) {
              cur_payoff[0][i] = save_values[j][0][i];
              cur_payoff[1][i] = total_payoff[i];
            }

            b[1] = 1.0;
          }

          /* Sort this vector by increasing index */
          num_f_sort_matrix(nb_paths, 2, 0, cur_payoff);

          /* Find the maximum */
          value_max = payoff = 0.0;
          index = -1;

          for (i = 0; i < nb_paths; i++) {
            payoff += cur_payoff[1][i];

            if (payoff > value_max) {
              value_max = payoff;
              index = i;
            }
          }

          params->dOneTimeCall[j] = (value_max - payoff) / nb_paths;
        }

        /* Regress the payoff */
        b = params->dCoefLin[j];
        save_valuesj = save_values[j];

        if (params->iMultiIndex) {
          for (i = 0; i < nb_paths; i++) {
            cur_payoff[1][i] =
                next_payoff[i] +
                save_values[j + add_one_current][params->iNbIndex][i];
            y[i] = -cur_payoff[1][i];

            if (params->iHasFees) {
              /* Substract Paid Fees */
              y[i] -=
                  save_valuesj[params->iNbIndex + params->iHasNumeraire + 1][i];
            }

            x[0][i] = 1.0;
          }

          for (k = 0; k < iNbIndex; k++) {
            memcpy(x[k + 1], save_valuesj[k], nb_paths * sizeof(double));
          }

          err = find_linear_coef(iNbIndex + 1, nb_paths, y, x, b, mat);

          while (err && iNbIndex > 1) {
            iNbIndex--;
            err = find_linear_coef(iNbIndex + 1, nb_paths, y, x, b, mat);
          }

          if (err)
            goto FREE_RETURN;

          /* Rescale the coeff */
          if (params->iRescaleLinCoefs) {
            res = 0.0;
            for (k = 0; k < iNbIndex; k++) {
              res += b[k + 1];
            }

            if (fabs(res) < 1.0E-08) {
              res = 1.0;
            } else {
              res = fabs(1.0 / res);
            }

            for (k = 0; k < iNbIndex; k++) {
              b[k + 1] *= res;
            }

            b[0] = 0.0;
          }

          /* Update the Index */
          for (i = 0; i < nb_paths; i++) {
            cur_payoff[0][i] = b[0];

            for (k = 0; k < iNbIndex; k++) {
              cur_payoff[0][i] += b[k + 1] * save_valuesj[k][i];
            }

            cur_payoff[2][i] = i;
          }
        } else {
          for (i = 0; i < nb_paths; i++) {
            cur_payoff[0][i] = save_valuesj[0][i];
            cur_payoff[1][i] =
                next_payoff[i] +
                save_values[j + add_one_current][params->iNbIndex][i];
            cur_payoff[2][i] = i;
          }

          b[1] = 1.0;
        }

        /* Calculate STD for Smoothing */
        if (iDoSmoothing) {
          dValue = 0.0;
          dValue2 = 0.0;

          for (i = 0; i < nb_paths; i++) {
            dValue += cur_payoff[0][i];
            dValue2 += cur_payoff[0][i] * cur_payoff[0][i];
          }

          dStdCallSpread[j] =
              dValue2 / nb_paths - (dValue / nb_paths * dValue / nb_paths);

          if (dStdCallSpread[j] > 1.0E-10) {
            dStdCallSpread[j] = sqrt(dStdCallSpread[j]);
          } else {
            dStdCallSpread[j] = 1.0E-05;
          }
        }

        /* Save into file */
        if (stream) {
          fprintf(stream, "Optimising Option %d\n", j + 1);

          for (i = 0; i < nb_paths; i++) {
            fprintf(stream, "%f	", cur_payoff[1][i]);

            for (k = 0; k < iNbIndex; k++) {
              fprintf(stream, "%f	", save_valuesj[k][i]);
            }

            fprintf(stream, "\n");
          }

          fprintf(stream, "\n");
        }

        /* Sort this vector by increasing index */
        num_f_sort_matrix(nb_paths, 3, 0, cur_payoff);

        /* Find the maximum */

        value_max = payoff = 0.0;
        index = -1;

        for (i = 0; i < nb_paths; i++) {
          payoff += cur_payoff[1][i];

          if (params->iHasFees) {
            payoff +=
                save_valuesj[params->iNbIndex + params->iHasNumeraire + 1][i];
          }

          if (payoff > value_max) {
            value_max = payoff;
            index = i;
          }
        }

        if (index >= nb_paths - 1) {
          params->dBarrier[j] = cur_payoff[0][nb_paths - 1] +
                                0.1 * fabs(cur_payoff[0][nb_paths - 1]);
        } else if (index >= 0) {
          params->dBarrier[j] =
              0.5 * (cur_payoff[0][index] + cur_payoff[0][index + 1]);
        } else {
          params->dBarrier[j] =
              cur_payoff[0][0] - 0.000001 * fabs(cur_payoff[0][0]);
        }

        /* Update payoff */
        if (params->iHasFees) {
          for (i = 0; i <= index; i++) {
            next_payoff[(int)cur_payoff[2][i]] = cur_payoff[1][i];
          }

          for (i = index + 1; i < nb_paths; i++) {
            next_payoff[(int)cur_payoff[2][i]] =
                -save_valuesj[params->iNbIndex + params->iHasNumeraire + 1][i];
          }
        } else {
          for (i = 0; i <= index; i++) {
            next_payoff[(int)cur_payoff[2][i]] = cur_payoff[1][i];
          }

          for (i = index + 1; i < nb_paths; i++) {
            next_payoff[(int)cur_payoff[2][i]] = 0.0;
          }
        }

        if (params->iCalcExeProba) {
          for (i = index + 1; i < nb_paths; i++) {
            call_col[(int)cur_payoff[2][i]] = j;
          }
        }
      } else {
        if (params->iCallCurrent) {
          for (i = 0; i < nb_paths; i++) {
            next_payoff[i] += save_values[j][params->iNbIndex][i];
          }
        } else if (j < nb_dates - 1) {
          for (i = 0; i < nb_paths; i++) {
            next_payoff[i] += save_values[j + 1][params->iNbIndex][i];
          }
        }

        params->dBarrier[j] = 0.0;

        if (params->iCalcOneTime) {
          params->dOneTimeCall[j] = 0.0;
        }

        if (params->iCalcOneTimePartial) {
          params->dOneTimePartial[j] = 0.0;
        }
      }
    }

    /* add the last PV if needed */
    if (!params->iCallCurrent) {
      for (i = 0; i < nb_paths; i++) {
        next_payoff[i] += save_values[0][params->iNbIndex][i];
      }
    }
  } else {
    /* Knock In case */

    /*	Perform the optimisation	*/

    sum_next = 0.0;

    for (j = nb_dates - 1; j >= 0; j--) {
      if (j == nb_dates - 1 && params->iNbIndex > 1 &&
          params->iRemoveLastOnLast == 1) {
        /* Remove last index for last date */
        iNbIndex = params->iNbIndex - 1;
      } else {
        iNbIndex = params->iNbIndex;
      }

      if (optimise[j]) {
        b = params->dCoefLin[j];

        if (params->iMultiIndex) {
          /* Regress the payoff */
          save_valuesj = save_values[j];

          for (i = 0; i < nb_paths; i++) {
            cur_payoff[1][i] =
                save_values[j + add_one_current][params->iNbIndex][i];
            cur_payoff[2][i] = next_payoff[i];

            if (params->iHasFees) {
              /* Substract Paid Fees */
              cur_payoff[1][i] -=
                  save_valuesj[params->iNbIndex + params->iHasNumeraire + 1][i];
            }

            y[i] = cur_payoff[1][i] - cur_payoff[2][i];
            x[0][i] = 1.0;
          }

          for (k = 0; k < iNbIndex; k++) {
            memcpy(x[k + 1], save_valuesj[k], nb_paths * sizeof(double));
          }

          err = find_linear_coef(iNbIndex + 1, nb_paths, y, x, b, mat);

          while (err && iNbIndex > 1) {
            iNbIndex--;
            err = find_linear_coef(iNbIndex + 1, nb_paths, y, x, b, mat);
          }

          if (err)
            goto FREE_RETURN;

          /* Rescale the coeff */
          if (params->iRescaleLinCoefs) {
            res = 0.0;
            for (k = 0; k < iNbIndex; k++) {
              res += b[k + 1];
            }

            if (fabs(res) < 1.0E-08) {
              res = 1.0;
            } else {
              res = fabs(1.0 / res);
            }

            for (k = 0; k < iNbIndex; k++) {
              b[k + 1] *= res;
            }
          }

          /* Update the payoff */

          /* add the current PV's */
          for (i = 0; i < nb_paths; i++) {
            cur_payoff[0][i] = 0.0;

            for (k = 0; k < iNbIndex; k++) {
              cur_payoff[0][i] += b[k + 1] * save_valuesj[k][i];
            }

            cur_payoff[3][i] = i;
          }
        } else {
          /* First update the payoff */
          for (i = 0; i < nb_paths; i++) {
            cur_payoff[0][i] = save_values[j][0][i];
            cur_payoff[1][i] =
                save_values[j + add_one_current][params->iNbIndex][i];
            cur_payoff[2][i] = next_payoff[i];
            cur_payoff[3][i] = i;

            if (params->iHasFees) {
              /* Substract Paid Fees */
              cur_payoff[1][i] -=
                  save_valuesj[params->iNbIndex + params->iHasNumeraire + 1][i];
            }
          }

          b[1] = 1.0;
        }

        /* Sort this vector by increasing index */

        num_f_sort_matrix(nb_paths, 4, 0, cur_payoff);

        /* Find the maximum */

        value_max = payoff = sum_next;
        index = nb_paths - 1;

        for (i = nb_paths - 1; i >= 0; i--) {
          payoff += cur_payoff[1][i] - cur_payoff[2][i];

          if (payoff > value_max) {
            value_max = payoff;
            index = i - 1;
          }
        }

        sum_next = value_max;

        if (index >= nb_paths - 1) {
          params->dBarrier[j] = cur_payoff[0][nb_paths - 1] +
                                0.1 * fabs(cur_payoff[0][nb_paths - 1]);
        } else if (index >= 0) {
          params->dBarrier[j] =
              0.5 * (cur_payoff[0][index] + cur_payoff[0][index + 1]);
        } else {
          params->dBarrier[j] =
              cur_payoff[0][0] - 0.000001 * fabs(cur_payoff[0][0]);
        }

        /* Update payoff */
        for (i = index + 1; i < nb_paths; i++) {
          next_payoff[(int)cur_payoff[3][i]] = cur_payoff[1][i];
        }

        if (params->iCalcExeProba) {
          for (i = index + 1; i < nb_paths; i++) {
            call_col[(int)cur_payoff[3][i]] = j;
          }
        }
      } else {
        if (params->iAddNonOptimisedForKI) {
          for (i = 0; i < nb_paths; i++) {
            next_payoff[i] += save_values[j][params->iNbIndex][i];
            sum_next += save_values[j][params->iNbIndex][i];
          }
        }
      }
    }
  }

  if (iDoSmoothing && params->iIsKO) {
    /* Calculate Call Spread */
    for (j = 0; j < nb_dates; j++) {
      dStdCallSpread[j] *= params->dCallSpread;
    }

    for (i = 0; i < nb_paths; i++) {
      dNotional = 1.0;
      dValue = 0.0;

      for (j = 0; j < nb_dates; j++) {
        if (optimise[j]) {
          /* reconstruct Index */
          dIndex = 0.0;

          for (k = 0; k < iNbIndex; k++) {
            dIndex += params->dCoefLin[j][k + 1] * save_values[j][k][i];
          }

          if (!params->iCallCurrent) {
            /* First Add the PV */
            dValue += dNotional * save_values[j][params->iNbIndex][i];
          }

          if (dIndex > params->dBarrier[j] + dStdCallSpread[j]) {
            dNotional = 0.0;
            j = nb_dates;
          } else if (dIndex > params->dBarrier[j] - dStdCallSpread[j]) {
            dNotional *= (params->dBarrier[j] + dStdCallSpread[j] - dIndex) /
                         (2.0 * dStdCallSpread[j]);

            if (params->iCallCurrent) {
              /* Eventually Add the PV */
              dValue += dNotional * save_values[j][params->iNbIndex][i];
            }
          } else {
            if (params->iCallCurrent) {
              /* Eventually Add the PV */
              dValue += dNotional * save_values[j][params->iNbIndex][i];
            }
          }
        } else {
          dValue += dNotional * save_values[j][params->iNbIndex][i];
        }
      }

      next_payoff[i] = dValue;
    }
  }

  res = 0.0;
  res2 = 0.0;

  /* Calculate final result */

  for (i = 0; i < nb_paths; i++) {
    res += next_payoff[i] / nb_paths;
    res2 += next_payoff[i] * next_payoff[i] / nb_paths;
  }

  if (params->iCalcExeProba) {
    for (i = 0; i < nb_paths; i++) {
      if (call_col[i] > 0) {
        params->dExeProba[call_col[i]] += inv_numpaths;
      }
    }
  }

  if (params->iDoInfos) {
    for (j = 0; j < nb_dates; j++) {
      params->dCoefLin[j][1] = 0.0;
    }

    for (j = 0; j < nb_dates; j++) {
      params->dCoefLin[j][params->iNbIndex] /= (nb_paths * 1.0);
    }
  }

  res2 = (res2 - res * res) / nb_paths;

  if (res2 < 1.0E-10) {
    res2 = 0.0;
  } else {
    res2 = sqrt(res2);
  }

  *value = res;
  *error = res2;

  if (stream)
    fclose(stream);

FREE_RETURN:

  if (next_payoff) {
    free_dvector(next_payoff, 0, nb_paths - 1);
  }

  if (cur_payoff) {
    free_dmatrix(cur_payoff, 0, 3 - params->iIsKO, 0, nb_paths - 1);
  }

  if (x) {
    free_dmatrix(x, 0, params->iNbIndex, 0, nb_paths - 1);
  }

  if (y) {
    free_dvector(y, 0, nb_paths - 1);
  }

  if (mat) {
    free_dmatrix(mat, 0, params->iNbIndex, 0, nb_paths - 1);
  }

  if (call_col)
    free(call_col);

  if (total_payoff)
    free(total_payoff);

  if (dStdCallSpread)
    free(dStdCallSpread);

  return err;
}

Err optimise_boundary_info(
    double ***save_values, long nb_dates, long nb_paths, int *optimise,
    int iCallCurrent, int iIsKO, double *boundary, double *value, double *error,
    double **infos) /* Col 1: One time Call        , Col2:
                       Exerc Proba        , Col3: New PV */
{
  double *next_payoff = NULL, **cur_payoff = NULL;

  long *call_date = NULL, *easy_date = NULL;

  double value_max, payoff, sum_next;
  double res, res2;
  double one_time;
  long i, j, index, pathnb;

  // FILE	*stream = fopen("C:\\optim.txt"        , "w+");

  Err err = NULL;

  /* Allocate new variables */
  next_payoff = dvector(0, nb_paths - 1);
  call_date = lvector(0, nb_paths - 1);
  easy_date = lvector(0, nb_paths - 1);
  cur_payoff = dmatrix(0, 3 - iIsKO, 0, nb_paths - 1);

  /*
  for (i=0; i<nb_paths; i++)
  {
          fprintf(stream        , "%f	%f\n"        , save_values[i][1][0] ,
  save_values[i][2][0]);
  }

  fclose(stream);
  */

  if (!next_payoff || !cur_payoff) {
    err = "Memory allocation failure in optimise_boundary";
    goto FREE_RETURN;
  }

  if (iIsKO) {

  } else {
    /* Knock In case */

    /*	Perform the optimisation	*/

    sum_next = 0.0;

    for (j = nb_dates - 1; j >= 0; j--) {
      if (optimise[j]) {
        /* First update the payoff */

        if (iCallCurrent) {
          /* add the current PV's */

          for (i = 0; i < nb_paths; i++) {
            cur_payoff[0][i] = save_values[j][0][i];
            cur_payoff[1][i] = save_values[j][1][i];
            cur_payoff[2][i] = next_payoff[i];
            cur_payoff[3][i] = i;
          }
        } else {
          /* recopy the PV */
          for (i = 0; i < nb_paths; i++) {
            cur_payoff[0][i] = save_values[j][0][i];
            cur_payoff[1][i] = save_values[j + 1][1][i];
            cur_payoff[2][i] = next_payoff[i];
            cur_payoff[3][i] = i;
          }
        }

        /* Sort this vector by increasing index */

        num_f_sort_matrix(nb_paths, 4, 0, cur_payoff);

        /* Find the maximum */

        value_max = payoff = sum_next;
        index = nb_paths - 1;
        one_time = 0.0;

        for (i = nb_paths - 1; i >= 0; i--) {
          payoff += cur_payoff[1][i] - cur_payoff[2][i];

          if (payoff > value_max) {
            value_max = payoff;
            index = i - 1;
          }

          one_time += max(cur_payoff[1][i], 0.0);

          pathnb = (int)cur_payoff[3][i];
        }

        infos[j][0] = one_time;
        infos[j][1] = nb_paths - index - 1;

        sum_next = value_max;

        if (index == nb_paths - 1) {
          boundary[j] = cur_payoff[0][index] * 1.1;
        } else if (index >= 0) {
          boundary[j] = 0.5 * (cur_payoff[0][index] + cur_payoff[0][index + 1]);
        } else {
          boundary[j] = cur_payoff[0][0] * 0.999999;
        }

        /* Update payoff */

        for (i = 0; i <= index; i++) {
          pathnb = (int)cur_payoff[3][i];
          if (cur_payoff[1][i] > 0.0) {
            if (call_date[pathnb] > 0) {
              infos[call_date[pathnb]][3 + nb_dates + j] +=
                  cur_payoff[2][i] - cur_payoff[1][i];

              if (easy_date[pathnb]) {
                infos[call_date[pathnb]][3 + nb_dates + call_date[pathnb]] -=
                    cur_payoff[2][i];
                infos[call_date[pathnb]][3 + call_date[pathnb]] -= 1.0;
                easy_date[pathnb] = 0;
              }
            } else {
              /* spoiled IV */
              infos[j + 1][3 + nb_dates + j] += -cur_payoff[1][i];
            }
          }
        }

        for (i = index + 1; i < nb_paths; i++) {
          pathnb = (int)cur_payoff[3][i];

          if (call_date[pathnb] > 0) {
            infos[call_date[pathnb]][1] -= 1;
            infos[j][3 + call_date[pathnb]] += 1;
            infos[j][3 + nb_dates + call_date[pathnb]] +=
                cur_payoff[1][i] - cur_payoff[2][i];

            if (easy_date[pathnb]) {
              infos[call_date[pathnb]][3 + nb_dates + call_date[pathnb]] -=
                  cur_payoff[2][i];
              infos[call_date[pathnb]][3 + call_date[pathnb]] -= 1.0;
              easy_date[pathnb] = 0;
            }

            infos[call_date[pathnb]][2] -= next_payoff[pathnb];
          } else {
            infos[j][3 + nb_dates + j] += cur_payoff[1][i];
            infos[j][3 + j] += 1.0;
            easy_date[pathnb] = 1;
          }

          call_date[pathnb] = j;
          next_payoff[pathnb] = cur_payoff[1][i];
          infos[j][2] += cur_payoff[1][i];
        }
      } else {
        for (i = 0; i < nb_paths; i++) {
          next_payoff[i] += save_values[j][1][i];
          sum_next += save_values[j][1][i];
        }

        infos[j][0] = 0.0;
        infos[j][1] = 0.0;
        infos[j][2] = 0.0;
      }
    }
  }

  res = 0.0;
  res2 = 0.0;

  /* Calculate final result */

  for (i = 0; i < nb_paths; i++) {
    res += next_payoff[i] / nb_paths;
    res2 += next_payoff[i] * next_payoff[i] / nb_paths;
  }

  res2 = (res2 - res * res) / nb_paths;

  if (res2 < 1.0E-10) {
    res2 = 0.0;
  } else {
    res2 = sqrt(res2);
  }

  for (j = nb_dates - 1; j >= 0; j--) {
    for (i = 2 + 2 * nb_dates; i >= 0; i--) {
      infos[j][i] /= nb_paths;
    }

    /*
    infos[j][3 + nb_dates + j] += infos[j][2];
    */
  }

  *value = res;
  *error = res2;

FREE_RETURN:

  if (next_payoff) {
    free_dvector(next_payoff, 0, nb_paths - 1);
  }

  if (cur_payoff) {
    free_dmatrix(cur_payoff, 0, 3 - iIsKO, 0, nb_paths - 1);
  }

  if (call_date) {
    free_lvector(call_date, 0, nb_paths - 1);
  }

  if (easy_date) {
    free_lvector(easy_date, 0, nb_paths - 1);
  }

  return err;
}

static void multiplic_matrix_vector(long nx, long ny, double **x, double *y,
                                    double *res) {
  int i, j;
  double sum;
  double *xi;

  for (i = 0; i < nx; i++) {
    xi = x[i];
    sum = 0.0;

    for (j = 0; j < ny; j++) {
      sum += xi[j] * y[j];
    }

    res[i] = sum;
  }
}

static void multiplic_with_special(long nx, long ny, double **x, double **y,
                                   double **res) {
  int i, j, k;
  double sum;
  double *xi, *resi;

  for (i = 0; i < nx; i++) {
    xi = x[i + 1];
    resi = res[i];

    for (j = 0; j < ny; j++) {
      sum = 0.0;

      for (k = 0; k < nx; k++) {
        sum += xi[k + 1] * y[k][j];
      }

      resi[j] = sum;
    }
  }
}

static void square_matrix(long n, long m, double **x, double **res) {
  int i, j, k;
  double sum;
  double *xi, *xj, *resi;

  for (i = 0; i < n; i++) {
    xi = x[i];
    resi = res[i];

    for (j = 0; j < n; j++) {
      xj = x[j];
      sum = 0.0;

      for (k = 0; k < m; k++) {
        sum += xi[k] * xj[k];
      }

      resi[j] = sum;
    }
  }
}

/* find b such that Y = Transp(X) * B */
Err find_linear_coef(
    long p, /* number of variables = number of row in Transp(X) */
    long n, /* number of observations */
    double *y, double **x, double *b, double **mat3 /* for memory saving */) {
  double **mat1 = NULL, **mat2 = NULL;

  Err err = NULL;

  mat1 = dmatrix(0, p - 1, 0, p - 1);

  if (!mat1) {
    err = "Memory allocation faillure in find_linear_coef";
    goto FREE_RETURN;
  }

  square_matrix(p, n, x, mat1);

  /* Allocation done inside... */
  mat2 = inverse_matrix(mat1, 0, p - 1);

  if (!mat2) {
    err = "Cannot invert matrix in find_linear_coef";
    goto FREE_RETURN;
  }

  multiplic_with_special(p, n, mat2, x, mat3);
  multiplic_matrix_vector(p, n, mat3, y, b);

FREE_RETURN:

  if (mat1)
    free_dmatrix(mat1, 0, p - 1, 0, p - 1);
  if (mat2)
    free_dmatrix(mat2, 1, p, 1, p);

  return err;
}

int twiddle(int *x, int *y, int *z, int *p) {
  static int i, j, k;

  j = 1;

  while (p[j] <= 0) {
    j++;
  }

  if (p[j - 1] == 0) {
    for (i = j - 1; i > 1; i--) {
      p[i] = -1;
    }

    p[j] = 0;
    *x = *z = 0;
    p[1] = 1;
    *y = j - 1;
  } else {
    if (j > 1)
      p[j - 1] = 0;
    do
      j++;
    while (p[j] > 0);
    k = j - 1;
    i = j;
    while (p[i] == 0)
      p[i++] = -1;
    if (p[i] == -1) {
      p[i] = p[k];
      *z = p[k] - 1;
      *x = i - 1;
      *y = k - 1;
      p[k] = -1;
    } else {
      if (i == p[0])
        return (1);
      else {
        p[j] = p[i];
        *z = p[i] - 1;
        p[i] = 0;
        *x = j - 1;
        *y = i - 1;
      }
    }
  }

  return (0);
}

void inittwiddle(m, n, p) int m, n, *p;
{
  int i;
  p[0] = n + 1;
  for (i = 1; i != n - m + 1; i++)
    p[i] = 0;
  while (i != n + 1) {
    p[i] = i + m - n;
    i++;
  }
  p[n + 1] = -2;
  if (m == 0)
    p[1] = 1;
}

/************************
  Here is a sample use of twiddle() and inittwiddle():
#define N 5
#define M 2
#include "stdio.h"
void main()
  {
  int i        , x        , y        , z        , p[N+2]        , b[N];
  inittwiddle(M        , N        , p);
  for(i = 0; i != N-M; i++)
    {
    b[i] = 0;
    putchar('0');
    }
  while(i != N)
    {
    b[i++] = 1;
    putchar('1');
    }
  putchar('\n');
  while(!twiddle(&x        , &y        , &z        , p))
    {
    b[x] = 1;
    b[y] = 0;
    for(i = 0; i != N; i++)
      putchar(b[i]? '1': '0');
    putchar('\n');
    }
  }
************************/

Err find_best_optim_dates(double ***save_values, long nb_dates, long nb_target,
                          long nb_paths, int *optimise, MCEBPARAMS params,
                          double *value, double *error) {
  Err err = NULL;

  int i, j, x, y, z;
  int *optim_test = NULL, *best_optim = NULL, *p = NULL;

  double res, best_res;

  /* We will first assume that all the dates are optimisation dates */

  optim_test = calloc(nb_dates, sizeof(int));
  best_optim = calloc(nb_dates, sizeof(int));
  p = calloc(nb_dates + 2, sizeof(double));

  if (!optim_test || !p) {
    err = "Memory allocation faillure in find_best_optim_dates";
    goto FREE_RETURN;
  }

  inittwiddle(nb_target, nb_dates, p);

  for (i = 0; i != nb_dates - nb_target; i++) {
    optim_test[i] = 0;
  }

  while (i != nb_dates) {
    optim_test[i++] = 1;
  }

  err = find_and_optimise_boundary(save_values, nb_dates, nb_paths, optim_test,
                                   params, &best_res, error);

  if (err)
    goto FREE_RETURN;

  memcpy(best_optim, optim_test, nb_dates * sizeof(int));

  while (!twiddle(&x, &y, &z, p)) {
    optim_test[x] = 1;
    optim_test[y] = 0;

    err = find_and_optimise_boundary(save_values, nb_dates, nb_paths,
                                     optim_test, params, &res, error);

    if (err)
      goto FREE_RETURN;

    if (res > best_res) {
      best_res = res;
      memcpy(best_optim, optim_test, nb_dates * sizeof(int));
    }
  }

  /* Last call with optim */
  for (i = 0; i < nb_dates; i++) {
    params->dBarrier[i] = 0.0;

    for (j = 0; j <= params->iNbIndex; j++) {
      params->dCoefLin[i][j] = 0.0;
    }
  }

  err = find_and_optimise_boundary(save_values, nb_dates, nb_paths, best_optim,
                                   params, value, error);

FREE_RETURN:

  if (optim_test)
    free(optim_test);
  if (best_optim)
    free(best_optim);
  if (p)
    free(p);

  return err;
}