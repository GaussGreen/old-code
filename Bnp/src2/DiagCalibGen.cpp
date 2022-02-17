
#include "DiagCalibGen.h"
#include "DiagCalibDLM.h"
#include "Math.h"
#include "opfnctns.h"

#define ONE_MONTH 0.083333333
#define MIN_CALIB_TIME 0.03

Err Initialise_CalibParams(int UseJumps, double Precision, int NbIterMax,
                           int RecalibAtEnd, int do_calib, int retry_after_fail,
                           int compute_sensi, double sensi_minshift,
                           CALIBGEN_PARAMS CalibParams) {
  CalibParams->UseJumps = UseJumps;
  CalibParams->NbIterMax = NbIterMax;
  CalibParams->Precision = Precision;
  CalibParams->RecalibAtEnd = RecalibAtEnd;

  CalibParams->res_iter = NULL;
  CalibParams->res_iter = dmatrix(0, NbIterMax - 1, 0, 1);

  if (!CalibParams->res_iter) {
    return "Memory allocation error in Initialise_CalibParams";
  }

  CalibParams->do_calib = do_calib;
  CalibParams->retry_after_fail = retry_after_fail;

  CalibParams->compute_sensi = compute_sensi;
  CalibParams->delta = sensi_minshift;
  CalibParams->sensi = 0.0;

  CalibParams->total_computed = 0;

  return NULL;
}

Free_CalibParams(CALIBGEN_PARAMS CalibParams) {
  if (CalibParams) {
    if (CalibParams->res_iter)
      free_dmatrix(CalibParams->res_iter, 0, CalibParams->NbIterMax - 1, 0, 1);
    CalibParams->res_iter = NULL;
  }
}

double solve_for_next_coef_gen(double **res_iter, int nb_iter,
                               double premium_tgt, double precision,
                               int method) /* 0: linear 1: quadratic */
{
  static double vega, coef, a, b, c, delta, res1, res2;
  static int i, j, k, type;

  if (nb_iter == 2 || method == 0) {
    /* linear interpolation */

    i = 0;

    if (res_iter[0][1] < res_iter[1][1]) {
      while (i < nb_iter && res_iter[i][1] < premium_tgt) {
        i++;
      }
    } else {
      while (i < nb_iter && res_iter[i][1] > premium_tgt) {
        i++;
      }
    }

    if (i > 0) {
      if (i == nb_iter) {
        i = i - 2;
      } else {
        i = i - 1;
      }
    }

    /* interpolation between i and i+1 */
    vega = (res_iter[i + 1][1] - res_iter[i][1]) /
           (res_iter[i + 1][0] - res_iter[i][0]);
    coef = res_iter[i][0] + (premium_tgt - res_iter[i][1]) / vega;

    /* check this is OK */
    for (j = 0; j < nb_iter; j++) {
      if (fabs(res_iter[j][0] - coef) < 1.0E-10 &&
          fabs(res_iter[j][1] - premium_tgt) > precision) {
        /* We have a problem */
        if (j > i + 1) {
          coef = 0.5 * (res_iter[i + 1][0] + res_iter[j][0]);

          for (k = i + 2; k < j; k++) {
            if (fabs(res_iter[k][0] - coef) < 1.0E-10) {
              coef = 0.5 * (res_iter[i][0] + res_iter[i + 1][0]);
              break;
            }
          }
        } else {
          coef = 0.5 * (res_iter[i][0] + res_iter[j][0]);

          for (k = j + 1; k < i; k++) {
            if (fabs(res_iter[k][0] - coef) < 1.0E-10) {
              coef = 0.5 * (res_iter[i][0] + res_iter[i + 1][0]);
              break;
            }
          }
        }

        break;
      }
    }
  } else if (nb_iter > 2) {
    /* quadratic interpolation */

    if (nb_iter == 3) {
      i = 1;
    } else {
      i = 0;

      while (i < nb_iter && res_iter[i][1] < premium_tgt) {
        i++;
      }

      if (i <= 1) {
        i = 1;
      } else if (i >= nb_iter - 1) {
        i = nb_iter - 2;
      } else {
        if ((premium_tgt - res_iter[i][1]) <
            0.5 * (res_iter[i + 1][1] - res_iter[i][1])) {
          i--;
        }
      }
    }

    /* interpolation between i-1  , i and i+1 */

    /* check that numbers are different */
    if (fabs(res_iter[i][0] - res_iter[i - 1][0]) < 1.0E-12 ||
        fabs(res_iter[i][0] - res_iter[i + 1][0]) < 1.0E-12) {
      coef = res_iter[i][0];
      return coef;
    }

    /* check that the function is really monotonic */
    if (fabs(res_iter[i][1] - res_iter[i - 1][1]) < 1.0E-08) {
      if ((premium_tgt - res_iter[i][1]) * (premium_tgt - res_iter[i - 1][1]) <
          0.0) {
        coef = 0.5 * (res_iter[i][0] + res_iter[i - 1][0]);
      } else {
        coef = 0.5 * (res_iter[i][0] + res_iter[i + 1][0]);
      }

      return coef;
    }

    if (fabs(res_iter[i][1] - res_iter[i + 1][1]) < 1.0E-08) {
      if ((premium_tgt - res_iter[i][1]) * (premium_tgt - res_iter[i + 1][1]) <
          0.0) {
        coef = 0.5 * (res_iter[i][0] + res_iter[i + 1][0]);
      } else {
        coef = 0.5 * (res_iter[i][0] + res_iter[i - 1][0]);
      }

      return coef;
    }

    /* solve quadratic equation */

    a = res_iter[i - 1][1] / (res_iter[i - 1][0] - res_iter[i][0]) /
            (res_iter[i - 1][0] - res_iter[i + 1][0]) +
        res_iter[i][1] / (res_iter[i][0] - res_iter[i - 1][0]) /
            (res_iter[i][0] - res_iter[i + 1][0]) +
        res_iter[i + 1][1] / (res_iter[i + 1][0] - res_iter[i - 1][0]) /
            (res_iter[i + 1][0] - res_iter[i][0]);

    if (fabs(a) < 1.0E-12) {
      /* we are linear */
      coef =
          solve_for_next_coef_gen(res_iter, nb_iter, premium_tgt, precision, 0);
      return coef;
    }

    b = res_iter[i - 1][1] / (res_iter[i - 1][0] - res_iter[i][0]) /
            (res_iter[i - 1][0] - res_iter[i + 1][0]) *
            (res_iter[i][0] + res_iter[i + 1][0]) +
        res_iter[i][1] / (res_iter[i][0] - res_iter[i - 1][0]) /
            (res_iter[i][0] - res_iter[i + 1][0]) *
            (res_iter[i - 1][0] + res_iter[i + 1][0]) +
        res_iter[i + 1][1] / (res_iter[i + 1][0] - res_iter[i - 1][0]) /
            (res_iter[i + 1][0] - res_iter[i][0]) *
            (res_iter[i][0] + res_iter[i - 1][0]);
    b = -b;

    c = res_iter[i - 1][1] / (res_iter[i - 1][0] - res_iter[i][0]) /
            (res_iter[i - 1][0] - res_iter[i + 1][0]) * res_iter[i][0] *
            res_iter[i + 1][0] +
        res_iter[i][1] / (res_iter[i][0] - res_iter[i - 1][0]) /
            (res_iter[i][0] - res_iter[i + 1][0]) * res_iter[i - 1][0] *
            res_iter[i + 1][0] +
        res_iter[i + 1][1] / (res_iter[i + 1][0] - res_iter[i - 1][0]) /
            (res_iter[i + 1][0] - res_iter[i][0]) * res_iter[i - 1][0] *
            res_iter[i][0] -
        premium_tgt;

    delta = b * b - 4 * a * c;

    if (delta > 0) {
      delta = sqrt(delta);
      res1 = (-b + delta) / (2 * a);
      res2 = (-b - delta) / (2 * a);

      /* choose the one */
      if (fabs(res2 - res_iter[i][0]) < fabs(res1 - res_iter[i][0])) {
        coef = res2;
      } else {
        coef = res1;
      }
    } else {
      coef =
          solve_for_next_coef_gen(res_iter, nb_iter, premium_tgt, precision, 0);
    }
  } else {
    /* only one point is available */
    if (res_iter[0][1] > premium_tgt) {
      coef = res_iter[0][0] * 0.9;
    } else {
      coef = res_iter[0][0] * 1.1;
    }
  }

  /* check that we are not back on a previous attempt */
  if (method > 0) {
    for (i = 0; i < nb_iter; i++) {
      if (fabs(res_iter[i][0] - coef) < 1.0E-10 &&
          fabs(res_iter[i][1] - premium_tgt) > precision) {
        /* BAD: need another guess */
        coef = solve_for_next_coef_gen(res_iter, nb_iter, premium_tgt,
                                       precision, 0);
        return coef;
      }
    }
  }

  return coef;
}

/* Find the new sensitivity */
void solve_for_next_sensi(double **res_iter, int nb_computed, double value,
                          double delta, double *sensi) {
  int i, index;
  double sens_down, sens_up;
  double has_down, has_up;

  /* Set the Sensitivity */
  if (nb_computed > 1 &&
      (res_iter[nb_computed - 1][0] - res_iter[0][0]) > delta) {

    for (index = 0;
         (index < nb_computed) && (res_iter[index][0] < value - 1.0E-10);
         index++)
      ;

    if (index >= nb_computed) {
      index = index - 1;
    }

    if (index > 0 && fabs(value - res_iter[index - 1][0]) <
                         fabs(value - res_iter[index][0])) {
      index--;
    }

    i = index + 1;

    while (i < nb_computed && (res_iter[i][0] - res_iter[index][0]) < delta) {
      i++;
    }

    if (i < nb_computed) {
      sens_up = (res_iter[i][1] - res_iter[index][1]) /
                (res_iter[i][0] - res_iter[index][0]);
      has_up = 1.0;
    } else {
      sens_up = 0.0;
      has_up = 0.0;
    }

    i = index - 1;

    while (i >= 0 && (res_iter[index][0] - res_iter[i][0]) < delta) {
      i--;
    }

    if (i >= 0) {
      sens_down = (res_iter[i][1] - res_iter[index][1]) /
                  (res_iter[i][0] - res_iter[index][0]);
      has_down = 1.0;
    } else {
      sens_down = 0.0;
      has_down = 0.0;
    }

    if (has_up + has_down > 0.5) {
      *sensi = (sens_down + sens_up) / (has_down + has_up);
    }
  }
}

Err CalibrateNextParam(void *Inst, void *InstConst, void *GlobalConst,
                       void *Model,

                       CALIBGEN_PARAMS CalibConsts,

                       int index_param, double last_param, double limit_down,
                       double limit_up,

                       CALIBFUNCTIONS AllFunctions,

                       double *res_param, int *success) {
  Err err = NULL;
  double param1, param2, price1, price2;
  int j, l, iter, nb_computed;
  double target;
  double **res_iter;

  /* Initialisation */
  if (CalibConsts->do_calib) {
    res_iter = CalibConsts->res_iter;
    CalibConsts->nb_computed = 0;

    AllFunctions->GetTarget(Inst, GlobalConst, Model, CalibConsts, &target);

    if (err)
      goto FREE_RETURN;

    *success = 1;

    /* First Pricing */
    err = AllFunctions->GetFirstGuess(Model, GlobalConst, index_param, target,
                                      &param1);

    if (err)
      goto FREE_RETURN;

    err = AllFunctions->BumpParam(Model, GlobalConst, index_param, param1);

    if (err)
      goto FREE_RETURN;

    err = AllFunctions->UpdateConstsAfterParam(Inst, InstConst, GlobalConst,
                                               Model, CalibConsts);

    if (err)
      goto FREE_RETURN;

    err = AllFunctions->PriceInst(Inst, InstConst, GlobalConst, Model, &price1);

    if (err)
      goto FREE_RETURN;

    res_iter[0][0] = param1;
    res_iter[0][1] = price1;
    CalibConsts->nb_computed = 1;
    nb_computed = 1;

    if (fabs(target - price1) > CalibConsts->Precision) {
      /* Guess for next param */
      if (CalibConsts->compute_sensi && fabs(CalibConsts->sensi) > 0.0) {
        param2 = param1 + (target - price1) / CalibConsts->sensi;
      } else {
        err = AllFunctions->GetSecondGuess(Model, GlobalConst, index_param,
                                           param1, price1, target, &param2);

        if (err)
          goto FREE_RETURN;
      }

      /* Check the bounds */
      if (param2 < limit_down) {
        param2 = limit_down;
      }

      if (param2 > limit_up) {
        param2 = limit_up;
      }

      /* Check that it was not already the case at first step !!! */
      if (fabs(param1 - param2) < 1.0E-10) {
        /* Problem  , stop calibration */
        if (!CalibConsts->UseJumps) {
          param2 = last_param;
        }

        smessage("Calibration: failed to calibrate instrument nb %d",
                 index_param);
        *success = 0;
      }
    } else {
      param2 = param1;
    }

    iter = 1;

    while (fabs(target - price1) > CalibConsts->Precision &&
           iter < CalibConsts->NbIterMax && *success) {
      iter++;

      err = AllFunctions->BumpParam(Model, GlobalConst, index_param, param2);

      if (err)
        goto FREE_RETURN;

      err = AllFunctions->UpdateConstsAfterParam(Inst, InstConst, GlobalConst,
                                                 Model, CalibConsts);

      if (err) {
        if (!CalibConsts->retry_after_fail) {
          goto FREE_RETURN;
        }

        /* try the average */
        param2 = 0.5 * (param1 + param2);
        err = NULL;
      } else {

        err = AllFunctions->PriceInst(Inst, InstConst, GlobalConst, Model,
                                      &price2);

        if (err)
          goto FREE_RETURN;

        nb_computed++;

        param1 = param2;
        price1 = price2;

        /* Save Res */

        l = 0;
        while (l < nb_computed - 1 && res_iter[l][0] < param2) {
          l++;
        }

        if (l < nb_computed - 1) {
          for (j = nb_computed - 2; j >= l; j--) {
            res_iter[j + 1][0] = res_iter[j][0];
            res_iter[j + 1][1] = res_iter[j][1];
          }
        }

        res_iter[l][0] = param2;
        res_iter[l][1] = price2;

        CalibConsts->nb_computed = nb_computed;

        param2 = solve_for_next_coef_gen(res_iter, nb_computed, target,
                                         CalibConsts->Precision, 2);
      }

      /* Check the bounds */
      if (param2 < limit_down) {
        param2 = limit_down;
      }

      if (param2 > limit_up) {
        param2 = limit_up;
      }

      /* Check if this is already the second time */
      if ((res_iter[0][0] < limit_down + 1.0E-12 &&
           param2 < limit_down + 1.0E-12) ||
          (res_iter[iter - 1][0] > limit_up - 1.0E-12 &&
           param2 > limit_up - 1.0E-12)) {
        /* stop this calibration */
        if (CalibConsts->UseJumps) {
          if (param2 - 1.0E-08 < limit_down) {
            param2 = limit_down;
          } else {
            param2 = limit_up;
          }
        } else {
          param2 = last_param;
        }

        smessage("Calibration: failed to calibrate instrument nb %d",
                 index_param);
        *success = 0;
      }
    }

    if (fabs(target - price1) > CalibConsts->Precision &&
        iter >= CalibConsts->NbIterMax) {
      smessage("Calibration: failed to calibrate instrument nb %d   , Max "
               "iteration reached",
               index_param);
      *success = 0;
    }

    if ((*success == 0 && !CalibConsts->UseJumps) ||
        (CalibConsts->RecalibAtEnd && nb_computed > 1)) {
      err = AllFunctions->BumpParam(Model, GlobalConst, index_param, param2);

      if (err)
        goto FREE_RETURN;

      err = AllFunctions->UpdateConstsAfterParam(Inst, InstConst, GlobalConst,
                                                 Model, CalibConsts);

      if (err)
        goto FREE_RETURN;
    }
  } else {
    /* Just update structures */
    param2 = 0.0;

    err = AllFunctions->UpdateConstsAfterParam(Inst, InstConst, GlobalConst,
                                               Model, CalibConsts);
    if (err)
      goto FREE_RETURN;
  }

FREE_RETURN:

  *res_param = param2;

  return err;
}

Err CalibrateParamTS(int StartIndexInst, int EndIndexInst, void **AllInst,
                     void **InstConst, void *GlobalConst, void *Model,

                     CALIBGEN_PARAMS CalibConsts,

                     CALIBFUNCTIONS AllFunctions) {
  Err err = NULL;
  int i, success;
  double limit_down, limit_up, last_param, new_param;

  new_param = 0.0;

  for (i = StartIndexInst; i <= EndIndexInst; i++) {
    err =
        AllFunctions->GetLimitAndLastParam(Model, CalibConsts, GlobalConst, i,
                                           &last_param, &limit_down, &limit_up);

    if (err)
      goto FREE_RETURN;

    err = CalibrateNextParam(AllInst[i], InstConst[i], GlobalConst, Model,
                             CalibConsts, i, last_param, limit_down, limit_up,
                             AllFunctions, &new_param, &success);

    if (err)
      goto FREE_RETURN;

    err = AllFunctions->SetParam(Model, CalibConsts, GlobalConst, i, new_param);

    if (err)
      goto FREE_RETURN;

    if (CalibConsts->compute_sensi) {
      solve_for_next_sensi(CalibConsts->res_iter, CalibConsts->nb_computed,
                           new_param, CalibConsts->delta,
                           &(CalibConsts->sensi));
    }

    CalibConsts->total_computed += CalibConsts->nb_computed;
  }

  AllFunctions->ExtrapolParam(Model, GlobalConst, EndIndexInst);

FREE_RETURN:

  return err;
}