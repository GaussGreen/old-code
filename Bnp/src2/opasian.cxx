/******************************************************************************/
/*                                                                            */
/*       MODULE      :	SRT_F_TMEOPTION.C
 */
/*		 AUTHOR		 :	AL (modif. from asian equity of CY) */
/*       CREATED     :	April 1998
 */
/*                                                                            */
/******************************************************************************/

#include "math.h"
#include "num_h_gamma.h"
#include "opfnctns.h"
#include "srt_h_all.h"

/* -------------------------------------------------------------------------- */
/* solve the system in (F        ,d        ,Sig) knowing (M1        ,M2 ,M3)
   s.t. F+d=M1 F^2*(exp(sig^2*T)-1)=M2-M1^2
                                F^3*(exp(3*sig^2*T)-1)+3*(M2-M1^2)*(M1-F)=M3-M1^3
 */

void find_d(double *ans, double spot, double A1, double A2, double A3) {
  double aux, auxder, v2, v3, x, y, forward_, disp_, tau_;
  int i;

  v2 = A2 - A1 * A1;
  v3 = A3 - 3.0 * v2 * A1 - A1 * A1 * A1;
  y = A2 / A1;

  for (i = 0; i < 10; i++) {
    aux = y * y * y + 3.0 * v2 * y - v3;
    auxder = 3.0 * (y * y + v2);
    y = y - aux / auxder;
  }

  x = 1.0 + y * y / v2;
  tau_ = sqrt((log(1.0 + y * y / v2)));
  forward_ = y / (x - 1.0);
  disp_ = A1 - forward_;

  ans[0] = disp_ * spot;
  ans[1] = forward_ * spot;
  ans[2] = tau_;
}

double thirdmomentformatchingmoments(double F, double M1, double M2) {
  double var, answer;

  var = log(1.0 + (M2 - M1 * M1) / (F * F));
  answer = F * F * F * (exp(3 * var) - 1) + 3 * (M2 - M1 * M1) * (M1 - F) +
           M1 * M1 * M1;
  return answer;
}

/* -------------------------------------------------------------------------- */
/* solve the system in (F        ,d        ,Sig) knowing (M1        ,M2 ,M3)
s.t. F+d=M1 F^2*(exp(sig^2*T)-1)=M2-M1^2
                                F^3*(exp(3*sig^2*T)-1)+3*(M2-M1^2)*(M1-F)=M3-M1^3
Use a dichotomy */
Err srt_f_matching_moments(double *shift, double *varshift, double M1,
                           double M2, double M3) {
  Err err = NULL;
  double Fmin, Fmax, Fmid;
  double M3min, M3max, M3mid;

  Fmin = M1;
  Fmax = 10 * M1;
  M3min = thirdmomentformatchingmoments(Fmin, M1, M2);
  M3max = thirdmomentformatchingmoments(Fmax, M1, M2);
  while (Fmax - Fmin > 5e-6) /* accuracy 0.05bp */
  {
    Fmid = 0.5 * (Fmin + Fmax);
    M3mid = thirdmomentformatchingmoments(Fmid, M1, M2);
    if (M3mid > M3) {
      Fmin = Fmid;
    } else
      Fmax = Fmid;
  }

  *shift = Fmid - M1;
  *varshift = log(1.0 + (M2 - M1 * M1) / (Fmid * Fmid));

  return err;
}

Err srt_f_sumshiftedlognormal(double *Forwards, double *Shifts,
                              double *VolShifts, double *Weights,
                              double **Correl, int Num_assets, double Maturity,
                              double *Forward_bskt, double *Shift_bskt,
                              double *VolShift_bskt) {
  Err err = NULL;

  double **SecondMomentsShifted, ***ThirdMomentsShifted, **SecondMoments,
      ***ThirdMoments;
  double Temp_Fwd_bskt, SecondMoment_bskt, ThirdMoment_bskt;
  double Tempshift_bskt, varshift_bskt;
  int i, j, k;

  /* Allocates memory */
  SecondMomentsShifted = dmatrix(0, Num_assets - 1, 0, Num_assets - 1);
  ThirdMomentsShifted =
      f3tensor(0, Num_assets - 1, 0, Num_assets - 1, 0, Num_assets - 1);
  SecondMoments = dmatrix(0, Num_assets - 1, 0, Num_assets - 1);
  ThirdMoments =
      f3tensor(0, Num_assets - 1, 0, Num_assets - 1, 0, Num_assets - 1);

  /* Calculates the 2 and 3 moments of the (assets+shifts) */
  for (i = 0; i < Num_assets; i++)
    for (j = i; j < Num_assets; j++) {
      SecondMomentsShifted[i][j] =
          (Forwards[i] + Shifts[i]) * (Forwards[j] + Shifts[j]) *
          exp(Correl[i][j] * VolShifts[i] * VolShifts[j] * Maturity);
      SecondMomentsShifted[j][i] = SecondMomentsShifted[i][j];
    }

  for (i = 0; i < Num_assets; i++)
    for (j = i; j < Num_assets; j++)
      for (k = j; k < Num_assets; k++) {
        ThirdMomentsShifted[i][j][k] =
            (Forwards[i] + Shifts[i]) * (Forwards[j] + Shifts[j]) *
            (Forwards[k] + Shifts[k]) *
            exp((Correl[i][j] * VolShifts[i] * VolShifts[j] +
                 Correl[i][k] * VolShifts[i] * VolShifts[k] +
                 Correl[j][k] * VolShifts[j] * VolShifts[k]) *
                Maturity);
        ThirdMomentsShifted[i][k][j] = ThirdMomentsShifted[i][j][k];
        ThirdMomentsShifted[j][i][k] = ThirdMomentsShifted[i][j][k];
        ThirdMomentsShifted[j][k][i] = ThirdMomentsShifted[i][j][k];
        ThirdMomentsShifted[k][i][j] = ThirdMomentsShifted[i][j][k];
        ThirdMomentsShifted[k][j][i] = ThirdMomentsShifted[i][j][k];
      }

  /* Calculates the 2 and 3 moments of the assets */
  for (i = 0; i < Num_assets; i++)
    for (j = i; j < Num_assets; j++) {
      SecondMoments[i][j] = SecondMomentsShifted[i][j] -
                            Shifts[i] * Forwards[j] - Shifts[j] * Forwards[i] -
                            Shifts[i] * Shifts[j];
      SecondMoments[j][i] = SecondMoments[i][j];
    }

  for (i = 0; i < Num_assets; i++)
    for (j = i; j < Num_assets; j++)
      for (k = i; k < Num_assets; k++) {
        ThirdMoments[i][j][k] = ThirdMomentsShifted[i][j][k] -
                                SecondMoments[i][j] * Shifts[k] -
                                SecondMoments[i][k] * Shifts[j] -
                                Forwards[i] * Shifts[j] * Shifts[k] -
                                Shifts[i] * SecondMoments[j][k] -
                                Shifts[i] * Shifts[k] * Forwards[j] -
                                Shifts[i] * Shifts[j] * Forwards[k] -
                                Shifts[i] * Shifts[j] * Shifts[k];
        ThirdMoments[i][k][j] = ThirdMoments[i][j][k];
        ThirdMoments[j][i][k] = ThirdMoments[i][j][k];
        ThirdMoments[j][k][i] = ThirdMoments[i][j][k];
        ThirdMoments[k][i][j] = ThirdMoments[i][j][k];
        ThirdMoments[k][j][i] = ThirdMoments[i][j][k];
      }

  /* Calculates the three first moments of the basket */
  Temp_Fwd_bskt = 0;
  for (i = 0; i < Num_assets; i++)
    Temp_Fwd_bskt += Weights[i] * Forwards[i];

  SecondMoment_bskt = 0;
  for (i = 0; i < Num_assets; i++)
    for (j = 0; j < Num_assets; j++) {
      SecondMoment_bskt += Weights[i] * Weights[j] * SecondMoments[i][j];
    }

  ThirdMoment_bskt = 0;
  for (i = 0; i < Num_assets; i++)
    for (j = 0; j < Num_assets; j++)
      for (k = 0; k < Num_assets; k++) {
        ThirdMoment_bskt +=
            Weights[i] * Weights[j] * Weights[k] * ThirdMoments[i][j][k];
      }

  /* Calculates the shift and the vol shift of the basket by matching moments */
  err = srt_f_matching_moments(&Tempshift_bskt, &varshift_bskt, Temp_Fwd_bskt,
                               SecondMoment_bskt, ThirdMoment_bskt);

  *Forward_bskt = Temp_Fwd_bskt;
  *Shift_bskt = Tempshift_bskt;
  *VolShift_bskt = sqrt(varshift_bskt / Maturity);

  /* Free memory */
  free_dmatrix(SecondMomentsShifted, 0, Num_assets - 1, 0, Num_assets - 1);
  free_f3tensor(ThirdMomentsShifted, 0, Num_assets - 1, 0, Num_assets - 1, 0,
                Num_assets - 1);
  free_dmatrix(SecondMoments, 0, Num_assets - 1, 0, Num_assets - 1);
  free_f3tensor(ThirdMoments, 0, Num_assets - 1, 0, Num_assets - 1, 0,
                Num_assets - 1);

  return err;
}

/* -------------------------------------------------------------------------- */
/* Asian price */
double eval_asian(int ndate, double *date, double *fwds, double *vols,
                  double spot, double maturity, double strike, int npast_fix,
                  double avg_cur, double disc, SrtCallPutType call_put) {
  long i;
  double *forwards, *fforwards, *ffforwards;
  double Mom1, Mom2, Mom3, extra, aux1, aux2;
  double ffMom1, fMom1, fMom2, sffMom1, sfMom1, sfMom2, wf;
  double ans, *sans;

  /*----------------------------*/
  /* computing the first 3 moments */
  forwards = (double *)calloc(ndate, sizeof(double));
  fforwards = (double *)calloc(ndate, sizeof(double));
  ffforwards = (double *)calloc(ndate, sizeof(double));
  sans = (double *)calloc(4, sizeof(double));

  Mom1 = fMom1 = ffMom1 = 0.0;
  Mom2 = fMom2 = 0.0;
  Mom3 = 0.0;
  for (i = 0; i < ndate; i++) {
    forwards[i] = fwds[i] / (npast_fix + ndate) / spot;
    fforwards[i] = forwards[i] * exp(vols[i] * vols[i] * date[i]);
    ffforwards[i] = forwards[i] * exp(2 * vols[i] * vols[i] * date[i]);
    sfMom1 = fMom1;
    sffMom1 = ffMom1;
    sfMom2 = fMom2;
    Mom1 += forwards[i];
    fMom1 += fforwards[i];
    ffMom1 += ffforwards[i];
    Mom2 += 2.0 * sfMom1 * forwards[i];
    Mom2 += forwards[i] * fforwards[i];
    fMom2 += 2.0 * sffMom1 * fforwards[i];
    fMom2 += fforwards[i] * ffforwards[i];
    Mom3 += 3.0 * sfMom2 * forwards[i];
    Mom3 += 3.0 * sffMom1 * forwards[i] * fforwards[i];
    wf = fforwards[i];
    Mom3 += (wf * wf * wf);
  }
  aux1 = Mom1;
  aux2 = Mom2;
  extra = (avg_cur * npast_fix) / (npast_fix + ndate) / spot;
  Mom1 += extra;
  Mom2 += extra * (2.0 * aux1 + extra);
  Mom3 += extra * (3.0 * aux2 + extra * (extra + 3.0 * aux1));

  /* solve the system in (F        ,d        ,Sig) knowing (M1        ,M2 ,M3)
   */
  find_d(sans, spot, Mom1, Mom2, Mom3);
  /* compute the price of the asian */
  ans = srt_f_optblksch(sans[1], strike - sans[0], sans[2], maturity, disc,
                        call_put, PREMIUM);

  free(sans);
  free(forwards);
  free(fforwards);
  free(ffforwards);

  return (ans);
}

Err srt_f_asian(int nforwards, double *forwards_date, double *forwards,
                int nfixing, double *fixing_date, int nvol, double *vol_date,
                double *volat, double spot, double maturity, double strike,
                int npast_fix, double avg_cur, double disc,
                SrtCallPutType call_put, SrtGreekType greek, double *answer) {
  Err err = NULL;
  double eps = 0.01, x;
  double *fwds, *vols;
  int i, nf, nv;

  /* interpolating data on fixing dates */
  fwds = (double *)calloc(nfixing, sizeof(double));
  vols = (double *)calloc(nfixing, sizeof(double));

  nf = nv = 0;
  for (i = 0; i < nfixing; i++) {
    while ((nf < nforwards) && (forwards_date[nf] < fixing_date[i]))
      nf++;
    if (nf == 0)
      fwds[i] = forwards[0];
    else if (nf == nforwards)
      fwds[i] = forwards[nforwards - 1];
    else {
      x = (fixing_date[i] - forwards_date[nf - 1]) /
          (forwards_date[nf] - forwards_date[nf - 1]);
      fwds[i] = x * forwards[nf] + (1.0 - x) * forwards[nf - 1];
    }

    while ((nv < nvol) && (vol_date[nv] < fixing_date[i]))
      nv++;
    if (nv == 0)
      vols[i] = volat[0];
    else if (nf == nforwards)
      vols[i] = volat[nforwards - 1];
    else {
      x = (fixing_date[i] - vol_date[nv - 1]) /
          (vol_date[nv] - vol_date[nv - 1]);
      vols[i] = x * volat[nv] + (1.0 - x) * volat[nv - 1];
    }
  }

  /*------------------------------------------*/
  switch (greek) {
  case PREMIUM:

    *answer = eval_asian(nfixing, fixing_date, fwds, vols, spot, maturity,
                         strike, npast_fix, avg_cur, disc, call_put);
    break;

  default:
    *answer = 0.0;
    break;
  }

  free(fwds);
  free(vols);

  return (err);
}

/* Call/Put price using Gamma Inverse function */
double srt_f_GammaInvOption(double dStrike, double dM1, double dM2,
                            SrtCallPutType call_put, SrtGreekType greek) {
  return srt_f_GammaInvOpt(dStrike, (2.0 * dM2 - dM1 * dM1) / (dM2 - dM1 * dM1),
                           (dM2 - dM1 * dM1) / dM2 / dM1, call_put, greek);
}

double srt_f_GammaInvOpt(double strike, double a, double b,
                         SrtCallPutType call_put, SrtGreekType greek) {
  double dInvStrike = .0;
  if ((a != 1.0) && (strike > 0) && (greek == PREMIUM)) {
    dInvStrike = 1.0 / strike / b; /* Gamma function scaling */
    if (call_put == SRT_CALL)
      return gammp(a - 1.0, dInvStrike) / b / (a - 1.0) -
             strike * gammp(a, dInvStrike);
    else
      return strike * (1.0 - gammp(a, dInvStrike)) -
             (1.0 - gammp(a - 1.0, dInvStrike)) / b / (a - 1.0);
  } else
    return .0;
}

/* Call/Put price using Gamma Inverse function */
double srt_f_GammaOption(double dStrike, double dM1, double dM2,
                         SrtCallPutType call_put, SrtGreekType greek) {
  return srt_f_GammaOpt(dStrike, dM1 * dM1 / (dM2 - dM1 * dM1),
                        (dM2 - dM1 * dM1) / dM1, call_put, greek);
}

double srt_f_GammaOpt(double strike, double a, double b,
                      SrtCallPutType call_put, SrtGreekType greek) {
  double dStrike = .0;
  if ((a != 1.0) && (strike > 0) && (greek == PREMIUM)) {
    dStrike = strike / b; /* Gamma function scaling */
    if (call_put == SRT_CALL)
      return a * b * (1.0 - gammp(a + 1.0, dStrike)) -
             strike * (1.0 - gammp(a, dStrike));
    else
      return strike * gammp(a, dStrike) - a * b * gammp(a + 1.0, dStrike);
  } else
    return .0;
}

/* Two moments and three moments matching code for Asian option
        In: Cms Start Dates (in the correct basis (Date - today)/365 )        ,
                Cms forward and Cms Vol        ,
                Correl as a vector (since the difference between two dates
   should be constant and Correl[0] = 1 and Correl[i] = <w1        ,w2> iCms the
   number of Cms        ,
                ...
                iMethod == 0 for the second moment matching method
                                == 1 for the third moment matching method
                                == 2 for the gamma price
                                == 3 for the inv gamma price
        The result should be disounted from Maturity to today */

Err srt_f_AsianCapPrice(double *dStartDates, double *dCmsForwards,
                        double *dCmsVols, double *dCorrels, int iCms,
                        double dStrike, SrtCallPutType SrtCallPut,
                        SrtGreekType SrtGreek, SrtDiffusionType SrtVolType,
                        int iMethod, double *dAnswer) {
  Err err = NULL;
  double dM1 = .0, dM2 = .0, dM21 = .0, dM22 = .0, dM3 = .0, *dNewton, dSub22,
         dSigma, dForward;
  int i, j, k;

  /* In Case there are no instruments ... */
  if (iCms == 0) {
    *dAnswer = .0;
    return err;
  }

  /* Greek Switch */
  switch (SrtGreek) {
  case PREMIUM:

    for (i = 0; i < iCms; i++) {
      /* M1 is the first moment */
      dM1 += dCmsForwards[i];
      dM21 += dCmsForwards[i] * dCmsForwards[i] *
              exp(dCmsVols[i] * dCmsVols[i] * dStartDates[i]);
      dSub22 = .0;
      for (j = i + 1; j < iCms; j++)
        dSub22 += dCmsForwards[j] * exp(dCorrels[j - i] * dCmsVols[i] *
                                        dCmsVols[j] * dStartDates[i]);

      dM22 += dCmsForwards[i] * dSub22;

      if (iMethod == 3)
        for (j = 0; j < iCms; j++)
          for (k = 0; k < iCms; k++)
            dM3 += dCmsForwards[i] * dCmsForwards[j] * dCmsForwards[k] *
                   exp(dCorrels[abs(j - i)] * dCmsVols[i] * dCmsVols[j] *
                           min(dStartDates[i], dStartDates[j]) +
                       dCorrels[abs(j - k)] * dCmsVols[k] * dCmsVols[j] *
                           min(dStartDates[k], dStartDates[j]) +
                       dCorrels[abs(i - k)] * dCmsVols[i] * dCmsVols[k] *
                           min(dStartDates[i], dStartDates[k]));
    }

    dM1 = dM1 / iCms;
    dM2 = (dM21 + 2.0 * dM22) / iCms / iCms;

    if (iMethod == 3) {
      dM3 /= iCms * iCms * iCms;
      dNewton = (double *)calloc(4, sizeof(double));
      find_d(dNewton, 1, dM1, dM2, dM3);
      dStrike -= dNewton[0];
      dForward = dNewton[1];
      dSigma = dNewton[2];
    } else {
      dSigma = sqrt(log(dM2) - 2.0 * log(dM1));
      dForward = dM1;
    }

    /* The option maturity is equal to 1 since it's include in the matching
     * moments */
    if (SrtVolType == SRT_NORMAL)
      *dAnswer = srt_f_optblknrm(dForward, dStrike, dSigma, 1.0, 1.0,
                                 SrtCallPut, SrtGreek);
    else {
      if (iMethod == 1)
        *dAnswer = srt_f_GammaOption(dStrike, dM1, dM2, SrtCallPut, PREMIUM);
      if (iMethod == 2)
        *dAnswer = srt_f_GammaInvOption(dStrike, dM1, dM2, SrtCallPut, PREMIUM);
      if ((iMethod == 0) || (iMethod == 3))
        *dAnswer = srt_f_optblksch(dForward, dStrike, dSigma, 1.0, 1.0,
                                   SrtCallPut, SrtGreek);
    }
    break;

  default:
    *dAnswer = 0.0;
    break;
  }
  return err;
}

/*--------------------------------- End of File
 * -------------------------------------*/