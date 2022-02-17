/* Calibrates the SABR parameters to a market smile */

// added by Pierre not so sure
#include <opsabrgeneric.h"

#include "utallhdr.h"
#include <FXSabrAdi.h"
#include <OPFNCTNS.H>
#include <OPSABRCALIB.H>
#include <OPSABRCALIBBVM.H>
#include <math.h"
#include <opsabrgeneric.H>

Err srt_f_optbvmsabrvol(double Forward, double Strike, double Maturity,
                        double Sigma, double Alpha, double Delay, double Rho,
                        SrtDiffusionType input, SrtDiffusionType output,
                        double *vol) {
  double res;
  double SigmaBeta;
  Err err = NULL;

  /* Control of passed arguments */
  if (Forward < 0.0) {
    err = "Fatal: Forward should be positive";
    goto FREE_RETURN;
  }

  if (Sigma < 0.0) {
    err = "Fatal: SigmaBeta should be positive";
    goto FREE_RETURN;
  }

  if (Alpha < 0.0) {
    err = "Fatal: Alpha should be positive";
    goto FREE_RETURN;
  }

  if (Maturity < 0.0) {
    err = "Error: Maturity should better be positive";
    goto FREE_RETURN;
  }

  if ((Rho < -1.0) || (Rho > 1.0)) {
    err = "Error: Rho should be > -1 and < 1";
    goto FREE_RETURN;
  }

  if (input == SRT_BETAVOL) {
    res = op_sabrgen(Forward, Strike, Maturity, Sigma, Alpha, Delay, 0, 0, Rho,
                     vol_bvm);

    if (output == SRT_NORMAL) {
      err = srt_f_optsarbvol(Forward, Strike, Maturity, res, 0, 1, 0,
                             SRT_BETAVOL, SRT_NORMAL, &res);
      if (err) {
        goto FREE_RETURN;
      }
    }
  } else {
    if (input == SRT_NORMAL) {
      /* first convert to lognormal */
      err = srt_f_optsarbvol(Forward, Strike, Maturity, Sigma, 0, 0, 0,
                             SRT_BETAVOL, SRT_LOGNORMAL, &Sigma);
      if (err) {
        goto FREE_RETURN;
      }
    }

    if (output == SRT_LOGNORMAL) {
      if (fabs(Strike - Forward) < 1.0E-10) {
        res = Sigma;
      } else {
        SigmaBeta = op_sabrgen_calib(Forward, Forward, Maturity, Sigma, Alpha,
                                     Delay, 0, 0, Rho, vol_bvm);

        res = op_sabrgen(Forward, Strike, Maturity, SigmaBeta, Alpha, Delay, 0,
                         0, Rho, vol_bvm);
      }
    } else {
      res = op_sabrgen_calib(Forward, Strike, Maturity, Sigma, Alpha, Delay, 0,
                             0, Rho, vol_bvm);
    }
  }

  *vol = res;

FREE_RETURN:

  return err;
}

Err srt_f_optbvmsabr_mr_vol(double Forward, double Strike, double Maturity,
                            double Sigma, double Alpha, double Delay,
                            double Rho, double Lambda, double *numericalparams,
                            int nparams, SrtDiffusionType input,
                            SrtDiffusionType output, double *vol) {
  Err err = NULL;
  double modAlpha;

  if (Lambda < 0.0) {
    return serror("Lambda has to be positive or null");
  } else if (Lambda == 0) {
    modAlpha = Alpha;
  } else {
    modAlpha = Alpha * sqrt((1 - exp(-2 * Lambda * Maturity)) /
                            (2 * Lambda * Maturity));
  }

  err = srt_f_optbvmsabrvol(Forward, Strike, Maturity, Sigma, modAlpha, Delay,
                            Rho, input, output, vol);

  return err;
}

Err HardLimitStructurePrice(double forward,
                            /*	Option specs	*/
                            double strike, double maturity,
                            double putSpreadWidth,
                            /*	Model specs	*/
                            double volinput, double Alpha, double Beta,
                            double Rho,

                            SrtDiffusionType input, SrtDiffusionType output,

                            double *price) {
  Err err = NULL;
  double volK, volKdK;
  double priceK, priceKdK;

  err = srt_f_optsarbvol(forward, strike, maturity, volinput, Alpha, Beta, Rho,
                         input, output, &volK);
  priceK = srt_f_optblksch(forward, strike, volK, maturity, 1.0, SRT_PUT,
                           SRT_PREMIUM);

  err = srt_f_optsarbvol(forward, strike + putSpreadWidth, maturity, volinput,
                         Alpha, Beta, Rho, input, output, &volKdK);
  priceKdK = srt_f_optblksch(forward, strike + putSpreadWidth, volKdK, maturity,
                             1.0, SRT_PUT, SRT_PREMIUM);

  *price = ((priceKdK - priceK) / putSpreadWidth - priceK / strike) / strike;

  return err;
}

Err HardLimitStructurePriceGen(
    double forward,
    /*	Option specs	*/
    double strike, double maturity, double putSpreadWidth,
    /*	Model specs	*/
    double volinput, double Alpha, double a, double b, double c, double Rho,
    double (*vol_local)(double x, double a, double b, double c, int type),

    SrtDiffusionType input, SrtDiffusionType output,

    double *price) {
  Err err = NULL;
  double volK, volKdK;
  double priceK, priceKdK;

  volK = op_sabrgen(forward, strike, maturity, volinput, Alpha, a, b, c, Rho,
                    vol_local);
  priceK = srt_f_optblksch(forward, strike, volK, maturity, 1.0, SRT_PUT,
                           SRT_PREMIUM);

  volKdK = op_sabrgen(forward, strike + putSpreadWidth, maturity, volinput,
                      Alpha, a, b, c, Rho, vol_local);
  priceKdK = srt_f_optblksch(forward, strike + putSpreadWidth, volKdK, maturity,
                             1.0, SRT_PUT, SRT_PREMIUM);

  *price = ((priceKdK - priceK) / putSpreadWidth - priceK / strike) / strike;

  return err;
}

Err alphaHardLimit(double forward,
                   /*	Option specs	*/
                   double strike, double maturity, double putSpreadWidth,
                   /*	Model specs	*/
                   double volinput, double Alpha, double Beta, double Rho,

                   SrtDiffusionType input,

                   double *AlphaHardLimit) {
  Err err = NULL;
  int compt;
  double Alpha_opt, Alpha1, Alpha2, Alpha_max;
  double price1, price2;

  Alpha_max = 1.0;
  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha_max, Beta, Rho, input,
                                SRT_LOGNORMAL, &price1);
  if (price1 > 0) {
    Alpha_opt = Alpha_max;
    *AlphaHardLimit = Alpha_max;
    goto FREE_RETURN;
  }

  Alpha1 = 0.2;
  Alpha_opt = Alpha1;
  Alpha2 = Alpha1 + 0.05;

  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha1, Beta, Rho, input,
                                SRT_LOGNORMAL, &price1);
  if (fabs(price1) < 1e-8) {
    Alpha_opt = Alpha1;
    *AlphaHardLimit = Alpha1;
    goto FREE_RETURN;
  }

  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha2, Beta, Rho, input,
                                SRT_LOGNORMAL, &price2);
  if (fabs(price2) < 1e-8) {
    Alpha_opt = Alpha2;
    *AlphaHardLimit = Alpha2;
    goto FREE_RETURN;
  }

  if (fabs(price2 - price1) < 1e-12) {
    Alpha_opt = Alpha1;
    *AlphaHardLimit = Alpha1;
    goto FREE_RETURN;
  }

  Alpha1 =
      DMIN(Alpha_max, Alpha1 - price1 * (Alpha2 - Alpha1) / (price2 - price1));
  if (Alpha1 < 0) {
    Alpha1 = Alpha_opt + 0.1;
  }
  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha1, Beta, Rho, input,
                                SRT_LOGNORMAL, &price1);
  if (fabs(price1) < 1e-8) {
    Alpha_opt = Alpha1;
    *AlphaHardLimit = Alpha1;
    goto FREE_RETURN;
  }

  compt = 0;
  while ((compt < 20) && (fabs(price1) > 1e-8)) {
    Alpha2 = Alpha1 + 0.05;

    err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                  volinput, Alpha2, Beta, Rho, input,
                                  SRT_LOGNORMAL, &price2);
    if (fabs(price2) < 1e-8) {
      Alpha_opt = Alpha2;
      compt = 20;
      *AlphaHardLimit = Alpha2;
      goto FREE_RETURN;
    }

    if (fabs(price2 - price1) < 1e-12) {
      Alpha_opt = Alpha1;
      compt = 20;
      *AlphaHardLimit = Alpha1;
      goto FREE_RETURN;
    }

    Alpha1 = DMIN(Alpha_max,
                  Alpha1 - price1 * (Alpha2 - Alpha1) / (price2 - price1));
    if (Alpha1 < 0) {
      Alpha1 = Alpha_opt + 0.1;
    }
    err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                  volinput, Alpha1, Beta, Rho, input,
                                  SRT_LOGNORMAL, &price1);
    if (fabs(price1) < 1e-8) {
      compt = 20;
      Alpha_opt = Alpha1;
      *AlphaHardLimit = Alpha1;
      goto FREE_RETURN;
    }

    if (fabs(Alpha_opt - Alpha1) < 1e-4) {
      Alpha_opt = Alpha1;
      compt = 20;
      *AlphaHardLimit = Alpha1;
      goto FREE_RETURN;
    }
    compt++;
    Alpha_opt = Alpha1;
  }

  *AlphaHardLimit = Alpha1;

FREE_RETURN:

  return err;
}

Err alphaHardLimitDicho(double forward,
                        /*	Option specs	*/
                        double strike, double maturity, double putSpreadWidth,
                        /*	Model specs	*/
                        double volinput, double Alpha, double Beta, double Rho,

                        SrtDiffusionType input,

                        double *AlphaHardLimit) {
  Err err = NULL;
  int compt;
  double Alpha_opt, Alpha1, Alpha2;
  double price, price1, price2;

  Alpha1 = 1.0;
  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha1, Beta, Rho, input,
                                SRT_LOGNORMAL, &price1);
  if (price1 > 0) {
    Alpha_opt = Alpha1;
    *AlphaHardLimit = Alpha1;
    goto FREE_RETURN;
  }

  Alpha2 = 0.01;
  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha2, Beta, Rho, input,
                                SRT_LOGNORMAL, &price2);
  if (price2 < 0) {
    Alpha_opt = Alpha2;
    *AlphaHardLimit = Alpha2;
    goto FREE_RETURN;
  }

  price = price1;
  Alpha_opt = Alpha1;
  compt = 0;
  while ((compt < 20) && (fabs(price) > 1e-8)) {
    Alpha_opt = 0.5 * (Alpha1 + Alpha2);

    err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                  volinput, Alpha_opt, Beta, Rho, input,
                                  SRT_LOGNORMAL, &price);
    if (fabs(price) < 1e-8) {
      compt = 20;
      *AlphaHardLimit = Alpha_opt;
      goto FREE_RETURN;
    }

    if (price > 0) {
      Alpha2 = Alpha_opt;
    } else {
      Alpha1 = Alpha_opt;
    }

    compt++;
  }

  *AlphaHardLimit = Alpha_opt;

FREE_RETURN:

  return err;
}

Err betaHardLimitDicho(double forward,
                       /*	Option specs	*/
                       double strike, double maturity, double putSpreadWidth,
                       /*	Model specs	*/
                       double volinput, double Alpha, double Beta, double Rho,

                       SrtDiffusionType input,

                       double *BetaHardLimit) {
  Err err = NULL;
  int compt;
  double Beta_opt, Beta1, Beta2;
  double price, price1, price2;

  Beta1 = 0.0;
  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha, Beta1, Rho, input,
                                SRT_LOGNORMAL, &price1);
  if (price1 > 0) {
    Beta_opt = Beta1;
    *BetaHardLimit = Beta1;
    goto FREE_RETURN;
  }

  Beta2 = 1.0;
  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha, Beta2, Rho, input,
                                SRT_LOGNORMAL, &price2);
  if (price2 < 0) {
    Beta_opt = Beta2;
    *BetaHardLimit = Beta2;
    goto FREE_RETURN;
  }

  price = price1;
  Beta_opt = Beta1;
  compt = 0;
  while ((compt < 20) && (fabs(price) > 1e-8)) {
    Beta_opt = 0.5 * (Beta1 + Beta2);

    err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                  volinput, Alpha, Beta_opt, Rho, input,
                                  SRT_LOGNORMAL, &price);
    if (fabs(price) < 1e-8) {
      compt = 20;
      *BetaHardLimit = Beta_opt;
      goto FREE_RETURN;
    }

    if (price > 0) {
      Beta2 = Beta_opt;
    } else {
      Beta1 = Beta_opt;
    }

    compt++;
  }

  *BetaHardLimit = Beta_opt;

FREE_RETURN:

  return err;
}

Err rhoHardLimitDicho(double forward,
                      /*	Option specs	*/
                      double strike, double maturity, double putSpreadWidth,
                      /*	Model specs	*/
                      double volinput, double Alpha, double Beta, double Rho,

                      SrtDiffusionType input,

                      double *RhoHardLimit) {
  Err err = NULL;
  int compt;
  double Rho_opt, Rho1, Rho2;
  double price, price1, price2;

  Rho1 = -1.0;
  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha, Beta, Rho1, input,
                                SRT_LOGNORMAL, &price1);
  if (price1 > 0) {
    Rho_opt = Rho1;
    *RhoHardLimit = Rho1;
    goto FREE_RETURN;
  }

  Rho2 = 1.0;
  err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                volinput, Alpha, Beta, Rho2, input,
                                SRT_LOGNORMAL, &price2);
  if (price2 < 0) {
    Rho_opt = Rho2;
    *RhoHardLimit = Rho2;
    goto FREE_RETURN;
  }

  price = price1;
  Rho_opt = Rho1;
  compt = 0;
  while ((compt < 20) && (fabs(price) > 1e-8)) {
    Rho_opt = 0.5 * (Rho1 + Rho2);

    err = HardLimitStructurePrice(forward, strike, maturity, putSpreadWidth,
                                  volinput, Alpha, Beta, Rho_opt, input,
                                  SRT_LOGNORMAL, &price);
    if (fabs(price) < 1e-8) {
      compt = 20;
      *RhoHardLimit = Rho_opt;
      goto FREE_RETURN;
    }

    if (price > 0) {
      Rho2 = Rho_opt;
    } else {
      Rho1 = Rho_opt;
    }

    compt++;
  }

  *RhoHardLimit = Rho_opt;

FREE_RETURN:

  return err;
}

// added by Pierre Trying to make it generic

Err SABRGenHardLimitDichoGen(
    double forward,
    /*	Option specs	*/
    double strike, double maturity, double putSpreadWidth,
    /*	Model specs	*/
    double volinput, double Alpha, double a, double b, double c, double Rho,
    double (*vol_local)(double x, double a, double b, double c, int type),
    SrtDiffusionType input,
    int HardLimType, // 0:Alpha  , 1:a  ,2:b  , 3:c  , 4:Rho
    double *HardLimit) {
  Err err = NULL;
  int compt;
  double Params1[5], Params2[5], ParamsOpt[5];
  double price, price1, price2;

  Params1[0] = Alpha;
  Params1[1] = a;
  Params1[2] = b;
  Params1[3] = c;
  Params1[4] = Rho;

  Params2[0] = Alpha;
  Params2[1] = a;
  Params2[2] = b;
  Params2[3] = c;
  Params2[4] = Rho;

  ParamsOpt[0] = Alpha;
  ParamsOpt[1] = a;
  ParamsOpt[2] = b;
  ParamsOpt[3] = c;
  ParamsOpt[4] = Rho;

  // assign the lower and upper bounf before starting the binomial search
  // reads lower and upper bounds from the vol_loc function
  //		Params1[HardLimType] = vol_local(0.01  ,Alpha  ,a  ,b  ,c  ,Rho
  //,10+HardLimType); 		Params2[HardLimType] = vol_local(0.01  ,Alpha  ,a  ,b  ,c
  //,Rho  ,20+HardLimType);

  // compute value at end of the range
  err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                   volinput, Params1[0], Params1[1], Params1[2],
                                   Params1[3], Params1[4], vol_local, input,
                                   SRT_LOGNORMAL, &price1);

  err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                   volinput, Params2[0], Params2[1], Params2[2],
                                   Params2[3], Params2[4], vol_local, input,
                                   SRT_LOGNORMAL, &price2);

  // Case Solution is on boundary and everything is working
  if ((price1 > 0) && (price2 > 0)) {
    *HardLimit = -1000;
    goto FREE_RETURN;
  }
  // Case Solution is on boundary and screwed
  else if ((price1 < 0) && (price2 < 0)) {
    *HardLimit = -1000000;
    goto FREE_RETURN;
  }
  // Case Solution is in-between
  else {

    ParamsOpt[HardLimType] = Params1[HardLimType];
    *HardLimit = ParamsOpt[HardLimType];

    price = -1; // force at least oen iteration

    while ((compt < 20) && (fabs(price) > 1e-8)) {
      ParamsOpt[HardLimType] =
          0.5 * (Params1[HardLimType] + Params2[HardLimType]);

      err = HardLimitStructurePriceGen(
          forward, strike, maturity, putSpreadWidth, volinput, ParamsOpt[0],
          ParamsOpt[1], ParamsOpt[2], ParamsOpt[3], ParamsOpt[4], vol_local,
          input, SRT_LOGNORMAL, &price);
      if (fabs(price) < 1e-8) {
        compt = 20;
        *HardLimit = ParamsOpt[HardLimType];
        goto FREE_RETURN;
      }

      if (price2 * price < 0)
        Params2[HardLimType] = ParamsOpt[HardLimType];
      else
        Params1[HardLimType] = ParamsOpt[HardLimType];

      compt++;
    }

    *HardLimit = ParamsOpt[HardLimType];
  }

FREE_RETURN:

  return err;
}

Err alphaHardLimitDichoGen(
    double forward,
    /*	Option specs	*/
    double strike, double maturity, double putSpreadWidth,
    /*	Model specs	*/
    double volinput, double Alpha, double delay, double b, double c, double Rho,
    double (*vol_local)(double x, double a, double b, double c, int type),

    SrtDiffusionType input,

    double *AlphaHardLimit) {
  Err err = NULL;
  int compt;
  double Alpha_opt, Alpha1, Alpha2;
  double price, price1, price2;

  Alpha1 = 1.0;
  err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                   volinput, Alpha1, delay, b, c, Rho,
                                   vol_local, input, SRT_LOGNORMAL, &price1);
  if (price1 > 0) {
    Alpha_opt = Alpha1;
    *AlphaHardLimit = Alpha1;
    goto FREE_RETURN;
  }

  Alpha2 = 0.01;
  err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                   volinput, Alpha2, delay, b, c, Rho,
                                   vol_local, input, SRT_LOGNORMAL, &price2);
  if (price2 < 0) {
    Alpha_opt = Alpha2;
    *AlphaHardLimit = Alpha2;
    goto FREE_RETURN;
  }

  price = price1;
  Alpha_opt = Alpha1;
  compt = 0;
  while ((compt < 20) && (fabs(price) > 1e-8)) {
    Alpha_opt = 0.5 * (Alpha1 + Alpha2);

    err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                     volinput, Alpha_opt, delay, b, c, Rho,
                                     vol_local, input, SRT_LOGNORMAL, &price);
    if (fabs(price) < 1e-8) {
      compt = 20;
      *AlphaHardLimit = Alpha_opt;
      goto FREE_RETURN;
    }

    if (price > 0) {
      Alpha2 = Alpha_opt;
    } else {
      Alpha1 = Alpha_opt;
    }

    compt++;
  }

  *AlphaHardLimit = Alpha_opt;

FREE_RETURN:

  return err;
}

Err delayHardLimitDichoGen(
    double forward,
    /*	Option specs	*/
    double strike, double maturity, double putSpreadWidth,
    /*	Model specs	*/
    double volinput, double Alpha, double delay, double b, double c, double Rho,
    double (*vol_local)(double x, double a, double b, double c, int type),

    SrtDiffusionType input,

    double *delayHardLimit) {
  Err err = NULL;
  int compt;
  double delay_opt, delay1, delay2;
  double price, price1, price2;

  delay1 = 10000;
  err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                   volinput, Alpha, delay1, b, c, Rho,
                                   vol_local, input, SRT_LOGNORMAL, &price1);
  if (price1 > 0) {
    delay_opt = delay1;
    *delayHardLimit = delay1;
    goto FREE_RETURN;
  }

  delay2 = 0.0001;
  err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                   volinput, Alpha, delay2, b, c, Rho,
                                   vol_local, input, SRT_LOGNORMAL, &price2);
  if (price2 < 0) {
    delay_opt = delay2;
    *delayHardLimit = delay2;
    goto FREE_RETURN;
  }

  price = price1;
  delay_opt = delay1;
  compt = 0;
  while ((compt < 20) && (fabs(price) > 1e-8)) {
    delay_opt = 0.5 * (delay1 + delay2);

    err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                     volinput, Alpha, delay_opt, b, c, Rho,
                                     vol_local, input, SRT_LOGNORMAL, &price);
    if (fabs(price) < 1e-8) {
      compt = 20;
      *delayHardLimit = delay_opt;
      goto FREE_RETURN;
    }

    if (price > 0) {
      delay2 = delay_opt;
    } else {
      delay1 = delay_opt;
    }

    compt++;
  }

  *delayHardLimit = delay_opt;

FREE_RETURN:

  return err;
}

Err rhoHardLimitDichoGen(double forward,
                         /*	Option specs	*/
                         double strike, double maturity, double putSpreadWidth,
                         /*	Model specs	*/
                         double volinput, double Alpha, double delay, double b,
                         double c, double Rho,
                         double (*vol_local)(double x, double a, double b,
                                             double c, int type),

                         SrtDiffusionType input,

                         double *RhoHardLimit) {
  Err err = NULL;
  int compt;
  double Rho_opt, Rho1, Rho2;
  double price, price1, price2;

  Rho1 = -1.0;
  err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                   volinput, Alpha, delay, b, c, Rho1,
                                   vol_local, input, SRT_LOGNORMAL, &price1);
  if (price1 > 0) {
    Rho_opt = Rho1;
    *RhoHardLimit = Rho1;
    goto FREE_RETURN;
  }

  Rho2 = 1.0;
  err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                   volinput, Alpha, delay, b, c, Rho2,
                                   vol_local, input, SRT_LOGNORMAL, &price2);
  if (price2 < 0) {
    Rho_opt = Rho2;
    *RhoHardLimit = Rho2;
    goto FREE_RETURN;
  }

  price = price1;
  Rho_opt = Rho1;
  compt = 0;
  while ((compt < 20) && (fabs(price) > 1e-8)) {
    Rho_opt = 0.5 * (Rho1 + Rho2);

    err = HardLimitStructurePriceGen(forward, strike, maturity, putSpreadWidth,
                                     volinput, Alpha, delay, b, c, Rho_opt,
                                     vol_local, input, SRT_LOGNORMAL, &price);
    if (fabs(price) < 1e-8) {
      compt = 20;
      *RhoHardLimit = Rho_opt;
      goto FREE_RETURN;
    }

    if (price > 0) {
      Rho2 = Rho_opt;
    } else {
      Rho1 = Rho_opt;
    }

    compt++;
  }

  *RhoHardLimit = Rho_opt;

FREE_RETURN:

  return err;
}

/* Computes the value and the gradient of the SABR price at a given strike */

Err op_bvmsabr_pricing(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double Alpha, double Delay, double Rho,
    /*	Vol input	*/
    double atmvol,
    /*	Premium	*/
    double *prem) {
  double Betavol;

  Betavol = op_sabrgen_calib(forward, forward, maturity, atmvol, Alpha, Delay,
                             0, 0, Rho, vol_bvm);

  *prem = op_sabrgen(forward, strike, maturity, Betavol, Alpha, Delay, 0, 0,
                     Rho, vol_bvm);

  *prem = srt_f_optblksch(forward, strike, *prem, maturity, 1, SRT_CALL,
                          SRT_PREMIUM);

  return NULL;
}

/*	Sensitivity to model parameters	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_bvmsabr_model_sens(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double Alpha, double delay, double Rho,
    /*	Vol input	*/
    double input_vol,
    /*	Options	*/
    /*	Size of bumps	*/
    double Alpha_shift, double delay_shift, double Rho_shift,
    /*	Factor that multiplies results	*/
    double Alpha_fac, double delay_fac, double Rho_fac,
    /*	Answers	*/
    double *Alpha_sens, double *delay_sens, double *Rho_sens) {
  double frozen_vol;
  double prem, Alpha_prem, delay_prem, Rho_prem;
  Err err;

  frozen_vol = input_vol;

  /*	Calculate premium before shift	*/
  if (err = op_bvmsabr_pricing(forward, strike, maturity, disc, call_put, Alpha,
                               delay, Rho, frozen_vol, &prem)) {
    return err;
  }

  /*	Calculate premium after shift of Alpha	*/
  if (err = op_bvmsabr_pricing(forward, strike, maturity, disc, call_put,
                               Alpha + Alpha_shift, delay, Rho, frozen_vol,
                               &Alpha_prem)) {
    return err;
  }

  /*	Calculate premium after shift of delay	*/
  if (err = op_bvmsabr_pricing(forward, strike, maturity, disc, call_put, Alpha,
                               delay + delay_shift, Rho, frozen_vol,
                               &delay_prem)) {
    return err;
  }

  /*	Calculate premium after shift of Rho	*/
  if (err = op_bvmsabr_pricing(forward, strike, maturity, disc, call_put, Alpha,
                               delay, Rho + Rho_shift, frozen_vol, &Rho_prem)) {
    return err;
  }

  /*	Calculate sensitivites	*/
  *Alpha_sens = Alpha_fac * (Alpha_prem - prem) / Alpha_shift;
  *delay_sens = delay_fac * (delay_prem - prem) / delay_shift;
  *Rho_sens = Rho_fac * (Rho_prem - prem) / Rho_shift;

  return NULL;
}

/*	Vega/Volga	*/
/*	All parameters checks are supposed to have been done before call
 */
Err op_bvmsabr_vega_volga(
    /*	Black-Scholes parameters	*/
    double forward, double strike, double maturity, double disc,
    SrtCallPutType call_put,
    /*	Model parameters	*/
    double Alpha, double Delay, double Rho,
    /*	Vol input	*/
    double input_vol,
    /*	Options	*/
    /*	Size of bumps	*/
    double vega_shift, double volga_shift,
    /*	Multiplicative or additive	*/
    /*	0: Additive  , 1: Multiplicative	*/
    int vega_mult, int volga_mult,
    /*	Calculate vega or volga	*/
    /*	0: Vega  , 1: Volga	*/
    int vega_volga,
    /*	Factor that multiplies result	*/
    double fac,
    /*	Answer	*/
    double *sens) {
  double bump_vol;
  double prem, prem2, vega1, vega2;
  double act_vega_shift, act_volga_shift;
  Err err;

  bump_vol = input_vol;

  /*	Calculate premium before shift	*/
  if (err = op_bvmsabr_pricing(forward, strike, maturity, disc, call_put, Alpha,
                               Delay, Rho, bump_vol, &prem)) {
    return err;
  }

  /*	Calculate actual shifts	*/

  if (vega_mult == 1) {
    act_vega_shift = vega_shift * bump_vol;
  } else {
    act_vega_shift = vega_shift;
  }

  if (volga_mult == 1) {
    act_volga_shift = volga_shift * bump_vol;
  } else {
    act_volga_shift = volga_shift;
  }

  /*	Calculate premium after shift	*/
  if (err = op_bvmsabr_pricing(forward, strike, maturity, disc, call_put, Alpha,
                               Delay, Rho, bump_vol + act_vega_shift, &prem2)) {
    return err;
  }

  /*	Calculate vega before shift	*/
  vega1 = (prem2 - prem) / vega_shift;
  *sens = fac * vega1;

  /*	If volga is to be calculated  , calculate delta after shift	*/
  if (vega_volga == 1) {
    if (err = op_bvmsabr_vega_volga(
            forward, strike, maturity, disc, call_put, Alpha, Delay, Rho,
            bump_vol + act_volga_shift, vega_shift, volga_shift, vega_mult,
            volga_mult, 0, 1.0, &vega2)) {
      return err;
    }

    /*	Calculate volga	*/
    *sens = fac * (vega2 - vega1) / volga_shift;
  }

  return NULL;
}

Err BVMSabrPriceGradient(double strike, double paramSABR[], double *price,
                         double *gradient, int nbr_paramSABR) {
  Err err = NULL;
  double Alpha, Delay, Rho, forward, maturity, ATMvol, Betavol;
  int out_range = 0, i;

  Alpha = paramSABR[1];
  Delay = paramSABR[2];
  Rho = paramSABR[3];
  forward = paramSABR[4];
  maturity = paramSABR[5];
  ATMvol = paramSABR[6];

  Betavol = op_sabrgen_calib(forward, strike, maturity, ATMvol, Alpha, Delay, 0,
                             0, Rho, vol_bvm);

  if ((Alpha < 0) || (Delay < 0) || (Rho < -1) || (Rho > 1))
    out_range = 1;

  if (out_range == 1) {
    *price = -1000000;
    for (i = 1; i <= 6; i++) {
      gradient[i] = 1000000;
    }
  } else {
    err = op_bvmsabr_pricing(forward, strike, maturity, 1, SRT_CALL, Alpha,
                             Delay, Rho, ATMvol, price);

    err = op_bvmsabr_model_sens(forward, strike, maturity, 1, SRT_CALL, Alpha,
                                Delay, Rho, ATMvol, 1e-4, 1e-4, 1e-4, 1, 1, 1,
                                &(gradient[1]), &(gradient[2]), &(gradient[3]));

    gradient[4] = 0; /* No need to calculate it */
    gradient[5] = 0; /* No need to calculate it */

    //		err=op_sabr_vega_volga(
    err = op_bvmsabr_vega_volga(forward, strike, maturity, 1, SRT_CALL, Alpha,
                                Delay, Rho, ATMvol, 1e-5, 1e-4, 1, 1, 0, 1,
                                &(gradient[6]));
  }

  return err;
}

Err init_BVMSABR_parameters(double forward, double *strikes,
                            double *market_vols, int nbr_strikes, double ATMVol,
                            double *Alpha, int freeze_Alpha, double *Delay,
                            int freeze_Delay, double *Rho, int freeze_Rho) {
  Err err = NULL;

  return err;
}

Err opBVMsabrcalib(double forward, double maturity, int nbr_strikes,
                   double *strikes,     /*ACHTUNG: The Vector starts at 1*/
                   double *market_vols, /*ACHTUNG: The Vector starts at 1*/
                   double *ATMVol, double *Alpha,
                   int freeze_Alpha, /* if 0  , Alpha is not calibrated */
                   double *Delay,
                   int freeze_Delay, /* If 0  , Delay is not calibrated */
                   double *Rho,
                   int freeze_Rho, /* if 0  , Rho is not calibrated */
                   double *fitting_error) {
  Err err = NULL;
  double *weights_of_strikes = NULL,
         *paramSABR = NULL, /* Parameters of SabrPriceGradient */
             *market_prices = NULL;
  long *use_paramSABR =
           NULL, /* 0 for frozen parameter  , 1 for parameter to optimize */
      nbr_paramSABR;
  int nbr_iter = 200; /* Maximum number of iterations in Levenberg-Marquardt */
  int freeze_ATMVol;
  int i;

  /* Memory allocation */
  nbr_paramSABR = 6; /* Alpha  ,Beta  ,Rho  ,Forward  ,Maturity and ATMVol */
  weights_of_strikes = dvector(1, nbr_strikes);
  paramSABR = dvector(1, nbr_paramSABR);
  market_prices = dvector(1, nbr_strikes);
  use_paramSABR = lngvector(1, nbr_paramSABR);

  /* Computes the ATM vol by interpolation */
  err = compute_ATMVol(forward, strikes, market_vols, nbr_strikes,
                       &freeze_ATMVol, ATMVol);

  /* Computes a first guess for the BVMSABR parameters */
  err = init_BVMSABR_parameters(forward, strikes, market_vols, nbr_strikes,
                                *ATMVol, Alpha, freeze_Alpha, Delay,
                                freeze_Delay, Rho, freeze_Rho);

  /* Assignment of weights to strikes for fitting */
  err = compute_weights(strikes, nbr_strikes, forward, *ATMVol, maturity,
                        weights_of_strikes);

  /* Initializes SABR parameters */
  paramSABR[1] = *Alpha;
  paramSABR[2] = *Delay;
  paramSABR[3] = *Rho;
  paramSABR[4] = forward;
  paramSABR[5] = maturity;
  paramSABR[6] = *ATMVol;

  /* Chooses the parameters to be calibrated */
  use_paramSABR[1] = freeze_Alpha;
  use_paramSABR[2] = freeze_Delay;
  use_paramSABR[3] = freeze_Rho;
  use_paramSABR[4] = 0; /* Not calibrating the forward */
  use_paramSABR[5] = 0; /* Not calibrating the maturity */
  use_paramSABR[6] = freeze_ATMVol;

  /* Compute the market prices */
  for (i = 1; i <= nbr_strikes; i++) {
    market_prices[i] = srt_f_optblksch(forward, strikes[i], market_vols[i],
                                       maturity, 1, SRT_CALL, PREMIUM);
  }

  /* Call Levenberg-Marcquardt */

  err = levenberg_marquardt_select(
      strikes,                /* From [1] to [nbr_strikes] */
      market_prices,          /* From [1] to [nbr_strikes] */
      weights_of_strikes,     /* From [1] to [nbr_strikes] */
      nbr_strikes, paramSABR, /* From [1] to [nparam] */
      use_paramSABR,          /* From [1] to [nparam] */
      nbr_paramSABR, nbr_iter, BVMSabrPriceGradient, fitting_error);

  /* Fills the calibrated parameters */
  *Alpha = paramSABR[1];
  *Delay = paramSABR[2];
  *Rho = paramSABR[3];
  *ATMVol = paramSABR[6];

  /* Free memory */

  free_dvector(weights_of_strikes, 1, nbr_strikes);
  free_dvector(market_prices, 1, nbr_strikes);
  free_dvector(paramSABR, 1, nbr_paramSABR);
  free_lngvector(use_paramSABR, 1, nbr_paramSABR);

  return err;
}