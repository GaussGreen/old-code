/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        Description	:	implementation of the Heston model: Volatility
conversions

                                (C) 2002 BNP Paribas.. All rights reserved.

        Author		:	 Stefano Galluccio

        Created		:	14.11.2002

        History		:

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/

#include "opfnctns.h"
#include <math.h"
#include <opHeston.h"
#include <utallhdr.h"

/*-----------------------------------------------------------------------------------------------
        This function converts among different types of volatilities by matching
the premiums of the corresponding calls/puts. It is the equivalent of the add-in
SABRVol A set of increasing strikes *Strikes must be provided as input  , and
the equivalent volatilities are returned as *NewVols
/*-----------------------------------------------------------------------------------------------*/

#define MAXIT 100
Err HestonVol(double Forward, double *Strikes, int nStrikes, double Maturity,
              double Vol, double VolInfty, double Alpha, double Beta,
              double Rho, double MeanRev, SrtVolConversion TypeConversion,
              double UpperBound, int nSteps, int IntegerType,
              SRT_Boolean isVolInfFix, double *NewVols) {

  Err err = NULL;
  int j; //  ,IntegrType=3; /* sets up the integration type to the optimised one
         //  for a single  */
         //	double df  ,dx  ,dxold  ,f  ,fh  ,fl  ,x1=0.0001  ,x2=2.0;
         //	double temp  ,xh  ,xl  ,rts  ,xacc=1.e-7;
  double BSPremium, /*HestonPremium  ,*/ ImplVol;
  double *result = NULL, *result_low = NULL, *result_high = NULL,
         *Strike = NULL;
  double BetaVol;

  result = dvector(1, nStrikes);
  result_low = dvector(1, nStrikes);
  result_high = dvector(1, nStrikes);
  Strike = dvector(1, 1);

  switch (TypeConversion) {

  case LOG_TO_LOG:
    err = HestonCalibrateSigmaInit(Forward, Maturity, Vol, Alpha, VolInfty,
                                   MeanRev, Beta, Rho, UpperBound, nSteps,
                                   IntegerType, isVolInfFix, &BetaVol);
    err = HestonPrice(Forward, Strikes, nStrikes, Maturity, BetaVol, Alpha,
                      VolInfty, MeanRev, Beta, Rho, 1.0, UpperBound, SRT_CALL,
                      PREMIUM, isVolInfFix, IntegerType, nSteps, result);

    for (j = 1; j <= nStrikes; j++) {

      err = srt_f_optimpvol(result[j], Forward, Strikes[j], Maturity, 1.0,
                            SRT_CALL, SRT_LOGNORMAL, &ImplVol);

      NewVols[j] = ImplVol;
    }

    break;

  case LOG_TO_HESTON:

    err = HestonCalibrateSigmaInit(Forward, Maturity, Vol, Alpha, VolInfty,
                                   MeanRev, Beta, Rho, UpperBound, nSteps,
                                   IntegerType, isVolInfFix, &Vol);

    for (j = 1; j <= nStrikes; j++) {
      NewVols[j] = Vol;
    }

    /* computes premium at the left boundary */
    /*			err =  HestonPrice(
                                                                    Forward  ,
                                                                    Strikes  ,
                                                                    1  ,
                                                                    Maturity  ,
                                                                    x1  ,
                                                                    Alpha  ,
                                                                    VolInfty  ,
                                                                    MeanRev  ,
                                                                    Beta  ,
                                                                    Rho  ,
                                                                    1.0  ,
                                                                    UpperBound ,
                                                                    SRT_CALL  ,
                                                                    PREMIUM  ,
                                                                    isVolInfFix
    , IntegerType  , nSteps  , result_low
                                                              );

    //			 computes premium at the right boundary
                            err =  HestonPrice(
                                                                    Forward  ,
                                                                    Strikes  ,
                                                                    1  ,
                                                                    Maturity  ,
                                                                    x2  ,
                                                                    Alpha  ,
                                                                    VolInfty  ,
                                                                    MeanRev  ,
                                                                    Beta  ,
                                                                    Rho  ,
                                                                    1.0  ,
                                                                    UpperBound ,
                                                                    SRT_CALL  ,
                                                                    PREMIUM  ,
                                                                    isVolInfFix
    , IntegerType  , nSteps  , result_high
                                                              );


                            for (i=1;i<=nStrikes;i++) {

                                    BSPremium = srt_f_optblksch(
                                                                                            Forward  ,
                                                                                            Strikes[i]  ,
                                                                                            Vol  ,
                                                                                            Maturity  ,
                                                                                            1.0  ,
                                                                                            SRT_CALL  ,
                                                                                            PREMIUM
                                                                            );


                                    fl=result_low[i]-BSPremium;
                                    fh=result_high[i]-BSPremium;

                                    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 &&
    fh < 0.0));  // Throw exception: not bracketed

                                    if (fl == 0.0) {
                                            NewVols[i]=x1;
                                            goto EndLoop;
                                            // return err;
                                    }
                                    if (fh == 0.0) {
                                            NewVols[i]=x2;
                                            goto EndLoop;
                                            // return err;
                                    }
                                    if (fl < 0.0) {
                                            xl=x1;
                                            xh=x2;
                                    } else {
                                            xh=x1;
                                            xl=x2;
                                    }

                                    rts=0.5*(x1+x2);
                                    dxold=fabs(x2-x1);
                                    dx=dxold;

                                    Strike[1]=Strikes[i];

                                    err =  HestonPrice(
                                                                            Forward
    , Strike  , 1  , Maturity  , rts  , Alpha  , VolInfty  , MeanRev  , Beta  ,
                                                                            Rho
    , 1.0  , UpperBound  , SRT_CALL  , PREMIUM  , isVolInfFix  , IntegerType  ,
                                                                            nSteps
    , result
                                                                      );

                                    f=result[1]-BSPremium;

                                    err =  HestonPrice(
                                                                            Forward
    , Strike  , 1  , Maturity  , rts  , Alpha  , VolInfty  , MeanRev  , Beta  ,
                                                                            Rho
    , 1.0  , UpperBound  , SRT_CALL  , VEGA  , isVolInfFix  , IntegerType  ,
                                                                            nSteps
    , result
                                                                      );
                                    df=result[1];

                                            for (j=1;j<=MAXIT;j++) {
                                                    if
    ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
                                                            || (fabs(2.0*f) >
    fabs(dxold*df))) { dxold=dx; dx=0.5*(xh-xl); rts=xl+dx;

                                                            if (xl == rts) {
                                                                    NewVols[i]=rts;
                                                                    goto
    EndLoop;
                                                                    // return
    err;
                                                            }

                                                    } else {
                                                            dxold=dx;
                                                            dx=f/df;
                                                            temp=rts;
                                                            rts -= dx;

                                                            if (temp == rts) {
                                                                    NewVols[i]=rts;
                                                                    goto
    EndLoop;
                                                                    // return
    err;
                                                            }

                                                    }
                                                    if (fabs(dx) < xacc) {
                                                            NewVols[i]=rts;
                                                            goto EndLoop;
                                                            // return err;
                                                    }

                                                    err =  HestonPrice(
                                                                                            Forward  ,
                                                                                            Strike  ,
                                                                                            1  ,
                                                                                            Maturity  ,
                                                                                            rts  ,
                                                                                            Alpha  ,
                                                                                            VolInfty  ,
                                                                                            MeanRev  ,
                                                                                            Beta  ,
                                                                                            Rho  ,
                                                                                            1.0  ,
                                                                                            UpperBound  ,
                                                                                            SRT_CALL  ,
                                                                                            PREMIUM  ,
                                                                                            isVolInfFix  ,
                                                                                            IntegerType  ,
                                                                                            nSteps  ,
                                                                                            result
                                                                                      );

                                                    f=result[1]-BSPremium;

                                                    if (f < 0.0)
                                                            xl=rts;
                                                    else
                                                            xh=rts;
                                            }

    EndLoop: ;
                    }

                      return err;
    */

    break;

  case HESTON_TO_LOG:

    err = HestonPrice(Forward, Strikes, nStrikes, Maturity, Vol, Alpha,
                      VolInfty, MeanRev, Beta, Rho, 1.0, UpperBound, SRT_CALL,
                      PREMIUM, isVolInfFix, IntegerType, nSteps, result);

    for (j = 1; j <= nStrikes; j++) {

      err = srt_f_optimpvol(result[j], Forward, Strikes[j], Maturity, 1.0,
                            SRT_CALL, SRT_LOGNORMAL, &ImplVol);

      NewVols[j] = ImplVol;
    }

    break;

  case HESTON_TO_HESTON:

    break;

  case HESTON_TO_NORMAL:

    err = HestonPrice(Forward, Strikes, nStrikes, Maturity, Vol, Alpha,
                      VolInfty, MeanRev, Beta, Rho, 1.0, UpperBound, SRT_CALL,
                      PREMIUM, isVolInfFix, IntegerType, nSteps, result);

    for (j = 1; j <= nStrikes; j++) {

      err = srt_f_optimpvol(result[j], Forward, Strikes[j], Maturity, 1.0,
                            SRT_CALL, SRT_NORMAL, &ImplVol);

      NewVols[j] = ImplVol;
    }

    break;

  case NORMAL_TO_HESTON:

    /*			err =  HestonPrice(
                                                                    Forward  ,
                                                                    Strikes  ,
                                                                    1  ,
                                                                    Maturity  ,
                                                                    x1  ,
                                                                    Alpha  ,
                                                                    VolInfty  ,
                                                                    MeanRev  ,
                                                                    Beta  ,
                                                                    Rho  ,
                                                                    1.0  ,
                                                                    UpperBound ,
                                                                    SRT_CALL  ,
                                                                    PREMIUM  ,
                                                                    isVolInfFix
    , IntegerType  , nSteps  , result_low
                                                              );

                            err =  HestonPrice(
                                                                    Forward  ,
                                                                    Strikes  ,
                                                                    1  ,
                                                                    Maturity  ,
                                                                    x2  ,
                                                                    Alpha  ,
                                                                    VolInfty  ,
                                                                    MeanRev  ,
                                                                    Beta  ,
                                                                    Rho  ,
                                                                    1.0  ,
                                                                    UpperBound ,
                                                                    SRT_CALL  ,
                                                                    PREMIUM  ,
                                                                    isVolInfFix
    , IntegerType  , nSteps  , result_high
                                                              );


                            for (i=1;i<=nStrikes;i++) {

                                    BSPremium = srt_f_optblknrm(
                                                                                            Forward  ,
                                                                                            Strikes[i]  ,
                                                                                            Vol  ,
                                                                                            Maturity  ,
                                                                                            1.0  ,
                                                                                            SRT_CALL  ,
                                                                                            PREMIUM
                                                                            );

                                    fl=result_low[i]-BSPremium;
                                    fh=result_high[i]-BSPremium;

                                    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 &&
    fh < 0.0))   ;

                                    if (fl == 0.0) {
                                            NewVols[i]=x1;
                                            goto EndLoop2;
                                            // return err;
                                    }
                                    if (fh == 0.0) {
                                            NewVols[i]=x2;
                                            goto EndLoop2;
                                            // return err;
                                    }
                                    if (fl < 0.0) {
                                            xl=x1;
                                            xh=x2;
                                    } else {
                                            xh=x1;
                                            xl=x2;
                                    }

                                    rts=0.5*(x1+x2);
                                    dxold=fabs(x2-x1);
                                    dx=dxold;

                                    Strike[1]=Strikes[i];

                                    err =  HestonPrice(
                                                                            Forward
    , Strike  , 1  , Maturity  , rts  , Alpha  , VolInfty  , MeanRev  , Beta  ,
                                                                            Rho
    , 1.0  , UpperBound  , SRT_CALL  , PREMIUM  , isVolInfFix  , IntegerType  ,
                                                                            nSteps
    , result
                                                                      );

                                    f=result[1]-BSPremium;

                                    err =  HestonPrice(
                                                                            Forward
    , Strike  , 1  , Maturity  , rts  , Alpha  , VolInfty  , MeanRev  , Beta  ,
                                                                            Rho
    , 1.0  , UpperBound  , SRT_CALL  , VEGA  , isVolInfFix  , IntegerType  ,
                                                                            nSteps
    , result
                                                                      );
                                    df=result[1];

                                            for (j=1;j<=MAXIT;j++) {
                                                    if
    ((((rts-xh)*df-f)*((rts-xl)*df-f) >= 0.0)
                                                            || (fabs(2.0*f) >
    fabs(dxold*df))) { dxold=dx; dx=0.5*(xh-xl); rts=xl+dx;

                                                            if (xl == rts) {
                                                                    NewVols[i]=rts;
                                                                    goto
    EndLoop;
                                                                    // return
    err;
                                                            }

                                                    } else {
                                                            dxold=dx;
                                                            dx=f/df;
                                                            temp=rts;
                                                            rts -= dx;

                                                            if (temp == rts) {
                                                                    NewVols[i]=rts;
                                                                    goto
    EndLoop2;
                                                                    // return
    err;
                                                            }

                                                    }
                                                    if (fabs(dx) < xacc) {
                                                            NewVols[i]=rts;
                                                            goto EndLoop2;
                                                            // return err;
                                                    }

                                                    err =  HestonPrice(
                                                                                            Forward  ,
                                                                                            Strike  ,
                                                                                            1  ,
                                                                                            Maturity  ,
                                                                                            rts  ,
                                                                                            Alpha  ,
                                                                                            VolInfty  ,
                                                                                            MeanRev  ,
                                                                                            Beta  ,
                                                                                            Rho  ,
                                                                                            1.0  ,
                                                                                            UpperBound  ,
                                                                                            SRT_CALL  ,
                                                                                            PREMIUM  ,
                                                                                            isVolInfFix  ,
                                                                                            IntegerType  ,
                                                                                            nSteps  ,
                                                                                            result
                                                                                      );

                                                    f=result[1]-BSPremium;

                                                    if (f < 0.0)
                                                            xl=rts;
                                                    else
                                                            xh=rts;
                                            }

    EndLoop2: ;
                    }

                      return err;
    */

    BSPremium = srt_f_optblknrm(Forward, Forward, Vol, Maturity, 1.0, SRT_CALL,
                                PREMIUM);

    err = srt_f_optimpvol(BSPremium, Forward, Forward, Maturity, 1.0, SRT_CALL,
                          SRT_LOGNORMAL, &ImplVol);

    err = HestonCalibrateSigmaInit(Forward, Maturity, ImplVol, Alpha, VolInfty,
                                   MeanRev, Beta, Rho, UpperBound, nSteps,
                                   IntegerType, isVolInfFix, &Vol);

    for (j = 1; j <= nStrikes; j++) {
      NewVols[j] = Vol;
    }

    break;

  case NORMAL_TO_LOG:

    for (j = 1; j <= nStrikes; j++) {

      BSPremium = srt_f_optblknrm(Forward, Strikes[j], Vol, Maturity, 1.0,
                                  SRT_CALL, PREMIUM);

      err = srt_f_optimpvol(BSPremium, Forward, Strikes[j], Maturity, 1.0,
                            SRT_CALL, SRT_LOGNORMAL, &ImplVol);

      NewVols[j] = ImplVol;
    }

    break;

  case LOG_TO_NORMAL:

    for (j = 1; j <= nStrikes; j++) {

      BSPremium = srt_f_optblksch(Forward, Strikes[j], Vol, Maturity, 1.0,
                                  SRT_CALL, PREMIUM);

      err = srt_f_optimpvol(BSPremium, Forward, Strikes[j], Maturity, 1.0,
                            SRT_CALL, SRT_NORMAL, &ImplVol);

      NewVols[j] = ImplVol;
    }

    break;
  }

  free_dvector(result, 1, nStrikes);
  free_dvector(result_low, 1, nStrikes);
  free_dvector(result_high, 1, nStrikes);
  free_dvector(Strike, 1, 1);

  return err;
}

Err srt_opthestonvol2(double Forward, double Strike, double Maturity,
                      double Sigma, double Alpha, double Beta, double Rho,
                      double Lambda, double *NumericalParams, int nParams,
                      SrtDiffusionType input, SrtDiffusionType output,
                      double *OutputVol) {
  Err err = NULL;
  double *StrikeVect = NULL;
  double *OutPutVolVect = NULL;

  StrikeVect = dvector(1, 1);
  OutPutVolVect = dvector(1, 1);
  StrikeVect[1] = Strike;

  err = srt_opthestonvol(Forward, StrikeVect, 1, Maturity, Sigma, Alpha, Beta,
                         Rho, Lambda, input, output, NumericalParams, nParams,
                         OutPutVolVect);

  *OutputVol = OutPutVolVect[1];

  if (StrikeVect)
    free_dvector(StrikeVect, 1, 1);
  if (OutPutVolVect)
    free_dvector(OutPutVolVect, 1, 1);

  return err;
}

Err srt_opthestonvol(double Forward, double *Strikes, int nStrikes,
                     double Maturity, double Sigma, double Alpha, double Beta,
                     double Rho, double Lambda, SrtDiffusionType input,
                     SrtDiffusionType output, double *NumericalParams,
                     int nParams, double *OutputVols) {

  Err err = NULL;
  int j;
  double ImplVol;
  double *result = NULL, *result_low = NULL, *result_high = NULL;
  double BetaVol, ATMLogVol;
  double price;
  double UpperBound;
  int nSteps;
  int IntegerType;

  if (nParams != 3) {
    err = "srt_opthestonvol needs 3-dim numerical params array";
    goto FREE_RETURN;
  }

  UpperBound = NumericalParams[0];
  nSteps = (int)(NumericalParams[1] + 1e-5);
  IntegerType = (int)(NumericalParams[2] + 1e-5);

  result = dvector(1, nStrikes);
  result_low = dvector(1, nStrikes);
  result_high = dvector(1, nStrikes);

  if (input == SRT_LOGNORMAL) {
    err = HestonCalibrateSigmaInit(Forward, Maturity, Sigma, Alpha, Sigma,
                                   Lambda, Beta, Rho, UpperBound, nSteps,
                                   IntegerType, 1, &BetaVol);

    err = HestonPrice(Forward, Strikes, nStrikes, Maturity, BetaVol, Alpha,
                      BetaVol, Lambda, Beta, Rho, 1.0, UpperBound, SRT_CALL,
                      PREMIUM, 1, IntegerType, nSteps, result);

    if (output == SRT_LOGNORMAL) {
      for (j = 1; j <= nStrikes; j++) {

        err = srt_f_optimpvol(result[j], Forward, Strikes[j], Maturity, 1.0,
                              SRT_CALL, SRT_LOGNORMAL, &ImplVol);
        OutputVols[j] = ImplVol;
      }
    } else if (output == SRT_BETAVOL) {
      for (j = 1; j <= nStrikes; j++) {
        OutputVols[j] = Sigma;
      }
    } else if (output == SRT_NORMAL) {
      for (j = 1; j <= nStrikes; j++) {

        err = srt_f_optimpvol(result[j], Forward, Strikes[j], Maturity, 1.0,
                              SRT_CALL, SRT_NORMAL, &ImplVol);
        OutputVols[j] = ImplVol;
      }
    } else {
      err = "OutPut Vol has to be LOG  , NORM or BETA";
      goto FREE_RETURN;
    }
  } else if (input == SRT_NORMAL) {
    price = srt_f_optblksch(Forward, Forward, Sigma, Maturity, 1.0, SRT_CALL,
                            PREMIUM);
    err = srt_f_optimpvol(price, Forward, Forward, Maturity, 1.0, SRT_CALL,
                          SRT_LOGNORMAL, &ATMLogVol);
    err = HestonCalibrateSigmaInit(Forward, Maturity, ATMLogVol, Alpha,
                                   ATMLogVol, Lambda, Beta, Rho, UpperBound,
                                   nSteps, IntegerType, 1, &BetaVol);
    err = HestonPrice(Forward, Strikes, nStrikes, Maturity, BetaVol, Alpha,
                      BetaVol, Lambda, Beta, Rho, 1.0, UpperBound, SRT_CALL,
                      PREMIUM, 1, IntegerType, nSteps, result);

    if (output == SRT_LOGNORMAL) {
      for (j = 1; j <= nStrikes; j++) {

        err = srt_f_optimpvol(result[j], Forward, Strikes[j], Maturity, 1.0,
                              SRT_CALL, SRT_LOGNORMAL, &ImplVol);
        OutputVols[j] = ImplVol;
      }
    } else if (output == SRT_BETAVOL) {
      for (j = 1; j <= nStrikes; j++) {
        OutputVols[j] = Sigma;
      }
    } else if (output == SRT_NORMAL) {
      for (j = 1; j <= nStrikes; j++) {

        err = srt_f_optimpvol(result[j], Forward, Strikes[j], Maturity, 1.0,
                              SRT_CALL, SRT_NORMAL, &ImplVol);
        OutputVols[j] = ImplVol;
      }
    } else {
      err = "OutPut Vol has to be LOG  , NORM or BETA";
      goto FREE_RETURN;
    }
  } else if (input == SRT_BETAVOL) {
    err = HestonPrice(Forward, Strikes, nStrikes, Maturity, Sigma, Alpha, Sigma,
                      Lambda, Beta, Rho, 1.0, UpperBound, SRT_CALL, PREMIUM, 1,
                      IntegerType, nSteps, result);

    if (output == SRT_LOGNORMAL) {
      for (j = 1; j <= nStrikes; j++) {

        err = srt_f_optimpvol(result[j], Forward, Strikes[j], Maturity, 1.0,
                              SRT_CALL, SRT_LOGNORMAL, &ImplVol);
        OutputVols[j] = ImplVol;
      }
    } else if (output == SRT_BETAVOL) {
      for (j = 1; j <= nStrikes; j++) {
        OutputVols[j] = Sigma;
      }
    } else if (output == SRT_NORMAL) {
      for (j = 1; j <= nStrikes; j++) {

        err = srt_f_optimpvol(result[j], Forward, Strikes[j], Maturity, 1.0,
                              SRT_CALL, SRT_NORMAL, &ImplVol);
        OutputVols[j] = ImplVol;
      }
    } else {
      err = "OutPut Vol has to be LOG  , NORM or BETA";
      goto FREE_RETURN;
    }
  } else {
    err = "InPut Vol has to be LOG  , NORM or BETA";
    goto FREE_RETURN;
  }

FREE_RETURN:

  free_dvector(result, 1, nStrikes);
  free_dvector(result_low, 1, nStrikes);
  free_dvector(result_high, 1, nStrikes);

  return err;
}

#undef MAXIT

/*-----------------------------------------------------------------------------------------------
  This function returns the equivalent BS volatility for a given set of Heston
parameters and a given strike. This overloads the function HestonVol for this
particular case
-----------------------------------------------------------------------------------------------*/

Err HestonATMVol(double Forward, double Maturity, double Sigma, double Alpha,
                 double SigmaIfty, double b, double Beta, double Rho,
                 double Disc, double UpperBound, int nSteps,
                 SRT_Boolean isVolInfFix, double *result) {
  Err err = NULL;
  int nStrikes = 1;
  double *Strike;
  double *Price;
  double vol;
  int IntegrType = 3;

  Strike = dvector(0, 1);
  Price = dvector(0, 1);
  Strike[1] = Forward;

  if (isVolInfFix)
    SigmaIfty = Sigma;

  err = HestonPrice(Forward, Strike, nStrikes, Maturity, Sigma, Alpha,
                    SigmaIfty, b, Beta, Rho, Disc, UpperBound, SRT_CALL,
                    PREMIUM, isVolInfFix, IntegrType, nSteps, Price);

  err = srt_f_optimpvol(Price[1], Forward, Forward, Maturity, 1, SRT_CALL,
                        SRT_LOGNORMAL, &vol);

  *result = vol;
  free_dvector(Strike, 0, 1);
  free_dvector(Price, 0, 1);

  return err;
}