/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        Description	:	implementation of the Heston model : Pricing

                                (C) 2002 BNP Paribas.. All rights reserved.

        Author		:	 Stefano Galluccio

        Created		:	14.11.2002

        History		:

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/

#include "opfnctns.h"
#include <math.h"
#include <opHeston.h"
#include <utallhdr.h"

/* 1) If the argument "greek" is set to "PREMIUM"  ,
   this function returns a pointer *result that contains the values
   of the shifted-Heston european options
   in correspondence to a string of increasing strikes to be input in *Strikes.

   2) If the argument "greek" is set to "DENSITY"  ,
   this function returns a pointer *result that contains the values of the
   probability density associated to the input parameters in correspondence to a
   string of increasing strikes to be input in *Strikes.

   3)If the argument "greek" is set to "CUMDENSITY"  ,
   this function returns a pointer *result that contains the values of the
   probability cumulative density associated to the input parameters in
   correspondence to a string of increasing strikes to be input in *Strikes.

  */

Err ConvertSABRParamInHestonParam(double Forward, double Maturity,
                                  double SABRBetaVol, double SABRAlpha,
                                  double SABRMeanReversion, double SABRBeta,
                                  double SABRRho, double *HESTONVol,
                                  double *HESTONAlpha,
                                  double *HESTONMeanReversion,
                                  double *HESTONBeta, double *HESTONRho) {
  Err err = NULL;
  double shift;

  shift = Forward * (1. - SABRBeta) / (1e-20 + SABRBeta);
  *HESTONVol = pow(Forward, SABRBeta) * SABRBetaVol / (shift + Forward);
  *HESTONAlpha = SABRAlpha * 2.0 * (*HESTONVol);
  //	*HESTONAlpha = SABRAlpha * 2.0 * (SABRBetaVol);
  *HESTONMeanReversion = SABRMeanReversion * 2.0;
  *HESTONBeta = SABRBeta;
  *HESTONRho = SABRRho;

  return err;
}

Err ConvertHESTONParamInSABRParam(double Forward, double Maturity,
                                  double HESTONVol, double HESTONAlpha,
                                  double HESTONMeanReversion, double HESTONBeta,
                                  double HESTONRho, double *SABRBetaVol,
                                  double *SABRAlpha, double *SABRMeanReversion,
                                  double *SABRBeta, double *SABRRho) {
  Err err = NULL;
  double shift;

  shift = Forward * (1. - HESTONBeta) / (1e-20 + HESTONBeta);
  *SABRBetaVol = HESTONVol * (shift + Forward) / pow(Forward, HESTONBeta);
  *SABRAlpha = HESTONAlpha / (2.0 * HESTONVol);
  *SABRMeanReversion = HESTONMeanReversion * 0.5;
  *SABRBeta = HESTONBeta;
  *SABRRho = HESTONRho;

  return err;
}

Err HestonPrice(double Forward, double *Strike, int nStrikes, double Maturity,
                double Sigma, double Alpha, double SigmaIfty, double b,
                double Beta, double Rho, double Disc, double UpperBound,
                SrtCallPutType call_put, SrtGreekType greek,
                SRT_Boolean isVolInfFix, int IntegrType, int nSteps,
                double *result) {

  double a, CosT1, SinT1, CosT2, SinT2;
  double *x, *w, *x1, *w1, G11, *G12, G21, *G22, *GKK, n1, n2, n3, RePsi1,
      ImPsi1, RePsi2, ImPsi2;
  double DRePsi1, DRePsi2, DImPsi1, DImPsi2, DDRePsi1, DDRePsi2, DDImPsi1,
      DDImPsi2, *premium;
  double esp, esp2, ImAppr1, ImAppr2, dStdDev, dApproxBSVol,
      nStdDev = 10.0; //  ,*Strike_used;
  int i, j, k, nStp = 5, iter, iThresh, nOptSteps1, nOptSteps2, nSteps_opt;
  double Int, Int1, NodeLeft, NodeRight, NodeLeft1, NodeRight1,
      dScaling = 1., dGear, Gamma, UpperBound_opt;
  Err err = NULL;
  double HSigma, HAlpha, Hb, HBeta, HSigmaIfty, HRho;

  err =
      ConvertSABRParamInHestonParam(Forward, Maturity, Sigma, Alpha, b, Beta,
                                    Rho, &HSigma, &HAlpha, &Hb, &HBeta, &HRho);

  HSigmaIfty = HSigma;
  isVolInfFix = 1;

  //	IntegrType=3;

  G12 = dvector(1, nStrikes);
  G22 = dvector(1, nStrikes);
  GKK = dvector(1, nStrikes);
  premium = dvector(1, nStrikes);

  // Transforms the input Beta parameter into a Lognormal shift

  //	Gamma = HESTON_CONST*(1.-Beta)/Beta;
  Gamma = Forward * (1. - HBeta) / HBeta;

  // here I define the max n. of Std. Dev.

  dApproxBSVol = (Forward + Gamma) / Forward * HSigma;
  dStdDev = Forward * sqrt(exp(dApproxBSVol * dApproxBSVol * Maturity) - 1.);

  // if Infinity Vol is fixed  , then a is set to the initial vol
  if (isVolInfFix)
    HSigmaIfty = HSigma;
  a = HSigmaIfty * HSigmaIfty * Hb;

  esp = exp(-Hb * Maturity);
  esp2 = exp(-Maturity * (Hb - HAlpha * HRho));

  /* Here I redefinee Strikes and Forward (to improve the convergence) by a
   * scaling factor */

  //	Forward+=Gamma;
  Forward = (Forward + Gamma) * dScaling;

  for (j = 1; j <= nStrikes; j++) {
    Strike[j] = (Strike[j] + Gamma) * dScaling;

    //		Strike[j] += Gamma;
  }

  switch (greek) {

  case PREMIUM:

    switch (IntegrType) {

    case 1:

      x = dvector(1, nStp);
      w = dvector(1, nStp);
      x1 = dvector(1, nStp);
      w1 = dvector(1, nStp);

      for (j = 1; j <= nStrikes; j++) {
        G12[j] = 0.0;
        G22[j] = 0.0;
      }

      ImAppr2 =
          -(log(Forward) +
            (esp2 - 1) * HSigma * HSigma / 2. / (Hb - HAlpha * HRho) +
            Maturity * a / 2. / (Hb - HAlpha * HRho) +
            a * (esp2 - 1.) / 2 / (Hb - HAlpha * HRho) / (Hb - HAlpha * HRho));
      ImAppr1 = -(log(Forward) + HSigma * HSigma / Hb * (1. - esp) -
                  Maturity * a / 2. / Hb + a / 2. / Hb / Hb * (1. - esp));

      for (j = 1; j <= nStrikes; j++) {

        iter = 0;
        NodeLeft = 0.0;
        NodeLeft1 = 0.0;
        do {

          /* find the next node on the axis */

          NodeRight = HestonFindRightNode(ImAppr1, Strike[j], NodeLeft,
                                          isVolInfFix, Maturity, a, Hb, HAlpha,
                                          HRho, Forward, HSigma, 1.0, greek);
          GaussLeg(NodeLeft, NodeRight, x, w, nStp);

          NodeRight1 = HestonFindRightNode(ImAppr2, Strike[j], NodeLeft,
                                           isVolInfFix, Maturity, a, Hb, HAlpha,
                                           HRho, Forward, HSigma, 0.0, greek);
          GaussLeg(NodeLeft1, NodeRight1, x1, w1, nStp);

          Int = 0.0;
          Int1 = 0.0;
          for (k = 1; k <= nStp; k++) {

            err =
                PsiFunction(isVolInfFix, Maturity, 1., -x[k], a, Hb, HAlpha,
                            HRho, log(Forward), HSigma * HSigma, greek, &RePsi1,
                            &ImPsi1, &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);

            G12[j] +=
                w[k] / x[k] * exp(RePsi1) * sin(ImPsi1 + x[k] * log(Strike[j]));

            err =
                PsiFunction(isVolInfFix, Maturity, 0., -x1[k], a, Hb, HAlpha,
                            HRho, log(Forward), HSigma * HSigma, greek, &RePsi2,
                            &ImPsi2, &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

            G22[j] += w1[k] / x1[k] * exp(RePsi2) *
                      sin(ImPsi2 + x1[k] * log(Strike[j]));
          }

          NodeLeft = NodeRight;
          NodeLeft1 = NodeRight1;

          iter++;

        } while ((NodeLeft <= UpperBound) && (iter <= 1000));
      }

      /* computes the remaining two terms that need no integration */

      err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      G11 = exp(RePsi1);

      err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
      G21 = exp(RePsi2);

      /* finally puts all together and evaluates the option */

      for (j = 1; j <= nStrikes; j++) {

        n1 = G11 / 2. - 1. / SRT_PI * G12[j];
        n2 = G21 / 2. - 1. / SRT_PI * G22[j];

        switch (call_put) {

        case SRT_CALL:
          result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc);

          if (result[j] < (Forward - Strike[j]) * Disc)
            result[j] = (Forward - Strike[j]) * Disc + 1.e-8;

          break;

        default:
          result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc -
                                       (Forward - Strike[j]) * Disc);

          if (result[j] < (-Forward + Strike[j]) * Disc)
            result[j] = (-Forward + Strike[j]) * Disc + 1.e-8;
          break;
        }
      }

      free_dvector(x, 1, nStp);
      free_dvector(w, 1, nStp);
      free_dvector(x1, 1, nStp);
      free_dvector(w1, 1, nStp);

      break;

    case 0:

      x = dvector(1, nSteps);
      w = dvector(1, nSteps);
      GaussLeg(0., UpperBound, x, w, nSteps);

      for (j = 1; j <= nStrikes; j++) {
        G12[j] = 0.0;
        G22[j] = 0.0;
      }

      /*
                              // computes the approximated Re(Psi)/v

                              PsiFunction (isVolInfFix  ,Maturity  ,1.
         ,-UpperBound  ,a  ,b  ,Alpha  ,Rho  ,log(Forward)  ,Sigma*Sigma  ,greek
         ,&ReUpBound  ,&ImUpBound  ,&DRePsi1  ,&DImPsi1  ,&DDRePsi1 ,&DDImPsi1);
                              PsiFunction (isVolInfFix  ,Maturity  ,1.  ,0.0  ,a
         ,b  ,Alpha  ,Rho  ,log(Forward)  ,Sigma*Sigma  ,greek  ,&ReLowBound
         ,&ImLowBound  ,&DRePsi1  ,&DImPsi1  ,&DDRePsi1  ,&DDImPsi1);

                              ReAppr1=(ReUpBound-ReLowBound)/(UpperBound-x[2]);
                              ReAppr1c=ReLowBound;

                              PsiFunction (isVolInfFix  ,Maturity  ,0.
         ,-UpperBound  ,a  ,b  ,Alpha  ,Rho  ,log(Forward)  ,Sigma*Sigma  ,greek
         ,&ReUpBound  ,&ImUpBound  ,&DRePsi1  ,&DImPsi1  ,&DDRePsi1 ,&DDImPsi1);
                              PsiFunction (isVolInfFix  ,Maturity  ,0.  ,0.0  ,a
         ,b  ,Alpha  ,Rho  ,log(Forward)  ,Sigma*Sigma  ,greek  ,&ReLowBound
         ,&ImLowBound  ,&DRePsi1  ,&DImPsi1  ,&DDRePsi1  ,&DDImPsi1);

                              ReAppr2=(ReUpBound-ReLowBound)/(UpperBound-x[2]);
                              ReAppr2c=ReLowBound;

                              // computes the approximated Im(Psi)/v

                              ImAppr2=-(log(Forward)+(esp2-1)*Sigma*Sigma/2./(b-Alpha*Rho)+Maturity*a/2./(b-Alpha*Rho)
                                                      +a*(esp2-1.)/2/(b-Alpha*Rho)/(b-Alpha*Rho));
                              ImAppr1=-(log(Forward)+Sigma*Sigma/b*(1.-esp)-Maturity*a/2./b+a/2./b/b*(1.-esp));

      */
      for (i = 1; i <= nSteps; i++) {

        err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                          log(Forward), HSigma * HSigma, greek, &RePsi1,
                          &ImPsi1, &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
        err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                          log(Forward), HSigma * HSigma, greek, &RePsi2,
                          &ImPsi2, &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

        // runs over the strikes

        for (j = 1; j <= nStrikes; j++) {

          G12[j] +=
              w[i] / x[i] * exp(RePsi1) * sin(ImPsi1 + x[i] * log(Strike[j]));
          G22[j] +=
              w[i] / x[i] * exp(RePsi2) * sin(ImPsi2 + x[i] * log(Strike[j]));

          /*
                                          G12[j]+=w[i]/x[i]*exp(RePsi1)*sin(ImPsi1+x[i]*log(Strike[j]))
             - w[i]/x[i]*exp(ReAppr1c+ReAppr1*x[i]) *
             sin(x[i]*log(Strike[j])+x[i]*ImAppr1);

                                          G22[j]+=w[i]/x[i]*exp(RePsi2)*sin(ImPsi2+x[i]*log(Strike[j]))-
                                              w[i]/x[i]*exp(ReAppr2c+ReAppr2*x[i])
             * sin(x[i]*log(Strike[j])+x[i]*ImAppr2);
          */
        }
      }

      err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      G11 = exp(RePsi1);

      err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
      G21 = exp(RePsi2);

      for (j = 1; j <= nStrikes; j++) {

        /*
                                        ContrVar =
           exp(ReAppr1c)*atan(-(log(Strike[j])+ImAppr1)/ReAppr1);
                                        n1=G11/2.-1./SRT_PI*(G12[j]+ContrVar);

                                        ContrVar =
           exp(ReAppr2c)*atan(-(log(Strike[j])+ImAppr2)/ReAppr2);
                                        n2=G21/2.-1./SRT_PI*(G22[j]+ContrVar);
        */
        n1 = G11 / 2. - 1. / SRT_PI * G12[j];
        n2 = G21 / 2. - 1. / SRT_PI * G22[j];

        switch (call_put) {

        case SRT_CALL:
          result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc) / dScaling;

          if (result[j] < (Forward - Strike[j]) * Disc)
            result[j] = (Forward - Strike[j]) * Disc / dScaling + 1.e-8;

          break;

        default:
          result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc / dScaling -
                                       (Forward - Strike[j]) * Disc) /
                      dScaling;

          if (result[j] < (-Forward + Strike[j]) * Disc)
            result[j] = (-Forward + Strike[j]) * Disc + 1.e-8;
          break;
        }
      }

      free_dvector(x, 1, nSteps);
      free_dvector(w, 1, nSteps);

      break;

    case 3:

      //		Strike_used=dvector(1  ,nStrikes);

      for (j = 1; j <= nStrikes; j++) {
        G12[j] = 0.0;
        G22[j] = 0.0;
      }

      ImAppr2 =
          -(log(Forward) +
            (esp2 - 1) * HSigma * HSigma / 2. / (Hb - HAlpha * HRho) +
            Maturity * a / 2. / (Hb - HAlpha * HRho) +
            a * (esp2 - 1.) / 2 / (Hb - HAlpha * HRho) / (Hb - HAlpha * HRho));
      ImAppr1 = -(log(Forward) + HSigma * HSigma / Hb * (1. - esp) -
                  Maturity * a / 2. / Hb + a / 2. / Hb / Hb * (1. - esp));

      UpperBound_opt = 60. + 20 * pow(DMAX(0.0, 5. - Maturity), 1.4);

      for (j = 1; j <= nStrikes; j++) {

        // computes the approximated n. of oscillations in a unit interval
        // and then the required n. integration steps. For very short maturities
        // , the precision is forced to increse  , as there is little time value

        // Here it branches following the rule: if Strike is larger than nStdDev
        // from the forward or smaller than -nStdDev  , then it simply evaluates
        // the price at +- nStdDev

        dGear = 2 + 3. / (Maturity * Maturity * Maturity);
        nOptSteps1 = IMAX(IMIN((int)(dGear * UpperBound *
                                     (fabs(ImAppr1 + log(Strike[j]))) / SRT_PI),
                               nSteps),
                          22);
        nOptSteps2 = IMAX(IMIN((int)(dGear * UpperBound *
                                     (fabs(ImAppr2 + log(Strike[j]))) / SRT_PI),
                               nSteps),
                          22);

        nSteps_opt = IMAX(nOptSteps1, nOptSteps2);
        x = dvector(1, nSteps_opt);
        w = dvector(1, nSteps_opt);

        GaussLeg(0., UpperBound_opt, x, w, nSteps_opt);

        //		Strike_used[j] = DMIN(Strike[j]  ,Forward +
        //nStdDev*dStdDev); 		Strike_used[j] = DMAX(Strike_used[j]  , DMAX(0.0015
        //,Forward - nStdDev*dStdDev));

        Strike[j] = DMIN(Strike[j], Forward + nStdDev * dStdDev);
        Strike[j] = DMAX(Strike[j], DMAX(0.00001, Forward - nStdDev * dStdDev));

        for (i = 1; i <= nSteps_opt; i++) {

          err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha,
                            HRho, log(Forward), HSigma * HSigma, greek, &RePsi1,
                            &ImPsi1, &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
          err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha,
                            HRho, log(Forward), HSigma * HSigma, greek, &RePsi2,
                            &ImPsi2, &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

          //				G12[j]+=w[i]/x[i]*exp(RePsi1)*sin(ImPsi1+x[i]*log(Strike_used[j]));
          //				G22[j]+=w[i]/x[i]*exp(RePsi2)*sin(ImPsi2+x[i]*log(Strike_used[j]));

          G12[j] +=
              w[i] / x[i] * exp(RePsi1) * sin(ImPsi1 + x[i] * log(Strike[j]));
          G22[j] +=
              w[i] / x[i] * exp(RePsi2) * sin(ImPsi2 + x[i] * log(Strike[j]));
        }

        err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                          log(Forward), HSigma * HSigma, greek, &RePsi1,
                          &ImPsi1, &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
        G11 = exp(RePsi1);

        err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                          log(Forward), HSigma * HSigma, greek, &RePsi2,
                          &ImPsi2, &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
        G21 = exp(RePsi2);

        n1 = G11 / 2. - 1. / SRT_PI * G12[j];
        n2 = G21 / 2. - 1. / SRT_PI * G22[j];

        switch (call_put) {

        case SRT_CALL:
          result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc) / dScaling;
          //					result[j]=DMAX(1.e-13
          //,(n1-Strike_used[j]*n2)*Disc)/dScaling;

          if (result[j] < (Forward - Strike[j]) * Disc)
            result[j] = (Forward - Strike[j]) * Disc / dScaling + 1.e-8;
          //					if (result[j]< (Forward-Strike_used[j])*Disc)
          //result[j]=(Forward-Strike_used[j])*Disc/dScaling+1.e-8;

          break;

        default:
          result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc / dScaling -
                                       (Forward - Strike[j]) * Disc) /
                      dScaling;
          //					result[j]= DMAX(1.e-13
          //,(n1-Strike_used[j]*n2)*Disc/dScaling -
          //(Forward-Strike_used[j])*Disc)/dScaling;

          if (result[j] < (-Forward + Strike[j]) * Disc)
            result[j] = (-Forward + Strike[j]) * Disc + 1.e-8;
          //					if (result[j]< (-Forward+Strike_used[j])*Disc)
          //result[j]=(-Forward+Strike_used[j])*Disc+1.e-8;
          break;
        }

        free_dvector(x, 1, nSteps_opt);
        free_dvector(w, 1, nSteps_opt);
      }

      //			free_dvector(Strike_used  ,1  ,nStrikes);

      break;

    default:

      /* runs over strikes to set the breakpoint between the two methods */

      iThresh = nStrikes;
      for (j = 0; j < nStrikes; j++) {
        if (Strike[j + 1] > Forward * 0.1 * dScaling) {
          iThresh = j;
          break;
        }
      }

      /* now I integrate with the optimised algorithm for strikes < breakpoint
         and with the parallel algorithm for strikes  > breakpoint */

      x = dvector(1, nStp);
      w = dvector(1, nStp);
      x1 = dvector(1, nStp);
      w1 = dvector(1, nStp);

      for (j = 1; j <= iThresh; j++) {
        G12[j] = 0.0;
        G22[j] = 0.0;
      }

      ImAppr2 =
          -(log(Forward) +
            (esp2 - 1) * HSigma * HSigma / 2. / (Hb - HAlpha * HRho) +
            Maturity * a / 2. / (Hb - HAlpha * HRho) +
            a * (esp2 - 1.) / 2 / (Hb - HAlpha * HRho) / (Hb - HAlpha * HRho));
      ImAppr1 = -(log(Forward) + HSigma * HSigma / Hb * (1. - esp) -
                  Maturity * a / 2. / Hb + a / 2. / Hb / Hb * (1. - esp));

      for (j = 1; j <= iThresh; j++) {

        iter = 0;
        NodeLeft = 0.0;
        NodeLeft1 = 0.0;
        do {

          /* find the next node on the axis */

          NodeRight = HestonFindRightNode(ImAppr1, Strike[j], NodeLeft,
                                          isVolInfFix, Maturity, a, Hb, HAlpha,
                                          HRho, Forward, HSigma, 1.0, greek);
          GaussLeg(NodeLeft, NodeRight, x, w, nStp);

          NodeRight1 = HestonFindRightNode(ImAppr2, Strike[j], NodeLeft,
                                           isVolInfFix, Maturity, a, Hb, HAlpha,
                                           HRho, Forward, HSigma, 0.0, greek);
          GaussLeg(NodeLeft1, NodeRight1, x1, w1, nStp);

          Int = 0.0;
          Int1 = 0.0;
          for (k = 1; k <= nStp; k++) {

            err =
                PsiFunction(isVolInfFix, Maturity, 1., -x[k], a, Hb, HAlpha,
                            HRho, log(Forward), HSigma * HSigma, greek, &RePsi1,
                            &ImPsi1, &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);

            G12[j] +=
                w[k] / x[k] * exp(RePsi1) * sin(ImPsi1 + x[k] * log(Strike[j]));

            err =
                PsiFunction(isVolInfFix, Maturity, 0., -x1[k], a, Hb, HAlpha,
                            HRho, log(Forward), HSigma * HSigma, greek, &RePsi2,
                            &ImPsi2, &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

            G22[j] += w1[k] / x1[k] * exp(RePsi2) *
                      sin(ImPsi2 + x1[k] * log(Strike[j]));
          }

          NodeLeft = NodeRight;
          NodeLeft1 = NodeRight1;

          iter++;

        } while ((NodeLeft <= UpperBound) && (iter <= 1000));
      }

      free_dvector(x, 1, nStp);
      free_dvector(w, 1, nStp);
      free_dvector(x1, 1, nStp);
      free_dvector(w1, 1, nStp);

      /* second part  , now I use the parallel integrator for strikes > 30 bps
       */

      x = dvector(1, nSteps);
      w = dvector(1, nSteps);
      GaussLeg(0., UpperBound, x, w, nSteps);

      for (j = iThresh + 1; j <= nStrikes; j++) {
        G12[j] = 0.0;
        G22[j] = 0.0;
      }

      for (i = 1; i <= nSteps; i++) {

        err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                          log(Forward), HSigma * HSigma, greek, &RePsi1,
                          &ImPsi1, &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
        err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                          log(Forward), HSigma * HSigma, greek, &RePsi2,
                          &ImPsi2, &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

        // runs over the strikes

        for (j = iThresh + 1; j <= nStrikes; j++) {

          G12[j] +=
              w[i] / x[i] * exp(RePsi1) * sin(ImPsi1 + x[i] * log(Strike[j]));
          G22[j] +=
              w[i] / x[i] * exp(RePsi2) * sin(ImPsi2 + x[i] * log(Strike[j]));
        }
      }

      err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      G11 = exp(RePsi1);

      err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
      G21 = exp(RePsi2);

      /* Now sums up all together */

      for (j = 1; j <= nStrikes; j++) {

        n1 = G11 / 2. - 1. / SRT_PI * G12[j];
        n2 = G21 / 2. - 1. / SRT_PI * G22[j];

        switch (call_put) {

        case SRT_CALL:
          result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc) / dScaling;

          if (result[j] < (Forward - Strike[j]) * Disc)
            result[j] = (Forward - Strike[j]) * Disc / dScaling + 1.e-8;

          break;

        default:
          result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc / dScaling -
                                       (Forward - Strike[j]) * Disc) /
                      dScaling;

          if (result[j] < (-Forward + Strike[j]) * Disc)
            result[j] = (-Forward + Strike[j]) * Disc + 1.e-8;
          break;
        }
      }

      free_dvector(x, 1, nSteps);
      free_dvector(w, 1, nSteps);

      break;
    }

    break;

  case DELTA:

    x = dvector(1, nSteps);
    w = dvector(1, nSteps);
    GaussLeg(0., UpperBound, x, w, nSteps);

    for (j = 1; j <= nStrikes; j++) {
      G12[j] = 0.0;
      G22[j] = 0.0;
    }

    for (i = 1; i <= nSteps; i++) {

      err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

      // runs over the strikes

      for (j = 1; j <= nStrikes; j++) {

        SinT1 = sin(ImPsi1 + x[i] * log(Strike[j]));
        CosT1 = cos(ImPsi1 + x[i] * log(Strike[j]));

        SinT2 = sin(ImPsi2 + x[i] * log(Strike[j]));
        CosT2 = cos(ImPsi2 + x[i] * log(Strike[j]));

        G12[j] +=
            w[i] / x[i] * exp(RePsi1) * (SinT1 * DRePsi1 + CosT1 * DImPsi1);
        G22[j] +=
            w[i] / x[i] * exp(RePsi2) * (SinT2 * DRePsi2 + CosT2 * DImPsi2);
      }
    }

    err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                      &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
    G11 = exp(RePsi1) * DRePsi1;

    err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                      &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
    G21 = exp(RePsi2) * DRePsi2;

    for (j = 1; j <= nStrikes; j++) {

      n1 = G11 / 2. - 1. / SRT_PI * G12[j];
      n2 = G21 / 2. - 1. / SRT_PI * G22[j];

      switch (call_put) {

      case SRT_CALL:
        result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc);

        break;

      default:
        result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc) - Disc;

        break;
      }
    }

    free_dvector(x, 1, nSteps);
    free_dvector(w, 1, nSteps);

    break;

  case GAMMA:

    x = dvector(1, nSteps);
    w = dvector(1, nSteps);
    GaussLeg(0., UpperBound, x, w, nSteps);

    for (j = 1; j <= nStrikes; j++) {
      G12[j] = 0.0;
      G22[j] = 0.0;
    }

    for (i = 1; i <= nSteps; i++) {

      err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

      // runs over the strikes

      for (j = 1; j <= nStrikes; j++) {

        SinT1 = sin(ImPsi1 + x[i] * log(Strike[j]));
        CosT1 = cos(ImPsi1 + x[i] * log(Strike[j]));

        SinT2 = sin(ImPsi2 + x[i] * log(Strike[j]));
        CosT2 = cos(ImPsi2 + x[i] * log(Strike[j]));

        G12[j] += w[i] / x[i] * exp(RePsi1) *
                  (SinT1 * DRePsi1 * DRePsi1 + CosT1 * DImPsi1 * DRePsi1 +
                   DDRePsi1 * SinT1 + DRePsi1 * DImPsi1 * CosT1 -
                   DImPsi1 * DImPsi1 * SinT1 + DDImPsi1 * CosT1);
        G22[j] += w[i] / x[i] * exp(RePsi2) *
                  (SinT2 * DRePsi2 * DRePsi2 + CosT2 * DImPsi2 * DRePsi2 +
                   DDRePsi2 * SinT2 + DRePsi2 * DImPsi2 * CosT2 -
                   DImPsi2 * DImPsi2 * SinT2 + DDImPsi2 * CosT2);
      }
    }

    err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                      &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
    G11 = exp(RePsi1) * (DRePsi1 * DRePsi1 + DDRePsi1);

    err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                      &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
    G21 = exp(RePsi2) * (DRePsi2 * DRePsi2 + DDRePsi2);

    for (j = 1; j <= nStrikes; j++) {

      n1 = G11 / 2. - 1. / SRT_PI * G12[j];
      n2 = G21 / 2. - 1. / SRT_PI * G22[j];

      result[j] = DMAX(1.e-13, (n1 - Strike[j] * n2) * Disc);
    }

    free_dvector(x, 1, nSteps);
    free_dvector(w, 1, nSteps);
    break;

  case VEGA:

    x = dvector(1, nSteps);
    w = dvector(1, nSteps);
    GaussLeg(0., UpperBound, x, w, nSteps);

    for (j = 1; j <= nStrikes; j++) {
      G12[j] = 0.0;
      G22[j] = 0.0;
    }

    for (i = 1; i <= nSteps; i++) {

      err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

      // runs over the strikes

      for (j = 1; j <= nStrikes; j++) {

        SinT1 = sin(ImPsi1 + x[i] * log(Strike[j]));
        CosT1 = cos(ImPsi1 + x[i] * log(Strike[j]));

        SinT2 = sin(ImPsi2 + x[i] * log(Strike[j]));
        CosT2 = cos(ImPsi2 + x[i] * log(Strike[j]));

        G12[j] +=
            w[i] / x[i] * exp(RePsi1) * (SinT1 * DRePsi1 + CosT1 * DImPsi1);
        G22[j] +=
            w[i] / x[i] * exp(RePsi2) * (SinT2 * DRePsi2 + CosT2 * DImPsi2);
      }
    }

    err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                      &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
    G11 = exp(RePsi1) * DRePsi1;

    err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                      &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
    G21 = exp(RePsi2) * DRePsi2;

    for (j = 1; j <= nStrikes; j++) {

      n1 = G11 / 2. - 1. / SRT_PI * G12[j];
      n2 = G21 / 2. - 1. / SRT_PI * G22[j];

      result[j] = (n1 - Strike[j] * n2) * Disc;
    }

    free_dvector(x, 1, nSteps);
    free_dvector(w, 1, nSteps);
    break;

  case VOLGA:

    x = dvector(1, nSteps);
    w = dvector(1, nSteps);
    GaussLeg(0., UpperBound, x, w, nSteps);

    for (j = 1; j <= nStrikes; j++) {
      G12[j] = 0.0;
      G22[j] = 0.0;
    }

    for (i = 1; i <= nSteps; i++) {

      err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

      // runs over the strikes

      for (j = 1; j <= nStrikes; j++) {

        SinT1 = sin(ImPsi1 + x[i] * log(Strike[j]));
        CosT1 = cos(ImPsi1 + x[i] * log(Strike[j]));

        SinT2 = sin(ImPsi2 + x[i] * log(Strike[j]));
        CosT2 = cos(ImPsi2 + x[i] * log(Strike[j]));

        G12[j] += w[i] / x[i] * exp(RePsi1) *
                  (SinT1 * DRePsi1 * DRePsi1 + CosT1 * DImPsi1 * DRePsi1 +
                   DDRePsi1 * SinT1 + DRePsi1 * DImPsi1 * CosT1 -
                   DImPsi1 * DImPsi1 * SinT1 + DDImPsi1 * CosT1);
        G22[j] += w[i] / x[i] * exp(RePsi2) *
                  (SinT2 * DRePsi2 * DRePsi2 + CosT2 * DImPsi2 * DRePsi2 +
                   DDRePsi2 * SinT2 + DRePsi2 * DImPsi2 * CosT2 -
                   DImPsi2 * DImPsi2 * SinT2 + DDImPsi2 * CosT2);
      }
    }

    err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                      &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
    G11 = exp(RePsi1) * (DRePsi1 * DRePsi1 + DDRePsi1);

    err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                      &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
    G21 = exp(RePsi2) * (DRePsi2 * DRePsi2 + DDRePsi2);

    for (j = 1; j <= nStrikes; j++) {

      n1 = G11 / 2. - 1. / SRT_PI * G12[j];
      n2 = G21 / 2. - 1. / SRT_PI * G22[j];

      result[j] = (n1 - Strike[j] * n2) * Disc;
    }

    free_dvector(x, 1, nSteps);
    free_dvector(w, 1, nSteps);
    break;

  case VANNA:

    x = dvector(1, nSteps);
    w = dvector(1, nSteps);
    GaussLeg(0., UpperBound, x, w, nSteps);

    for (j = 1; j <= nStrikes; j++) {
      G12[j] = 0.0;
      G22[j] = 0.0;
    }

    for (i = 1; i <= nSteps; i++) {

      err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

      // runs over the strikes

      for (j = 1; j <= nStrikes; j++) {

        SinT1 = sin(ImPsi1 + x[i] * log(Strike[j]));
        CosT1 = cos(ImPsi1 + x[i] * log(Strike[j]));

        SinT2 = sin(ImPsi2 + x[i] * log(Strike[j]));
        CosT2 = cos(ImPsi2 + x[i] * log(Strike[j]));

        G12[j] += w[i] / x[i] * exp(RePsi1) *
                  (SinT1 * DDRePsi1 * DRePsi1 + CosT1 * DRePsi1 * DDImPsi1 +
                   DImPsi1 * DDRePsi1 * CosT1 - DDImPsi1 * DImPsi1 * SinT1);
        G22[j] += w[i] / x[i] * exp(RePsi2) *
                  (SinT2 * DDRePsi2 * DRePsi2 + CosT2 * DRePsi2 * DDImPsi2 +
                   DImPsi2 * DDRePsi2 * CosT2 - DDImPsi2 * DImPsi2 * SinT2);
      }
    }

    err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                      &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
    G11 = exp(RePsi1) * (DRePsi1 * DDRePsi1);

    err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                      &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
    G21 = exp(RePsi2) * (DRePsi2 * DDRePsi2);

    for (j = 1; j <= nStrikes; j++) {

      n1 = G11 / 2. - 1. / SRT_PI * G12[j];
      n2 = G21 / 2. - 1. / SRT_PI * G22[j];

      result[j] = (n1 - Strike[j] * n2) * Disc;
    }

    free_dvector(x, 1, nSteps);
    free_dvector(w, 1, nSteps);
    break;

  case THETA:

    x = dvector(1, nSteps);
    w = dvector(1, nSteps);
    GaussLeg(0., UpperBound, x, w, nSteps);

    /* I first compute the derivative of the term containing the DF */

    greek = PREMIUM;
    for (j = 1; j <= nStrikes; j++) {
      G12[j] = 0.0;
      G22[j] = 0.0;
    }

    for (i = 1; i <= nSteps; i++) {

      err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

      // runs over the strikes

      for (j = 1; j <= nStrikes; j++) {

        G12[j] +=
            w[i] / x[i] * exp(RePsi1) * sin(ImPsi1 + x[i] * log(Strike[j]));
        G22[j] +=
            w[i] / x[i] * exp(RePsi2) * sin(ImPsi2 + x[i] * log(Strike[j]));
      }
    }

    err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                      &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
    G11 = exp(RePsi1);

    err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                      &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
    G21 = exp(RePsi2);

    for (j = 1; j <= nStrikes; j++) {

      n1 = G11 / 2. - 1. / SRT_PI * G12[j];
      n2 = G21 / 2. - 1. / SRT_PI * G22[j];

      switch (call_put) {

      case SRT_CALL:
        result[j] = -(n1 - Strike[j] * n2) * Disc * log(Disc) / Maturity;

        break;

      default:
        result[j] = -(n1 - Strike[j] * n2) * Disc * log(Disc) / Maturity +
                    (Forward - Strike[j]) * Disc * log(Disc) / Maturity;

        break;
      }
    }

    /* I now compute the derivative of the term equivalent to [Fwd*N(d1)-KN(d2)]
     */

    greek = THETA;
    for (j = 1; j <= nStrikes; j++) {
      G12[j] = 0.0;
      G22[j] = 0.0;
    }

    for (i = 1; i <= nSteps; i++) {

      err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

      // runs over the strikes

      for (j = 1; j <= nStrikes; j++) {

        SinT1 = sin(ImPsi1 + x[i] * log(Strike[j]));
        CosT1 = cos(ImPsi1 + x[i] * log(Strike[j]));

        SinT2 = sin(ImPsi2 + x[i] * log(Strike[j]));
        CosT2 = cos(ImPsi2 + x[i] * log(Strike[j]));

        G12[j] +=
            w[i] / x[i] * exp(RePsi1) * (SinT1 * DRePsi1 + CosT1 * DImPsi1);
        G22[j] +=
            w[i] / x[i] * exp(RePsi2) * (SinT2 * DRePsi2 + CosT2 * DImPsi2);
      }
    }

    err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                      &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
    G11 = exp(RePsi1) * DRePsi1;

    err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                      log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                      &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
    G21 = exp(RePsi2) * DRePsi2;

    for (j = 1; j <= nStrikes; j++) {

      n1 = G11 / 2. - 1. / SRT_PI * G12[j];
      n2 = G21 / 2. - 1. / SRT_PI * G22[j];

      result[j] += -(n1 - Strike[j] * n2) * Disc;
    }

    free_dvector(x, 1, nSteps);
    free_dvector(w, 1, nSteps);
    break;

  case DENSITY:

    x = dvector(1, nSteps);
    w = dvector(1, nSteps);
    GaussLeg(0., UpperBound, x, w, nSteps);

    /* now I split the integration in two parts */

    for (i = 1; i <= nSteps; i++) {

      err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

      // runs over the strikes

      for (j = 1; j <= nStrikes; j++) {

        G12[j] += w[i] * exp(RePsi1) *
                  (-cos(ImPsi1 + x[i] * log(Strike[j])) -
                   x[i] * sin(ImPsi1 + x[i] * log(Strike[j]))) /
                  Strike[j] / Strike[j];
        G22[j] += w[i] * exp(RePsi2) *
                  (-cos(ImPsi2 + x[i] * log(Strike[j])) -
                   x[i] * sin(ImPsi2 + x[i] * log(Strike[j]))) /
                  Strike[j] / Strike[j];
        GKK[j] += w[i] * exp(RePsi2) * cos(ImPsi2 + x[i] * log(Strike[j])) /
                  Strike[j];
      }
    }

    for (j = nStrikes; j >= 1; j--) {

      n1 = -1. / SRT_PI * G12[j];
      n2 = -1. / SRT_PI * G22[j];
      n3 = -1. / SRT_PI * GKK[j];

      result[j] = DMAX(0.0, n1 - 2. * n3 - Strike[j] * n2);
    }

    free_dvector(x, 1, nSteps);
    free_dvector(w, 1, nSteps);
    break;

  case CUMDENSITY:

    switch (IntegrType) {

    case 1:

      x = dvector(1, nStp);
      w = dvector(1, nStp);
      x1 = dvector(1, nStp);
      w1 = dvector(1, nStp);

      for (j = 1; j <= nStrikes; j++) {

        G12[j] = 0.0;
        G22[j] = 0.0;
        GKK[j] = 0.0;
      }

      ImAppr2 =
          -(log(Forward) +
            (esp2 - 1) * HSigma * HSigma / 2. / (Hb - HAlpha * HRho) +
            Maturity * a / 2. / (Hb - HAlpha * HRho) +
            a * (esp2 - 1.) / 2 / (Hb - HAlpha * HRho) / (Hb - HAlpha * HRho));
      ImAppr1 = -(log(Forward) + HSigma * HSigma / Hb * (1. - esp) -
                  Maturity * a / 2. / Hb + a / 2. / Hb / Hb * (1. - esp));

      for (j = 1; j <= nStrikes; j++) {

        iter = 0;
        NodeLeft = 0.0;
        NodeLeft1 = 0.0;
        do {

          /* find the next node on the axis */

          NodeRight = HestonFindRightNode(ImAppr1, Strike[j], NodeLeft,
                                          isVolInfFix, Maturity, a, Hb, HAlpha,
                                          HRho, Forward, HSigma, 1.0, greek);
          GaussLeg(NodeLeft, NodeRight, x, w, nStp);

          NodeRight1 = HestonFindRightNode(ImAppr2, Strike[j], NodeLeft,
                                           isVolInfFix, Maturity, a, Hb, HAlpha,
                                           HRho, Forward, HSigma, 0.0, greek);
          GaussLeg(NodeLeft1, NodeRight1, x1, w1, nStp);

          Int = 0.0;
          Int1 = 0.0;
          for (k = 1; k <= nStp; k++) {

            err =
                PsiFunction(isVolInfFix, Maturity, 1., -x[k], a, Hb, HAlpha,
                            HRho, log(Forward), HSigma * HSigma, greek, &RePsi1,
                            &ImPsi1, &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);

            G12[j] += w[k] * exp(RePsi1) * cos(ImPsi1 + x[k] * log(Strike[j])) /
                      Strike[j];

            err =
                PsiFunction(isVolInfFix, Maturity, 0., -x1[k], a, Hb, HAlpha,
                            HRho, log(Forward), HSigma * HSigma, greek, &RePsi2,
                            &ImPsi2, &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

            G22[j] += w1[k] * exp(RePsi2) *
                      cos(ImPsi2 + x1[k] * log(Strike[j])) / Strike[j];
            GKK[j] += w1[k] / x1[k] * exp(RePsi2) *
                      sin(ImPsi2 + x1[k] * log(Strike[j]));
          }

          NodeLeft = NodeRight;
          NodeLeft1 = NodeRight1;

          iter++;

        } while ((NodeLeft <= UpperBound) && (iter <= 1000));
      }

      /* computes the remaining two terms that need no integration */

      err = PsiFunction(isVolInfFix, Maturity, 1., 0., a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi1, &ImPsi1,
                        &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
      G11 = exp(RePsi1);

      err = PsiFunction(isVolInfFix, Maturity, 0., 0., a, Hb, HAlpha, HRho,
                        log(Forward), HSigma * HSigma, greek, &RePsi2, &ImPsi2,
                        &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);
      G21 = exp(RePsi2);

      /* finally puts all together and evaluates the option */

      for (j = 1; j <= nStrikes; j++) {

        n1 = -1. / SRT_PI * G12[j];
        n2 = -1. / SRT_PI * G22[j];
        n3 = -1. / SRT_PI * GKK[j];

        result[j] = DMAX(DMIN(1.0, 0.5 + (n1 - n3 - Strike[j] * n2)), 0.0);
      }

      free_dvector(x, 1, nStp);
      free_dvector(w, 1, nStp);
      free_dvector(x1, 1, nStp);
      free_dvector(w1, 1, nStp);

      break;

    default:

      x = dvector(1, nSteps);
      w = dvector(1, nSteps);
      GaussLeg(0., UpperBound, x, w, nSteps);

      for (j = 1; j <= nStrikes; j++) {

        G12[j] = 0.0;
        G22[j] = 0.0;
        GKK[j] = 0.0;
      }

      for (i = 1; i <= nSteps; i++) {

        err = PsiFunction(isVolInfFix, Maturity, 1., -x[i], a, Hb, HAlpha, HRho,
                          log(Forward), HSigma * HSigma, greek, &RePsi1,
                          &ImPsi1, &DRePsi1, &DImPsi1, &DDRePsi1, &DDImPsi1);
        err = PsiFunction(isVolInfFix, Maturity, 0., -x[i], a, Hb, HAlpha, HRho,
                          log(Forward), HSigma * HSigma, greek, &RePsi2,
                          &ImPsi2, &DRePsi2, &DImPsi2, &DDRePsi2, &DDImPsi2);

        // runs over the strikes

        for (j = 1; j <= nStrikes; j++) {

          G12[j] += w[i] * exp(RePsi1) * cos(ImPsi1 + x[i] * log(Strike[j])) /
                    Strike[j];
          G22[j] += w[i] * exp(RePsi2) * cos(ImPsi2 + x[i] * log(Strike[j])) /
                    Strike[j];
          GKK[j] +=
              w[i] / x[i] * exp(RePsi2) * sin(ImPsi2 + x[i] * log(Strike[j]));
        }
      }

      for (j = 1; j <= nStrikes; j++) {

        n1 = -1. / SRT_PI * G12[j];
        n2 = -1. / SRT_PI * G22[j];
        n3 = -1. / SRT_PI * GKK[j];

        result[j] = DMAX(DMIN(1.0, 0.5 + (n1 - n3 - Strike[j] * n2)), 0.0);
      }

      free_dvector(x, 1, nSteps);
      free_dvector(w, 1, nSteps);
    }

    break;
  }

  Forward -= Gamma;

  for (j = 1; j <= nStrikes; j++) {
    Strike[j] -= Gamma;
  }

  free_dvector(G12, 1, nStrikes);
  free_dvector(G22, 1, nStrikes);
  free_dvector(GKK, 1, nStrikes);
  free_dvector(premium, 1, nStrikes);

  return NULL;
}
