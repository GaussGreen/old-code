/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        Description	:	implementation of the Heston model: Generic
mathematical auxiliary routines

                                (C) 2002 BNP Paribas.. All rights reserved.

        Author		:	 Stefano Galluccio

        Created		:	14.11.2002

        History		:

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/

#include "opfnctns.h"
#include <math.h"
#include <opHeston.h"
#include <utallhdr.h"

/*---------------------------------------------------------------------------------------
  This function returns the Real and Imaginary parts of the PsiFunction
(Arrow-Debreau) and their derivatives that form the integrand of the option
price.
---------------------------------------------------------------------------------------*/

Err PsiFunction(SRT_Boolean isVolInfFix, double T, double U, double V, double a,
                double b, double alpha, double rho, double x1, double x2,
                SrtGreekType greek, double *RePsi, double *ImPsi,
                double *DRePsi, double *DImPsi, double *DDRePsi,
                double *DDImPsi)

{

  double Phi1, Phi2, Phi3, Phi4, Psi1, Psi2, Psi3, Psi4, NumRe, NumIm, DenRe,
      DenIm, Chi, Den, DDen, PhiG, PsiG, PhiG2, PsiG2, Phi3G, Psi3G, SinT,
      Sin2T, CosT, Cos2T, U2 = U * U, U3 = U * U * U, V2 = V * V,
                          V3 = V * V * V, r = 0.0, NumRe2, NumIm2, Phi2G, Psi2G,
                          Esp, Esp2, DNumRe, DNumIm, DLogArg;
  double alpha2 = alpha * alpha, alpharho = alpha * rho,
         alpharho2 = alpharho * alpharho, Talpha2 = T / alpha2,
         ralpha2 = r * alpha2, aalpha2 = a / alpha2, Ualpharho = U * alpharho,
         Valpharho = V * alpharho, PowG, ATanG, VPhiG, VPsiG, Den2ReDen2Im,
         SinTEsp, CosTEsp, Sin2TEsp2, Cos2TEsp2;

  Err err = NULL;

  //  Defines the Re and Im part of the square of the function gamma

  PhiG2 = b * b + alpha2 * (U - U2 + V2) - 2 * b * Ualpharho +
          alpharho2 * (U2 - V2);

  PsiG2 = V * (alpha2 * (1 - 2 * U) - 2 * (b * alpharho - U * alpharho2));

  //  Defines the Re and Im part of the function gamma

  PowG = pow(PhiG2 * PhiG2 + PsiG2 * PsiG2, 0.25);
  ATanG = atan(PsiG2 / PhiG2) / 2.0;
  PhiG = PowG * cos(ATanG);
  PsiG = PowG * sin(ATanG);
  VPhiG = V * PhiG;
  VPsiG = V * PsiG;

  Phi2G = PhiG * PhiG;
  Psi2G = PsiG * PsiG;
  Phi3G = Phi2G * PhiG;
  Psi3G = Psi2G * PsiG;

  SinT = sin(T * PsiG);
  Sin2T = SinT * SinT;
  CosT = cos(T * PsiG);
  Cos2T = CosT * CosT;
  Esp = exp(-T * PhiG);
  Esp2 = Esp * Esp;
  SinTEsp = SinT * Esp;
  CosTEsp = CosT * Esp;
  Sin2TEsp2 = Sin2T * Esp2;
  Cos2TEsp2 = Cos2T * Esp2;

  //  Defines the L1 function
  Phi1 = x1 * U;
  Psi1 = x1 * V;

  //  Defines the L2 function
  NumRe = -x2 * ((1 - 2 * CosTEsp + Cos2TEsp2 + Sin2TEsp2) *
                     (b * (U - U2 + V2) - alpharho * (U2 - U3 + V2 - U * V2)) +
                 (1 - Cos2TEsp2 - Sin2TEsp2) *
                     (PhiG * (U - U2 + V2) + VPsiG * (1 - 2 * U)) -
                 2 * SinTEsp * (VPhiG * (1 - 2 * U) - PsiG * (U - U2 + V2)));

  NumIm = -x2 * ((1 - 2 * CosTEsp + Cos2TEsp2 + Sin2TEsp2) *
                     (b * V * (1 - 2 * U) + alpharho * (U2 * V + V3)) +
                 (1 - Cos2TEsp2 - Sin2TEsp2) *
                     (VPhiG * (1 - 2 * U) - PsiG * (U - U2 + V2)) +
                 2 * SinTEsp * (PhiG * (U - U2 + V2) + VPsiG * (1 - 2 * U)));

  DenRe = b - Ualpharho + PhiG -
          Esp * (CosT * (b - Ualpharho - PhiG) - SinT * (Valpharho + PsiG));

  DenIm = -Valpharho + PsiG +
          Esp * (CosT * (Valpharho + PsiG) + SinT * (b - Ualpharho - PhiG));

  Den2ReDen2Im = 1.0 / (DenRe * DenRe + DenIm * DenIm);
  Phi2 = NumRe * Den2ReDen2Im;
  Psi2 = NumIm * Den2ReDen2Im;

  //  Defines the L3 function

  Phi3 = Talpha2 * (ralpha2 * (U - 1) + a * (b - Ualpharho - PhiG));
  Psi3 = Talpha2 * (ralpha2 * V - a * (Valpharho + PsiG));

  //  Defines the L4 function

  Chi = Phi2G + Psi2G;
  Den = 2 * Chi;

  NumRe =
      PhiG * (b - Ualpharho) - Phi2G + 2 * Chi - Valpharho * PsiG - Psi2G -
      Esp *
          (CosT * (PhiG * (b - Ualpharho) - Phi2G - Valpharho * PsiG - Psi2G) -
           SinT * (Valpharho * PhiG + PsiG * (b - Ualpharho)));

  NumIm =
      -Valpharho * PhiG - PsiG * (b - Ualpharho) +
      Esp *
          (CosT * (Valpharho * PhiG + PsiG * (b - Ualpharho)) +
           SinT * (PhiG * (b - Ualpharho) - Phi2G - Valpharho * PsiG - Psi2G));

  Phi4 = -aalpha2 * log((NumRe * NumRe + NumIm * NumIm) / (Den * Den));
  Psi4 = -2 * aalpha2 * atan(NumIm / NumRe);

  // computes Re and Im part of the Psi function

  *RePsi = Phi1 + Phi2 + Phi3 + Phi4;
  *ImPsi = Psi1 + Psi2 + Psi3 + Psi4;

  // Now I compute the greeks

  switch (greek) {

  case PREMIUM:

    *DRePsi = 0.0;
    *DImPsi = 0.0;
    *DDRePsi = 0.0;
    *DDImPsi = 0.0;
    break;

  case DELTA:

    x1 = exp(x1); // exp(x1) is the fwd.
    *DRePsi = 1. / x1 * U;
    *DImPsi = 1. / x1 * V;
    *DDRePsi = 0.0;
    *DDImPsi = 0.0;

    break;

  case GAMMA:

    x1 = exp(x1); // exp(x1) is the fwd.
    *DRePsi = 1. / x1 * U;
    *DImPsi = 1. / x1 * V;
    *DDRePsi = -1. / x1 / x1 * U;
    *DDImPsi = -1. / x1 / x1 * V;

    break;

  case VEGA:

    //		a=2.*b*sqrt(x2);
    x2 = 2. * sqrt(x2);

    NumRe =
        -(b * U * x2) + b * U2 * x2 - b * V2 * x2 + U2 * x2 * alpha * rho -
        U3 * x2 * alpha * rho + V2 * x2 * alpha * rho -
        U * V2 * x2 * alpha * rho - U * x2 * PhiG + U2 * x2 * PhiG -
        V2 * x2 * PhiG - V * x2 * PsiG + 2 * U * V * x2 * PsiG +
        (2 * b * U * x2 * CosT) * Esp - (2 * b * U2 * x2 * CosT) * Esp +
        (2 * b * V2 * x2 * CosT) * Esp -
        (2 * U2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U3 * x2 * alpha * rho * CosT) * Esp -
        (2 * V2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U * V2 * x2 * alpha * rho * CosT) * Esp -
        (b * U * x2 * Cos2T) * Esp2 + (b * U2 * x2 * Cos2T) * Esp2 -
        (b * V2 * x2 * Cos2T) * Esp2 + (U2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U3 * x2 * alpha * rho * Cos2T) * Esp2 +
        (V2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Cos2T) * Esp2 +
        (U * x2 * PhiG * Cos2T) * Esp2 - (U2 * x2 * PhiG * Cos2T) * Esp2 +
        (V2 * x2 * PhiG * Cos2T) * Esp2 + (V * x2 * PsiG * Cos2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Cos2T) * Esp2 +
        (2 * V * x2 * PhiG * SinT) * Esp -
        (4 * U * V * x2 * PhiG * SinT) * Esp -
        (2 * U * x2 * PsiG * SinT) * Esp + (2 * U2 * x2 * PsiG * SinT) * Esp -
        (2 * V2 * x2 * PsiG * SinT) * Esp - (b * U * x2 * Sin2T) * Esp2 +
        (b * U2 * x2 * Sin2T) * Esp2 - (b * V2 * x2 * Sin2T) * Esp2 +
        (U2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U3 * x2 * alpha * rho * Sin2T) * Esp2 +
        (V2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Sin2T) * Esp2 +
        (U * x2 * PhiG * Sin2T) * Esp2 - (U2 * x2 * PhiG * Sin2T) * Esp2 +
        (V2 * x2 * PhiG * Sin2T) * Esp2 + (V * x2 * PsiG * Sin2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Sin2T) * Esp2;

    NumIm = -(b * V * x2) + 2 * b * U * V * x2 - U2 * V * x2 * alpha * rho -
            V3 * x2 * alpha * rho - V * x2 * PhiG + 2 * U * V * x2 * PhiG +
            U * x2 * PsiG - U2 * x2 * PsiG + V2 * x2 * PsiG +
            (2 * b * V * x2 * CosT) * Esp - (4 * b * U * V * x2 * CosT) * Esp +
            (2 * U2 * V * x2 * alpha * rho * CosT) * Esp +
            (2 * V3 * x2 * alpha * rho * CosT) * Esp -
            (b * V * x2 * Cos2T) * Esp2 + (2 * b * U * V * x2 * Cos2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Cos2T) * Esp2 -
            (V3 * x2 * alpha * rho * Cos2T) * Esp2 +
            (V * x2 * PhiG * Cos2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Cos2T) * Esp2 -
            (U * x2 * PsiG * Cos2T) * Esp2 + (U2 * x2 * PsiG * Cos2T) * Esp2 -
            (V2 * x2 * PsiG * Cos2T) * Esp2 - (2 * U * x2 * PhiG * SinT) * Esp +
            (2 * U2 * x2 * PhiG * SinT) * Esp -
            (2 * V2 * x2 * PhiG * SinT) * Esp -
            (2 * V * x2 * PsiG * SinT) * Esp +
            (4 * U * V * x2 * PsiG * SinT) * Esp - (b * V * x2 * Sin2T) * Esp2 +
            (2 * b * U * V * x2 * Sin2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Sin2T) * Esp2 -
            (V3 * x2 * alpha * rho * Sin2T) * Esp2 +
            (V * x2 * PhiG * Sin2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Sin2T) * Esp2 -
            (U * x2 * PsiG * Sin2T) * Esp2 + (U2 * x2 * PsiG * Sin2T) * Esp2 -
            (V2 * x2 * PsiG * Sin2T) * Esp2;

    DenRe = b - U * alpha * rho + PhiG - (b * CosT) * Esp +
            (U * alpha * rho * CosT) * Esp + (PhiG * CosT) * Esp +
            (V * alpha * rho * SinT) * Esp + (PsiG * SinT) * Esp;

    DenIm = -(V * alpha * rho) + PsiG + (V * alpha * rho * CosT) * Esp +
            (PsiG * CosT) * Esp + (b * SinT) * Esp -
            (U * alpha * rho * SinT) * Esp - (PhiG * SinT) * Esp;

    *DRePsi = NumRe / (DenRe * DenRe + DenIm * DenIm);
    *DImPsi = NumIm / (DenRe * DenRe + DenIm * DenIm);
    *DDRePsi = 0.0;
    *DDImPsi = 0.0;

    /*
    This part is needed if one is interested in the case when Sigma(0) =
    Sigma(infinity)
*/

    if (isVolInfFix) {

      a = b * x2;
      Phi3 = T / alpha / alpha * (a * b - a * U * alpha * rho - a * PhiG);
      Psi3 = T / alpha / alpha * (-a * V * alpha * rho - a * PsiG);
      Chi = PhiG * PhiG + PsiG * PsiG;
      Den = 2 * Chi;

      NumRe2 = b * PhiG - U * alpha * rho * PhiG - Phi2G + 2 * Chi -
               V * alpha * rho * PsiG - Psi2G - (b * PhiG * CosT) * Esp +
               (U * alpha * rho * PhiG * CosT) * Esp + (Phi2G * CosT) * Esp +
               (V * alpha * rho * PsiG * CosT) * Esp + (Psi2G * CosT) * Esp +
               (V * alpha * rho * PhiG * SinT) * Esp + (b * PsiG * SinT) * Esp -
               (U * alpha * rho * PsiG * SinT) * Esp;

      NumIm2 = -(V * alpha * rho * PhiG) - b * PsiG + U * alpha * rho * PsiG +
               (V * alpha * rho * PhiG * CosT) * Esp + (b * PsiG * CosT) * Esp -
               (U * alpha * rho * PsiG * CosT) * Esp + (b * PhiG * SinT) * Esp -
               (U * alpha * rho * PhiG * SinT) * Esp - (Phi2G * SinT) * Esp -
               (V * alpha * rho * PsiG * SinT) * Esp - (Psi2G * SinT) * Esp;

      Phi4 = -a / alpha / alpha *
             log(NumRe2 * NumRe2 / Den / Den + NumIm2 * NumIm2 / Den / Den);
      Psi4 = -2. * a / alpha / alpha * atan(NumIm2 / NumRe2);

      *DRePsi += Phi3 + Phi4;
      *DImPsi += Psi3 + Psi4;
      *DDRePsi = 0.0;
      *DDImPsi = 0.0;
    }

    break;

  case VANNA:

    x1 = exp(x1); // exp(x1) is the fwd.
    *DRePsi = 1. / x1 * U;
    *DImPsi = 1. / x1 * V;

    x2 = 2. * sqrt(x2);

    NumRe =
        -(b * U * x2) + b * U2 * x2 - b * V2 * x2 + U2 * x2 * alpha * rho -
        U3 * x2 * alpha * rho + V2 * x2 * alpha * rho -
        U * V2 * x2 * alpha * rho - U * x2 * PhiG + U2 * x2 * PhiG -
        V2 * x2 * PhiG - V * x2 * PsiG + 2 * U * V * x2 * PsiG +
        (2 * b * U * x2 * CosT) * Esp - (2 * b * U2 * x2 * CosT) * Esp +
        (2 * b * V2 * x2 * CosT) * Esp -
        (2 * U2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U3 * x2 * alpha * rho * CosT) * Esp -
        (2 * V2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U * V2 * x2 * alpha * rho * CosT) * Esp -
        (b * U * x2 * Cos2T) * Esp2 + (b * U2 * x2 * Cos2T) * Esp2 -
        (b * V2 * x2 * Cos2T) * Esp2 + (U2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U3 * x2 * alpha * rho * Cos2T) * Esp2 +
        (V2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Cos2T) * Esp2 +
        (U * x2 * PhiG * Cos2T) * Esp2 - (U2 * x2 * PhiG * Cos2T) * Esp2 +
        (V2 * x2 * PhiG * Cos2T) * Esp2 + (V * x2 * PsiG * Cos2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Cos2T) * Esp2 +
        (2 * V * x2 * PhiG * SinT) * Esp -
        (4 * U * V * x2 * PhiG * SinT) * Esp -
        (2 * U * x2 * PsiG * SinT) * Esp + (2 * U2 * x2 * PsiG * SinT) * Esp -
        (2 * V2 * x2 * PsiG * SinT) * Esp - (b * U * x2 * Sin2T) * Esp2 +
        (b * U2 * x2 * Sin2T) * Esp2 - (b * V2 * x2 * Sin2T) * Esp2 +
        (U2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U3 * x2 * alpha * rho * Sin2T) * Esp2 +
        (V2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Sin2T) * Esp2 +
        (U * x2 * PhiG * Sin2T) * Esp2 - (U2 * x2 * PhiG * Sin2T) * Esp2 +
        (V2 * x2 * PhiG * Sin2T) * Esp2 + (V * x2 * PsiG * Sin2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Sin2T) * Esp2;

    NumIm = -(b * V * x2) + 2 * b * U * V * x2 - U2 * V * x2 * alpha * rho -
            V3 * x2 * alpha * rho - V * x2 * PhiG + 2 * U * V * x2 * PhiG +
            U * x2 * PsiG - U2 * x2 * PsiG + V2 * x2 * PsiG +
            (2 * b * V * x2 * CosT) * Esp - (4 * b * U * V * x2 * CosT) * Esp +
            (2 * U2 * V * x2 * alpha * rho * CosT) * Esp +
            (2 * V3 * x2 * alpha * rho * CosT) * Esp -
            (b * V * x2 * Cos2T) * Esp2 + (2 * b * U * V * x2 * Cos2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Cos2T) * Esp2 -
            (V3 * x2 * alpha * rho * Cos2T) * Esp2 +
            (V * x2 * PhiG * Cos2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Cos2T) * Esp2 -
            (U * x2 * PsiG * Cos2T) * Esp2 + (U2 * x2 * PsiG * Cos2T) * Esp2 -
            (V2 * x2 * PsiG * Cos2T) * Esp2 - (2 * U * x2 * PhiG * SinT) * Esp +
            (2 * U2 * x2 * PhiG * SinT) * Esp -
            (2 * V2 * x2 * PhiG * SinT) * Esp -
            (2 * V * x2 * PsiG * SinT) * Esp +
            (4 * U * V * x2 * PsiG * SinT) * Esp - (b * V * x2 * Sin2T) * Esp2 +
            (2 * b * U * V * x2 * Sin2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Sin2T) * Esp2 -
            (V3 * x2 * alpha * rho * Sin2T) * Esp2 +
            (V * x2 * PhiG * Sin2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Sin2T) * Esp2 -
            (U * x2 * PsiG * Sin2T) * Esp2 + (U2 * x2 * PsiG * Sin2T) * Esp2 -
            (V2 * x2 * PsiG * Sin2T) * Esp2;

    DenRe = b - U * alpha * rho + PhiG - (b * CosT) * Esp +
            (U * alpha * rho * CosT) * Esp + (PhiG * CosT) * Esp +
            (V * alpha * rho * SinT) * Esp + (PsiG * SinT) * Esp;

    DenIm = -(V * alpha * rho) + PsiG + (V * alpha * rho * CosT) * Esp +
            (PsiG * CosT) * Esp + (b * SinT) * Esp -
            (U * alpha * rho * SinT) * Esp - (PhiG * SinT) * Esp;

    *DDRePsi = NumRe / (DenRe * DenRe + DenIm * DenIm);
    *DDImPsi = NumIm / (DenRe * DenRe + DenIm * DenIm);

    /*

    This part is needed if one is interested in the case when Sigma(0) =
    Sigma(infinity)
    */
    if (isVolInfFix) {

      a = b * x2;
      Phi3 = T / alpha / alpha * (a * b - a * U * alpha * rho - a * PhiG);
      Psi3 = T / alpha / alpha * (-a * V * alpha * rho - a * PsiG);
      Chi = PhiG * PhiG + PsiG * PsiG;
      Den = 2 * Chi;

      NumRe2 = b * PhiG - U * alpha * rho * PhiG - Phi2G + 2 * Chi -
               V * alpha * rho * PsiG - Psi2G - (b * PhiG * CosT) * Esp +
               (U * alpha * rho * PhiG * CosT) * Esp + (Phi2G * CosT) * Esp +
               (V * alpha * rho * PsiG * CosT) * Esp + (Psi2G * CosT) * Esp +
               (V * alpha * rho * PhiG * SinT) * Esp + (b * PsiG * SinT) * Esp -
               (U * alpha * rho * PsiG * SinT) * Esp;

      NumIm2 = -(V * alpha * rho * PhiG) - b * PsiG + U * alpha * rho * PsiG +
               (V * alpha * rho * PhiG * CosT) * Esp + (b * PsiG * CosT) * Esp -
               (U * alpha * rho * PsiG * CosT) * Esp + (b * PhiG * SinT) * Esp -
               (U * alpha * rho * PhiG * SinT) * Esp - (Phi2G * SinT) * Esp -
               (V * alpha * rho * PsiG * SinT) * Esp - (Psi2G * SinT) * Esp;

      Phi4 = -a / alpha / alpha *
             log(NumRe2 * NumRe2 / Den / Den + NumIm2 * NumIm2 / Den / Den);
      Psi4 = -2. * a / alpha / alpha * atan(NumIm2 / NumRe2);

      *DDRePsi += Phi3 + Phi4;
      *DDImPsi += Psi3 + Psi4;
    }

    break;

  case VOLGA:

    x2 = 2. * sqrt(x2);

    NumRe =
        -(b * U * x2) + b * U2 * x2 - b * V2 * x2 + U2 * x2 * alpha * rho -
        U3 * x2 * alpha * rho + V2 * x2 * alpha * rho -
        U * V2 * x2 * alpha * rho - U * x2 * PhiG + U2 * x2 * PhiG -
        V2 * x2 * PhiG - V * x2 * PsiG + 2 * U * V * x2 * PsiG +
        (2 * b * U * x2 * CosT) * Esp - (2 * b * U2 * x2 * CosT) * Esp +
        (2 * b * V2 * x2 * CosT) * Esp -
        (2 * U2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U3 * x2 * alpha * rho * CosT) * Esp -
        (2 * V2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U * V2 * x2 * alpha * rho * CosT) * Esp -
        (b * U * x2 * Cos2T) * Esp2 + (b * U2 * x2 * Cos2T) * Esp2 -
        (b * V2 * x2 * Cos2T) * Esp2 + (U2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U3 * x2 * alpha * rho * Cos2T) * Esp2 +
        (V2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Cos2T) * Esp2 +
        (U * x2 * PhiG * Cos2T) * Esp2 - (U2 * x2 * PhiG * Cos2T) * Esp2 +
        (V2 * x2 * PhiG * Cos2T) * Esp2 + (V * x2 * PsiG * Cos2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Cos2T) * Esp2 +
        (2 * V * x2 * PhiG * SinT) * Esp -
        (4 * U * V * x2 * PhiG * SinT) * Esp -
        (2 * U * x2 * PsiG * SinT) * Esp + (2 * U2 * x2 * PsiG * SinT) * Esp -
        (2 * V2 * x2 * PsiG * SinT) * Esp - (b * U * x2 * Sin2T) * Esp2 +
        (b * U2 * x2 * Sin2T) * Esp2 - (b * V2 * x2 * Sin2T) * Esp2 +
        (U2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U3 * x2 * alpha * rho * Sin2T) * Esp2 +
        (V2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Sin2T) * Esp2 +
        (U * x2 * PhiG * Sin2T) * Esp2 - (U2 * x2 * PhiG * Sin2T) * Esp2 +
        (V2 * x2 * PhiG * Sin2T) * Esp2 + (V * x2 * PsiG * Sin2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Sin2T) * Esp2;

    NumIm = -(b * V * x2) + 2 * b * U * V * x2 - U2 * V * x2 * alpha * rho -
            V3 * x2 * alpha * rho - V * x2 * PhiG + 2 * U * V * x2 * PhiG +
            U * x2 * PsiG - U2 * x2 * PsiG + V2 * x2 * PsiG +
            (2 * b * V * x2 * CosT) * Esp - (4 * b * U * V * x2 * CosT) * Esp +
            (2 * U2 * V * x2 * alpha * rho * CosT) * Esp +
            (2 * V3 * x2 * alpha * rho * CosT) * Esp -
            (b * V * x2 * Cos2T) * Esp2 + (2 * b * U * V * x2 * Cos2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Cos2T) * Esp2 -
            (V3 * x2 * alpha * rho * Cos2T) * Esp2 +
            (V * x2 * PhiG * Cos2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Cos2T) * Esp2 -
            (U * x2 * PsiG * Cos2T) * Esp2 + (U2 * x2 * PsiG * Cos2T) * Esp2 -
            (V2 * x2 * PsiG * Cos2T) * Esp2 - (2 * U * x2 * PhiG * SinT) * Esp +
            (2 * U2 * x2 * PhiG * SinT) * Esp -
            (2 * V2 * x2 * PhiG * SinT) * Esp -
            (2 * V * x2 * PsiG * SinT) * Esp +
            (4 * U * V * x2 * PsiG * SinT) * Esp - (b * V * x2 * Sin2T) * Esp2 +
            (2 * b * U * V * x2 * Sin2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Sin2T) * Esp2 -
            (V3 * x2 * alpha * rho * Sin2T) * Esp2 +
            (V * x2 * PhiG * Sin2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Sin2T) * Esp2 -
            (U * x2 * PsiG * Sin2T) * Esp2 + (U2 * x2 * PsiG * Sin2T) * Esp2 -
            (V2 * x2 * PsiG * Sin2T) * Esp2;

    DenRe = b - U * alpha * rho + PhiG - (b * CosT) * Esp +
            (U * alpha * rho * CosT) * Esp + (PhiG * CosT) * Esp +
            (V * alpha * rho * SinT) * Esp + (PsiG * SinT) * Esp;

    DenIm = -(V * alpha * rho) + PsiG + (V * alpha * rho * CosT) * Esp +
            (PsiG * CosT) * Esp + (b * SinT) * Esp -
            (U * alpha * rho * SinT) * Esp - (PhiG * SinT) * Esp;

    *DRePsi = NumRe / (DenRe * DenRe + DenIm * DenIm);
    *DImPsi = NumIm / (DenRe * DenRe + DenIm * DenIm);

    /*
    This part is needed if one is interested in the case when Sigma(0) =
    Sigma(infinity)
    */
    if (isVolInfFix) {

      a = b * x2;
      Phi3 = T / alpha / alpha * (a * b - a * U * alpha * rho - a * PhiG);
      Psi3 = T / alpha / alpha * (-a * V * alpha * rho - a * PsiG);

      Chi = PhiG * PhiG + PsiG * PsiG;
      Den = 2 * Chi;

      NumRe2 = b * PhiG - U * alpha * rho * PhiG - Phi2G + 2 * Chi -
               V * alpha * rho * PsiG - Psi2G - (b * PhiG * CosT) * Esp +
               (U * alpha * rho * PhiG * CosT) * Esp + (Phi2G * CosT) * Esp +
               (V * alpha * rho * PsiG * CosT) * Esp + (Psi2G * CosT) * Esp +
               (V * alpha * rho * PhiG * SinT) * Esp + (b * PsiG * SinT) * Esp -
               (U * alpha * rho * PsiG * SinT) * Esp;

      NumIm2 = -(V * alpha * rho * PhiG) - b * PsiG + U * alpha * rho * PsiG +
               (V * alpha * rho * PhiG * CosT) * Esp + (b * PsiG * CosT) * Esp -
               (U * alpha * rho * PsiG * CosT) * Esp + (b * PhiG * SinT) * Esp -
               (U * alpha * rho * PhiG * SinT) * Esp - (Phi2G * SinT) * Esp -
               (V * alpha * rho * PsiG * SinT) * Esp - (Psi2G * SinT) * Esp;

      Phi4 = -a / alpha / alpha *
             log(NumRe2 * NumRe2 / Den / Den + NumIm2 * NumIm2 / Den / Den);
      Psi4 = -2. * a / alpha / alpha * atan(NumIm2 / NumRe2);

      *DRePsi += Phi3 + Phi4;
      *DImPsi += Psi3 + Psi4;
    }

    x2 = 2.;
    NumRe =
        -(b * U * x2) + b * U2 * x2 - b * V2 * x2 + U2 * x2 * alpha * rho -
        U3 * x2 * alpha * rho + V2 * x2 * alpha * rho -
        U * V2 * x2 * alpha * rho - U * x2 * PhiG + U2 * x2 * PhiG -
        V2 * x2 * PhiG - V * x2 * PsiG + 2 * U * V * x2 * PsiG +
        (2 * b * U * x2 * CosT) * Esp - (2 * b * U2 * x2 * CosT) * Esp +
        (2 * b * V2 * x2 * CosT) * Esp -
        (2 * U2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U3 * x2 * alpha * rho * CosT) * Esp -
        (2 * V2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U * V2 * x2 * alpha * rho * CosT) * Esp -
        (b * U * x2 * Cos2T) * Esp2 + (b * U2 * x2 * Cos2T) * Esp2 -
        (b * V2 * x2 * Cos2T) * Esp2 + (U2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U3 * x2 * alpha * rho * Cos2T) * Esp2 +
        (V2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Cos2T) * Esp2 +
        (U * x2 * PhiG * Cos2T) * Esp2 - (U2 * x2 * PhiG * Cos2T) * Esp2 +
        (V2 * x2 * PhiG * Cos2T) * Esp2 + (V * x2 * PsiG * Cos2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Cos2T) * Esp2 +
        (2 * V * x2 * PhiG * SinT) * Esp -
        (4 * U * V * x2 * PhiG * SinT) * Esp -
        (2 * U * x2 * PsiG * SinT) * Esp + (2 * U2 * x2 * PsiG * SinT) * Esp -
        (2 * V2 * x2 * PsiG * SinT) * Esp - (b * U * x2 * Sin2T) * Esp2 +
        (b * U2 * x2 * Sin2T) * Esp2 - (b * V2 * x2 * Sin2T) * Esp2 +
        (U2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U3 * x2 * alpha * rho * Sin2T) * Esp2 +
        (V2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Sin2T) * Esp2 +
        (U * x2 * PhiG * Sin2T) * Esp2 - (U2 * x2 * PhiG * Sin2T) * Esp2 +
        (V2 * x2 * PhiG * Sin2T) * Esp2 + (V * x2 * PsiG * Sin2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Sin2T) * Esp2;

    NumIm = -(b * V * x2) + 2 * b * U * V * x2 - U2 * V * x2 * alpha * rho -
            V3 * x2 * alpha * rho - V * x2 * PhiG + 2 * U * V * x2 * PhiG +
            U * x2 * PsiG - U2 * x2 * PsiG + V2 * x2 * PsiG +
            (2 * b * V * x2 * CosT) * Esp - (4 * b * U * V * x2 * CosT) * Esp +
            (2 * U2 * V * x2 * alpha * rho * CosT) * Esp +
            (2 * V3 * x2 * alpha * rho * CosT) * Esp -
            (b * V * x2 * Cos2T) * Esp2 + (2 * b * U * V * x2 * Cos2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Cos2T) * Esp2 -
            (V3 * x2 * alpha * rho * Cos2T) * Esp2 +
            (V * x2 * PhiG * Cos2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Cos2T) * Esp2 -
            (U * x2 * PsiG * Cos2T) * Esp2 + (U2 * x2 * PsiG * Cos2T) * Esp2 -
            (V2 * x2 * PsiG * Cos2T) * Esp2 - (2 * U * x2 * PhiG * SinT) * Esp +
            (2 * U2 * x2 * PhiG * SinT) * Esp -
            (2 * V2 * x2 * PhiG * SinT) * Esp -
            (2 * V * x2 * PsiG * SinT) * Esp +
            (4 * U * V * x2 * PsiG * SinT) * Esp - (b * V * x2 * Sin2T) * Esp2 +
            (2 * b * U * V * x2 * Sin2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Sin2T) * Esp2 -
            (V3 * x2 * alpha * rho * Sin2T) * Esp2 +
            (V * x2 * PhiG * Sin2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Sin2T) * Esp2 -
            (U * x2 * PsiG * Sin2T) * Esp2 + (U2 * x2 * PsiG * Sin2T) * Esp2 -
            (V2 * x2 * PsiG * Sin2T) * Esp2;

    DenRe = b - U * alpha * rho + PhiG - (b * CosT) * Esp +
            (U * alpha * rho * CosT) * Esp + (PhiG * CosT) * Esp +
            (V * alpha * rho * SinT) * Esp + (PsiG * SinT) * Esp;

    DenIm = -(V * alpha * rho) + PsiG + (V * alpha * rho * CosT) * Esp +
            (PsiG * CosT) * Esp + (b * SinT) * Esp -
            (U * alpha * rho * SinT) * Esp - (PhiG * SinT) * Esp;

    *DDRePsi = NumRe / (DenRe * DenRe + DenIm * DenIm);
    *DDImPsi = NumIm / (DenRe * DenRe + DenIm * DenIm);

    /*
    This part is needed if one is interested in the case when Sigma(0) =
    Sigma(infinity)
    */
    if (isVolInfFix) {

      a = 2. * b;
      Phi3 = T / alpha / alpha * (a * b - a * U * alpha * rho - a * PhiG);
      Psi3 = T / alpha / alpha * (-a * V * alpha * rho - a * PsiG);
      Chi = PhiG * PhiG + PsiG * PsiG;
      Den = 2 * Chi;

      NumRe2 = b * PhiG - U * alpha * rho * PhiG - Phi2G + 2 * Chi -
               V * alpha * rho * PsiG - Psi2G - (b * PhiG * CosT) * Esp +
               (U * alpha * rho * PhiG * CosT) * Esp + (Phi2G * CosT) * Esp +
               (V * alpha * rho * PsiG * CosT) * Esp + (Psi2G * CosT) * Esp +
               (V * alpha * rho * PhiG * SinT) * Esp + (b * PsiG * SinT) * Esp -
               (U * alpha * rho * PsiG * SinT) * Esp;

      NumIm2 = -(V * alpha * rho * PhiG) - b * PsiG + U * alpha * rho * PsiG +
               (V * alpha * rho * PhiG * CosT) * Esp + (b * PsiG * CosT) * Esp -
               (U * alpha * rho * PsiG * CosT) * Esp + (b * PhiG * SinT) * Esp -
               (U * alpha * rho * PhiG * SinT) * Esp - (Phi2G * SinT) * Esp -
               (V * alpha * rho * PsiG * SinT) * Esp - (Psi2G * SinT) * Esp;

      Phi4 = -a / alpha / alpha *
             log(NumRe2 * NumRe2 / Den / Den + NumIm2 * NumIm2 / Den / Den);
      Psi4 = -2. * a / alpha / alpha * atan(NumIm2 / NumRe2);

      *DDRePsi += Phi3 + Phi4;
      *DDImPsi += Psi3 + Psi4;
    }

    break;

  case THETA:

    /* L2 term contribution */

    NumRe =
        -(b * U * x2) + b * U2 * x2 - b * V2 * x2 + U2 * x2 * alpha * rho -
        U3 * x2 * alpha * rho + V2 * x2 * alpha * rho -
        U * V2 * x2 * alpha * rho - U * x2 * PhiG + U2 * x2 * PhiG -
        V2 * x2 * PhiG - V * x2 * PsiG + 2 * U * V * x2 * PsiG +
        (2 * b * U * x2 * CosT) * Esp - (2 * b * U2 * x2 * CosT) * Esp +
        (2 * b * V2 * x2 * CosT) * Esp -
        (2 * U2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U3 * x2 * alpha * rho * CosT) * Esp -
        (2 * V2 * x2 * alpha * rho * CosT) * Esp +
        (2 * U * V2 * x2 * alpha * rho * CosT) * Esp -
        (b * U * x2 * Cos2T) * Esp2 + (b * U2 * x2 * Cos2T) * Esp2 -
        (b * V2 * x2 * Cos2T) * Esp2 + (U2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U3 * x2 * alpha * rho * Cos2T) * Esp2 +
        (V2 * x2 * alpha * rho * Cos2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Cos2T) * Esp2 +
        (U * x2 * PhiG * Cos2T) * Esp2 - (U2 * x2 * PhiG * Cos2T) * Esp2 +
        (V2 * x2 * PhiG * Cos2T) * Esp2 + (V * x2 * PsiG * Cos2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Cos2T) * Esp2 +
        (2 * V * x2 * PhiG * SinT) * Esp -
        (4 * U * V * x2 * PhiG * SinT) * Esp -
        (2 * U * x2 * PsiG * SinT) * Esp + (2 * U2 * x2 * PsiG * SinT) * Esp -
        (2 * V2 * x2 * PsiG * SinT) * Esp - (b * U * x2 * Sin2T) * Esp2 +
        (b * U2 * x2 * Sin2T) * Esp2 - (b * V2 * x2 * Sin2T) * Esp2 +
        (U2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U3 * x2 * alpha * rho * Sin2T) * Esp2 +
        (V2 * x2 * alpha * rho * Sin2T) * Esp2 -
        (U * V2 * x2 * alpha * rho * Sin2T) * Esp2 +
        (U * x2 * PhiG * Sin2T) * Esp2 - (U2 * x2 * PhiG * Sin2T) * Esp2 +
        (V2 * x2 * PhiG * Sin2T) * Esp2 + (V * x2 * PsiG * Sin2T) * Esp2 -
        (2 * U * V * x2 * PsiG * Sin2T) * Esp2;

    NumIm = -(b * V * x2) + 2 * b * U * V * x2 - U2 * V * x2 * alpha * rho -
            V3 * x2 * alpha * rho - V * x2 * PhiG + 2 * U * V * x2 * PhiG +
            U * x2 * PsiG - U2 * x2 * PsiG + V2 * x2 * PsiG +
            (2 * b * V * x2 * CosT) * Esp - (4 * b * U * V * x2 * CosT) * Esp +
            (2 * U2 * V * x2 * alpha * rho * CosT) * Esp +
            (2 * V3 * x2 * alpha * rho * CosT) * Esp -
            (b * V * x2 * Cos2T) * Esp2 + (2 * b * U * V * x2 * Cos2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Cos2T) * Esp2 -
            (V3 * x2 * alpha * rho * Cos2T) * Esp2 +
            (V * x2 * PhiG * Cos2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Cos2T) * Esp2 -
            (U * x2 * PsiG * Cos2T) * Esp2 + (U2 * x2 * PsiG * Cos2T) * Esp2 -
            (V2 * x2 * PsiG * Cos2T) * Esp2 - (2 * U * x2 * PhiG * SinT) * Esp +
            (2 * U2 * x2 * PhiG * SinT) * Esp -
            (2 * V2 * x2 * PhiG * SinT) * Esp -
            (2 * V * x2 * PsiG * SinT) * Esp +
            (4 * U * V * x2 * PsiG * SinT) * Esp - (b * V * x2 * Sin2T) * Esp2 +
            (2 * b * U * V * x2 * Sin2T) * Esp2 -
            (U2 * V * x2 * alpha * rho * Sin2T) * Esp2 -
            (V3 * x2 * alpha * rho * Sin2T) * Esp2 +
            (V * x2 * PhiG * Sin2T) * Esp2 -
            (2 * U * V * x2 * PhiG * Sin2T) * Esp2 -
            (U * x2 * PsiG * Sin2T) * Esp2 + (U2 * x2 * PsiG * Sin2T) * Esp2 -
            (V2 * x2 * PsiG * Sin2T) * Esp2;

    DNumRe =
        (2 * x2 *
         (PhiG *
              (b * U - b * U2 + b * V2 - U2 * alpha * rho + U3 * alpha * rho -
               V2 * alpha * rho + U * V2 * alpha * rho - U * PhiG + U2 * PhiG -
               V2 * PhiG - V * PsiG + 2 * U * V * PsiG) *
              Esp2 +

          (b * (-U + U2 - V2) * PhiG + U2 * alpha * rho * PhiG -
           U3 * alpha * rho * PhiG + V2 * alpha * rho * PhiG -
           U * V2 * alpha * rho * PhiG + V * PhiG * PsiG -
           2 * U * V * PhiG * PsiG - U * Psi2G + U2 * Psi2G - V2 * Psi2G) *
              Esp * CosT +

          Esp *
              (-(V * Phi2G) + 2 * U * V * Phi2G - b * U * PsiG + b * U2 * PsiG -
               b * V2 * PsiG + U2 * alpha * rho * PsiG -
               U3 * alpha * rho * PsiG + V2 * alpha * rho * PsiG -
               U * V2 * alpha * rho * PsiG + U * PhiG * PsiG -
               U2 * PhiG * PsiG + V2 * PhiG * PsiG) *
              SinT)); /* OK*/

    DNumIm =
        (2 * x2 *
         (PhiG *
              (b * V - 2 * b * U * V + U2 * V * alpha * rho + V3 * alpha * rho -
               V * PhiG + 2 * U * V * PhiG + U * PsiG - U2 * PsiG + V2 * PsiG) *
              Esp2 +

          Esp *
              (b * (-1 + 2 * U) * V * PhiG - U2 * V * alpha * rho * PhiG -
               V3 * alpha * rho * PhiG - U * PhiG * PsiG + U2 * PhiG * PsiG -
               V2 * PhiG * PsiG - V * Psi2G + 2 * U * V * Psi2G) *
              CosT -

          Esp *
              (-(U * Phi2G) + U2 * Phi2G - V2 * Phi2G + b * V * PsiG -
               2 * b * U * V * PsiG + U2 * V * alpha * rho * PsiG +
               V3 * alpha * rho * PsiG - V * PhiG * PsiG +
               2 * U * V * PhiG * PsiG) *
              SinT)); /* OK */

    DDen =
        (-2 *
         (PhiG *
              (b * b - 2 * b * U * alpha * rho +
               U2 * alpha * alpha * rho * rho + V2 * alpha * alpha * rho * rho -
               2 * b * PhiG + 2 * U * alpha * rho * PhiG + Phi2G +
               2 * V * alpha * rho * PsiG + Psi2G) *
              Esp2 +

          Esp *
              (-(b * b * PhiG) + 2 * b * U * alpha * rho * PhiG -
               U2 * alpha * alpha * rho * rho * PhiG -
               V2 * alpha * alpha * rho * rho * PhiG + Phi3G -
               2 * V * alpha * rho * PhiG * PsiG - 2 * b * Psi2G +
               2 * U * alpha * rho * Psi2G + PhiG * Psi2G) *
              CosT +

          Esp *
              (2 * V * alpha * rho * Phi2G - b * b * PsiG +
               2 * b * U * alpha * rho * PsiG -
               U2 * alpha * alpha * rho * rho * PsiG -
               V2 * alpha * alpha * rho * rho * PsiG + 2 * b * PhiG * PsiG -
               2 * U * alpha * rho * PhiG * PsiG + Phi2G * PsiG + Psi3G) *
              SinT)); /* OK */

    DenRe = b - U * alpha * rho + PhiG - (b * CosT) * Esp +
            (U * alpha * rho * CosT) * Esp + (PhiG * CosT) * Esp +
            (V * alpha * rho * SinT) * Esp + (PsiG * SinT) * Esp;

    DenIm = -(V * alpha * rho) + PsiG + (V * alpha * rho * CosT) * Esp +
            (PsiG * CosT) * Esp + (b * SinT) * Esp -
            (U * alpha * rho * SinT) * Esp - (PhiG * SinT) * Esp;

    Den = DenRe * DenRe + DenIm * DenIm;

    Phi2 = (DNumRe * Den - DDen * NumRe) / Den / Den; /* OK */
    Psi2 = (DNumIm * Den - DDen * NumIm) / Den / Den; /* OK */

    /* L3 term contribution */

    Phi3 =
        1. / alpha / alpha *
        (r * alpha * alpha * (U - 1) + a * b - a * U * alpha * rho - a * PhiG);
    Psi3 = 1. / alpha / alpha *
           (r * V * alpha * alpha - a * V * alpha * rho - a * PsiG);

    /* L4 term contribution */

    Chi = PhiG * PhiG + PsiG * PsiG;
    Den = 2 * Chi;

    NumRe = b * PhiG - U * alpha * rho * PhiG - Phi2G + 2 * Chi -
            V * alpha * rho * PsiG - Psi2G - (b * PhiG * CosT) * Esp +
            (U * alpha * rho * PhiG * CosT) * Esp + (Phi2G * CosT) * Esp +
            (V * alpha * rho * PsiG * CosT) * Esp + (Psi2G * CosT) * Esp +
            (V * alpha * rho * PhiG * SinT) * Esp + (b * PsiG * SinT) * Esp -
            (U * alpha * rho * PsiG * SinT) * Esp; /* OK */

    NumIm = -(V * alpha * rho * PhiG) - b * PsiG + U * alpha * rho * PsiG +
            (V * alpha * rho * PhiG * CosT) * Esp + (b * PsiG * CosT) * Esp -
            (U * alpha * rho * PsiG * CosT) * Esp + (b * PhiG * SinT) * Esp -
            (U * alpha * rho * PhiG * SinT) * Esp - (Phi2G * SinT) * Esp -
            (V * alpha * rho * PsiG * SinT) * Esp -
            (Psi2G * SinT) * Esp; /* OK */

    DLogArg =
        -(PhiG *
              (b * b - 2 * b * U * alpha * rho +
               U2 * alpha * alpha * rho * rho + V2 * alpha * alpha * rho * rho -
               2 * b * PhiG + 2 * U * alpha * rho * PhiG + Phi2G +
               2 * V * alpha * rho * PsiG + Psi2G) *
              Esp2 +
          Esp *
              (-(b * b * PhiG) + 2 * b * U * alpha * rho * PhiG -
               U2 * alpha * alpha * rho * rho * PhiG -
               V2 * alpha * alpha * rho * rho * PhiG + Phi3G -
               2 * V * alpha * rho * PhiG * PsiG - 2 * b * Psi2G +
               2 * U * alpha * rho * Psi2G + PhiG * Psi2G) *
              CosT +
          Esp *
              (2 * V * alpha * rho * Phi2G - b * b * PsiG +
               2 * b * U * alpha * rho * PsiG -
               U2 * alpha * alpha * rho * rho * PsiG -
               V2 * alpha * alpha * rho * rho * PsiG + 2 * b * PhiG * PsiG -
               2 * U * alpha * rho * PhiG * PsiG + Phi2G * PsiG + Psi3G) *
              SinT) /
        (2. * (Phi2G + Psi2G)); /* OK */

    Phi4 = -a / alpha / alpha /
           (NumRe * NumRe / Den / Den + NumIm * NumIm / Den / Den) * DLogArg;

    DNumRe = -((Phi2G + Psi2G) * ((-b + U * alpha * rho + PhiG) * CosT +
                                  (V * alpha * rho + PsiG) * SinT)) *
             Esp; /*OK */
    DNumIm = ((Phi2G + Psi2G) * (-((V * alpha * rho + PsiG) * CosT) +
                                 (-b + U * alpha * rho + PhiG) * SinT)) *
             Esp; /* OK */
    Psi4 = (DNumIm * NumRe - DNumRe * NumIm) / NumRe / NumRe;
    Psi4 = -2 * a / alpha / alpha * Psi4 * 1. /
           (1 + NumIm * NumIm / NumRe / NumRe);

    /* computes the output */

    *DRePsi = Phi2 + Phi3 + Phi4;
    *DImPsi = Psi2 + Psi3 + Psi4;
    *DDRePsi = 0.0;
    *DDImPsi = 0.0;

    break;

  case DENSITY:

    *DRePsi = 0.0;
    *DImPsi = 0.0;
    *DDRePsi = 0.0;
    *DDImPsi = 0.0;

    break;
  }

  return err;
}

/*---------------------------------------------------------------------------------------------------

  This function returns the location of the node of the PsiFunction
(Arrow-Debreau function) at the right of the last found onde that has to be
input as  "NodeLeft"

----------------------------------------------------------------------------------------------------*/

#define NRANSI
#define MAXIT 60
#define UNUSED (-1.11e30)

double HestonFindRightNode(double ImAppr, double Strike, double NodeLeft,
                           SRT_Boolean isVolInfFix, double Maturity, double a,
                           double b, double Alpha, double Rho, double Forward,
                           double Sigma, double U, SrtGreekType greek) {
  double Node = 0.0, x1, x2, RePsi, ImPsi, DRePsi, DImPsi, DDRePsi, DDImPsi;
  int j;
  double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew, xacc = 1.e-7;
  Err err = NULL;

  /* finds x1  , the left boundary */

  x1 = NodeLeft + fabs(SRT_PI / 2. / (ImAppr + log(Strike)));
  x2 = NodeLeft + fabs(SRT_PI * 3. / 2. / (ImAppr + log(Strike)));

  err = PsiFunction(isVolInfFix, Maturity, U, -x1, a, b, Alpha, Rho,
                    log(Forward), Sigma * Sigma, greek, &RePsi, &ImPsi, &DRePsi,
                    &DImPsi, &DDRePsi, &DDImPsi);
  fl = sin(fabs(ImPsi + x1 * log(Strike)));

  err = PsiFunction(isVolInfFix, Maturity, U, -x2, a, b, Alpha, Rho,
                    log(Forward), Sigma * Sigma, greek, &RePsi, &ImPsi, &DRePsi,
                    &DImPsi, &DDRePsi, &DDImPsi);
  fh = sin(fabs(ImPsi + x2 * log(Strike)));

  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl = x1;
    xh = x2;
    ans = UNUSED;
    for (j = 1; j <= MAXIT; j++) {
      xm = 0.5 * (xl + xh);

      err = PsiFunction(isVolInfFix, Maturity, U, -xm, a, b, Alpha, Rho,
                        log(Forward), Sigma * Sigma, greek, &RePsi, &ImPsi,
                        &DRePsi, &DImPsi, &DDRePsi, &DDImPsi);
      fm = sin(fabs(ImPsi + xm * log(Strike)));

      s = sqrt(fm * fm - fl * fh);
      if (s == 0.0)
        return ans;
      xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
      if (fabs(xnew - ans) <= xacc)
        return ans;
      ans = xnew;

      err = PsiFunction(isVolInfFix, Maturity, U, -ans, a, b, Alpha, Rho,
                        log(Forward), Sigma * Sigma, greek, &RePsi, &ImPsi,
                        &DRePsi, &DImPsi, &DDRePsi, &DDImPsi);
      fnew = sin(fabs(ImPsi + ans * log(Strike)));

      if (fnew == 0.0)
        return ans;
      if (SIGN(fm, fnew) != fm) {
        xl = xm;
        fl = fm;
        xh = ans;
        fh = fnew;
      } else if (SIGN(fl, fnew) != fl) {
        xh = ans;
        fh = fnew;
      } else if (SIGN(fh, fnew) != fh) {
        xl = ans;
        fl = fnew;
      } else {
        return 0.0;
      }
      if (fabs(xh - xl) <= xacc)
        return ans;
    }
    return 0.0;
    // nrerror("zriddr exceed maximum iterations");
  } else {
    if (fl == 0.0)
      return x1;
    if (fh == 0.0)
      return x2;
    // nrerror("root must be bracketed in zriddr.");
  }

  return 0.0;
}

#undef MAXIT
#undef UNUSED
#undef NRANSI