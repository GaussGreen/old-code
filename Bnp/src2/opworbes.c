
#include <stdio.h>

#include "math.h"
#include "num_h_gausslegendre.h"
#include "num_h_proba.h"
#include "opfnctns.h"
#include "spfnctns.h"
#include "utallhdr.h"

#define GAUSS_LIMIT 6.0
#define precision 1e-6
#define DOUBLEMAX 300
#define DOUBLEMIN 1e-12;

/* Some Static Declarations	*/

static double         s_FirstStrike, s_SecondStrike;
static double         s_FirstBSVol, s_SecondBSVol;
static double         s_Gearing;
static double         s_MeanLnFirstForward, s_MeanLnSecondForward, s_Rho, s_Maturity;
static double         s_VlimitInf, s_VlimitSup;
static long           s_NumPoints;
static SrtCallPutType s_Call_PutOne, s_Call_PutTwo;

/* Computaion of the Simpson integrand. The Simpson integrand is a U variable function.
 It is the result of a Gauss Legendre integration.  */

static double simp_func(
    double         U,
    double         FirstBSVol,
    double         SecondBSVol,
    double         Gearing,
    double         MeanLnFirstForward,
    double         MeanLnSecondForward,
    double         Rho,
    double         Maturity,
    long           NumPoints,
    SrtCallPutType Call_PutOne,
    SrtCallPutType Call_PutTwo)

{
    double *GaussLegWeights, *GaussLegAbscissas;
    int     i;
    double  X, YlimitSup, YlimitInf, VlimitSup, VlimitInf;
    double  term, exp_term;
    double  result;

    GaussLegWeights   = dvector(1, NumPoints);
    GaussLegAbscissas = dvector(1, NumPoints);

    if ((Call_PutOne == SRT_CALL) && (Call_PutTwo == SRT_CALL))
    {
        /* Computation of the boundaries of the Gauss Legendre integral */

        X         = exp(U * sqrt(Maturity) * FirstBSVol + MeanLnFirstForward);
        YlimitSup = (X + s_Gearing * s_SecondStrike - s_FirstStrike) / s_Gearing;
        if (YlimitSup < 0)
            VlimitSup = -GAUSS_LIMIT;
        else
        {
            VlimitSup = (log(YlimitSup) - MeanLnSecondForward) / (sqrt(Maturity) * SecondBSVol);
            VlimitSup = DMIN(GAUSS_LIMIT, VlimitSup);
        }
        /* Some initialisation */

        result   = 0.0;
        term     = 0.0;
        exp_term = 1.0;

        /* Computation of the Gauss Legendre integral */

        if (VlimitSup > s_VlimitInf)
        {
            gauleg(s_VlimitInf, VlimitSup, GaussLegAbscissas, GaussLegWeights, NumPoints);

            for (i = 1; i <= NumPoints; i++)
            {
                term = -0.5 * GaussLegAbscissas[i] * GaussLegAbscissas[i] / (1 - Rho * Rho) +
                       (Rho / (1 - Rho * Rho)) * U * GaussLegAbscissas[i];

                if (fabs(term) < DOUBLEMAX)
                {
                    exp_term = exp(term);
                    result += GaussLegWeights[i] * exp_term;
                }

                else
                    result += DOUBLEMIN;
            }
        }

        else
            result = DOUBLEMIN;
    }

    else if ((Call_PutOne == SRT_PUT) && (Call_PutTwo == SRT_PUT))

    {
        /* Computation of the boundaries of the Gauss Legendre integral */

        X         = exp(U * sqrt(Maturity) * FirstBSVol + MeanLnFirstForward);
        YlimitInf = (X + s_Gearing * s_SecondStrike - s_FirstStrike) / s_Gearing;
        if (YlimitInf < 0)
            VlimitInf = -GAUSS_LIMIT;

        else
        {
            VlimitInf = (log(YlimitInf) - MeanLnSecondForward) / (sqrt(Maturity) * SecondBSVol);
            VlimitInf = DMAX(-GAUSS_LIMIT, VlimitInf);
        }
        /* Some initialisation */
        result   = 0.0;
        term     = 0;
        exp_term = 1.0;

        /* Computation of the Gauss Legendre integral */

        if (s_VlimitSup > VlimitInf)
        {
            gauleg(VlimitInf, s_VlimitSup, GaussLegAbscissas, GaussLegWeights, NumPoints);

            for (i = 1; i <= NumPoints; i++)
            {
                term = -0.5 * GaussLegAbscissas[i] * GaussLegAbscissas[i] / (1 - Rho * Rho) +
                       (Rho / (1 - Rho * Rho)) * U * GaussLegAbscissas[i];

                if (fabs(term) < DOUBLEMAX)
                {
                    exp_term = exp(term);
                    result += GaussLegWeights[i] * exp_term;
                }

                else
                    result += DOUBLEMIN;
            }
        }

        else
            result = DOUBLEMIN;
    }

    else if ((Call_PutOne == SRT_CALL) && (Call_PutTwo == SRT_PUT))
    {
        /* Computation of the boundaries of the Gauss Legendre integral */

        X         = exp(U * sqrt(Maturity) * FirstBSVol + MeanLnFirstForward);
        YlimitInf = (s_Gearing * s_SecondStrike + s_FirstStrike - X) / s_Gearing;
        if (YlimitInf < 0)
            VlimitInf = -GAUSS_LIMIT;

        else
        {
            VlimitInf = (log(YlimitInf) - MeanLnSecondForward) / (sqrt(Maturity) * SecondBSVol);
            VlimitInf = DMAX(-GAUSS_LIMIT, VlimitInf);
        }
        /* Some initialisation */

        result   = 0.0;
        term     = 0.0;
        exp_term = 1.0;

        /* Computation of the Gauss Legendre integral */

        if (s_VlimitSup > VlimitInf)
        {
            gauleg(VlimitInf, s_VlimitSup, GaussLegAbscissas, GaussLegWeights, NumPoints);

            for (i = 1; i <= NumPoints; i++)
            {
                term = -0.5 * GaussLegAbscissas[i] * GaussLegAbscissas[i] / (1 - Rho * Rho) +
                       (Rho / (1 - Rho * Rho)) * U * GaussLegAbscissas[i];

                if (fabs(term) < DOUBLEMAX)
                {
                    exp_term = exp(term);
                    result += GaussLegWeights[i] * exp_term;
                }

                else
                    result += DOUBLEMIN;
            }
        }

        else
            result = DOUBLEMIN;
    }

    result = result * exp(-(0.5 / (1 - Rho * Rho)) * U * U) / (2 * SRT_PI * sqrt(1 - Rho * Rho));

    return result;
}

static double simp_integrand(double U)

{
    double result;
    result = simp_func(
        U,
        s_FirstBSVol,
        s_SecondBSVol,
        s_Gearing,
        s_MeanLnFirstForward,
        s_MeanLnSecondForward,
        s_Rho,
        s_Maturity,
        s_NumPoints,
        s_Call_PutOne,
        s_Call_PutTwo);

    return result;
}

/* The main fonction. Return the min (max) of two Call or two Put Options */

double srt_f_optworstbest(
    double           FirstForward,
    double           SecondForward,
    double           FirstStrike,
    double           SecondStrike,
    double           Gearing,
    double           FirstBSVol,
    double           SecondBSVol,
    double           Rho,
    double           Df,
    double           Maturity,
    SrtCallPutType   Call_PutOne,
    SrtCallPutType   Call_PutTwo,
    SrtBestWorstType Best_Worst,
    SrtGreekType     Greek)

{
    double MeanLnFirstForwardp, MeanLnSecondForwardp;
    double MeanLnFirstForwardpp, MeanLnSecondForwardpp;
    double SqFirstBSVol, SqSecondBSVol;
    double SqrtMaturity;
    double A, B, Ap, Bp, App, Bpp;
    double UlimitInf, UlimitSup;
    double Premium;
    double FirstOptionPrice, SecondOptionPrice;
    double answer, shift;

    s_FirstBSVol   = FirstBSVol;
    s_SecondBSVol  = SecondBSVol;
    s_Gearing      = Gearing;
    s_NumPoints    = 10;
    s_FirstStrike  = FirstStrike;
    s_SecondStrike = SecondStrike;
    SqFirstBSVol   = FirstBSVol * FirstBSVol;
    SqSecondBSVol  = SecondBSVol * SecondBSVol;

    /* computation of the mean of the logarithm of the forward FX rate under QT (see paper) */

    s_MeanLnFirstForward  = log(FirstForward) - 0.5 * Maturity * SqFirstBSVol;
    s_MeanLnSecondForward = log(SecondForward) - 0.5 * Maturity * SqSecondBSVol;

    s_Rho         = Rho;
    s_Maturity    = Maturity;
    s_Call_PutOne = Call_PutOne;
    s_Call_PutTwo = Call_PutTwo;

    SqrtMaturity = sqrt(Maturity);

    /* computation of the mean of the logarithm of the forward FX rate under QX1 (see paper) */

    MeanLnFirstForwardp  = log(FirstForward) + 0.5 * Maturity * SqFirstBSVol;
    MeanLnSecondForwardp = log(SecondForward) - 0.5 * Maturity * SqSecondBSVol +
                           Rho * FirstBSVol * SecondBSVol * Maturity;

    /* computation of the mean of the logarithm of the forward FX rate under QX2 (see paper) */

    MeanLnFirstForwardpp = log(FirstForward) - 0.5 * Maturity * SqFirstBSVol +
                           Rho * FirstBSVol * SecondBSVol * Maturity;
    MeanLnSecondForwardpp = log(SecondForward) + 0.5 * Maturity * SqSecondBSVol;

    A = (log(FirstStrike) - s_MeanLnFirstForward) / (SqrtMaturity * s_FirstBSVol);
    B = (log(SecondStrike) - s_MeanLnSecondForward) / (SqrtMaturity * s_SecondBSVol);

    Ap = (log(FirstStrike) - MeanLnFirstForwardp) / (SqrtMaturity * s_FirstBSVol);
    Bp = (log(SecondStrike) - MeanLnSecondForwardp) / (SqrtMaturity * s_SecondBSVol);

    App = (log(FirstStrike) - MeanLnFirstForwardpp) / (SqrtMaturity * s_FirstBSVol);
    Bpp = (log(SecondStrike) - MeanLnSecondForwardpp) / (SqrtMaturity * s_SecondBSVol);

    Premium = 0.0;

    /* The first part of the price */

    if (Call_PutOne == SRT_CALL)
    {
        Premium = Df * srt_f_optblksch(
                           FirstForward, FirstStrike, FirstBSVol, Maturity, 1.0, SRT_CALL, PREMIUM);
    }

    if (Call_PutOne == SRT_PUT)
    {
        Premium = Df * srt_f_optblksch(
                           FirstForward, FirstStrike, FirstBSVol, Maturity, 1.0, SRT_PUT, PREMIUM);
    }
    /* The second part of the price */

    if ((Call_PutOne == SRT_CALL) && (Call_PutTwo == SRT_CALL))
    {
        Premium += FirstStrike * Df * (norm_accurate(B) - bivar(A, B, Rho));
        Premium -= FirstForward * Df * (norm_accurate(Bp) - bivar(Ap, Bp, Rho));
    }

    else if ((Call_PutOne == SRT_PUT) && (Call_PutTwo == SRT_PUT))
    {
        Premium -= FirstStrike * Df * (norm_accurate(A) - bivar(B, A, Rho));
        Premium += FirstForward * Df * (norm_accurate(Ap) - bivar(Bp, Ap, Rho));
    }

    else if ((Call_PutOne == SRT_CALL) && (Call_PutTwo == SRT_PUT))
    {
        Premium += FirstStrike * Df * (1 + bivar(A, B, Rho) - norm(A) - norm(B));
        Premium -= FirstForward * Df * (1 + bivar(Ap, Bp, Rho) - norm(Ap) - norm(Bp));
    }

    /* The integral part of the pricing  */

    if ((Call_PutOne == SRT_CALL) && (Call_PutTwo == SRT_CALL))
    {
        UlimitSup   = GAUSS_LIMIT;
        UlimitInf   = A;
        s_VlimitInf = B;

        /* Computaion of the probability under QT */
        Premium -= Df * (s_Gearing * SecondStrike - FirstStrike) *
                   sm_qsimp(simp_integrand, UlimitInf, UlimitSup, precision);

        UlimitInf   = Ap;
        s_VlimitInf = Bp;

        s_MeanLnFirstForward  = MeanLnFirstForwardp;
        s_MeanLnSecondForward = MeanLnSecondForwardp;

        /* Computaion of the probability under QX1 */

        Premium -= Df * FirstForward * sm_qsimp(simp_integrand, UlimitInf, UlimitSup, precision);

        UlimitInf   = App;
        s_VlimitInf = Bpp;

        s_MeanLnFirstForward  = MeanLnFirstForwardpp;
        s_MeanLnSecondForward = MeanLnSecondForwardpp;

        /* Computation of the probability under QX2 */
        Premium += Df * s_Gearing * SecondForward *
                   sm_qsimp(simp_integrand, UlimitInf, UlimitSup, precision);
    }

    else if ((Call_PutOne == SRT_PUT) && (Call_PutTwo == SRT_PUT))
    {
        UlimitInf   = -GAUSS_LIMIT;
        UlimitSup   = A;
        s_VlimitSup = B;

        Premium += Df * (s_Gearing * SecondStrike - FirstStrike) *
                   sm_qsimp(simp_integrand, UlimitInf, UlimitSup, precision);

        UlimitSup   = Ap;
        s_VlimitSup = Bp;

        s_MeanLnFirstForward  = MeanLnFirstForwardp;
        s_MeanLnSecondForward = MeanLnSecondForwardp;

        Premium += Df * FirstForward * sm_qsimp(simp_integrand, UlimitInf, UlimitSup, precision);

        UlimitSup   = App;
        s_VlimitSup = Bpp;

        s_MeanLnFirstForward  = MeanLnFirstForwardpp;
        s_MeanLnSecondForward = MeanLnSecondForwardpp;

        Premium -= Df * s_Gearing * SecondForward *
                   sm_qsimp(simp_integrand, UlimitInf, UlimitSup, precision);
    }

    else if ((Call_PutOne == SRT_CALL) && (Call_PutTwo == SRT_PUT))
    {
        UlimitSup   = GAUSS_LIMIT;
        UlimitInf   = A;
        s_VlimitSup = B;

        Premium += Df * (s_Gearing * SecondStrike + FirstStrike) *
                   sm_qsimp(simp_integrand, UlimitInf, UlimitSup, precision);

        UlimitSup   = GAUSS_LIMIT;
        UlimitInf   = Ap;
        s_VlimitSup = Bp;

        s_MeanLnFirstForward  = MeanLnFirstForwardp;
        s_MeanLnSecondForward = MeanLnSecondForwardp;

        Premium -= Df * FirstForward * sm_qsimp(simp_integrand, UlimitInf, UlimitSup, precision);

        UlimitSup   = GAUSS_LIMIT;
        UlimitInf   = App;
        s_VlimitSup = Bpp;

        s_MeanLnFirstForward  = MeanLnFirstForwardpp;
        s_MeanLnSecondForward = MeanLnSecondForwardpp;

        Premium -= Df * s_Gearing * SecondForward *
                   sm_qsimp(simp_integrand, UlimitInf, UlimitSup, precision);
    }

    switch (Greek)
    {
    case PREMIUM:

        switch (Best_Worst)
        {
        case SRT_BESTOFF:

            FirstOptionPrice =
                Df *
                srt_f_optblksch(
                    FirstForward, FirstStrike, FirstBSVol, Maturity, 1.0, Call_PutOne, PREMIUM);

            SecondOptionPrice =
                Df * s_Gearing *
                srt_f_optblksch(
                    SecondForward, SecondStrike, SecondBSVol, Maturity, 1.0, Call_PutTwo, PREMIUM);

            Premium = FirstOptionPrice + SecondOptionPrice - Premium;

            return (Premium);
            break;

        case SRT_WORSTOFF:
            return (Premium);
            break;
        }

        return (Premium);
        break;

    case DELTAX:
        shift  = FirstForward / 10000;
        answer = (srt_f_optworstbest(
                      FirstForward + shift,
                      SecondForward,
                      FirstStrike,
                      SecondStrike,
                      Gearing,
                      FirstBSVol,
                      SecondBSVol,
                      Rho,
                      Df,
                      Maturity,
                      Call_PutOne,
                      Call_PutTwo,
                      Best_Worst,
                      PREMIUM) -
                  Premium) /
                 shift;
        return answer;
        break;

    case DELTAY:
        shift  = SecondForward / 10000;
        answer = (srt_f_optworstbest(
                      FirstForward,
                      SecondForward + shift,
                      FirstStrike,
                      SecondStrike,
                      Gearing,
                      FirstBSVol,
                      SecondBSVol,
                      Rho,
                      Df,
                      Maturity,
                      Call_PutOne,
                      Call_PutTwo,
                      Best_Worst,
                      PREMIUM) -
                  Premium) /
                 shift;

        return answer;
        break;

    case GAMMAX:
        shift  = FirstForward / 1000;
        answer = srt_f_optworstbest(
            FirstForward + shift,
            SecondForward,
            FirstStrike,
            SecondStrike,
            Gearing,
            FirstBSVol,
            SecondBSVol,
            Rho,
            Df,
            Maturity,
            Call_PutOne,
            Call_PutTwo,
            Best_Worst,
            PREMIUM);

        answer += srt_f_optworstbest(
            FirstForward - shift,
            SecondForward,
            FirstStrike,
            SecondStrike,
            Gearing,
            FirstBSVol,
            SecondBSVol,
            Rho,
            Df,
            Maturity,
            Call_PutOne,
            Call_PutTwo,
            Best_Worst,
            PREMIUM);

        answer -= 2 * Premium;
        return (answer);
        break;

    case GAMMAY:
        shift  = SecondForward / 1000;
        answer = srt_f_optworstbest(
            FirstForward,
            SecondForward + shift,
            FirstStrike,
            SecondStrike,
            Gearing,
            FirstBSVol,
            SecondBSVol,
            Rho,
            Df,
            Maturity,
            Call_PutOne,
            Call_PutTwo,
            Best_Worst,
            PREMIUM);

        answer += srt_f_optworstbest(
            FirstForward,
            SecondForward - shift,
            FirstStrike,
            SecondStrike,
            Gearing,
            FirstBSVol,
            SecondBSVol,
            Rho,
            Df,
            Maturity,
            Call_PutOne,
            Call_PutTwo,
            Best_Worst,
            PREMIUM);

        answer -= 2 * Premium;
        return (answer);
        break;

    case VEGAX:

        shift  = GVOPT.vol_add;
        answer = (srt_f_optworstbest(
                      FirstForward,
                      SecondForward,
                      FirstStrike,
                      SecondStrike,
                      Gearing,
                      FirstBSVol + shift,
                      SecondBSVol,
                      Rho,
                      Df,
                      Maturity,
                      Call_PutOne,
                      Call_PutTwo,
                      Best_Worst,
                      PREMIUM) -
                  Premium) /
                 shift;
        return (answer);
        break;

    case VEGAY:
        shift  = GVOPT.vol_add;
        answer = (srt_f_optworstbest(
                      FirstForward,
                      SecondForward,
                      FirstStrike,
                      SecondStrike,
                      Gearing,
                      FirstBSVol,
                      SecondBSVol + shift,
                      Rho,
                      Df,
                      Maturity,
                      Call_PutOne,
                      Call_PutTwo,
                      Best_Worst,
                      PREMIUM) -
                  Premium) /
                 shift;
        return (answer);
        break;

    case THETA:
        shift  = YEARS_IN_DAY;
        answer = (srt_f_optworstbest(
                      FirstForward,
                      SecondForward,
                      FirstStrike,
                      SecondStrike,
                      Gearing,
                      FirstBSVol,
                      SecondBSVol + shift,
                      Rho,
                      Df * exp(-shift * log(Df) / Maturity),
                      Maturity - shift,
                      Call_PutOne,
                      Call_PutTwo,
                      Best_Worst,
                      PREMIUM) -
                  Premium) /
                 shift;

        return (answer);
        break;

    default:
        return UNKNOWN_GREEK;
    }

    return Premium;
}
