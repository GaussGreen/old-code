

#include "srt_h_resetable.h"

#include "grf_h_public.h"
#include "math.h"
#include "num_h_gamma.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_calib.h"
#include "srt_h_closedform.h"
#include "srt_h_grfclsdfrm.h"
#include "swp_h_all.h"
#include "utallhdr.h"

#define prec                                                                     \
    0.0001 /* Related to the approximation made when considering a finite number \
                                            of jumps*/

/* MODELIZES THE JUMP PART OF THE DYNAMIC */
struct Sjumps
{
    double  nsizes; /* Number of possible events */
    double* sizes;  /* For each event, value of Prod((1+U_i)^n_i) */
    double* probas; /* Conditional probabilities */
};

/* CALCULATES NUMBER OF p-UPLETS WHOSE SUM IS i, FOR
   i=0...max */
Err CalculeNbrespuplets(double* nbrespuplets, double max, int p)
/* The output nbrepuplets[i] is the number of p-uplets of integers whose sum is  i-1 */
{
    double** u;
    int      i, j, k;
    double   s;

    /* Memory allocation */
    u = dmatrix(1, (long)max + 1, 1, (long)p);

    /* Initialisation. */
    for (i = 1; i < max + 2; i++)
        u[i][1] = 1;

    /* Induction relation */
    for (j = 2; j < p + 1; j++)
    {
        u[1][j] = 1;
        for (i = 2; i < max + 2; i++)
        {
            s = 0;
            for (k = 1; k < i + 1; k++)
                s += u[k][j - 1];
            u[i][j] = s;
        }
    }
    for (i = 1; i < max + 2; i++)
        nbrespuplets[i] = u[i][p];

    /* Free memory */
    free_dmatrix(u, 1, (long)max + 1, 1, p);

    return NULL;
}

/* CALCULATE AN UPPER BOUND OF THE TOTAL NUMBER OF JUMPS ON THE MATURITY OF OPTION*/
Err FSautsMax(
    double** Poisson, /* Array of p Poissons */
    int      p,
    double*  NSautsMax)
{
    int    i;
    double intensmax, s;

    /* Calculates the largest intensity of Poissons */
    intensmax = Poisson[1][1];
    for (i = 2; i <= p; i++)
        if (Poisson[i][1] > intensmax)
            intensmax = Poisson[i][1];

    /* Calculates the upper bound */
    s = 0;
    i = 0;
    while (s < (1 - prec) * intensmax) /* Guaranteed that the relative error on the expectation
                                          is smaller than prec  */
    {
        s += exp(-intensmax) * i * exp(i * log(intensmax)) / fact((int)i);
        i += 1;
    }
    *NSautsMax = i;

    return NULL;
}

/* LOI DE PROBA D'UN p-UPLET DE POISSONS INDEPENDENTS  */
Err Probas_Poissons(
    double** Poissons, /* For i=1..p, Poisson[i]
has two parameters: intensity and amplitude */
    int     p,
    double* n_sauts,
    double* proba,           /* Conditional proba of the following event: for all i,
the i=th Poisson jumps n_sauts[i] times */
    double* amplitude_totale /* For each event, value of Prod((1+U_i)^n_i) */
)
{
    int i;

    *proba            = 1;
    *amplitude_totale = 1;
    for (i = 1; i < p + 1; i++)
    {
        *proba = *proba * exp(-Poissons[i][1] + n_sauts[i] * log(Poissons[i][1])) /
                 fact((int)n_sauts[i]);
        *amplitude_totale = *amplitude_totale * pow(1 + Poissons[i][2], (double)n_sauts[i]);
    }

    return NULL;
}

/* ITERATIVE FUNCTION*/
Err FLoop(double** enspuplets, double* pupletloc, int levelloc, int maxloc, int p, int* n)
{
    int i, j;

    for (i = 0; i <= maxloc; i++)
    {
        pupletloc[levelloc] = i;
        if (levelloc == p)
        {
            *n += 1;
            for (j = 1; j <= p; j++)
                enspuplets[*n][j] = pupletloc[j];
        }
        else
            FLoop(enspuplets, pupletloc, levelloc + 1, maxloc - i, p, n);
    }

    return NULL;
}

/* CALCULATES AN APPROXIMATION OF THE TIME DEPENDENT JUMP PROCESS */
Err srt_f_tooljumps(
    double** Poissons, /*  For i=1..p, Poisson[i]
has two parameters: intensity and amplitude*/
    int            p,
    double         NSautsMax, /* Upper bound ofthe total number of jumps*/
    struct Sjumps* sauts,     /*result of the function*/
    int*           nbretotalpuplets)
{
    double** enspuplets; /* Set of all p-uplets (n1,..,np) */
    double*  puplet;
    double*  nbrespuplets;
    int      n;
    /* Memory allocation */
    nbrespuplets = dvector(1, (long)NSautsMax + 1);
    CalculeNbrespuplets(nbrespuplets, NSautsMax, p);
    *nbretotalpuplets = 0;
    for (n = 1; n < (NSautsMax + 2); n++)
        *nbretotalpuplets += (int)nbrespuplets[n];
    enspuplets = dmatrix(1, *nbretotalpuplets, 1, p);
    puplet     = dvector(1, p);

    (*sauts).nsizes = *nbretotalpuplets;
    (*sauts).sizes  = dvector(1, *nbretotalpuplets);
    (*sauts).probas = dvector(1, *nbretotalpuplets);

    /* Call the iterative function which generates  p-uplets whose sum is <= max */
    n = 0;
    FLoop(enspuplets, puplet, 1, (int)NSautsMax, p, &n);

    /* Calculates conditional probas and resulting size of jumps for each event.*/
    for (n = 1; n <= *nbretotalpuplets; n++)
        Probas_Poissons(Poissons, p, enspuplets[n], &((*sauts).probas[n]), &((*sauts).sizes[n]));

    /* Free memory */
    free_dvector(nbrespuplets, 1, (long)NSautsMax + 1);
    free_dmatrix(enspuplets, 1, *nbretotalpuplets, 1, p);
    free_dvector(puplet, 1, p);

    return NULL;
}

Err optmertonpremiumtimedependent(
    double  dFwd,
    double  Strike,
    int     n_periods, /* The 6 following variables are vectors of length nperiods */
    double* sigma,
    double* U1,
    double* lambda1,
    double* U2,
    double* lambda2,
    double* mat,
    char*   Logornorm,
    double* answer)
{
    double**      Poissons;
    int           i;
    double        NSautsMax; /* Total number of jumps during the life of the option */
    int           nevents;   /* Total number of possible jump events */
    struct Sjumps sauts;     /* Describes the outcome of all jumps events.*/
    double        TermCompens;
    double        dFwdJump; /* Fwd used in a term of the sum in Merton's formula */
    double        maturity;
    double        sigma_eq;
    double        premium; /* Price of the option */

    /* Calculates the maturity of the option and the gaussian vol */
    maturity = 0;
    sigma_eq = 0;
    for (i = 1; i <= n_periods; i++)
    {
        maturity += mat[i];
        sigma_eq += sigma[i] * sigma[i] * mat[i];
    }
    sigma_eq = sqrt(sigma_eq / maturity);

    /* Creates 2*n_periods Poisson processes */
    Poissons = dmatrix(1, 2 * n_periods, 1, 2);
    for (i = 1; i <= n_periods; i++)
    {
        Poissons[2 * i - 1][1] = lambda1[i] * mat[i];
        Poissons[2 * i - 1][2] = U1[i];
        Poissons[2 * i][1]     = lambda2[i] * mat[i];
        Poissons[2 * i][2]     = U2[i];
    }

    /* Calculates the maximum number of jumps */
    FSautsMax(Poissons, 2 * n_periods, &NSautsMax);

    /* Calculates the different events and their outcomes */
    srt_f_tooljumps(Poissons, 2 * n_periods, NSautsMax, &sauts, &nevents);

    /* Calculates the term due to the compensation of jumps */
    TermCompens = 1;
    for (i = 1; i <= 2 * n_periods; i++)
        TermCompens = TermCompens * exp(-Poissons[i][2] * Poissons[i][1]);

    /* Calculates Merton's formula */
    premium = 0;
    for (i = 1; i <= nevents; i++)
    {
        dFwdJump = dFwd * TermCompens * sauts.sizes[i];
        premium += sauts.probas[i] *
                   srt_f_optblksch(dFwdJump, Strike, sigma_eq, maturity, 1, SRT_CALL, PREMIUM);
    }

    *answer = premium;
    /* Free memory */
    free_dmatrix(Poissons, 1, 2 * n_periods, 1, 2);
    free_dvector(sauts.sizes, 1, nevents);
    free_dvector(sauts.probas, 1, nevents);

    return NULL;
}

Err optmertonsmiletimedependent(
    double  dFwd,
    int     n_strikes,
    double* strikes,
    int     n_periods, /* The 6 following variables are vectors of length nperiods */
    double* param,
    double* mat,
    char*   Logornorm,
    double* impvols)

{
    double**      Poissons;
    int           i, j;
    double        NSautsMax; /* Total number of jumps during the life of the option */
    int           nevents;   /* Total number of possible jump events */
    struct Sjumps sauts;     /* Describes the outcome of all jumps events.*/
    double        TermCompens;
    double        dFwdJump; /* Fwd used in a term of the sum in Merton's formula */
    double        maturity;
    double        sigma_eq;
    double*       premium; /* Price of the option */
    double*       sigma;
    double*       U1;
    double*       lambda1;
    double*       U2;
    double*       lambda2;

    /*Memory allocation*/
    premium = dvector(1, n_strikes);
    sigma   = dvector(1, n_periods);
    U1      = dvector(1, n_periods);
    lambda1 = dvector(1, n_periods);
    U2      = dvector(1, n_periods);
    lambda2 = dvector(1, n_periods);

    /* puts the parameters in their respective vectors */
    for (i = 1; i <= n_periods; i++)
    {
        sigma[i]   = param[5 * (i - 1) + 1];
        U1[i]      = param[5 * (i - 1) + 2];
        lambda1[i] = param[5 * (i - 1) + 3];
        U2[i]      = param[5 * (i - 1) + 4];
        lambda2[i] = param[5 * i];
    }

    /* Calculates the maturity of the option and the gaussian vol */
    maturity = 0;
    sigma_eq = 0;
    for (i = 1; i <= n_periods; i++)
    {
        maturity += mat[i];
        sigma_eq += sigma[i] * sigma[i] * mat[i];
    }
    sigma_eq = sqrt(sigma_eq / maturity);

    /* Creates 2*n_periods Poisson processes */
    Poissons = dmatrix(1, 2 * n_periods, 1, 2);

    for (i = 1; i <= n_periods; i++)
    {
        Poissons[2 * i - 1][1] = lambda1[i] * mat[i];
        Poissons[2 * i - 1][2] = U1[i];
        Poissons[2 * i][1]     = lambda2[i] * mat[i];
        Poissons[2 * i][2]     = U2[i];
    }

    /* Calculates the maximum number of jumps */
    FSautsMax(Poissons, 2 * n_periods, &NSautsMax);

    /* Calculates the different events and their outcomes */
    srt_f_tooljumps(Poissons, 2 * n_periods, NSautsMax, &sauts, &nevents);

    /* Calculates the term due to the compensation of jumps */
    TermCompens = 1;
    for (i = 1; i <= 2 * n_periods; i++)
        TermCompens = TermCompens * exp(-Poissons[i][2] * Poissons[i][1]);

    /* Calculates Merton's formula */
    for (i = 1; i <= n_strikes; i++)
    {
        premium[i] = 0;
    }

    for (i = 1; i <= n_strikes; i++)
    {
        for (j = 1; j <= nevents; j++)
        {
            dFwdJump = dFwd * TermCompens * sauts.sizes[j];
            premium[i] +=
                sauts.probas[j] *
                srt_f_optblksch(dFwdJump, strikes[i], sigma_eq, maturity, 1, SRT_CALL, PREMIUM);
        }
    }

    for (i = 1; i <= n_strikes; i++)
    {
        srt_f_optimpvol(
            premium[i], dFwd, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &(impvols[i]));
    }

    /* Free memory */
    free_dmatrix(Poissons, 1, 2 * n_periods, 1, 2);
    free_dvector(sauts.sizes, 1, nevents);
    free_dvector(sauts.probas, 1, nevents);
    free_dvector(premium, 1, n_strikes);
    free_dvector(sigma, 1, n_periods);
    free_dvector(U1, 1, n_periods);
    free_dvector(lambda1, 1, n_periods);
    free_dvector(U2, 1, n_periods);
    free_dvector(lambda2, 1, n_periods);
    return NULL;
}

double optmertonpremium(
    double dFwd,
    double Strike,
    double sigma,
    double U1,
    double lambda1,
    double U2,
    double lambda2,
    double T,
    char*  Logornorm)
{
    double lambda1p;
    double lambda2p;
    double N1max;
    double N2max;
    double d;
    double d1;
    double d2;
    int    n1, n2;
    double premium1 = 0;
    double premium2 = 0;
    double premium;
    double proba1, proba1p, proba2, proba2p;
    double probainf;

    probainf = 0.000001;

    if (strcmp(Logornorm, "Normal") == 0)
    {
        N1max = lambda1 * T + 10 * sqrt(lambda1 * T);
        N2max = lambda2 * T + 10 * sqrt(lambda2 * T);

        for (n1 = 0; n1 < N1max; n1++)
        {
            for (n2 = 0; n2 < N2max; n2++)
            {
                proba1 = exp(-lambda1 * T) * pow(lambda1 * T, n1) / fact(n1);
                proba2 = exp(-lambda2 * T) * pow(lambda2 * T, n2) / fact(n2);

                if (proba1 * proba2 > probainf)
                {
                    d = Strike - dFwd - U1 * (n1 - lambda1 * T) - U2 * (n2 - lambda2 * T);

                    premium1 += proba1 * proba2 * sigma * T *
                                exp(-d * d / (2 * sigma * sigma * T)) / sqrt(2 * SRT_PI * T);

                    premium2 += proba1 * proba2 * d * norm(-d / sigma);
                }
            }
        }

        premium = premium1 - premium2;
    }

    else
    {
        lambda1p = lambda1 * (1 + U1);
        lambda2p = lambda2 * (1 + U2);

        N1max = lambda1p * T + 10 * sqrt(lambda1p * T);
        N2max = lambda2 * T + 10 * sqrt(lambda2 * T);

        for (n1 = 0; n1 < N1max; n1++)
        {
            for (n2 = 0; n2 < N2max; n2++)
            {
                proba1p = exp(-lambda1p * T) * pow(lambda1p * T, n1) / fact(n1);
                proba2p = exp(-lambda2p * T) * pow(lambda2p * T, n2) / fact(n2);

                proba1 = exp(-lambda1 * T) * pow(lambda1 * T, n1) / fact(n1);
                proba2 = exp(-lambda2 * T) * pow(lambda2 * T, n2) / fact(n2);

                if ((proba1p * proba2p > probainf) || (proba1 * proba2 > probainf))
                {
                    d1 = (log(dFwd / Strike) - U1 * lambda1 * T - U2 * lambda2 * T +
                          n1 * log(1 + U1) + n2 * log(1 + U2) + sigma * sigma * T / 2) /
                         (sigma * sqrt(T));
                    d2 = d1 - sigma * sqrt(T);

                    premium1 += proba1p * proba2p * norm(d1);

                    premium2 += proba1 * proba2 * norm(d2);
                }
            }
        }

        premium = dFwd * premium1 - Strike * premium2;
    }

    return premium;
}

Err optmerton(
    double       dFwd,
    double       Strike,
    double       sigma,
    double       U1,
    double       lambda1,
    double       U2,
    double       lambda2,
    double       maturity,
    char*        Logornorm,
    SrtGreekType Greek,
    double*      answer)

{
    double shift;

    switch (Greek)
    {
    case PREMIUM:
        *answer =
            optmertonpremium(dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm);

        break;
    case DELTA:
        shift = dFwd / 10000;
        *answer =
            (optmertonpremium(
                 dFwd + shift, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm) -
             optmertonpremium(dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm)) /
            shift;
        break;
    case GAMMA:
        shift = dFwd / 10000;

        *answer = optmertonpremium(
            dFwd + shift, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm);
        *answer += optmertonpremium(
            dFwd - shift, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm);
        *answer -= 2 * optmertonpremium(
                           dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm);
        *answer /= shift * shift;

        break;

    case VEGA:

        shift = sigma / 100;

        *answer =
            (optmertonpremium(
                 dFwd, Strike, sigma + shift, U1, lambda1, U2, lambda2, maturity, Logornorm) -
             optmertonpremium(dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm)) /
            shift;

        break;

    case THETA:

        shift = YEARS_IN_DAY;

        *answer =
            (optmertonpremium(
                 dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity - shift, Logornorm) -
             optmertonpremium(dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm)) /
            shift;

        break;
    case DU1:

        shift = U1 / 100;

        *answer =
            (optmertonpremium(
                 dFwd, Strike, sigma, U1 + shift, lambda1, U2, lambda2, maturity, Logornorm) -
             optmertonpremium(dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm)) /
            shift;
        break;

    case DLAMBDA1:

        shift = lambda1 / 100;

        *answer =
            (optmertonpremium(
                 dFwd, Strike, sigma, U1, lambda1 + shift, U2, lambda2, maturity, Logornorm) -
             optmertonpremium(dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm)) /
            shift;
        break;

    case DU2:

        shift = U2 / 100;

        *answer =
            (optmertonpremium(
                 dFwd, Strike, sigma, U1, lambda1, U2 - shift, lambda2, maturity, Logornorm) -
             optmertonpremium(dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm)) /
            shift;
        break;

    case DLAMBDA2:

        shift = lambda2 / 100;

        *answer =
            (optmertonpremium(
                 dFwd, Strike, sigma, U1, lambda1, U2, lambda2 + shift, maturity, Logornorm) -
             optmertonpremium(dFwd, Strike, sigma, U1, lambda1, U2, lambda2, maturity, Logornorm)) /
            shift;
        break;

    default:
        *answer = UNKNOWN_GREEK;
    }

    return NULL;
}

Err optmertonsmile(
    double  dFwd,
    int     n_strikes,
    double* strikes,
    double* parameters,
    double  maturity,
    char*   Logornorm,
    double* impvols)
{
    int    i;
    double price;

    for (i = 1; i <= n_strikes; i++)
    {
        price = optmertonpremium(
            dFwd,
            strikes[i],
            parameters[1],
            parameters[2],
            parameters[3],
            parameters[4],
            parameters[5],
            maturity,
            Logornorm);

        if (strcmp(Logornorm, "Normal") == 0)
        {
            srt_f_optimpvol(
                price, dFwd, strikes[i], maturity, 1.0, SRT_CALL, SRT_NORMAL, &(impvols[i]));
        }
        else
        {
            srt_f_optimpvol(
                price, dFwd, strikes[i], maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &(impvols[i]));
        }
    }

    return NULL;
}
