/*	RainbowOpt.c
        Author: D. Mayevski
        Purpose: calculate E[min(nx*Xt - kX, ny*Yt - kY, 0)] with Xt, Yt lognormal (using 1D
   numerical integration)
*/

#include "RainbowOpt.h"

#include "Math.h"
#include "RandomGen.h"
#include "num_h_gausslegendre.h"
#include "num_h_hermite.h"
#include "uterror.h"

#define NPTS 16
#define NSTD 6.0
#define INV_SQRT2PI 0.398942280401433

//	Calculates E[min(nx*Xt - kX, ny*Yt - kY, 0)] with Xt, Yt lognormal (using 1D numerical
//integration)

Err OptRainbow(
    double  fwdx,
    double  nx,
    double  kx,
    double  sigx,
    double  fwdy,
    double  ny,
    double  ky,
    double  sigy,
    double  mat,
    double  rho,
    double* res)
{
    Err    err = NULL;
    double w[NPTS], u[NPTS];
    double sqrtT = sqrt(mat);
    double umin = -NSTD, umax = NSTD;
    double vx = sigx * sqrtT, vy = sigy * sqrtT, volBS = vy * sqrt(1.0 - rho * rho);
    double d0  = (log(kx / (nx * fwdx)) + vx * vx / 2.0) / vx;
    double sum = 0.0, fwdBS, strikeBS, BS, xt;
    int    i;

    if (fabs(ny) < 1e-8)
        return serror("nY is too close to zero");

    if (d0 > umin)
    {
        gauleg(umin, (d0 < umax ? d0 : umax), u - 1, w - 1, NPTS);

        for (i = 0; i < NPTS; i++)
        {
            xt       = fwdx * exp(vx * u[i] - vx * vx / 2.0);
            fwdBS    = fwdy * exp(vy * rho * u[i] - vy * vy * rho * rho / 2.0);
            strikeBS = (ky - kx + nx * xt) / ny;
            BS       = srt_f_optblksch(
                fwdBS, strikeBS, volBS, 1.0, 1.0, (ny > 0.0 ? SRT_PUT : SRT_CALL), PREMIUM);
            sum += w[i] * exp(-u[i] * u[i] / 2.0) * (kx - nx * xt + BS * fabs(ny));
        }
    }
    if (d0 < umax)
    {
        gauleg((d0 > umin ? d0 : umin), umax, u - 1, w - 1, NPTS);

        for (i = 0; i < NPTS; i++)
        {
            fwdBS    = fwdy * exp(vy * rho * u[i] - vy * vy * rho * rho / 2.0);
            strikeBS = ky / ny;
            BS       = srt_f_optblksch(
                fwdBS, strikeBS, volBS, 1.0, 1.0, (ny > 0.0 ? SRT_PUT : SRT_CALL), PREMIUM);
            sum += w[i] * exp(-u[i] * u[i] / 2.0) * BS * fabs(ny);
        }
    }
    *res = -INV_SQRT2PI * sum;

    return err;
}

#undef NPTS
#define NPTS 16
#define INV_SQRTPI 0.564189583547756
#define SQRT2 1.4142135623731

//	Calculates E[max(nx*Xt - ny*Yt - K, 0)] with Xt, Yt lognormal (using 1D numerical
//integration)

Err OptSpread(
    double  fwdx,
    double  nx,
    double  sigx,
    double  fwdy,
    double  ny,
    double  sigy,
    double  K,
    double  mat,
    double  rho,
    double* res)
{
    Err    err = NULL;
    double w[NPTS], u[NPTS];
    double sqrtT = sqrt(mat);
    double vx = sigx * sqrtT, vy = sigy * sqrtT, volBS = vy * sqrt(1.0 - rho * rho);
    double sum = 0.0, fwdBS, strikeBS, BS, xt;
    int    i;

    if (fabs(ny) < 1e-8)
        return serror("nY is too close to zero");

    err = gauss_hermite(u - 1, w - 1, NPTS);
    if (err)
        return err;
    for (i = 0; i < NPTS; i++)
        u[i] *= SQRT2;

    for (i = 0; i < NPTS; i++)
    {
        xt       = fwdx * exp(vx * u[i] - vx * vx / 2.0);
        fwdBS    = fwdy * exp(vy * rho * u[i] - vy * vy * rho * rho / 2.0);
        strikeBS = (nx * xt - K) / ny;
        BS       = srt_f_optblksch(
            fwdBS, strikeBS, volBS, 1.0, 1.0, (ny > 0.0 ? SRT_PUT : SRT_CALL), PREMIUM);
        sum += w[i] * BS * fabs(ny);
    }
    *res = INV_SQRTPI * sum;

    return err;
}

// Calculates E[max(nx*Xt - ny*Yt - K, 0)] with Xt, Yt following SABR distributions using MC
// Note: nsteps - number of time discretization steps PER YEAR

Err OptSpreadSabrMC(
    double  fwdx,
    double  nx,
    double  sigx,
    double  alphax,
    double  betax,
    double  rhox,
    double  fwdy,
    double  ny,
    double  sigy,
    double  alphay,
    double  betay,
    double  rhoy,
    double  K,
    double  mat,
    double  rho,
    long    npaths,
    long    nsteps,
    double* res,
    double* std)
{
    Err        err = NULL;
    long       i, j, k, m, N = (long)(nsteps * mat + 1e-5), seed = -12345678;
    double     dt = mat / N, sqrtdt = sqrt(dt);
    SRandomGen rg;
    double **  sv = NULL, **corr = NULL, **chol = NULL, br[4], cbr[4], *pvs = NULL;

    memset(&rg, 0, sizeof(SRandomGen));
    err = ABS_Init(&rg, seed, npaths, 4, 0);
    if (err)
        goto FREE_RETURN;

    sv   = dmatrix(0, npaths - 1, 0, 3);
    corr = dmatrix(0, 3, 0, 3);
    chol = dmatrix(0, 3, 0, 3);
    pvs  = dvector(0, npaths - 1);
    if (!sv || !corr || !chol || !pvs)
    {
        err = serror("Memory failure in OptSpreadSabrMC");
        goto FREE_RETURN;
    }

    // Initialize state variables at time 0 for all paths
    for (j = 0; j < npaths; j++)
    {
        sv[j][0] = fwdx;
        sv[j][1] = sigx;
        sv[j][2] = fwdy;
        sv[j][3] = sigy;
    }

    memset(&corr[0][0], 0, 16 * sizeof(double));
    for (k = 0; k < 4; k++)
        corr[k][k] = 1.0;
    corr[1][0] = corr[0][1] = rhox;
    corr[2][0] = corr[0][2] = rho;
    corr[3][2] = corr[2][3] = rhoy;

    // Cholesky decomposition
    err = choldc(4, corr, chol);
    if (err)
        goto FREE_RETURN;

    // Move forward through time
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < npaths; j++)
        {
            // Generate independent brownian increments
            for (k = 0; k < 4; k++)
            {
                err = rg.Gauss(&rg, &br[k]);
                if (err)
                    goto FREE_RETURN;
                br[k] *= sqrtdt;
            }
            // Correlate the increments
            memset(cbr, 0, 4 * sizeof(double));
            for (k = 0; k < 4; k++)
                for (m = 0; m <= k; m++)
                    cbr[k] += chol[k][m] * br[m];

            // Calculate state variables at the next time step
            sv[j][0] += sv[j][1] * pow(fabs(sv[j][0]), betax) * cbr[0];
            sv[j][1] += sv[j][1] * alphax * cbr[1];
            sv[j][2] += sv[j][3] * pow(fabs(sv[j][2]), betay) * cbr[2];
            sv[j][3] += sv[j][3] * alphay * cbr[3];

#ifdef _DEBUG
            if (sv[j][0] < 0.0)
                smessage("SV[0] < 0.0. Please check !!!");
            if (sv[j][1] < 0.0)
                smessage("SV[1] < 0.0. Please check !!!");
            if (sv[j][2] < 0.0)
                smessage("SV[2] < 0.0. Please check !!!");
            if (sv[j][3] < 0.0)
                smessage("SV[3] < 0.0. Please check !!!");
#endif  // #ifdef _DEBUG
        }
    }

    // Calculate PVs for all paths and the expectation
    *res = 0.0;
    for (j = 0; j < npaths; j++)
    {
        pvs[j] = nx * sv[j][0] - ny * sv[j][2] - K;
        if (pvs[j] < 0.0)
            pvs[j] = 0.0;
        *res += pvs[j] / npaths;
    }
    // Calculate std
    *std = 0.0;
    for (j = 0; j < npaths; j++)
        *std += (pvs[j] - *res) * (pvs[j] - *res) / (npaths - 1);
    *std = sqrt(*std / npaths);

FREE_RETURN:

    ABS_Free(&rg);
    if (sv)
        free_dmatrix(sv, 0, npaths - 1, 0, 3);
    if (corr)
        free_dmatrix(corr, 0, 3, 0, 3);
    if (chol)
        free_dmatrix(chol, 0, 3, 0, 3);
    if (pvs)
        free_dvector(pvs, 0, npaths - 1);
    return err;
}

/*	----------------------------------------------------------------------------------------------------
        if call
                calculates E[max((nx*Xt + ny*Yt) - K, 0)] with Xt, Yt lognormal (using 1D numerical
   integration) if put calculates E[max(K - (nx*Xt + ny*Yt), 0)] with Xt, Yt lognormal (using 1D
   numerical integration)

        ----------------------------------------------------------------------------------------------------
 */

#define SwapVariables(Type, A, B) \
    {                             \
        Type C;                   \
        C = A;                    \
        A = B;                    \
        B = C;                    \
    }

Err OptSpreadNew(
    double         fwdx,
    double         nx,
    double         sigx,
    double         fwdy,
    double         ny,
    double         sigy,
    double         K,
    double         mat,
    double         rho,
    SrtCallPutType call_put,
    double*        res)
{
    /* Declaration of local variables */
    Err    err = NULL;
    double w[NPTS], u[NPTS];

    double sum, sqrtT, vx, vy;
    double volBS, fwdBS, strikeBS;
    int    i;

    /* Check of inputs */
    if ((rho > 1) || (rho < -1))
        return err = " rho should be between -1 and 1";

    /* Case  mat <0 */
    if (mat < 0)
    {
        *res = 0;
    }
    else
    {
        if ((nx != 0) && (ny != 0))
        {
            /* Integration on the underlying with the less normal vol		*/
            /* By default x is the underlying on which we integrate			*/
            /* swap x and y if necessary
             */
            /* if fwd<0 then srt_f_optblksch return value with 0 normal vol	*/
            if (fabs(nx * max(fwdx, 0) * sigx) > fabs(ny * max(fwdy, 0) * sigy))
            {
                /* swap x and y */
                SwapVariables(double, nx, ny);
                SwapVariables(double, fwdx, fwdy);
                SwapVariables(double, sigx, sigy);
            }

            /* init parameters*/
            sqrtT = sqrt(mat);
            vx    = sigx * sqrtT;
            vy    = sigy * sqrtT;
            volBS = vy * sqrt(1.0 - rho * rho);

            /* Calculate the abscissas and weights of the n-point Gauss-Hermite quadrature */
            if (err = hermite_gauss_quick(NPTS, u - 1, w - 1))
                return err;

            for (i = 0; i < NPTS; i++)
                u[i] *= SQRT2;

            sum = 0;
            for (i = 0; i < NPTS; i++)
            {
                fwdBS    = fwdy * exp(vy * rho * u[i] - vy * vy * rho * rho / 2.0);
                strikeBS = (K - nx * fwdx * exp(vx * u[i] - vx * vx / 2.0)) / ny;
                sum += w[i] * fabs(ny) *
                       srt_f_optblksch(
                           fwdBS,
                           strikeBS,
                           volBS,
                           1.0,
                           1.0,
                           (ny > 0.0 ? call_put : (call_put == SRT_CALL ? SRT_PUT : SRT_CALL)),
                           PREMIUM);
            }
            *res = INV_SQRTPI * sum;
        }
        else
        {
            if ((nx == 0) && (ny == 0))
            {
                *res = (call_put == SRT_CALL ? max(-K, 0) : max(K, 0));
            }
            else
            {
                if (nx == 0)
                {
                    strikeBS = K / ny;
                    *res     = fabs(ny) *
                           srt_f_optblksch(
                               fwdy,
                               strikeBS,
                               sigy,
                               mat,
                               1.0,
                               (ny > 0.0 ? call_put : (call_put == SRT_CALL ? SRT_PUT : SRT_CALL)),
                               PREMIUM);
                }
                else
                {
                    /* ny = 0 */
                    strikeBS = K / nx;
                    *res     = fabs(nx) *
                           srt_f_optblksch(
                               fwdx,
                               strikeBS,
                               sigx,
                               mat,
                               1.0,
                               (nx > 0.0 ? call_put : (call_put == SRT_CALL ? SRT_PUT : SRT_CALL)),
                               PREMIUM);
                }
            }
        }
    }

    return err;
}