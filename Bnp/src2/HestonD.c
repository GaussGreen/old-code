// HestonD.c : New Heston implementation

#include "srt_h_all.h"
#undef SIGN
#include "math.h"
#include "nag.h"
#include "nag_stdlib.h"
#include "nagd01.h"
#include "nagd02.h"
//#include "nage04.h"
//#include "d_amoeba.h"
#include "HestonD.h"
#include "LGMSVClosedForm.h"
#include "intde2.h"
#include "opfnctns.h"

// Structure for communication with NAG Runge-Kutta integrator:
typedef struct _SHestonComm_ODE
{
    double sigma, alpha, rho, lam;
    double u_re, u_im;
    int    want[4];  // Flags "want derivatives" : { sigma, alpha, rho, lam }
} SHestonComm_ODE;

// Cache to store calculated integrand values:
typedef struct _SHestonCache SHestonCache;
struct _SHestonCache
{
    int           maxpts;  // Max number of stored points in cache
    int           cur_pt;  // Currently requested point
    double*       v;       // Stored points
    double***     fn;      // Stored integrand values: per v, per k, per fn
    SHestonCache* next;    // Linked list
};

// Structure for communication with a quadrature integrator
typedef struct _SHestonComm_InvFT
{
    double mat;
    double sigma, alpha, rho, lam;
    int    want[4];  // Flags "want derivatives" : { sigma, alpha, rho, lam }
    double v_max;    // Upper limit of integration

    int     nk;       // Number of strikes
    int     want_k;   // Current strike
    int     want_fn;  // Current function (integrand or its derivative)
    double* k;        // Strikes ( log(K+sh/FFX0+sh) )

    // Integrand cache:
    SHestonCache *c_head, *c_cur;

    // Counter for function calls:
    int count;

} SHestonComm_InvFT;

typedef struct _SHestonComm_Calib
{
    double mat, fwd;
    double sigma, alpha, rho, lam, shift;
    double k[2], K[2];          // Strikes ( log(K+sh/FFX0+sh) and original)
    double volATM, RR, BF;      // Market data
    double ar, lastRR, lastBF;  // temporary calibration variables
    int    count;               // Counter for function calls
} SHestonComm_Calib;

static void NAG_CALL EvalDerivatives(Integer neq, double t, double y[], double yp[], Nag_User* comm)
{
    SHestonComm_ODE* p = (SHestonComm_ODE*)comm->p;

    double *yg[4], *ypg[4];
    int     i, offs;
    double  c1, c2, c3;

    for (i = 0, offs = 4; i < 4; i++)
        if (p->want[i])
        {
            yg[i]  = y + offs;
            ypg[i] = yp + offs;
            offs += 4;
        }

    // dy1/dt, dy2/dt

    yp[0] = -p->lam * y[2];
    yp[1] = -p->lam * y[3];

    for (i = 0; i < 4; i++)
        if (p->want[i])
        {
            ypg[i][0] = -p->lam * yg[i][2];
            ypg[i][1] = -p->lam * yg[i][3];
        }
    if (p->want[3])
    {
        ypg[3][0] += -y[2];
        ypg[3][1] += -y[3];
    }

    // dy3/dt, dy4/dt

    c1 = p->alpha * p->alpha;
    c2 = p->rho * p->sigma * p->alpha;
    c3 = p->sigma * p->sigma;

    yp[2] = -0.5 * c1 * (y[2] * y[2] - y[3] * y[3]) + p->lam * y[2] +
            c2 * (p->u_re * y[3] + p->u_im * y[2]) +
            0.5 * c3 * (p->u_re * p->u_re - p->u_im * p->u_im - p->u_im);

    yp[3] = -c1 * y[2] * y[3] + p->lam * y[3] - c2 * (p->u_re * y[2] - p->u_im * y[3]) +
            0.5 * c3 * (p->u_re + 2.0 * p->u_re * p->u_im);

    for (i = 0; i < 4; i++)
        if (p->want[i])
        {
            ypg[i][2] = -c1 * (yg[i][2] * y[2] - yg[i][3] * y[3]) + p->lam * yg[i][2] +
                        c2 * (p->u_re * yg[i][3] + p->u_im * yg[i][2]);

            ypg[i][3] = -c1 * (yg[i][2] * y[3] + yg[i][3] * y[2]) + p->lam * yg[i][3] -
                        c2 * (p->u_re * yg[i][2] - p->u_im * yg[i][3]);
        }

    if (p->want[0])
    {
        ypg[0][2] += p->rho * p->alpha * (p->u_re * y[3] + p->u_im * y[2]) +
                     p->sigma * (p->u_re * p->u_re - p->u_im * p->u_im - p->u_im);
        ypg[0][3] += -p->rho * p->alpha * (p->u_re * y[2] - p->u_im * y[3]) +
                     p->sigma * (p->u_re + 2.0 * p->u_re * p->u_im);
    }
    if (p->want[1])
    {
        ypg[1][2] += -p->alpha * (y[2] * y[2] - y[3] * y[3]) +
                     p->rho * p->sigma * (p->u_re * y[3] + p->u_im * y[2]);
        ypg[1][3] +=
            -2.0 * p->alpha * y[2] * y[3] - p->rho * p->sigma * (p->u_re * y[2] - p->u_im * y[3]);
    }
    if (p->want[2])
    {
        ypg[2][2] += p->alpha * p->sigma * (p->u_re * y[3] + p->u_im * y[2]);
        ypg[2][3] += -p->alpha * p->sigma * (p->u_re * y[2] - p->u_im * y[3]);
    }
    if (p->want[3])
    {
        ypg[3][2] += y[2];
        ypg[3][3] += y[3];
    }
}

static Err HestonCalcPhi(SHestonComm_InvFT* q, double u_re, double u_im, double* h_re, double* h_im)
{
    Err             err = NULL;
    int             j, k, neq;
    SHestonComm_ODE comm_RK;
    Nag_User        comm_Nag;
    NagError        fail;
    Nag_ODE_RK      opt;
    double          y[20], yp[20], ymax[20], RKthres[20], tgot, norm;
    double *        yg[4], *ypg[4], yg_re, yg_im;
    const double    RKtol = 1e-5, thres_min = 1e-8, thres_coef = 1e-7;  // adjustable

    memset(&fail, 0, sizeof(NagError));
    memset(&opt, 0, sizeof(Nag_ODE_RK));

    for (k = 0, neq = 4; k < 4; k++)
        if (q->want[k])
        {
            yg[k]  = y + neq;
            ypg[k] = yp + neq;
            neq += 4;
        }

    // Preinitialize comm structure for NAG

    comm_Nag.p    = &comm_RK;
    comm_RK.u_re  = u_re;
    comm_RK.u_im  = u_im;
    comm_RK.sigma = q->sigma;
    comm_RK.alpha = q->alpha;
    comm_RK.rho   = q->rho;
    comm_RK.lam   = q->lam;
    memcpy(comm_RK.want, q->want, 4 * sizeof(int));

    memset(y, 0, neq * sizeof(double));  // Final y values are all zeros.
    for (j = 0; j < neq; j++)
        RKthres[j] = thres_min;

    // Calculate solution at time 0:

    nag_ode_ivp_rk_setup(
        neq,
        q->mat,
        y,
        0.0,
        RKtol,
        RKthres,
        Nag_RK_4_5,
        Nag_RK_range,
        Nag_ErrorAssess_off,
        0.0,
        &opt,
        &fail);
    if (fail.code != NE_NOERROR)
    {
        err = serror(fail.message);
        goto FREE_RETURN;
    }

    nag_ode_ivp_rk_range(neq, EvalDerivatives, 0.0, &tgot, y, yp, ymax, &opt, &comm_Nag, &fail);
    if (fail.code != NE_NOERROR)
    {
        err = serror(fail.message);
        goto FREE_RETURN;
    }

    norm    = exp(y[0] + y[2]);
    h_re[0] = norm * cos(y[1] + y[3]);
    h_im[0] = norm * sin(y[1] + y[3]);

    for (k = 0; k < 4; k++)
        if (q->want[k])
        {
            yg_re       = yg[k][0] + yg[k][2];
            yg_im       = yg[k][1] + yg[k][3];
            h_re[k + 1] = h_re[0] * yg_re - h_im[0] * yg_im;
            h_im[k + 1] = h_re[0] * yg_im + h_im[0] * yg_re;
        }

FREE_RETURN:
    nag_ode_ivp_rk_free(&opt);

    return err;
}

static Err HestonCalcPhiClosedForm(
    SHestonComm_InvFT* q, double u_re, double u_im, double* h_re, double* h_im)
{
    double A_re = 0.0, A_im = 0.0, C_re = 0.0, C_im = 0.0;
    double c1, c2, norm;

    c1 = q->rho * q->alpha * q->sigma;
    c2 = 0.5 * q->sigma * q->sigma;

    LGMSVUpdateADClosedForm(
        q->mat,
        q->lam,
        -0.5 * q->alpha * q->alpha,
        q->lam + c1 * u_im,
        -c1 * u_re,
        c2 * (u_re * u_re - u_im * u_im - u_im),
        c2 * (u_re + 2.0 * u_re * u_im),
        &C_re,
        &C_im,
        &A_re,
        &A_im);

    norm    = exp(A_re + C_re);
    h_re[0] = norm * cos(A_im + C_im);
    h_im[0] = norm * sin(A_im + C_im);

    return NULL;
}

static double TailIntegrand(double v, void* comm)
{
    SHestonComm_InvFT* q = (SHestonComm_InvFT*)comm;
    double             dd1, dd2, vk = v * q->k[q->want_k];

    dd1 = 1.0 + v * v;
    dd2 = dd1 * v;
    return cos(vk) / dd1 + sin(vk) / dd2;
}

static double IntegrateTail(SHestonComm_InvFT* q)
{
    const double tiny = 1.0e-307, halfpi = 1.5707963267949, tol = 1e-6;  // adjustable
    const int    lenaw = 8000;
    double       aw[8000], z, int_err, freq = fabs(q->k[q->want_k]);

    if (freq < 1e-16)
        return halfpi - atan(q->v_max);

    intdeoini(lenaw, tiny, tol, aw);
    intdeo(TailIntegrand, q->v_max, freq, aw, &z, &int_err, q);
    if (int_err < 0.0)
    {
        smessage("Tail integration failed");
        return log(-1.0);
    }

    return z;
}

static Err HestonInitCache(SHestonCache* cache, int maxpts, int nk)
{
    cache->maxpts = maxpts;
    cache->v      = (double*)calloc(maxpts, sizeof(double));
    cache->fn     = f3tensor(0, maxpts - 1, 0, nk - 1, 0, 4);
    if (!cache->v || !cache->fn)
        return serror("Memory failure");
    memset(cache->v, 0, maxpts * sizeof(double));
    return NULL;
}

static Err HestonFreeCache(SHestonCache* cache, int nk)
{
    if (cache->next)
        HestonFreeCache(cache->next, nk);
    free(cache->next);
    free(cache->v);
    if (cache->fn)
        free_f3tensor(cache->fn, 0, cache->maxpts - 1, 0, nk - 1, 0, 4);
    return NULL;
}

// Calculation of cos(vk)Re[Zeta(v)] + sin(vk)Im[Zeta(v)]

static double NAG_CALL Z_func(double v, Nag_User* comm)
{
    Err                err = NULL;
    SHestonComm_InvFT* q   = (SHestonComm_InvFT*)comm->p;
    double phi_re[5], phi_im[5], zeta_re[5], zeta_im[5];  //, phi_re_test[1], phi_im_test[1];
    double dd1, dd2;
    int    i, j;

    // Retrieve integrand value from cache or calculate it:

    if (fabs(v - q->c_head->v[0]) < 1e-16)
    {
        q->c_cur         = q->c_head;
        q->c_cur->cur_pt = 0;  // reset the counter if restarted from 0
    }
    else if (fabs(v - q->c_cur->v[q->c_cur->cur_pt]) > 1e-16)  // not yet calculated or no match
    {
        if (q->c_cur->v[q->c_cur->cur_pt] != 0.0)  // no match (branch point)
        {
            while (q->c_cur->next && fabs(v - q->c_cur->next->v[0]) > 1e-16)
                q->c_cur = q->c_cur->next;  // search in existing branches

            if (!q->c_cur->next)  // if not found - create a new branch
            {
                q->c_cur->next = (SHestonCache*)malloc(sizeof(SHestonCache));
                if (!q->c_cur->next)
                {
                    smessage("Memory failure");
                    return log(-1.0);
                }
                memset(q->c_cur->next, 0, sizeof(SHestonCache));
                err = HestonInitCache(q->c_cur->next, q->c_cur->maxpts, q->nk);
                if (err)
                {
                    smessage(err);
                    return log(-1.0);
                }
            }
            q->c_cur         = q->c_cur->next;
            q->c_cur->cur_pt = 0;
        }

        if (q->c_cur->v[q->c_cur->cur_pt] == 0.0)  // point not yet calculated -> calculate it
        {
            for (j = 0; j < 4 && !q->want[j]; j++)
                ;
            if (j < 4)
                err = HestonCalcPhi(q, v, -1.0, phi_re, phi_im);
            else
                err = HestonCalcPhiClosedForm(q, v, -1.0, phi_re, phi_im);
            if (err)
            {
                smessage(err);
                return log(-1.0);
            }            // return NaN if error
            q->count++;  // only increase the function call counter if not retrieving from cache
            q->c_cur->v[q->c_cur->cur_pt] = v;

            // calculate zeta for all functions:
            dd1 = 1.0 + v * v;
            dd2 = dd1 * v;

            for (j = 0; j < 5; j++)
                if (j == 0 || q->want[j - 1])
                {
                    zeta_re[j] = ((j == 0) - phi_re[j]) / dd1 + phi_im[j] / dd2;
                    zeta_im[j] = -phi_im[j] / dd1 + ((j == 0) - phi_re[j]) / dd2;
                }

            // calculate integrand for all strikes and all functions:
            for (i = 0; i < q->nk; i++)
            {
                dd1 = cos(v * q->k[i]);
                dd2 = sin(v * q->k[i]);

                for (j = 0; j < 5; j++)
                    if (j == 0 || q->want[j - 1])
                        q->c_cur->fn[q->c_cur->cur_pt][i][j] = dd1 * zeta_re[j] + dd2 * zeta_im[j];
            }
        }
    }

    if (++q->c_cur->cur_pt >= q->c_cur->maxpts)
    {
        smessage("maxpts exceeded in Z_func. Contact FIRST");
        return log(-1.0);
    }

    return q->c_cur->fn[q->c_cur->cur_pt - 1][q->want_k][q->want_fn];
}

static Err HestonDoIntegration(SHestonComm_InvFT* q, int idx_k, int idx_fn, double* res)
{
    Err              err    = NULL;
    const double     epsabs = 1e-7, epsrel = 1e-4;  // adjustable
    double           int_err;
    Nag_User         comm_Nag;
    NagError         fail;
    Nag_QuadProgress qp;

    memset(&fail, 0, sizeof(NagError));
    memset(&qp, 0, sizeof(Nag_QuadProgress));
    comm_Nag.p = q;
    q->want_k  = idx_k;
    q->want_fn = idx_fn;

    nag_1d_quad_gen_1(
        Z_func, 0.0, q->v_max, epsabs, epsrel, 200, res, &int_err, &qp, &comm_Nag, &fail);
    if (fail.code != NE_NOERROR)
    {
        err = serror(fail.message);
        goto FREE_RETURN;
    }

    if (idx_fn == 0)
        *res += IntegrateTail(q);

FREE_RETURN:
    NAG_FREE(qp.sub_int_beg_pts);
    NAG_FREE(qp.sub_int_end_pts);
    NAG_FREE(qp.sub_int_result);
    NAG_FREE(qp.sub_int_error);
    return err;
}

// Function called by NAG while calculating moments of Xt:

static void NAG_CALL
EvalDerivativesU(Integer neq, double t, double y[], double yp[], Nag_User* comm)
{
    SHestonComm_ODE* p = (SHestonComm_ODE*)comm->p;

    // dy1/dt, dy3/dt

    yp[0] = -p->lam * y[1];
    yp[2] = -p->lam * y[3];

    // dy2/dt, dy4/dt

    yp[1] = p->lam * y[1] + 0.5 * p->sigma * p->sigma;
    yp[3] = p->alpha * p->alpha * y[1] * y[1] + p->lam * y[3] +
            2.0 * p->rho * p->sigma * p->alpha * y[1] + p->sigma * p->sigma;
}

// Calculate mean and std of log(F_T/F_0)

static Err HestonCalcMoments(SHestonComm_InvFT* q, double* mean, double* std)
{
    Err             err = NULL;
    int             j;
    SHestonComm_ODE comm_RK;
    Nag_User        comm_Nag;
    NagError        fail;
    Nag_ODE_RK      opt;
    double          y[4], yp[4], ymax[4], RKthres[4], tgot;
    const double    RKtol = 1e-5, thres_min = 1e-8, thres_coef = 1e-7;  // adjustable

    memset(&fail, 0, sizeof(NagError));
    memset(&opt, 0, sizeof(Nag_ODE_RK));

    // Preinitialize comm structure for NAG

    comm_Nag.p   = &comm_RK;
    comm_RK.u_re = comm_RK.u_im = 0.0;
    comm_RK.sigma               = q->sigma;
    comm_RK.alpha               = q->alpha;
    comm_RK.rho                 = q->rho;
    comm_RK.lam                 = q->lam;

    memset(y, 0, 4 * sizeof(double));  // Final y values are all zeros.
    for (j = 0; j < 4; j++)
        RKthres[j] = thres_min;

    // Calculate solution at time 0

    nag_ode_ivp_rk_setup(
        4,
        q->mat,
        y,
        0.0,
        RKtol,
        RKthres,
        Nag_RK_4_5,
        Nag_RK_range,
        Nag_ErrorAssess_off,
        0.0,
        &opt,
        &fail);
    if (fail.code != NE_NOERROR)
    {
        err = serror(fail.message);
        goto FREE_RETURN;
    }

    nag_ode_ivp_rk_range(4, EvalDerivativesU, 0.0, &tgot, y, yp, ymax, &opt, &comm_Nag, &fail);
    if (fail.code != NE_NOERROR)
    {
        err = serror(fail.message);
        goto FREE_RETURN;
    }

    *mean = y[0] + y[1];
    *std  = sqrt(-y[2] - y[3]);

FREE_RETURN:
    nag_ode_ivp_rk_free(&opt);

    return err;
}

Err HestonDOptions(
    double   f0,
    double   sigma,
    double   alpha,
    double   rho,
    double   lam,
    double   shift,
    double   mat,
    int      nK,
    double*  K,
    char**   rec_pay_str,
    int*     want,
    double** res)
{
    Err               err = NULL;
    SHestonComm_InvFT comm_InvFT;
    SHestonCache      cache;
    const double      pi = 3.14159265358979, nstd = 20.0;  // adjustable
    const int         maxpts = 8000;
    SrtReceiverType   rec_pay;
    double            z, std;
    int               i, j;

    memset(&comm_InvFT, 0, sizeof(SHestonComm_InvFT));
    memset(&cache, 0, sizeof(SHestonCache));

    comm_InvFT.mat   = mat;
    comm_InvFT.sigma = sigma;
    comm_InvFT.alpha = alpha;
    comm_InvFT.rho   = rho;
    comm_InvFT.lam   = lam;
    memcpy(comm_InvFT.want, want, 4 * sizeof(int));

    // Calculate upper limit of integration

    err = HestonCalcMoments(&comm_InvFT, &z, &std);
    if (err)
        goto FREE_RETURN;

    comm_InvFT.v_max = nstd / std;

    comm_InvFT.nk = nK;
    comm_InvFT.k  = (double*)calloc(nK, sizeof(double));
    if (!comm_InvFT.k)
    {
        err = serror("Memory failure");
        goto FREE_RETURN;
    }

    err = HestonInitCache(&cache, maxpts, nK);
    if (err)
        goto FREE_RETURN;

    comm_InvFT.c_head = comm_InvFT.c_cur = &cache;

    for (i = 0; i < nK; i++)
        comm_InvFT.k[i] = log((K[i] + shift) / (f0 + shift));

    // Integrate for all strikes and all functions:

    for (i = 0; i < nK; i++)
    {
        err = interp_rec_pay(rec_pay_str[i], &rec_pay);
        if (err)
            goto FREE_RETURN;

        for (j = 0; j < 5; j++)
            if (j == 0 || want[j - 1])
            {
                err = HestonDoIntegration(&comm_InvFT, i, j, &z);
                if (err)
                    goto FREE_RETURN;

                res[i][j] = z / pi * (f0 + shift);
                if (j == 0 && rec_pay == SRT_PUT && K[i] > f0)
                    res[i][j] += K[i] - f0;
                if (j == 0 && rec_pay == SRT_CALL && K[i] < f0)
                    res[i][j] += f0 - K[i];
            }
    }

FREE_RETURN:
    free(comm_InvFT.k);
    HestonFreeCache(&cache, nK);

    return err;
}

Err NewtonD(
    double x,
    double x_min,
    double x_max,
    Err (*func)(double, double*, void*),
    double tol,
    int    maxiter,
    void*  comm)
{
    Err    err = NULL;
    double xx[3], ff[3];
    double a, b, c, delta, x1, x2;
    int    npts = 0, iter = 0;

    while (1)
    {
        if (iter++ > maxiter)
            return NULL;

        if (x < x_min)
            x = x_min;
        if (x > x_max)
            x = x_max;
        err = func(xx[0] = x, &ff[0], comm);
        if (err)
            return err;
        if (npts < 3)
            npts++;

        if (fabs(ff[0]) < tol)
            return NULL;

        if (npts == 1)
            x = xx[0] + 0.1 * ff[0];
        else if (fabs(ff[1] - ff[0]) < 1e-16)
            return NULL;
        else if (npts == 2)
            x = (xx[0] * ff[1] - xx[1] * ff[0]) / (ff[1] - ff[0]);
        else
        {
            a = ((b = (ff[2] - ff[1]) / (xx[2] - xx[1])) - (ff[1] - ff[0]) / (xx[1] - xx[0])) /
                (xx[2] - xx[0]);
            b -= a * (xx[1] + xx[2]);
            c = ff[0] - a * xx[0] * xx[0] - b * xx[0];

            delta = b * b - 4.0 * a * c;

            if (delta < 0.0)  // Switch back to linear interpolation
                x = (xx[0] * ff[1] - xx[1] * ff[0]) / (ff[1] - ff[0]);
            else
            {
                x1 = 0.5 * (-b + sqrt(delta)) / a;
                x2 = -b / a - x1;

                x = (fabs(xx[0] - x1) < fabs(xx[0] - x2) ? x1 : x2);
            }
        }

        xx[2] = xx[1];
        ff[2] = ff[1];
        xx[1] = xx[0];
        ff[1] = ff[0];
    }
    return serror("Should never get this message");
}

static Err EvalATMdiff(double x, double* diff, void* comm)
{
    Err                err   = NULL;
    SHestonComm_Calib* Calib = (SHestonComm_Calib*)comm;
    SHestonComm_InvFT  comm_InvFT;
    SHestonCache       cache;
    const double       pi = 3.14159265358979, nstd = 20.0;  // adjustable
    const int          maxpts = 8000;
    double             z, vol, std, zero = 0.0;

    memset(&comm_InvFT, 0, sizeof(SHestonComm_InvFT));
    memset(&cache, 0, sizeof(SHestonCache));

    comm_InvFT.mat   = Calib->mat;
    comm_InvFT.sigma = Calib->sigma = x;
    comm_InvFT.alpha                = Calib->alpha;
    comm_InvFT.rho                  = Calib->rho;
    comm_InvFT.lam                  = Calib->lam;
    memset(comm_InvFT.want, 0, 4 * sizeof(int));

    // Calculate upper limit of integration:

    err = HestonCalcMoments(&comm_InvFT, &z, &std);
    if (err)
        goto FREE_RETURN;

    comm_InvFT.v_max = nstd / std;

    comm_InvFT.nk = 1;
    comm_InvFT.k  = &zero;

    err = HestonInitCache(&cache, maxpts, 1);
    if (err)
        goto FREE_RETURN;

    comm_InvFT.c_head = comm_InvFT.c_cur = &cache;

    // Integrate:

    if (Calib->alpha / Calib->lam < 1e-3)
        z = srt_f_optblksch(
            Calib->fwd + Calib->shift,
            Calib->fwd + Calib->shift,
            Calib->sigma,
            Calib->mat,
            1.0,
            SRT_CALL,
            PREMIUM);
    else
    {
        err = HestonDoIntegration(&comm_InvFT, 0, 0, &z);
        if (err)
            goto FREE_RETURN;

        z *= (Calib->fwd + Calib->shift) / pi;
    }
    err =
        srt_f_optimpvol(z, Calib->fwd, Calib->fwd, Calib->mat, 1.0, SRT_CALL, SRT_LOGNORMAL, &vol);
    if (err)
        goto FREE_RETURN;

    *diff = vol - Calib->volATM;

    Calib->count++;

FREE_RETURN:
    HestonFreeCache(&cache, 1);

    return err;
}

static Err Get_RR_BF(SHestonComm_Calib* Calib, double* RR, double* BF)
{
    Err               err = NULL;
    int               i;
    SHestonComm_InvFT comm_InvFT;
    SHestonCache      cache;
    const double      pi = 3.14159265358979, nstd = 20.0, kmaxstd = 6.0;  // adjustable
    const int         maxpts = 8000;
    double            z, vols[2], std;

    memset(&comm_InvFT, 0, sizeof(SHestonComm_InvFT));
    memset(&cache, 0, sizeof(SHestonCache));

    comm_InvFT.mat   = Calib->mat;
    comm_InvFT.sigma = Calib->sigma;
    comm_InvFT.alpha = Calib->alpha;
    comm_InvFT.rho   = Calib->rho;
    comm_InvFT.lam   = Calib->lam;
    memset(comm_InvFT.want, 0, 4 * sizeof(int));

    // Calculate upper limit of integration:

    err = HestonCalcMoments(&comm_InvFT, &z, &std);
    if (err)
        goto FREE_RETURN;

    comm_InvFT.v_max = nstd / std;

    comm_InvFT.nk = 2;
    comm_InvFT.k  = Calib->k;

    err = HestonInitCache(&cache, maxpts, 2);
    if (err)
        goto FREE_RETURN;

    comm_InvFT.c_head = comm_InvFT.c_cur = &cache;

    // Integrate for all strikes:

    for (i = 0; i < 2; i++)
    {
        // Check that the strike is not too far from the money:

        if (fabs(Calib->k[i]) > kmaxstd * std)
        {
            err = serror("Calibration strike is too far from the money");
            goto FREE_RETURN;
        }

        if (Calib->alpha / Calib->lam < 1e-3)
            z = srt_f_optblksch(
                Calib->fwd + Calib->shift,
                Calib->K[i] + Calib->shift,
                Calib->sigma,
                Calib->mat,
                1.0,
                (Calib->K[i] < Calib->fwd ? SRT_PUT : SRT_CALL),
                PREMIUM);
        else
        {
            err = HestonDoIntegration(&comm_InvFT, i, 0, &z);
            if (err)
                goto FREE_RETURN;
            z *= (Calib->fwd + Calib->shift) / pi;
        }
        err = srt_f_optimpvol(
            z,
            Calib->fwd,
            Calib->K[i],
            Calib->mat,
            1.0,
            (Calib->K[i] < Calib->fwd ? SRT_PUT : SRT_CALL),
            SRT_LOGNORMAL,
            &vols[i]);
        if (err)
            goto FREE_RETURN;
    }

    *RR = vols[1] - vols[0];
    *BF = vols[1] + vols[0];

FREE_RETURN:
    HestonFreeCache(&cache, 2);
    return err;
}

static Err EvalBFdiff(double x, double* diff, void* comm)
{
    Err                err   = NULL;
    SHestonComm_Calib* Calib = (SHestonComm_Calib*)comm;

    Calib->alpha = sqrt(x);
    Calib->rho   = Calib->ar / (Calib->alpha + 1e-5);

    err = NewtonD(Calib->sigma, 1e-8, 2.0 / sqrt(Calib->mat), EvalATMdiff, 1e-4, 7, comm);
    if (err)
        return err;

    err = Get_RR_BF(Calib, &Calib->lastRR, &Calib->lastBF);
    if (err)
        return err;

    *diff = Calib->lastBF - Calib->BF;
    return NULL;
}

static Err EvalRRdiff(double x, double* diff, void* comm)
{
    Err                err   = NULL;
    SHestonComm_Calib* Calib = (SHestonComm_Calib*)comm;
    double             a_min, a_max;

    Calib->ar = x;
    a_min     = fabs(x) / 0.99;
    a_max     = 5.0 / sqrt(Calib->mat);
    if (a_max < a_min)
        a_max = a_min;
    if (Calib->alpha <= a_min)
        Calib->alpha = a_min + 0.05 * (a_max - a_min);

    err = NewtonD(
        Calib->alpha * Calib->alpha, a_min * a_min, a_max * a_max, EvalBFdiff, 1e-4, 7, comm);
    if (err)
        return err;

    *diff = Calib->lastRR - Calib->RR;
    return NULL;
}

Err HestonDCalibrate(
    double  f0,
    double  volf0,
    double  K1,
    double  vol1,
    double  K2,
    double  vol2,
    double  mat,
    double* sigma,
    double* alpha,
    double* rho,
    double  lam,
    double  shift,
    int     calib_smile)
{
    Err               err = NULL;
    SHestonComm_Calib Calib;
    int               i;
    double            ra_min, ra_max;

    memset(&Calib, 0, sizeof(SHestonComm_Calib));

    Calib.mat   = mat;
    Calib.fwd   = f0;
    Calib.sigma = volf0 * f0 / (f0 + shift);
    Calib.alpha = *alpha;
    Calib.rho   = *rho;
    Calib.lam   = lam;
    Calib.shift = shift;

    Calib.K[0] = K1;
    Calib.K[1] = K2;

    for (i = 0; i < 2; i++)
        if (Calib.K[i] > 0.0)
            Calib.k[i] = log((Calib.K[i] + shift) / (f0 + shift));

    Calib.volATM = volf0;
    Calib.RR     = vol2 - vol1;
    Calib.BF     = vol2 + vol1;

    Calib.count = 0;

    ra_max = 0.99 * 5.0 / sqrt(Calib.mat);
    ra_min = -ra_max;

    if (!calib_smile)
        err = NewtonD(Calib.sigma, 1e-8, 2.0 / sqrt(Calib.mat), EvalATMdiff, 1e-4, 7, &Calib);
    else
        err =
            NewtonD(Calib.rho * (Calib.alpha + 1e-5), ra_min, ra_max, EvalRRdiff, 1e-4, 7, &Calib);

    if (err)
        return err;

    // Return results:

    *sigma = Calib.sigma;
    *alpha = Calib.alpha;
    *rho   = Calib.rho;

    return NULL;
}

Err HestonDDensityFT(
    double  sigma,
    double  alpha,
    double  rho,
    double  lam,
    double  mat,
    double  u_re,
    double  u_im,
    double* h_re,
    double* h_im)
{
    double A_re = 0.0, A_im = 0.0, C_re = 0.0, C_im = 0.0;
    double c1, c2;

    c1 = rho * alpha * sigma;
    c2 = 0.5 * sigma * sigma;

    LGMSVUpdateADClosedForm(
        mat,
        lam,
        -0.5 * alpha * alpha,
        lam + c1 * u_im,
        -c1 * u_re,
        c2 * (u_re * u_re - u_im * u_im - u_im),
        c2 * (u_re + 2.0 * u_re * u_im),
        &C_re,
        &C_im,
        &A_re,
        &A_im);

    *h_re = A_re + C_re;
    *h_im = A_im + C_im;

    return NULL;
}
