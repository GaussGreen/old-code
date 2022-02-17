// FxHestonD.cxx : FX Heston model

#include "srt_h_all.h"
#undef SIGN
#include "FxHestonD.h"
#include "LGMSVClosedForm.h"
#include "intde2.h"
#include "math.h"
#include "nag.h"
#include "nag_stdlib.h"
#include "nagd01.h"
#include "nagd02.h"
#include "opfnctns.h"

// Structure for communication with NAG Runge-Kutta integrator:
typedef struct _SFxHestonComm_ODE {
  double sigma, alpha, rho, lam;
  double u_re, u_im;
} SFxHestonComm_ODE;

// Cache to store calculated integrand values:
typedef struct _SFxHestonCache SFxHestonCache;
struct _SFxHestonCache {
  int maxpts;           // Max number of stored points in cache
  int cur_pt;           // Currently requested point
  double *v;            // Stored points
  double **fn;          // Stored integrand values: per v        , per k
  SFxHestonCache *next; // Linked list
};

// Structure for communication with a quadrature integrator
typedef struct _SFxHestonComm_InvFT {
  SFxHestonUnd *o;
  double mat;   // Options common maturity (t_fix)
  double v_max; // Upper limit of integration

  int nk;     // Number of strikes
  int want_k; // Current strike
  double *k;  // Strikes ( log(K'+sh/FFX0+sh) )	K'= K*df_for(T_pay
              //       ,T_star)/df_dom(T_pay        ,T_star)

  // Integrand cache:
  SFxHestonCache *c_head, *c_cur;

  // Counter for function calls:
  int count;

} SFxHestonComm_InvFT;

typedef struct _SFxHestonComm_Calib {
  SFxHestonUnd *und;
  double *k1, *K1, *k2, *K2; // Strikes ( log(K'+sh/FFX0+sh) and original)
  double *volATM, *RR, *BF;  // Market data
  double *f0,
      *pay_to_star; // forwards and T_star adjustments at each exercise date
  int cal_from, cal_to, cal_cur; // Pieces currently being calibrated
  double ar, lastRR, lastBF;     // temporary calibration variables
  int count;                     // Counter for function calls
} SFxHestonComm_Calib;

static void NAG_CALL EvalDerivatives(Integer neq, double t, double y[],
                                     double yp[], Nag_User *comm) {
  SFxHestonComm_ODE *p = (SFxHestonComm_ODE *)comm->p;

  double c1, c2, c3;

  // dy1/dt        , dy2/dt

  yp[0] = -p->lam * y[2];
  yp[1] = -p->lam * y[3];

  // dy3/dt        , dy4/dt

  c1 = p->alpha * p->alpha;
  c2 = p->rho * p->sigma * p->alpha;
  c3 = p->sigma * p->sigma;

  yp[2] = -0.5 * c1 * (y[2] * y[2] - y[3] * y[3]) + p->lam * y[2] +
          c2 * (p->u_re * y[3] + p->u_im * y[2]) +
          0.5 * c3 * (p->u_re * p->u_re - p->u_im * p->u_im - p->u_im);

  yp[3] = -c1 * y[2] * y[3] + p->lam * y[3] -
          c2 * (p->u_re * y[2] - p->u_im * y[3]) +
          0.5 * c3 * (p->u_re + 2.0 * p->u_re * p->u_im);
}

static Err FxHestonCalcPhi(SFxHestonComm_InvFT *q, double u_re, double u_im,
                           double *h_re, double *h_im) {
  Err err = NULL;
  int i, j, neq = 4;
  SFxHestonUnd *o = q->o;
  SFxHestonComm_ODE comm_RK;
  Nag_User comm_Nag;
  NagError fail;
  Nag_ODE_RK opt;
  double y[4], yp[4], ymax[4], RKthres[4], tgot, norm;
  const double RKtol = 1e-5, thres_min = 1e-8, thres_coef = 1e-7; // adjustable

  memset(&fail, 0, sizeof(NagError));
  memset(&opt, 0, sizeof(Nag_ODE_RK));

  // Preinitialize comm structure for NAG

  comm_Nag.p = &comm_RK;
  comm_RK.u_re = u_re;
  comm_RK.u_im = u_im;

  memset(y, 0, neq * sizeof(double)); // Final y values are all zeros.
  for (j = 0; j < neq; j++)
    RKthres[j] = thres_min;

  // Proceed backwards integrating the system of ODE using Runge-Kutta

  tgot = q->mat;
  for (i = 0; i < o->nt - 1 && o->t[i] < tgot - 1e-5; i++)
    ;

  for (; i >= 0; i--) {
    comm_RK.sigma = o->sigma[i];
    comm_RK.alpha = o->alpha[i];
    comm_RK.rho = o->rho[i];
    comm_RK.lam = o->lam[i];

    // Calculate solution at time i-1:

    nag_ode_ivp_rk_setup(neq, tgot, y, (i ? o->t[i - 1] : 0.0), RKtol, RKthres,
                         Nag_RK_4_5, Nag_RK_range, Nag_ErrorAssess_off, 0.0,
                         &opt, &fail);
    if (fail.cxxode != NE_NOERROR) {
      err = serror(fail.message);
      goto FREE_RETURN;
    }

    nag_ode_ivp_rk_range(neq, EvalDerivatives, (i ? o->t[i - 1] : 0.0), &tgot,
                         y, yp, ymax, &opt, &comm_Nag, &fail);
    if (fail.cxxode != NE_NOERROR) {
      err = serror(fail.message);
      goto FREE_RETURN;
    }

    nag_ode_ivp_rk_free(&opt);

    for (j = 0; j < neq; j++) {
      RKthres[j] = thres_coef * ymax[j];
      if (RKthres[j] < thres_min)
        RKthres[j] = thres_min;
    }
  }

  norm = exp(y[0] + y[2]);
  h_re[0] = norm * cos(y[1] + y[3]);
  h_im[0] = norm * sin(y[1] + y[3]);

FREE_RETURN:
  nag_ode_ivp_rk_free(&opt);

  return err;
}

static Err FxHestonCalcPhiClosedForm(SFxHestonComm_InvFT *q, double u_re,
                                     double u_im, double *h_re, double *h_im) {
  SFxHestonUnd *o = q->o;
  int i;
  double A_re = 0.0, A_im = 0.0, C_re = 0.0, C_im = 0.0;
  double t, c1, c2, norm;

  t = q->mat;
  for (i = 0; i < o->nt - 1 && o->t[i] < t - 1e-5; i++)
    ;

  for (; i >= 0; i--) {
    c1 = o->rho[i] * o->alpha[i] * o->sigma[i];
    c2 = 0.5 * o->sigma[i] * o->sigma[i];

    LGMSVUpdateADClosedForm(
        (i ? t - o->t[i - 1] : t), o->lam[i], -0.5 * o->alpha[i] * o->alpha[i],
        o->lam[i] + c1 * u_im, -c1 * u_re,
        c2 * (u_re * u_re - u_im * u_im - u_im),
        c2 * (u_re + 2.0 * u_re * u_im), &C_re, &C_im, &A_re, &A_im);

    if (i > 0)
      t = o->t[i - 1];
  }

  norm = exp(A_re + C_re);
  h_re[0] = norm * cos(A_im + C_im);
  h_im[0] = norm * sin(A_im + C_im);

  return NULL;
}

static double TailIntegrand(double v, void *comm) {
  SFxHestonComm_InvFT *q = (SFxHestonComm_InvFT *)comm;
  double dd1, dd2, vk = v * q->k[q->want_k];

  dd1 = 1.0 + v * v;
  dd2 = dd1 * v;
  return cos(vk) / dd1 + sin(vk) / dd2;
}

static double IntegrateTail(SFxHestonComm_InvFT *q) {
  const double tiny = 1.0e-307, halfpi = 1.5707963267949,
               tol = 1e-6; // adjustable
  const int lenaw = 8000;
  double aw[8000], z, int_err, freq = fabs(q->k[q->want_k]);

  if (freq < 1e-16)
    return halfpi - atan(q->v_max);

  intdeoini(lenaw, tiny, tol, aw);
  intdeo(TailIntegrand, q->v_max, freq, aw, &z, &int_err, q);
  if (int_err < 0.0) {
    smessage("Tail integration failed");
    return log(-1.0);
  }

  return z;
}

static Err FxHestonInitCache(SFxHestonCache *cache, int maxpts, int nk) {
  cache->maxpts = maxpts;
  cache->v = (double *)calloc(maxpts, sizeof(double));
  cache->fn = dmatrix(0, maxpts - 1, 0, nk - 1);
  if (!cache->v || !cache->fn)
    return serror("Memory failure");
  memset(cache->v, 0, maxpts * sizeof(double));
  return NULL;
}

static Err FxHestonFreeCache(SFxHestonCache *cache, int nk) {
  if (cache->next)
    FxHestonFreeCache(cache->next, nk);
  free(cache->next);
  free(cache->v);
  if (cache->fn)
    free_dmatrix(cache->fn, 0, cache->maxpts - 1, 0, nk - 1);
  return NULL;
}

// Calculation of cos(vk)Re[Zeta(v)] + sin(vk)Im[Zeta(v)]

static double NAG_CALL Z_func(double v, Nag_User *comm) {
  Err err = NULL;
  SFxHestonComm_InvFT *q = (SFxHestonComm_InvFT *)comm->p;
  double phi_re, phi_im, zeta_re,
      zeta_im; //        , phi_re_test        , phi_im_test;
  double dd1, dd2;
  int i;

  // Retrieve integrand value from cache or calculate it:

  if (fabs(v - q->c_head->v[0]) < 1e-16) {
    q->c_cur = q->c_head;
    q->c_cur->cur_pt = 0; // reset the counter if restarted from 0
  } else if (fabs(v - q->c_cur->v[q->c_cur->cur_pt]) >
             1e-16) // not yet calculated or no match
  {
    if (q->c_cur->v[q->c_cur->cur_pt] != 0.0) // no match (branch point)
    {
      while (q->c_cur->next && fabs(v - q->c_cur->next->v[0]) > 1e-16)
        q->c_cur = q->c_cur->next; // search in existing branches

      if (!q->c_cur->next) // if not found - create a new branch
      {
        q->c_cur->next = (SFxHestonCache *)malloc(sizeof(SFxHestonCache));
        if (!q->c_cur->next) {
          smessage("Memory failure");
          return log(-1.0);
        }
        memset(q->c_cur->next, 0, sizeof(SFxHestonCache));
        err = FxHestonInitCache(q->c_cur->next, q->c_cur->maxpts, q->nk);
        if (err) {
          smessage(err);
          return log(-1.0);
        }
      }
      q->c_cur = q->c_cur->next;
      q->c_cur->cur_pt = 0;
    }

    if (q->c_cur->v[q->c_cur->cur_pt] ==
        0.0) // point not yet calculated -> calculate it
    {
      //			err = FxHestonCalcPhi(q        , v        , -1.0
      //, &phi_re
      //      , &phi_im);
      err = FxHestonCalcPhiClosedForm(q, v, -1.0, &phi_re, &phi_im);
      if (err) {
        smessage(err);
        return log(-1.0);
      }           // return NaN if error
      q->count++; // only increase the function call counter if not retrieving
                  // from cache
      q->c_cur->v[q->c_cur->cur_pt] = v;

      // calculate zeta:
      dd1 = 1.0 + v * v;
      dd2 = dd1 * v;

      zeta_re = (1.0 - phi_re) / dd1 + phi_im / dd2;
      zeta_im = -phi_im / dd1 + (1.0 - phi_re) / dd2;

      // calculate integrand for all strikes:
      for (i = 0; i < q->nk; i++) {
        dd1 = cos(v * q->k[i]);
        dd2 = sin(v * q->k[i]);

        q->c_cur->fn[q->c_cur->cur_pt][i] = dd1 * zeta_re + dd2 * zeta_im;
      }
    }
  }

  if (++q->c_cur->cur_pt >= q->c_cur->maxpts) {
    smessage("maxpts exceeded in Z_func. Contact FIRST");
    return log(-1.0);
  }

  return q->c_cur->fn[q->c_cur->cur_pt - 1][q->want_k];
}

static Err FxHestonDoIntegration(SFxHestonComm_InvFT *q, int idx_k,
                                 double *res) {
  Err err = NULL;
  const double epsabs = 1e-7, epsrel = 1e-4; // adjustable
  double int_err;
  Nag_User comm_Nag;
  NagError fail;
  Nag_QuadProgress qp;

  memset(&fail, 0, sizeof(NagError));
  memset(&qp, 0, sizeof(Nag_QuadProgress));
  comm_Nag.p = q;
  q->want_k = idx_k;

  nag_1d_quad_gen_1(Z_func, 0.0, q->v_max, epsabs, epsrel, 200, res, &int_err,
                    &qp, &comm_Nag, &fail);
  if (fail.cxxode != NE_NOERROR) {
    err = serror(fail.message);
    goto FREE_RETURN;
  }

  *res += IntegrateTail(q);

FREE_RETURN:
  NAG_FREE(qp.sub_int_beg_pts);
  NAG_FREE(qp.sub_int_end_pts);
  NAG_FREE(qp.sub_int_result);
  NAG_FREE(qp.sub_int_error);
  return err;
}

// Function called by NAG while calculating moments of Xt:

static void NAG_CALL EvalDerivativesU(Integer neq, double t, double y[],
                                      double yp[], Nag_User *comm) {
  SFxHestonComm_ODE *p = (SFxHestonComm_ODE *)comm->p;

  // dy1/dt        , dy3/dt

  yp[0] = -p->lam * y[1];
  yp[2] = -p->lam * y[3];

  // dy2/dt        , dy4/dt

  yp[1] = p->lam * y[1] + 0.5 * p->sigma * p->sigma;
  yp[3] = p->alpha * p->alpha * y[1] * y[1] + p->lam * y[3] +
          2.0 * p->rho * p->sigma * p->alpha * y[1] + p->sigma * p->sigma;
}

// Calculate mean and std of log(F_T/F_0)

static Err FxHestonCalcMoments(SFxHestonComm_InvFT *q, double *mean,
                               double *std) {
  Err err = NULL;
  SFxHestonUnd *o = q->o;
  int i, j;
  SFxHestonComm_ODE comm_RK;
  Nag_User comm_Nag;
  NagError fail;
  Nag_ODE_RK opt;
  double y[4], yp[4], ymax[4], RKthres[4], tgot;
  const double RKtol = 1e-5, thres_min = 1e-8, thres_coef = 1e-7; // adjustable

  memset(&fail, 0, sizeof(NagError));
  memset(&opt, 0, sizeof(Nag_ODE_RK));

  // Preinitialize comm structure for NAG

  comm_Nag.p = &comm_RK;
  comm_RK.u_re = comm_RK.u_im = 0.0;

  memset(y, 0, 4 * sizeof(double)); // Final y values are all zeros.
  for (j = 0; j < 4; j++)
    RKthres[j] = thres_min;

  // Proceed backwards integrating the system of ODE using Runge-Kutta

  tgot = q->mat;
  for (i = 0; i < o->nt - 1 && o->t[i] < tgot - 1e-5; i++)
    ;

  for (; i >= 0; i--) {
    comm_RK.sigma = o->sigma[i];
    comm_RK.alpha = o->alpha[i];
    comm_RK.rho = o->rho[i];
    comm_RK.lam = o->lam[i];

    // Calculate solution at time i-1:

    nag_ode_ivp_rk_setup(4, tgot, y, (i ? o->t[i - 1] : 0.0), RKtol, RKthres,
                         Nag_RK_4_5, Nag_RK_range, Nag_ErrorAssess_off, 0.0,
                         &opt, &fail);
    if (fail.cxxode != NE_NOERROR) {
      err = serror(fail.message);
      goto FREE_RETURN;
    }

    nag_ode_ivp_rk_range(4, EvalDerivativesU, (i ? o->t[i - 1] : 0.0), &tgot, y,
                         yp, ymax, &opt, &comm_Nag, &fail);
    if (fail.cxxode != NE_NOERROR) {
      err = serror(fail.message);
      goto FREE_RETURN;
    }

    nag_ode_ivp_rk_free(&opt);

    for (j = 0; j < 4; j++) {
      RKthres[j] = thres_coef * ymax[j];
      if (RKthres[j] < thres_min)
        RKthres[j] = thres_min;
    }
  }

  *mean = y[0] + y[1];
  *std = sqrt(-y[2] - y[3]);

FREE_RETURN:
  nag_ode_ivp_rk_free(&opt);

  return err;
}

Err FxHestonDOptions(SFxHestonUnd *und, long fix_date, long pay_date, int nK,
                     double *K, char **rec_pay_str, double *res) {
  Err err = NULL;
  SFxHestonComm_InvFT comm_InvFT;
  SFxHestonCache cache;
  double df_pay_dom, df_pay_for, df_star_dom, df_star_for, f0, pay_to_star;
  const double pi = 3.14159265358979, nstd = 20.0; // adjustable
  const int maxpts = 8000;
  SrtReceiverType rec_pay;
  double z, std;
  int i;

  memset(&comm_InvFT, 0, sizeof(SFxHestonComm_InvFT));
  memset(&cache, 0, sizeof(SFxHestonCache));

  comm_InvFT.o = und;
  comm_InvFT.mat = (fix_date - und->today) * YEARS_IN_DAY;

  // Calculate upper limit of integration

  err = FxHestonCalcMoments(&comm_InvFT, &z, &std);
  if (err)
    goto FREE_RETURN;

  comm_InvFT.v_max = nstd / std;

  comm_InvFT.nk = nK;
  comm_InvFT.k = (double *)calloc(nK, sizeof(double));
  if (!comm_InvFT.k) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }

  err = FxHestonInitCache(&cache, maxpts, nK);
  if (err)
    goto FREE_RETURN;

  comm_InvFT.cxx_head = comm_InvFT.cxx_cur = &cache;

  df_pay_dom = swp_f_df(und->today, pay_date, und->yc_dom);
  df_pay_for = swp_f_df(und->today, pay_date, und->yc_for);
  df_star_dom = swp_f_df(und->today, und->D_star, und->yc_dom);
  df_star_for = swp_f_df(und->today, und->D_star, und->yc_for);
  f0 = und->S0 * df_pay_for / df_pay_dom;
  pay_to_star = (df_star_for / df_pay_for) / (df_star_dom / df_pay_dom);

  for (i = 0; i < nK; i++)
    comm_InvFT.k[i] = log((K[i] * pay_to_star + und->shift) /
                          (f0 * pay_to_star + und->shift));

  // Integrate for all strikes

  for (i = 0; i < nK; i++) {
    err = interp_rec_pay(rec_pay_str[i], &rec_pay);
    if (err)
      goto FREE_RETURN;

    err = FxHestonDoIntegration(&comm_InvFT, i, &z);
    if (err)
      goto FREE_RETURN;

    res[i] = z / pi * (f0 + und->shift / pay_to_star);
    if (rec_pay == SRT_PUT && K[i] > f0)
      res[i] += K[i] - f0;
    if (rec_pay == SRT_CALL && K[i] < f0)
      res[i] += f0 - K[i];

    res[i] *= df_pay_dom;
  }

FREE_RETURN:
  free(comm_InvFT.k);
  FxHestonFreeCache(&cache, nK);

  return err;
}

char *FreeFxHestonUnd(void *ptr) {
  SrtUndPtr undptr = (SrtUndPtr)ptr;
  SFxHestonUnd *und = NULL;

  if (!ptr)
    return NULL;
  und = (SFxHestonUnd *)undptr->spec_desc;
  if (!und)
    goto FREE_RETURN;

  free(und->t);
  free(und->sigma);
  free(und->alpha);
  free(und->lam);
  free(und->rho);

FREE_RETURN:
  free(und);
  free(ptr);
  return NULL;
}

static Err EvalATMdiff(double x, double *diff, void *comm) {
  Err err = NULL;
  SFxHestonComm_Calib *Calib = (SFxHestonComm_Calib *)comm;
  SFxHestonUnd *und = Calib->und;
  SFxHestonComm_InvFT comm_InvFT;
  SFxHestonCache cache;
  const double pi = 3.14159265358979, nstd = 20.0; // adjustable
  const int maxpts = 8000;
  double z, vol, std, zero = 0.0;
  int idx = Calib->cal_cur;

  memset(&comm_InvFT, 0, sizeof(SFxHestonComm_InvFT));
  memset(&cache, 0, sizeof(SFxHestonCache));

  comm_InvFT.o = und;
  comm_InvFT.mat = und->t[idx];

  und->sigma[idx] = x;

  // Calculate upper limit of integration:

  err = FxHestonCalcMoments(&comm_InvFT, &z, &std);
  if (err)
    goto FREE_RETURN;

  comm_InvFT.v_max = nstd / std;

  comm_InvFT.nk = 1;
  comm_InvFT.k = &zero;

  err = FxHestonInitCache(&cache, maxpts, 1);
  if (err)
    goto FREE_RETURN;

  comm_InvFT.cxx_head = comm_InvFT.cxx_cur = &cache;

  // Integrate:

  err = FxHestonDoIntegration(&comm_InvFT, 0, &z);
  if (err)
    goto FREE_RETURN;

  z *= (Calib->f0[idx] + und->shift / Calib->pay_to_star[idx]) / pi;

  err = srt_f_optimpvol(z, Calib->f0[idx], Calib->f0[idx], und->t[idx], 1.0,
                        SRT_CALL, SRT_LOGNORMAL, &vol);
  if (err)
    goto FREE_RETURN;

  *diff = vol - Calib->volATM[idx];

  Calib->count++;

FREE_RETURN:
  FxHestonFreeCache(&cache, 1);

  return err;
}

static Err CalibSigma(SFxHestonComm_Calib *Calib) {
  Err err = NULL;
  SFxHestonUnd *und = Calib->und;
  int i;
  double sig_1st_guess;

  for (i = Calib->cal_from; i <= Calib->cal_to; i++) {
    Calib->cal_cur = i;

    if (i > 0)
      sig_1st_guess = und->sigma[i - 1];
    else
      sig_1st_guess = Calib->volATM[0] * Calib->f0[0] /
                      (Calib->f0[0] + und->shift / Calib->pay_to_star[0]);

    err = NewtonD(sig_1st_guess, 1e-8, 2.0 / sqrt(und->t[i]), EvalATMdiff, 1e-4,
                  7, Calib);
    if (err)
      return err;
  }
  return NULL;
}

static Err Get_RR_BF(SFxHestonComm_Calib *Calib, double *RR, double *BF) {
  Err err = NULL;
  SFxHestonUnd *und = Calib->und;
  SFxHestonComm_InvFT comm_InvFT;
  SFxHestonCache cache;
  const double pi = 3.14159265358979, nstd = 20.0, kmaxstd = 6.0; // adjustable
  const int maxpts = 8000;
  int i, idx = Calib->cal_to;
  double z, std, vols[2], kk[2] = {Calib->k1[idx], Calib->k2[idx]};
  double KK[2] = {Calib->K1[idx], Calib->K2[idx]};

  memset(&comm_InvFT, 0, sizeof(SFxHestonComm_InvFT));
  memset(&cache, 0, sizeof(SFxHestonCache));

  comm_InvFT.o = und;
  comm_InvFT.mat = und->t[idx];

  // Calculate upper limit of integration:

  err = FxHestonCalcMoments(&comm_InvFT, &z, &std);
  if (err)
    goto FREE_RETURN;

  comm_InvFT.v_max = nstd / std;

  comm_InvFT.nk = 2;
  comm_InvFT.k = kk;

  err = FxHestonInitCache(&cache, maxpts, 2);
  if (err)
    goto FREE_RETURN;

  comm_InvFT.cxx_head = comm_InvFT.cxx_cur = &cache;

  // Integrate for all strikes:

  for (i = 0; i < 2; i++) {
    // Check that the strike is not too far from the money:

    if (fabs(kk[i]) > kmaxstd * std) {
      err = serror("Calibration strike is too far from the money");
      goto FREE_RETURN;
    }

    err = FxHestonDoIntegration(&comm_InvFT, i, &z);
    if (err)
      goto FREE_RETURN;

    z *= (Calib->f0[idx] + und->shift / Calib->pay_to_star[idx]) / pi;

    err = srt_f_optimpvol(z, Calib->f0[idx], KK[i], und->t[idx], 1.0,
                          (KK[i] < Calib->f0[idx] ? SRT_PUT : SRT_CALL),
                          SRT_LOGNORMAL, &vols[i]);
    if (err)
      goto FREE_RETURN;
  }

  *RR = vols[1] - vols[0];
  *BF = vols[1] + vols[0];

FREE_RETURN:
  FxHestonFreeCache(&cache, 2);
  return err;
}

static Err EvalBFdiff(double x, double *diff, void *comm) {
  Err err = NULL;
  SFxHestonComm_Calib *Calib = (SFxHestonComm_Calib *)comm;
  SFxHestonUnd *und = Calib->und;
  int i;
  double alpha, rho;

  alpha = sqrt(x);
  rho = Calib->ar / (alpha + 1e-5);

  for (i = Calib->cal_from; i <= Calib->cal_to; i++) {
    und->alpha[i] = alpha;
    und->rho[i] = rho;
  }
  err = CalibSigma(Calib);
  if (err)
    return err;

  err = Get_RR_BF(Calib, &Calib->lastRR, &Calib->lastBF);
  if (err)
    return err;

  *diff = Calib->lastBF - Calib->BF[Calib->cal_to];
  return NULL;
}

static Err EvalRRdiff(double x, double *diff, void *comm) {
  Err err = NULL;
  SFxHestonComm_Calib *Calib = (SFxHestonComm_Calib *)comm;
  SFxHestonUnd *und = Calib->und;
  double a_min, a_max, alpha_1st_guess;
  int i;

  Calib->ar = x;
  a_min = fabs(x) / 0.99;
  a_max = 20.0 / sqrt(und->t[Calib->cal_to]); // alpha limit

  for (i = Calib->cal_from; i <= Calib->cal_to; i++) {
    if (a_min < und->lam[i] * 1e-1) // avoid BS - currently numerically unstable
      a_min = und->lam[i] * 1e-1;
  }

  if (a_max < a_min)
    a_max = a_min;

  if (Calib->cal_from > 0)
    alpha_1st_guess = und->alpha[Calib->cal_from - 1];
  if (Calib->cal_from == 0 || alpha_1st_guess <= a_min)
    alpha_1st_guess = a_min + 0.05 * (a_max - a_min);

  err = NewtonD(alpha_1st_guess * alpha_1st_guess, a_min * a_min, a_max * a_max,
                EvalBFdiff, 1e-4, 7, comm);
  if (err)
    return err;

  *diff = Calib->lastRR - Calib->RR[Calib->cal_to];
  return NULL;
}

Err FxHestonDCalibrate(char *yc_dom, char *yc_for, double spot_fx, long D_star,
                       double beta, int ndates, long *dates, double *lam,
                       double *volf0, double *K1, double *volK1, double *K2,
                       double *volK2, int *calib_smile,
                       // Output:
                       double *sigma, double *alpha, double *rho) {
  Err err = NULL;
  SFxHestonUnd und;
  SFxHestonComm_Calib Calib;
  SrtCurvePtr crv_dom;
  long today, spotdate, paydate;
  double df_pay_dom, df_pay_for, df_star_dom, df_star_for;
  int i;
  double ra_min, ra_max, t1, t2;

  memset(&und, 0, sizeof(SFxHestonUnd));
  memset(&Calib, 0, sizeof(SFxHestonComm_Calib));

  crv_dom = lookup_curve(yc_dom);
  today = get_today_from_curve(crv_dom);
  spotdate = add_unit(today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

  strcpy(und.yc_dom, yc_dom);
  strcpy(und.yc_for, yc_for);
  und.S0 = spot_fx * swp_f_df(today, spotdate, yc_dom) /
           swp_f_df(today, spotdate, yc_for);
  und.D_star = D_star;
  und.today = today;
  und.shift = und.S0 * swp_f_df(today, D_star, yc_for) /
              swp_f_df(today, D_star, yc_dom) * (1.0 - beta) / beta;
  und.nt = ndates;
  und.sigma = sigma;
  und.alpha = alpha;
  und.rho = rho;
  und.lam = lam;

  Calib.und = &und;
  Calib.K1 = K1;
  Calib.K2 = K2;
  Calib.volATM = volf0;

  und.t = calloc(ndates, sizeof(double));
  Calib.k1 = calloc(ndates, sizeof(double));
  Calib.k2 = calloc(ndates, sizeof(double));
  Calib.RR = calloc(ndates, sizeof(double));
  Calib.BF = calloc(ndates, sizeof(double));
  Calib.f0 = calloc(ndates, sizeof(double));
  Calib.pay_to_star = calloc(ndates, sizeof(double));

  if (!und.t || !Calib.k1 || !Calib.k2 || !Calib.RR || !Calib.BF || !Calib.f0 ||
      !Calib.pay_to_star) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }

  for (i = 0; i < ndates; i++) {
    und.t[i] = (dates[i] - today) * YEARS_IN_DAY;

    paydate = add_unit(dates[i], 2, SRT_BDAY, MODIFIED_SUCCEEDING);
    df_pay_dom = swp_f_df(today, paydate, yc_dom);
    df_pay_for = swp_f_df(today, paydate, yc_for);
    df_star_dom = swp_f_df(today, D_star, yc_dom);
    df_star_for = swp_f_df(today, D_star, yc_for);
    Calib.f0[i] = und.S0 * df_pay_for / df_pay_dom;
    Calib.pay_to_star[i] =
        (df_star_for / df_pay_for) / (df_star_dom / df_pay_dom);

    if (calib_smile && calib_smile[i]) {
      Calib.k1[i] = log((K1[i] * Calib.pay_to_star[i] + und.shift) /
                        (Calib.f0[i] * Calib.pay_to_star[i] + und.shift));
      Calib.k2[i] = log((K2[i] * Calib.pay_to_star[i] + und.shift) /
                        (Calib.f0[i] * Calib.pay_to_star[i] + und.shift));

      Calib.RR[i] = volK2[i] - volK1[i];
      Calib.BF[i] = volK2[i] + volK1[i];
    }
  }

  for (Calib.cxxal_from = 0; Calib.cxxal_from < ndates;
       Calib.cxxal_from = Calib.cxxal_to + 1) {
    for (Calib.cxxal_to = Calib.cxxal_from;
         Calib.cxxal_to < ndates - 1 &&
         (!calib_smile || !calib_smile[Calib.cxxal_to]);
         Calib.cxxal_to++)
      ;

    if (calib_smile && calib_smile[Calib.cxxal_to]) {
      ra_max = 0.99 * 20.0 / sqrt(und.t[Calib.cxxal_to]); // alpha limit
      ra_min = -ra_max;

      t1 = clock();
      err = NewtonD(0.0, ra_min, ra_max, EvalRRdiff, 1e-4, 7, &Calib);
      if (err)
        goto FREE_RETURN;
      t2 = clock();
      smessage("Smile calibration at step %d        , time in sec: %.2f",
               Calib.cxxal_to, (double)(t2 - t1) / CLOCKS_PER_SEC);
    } else {
      err = CalibSigma(&Calib);
      if (err)
        goto FREE_RETURN;
    }
  }

FREE_RETURN:
  free(und.t);
  free(Calib.k1);
  free(Calib.k2);
  free(Calib.RR);
  free(Calib.BF);
  free(Calib.f0);
  free(Calib.pay_to_star);

  return err;
}
