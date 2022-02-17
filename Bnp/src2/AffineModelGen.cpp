// AffineModelGen.c : generic affine model option pricing

#include "math.h"
#include "srt_h_all.h"
#undef SIGN
#include "AffineModelGen.h"
#include "FxSVCalib.h"
#include "nag.h"
#include "nagd02.h"

static void FillAffineCoefFxSV(SAffinCoef *p, double t) {
  double Gamma_d, Gamma_f;
  SFxSVUndFull *o = p->o;
  SFxSVUndDesc *pmdl = o->pmdl;
  int i, idx_d = p->idx_d, idx_f = p->idx_f, idx_x = p->idx_x;
  double T_t = p->T_pay - t;

  Gamma_d = -o->sig_d[idx_d] * (1.0 - exp(-o->lam_d * T_t)) / o->lam_d;
  Gamma_f = -o->sig_f[idx_f] * (1.0 - exp(-o->lam_f * T_t)) / o->lam_f;

  // Drift:
  p->B[0][0] = -(Gamma_d * Gamma_d / 2.0 + Gamma_f * Gamma_f / 2.0 -
                 pmdl->rho[idx_x][0][1] * Gamma_d * Gamma_f);
  p->B[0][1] = 0.0;
  p->B[0][2] = -pmdl->beta[idx_x] * (pmdl->rho[idx_x][1][2] * Gamma_f -
                                     pmdl->rho[idx_x][0][2] * Gamma_d);
  p->B[0][3] = -0.5 * pmdl->beta[idx_x] * pmdl->beta[idx_x];

  p->B[1][0] = pmdl->rho[idx_x][0][3] * pmdl->alpha[idx_x] * Gamma_d +
               pmdl->gamma[idx_x];
  p->B[1][1] = 0.0;
  p->B[1][2] = -pmdl->gamma[idx_x];
  p->B[1][3] = 0.0;

  p->B[2][0] = pmdl->alpha[idx_x] * pmdl->alpha[idx_x];
  p->B[2][1] = 0.0;
  p->B[2][2] = 2.0 * p->B[1][0];
  p->B[2][3] = 2.0 * p->B[1][2];

  // Cross products:

  for (i = 0; i <= 3; i++)
    p->D[0][0][i] = -2.0 * p->B[0][i];

  p->D[0][1][0] = p->D[1][0][0] =
      pmdl->alpha[idx_x] *
      (pmdl->rho[idx_x][1][3] * Gamma_f - pmdl->rho[idx_x][0][3] * Gamma_d);
  p->D[0][1][1] = p->D[1][0][1] = 0.0;
  p->D[0][1][2] = p->D[1][0][2] =
      pmdl->alpha[idx_x] * pmdl->rho[idx_x][2][3] * pmdl->beta[idx_x];
  p->D[0][1][3] = p->D[1][0][3] = 0.0;

  p->D[0][2][0] = p->D[2][0][0] = 0.0;
  p->D[0][2][1] = p->D[2][0][1] = 0.0;
  p->D[0][2][2] = p->D[2][0][2] = 2.0 * p->D[0][1][0];
  p->D[0][2][3] = p->D[2][0][3] = 2.0 * p->D[0][1][2];

  p->D[1][1][0] = p->B[2][0];
  p->D[1][1][1] = p->D[1][1][2] = p->D[1][1][3] = 0.0;

  p->D[1][2][0] = p->D[2][1][0] = 0.0;
  p->D[1][2][1] = p->D[2][1][1] = 0.0;
  p->D[1][2][2] = p->D[2][1][2] = 2.0 * p->B[2][0];
  p->D[1][2][3] = p->D[2][1][3] = 0.0;

  p->D[2][2][0] = p->D[2][2][1] = p->D[2][2][2] = 0.0;
  p->D[2][2][3] = 4.0 * p->B[2][0];
}

static void FillAffineCoefHeston(SAffinCoef *p, double t) {
  SFxSVUndDesc *pmdl = p->o->pmdl;
  int idx = p->idx_x;

  // Drift:
  p->B[0][0] = 0.0;
  p->B[0][1] = 0.0;
  p->B[0][2] = -0.5 * pmdl->beta[idx] * pmdl->beta[idx];

  p->B[1][0] = pmdl->gamma[idx];
  p->B[1][1] = 0.0;
  p->B[1][2] = -pmdl->gamma[idx];

  // Cross products:
  p->D[0][0][0] = 0.0;
  p->D[0][0][1] = 0.0;
  p->D[0][0][2] = pmdl->beta[idx] * pmdl->beta[idx];

  p->D[0][1][0] = p->D[1][0][0] = 0.0;
  p->D[0][1][1] = p->D[1][0][1] = 0.0;
  p->D[0][1][2] = p->D[1][0][2] =
      pmdl->alpha[idx] * pmdl->rho[idx][2][3] * pmdl->beta[idx];

  p->D[1][1][0] = 0.0;
  p->D[1][1][1] = 0.0;
  p->D[1][1][2] = pmdl->alpha[idx] * pmdl->alpha[idx];
}

static void NAG_CALL EvalDerivatives(Integer neq, double t, double y[],
                                     double yp[], Nag_User *comm) {
  SAffinCoef *p = (SAffinCoef *)comm->p;
  int i, j, k;

  memset(yp, 0, neq * sizeof(double));
  //	FillAffineCoefFxSV(p  , t);
  FillAffineCoefHeston(p, t);

  // Drift:
  for (i = 1; i <= p->n; i++) {
    for (j = 0; j <= p->n; j++) {
      yp[2 * j] -= y[2 * i] * p->B[i - 1][j];
      yp[2 * j + 1] -= y[2 * i + 1] * p->B[i - 1][j];
    }
  }

  // Cross products:
  for (i = 1; i <= p->n; i++) {
    for (k = 0; k <= p->n; k++) {
      yp[2 * k] -= 0.5 * (y[2 * i] * y[2 * i] - y[2 * i + 1] * y[2 * i + 1]) *
                   p->D[i - 1][i - 1][k];
      yp[2 * k + 1] -= y[2 * i] * y[2 * i + 1] * p->D[i - 1][i - 1][k];
    }

    for (j = i + 1; j <= p->n; j++) {
      for (k = 0; k <= p->n; k++) {
        yp[2 * k] -= (y[2 * i] * y[2 * j] - y[2 * i + 1] * y[2 * j + 1]) *
                     p->D[i - 1][j - 1][k];
        yp[2 * k + 1] -= (y[2 * i] * y[2 * j + 1] + y[2 * i + 1] * y[2 * j]) *
                         p->D[i - 1][j - 1][k];
      }
    }
  }
}

#define NFACTORS 2 // 3
#define NEQ (NFACTORS * 2 + 2)

Err AffineCalcPhi(SFxSVComm_InvFT *q, double u_re, double u_im, double *h_re,
                  double *h_im) {
  Err err = NULL;
  int idx_d, idx_f, idx_x, i, j, neq;
  SFxSVUndFull *o = q->o;
  SAffinCoef comm_RK;
  Nag_User comm_Nag;
  NagError fail;
  Nag_ODE_RK opt;
  double y[NEQ], yp[NEQ], ymax[NEQ], RKthres[NEQ], tgot, norm;
  const double RKtol = 1e-5, thres_min = 1e-8, thres_coef = 1e-7; // adjustable

  memset(&fail, 0, sizeof(NagError));
  memset(&opt, 0, sizeof(Nag_ODE_RK));
  memset(&comm_RK, 0, sizeof(SAffinCoef));

  neq = NEQ;

  // Preinitialize comm structure for NAG

  comm_Nag.p = &comm_RK;
  comm_RK.o = o;
  comm_RK.T_pay = q->T_pay;
  comm_RK.n = NFACTORS;
  comm_RK.B = dmatrix(0, comm_RK.n - 1, 0, comm_RK.n);
  comm_RK.D = f3tensor(0, comm_RK.n - 1, 0, comm_RK.n - 1, 0, comm_RK.n);

  memset(y, 0, neq * sizeof(double)); // Final y values are almost all zeros.
  y[2] = -u_im;
  y[3] = u_re;

  for (j = 0; j < neq; j++)
    RKthres[j] = thres_min;

  idx_d = o->nsig_d - 1;
  idx_f = o->nsig_f - 1;
  idx_x = o->pmdl->ntimes - 1;

  // Proceed backwards integrating the system of ODE using Runge-Kutta

  for (i = o->ntimes - 2; i >= 0; i--) {
    while (idx_d > 0 && o->sigtms_d[idx_d - 1] > o->times[i] + 1e-5)
      idx_d--;
    while (idx_f > 0 && o->sigtms_f[idx_f - 1] > o->times[i] + 1e-5)
      idx_f--;
    while (idx_x > 0 && o->pmdl->times[idx_x - 1] > o->times[i] + 1e-5)
      idx_x--;

    // Fill in local constant coefficients in comm_Nag
    comm_RK.idx_d = idx_d;
    comm_RK.idx_f = idx_f;
    comm_RK.idx_x = idx_x;

    // Calculate solution at time i

    nag_ode_ivp_rk_setup(neq, o->times[i + 1], y, o->times[i], RKtol, RKthres,
                         Nag_RK_4_5, Nag_RK_range, Nag_ErrorAssess_off, 0.0,
                         &opt, &fail);
    if (fail.code != NE_NOERROR) {
      err = serror(fail.message);
      goto FREE_RETURN;
    }

    nag_ode_ivp_rk_range(neq, EvalDerivatives, o->times[i], &tgot, y, yp, ymax,
                         &opt, &comm_Nag, &fail);
    if (fail.code != NE_NOERROR) {
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
  //	norm = exp(y[0] + y[4] + y[6]);
  //	h_re[0] = norm * cos(y[1] + y[5] + y[7]);
  //	h_im[0] = norm * sin(y[1] + y[5] + y[7]);
  norm = exp(y[0] + y[4]);
  h_re[0] = norm * cos(y[1] + y[5]);
  h_im[0] = norm * sin(y[1] + y[5]);

FREE_RETURN:
  nag_ode_ivp_rk_free(&opt);
  if (comm_RK.B)
    free_dmatrix(comm_RK.B, 0, comm_RK.n - 1, 0, comm_RK.n);
  if (comm_RK.D)
    free_f3tensor(comm_RK.D, 0, comm_RK.n - 1, 0, comm_RK.n - 1, 0, comm_RK.n);

  return err;
}

static void NAG_CALL EvalDerivativesU(Integer neq, double t, double y[],
                                      double yp[], Nag_User *comm) {
  SAffinCoef *p = (SAffinCoef *)comm->p;
  int i, j, k;

  memset(yp, 0, neq * sizeof(double));
  //	FillAffineCoefFxSV(p  , t);
  FillAffineCoefHeston(p, t);

  // Drift:
  for (i = 1; i <= p->n; i++) {
    for (j = 0; j <= p->n; j++) {
      yp[2 * j] -= y[2 * i] * p->B[i - 1][j];
      yp[2 * j + 1] -= y[2 * i + 1] * p->B[i - 1][j];
    }
  }

  // Cross products:
  for (i = 1; i <= p->n; i++) {
    for (k = 0; k <= p->n; k++) {
      yp[2 * k] -= -y[2 * i + 1] * y[2 * i + 1] * p->D[i - 1][i - 1][k];
    }

    for (j = i + 1; j <= p->n; j++) {
      for (k = 0; k <= p->n; k++) {
        yp[2 * k] -= -2.0 * y[2 * i + 1] * y[2 * j + 1] * p->D[i - 1][j - 1][k];
      }
    }
  }
}

Err AffineCalcMoments(SFxSVComm_InvFT *q, double *mean, double *std) {
  Err err = NULL;
  int idx_d, idx_f, idx_x, i, j;
  SFxSVUndFull *o = q->o;
  SAffinCoef comm_RK;
  Nag_User comm_Nag;
  NagError fail;
  Nag_ODE_RK opt;
  double y[NEQ], yp[NEQ], ymax[NEQ], RKthres[NEQ], tgot;
  const double RKtol = 1e-5, thres_min = 1e-8, thres_coef = 1e-7; // adjustable

  memset(&fail, 0, sizeof(NagError));
  memset(&opt, 0, sizeof(Nag_ODE_RK));
  memset(&comm_RK, 0, sizeof(SAffinCoef));

  // Preinitialize comm structure for NAG

  comm_Nag.p = &comm_RK;
  comm_RK.o = o;
  comm_RK.T_pay = q->T_pay;
  comm_RK.n = NFACTORS;
  comm_RK.B = dmatrix(0, comm_RK.n - 1, 0, comm_RK.n);
  comm_RK.D = f3tensor(0, comm_RK.n - 1, 0, comm_RK.n - 1, 0, comm_RK.n);

  memset(y, 0, NEQ * sizeof(double)); // Final y values are almost all zeros.
  y[3] = 1.0;

  for (j = 0; j < NEQ; j++)
    RKthres[j] = thres_min;

  idx_d = o->nsig_d - 1;
  idx_f = o->nsig_f - 1;
  idx_x = o->pmdl->ntimes - 1;

  // Proceed backwards integrating the system of 6 ODE using Runge-Kutta

  for (i = o->ntimes - 2; i >= 0; i--) {
    while (idx_d > 0 && o->sigtms_d[idx_d - 1] > o->times[i] + 1e-5)
      idx_d--;
    while (idx_f > 0 && o->sigtms_f[idx_f - 1] > o->times[i] + 1e-5)
      idx_f--;
    while (idx_x > 0 && o->pmdl->times[idx_x - 1] > o->times[i] + 1e-5)
      idx_x--;

    // Fill in local constant coefficients in comm_Nag
    comm_RK.idx_d = idx_d;
    comm_RK.idx_f = idx_f;
    comm_RK.idx_x = idx_x;

    // Calculate solution at time i

    nag_ode_ivp_rk_setup(NEQ, o->times[i + 1], y, o->times[i], RKtol, RKthres,
                         Nag_RK_4_5, Nag_RK_range, Nag_ErrorAssess_off, 0.0,
                         &opt, &fail);
    if (fail.code != NE_NOERROR) {
      err = serror(fail.message);
      goto FREE_RETURN;
    }

    nag_ode_ivp_rk_range(NEQ, EvalDerivativesU, o->times[i], &tgot, y, yp, ymax,
                         &opt, &comm_Nag, &fail);
    if (fail.code != NE_NOERROR) {
      err = serror(fail.message);
      goto FREE_RETURN;
    }

    nag_ode_ivp_rk_free(&opt);

    for (j = 0; j < NEQ; j++) {
      RKthres[j] = thres_coef * ymax[j];
      if (RKthres[j] < thres_min)
        RKthres[j] = thres_min;
    }
  }
  //	*mean = y[1] + y[5] + y[7];
  //	*std = sqrt( -(y[0] + y[4] + y[6]) );
  *mean = y[1] + y[5];
  *std = sqrt(-(y[0] + y[4]));

FREE_RETURN:
  nag_ode_ivp_rk_free(&opt);

  return err;
}
