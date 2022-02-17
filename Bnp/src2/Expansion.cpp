// Expansion.c : Automatic expansion pricing for stoch vol models

#include "Expansion.h"
#include "srt_h_all.h"
#undef SIGN
#include "nag.h"
#include "nagd02.h"

#define Expansion_MAXNBR 10
#define Expansion_MAXPOW 10

typedef struct _OldSPolynomial {
  int n_mono;
  long *mono;
  double **coef;
} OldSPolynomial;

typedef struct _OldSstructZ {
  OldSPolynomial Z[2], volZ[Expansion_MAXNBR];
  int d_idx, nd;
} OldSstructZ;

typedef struct _OldSHashTree OldSHashTree;
struct _OldSHashTree {
  OldSHashTree *facpow[Expansion_MAXPOW];
  long idx;
};

typedef struct _OldSControl {
  OldSstructZ **pool;
  int *local_pool;
  long pool_size, pool_pos, h_y0;
  OldSHashTree hash_tree;
  int nfac, nbit, nbr, nstp, maxmono;
  double ***rho; // [Bi][Bj][stp]
  OldSPolynomial *d_ODE;
  int d_pos, d_size;
} OldSControl;

Err OldInitControl(OldSControl *ctrl, int nfac, int nbit, int nbr, int nstp,
                   int maxmono, int d_size);
void OldFreeControl(OldSControl *ctrl);

#define Expansion_MAXORDER 40
#define Expansion_ZEROPOS 20

typedef struct _SOrdSeries {
  double val[Expansion_MAXORDER];
  int min, max;
} SOrdSeries;

Err ord_mult_const(SOrdSeries *s, double c);
Err ord_sum(SOrdSeries *s, SOrdSeries *s1);
Err ord_product(SOrdSeries *s, SOrdSeries *s1);
Err ord_shift(SOrdSeries *s, int shift);

#define MAXFAC 10

static Err OldHash(OldSControl *ctrl, long imono, long *hmono) {
  OldSHashTree *p = &ctrl->hash_tree;
  int i, mask = (1 << ctrl->nbit) - 1;
  long j;

  for (j = imono; j != 0; j >>= ctrl->nbit) {
    i = (j & mask);
    if (!p->facpow[i]) {
      p->facpow[i] = (OldSHashTree *)malloc(sizeof(OldSHashTree));
      if (!p->facpow[i])
        return serror("Memory failure");
      memset(p->facpow[i], 0, sizeof(OldSHashTree));
      if (ctrl->pool_pos >= ctrl->pool_size)
        return serror("Pool size exceeded");
      p->facpow[i]->idx = ctrl->pool_pos++;
    }
    p = p->facpow[i];
  }
  *hmono = p->idx;

  return NULL;
}

static void OldFreeHashTree(OldSHashTree *p) {
  int i;
  for (i = 0; i < Expansion_MAXPOW; i++) {
    if (p->facpow[i]) {
      OldFreeHashTree(p->facpow[i]);
      free(p->facpow[i]);
    }
  }
}

static Err AddProdMonoPoly(OldSControl *ctrl, OldSPolynomial *poly, double c,
                           long mono_idx, OldSPolynomial *p) {
  Err err;
  int j, k;
  long imono, hmono;

  for (j = 0; j < p->n_mono; j++) {
    imono = mono_idx + p->mono[j];
    err = OldHash(ctrl, imono, &hmono);
    if (err)
      return err;

    if (ctrl->local_pool[hmono] < 0) {
      ctrl->local_pool[hmono] = poly->n_mono++;
      if (poly->n_mono > ctrl->maxmono)
        return serror("ctrl->maxmono exceeded");
      poly->mono[poly->n_mono - 1] = imono;
    }
    for (k = 0; k < ctrl->nstp; k++)
      poly->coef[ctrl->local_pool[hmono]][k] += c * p->coef[j][k];
  }
  return NULL;
}

static Err AddProdMonoPolyPoly(OldSControl *ctrl, OldSPolynomial *poly,
                               double c, double *pc, long mono_idx,
                               OldSPolynomial *p1, OldSPolynomial *p2) {
  Err err;
  int i, j, k;
  long imono, hmono;

  for (i = 0; i < p1->n_mono; i++) {
    for (j = 0; j < p2->n_mono; j++) {
      imono = mono_idx + p1->mono[i] + p2->mono[j];
      err = OldHash(ctrl, imono, &hmono);
      if (err)
        return err;

      if (ctrl->local_pool[hmono] < 0) {
        ctrl->local_pool[hmono] = poly->n_mono++;
        if (poly->n_mono > ctrl->maxmono)
          return serror("ctrl->maxmono exceeded");
        poly->mono[poly->n_mono - 1] = imono;
      }
      for (k = 0; k < ctrl->nstp; k++)
        poly->coef[ctrl->local_pool[hmono]][k] +=
            (pc ? c * pc[k] : c) * p1->coef[i][k] * p2->coef[j][k];
    }
  }
  return NULL;
}

static Err OldFillStructZ(OldSControl *ctrl, long iZ) {
  Err err = NULL;
  int mask = (1 << ctrl->nbit) - 1;
  int facpow[MAXFAC], nfac = 0;
  long ifac[MAXFAC], hfac[MAXFAC], imono, hmono, hZ;
  int i, j, ki, kj, di, k;
  OldSstructZ *pZ;
  OldSPolynomial *dZ;

  err = OldHash(ctrl, iZ, &hZ);
  if (err)
    return err;

  if ((pZ = ctrl->pool[hZ]) == NULL) {
    pZ = ctrl->pool[hZ] = (OldSstructZ *)malloc(sizeof(OldSstructZ));
    if (!pZ)
      return serror("Memory failure");

    memset(pZ, 0, sizeof(OldSstructZ));

    for (i = 0; i < ctrl->nfac; i++) {
      facpow[nfac] = (iZ >> (i * ctrl->nbit)) & mask;
      if (facpow[nfac] > 0) {
        ifac[nfac] = (1 << (i * ctrl->nbit));
        err = OldHash(ctrl, ifac[nfac], &hfac[nfac]);
        if (err)
          return err;
        nfac++;
      }
    }

    for (i = 0; i < 2; i++) {
      pZ->Z[i].mono = (long *)calloc(ctrl->maxmono, sizeof(long));
      pZ->Z[i].coef = dmatrix(0, ctrl->maxmono - 1, 0, ctrl->nstp - 1);

      if (!pZ->Z[i].mono || !pZ->Z[i].coef)
        return serror("Memory failure");

      memset(&pZ->Z[i].coef[0][0], 0,
             ctrl->maxmono * ctrl->nstp * sizeof(double));
    }

    // Calculate drift:
    memset(ctrl->local_pool, -1, ctrl->pool_size * sizeof(int));

    for (i = 0; i < nfac; i++) {
      err = AddProdMonoPoly(ctrl, &pZ->Z[0], facpow[i], iZ - ifac[i],
                            &ctrl->pool[hfac[i]]->Z[0]);
      if (err)
        return err;

      if (facpow[i] >= 2) {
        for (ki = 0; ki < ctrl->nbr; ki++) {
          err = AddProdMonoPolyPoly(
              ctrl, &pZ->Z[0], 0.5 * facpow[i] * (facpow[i] - 1), NULL,
              iZ - 2 * ifac[i], &ctrl->pool[hfac[i]]->volZ[ki],
              &ctrl->pool[hfac[i]]->volZ[ki]);
          if (err)
            return err;

          for (kj = ki + 1; kj < ctrl->nbr; kj++) {
            err = AddProdMonoPolyPoly(
                ctrl, &pZ->Z[0], facpow[i] * (facpow[i] - 1), ctrl->rho[ki][kj],
                iZ - 2 * ifac[i], &ctrl->pool[hfac[i]]->volZ[ki],
                &ctrl->pool[hfac[i]]->volZ[kj]);
            if (err)
              return err;
          }
        }
      }

      for (j = i + 1; j < nfac; j++) {
        for (ki = 0; ki < ctrl->nbr; ki++) {
          for (kj = 0; kj < ctrl->nbr; kj++) {
            err = AddProdMonoPolyPoly(ctrl, &pZ->Z[0], facpow[i] * facpow[j],
                                      ctrl->rho[ki][kj], iZ - ifac[i] - ifac[j],
                                      &ctrl->pool[hfac[i]]->volZ[ki],
                                      &ctrl->pool[hfac[j]]->volZ[kj]);
            if (err)
              return err;
          }
        }
      }
    }

    // Calculate vol:

    for (ki = 0; ki < ctrl->nbr; ki++) {
      for (i = 0; i < nfac; i++) {
        if (ctrl->pool[hfac[i]]->volZ[ki].n_mono > 0) {
          if (!pZ->volZ[ki].mono) {
            pZ->volZ[ki].mono = (long *)calloc(ctrl->maxmono, sizeof(long));
            pZ->volZ[ki].coef =
                dmatrix(0, ctrl->maxmono - 1, 0, ctrl->nstp - 1);

            if (!pZ->volZ[ki].mono || !pZ->volZ[ki].coef)
              return serror("Memory failure");

            memset(&pZ->volZ[ki].coef[0][0], 0,
                   ctrl->maxmono * ctrl->nstp * sizeof(double));
            memset(ctrl->local_pool, -1, ctrl->pool_size * sizeof(int));
          }

          err = AddProdMonoPoly(ctrl, &pZ->volZ[ki], facpow[i], iZ - ifac[i],
                                &ctrl->pool[hfac[i]]->volZ[ki]);
          if (err)
            return err;
        }
      }
    }

    // Calculate <Z  ,Y0>:
    memset(ctrl->local_pool, -1, ctrl->pool_size * sizeof(int));

    for (ki = 0; ki < ctrl->nbr; ki++) {
      for (kj = 0; kj < ctrl->nbr; kj++) {
        err = AddProdMonoPolyPoly(ctrl, &pZ->Z[1], 1.0, ctrl->rho[ki][kj], 0,
                                  &pZ->volZ[ki],
                                  &ctrl->pool[ctrl->h_y0]->volZ[kj]);
        if (err)
          return err;
      }
    }
  } // if ( (pZ = ctrl->pool[iZ]) == NULL )

  //=====================================================================
  // Proceed recursively until all the monomials in Z are initialized:

  for (i = 0; i < 2; i++) {
    for (j = 0; j < pZ->Z[i].n_mono; j++) {
      imono = pZ->Z[i].mono[j];
      if (imono == iZ && i == 1)
        return serror("System unsolvable: <Z  ,Y0> depends on Z");
      err = OldHash(ctrl, imono, &hmono);
      if (err)
        return err;

      if (imono != 0 && imono != iZ &&
          (!ctrl->pool[hmono] || ctrl->pool[hmono]->nd == 0)) {
        err = OldFillStructZ(ctrl, imono);
        if (err)
          return err;
      }
      di = (imono ? ctrl->pool[hmono]->nd : 1) + i;
      if (di > pZ->nd)
        pZ->nd = di;
    }
  }

  // Build the ODE for iZ:

  pZ->d_idx = ctrl->d_pos;
  dZ = ctrl->d_ODE + pZ->d_idx;
  ctrl->d_pos += pZ->nd;
  if (ctrl->d_pos > ctrl->d_size)
    return serror("d_size exceeded");

  for (di = 0; di < pZ->nd;
       di++) // allocate memory for equations for each power of u
  {
    dZ[di].mono = (long *)calloc(ctrl->maxmono, sizeof(long));
    dZ[di].coef = dmatrix(0, ctrl->maxmono - 1, 0, ctrl->nstp - 1);

    if (!dZ[di].mono || !dZ[di].coef)
      return serror("Memory failure");
  }

  for (i = 0; i < 2; i++) {
    for (j = 0; j < pZ->Z[i].n_mono; j++) {
      imono = pZ->Z[i].mono[j];
      err = OldHash(ctrl, imono, &hmono);
      if (err)
        return err;

      if (imono == 0) {
        dZ[i].mono[dZ[i].n_mono] = -1;
        for (k = 0; k < ctrl->nstp; k++)
          dZ[i].coef[dZ[i].n_mono][k] = pZ->Z[i].coef[j][k];
        dZ[i].n_mono++;
        if (dZ[i].n_mono >= ctrl->maxmono)
          return serror("ctrl->maxmono exceeded (2)");
      } else {
        for (di = 0; di < ctrl->pool[hmono]->nd; di++) {
          dZ[i + di].mono[dZ[i + di].n_mono] = ctrl->pool[hmono]->d_idx + di;
          for (k = 0; k < ctrl->nstp; k++)
            dZ[i + di].coef[dZ[i + di].n_mono][k] = pZ->Z[i].coef[j][k];
          dZ[i + di].n_mono++;
          if (dZ[i + di].n_mono >= ctrl->maxmono)
            return serror("ctrl->maxmono exceeded (2)");
        }
      }
    }
  }
  return NULL;
}

typedef struct _SCommODE {
  int istp;
  OldSControl *ctrl;
} SCommODE;

static void NAG_CALL OldEvalDerivatives(Integer neq, double t, double y[],
                                        double yp[], Nag_User *comm) {
  OldSControl *ctrl = ((SCommODE *)comm->p)->ctrl;
  int i, j, k = ((SCommODE *)comm->p)->istp;
  double c;

  for (i = 0; i < ctrl->d_pos; i++) {
    yp[i] = 0.0;
    for (j = 0; j < ctrl->d_ODE[i].n_mono; j++) {
      c = ctrl->d_ODE[i].coef[j][k];
      if (ctrl->d_ODE[i].mono[j] >= 0)
        c *= y[ctrl->d_ODE[i].mono[j]];
      yp[i] += c;
    }
  }
}

static Err OldIntegrateD(OldSControl *ctrl, double *times, double *y) {
  Err err = NULL;
  int i, j, neq = ctrl->d_pos;
  SCommODE comm_RK;
  Nag_User comm_Nag;
  NagError fail;
  Nag_ODE_RK opt;
  double *yp = NULL, *ymax = NULL, *RKthres = NULL, t, tgot;
  const double RKtol = 1e-5, thres_min = 1e-8, thres_coef = 1e-7; // adjustable

  memset(&fail, 0, sizeof(NagError));
  memset(&opt, 0, sizeof(Nag_ODE_RK));

  yp = (double *)calloc(neq, sizeof(double));
  ymax = (double *)calloc(neq, sizeof(double));
  RKthres = (double *)calloc(neq, sizeof(double));

  if (!yp || !ymax || !RKthres) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }

  comm_Nag.p = &comm_RK;
  comm_RK.ctrl = ctrl;

  memset(y, 0, neq * sizeof(double)); // Initial y values are all zeros
  for (j = 0; j < neq; j++)
    RKthres[j] = thres_min;

  for (i = 0, t = 0.0; i < ctrl->nstp; t = times[i], i++) {
    comm_RK.istp = i;

    // Calculate solution at time i

    nag_ode_ivp_rk_setup(neq, t, y, times[i], RKtol, RKthres, Nag_RK_4_5,
                         Nag_RK_range, Nag_ErrorAssess_off, 0.0, &opt, &fail);
    if (fail.code != NE_NOERROR) {
      err = serror(fail.message);
      goto FREE_RETURN;
    }

    nag_ode_ivp_rk_range(neq, OldEvalDerivatives, times[i], &tgot, y, yp, ymax,
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

FREE_RETURN:
  nag_ode_ivp_rk_free(&opt);

  free(yp);
  free(ymax);
  free(RKthres);

  return err;
}

Err OldInitControl(OldSControl *ctrl, int nfac, int nbit, int nbr, int nstp,
                   int maxmono, int d_size) {
  Err err;
  int i, j;
  OldSstructZ *pZ;
  long hZ;

  ctrl->nfac = nfac;
  ctrl->nbit = nbit;
  ctrl->nbr = nbr;
  ctrl->nstp = nstp;
  ctrl->maxmono = maxmono;

  ctrl->pool_size = maxmono;
  ctrl->pool_pos = 1;

  ctrl->pool = (OldSstructZ **)calloc(ctrl->pool_size, sizeof(OldSstructZ *));
  if (!ctrl->pool)
    return serror("Memory failure");
  memset(ctrl->pool, 0, ctrl->pool_size * sizeof(OldSstructZ *));

  ctrl->local_pool = (int *)calloc(ctrl->pool_size, sizeof(int));
  if (!ctrl->local_pool)
    return serror("Memory failure");

  ctrl->d_size = d_size;

  ctrl->d_ODE = (OldSPolynomial *)calloc(d_size, sizeof(OldSPolynomial));
  if (!ctrl->d_ODE)
    return serror("Memory failure");
  memset(ctrl->d_ODE, 0, d_size * sizeof(OldSPolynomial));

  ctrl->rho = f3tensor(0, nbr - 1, 0, nbr - 1, 0, nstp - 1);
  if (!ctrl->rho)
    return serror("Memory failure");

  // Allocate memory for drifts and vols of basic factors

  for (i = 0; i < nfac; i++) {
    err = OldHash(ctrl, (1 << (i * nbit)), &hZ);
    if (err)
      return err;
    pZ = ctrl->pool[hZ] = (OldSstructZ *)malloc(sizeof(OldSstructZ));
    if (!pZ)
      return serror("Memory failure");
    memset(pZ, 0, sizeof(OldSstructZ));

    for (j = 0; j < 2; j++) {
      pZ->Z[j].mono = (long *)calloc(maxmono, sizeof(long));
      pZ->Z[j].coef = dmatrix(0, maxmono - 1, 0, nstp - 1);
      if (!pZ->Z[j].mono || !pZ->Z[j].coef)
        return serror("Memory failure");

      memset(&pZ->Z[j].coef[0][0], 0, maxmono * nstp * sizeof(double));
    }

    for (j = 0; j < nbr; j++) {
      pZ->volZ[j].mono = (long *)calloc(maxmono, sizeof(long));
      pZ->volZ[j].coef = dmatrix(0, maxmono - 1, 0, nstp - 1);
      if (!pZ->volZ[j].mono || !pZ->volZ[j].coef)
        return serror("Memory failure");

      memset(&pZ->volZ[j].coef[0][0], 0, maxmono * nstp * sizeof(double));
    }
  }

  return NULL;
}

void OldFreeControl(OldSControl *ctrl) {
  long i;
  int j;

  if (ctrl->pool)
    for (i = 0; i < ctrl->pool_pos; i++)
      if (ctrl->pool[i]) {
        for (j = 0; j < 2; j++) {
          free(ctrl->pool[i]->Z[j].mono);
          if (ctrl->pool[i]->Z[j].coef)
            free_dmatrix(ctrl->pool[i]->Z[j].coef, 0, ctrl->maxmono - 1, 0,
                         ctrl->nstp - 1);
        }
        for (j = 0; j < ctrl->nbr; j++) {
          free(ctrl->pool[i]->volZ[j].mono);
          if (ctrl->pool[i]->volZ[j].coef)
            free_dmatrix(ctrl->pool[i]->volZ[j].coef, 0, ctrl->maxmono - 1, 0,
                         ctrl->nstp - 1);
        }
        free(ctrl->pool[i]);
      }
  free(ctrl->pool);
  free(ctrl->local_pool);

  if (ctrl->d_ODE)
    for (i = 0; i < ctrl->d_pos; i++) {
      free(ctrl->d_ODE[i].mono);
      if (ctrl->d_ODE[i].coef)
        free_dmatrix(ctrl->d_ODE[i].coef, 0, ctrl->maxmono - 1, 0,
                     ctrl->nstp - 1);
    }
  free(ctrl->d_ODE);

  if (ctrl->rho)
    free_f3tensor(ctrl->rho, 0, ctrl->nbr - 1, 0, ctrl->nbr - 1, 0,
                  ctrl->nstp - 1);

  OldFreeHashTree(&ctrl->hash_tree);
}

#undef X1
#undef X2
#define Y0 (1 << 0)
#define X1 (1 << 2)
#define X2 (1 << 4)
#define Y1 (1 << 6)
#define Y2 (1 << 8)
#define X1X1 (X1 + X1)
#define Y1Y1 (Y1 + Y1)
#define MAXNEQ 500

Err sabr_test_get_d(double sigma, double alpha, double rho, double t,
                    double **d, int nd) {
  Err err = NULL;
  OldSControl C;
  OldSstructZ *pZ;
  double y[MAXNEQ];
  int i, ki, kj;
  long hY0, hX1, hX2, hY1, hY2, hY1Y1, hZ;

  memset(&C, 0, sizeof(OldSControl));

  err = OldInitControl(&C, 5, 2, 2, 1, 100, MAXNEQ);
  if (err)
    goto FREE_RETURN;

  C.rho[0][0][0] = C.rho[1][1][0] = 1.0;
  C.rho[0][1][0] = C.rho[1][0][0] = rho;

  // Get basic factors hash indices:

  if ((err = OldHash(&C, Y0, &hY0)) || (err = OldHash(&C, X1, &hX1)) ||
      (err = OldHash(&C, X2, &hX2)) || (err = OldHash(&C, Y1, &hY1)) ||
      (err = OldHash(&C, Y2, &hY2)) || (err = OldHash(&C, Y1Y1, &hY1Y1)))
    goto FREE_RETURN;

  // Initialize basic factors diffusion:

  // Y0 ==================================
  C.h_y0 = hY0;
  pZ = C.pool[hY0];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = 0; // const
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = 0; // const
  pZ->volZ[0].coef[0][0] = sigma;

  // X1 ==================================
  pZ = C.pool[hX1];

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = 0; // const
  pZ->volZ[1].coef[0][0] = alpha;

  // X2 ==================================
  pZ = C.pool[hX2];

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = X1;
  pZ->volZ[1].coef[0][0] = alpha;

  // Y1 ==================================
  pZ = C.pool[hY1];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X1;
  pZ->Z[0].coef[0][0] = -sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = X1;
  pZ->volZ[0].coef[0][0] = sigma;

  // Y2 ==================================
  pZ = C.pool[hY2];

  pZ->Z[0].n_mono = 2;
  pZ->Z[0].mono[0] = X1X1;
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;
  pZ->Z[0].mono[1] = X2;
  pZ->Z[0].coef[1][0] = -sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = X2;
  pZ->volZ[0].coef[0][0] = sigma;

  // Calculate <Z  ,Y0> for all basic factors:

  for (i = 0; i < C.nfac; i++) {
    err = OldHash(&C, (1 << (i * C.nbit)), &hZ);
    if (err)
      goto FREE_RETURN;

    pZ = C.pool[hZ];
    memset(C.local_pool, -1, C.pool_size * sizeof(int));

    for (ki = 0; ki < C.nbr; ki++) {
      for (kj = 0; kj < C.nbr; kj++) {
        err = AddProdMonoPolyPoly(&C, &pZ->Z[1], 1.0, C.rho[ki][kj], 0,
                                  &pZ->volZ[ki], &C.pool[hY0]->volZ[kj]);
        if (err)
          goto FREE_RETURN;
      }
    }
  }

  // Initial factors initialized  , now build the system of ODE:

  err = OldFillStructZ(&C, Y1);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y2);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1Y1);
  if (err)
    goto FREE_RETURN;

  // Integrate the system:

  err = OldIntegrateD(&C, &t, y);
  if (err)
    goto FREE_RETURN;

  // Return results:

  memset(d[0], 0, nd * sizeof(double));
  memset(d[1], 0, nd * sizeof(double));

  if (C.pool[hY1]->nd > nd || C.pool[hY2]->nd > nd ||
      C.pool[hY1Y1]->nd + 1 > nd) {
    err = serror("nd too small");
    goto FREE_RETURN;
  }

  for (i = 0; i < C.pool[hY1]->nd; i++)
    d[0][i] += y[C.pool[hY1]->d_idx + i];
  for (i = 0; i < C.pool[hY2]->nd; i++)
    d[1][i] += y[C.pool[hY2]->d_idx + i];
  for (i = 0; i < C.pool[hY1Y1]->nd; i++)
    d[1][i + 1] += 0.5 * y[C.pool[hY1Y1]->d_idx + i];

FREE_RETURN:
  OldFreeControl(&C);

  return err;
}

Err heston_test_get_d(double sigma, double gamma, double alpha, double rho,
                      double t, double *d, int nd) {
  Err err = NULL;
  OldSControl C;
  OldSstructZ *pZ;
  double y[MAXNEQ];
  int i, ki, kj;
  long hY0, hX1, hX2, hY1, hY2, hY1Y1, hZ;

  memset(&C, 0, sizeof(OldSControl));

  err = OldInitControl(&C, 5, 2, 2, 1, 100, MAXNEQ);
  if (err)
    goto FREE_RETURN;

  C.rho[0][0][0] = C.rho[1][1][0] = 1.0;
  C.rho[0][1][0] = C.rho[1][0][0] = rho;

  // Get basic factors hash indices:

  if ((err = OldHash(&C, Y0, &hY0)) || (err = OldHash(&C, X1, &hX1)) ||
      (err = OldHash(&C, X2, &hX2)) || (err = OldHash(&C, Y1, &hY1)) ||
      (err = OldHash(&C, Y2, &hY2)) || (err = OldHash(&C, Y1Y1, &hY1Y1)))
    goto FREE_RETURN;

  // Initialize basic factors diffusion:

  // Y0 ==================================
  C.h_y0 = hY0;
  pZ = C.pool[hY0];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = 0; // const
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = 0; // const
  pZ->volZ[0].coef[0][0] = sigma;

  // X1 ==================================
  pZ = C.pool[hX1];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X1;
  pZ->Z[0].coef[0][0] = -gamma;

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = 0; // const
  pZ->volZ[1].coef[0][0] = alpha;

  // X2 ==================================
  pZ = C.pool[hX2];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X2;
  pZ->Z[0].coef[0][0] = -gamma;

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = X1;
  pZ->volZ[1].coef[0][0] = 0.5 * alpha;

  // Y1 ==================================
  pZ = C.pool[hY1];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X1;
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = X1;
  pZ->volZ[0].coef[0][0] = 0.5 * sigma;

  // Y2 ==================================
  pZ = C.pool[hY2];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X2;
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;

  pZ->volZ[0].n_mono = 2;
  pZ->volZ[0].mono[0] = X2;
  pZ->volZ[0].coef[0][0] = 0.5 * sigma;
  pZ->volZ[0].mono[1] = X1X1;
  pZ->volZ[0].coef[1][0] = -sigma / 8.0;

  // Calculate <Z  ,Y0> for all basic factors:

  for (i = 0; i < C.nfac; i++) {
    err = OldHash(&C, (1 << (i * C.nbit)), &hZ);
    if (err)
      goto FREE_RETURN;

    pZ = C.pool[hZ];
    memset(C.local_pool, -1, C.pool_size * sizeof(int));

    for (ki = 0; ki < C.nbr; ki++) {
      for (kj = 0; kj < C.nbr; kj++) {
        err = AddProdMonoPolyPoly(&C, &pZ->Z[1], 1.0, C.rho[ki][kj], 0,
                                  &pZ->volZ[ki], &C.pool[hY0]->volZ[kj]);
        if (err)
          goto FREE_RETURN;
      }
    }
  }

  // Initial factors initialized  , now build the system of ODE:

  err = OldFillStructZ(&C, Y1);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y2);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1Y1);
  if (err)
    goto FREE_RETURN;

  // Integrate the system:

  err = OldIntegrateD(&C, &t, y);
  if (err)
    goto FREE_RETURN;

  // Return results:

  memset(d, 0, nd * sizeof(double));

  if (C.pool[hY1]->nd > nd || C.pool[hY2]->nd > nd ||
      C.pool[hY1Y1]->nd + 1 > nd) {
    err = serror("nd too small");
    goto FREE_RETURN;
  }

  for (i = 0; i < C.pool[hY1]->nd; i++)
    d[i] += y[C.pool[hY1]->d_idx + i];
  for (i = 0; i < C.pool[hY2]->nd; i++)
    d[i] += y[C.pool[hY2]->d_idx + i];
  for (i = 0; i < C.pool[hY1Y1]->nd; i++)
    d[i + 1] += 0.5 * y[C.pool[hY1Y1]->d_idx + i];

FREE_RETURN:
  OldFreeControl(&C);

  return err;
}

#undef Y0
#undef X1
#undef X2
#undef Y1
#undef Y2
#undef X1X1
#undef Y1Y1

#define Y0 (1 << 0)
#define X1 (1 << 2)
#define X2 (1 << 4)
#define Y1V (1 << 6)
#define Y1D (1 << 8)
#define Y2V (1 << 10)
#define Y2D (1 << 12)
#define X1X1 (X1 + X1)
#define Y1VY1V (Y1V + Y1V)
#define Y1VY1D (Y1V + Y1D)
#define Y1DY1D (Y1D + Y1D)
/*
Err sabr_test_get_all_d(double sigma  , double alpha  , double rho  , double t
, double **d  , int nd)
{
        Err				err = NULL;
        OldSControl		C;
        OldSstructZ		*pZ;
        double			y[MAXNEQ];
        int				i  , ki  , kj;
        long			hY0  , hX1  , hX2  , hY1V  , hY1D  , hY2V  ,
hY2D  , hY1VY1V  , hY1VY1D  , hY1DY1D  , hZ;

        memset(&C  , 0  , sizeof(OldSControl));

        err = OldInitControl(&C  , 7  , 2  , 2  , 1  , 100  , MAXNEQ);
        if (err) goto FREE_RETURN;

        C.rho[0][0][0] = C.rho[1][1][0] = 1.0;
        C.rho[0][1][0] = C.rho[1][0][0] = rho;

        // Get basic factors hash indices:

        if ( (err = OldHash(&C  , Y0  , &hY0)) ||
                 (err = OldHash(&C  , X1  , &hX1)) ||
                 (err = OldHash(&C  , X2  , &hX2)) ||
                 (err = OldHash(&C  , Y1V  , &hY1V)) ||
                 (err = OldHash(&C  , Y1D  , &hY1D)) ||
                 (err = OldHash(&C  , Y2V  , &hY2V)) ||
                 (err = OldHash(&C  , Y2D  , &hY2D)) ||
                 (err = OldHash(&C  , Y1VY1V  , &hY1VY1V)) ||
                 (err = OldHash(&C  , Y1VY1D  , &hY1VY1D)) ||
                 (err = OldHash(&C  , Y1DY1D  , &hY1DY1D)) )
                goto FREE_RETURN;

        // Initialize basic factors diffusion:

        // Y0 ==================================
        C.h_y0 = hY0;
        pZ = C.pool[hY0];

        pZ->Z[0].n_mono = 1;
        pZ->Z[0].mono[0] = 0;			// const
        pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;

        pZ->volZ[0].n_mono = 1;
        pZ->volZ[0].mono[0] = 0;		// const
        pZ->volZ[0].coef[0][0] = sigma;

        // X1 ==================================
        pZ = C.pool[hX1];

        pZ->volZ[1].n_mono = 1;
        pZ->volZ[1].mono[0] = 0;		// const
        pZ->volZ[1].coef[0][0] = alpha;

        // X2 ==================================
        pZ = C.pool[hX2];

        pZ->volZ[1].n_mono = 1;
        pZ->volZ[1].mono[0] = X1;
        pZ->volZ[1].coef[0][0] = alpha;

        // Y1V ==================================
        pZ = C.pool[hY1V];

        pZ->volZ[0].n_mono = 1;
        pZ->volZ[0].mono[0] = X1;
        pZ->volZ[0].coef[0][0] = sigma;

        // Y1D ==================================
        pZ = C.pool[hY1D];

        pZ->Z[0].n_mono = 1;
        pZ->Z[0].mono[0] = X1;
        pZ->Z[0].coef[0][0] = -sigma * sigma;

        // Y2V ==================================
        pZ = C.pool[hY2V];

        pZ->volZ[0].n_mono = 1;
        pZ->volZ[0].mono[0] = X2;
        pZ->volZ[0].coef[0][0] = sigma;

        // Y2D ==================================
        pZ = C.pool[hY2D];

        pZ->Z[0].n_mono = 2;
        pZ->Z[0].mono[0] = X1X1;
        pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;
        pZ->Z[0].mono[1] = X2;
        pZ->Z[0].coef[1][0] = -sigma * sigma;


        // Calculate <Z  ,Y0> for all basic factors:

        for (i=0; i < C.nfac; i++)
        {
                err = OldHash(&C  , (1 << (i * C.nbit))  , &hZ);
                if (err) goto FREE_RETURN;

                pZ = C.pool[hZ];
                memset(C.local_pool  , -1  , C.pool_size * sizeof(int));

                for (ki=0; ki < C.nbr; ki++)
                {
                        for (kj=0; kj < C.nbr; kj++)
                        {
                                err = AddProdMonoPolyPoly( &C  , &pZ->Z[1]  ,
                                        1.0  , C.rho[ki][kj]  , 0  ,
&pZ->volZ[ki]  , &C.pool[hY0]->volZ[kj] ); if (err) goto FREE_RETURN;
                        }
                }
        }

        // Initial factors initialized  , now build the system of ODE:

        err = OldFillStructZ(&C  , Y1V);
        if (err) goto FREE_RETURN;

        err = OldFillStructZ(&C  , Y1D);
        if (err) goto FREE_RETURN;

        err = OldFillStructZ(&C  , Y2V);
        if (err) goto FREE_RETURN;

        err = OldFillStructZ(&C  , Y2D);
        if (err) goto FREE_RETURN;

        err = OldFillStructZ(&C  , Y1VY1V);
        if (err) goto FREE_RETURN;

        err = OldFillStructZ(&C  , Y1VY1D);
        if (err) goto FREE_RETURN;

        err = OldFillStructZ(&C  , Y1DY1D);
        if (err) goto FREE_RETURN;

        // Integrate the system:

        err = OldIntegrateD(&C  , &t  , y);
        if (err) goto FREE_RETURN;

        // Return results:

        for (i=0; i < 7; i++) memset(d[i]  , 0  , nd * sizeof(double));

        if (C.pool[hY1VY1V]->nd > nd || C.pool[hY1VY1D]->nd > nd ||
C.pool[hY1DY1D]->nd > nd) { err = serror("nd too small");  goto FREE_RETURN; }

        for (i=0; i < C.pool[hY1V]->nd; i++) d[0][i] = y[C.pool[hY1V]->d_idx +
i]; for (i=0; i < C.pool[hY1D]->nd; i++) d[1][i] = y[C.pool[hY1D]->d_idx + i];
        for (i=0; i < C.pool[hY2V]->nd; i++) d[2][i] = y[C.pool[hY2V]->d_idx +
i]; for (i=0; i < C.pool[hY2D]->nd; i++) d[3][i] = y[C.pool[hY2D]->d_idx + i];
        for (i=0; i < C.pool[hY1VY1V]->nd; i++) d[4][i] =
y[C.pool[hY1VY1V]->d_idx + i] * 0.5; for (i=0; i < C.pool[hY1VY1D]->nd; i++)
d[5][i] = y[C.pool[hY1VY1D]->d_idx + i]; for (i=0; i < C.pool[hY1DY1D]->nd; i++)
d[6][i] = y[C.pool[hY1DY1D]->d_idx + i] * 0.5;

FREE_RETURN:
        OldFreeControl(&C);

        return err;
}
*/
Err heston_test_get_all_d2(double sigma, double gamma, double alpha, double rho,
                           double t, double **d, int nd) {
  Err err = NULL;
  OldSControl C;
  OldSstructZ *pZ;
  double y[MAXNEQ];
  int i, ki, kj;
  long hY0, hX1, hX2, hY1V, hY1D, hY2V, hY2D, hY1VY1V, hY1VY1D, hY1DY1D, hZ;

  memset(&C, 0, sizeof(OldSControl));

  err = OldInitControl(&C, 7, 2, 2, 1, 100, MAXNEQ);
  if (err)
    goto FREE_RETURN;

  C.rho[0][0][0] = C.rho[1][1][0] = 1.0;
  C.rho[0][1][0] = C.rho[1][0][0] = rho;

  // Get basic factors hash indices:

  if ((err = OldHash(&C, Y0, &hY0)) || (err = OldHash(&C, X1, &hX1)) ||
      (err = OldHash(&C, X2, &hX2)) || (err = OldHash(&C, Y1V, &hY1V)) ||
      (err = OldHash(&C, Y1D, &hY1D)) || (err = OldHash(&C, Y2V, &hY2V)) ||
      (err = OldHash(&C, Y2D, &hY2D)) ||
      (err = OldHash(&C, Y1VY1V, &hY1VY1V)) ||
      (err = OldHash(&C, Y1VY1D, &hY1VY1D)) ||
      (err = OldHash(&C, Y1DY1D, &hY1DY1D)))
    goto FREE_RETURN;

  // Initialize basic factors diffusion:

  // Y0 ==================================
  C.h_y0 = hY0;
  pZ = C.pool[hY0];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = 0; // const
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = 0; // const
  pZ->volZ[0].coef[0][0] = sigma;

  // X1 ==================================
  pZ = C.pool[hX1];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X1;
  pZ->Z[0].coef[0][0] = -gamma;

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = 0; // const
  pZ->volZ[1].coef[0][0] = alpha;

  // X2 ==================================
  pZ = C.pool[hX2];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X2;
  pZ->Z[0].coef[0][0] = -gamma;

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = X1;
  pZ->volZ[1].coef[0][0] = 0.5 * alpha;

  // Y1V ==================================
  pZ = C.pool[hY1V];

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = X1;
  pZ->volZ[0].coef[0][0] = 0.5 * sigma;

  // Y1D ==================================
  pZ = C.pool[hY1D];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X1;
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;

  // Y2V ==================================
  pZ = C.pool[hY2V];

  pZ->volZ[0].n_mono = 2;
  pZ->volZ[0].mono[0] = X2;
  pZ->volZ[0].coef[0][0] = 0.5 * sigma;
  pZ->volZ[0].mono[1] = X1X1;
  pZ->volZ[0].coef[1][0] = -sigma / 8.0;

  // Y2D ==================================
  pZ = C.pool[hY2D];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X2;
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;

  // Calculate <Z  ,Y0> for all basic factors:

  for (i = 0; i < C.nfac; i++) {
    err = OldHash(&C, (1 << (i * C.nbit)), &hZ);
    if (err)
      goto FREE_RETURN;

    pZ = C.pool[hZ];
    memset(C.local_pool, -1, C.pool_size * sizeof(int));

    for (ki = 0; ki < C.nbr; ki++) {
      for (kj = 0; kj < C.nbr; kj++) {
        err = AddProdMonoPolyPoly(&C, &pZ->Z[1], 1.0, C.rho[ki][kj], 0,
                                  &pZ->volZ[ki], &C.pool[hY0]->volZ[kj]);
        if (err)
          goto FREE_RETURN;
      }
    }
  }

  // Initial factors initialized  , now build the system of ODE:

  err = OldFillStructZ(&C, Y1V);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1D);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y2V);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y2D);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1VY1V);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1VY1D);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1DY1D);
  if (err)
    goto FREE_RETURN;

  // Integrate the system:

  err = OldIntegrateD(&C, &t, y);
  if (err)
    goto FREE_RETURN;

  // Return results:

  for (i = 0; i < 7; i++)
    memset(d[i], 0, nd * sizeof(double));

  if (C.pool[hY1VY1V]->nd > nd || C.pool[hY1VY1D]->nd > nd ||
      C.pool[hY1DY1D]->nd > nd) {
    err = serror("nd too small");
    goto FREE_RETURN;
  }

  for (i = 0; i < C.pool[hY1V]->nd; i++)
    d[0][i] = y[C.pool[hY1V]->d_idx + i];
  for (i = 0; i < C.pool[hY1D]->nd; i++)
    d[1][i] = y[C.pool[hY1D]->d_idx + i];
  for (i = 0; i < C.pool[hY2V]->nd; i++)
    d[2][i] = y[C.pool[hY2V]->d_idx + i];
  for (i = 0; i < C.pool[hY2D]->nd; i++)
    d[3][i] = y[C.pool[hY2D]->d_idx + i];
  for (i = 0; i < C.pool[hY1VY1V]->nd; i++)
    d[4][i] = y[C.pool[hY1VY1V]->d_idx + i] * 0.5;
  for (i = 0; i < C.pool[hY1VY1D]->nd; i++)
    d[5][i] = y[C.pool[hY1VY1D]->d_idx + i];
  for (i = 0; i < C.pool[hY1DY1D]->nd; i++)
    d[6][i] = y[C.pool[hY1DY1D]->d_idx + i] * 0.5;

FREE_RETURN:
  OldFreeControl(&C);

  return err;
}

#undef Y0
#undef X1
#undef X2
#undef Y1
#undef Y2
#undef X1X1
#undef Y1Y1

#define Y0 (1 << 0)
#define X1 (1 << 3)
#define X2 (1 << 6)
#define Y1 (1 << 9)
#define Y2 (1 << 12)
#define X3 (1 << 15)
#define Y3 (1 << 18)
#define X4 (1 << 21)
#define Y4 (1 << 24)
#define X1X1 (X1 + X1)
#define Y1Y1 (Y1 + Y1)
#define X1X2 (X1 + X2)
#define Y1Y2 (Y1 + Y2)
#define Y1Y1Y1 (Y1 + Y1 + Y1)
#define Y2Y2 (Y2 + Y2)
#define X1X3 (X1 + X3)
#define X2X2 (X2 + X2)
#define Y1Y3 (Y1 + Y3)
#define Y1Y1Y2 (Y1 + Y1 + Y2)
#define Y1Y1Y1Y1 (Y1 + Y1 + Y1 + Y1)

Err sabr_test4_get_d(double sigma, double alpha, double rho, double t,
                     double **d, int nd) {
  Err err = NULL;
  OldSControl C;
  OldSstructZ *pZ;
  double y[MAXNEQ];
  int i, ki, kj;
  long hY0, hX1, hX2, hY1, hY2, hY1Y1, hX3, hY3, hY1Y2, hY1Y1Y1;
  long hX4, hY4, hY2Y2, hY1Y3, hY1Y1Y2, hY1Y1Y1Y1, hZ;

  memset(&C, 0, sizeof(OldSControl));

  err = OldInitControl(&C, 9, 3, 2, 1, 100, MAXNEQ);
  if (err)
    goto FREE_RETURN;

  C.rho[0][0][0] = C.rho[1][1][0] = 1.0;
  C.rho[0][1][0] = C.rho[1][0][0] = rho;

  // Get basic factors hash indices:

  if ((err = OldHash(&C, Y0, &hY0)) || (err = OldHash(&C, X1, &hX1)) ||
      (err = OldHash(&C, X2, &hX2)) || (err = OldHash(&C, Y1, &hY1)) ||
      (err = OldHash(&C, Y2, &hY2)) || (err = OldHash(&C, Y1Y1, &hY1Y1)) ||
      (err = OldHash(&C, X3, &hX3)) || (err = OldHash(&C, Y3, &hY3)) ||
      (err = OldHash(&C, Y1Y2, &hY1Y2)) ||
      (err = OldHash(&C, Y1Y1Y1, &hY1Y1Y1)) || (err = OldHash(&C, X4, &hX4)) ||
      (err = OldHash(&C, Y4, &hY4)) || (err = OldHash(&C, Y2Y2, &hY2Y2)) ||
      (err = OldHash(&C, Y1Y3, &hY1Y3)) ||
      (err = OldHash(&C, Y1Y1Y2, &hY1Y1Y2)) ||
      (err = OldHash(&C, Y1Y1Y1Y1, &hY1Y1Y1Y1)))
    goto FREE_RETURN;

  // Initialize basic factors diffusion:

  // Y0 ==================================
  C.h_y0 = hY0;
  pZ = C.pool[hY0];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = 0; // const
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = 0; // const
  pZ->volZ[0].coef[0][0] = sigma;

  // X1 ==================================
  pZ = C.pool[hX1];

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = 0; // const
  pZ->volZ[1].coef[0][0] = alpha;

  // X2 ==================================
  pZ = C.pool[hX2];

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = X1;
  pZ->volZ[1].coef[0][0] = alpha;

  // X3 ==================================
  pZ = C.pool[hX3];

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = X2;
  pZ->volZ[1].coef[0][0] = alpha;

  // X4 ==================================
  pZ = C.pool[hX4];

  pZ->volZ[1].n_mono = 1;
  pZ->volZ[1].mono[0] = X3;
  pZ->volZ[1].coef[0][0] = alpha;

  // Y1 ==================================
  pZ = C.pool[hY1];

  pZ->Z[0].n_mono = 1;
  pZ->Z[0].mono[0] = X1;
  pZ->Z[0].coef[0][0] = -sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = X1;
  pZ->volZ[0].coef[0][0] = sigma;

  // Y2 ==================================
  pZ = C.pool[hY2];

  pZ->Z[0].n_mono = 2;
  pZ->Z[0].mono[0] = X1X1;
  pZ->Z[0].coef[0][0] = -0.5 * sigma * sigma;
  pZ->Z[0].mono[1] = X2;
  pZ->Z[0].coef[1][0] = -sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = X2;
  pZ->volZ[0].coef[0][0] = sigma;

  // Y3 ==================================
  pZ = C.pool[hY3];

  pZ->Z[0].n_mono = 2;
  pZ->Z[0].mono[0] = X3;
  pZ->Z[0].coef[0][0] = -sigma * sigma;
  pZ->Z[0].mono[1] = X1X2;
  pZ->Z[0].coef[1][0] = -sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = X3;
  pZ->volZ[0].coef[0][0] = sigma;

  // Y4 ==================================
  pZ = C.pool[hY4];

  pZ->Z[0].n_mono = 3;
  pZ->Z[0].mono[0] = X4;
  pZ->Z[0].coef[0][0] = -sigma * sigma;
  pZ->Z[0].mono[1] = X1X3;
  pZ->Z[0].coef[1][0] = -sigma * sigma;
  pZ->Z[0].mono[2] = X2X2;
  pZ->Z[0].coef[2][0] = -0.5 * sigma * sigma;

  pZ->volZ[0].n_mono = 1;
  pZ->volZ[0].mono[0] = X4;
  pZ->volZ[0].coef[0][0] = sigma;

  // Calculate <Z  ,Y0> for all basic factors:

  for (i = 0; i < C.nfac; i++) {
    err = OldHash(&C, (1 << (i * C.nbit)), &hZ);
    if (err)
      goto FREE_RETURN;

    pZ = C.pool[hZ];
    memset(C.local_pool, -1, C.pool_size * sizeof(int));

    for (ki = 0; ki < C.nbr; ki++) {
      for (kj = 0; kj < C.nbr; kj++) {
        err = AddProdMonoPolyPoly(&C, &pZ->Z[1], 1.0, C.rho[ki][kj], 0,
                                  &pZ->volZ[ki], &C.pool[hY0]->volZ[kj]);
        if (err)
          goto FREE_RETURN;
      }
    }
  }

  // Initial factors initialized  , now build the system of ODE:

  err = OldFillStructZ(&C, Y1);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y2);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y3);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1Y1);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1Y2);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1Y1Y1);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y4);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y2Y2);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1Y3);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1Y1Y2);
  if (err)
    goto FREE_RETURN;

  err = OldFillStructZ(&C, Y1Y1Y1Y1);
  if (err)
    goto FREE_RETURN;

  // Integrate the system:

  err = OldIntegrateD(&C, &t, y);
  if (err)
    goto FREE_RETURN;

  // Return results:

  for (i = 0; i < 4; i++)
    memset(d[i], 0, nd * sizeof(double));

  if (C.pool[hY1]->nd > nd || C.pool[hY2]->nd > nd || C.pool[hY3]->nd > nd ||
      C.pool[hY1Y1]->nd + 1 > nd || C.pool[hY1Y2]->nd + 1 > nd ||
      C.pool[hY1Y1Y1]->nd + 2 > nd ||
      C.pool[hY1Y1Y1Y1]->nd + 3 >
          nd) // Anyway  , there are the most of equations
  {
    err = serror("nd too small");
    goto FREE_RETURN;
  }

  for (i = 0; i < C.pool[hY1]->nd; i++)
    d[0][i] += y[C.pool[hY1]->d_idx + i];
  for (i = 0; i < C.pool[hY2]->nd; i++)
    d[1][i] += y[C.pool[hY2]->d_idx + i];
  for (i = 0; i < C.pool[hY3]->nd; i++)
    d[2][i] += y[C.pool[hY3]->d_idx + i];
  for (i = 0; i < C.pool[hY1Y1]->nd; i++)
    d[1][i + 1] += 0.5 * y[C.pool[hY1Y1]->d_idx + i];
  for (i = 0; i < C.pool[hY1Y2]->nd; i++)
    d[2][i + 1] += y[C.pool[hY1Y2]->d_idx + i];
  for (i = 0; i < C.pool[hY1Y1Y1]->nd; i++)
    d[2][i + 2] += y[C.pool[hY1Y1Y1]->d_idx + i] / 6.0;
  for (i = 0; i < C.pool[hY4]->nd; i++)
    d[3][i] += y[C.pool[hY4]->d_idx + i];
  for (i = 0; i < C.pool[hY2Y2]->nd; i++)
    d[3][i + 1] += 0.5 * y[C.pool[hY2Y2]->d_idx + i];
  for (i = 0; i < C.pool[hY1Y3]->nd; i++)
    d[3][i + 1] += y[C.pool[hY1Y3]->d_idx + i];
  for (i = 0; i < C.pool[hY1Y1Y2]->nd; i++)
    d[3][i + 2] += 0.5 * y[C.pool[hY1Y1Y2]->d_idx + i];
  for (i = 0; i < C.pool[hY1Y1Y1Y1]->nd; i++)
    d[3][i + 3] += y[C.pool[hY1Y1Y1Y1]->d_idx + i] / 24.0;

FREE_RETURN:
  OldFreeControl(&C);

  return err;
}

Err ord_mult_const(SOrdSeries *s, double c) {
  int i;
  double *p = s->val + Expansion_ZEROPOS;
  for (i = s->min; i <= s->max; i++)
    p[i] *= c;
  return NULL;
}

Err ord_sum(SOrdSeries *s, SOrdSeries *s1) {
  int i;
  double *p = s->val + Expansion_ZEROPOS, *p1 = s1->val + Expansion_ZEROPOS;

  if (s1->min < s->min)
    s->min = s1->min;
  if (s1->max > s->max)
    s->max = s1->max;

  for (i = s->min; i <= s->max; i++)
    p[i] += p1[i];
  return NULL;
}

Err ord_product(SOrdSeries *s, SOrdSeries *s1) {
  int i, j;
  SOrdSeries tmp;
  double *p = s->val + Expansion_ZEROPOS, *p1 = s1->val + Expansion_ZEROPOS,
         *pt = tmp.val + Expansion_ZEROPOS;

  memset(&tmp, 0, sizeof(SOrdSeries));

  tmp.min = s->min + s1->min;
  tmp.max = s->max + s1->max;

  if (tmp.min + Expansion_ZEROPOS < 0 ||
      tmp.max + Expansion_ZEROPOS >= Expansion_MAXORDER)
    return serror("Out of bounds in ord_product");

  for (i = s->min; i <= s->max; i++)
    for (j = s1->min; j <= s1->max; j++)
      pt[i + j] += p[i] * p1[j];

  memcpy(s, &tmp, sizeof(SOrdSeries));
  return NULL;
}

Err ord_shift(SOrdSeries *s, int shift) {
  SOrdSeries tmp;
  double *p = s->val + Expansion_ZEROPOS, *pt = tmp.val + Expansion_ZEROPOS;

  memset(&tmp, 0, sizeof(SOrdSeries));

  tmp.min = s->min + shift;
  tmp.max = s->max + shift;

  if (tmp.min + Expansion_ZEROPOS < 0 ||
      tmp.max + Expansion_ZEROPOS >= Expansion_MAXORDER)
    return serror("Out of bounds in ord_shift");

  memcpy(pt + tmp.min, p + s->min, (s->max - s->min + 1) * sizeof(double));
  memcpy(s, &tmp, sizeof(SOrdSeries));
  return NULL;
}

#define ORDER 2

Err calc_impvols(double **d, int *uplus, int *bplus, int *splus, int m, int n,
                 double eps, double *yk, int nyk, double **vols) {
  Err err = NULL;
  SOrdSeries **f = NULL, *D[ORDER] = {NULL, NULL}, ftmp, s[ORDER], s1s1;
  int i, j, k, ord;

  // Allocate memory:

  f = (SOrdSeries **)calloc(n - 1, sizeof(SOrdSeries *));
  if (!f) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }
  memset(f, 0, (n - 1) * sizeof(SOrdSeries *));

  for (i = 0; i < n - 1; i++) {
    f[i] = (SOrdSeries *)calloc(nyk, sizeof(SOrdSeries));
    if (!f[i]) {
      err = serror("Memory failure");
      goto FREE_RETURN;
    }
    memset(f[i], 0, sizeof(SOrdSeries));
  }

  for (i = 0; i < ORDER; i++) {
    D[i] = (SOrdSeries *)calloc(n - 1, sizeof(SOrdSeries));
    if (!D[i]) {
      err = serror("Memory failure");
      goto FREE_RETURN;
    }
    memset(D[i], 0, sizeof(SOrdSeries));
  }

  // Calculate all D's:

  for (i = n - 1; i >= 1; i--) {
    if (i < n - 1)
      for (j = 0; j < ORDER; j++)
        memcpy(&D[j][i - 1], &D[j][i], sizeof(SOrdSeries));

    for (j = 0; j < m; j++) {
      if (bplus[j] >= ORDER) {
        err = serror("Invalid bplus");
        goto FREE_RETURN;
      }

      if (i - uplus[j] >= 0) {
        ord = i - uplus[j] + splus[j];
        D[bplus[j]][i - 1].val[Expansion_ZEROPOS + ord] += d[j][i - uplus[j]];
        if (D[bplus[j]][i - 1].min > ord)
          D[bplus[j]][i - 1].min = ord;
        if (D[bplus[j]][i - 1].max < ord)
          D[bplus[j]][i - 1].max = ord;
      }
    }
  }

  // For all strikes:

  for (j = 0; j < nyk; j++) {
    // Calculate all f's:

    f[0][j].val[Expansion_ZEROPOS] = 1.0;

    for (i = 1; i < n - 1; i++) {
      memcpy(&f[i][j], &f[i - 1][j], sizeof(SOrdSeries));
      err = ord_mult_const(&f[i][j], -0.5);
      if (err)
        goto FREE_RETURN;

      memcpy(&ftmp, &f[i - 1][j], sizeof(SOrdSeries));
      if ((err = ord_mult_const(&ftmp, -yk[j] / (eps * eps))) ||
          (err = ord_shift(&ftmp, -1)) || (err = ord_sum(&f[i][j], &ftmp)))
        goto FREE_RETURN;

      if (i >= 2) {
        memcpy(&ftmp, &f[i - 2][j], sizeof(SOrdSeries));
        if ((err = ord_mult_const(&ftmp, -(i - 1) / (eps * eps))) ||
            (err = ord_shift(&ftmp, -2)) || (err = ord_sum(&f[i][j], &ftmp)))
          goto FREE_RETURN;
      }
    }

    // Calculate sigmas:

    for (k = 0; k < ORDER; k++) {
      memcpy(&s[k], &D[k][0], sizeof(SOrdSeries));

      for (i = 1; i < n - 1; i++) {
        memcpy(&ftmp, &f[i][j], sizeof(SOrdSeries));
        if ((err = ord_product(&ftmp, &D[k][i])) ||
            ((i % 2) && (err = ord_mult_const(&ftmp, -1.0))) ||
            (err = ord_sum(&s[k], &ftmp)))
          goto FREE_RETURN;
      }
    }

    if ((err = ord_mult_const(&s[0], 1.0 / eps)) ||
        (err = ord_shift(&s[0], -1)))
      goto FREE_RETURN;

    memcpy(&s1s1, &s[0], sizeof(SOrdSeries));
    err = ord_product(&s1s1, &s[0]);
    if (err)
      goto FREE_RETURN;

    memcpy(&ftmp, &s1s1, sizeof(SOrdSeries));
    if ((err = ord_mult_const(&ftmp, -0.5 * yk[j] * yk[j] / (eps * eps))) ||
        (err = ord_sum(&s[1], &ftmp)))
      goto FREE_RETURN;

    memcpy(&ftmp, &s1s1, sizeof(SOrdSeries));
    if ((err = ord_mult_const(&ftmp, eps * eps / 8.0)) ||
        (err = ord_shift(&ftmp, 2)) || (err = ord_sum(&s[1], &ftmp)))
      goto FREE_RETURN;

    if ((err = ord_mult_const(&s[1], 1.0 / eps)) ||
        (err = ord_shift(&s[1], -1)))
      goto FREE_RETURN;

    // Return results:

    for (k = 0; k < ORDER; k++) {
      vols[k][j] = 0.0;
      for (i = s[k].min; i <= ORDER - k; i++)
        vols[k][j] += s[k].val[Expansion_ZEROPOS + i];
    }
  }

FREE_RETURN:
  if (f)
    for (i = 0; i < n - 1; i++)
      free(f[i]);
  free(f);
  free(D[0]);
  free(D[1]);

  return err;
}
