// AffineCoFi.c:  Cornish-Fisher expansion for affine models

#include "math.h"
#include "srt_h_all.h"
#undef SIGN
#include "AffineCoFi.h"
#include "RandomGen.h"
#include "nag.h"
#include "nagd02.h"
#include "srt_h_hash.h"

typedef struct _SAffiCoef {
  double **B, ***D;
  int n, M;
  SHash *hash;
} SAffiCoef;

static void NAG_CALL EvalDerivatives(Integer neq, double t, double y[],
                                     double yp[], Nag_User *comm) {
  SAffiCoef *p = (SAffiCoef *)comm->p;
  int i, j, k, iim, ii[10], iisum[10], m_idx, jj[10], jsum, jjcoef, mjj_idx,
      cjj[10], cmjj_idx;
  double coef, coef2;

  memset(yp, 0, neq * sizeof(double));

  for (iim = 1; iim <= p->M; iim++) // order of the moment
  {
    memset(ii, 0,
           p->n * sizeof(int)); // order of differentiation on each variable
    memset(iisum, 0, p->n * sizeof(int));
    ii[p->n - 1] = iim;

    for (;;) {
      Hash(p->hash, ii, &m_idx);
      m_idx = (m_idx - 1) * (p->n + 1);

      // Drift:
      for (i = 1; i <= p->n; i++)
        for (j = 0; j <= p->n; j++)
          yp[m_idx + j] -= y[m_idx + i] * p->B[i - 1][j];

      // Cross products:

      memset(jj, 0, p->n * sizeof(int));
      for (k = p->n - 1; ii[k] == 0; k--)
        ;
      jj[k] = jsum = 1;
      jjcoef = ii[k];

      if (iim > 1)
        for (;;) {
          Hash(p->hash, jj, &mjj_idx);
          mjj_idx = (mjj_idx - 1) * (p->n + 1);

          for (k = 0; k < p->n; k++)
            cjj[k] = ii[k] - jj[k];
          Hash(p->hash, cjj, &cmjj_idx);
          cmjj_idx = (cmjj_idx - 1) * (p->n + 1);

          coef = (jsum % 2 - 0.5 * (iim % 2 + 1)) * jjcoef;
          coef2 = ((1 - iim % 2) * (jsum % 2) * 2 - 1) * jjcoef;

          for (i = 1; i <= p->n; i++) {
            for (k = 0; k <= p->n; k++)
              yp[m_idx + k] += coef * y[mjj_idx + i] * y[cmjj_idx + i] *
                               p->D[i - 1][i - 1][k];

            for (j = i + 1; j <= p->n; j++)
              for (k = 0; k <= p->n; k++)
                yp[m_idx + k] += coef2 * y[mjj_idx + i] * y[cmjj_idx + j] *
                                 p->D[i - 1][j - 1][k];
          }

          // increment jj:

          for (k = p->n - 1; k > 0 && jj[k] == ii[k]; k--)
            ;
          if (k == 0 && jj[0] >= ii[0] - 1)
            break;

          jjcoef = jjcoef * (ii[k] - jj[k]) / (jj[k] + 1);
          jj[k]++;
          jsum++;
          for (k++; k < p->n; k++) {
            jsum -= jj[k];
            jj[k] = 0;
          }
        }

      // increment ii:
      for (k = p->n - 2; k > 0 && ii[k] == iim - iisum[k - 1]; k--)
        ;
      if (k <= 0 && ii[0] == iim)
        break;

      ii[k]++;
      iisum[k]++;
      for (k++; k < p->n - 1; k++) {
        ii[k] = 0;
        iisum[k] = iisum[k - 1];
      }
      ii[k] = iim - iisum[k - 1];
    }
  }
}

static void FillAffiCoefHeston(SAffiCoef *p, double sigma, double lambda,
                               double alpha, double rho) {
  // Drift:
  p->B[0][0] = 0.0;
  p->B[0][1] = 0.0;
  p->B[0][2] = -0.5 * sigma * sigma;

  p->B[1][0] = lambda;
  p->B[1][1] = 0.0;
  p->B[1][2] = -lambda;

  // Cross products:
  p->D[0][0][0] = 0.0;
  p->D[0][0][1] = 0.0;
  p->D[0][0][2] = sigma * sigma;

  p->D[0][1][0] = p->D[1][0][0] = 0.0;
  p->D[0][1][1] = p->D[1][0][1] = 0.0;
  p->D[0][1][2] = p->D[1][0][2] = alpha * rho * sigma;

  p->D[1][1][0] = 0.0;
  p->D[1][1][1] = 0.0;
  p->D[1][1][2] = alpha * alpha;
}

#define MAXNEQ 1000

static Err HestonDMoments(double sigma, double lambda, double alpha, double rho,
                          double mat, int max_order, SPolynomial *poly) {
  Err err = NULL;
  int k, neq;
  SAffiCoef affi;
  int iim, ii[10], iisum[10], m_idx;
  long fact[10];
  double coef;
  Nag_User comm_Nag;
  NagError fail;
  Nag_ODE_RK opt;
  double y[MAXNEQ], yp[MAXNEQ], ymax[MAXNEQ], RKthres[MAXNEQ], tgot;
  const double RKtol = 1e-5,
               thres_min =
                   1e-8; //  , thres_coef = 1e-7;		// adjustable

  memset(&fail, 0, sizeof(NagError));
  memset(&opt, 0, sizeof(Nag_ODE_RK));
  memset(&affi, 0, sizeof(SAffiCoef));

  if (poly->hash->pos > 1)
    return serror("Second hand hash - oops !");

  comm_Nag.p = &affi;
  affi.hash = poly->hash;
  affi.n = 2; // Heston specific
  affi.M = max_order;
  affi.B = dmatrix(0, affi.n - 1, 0, affi.n);
  affi.D = f3tensor(0, affi.n - 1, 0, affi.n - 1, 0, affi.n);

  if (!affi.B || !affi.D) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }

  FillAffiCoefHeston(&affi, sigma, lambda, alpha, rho); // Heston specific

  // Calculate number of equations and initialize the hash:

  for (iim = 1; iim <= affi.M; iim++) // order of the moment
  {
    memset(ii, 0,
           affi.n * sizeof(int)); // order of differentiation on each variable
    memset(iisum, 0, affi.n * sizeof(int));
    ii[affi.n - 1] = iim;

    for (;;) {
      err = Hash(affi.hash, ii, &m_idx);
      if (err)
        goto FREE_RETURN;

      // increment ii:
      for (k = affi.n - 2; k > 0 && ii[k] == iim - iisum[k - 1]; k--)
        ;
      if (k <= 0 && ii[0] == iim)
        break;

      ii[k]++;
      iisum[k]++;
      for (k++; k < affi.n - 1; k++) {
        ii[k] = 0;
        iisum[k] = iisum[k - 1];
      }
      ii[k] = iim - iisum[k - 1];
    }
  }

  neq = (affi.hash->pos - 1) * (affi.n + 1);
  memset(y, 0, neq * sizeof(double)); // Final y values are almost all zeros
  memset(ii, 0, affi.n * sizeof(int));

  for (k = 0; k < affi.n; k++) {
    ii[k]++;
    Hash(affi.hash, ii, &m_idx);
    m_idx = (m_idx - 1) * (affi.n + 1);
    y[m_idx + k + 1] = 1.0;
    ii[k]--;
  }

  for (k = 0; k < neq; k++)
    RKthres[k] = thres_min;

  // Calculate solution at time 0

  nag_ode_ivp_rk_setup(neq, mat, y, 0.0, RKtol, RKthres, Nag_RK_4_5,
                       Nag_RK_range, Nag_ErrorAssess_off, 0.0, &opt, &fail);
  if (fail.code != NE_NOERROR) {
    err = serror(fail.message);
    goto FREE_RETURN;
  }

  nag_ode_ivp_rk_range(neq, EvalDerivatives, 0.0, &tgot, y, yp, ymax, &opt,
                       &comm_Nag, &fail);
  if (fail.code != NE_NOERROR) {
    err = serror(fail.message);
    goto FREE_RETURN;
  }

  /*	nag_ode_ivp_rk_free(&opt);

          for (k=0; k < neq; k++)
          {
                  RKthres[k] = thres_coef * ymax[k];
                  if (RKthres[k] < thres_min) RKthres[k] = thres_min;
          }
  */

  // Calculate the coefficients of the resulting polynomial:

  fact[0] = fact[1] = 1;
  for (k = 2; k <= max_order; k++)
    fact[k] = fact[k - 1] * k;

  for (iim = 1; iim <= affi.M; iim++) // order of the moment
  {
    memset(ii, 0,
           affi.n * sizeof(int)); // order of differentiation on each variable
    memset(iisum, 0, affi.n * sizeof(int));
    ii[affi.n - 1] = iim;

    for (;;) {
      Hash(affi.hash, ii, &m_idx);
      m_idx = (m_idx - 1) * (affi.n + 1);

      coef = y[m_idx] + y[m_idx + 2]; // Heston specific	(later: force all
                                      // factors initial values = 0)
      for (k = 0; k < affi.n; k++)
        coef /= fact[ii[k]];

      // Add the term to the polynomial:

      if (fabs(coef) > 1e-16) {
        err = Poly_AddConstMono(poly, coef, NULL, ii);
        if (err)
          goto FREE_RETURN;
      }

      // increment ii:
      for (k = affi.n - 2; k > 0 && ii[k] == iim - iisum[k - 1]; k--)
        ;
      if (k <= 0 && ii[0] == iim)
        break;

      ii[k]++;
      iisum[k]++;
      for (k++; k < affi.n - 1; k++) {
        ii[k] = 0;
        iisum[k] = iisum[k - 1];
      }
      ii[k] = iim - iisum[k - 1];
    }
  }

FREE_RETURN:
  nag_ode_ivp_rk_free(&opt);

  if (affi.B)
    free_dmatrix(affi.B, 0, affi.n - 1, 0, affi.n);
  if (affi.D)
    free_f3tensor(affi.D, 0, affi.n - 1, 0, affi.n - 1, 0, affi.n);

  return err;
}

Err HestonDDensityFTExpansion(double sigma, double alpha, double rho,
                              double lam, double mat, int order, int nu,
                              double *u_re, double *u_im, double *h_re,
                              double *h_im) {
  Err err = NULL;
  SHash hash;
  SPolynomial pp, poly;
  double uu_re, uu_im, tmp_re, tmp_im;
  int i, j, ii[2] = {0, 0}, hii;

  memset(&hash, 0, sizeof(SHash));
  memset(&poly, 0, sizeof(SPolynomial));
  memset(&pp, 0, sizeof(SPolynomial));

  err = InitHash(&hash, 2, MAXNEQ);
  if (err)
    goto FREE_RETURN;

  err = Poly_Init(&poly, &hash, MAXNEQ, 1);
  if (err)
    goto FREE_RETURN;

  err = Poly_Init(&pp, &hash, MAXNEQ, 1);
  if (err)
    goto FREE_RETURN;

  err = HestonDMoments(sigma, lam, alpha, rho, mat, order, &poly);
  if (err)
    goto FREE_RETURN;

  for (ii[0] = 1; ii[0] <= order; ii[0]++) {
    Hash(&hash, ii, &hii);

    if (poly.pool[hii] >= 0) {
      err = Poly_AddConstMono(
          &pp, (ii[0] % 4 >= 2 ? -1.0 : 1.0) * poly.coef[poly.pool[hii]][0],
          NULL, poly.mono[poly.pool[hii]]);
      if (err)
        goto FREE_RETURN;
    }
  }

  for (i = 0; i < nu; i++) {
    h_re[i] = h_im[i] = 0.0;
    for (ii[0] = 1; ii[0] <= order; ii[0]++) {
      uu_re = -u_im[i];
      uu_im = u_re[i];
      for (j = 1; j < ii[0]; j++) {
        tmp_re = -uu_re * u_im[i] - uu_im * u_re[i];
        tmp_im = uu_re * u_re[i] - uu_im * u_im[i];
        uu_re = tmp_re;
        uu_im = tmp_im;
      }
      Hash(&hash, ii, &hii);

      h_re[i] += uu_re * pp.coef[pp.pool[hii]][0];
      h_im[i] += uu_im * pp.coef[pp.pool[hii]][0];
    }
  }

FREE_RETURN:
  Poly_Free(&poly);
  Poly_Free(&pp);
  FreeHashTree(hash.hash_tree);

  return err;
}

static Err calc_g_deriv(int *iuv, SHash *uvHash, SPolynomial *g_pool,
                        SPolynomial *g) {
  Err err = NULL;
  int i, huv, ii[10], i1[10], hi, h1;

  err = Hash(uvHash, iuv, &huv);
  if (err)
    return err;

  if (g_pool[huv].mono)
    return NULL;

  err = Poly_Init(&g_pool[huv], g->hash, 5000, 1);
  if (err)
    return err;

  if (huv == 0) {
    g_pool[huv].n_mono = 1;
    memset(g_pool[huv].mono[0], 0, g->hash->nfac * sizeof(int)); // const
    g_pool[huv].coef[0][0] = 1.0;
    return NULL;
  }

  memset(i1, 0, uvHash->nfac * sizeof(int));
  for (i = 0; iuv[i] == 0; i++)
    ;
  i1[i]++;
  err = Hash(uvHash, i1, &h1);
  if (err)
    return err;

  if (huv == h1) {
    err = Poly_AddDPoly(&g_pool[huv], g, i);
    return err;
  }

  memcpy(ii, iuv, uvHash->nfac * sizeof(int));
  ii[i]--;
  err = Hash(uvHash, ii, &hi);
  if (err)
    return err;

  if (!g_pool[hi].mono) {
    err = calc_g_deriv(ii, uvHash, g_pool, g);
    if (err)
      return err;
  }

  if (!g_pool[h1].mono) {
    err = calc_g_deriv(i1, uvHash, g_pool, g);
    if (err)
      return err;
  }

  err = Poly_AddDPoly(&g_pool[huv], &g_pool[hi], i);
  if (err)
    return err;

  err = Poly_AddConstMonoPolyPoly(&g_pool[huv], 1.0, NULL, NULL, &g_pool[h1],
                                  &g_pool[hi]);
  if (err)
    return err;

  return NULL;
}

Err HestonDCoFiMC(double sigma, double alpha, double rho, double lam,
                  double mat, int order, int cf_ord, long npaths, double fwd,
                  int nK, double *K, double *pv, double *stddev) {
  Err err = NULL;
  SHash hash;
  SPolynomial poly, pp, dd, g_pool[MAXNEQ], phi, herm[MAXNEQ], a[10], corr,
      corr_cut, tmp[2];
  int i, j, r, ii[2] = {0, 0}, hii, hii1;
  int *reject_path = NULL, npaths_used = npaths;
  long fact[10], seed = -12345678;
  double d, *YY = NULL, **pvs = NULL, WW, dY, stdy, meany;
  SRandomGen rg;

  memset(&hash, 0, sizeof(SHash));
  memset(&poly, 0, sizeof(SPolynomial));
  memset(&pp, 0, sizeof(SPolynomial));
  memset(&dd, 0, sizeof(SPolynomial));
  memset(g_pool, 0, MAXNEQ * sizeof(SPolynomial));
  memset(&phi, 0, sizeof(SPolynomial));
  memset(herm, 0, MAXNEQ * sizeof(SPolynomial));
  memset(a, 0, 10 * sizeof(SPolynomial));
  memset(&corr, 0, sizeof(SPolynomial));
  memset(&corr_cut, 0, sizeof(SPolynomial));
  memset(tmp, 0, 2 * sizeof(SPolynomial));
  memset(&rg, 0, sizeof(SRandomGen));

  err = InitHash(&hash, 2, MAXNEQ);
  if (err)
    goto FREE_RETURN;

  err = Poly_Init(&poly, &hash, MAXNEQ, 1);
  if (err)
    goto FREE_RETURN;

  err = Poly_Init(&pp, &hash, MAXNEQ, 1);
  if (err)
    goto FREE_RETURN;

  err = HestonDMoments(sigma, lam, alpha, rho, mat, order, &poly);
  if (err)
    goto FREE_RETURN;

  for (ii[0] = 1; ii[0] <= order; ii[0]++) {
    Hash(&hash, ii, &hii);

    if (poly.pool[hii] >= 0) {
      err = Poly_AddConstMono(
          &pp, (ii[0] % 4 >= 2 ? -1.0 : 1.0) * poly.coef[poly.pool[hii]][0],
          NULL, poly.mono[poly.pool[hii]]);
      if (err)
        goto FREE_RETURN;
    }
  }
  ii[0] = 2;
  err = Poly_AddConstMono(&pp, -0.5, NULL, ii);
  if (err)
    goto FREE_RETURN;

  err = Poly_Init(&phi, &hash, 1, 1);
  if (err)
    goto FREE_RETURN;

  err = Poly_AddConstMono(&phi, -0.5, NULL, ii);
  if (err)
    goto FREE_RETURN;

  // init fact:
  fact[0] = fact[1] = 1;
  for (i = 2; i <= order; i++)
    fact[i] = fact[i - 1] * i;

  err = Poly_Init(&dd, &hash, 10, 1);
  if (err)
    goto FREE_RETURN;

  // calculate the dd coefficients:
  for (ii[0] = 1; ii[0] <= order; ii[0]++) {
    Hash(&hash, ii, &hii);

    err = calc_g_deriv(ii, &hash, g_pool, &pp);
    if (err)
      goto FREE_RETURN;

    err =
        Poly_AddConstMono(&dd, g_pool[hii].coef[0][0] / fact[ii[0]], NULL, ii);
    if (err)
      goto FREE_RETURN;
  }

  // Calculate a:
  for (i = 0; i < cf_ord; i++) {
    err = Poly_Init(&a[i], &hash, 5000, 1);
    if (err)
      goto FREE_RETURN;
  }

  for (ii[0] = 1; ii[0] <= order; ii[0]++) {
    Hash(&hash, ii, &hii);
    ii[0]--;
    Hash(&hash, ii, &hii1);

    err = calc_g_deriv(ii, &hash, herm, &phi);
    if (err)
      goto FREE_RETURN;
    ii[0]++;

    err = Poly_AddConstMonoPoly(
        &a[0], (ii[0] % 2 ? -1.0 : 1.0) * dd.coef[dd.pool[hii]][0], NULL, NULL,
        &herm[hii1]);
    if (err)
      goto FREE_RETURN;
  }

  // Calculate all the powers of a:

  for (r = 1; r < cf_ord; r++) {
    err = Poly_AddConstMonoPolyPoly(&a[r], 1.0, NULL, NULL, &a[0], &a[r - 1]);
    if (err)
      goto FREE_RETURN;
  }

  // Now compute the Cornish-Fisher expansion:

  err = Poly_Init(&corr, &hash, 50000, 1);
  if (err)
    goto FREE_RETURN;

  for (r = 0; r < 2; r++) {
    err = Poly_Init(&tmp[r], &hash, 50000, 1);
    if (err)
      goto FREE_RETURN;
  }

  ii[0] = 1;
  Hash(&hash, ii, &hii1);

  for (r = 0, d = 1.0; r < cf_ord; r++) {
    d = -d / (r + 1);

    Poly_Copy(&tmp[0], &a[r]);
    for (i = 0; i < r; i++) {
      Poly_Clear(&tmp[1]);
      err = Poly_AddDPoly(&tmp[1], &tmp[0], 0);
      if (err)
        goto FREE_RETURN;

      err = Poly_AddConstMonoPolyPoly(&tmp[1], r - i, NULL, NULL, &herm[hii1],
                                      &tmp[0]);
      if (err)
        goto FREE_RETURN;

      Poly_Copy(&tmp[0], &tmp[1]);
    }
    err = Poly_AddConstMonoPoly(&corr, d, NULL, NULL, &tmp[0]);
    if (err)
      goto FREE_RETURN;
  }

  // Now cut higher orders of eps and eps_b in corr:

  err = Poly_Init(&corr_cut, &hash, 2000, 1);
  if (err)
    goto FREE_RETURN;

  for (i = 0; i < corr.n_mono; i++) {
    if (corr.mono[i][0] <= order) {
      if (corr_cut.n_mono >= corr_cut.maxmono) {
        err = serror("Max mono exceeded");
        goto FREE_RETURN;
      }

      memcpy(corr_cut.mono[corr_cut.n_mono], corr.mono[i],
             hash.nfac * sizeof(int));
      corr_cut.coef[corr_cut.n_mono][0] = corr.coef[i][0];
      corr_cut.n_mono++;
    }
  }

  // ========================================================
  // Correction polynomial now calculated - start simulation:

  YY = (double *)calloc(npaths, sizeof(double));
  pvs = dmatrix(0, npaths - 1, 0, nK - 1);
  reject_path = (int *)calloc(npaths, sizeof(int));

  if (!YY || !reject_path || !pvs) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }

  err = ABS_Init(&rg, seed, npaths, 1, 0);
  if (err)
    goto FREE_RETURN;

  // Generate all paths simultaneously:

  ii[0] = 1;
  Hash(&hash, ii, &hii);
  meany = poly.coef[poly.pool[hii]][0];
  ii[0]++;
  Hash(&hash, ii, &hii);
  stdy = sqrt(-poly.coef[poly.pool[hii]][0] * 2.0);

  for (r = 0; r < npaths; r++) {
    err = rg.Gauss(&rg, &WW);
    if (err)
      goto FREE_RETURN;

    if (!reject_path[r]) {
      // Correct dY:

      dY = WW;
      for (i = 0; i < corr_cut.n_mono; i++) {
        d = 1.0;
        for (j = 0; j < corr_cut.mono[i][0]; j++)
          d *= WW;
        dY += d * corr_cut.coef[i][0];
      }
      YY[r] += dY;

      if (fabs(YY[r] - meany) > 10.0 * stdy) {
        reject_path[r] = 1;
        npaths_used--;
        if (npaths / npaths_used >= 2) {
          err = serror("At least half of paths rejected - check parameters");
          goto FREE_RETURN;
        }
      }
    }
  }

  // Calculate pv and stddev:
  memset(pv, 0, nK * sizeof(double));
  memset(stddev, 0, nK * sizeof(double));

  for (i = 0; i < nK; i++) {
    for (j = 0; j < npaths; j++)
      if (!reject_path[j]) {
        pvs[j][i] = fwd * exp(YY[j]) - K[i];
        if (pvs[j][i] < 0.0)
          pvs[j][i] = 0.0;

        pv[i] += pvs[j][i] / npaths_used;
      }
    for (j = 0; j < npaths; j++)
      if (!reject_path[j])
        stddev[i] +=
            (pvs[j][i] - pv[i]) * (pvs[j][i] - pv[i]) / (npaths_used - 1);
    stddev[i] = sqrt(stddev[i] / npaths_used);
  }

FREE_RETURN:
  Poly_Free(&poly);
  Poly_Free(&pp);
  Poly_Free(&dd);
  Poly_Free(&phi);
  Poly_Free(&corr);
  Poly_Free(&corr_cut);
  for (i = 0; i < hash.pos; i++) {
    Poly_Free(&g_pool[i]);
    Poly_Free(&herm[i]);
  }
  for (i = 0; i < cf_ord; i++)
    Poly_Free(&a[i]);
  for (i = 0; i < 2; i++)
    Poly_Free(&tmp[i]);
  FreeHashTree(hash.hash_tree);

  ABS_Free(&rg);
  free(YY);
  free(reject_path);
  if (pvs)
    free_dmatrix(pvs, 0, npaths - 1, 0, nK - 1);

  return err;
}
