/*--------------------------------------------------------------
        FILE: CheyBetaPDE.cxx
        PURPOSE: Cheyette beta PDE (forward & backward)
        AUTHOR: Dimitri Mayevski
        DATE: 22/11/2002
  --------------------------------------------------------------*/

#include "CheyBetaPDE.h"
#include "math.h"
#include "srt_h_all.h"

static void FreeGrid(SCheyBetaPDEGrid *g) {
  free(g->phi);
  if (g->ad)
    free_dmatrix(g->ad, 0, g->nphi - 1, 0, g->nx - 1);
  if (g->pv)
    free_f3tensor(g->pv, 0, g->nphi - 1, 0, g->nx - 1, 0, g->ninst - 1);
  memset(g, 0, sizeof(SCheyBetaPDEGrid));
}

static Err CopyGrid(SCheyBetaPDEGrid *g1, SCheyBetaPDEGrid *g2) {
  FreeGrid(g2);
  g2->nphi = g1->nphi;
  g2->nx = g1->nx;
  g2->ninst = g1->ninst;

  if (g1->phi) {
    g2->phi = (double *)calloc(g2->nphi, sizeof(double));
    if (!g2->phi)
      return serror("Memory failure");
    memcpy(g2->phi, g1->phi, g2->nphi * sizeof(double));
  }
  if (g1->ad) {
    g2->ad = dmatrix(0, g2->nphi - 1, 0, g2->nx - 1);
    if (!g2->ad)
      return serror("Memory failure");
    memcpy(&(g2->ad[0][0]), &(g1->ad[0][0]),
           g2->nphi * g2->nx * sizeof(double));
  }
  if (g1->pv) {
    g2->pv = f3tensor(0, g2->nphi - 1, 0, g2->nx - 1, 0, g2->ninst - 1);
    if (!g2->pv)
      return serror("Memory failure");
    memcpy(&(g2->pv[0][0][0]), &(g1->pv[0][0][0]),
           g2->nphi * g2->nx * g2->ninst * sizeof(double));
  }
  return NULL;
}

Err CheyBetaPDE_Init(SCheyBetaPDE *pde, SProductDesc *g, int nt, int nx,
                     int nphi, double cutcoef, int isfwd, SCheyBeta *pmdl) {
  int ntimes, i;
  long d1, d2;

  pde->g = g;
  pde->pmdl = pmdl;
  pde->cutcoef = (cutcoef < 50.0 ? cutcoef : 50.0);
  pde->isfwd = isfwd;

  if (g->nccy > 1)
    return serror("Error (CheyBeta PDE): Product is multi-currency");

  ntimes = g->nex + pmdl->nsig + 1;
  pde->times = (double *)calloc(ntimes, sizeof(double));
  pde->times[0] = 0.0;
  memcpy(pde->times + 1, g->ex, g->nex * sizeof(double));
  memcpy(pde->times + g->nex + 1, pmdl->sigtms, pmdl->nsig * sizeof(double));

  num_f_sort_vector(ntimes, pde->times);
  num_f_unique_vector(&ntimes, pde->times);
  while (pde->times[ntimes - 1] > g->ex[g->nex - 1] + 1e-5)
    ntimes--;
  num_f_fill_vector_newalgo(&ntimes, &pde->times,
                            (int)(nt * pde->times[ntimes - 1] + 1e-5));
  pde->nt = ntimes;

  // Calculate deterministic diffusion parameters
  pde->fr = (double *)calloc(ntimes - 1, sizeof(double));
  pde->gamT = (double *)calloc(ntimes, sizeof(double));
  pde->df = (double *)calloc(ntimes, sizeof(double));
  pde->phi_max = (double *)calloc(ntimes, sizeof(double));
  pde->phi_min = (double *)calloc(ntimes, sizeof(double));
  pde->phi_mid = (double *)calloc(ntimes, sizeof(double));
  pde->cumxvar = (double *)calloc(ntimes, sizeof(double));
  if (!pde->fr || !pde->gamT || !pde->df || !pde->phi_max || !pde->phi_min ||
      !pde->phi_mid || !pde->cumxvar)
    return serror("Memory failure");

  d1 = pmdl->today;
  pde->df[0] = 1.0;
  for (i = 0; i < ntimes - 1; i++) {
    d2 = (long)(pmdl->today + pde->times[i + 1] * DAYS_IN_YEAR + 1e-5);
    pde->fr[i] = swp_f_zr(d1, d2, pmdl->ycname);
    pde->df[i + 1] = swp_f_df(pmdl->today, d2, pmdl->ycname);
    pde->gamT[i] =
        (1.0 - exp(-pmdl->lambda * (pde->times[ntimes - 1] - pde->times[i]))) /
        pmdl->lambda;
    d1 = d2;
  }
  pde->gamT[ntimes - 1] = 0.0;
  pde->phi_max[0] = pde->phi_min[0] = pde->phi_mid[0] = pde->cumxvar[0] = 0.0;

  pde->last = 0;
  pde->nx = (nx >> 1) * 2 + 1; // Force pde->nx to be an odd number
  pde->nphi = (nphi > 3 ? nphi : 3);
  pde->dir = 0;
  pde->savetime = 0;

  pde->x = (double *)calloc(pde->nx, sizeof(double));
  if (!pde->x)
    return serror("Memory failure");

  pde->tmppde = (CNPDE_TEMP *)malloc(sizeof(CNPDE_TEMP));
  if (!pde->tmppde)
    return serror("Memory failure");

  num_f_pde_init(pde->tmppde, pde->nx, isfwd ? 1 : g->ninst);

  return NULL;
}

Err CheyBetaPDE_Free(SCheyBetaPDE *pde) {
  int i;
  for (i = 0; i < 3; i++)
    FreeGrid(&pde->grids[i]);
  free(pde->times);
  free(pde->df);
  free(pde->gamT);
  free(pde->fr);
  free(pde->phi_max);
  free(pde->phi_min);
  free(pde->phi_mid);
  free(pde->x);
  free(pde->cumxvar);

  if (pde->tmppde)
    num_f_pde_free(pde->tmppde, pde->nx, pde->isfwd ? 1 : pde->g->ninst);
  free(pde->tmppde);

  memset(pde, 0, sizeof(SCheyBetaPDE));
  return NULL;
}

static void MakeGrid(SCheyBetaPDE *pde) {
  SCheyBeta *pmdl = pde->pmdl;
  double std = 0.0, sig, t = 0.0, x_min, x_max, dx;
  int i;

  for (i = 0; i < pmdl->nsig && t < pde->times[pde->nt - 1] - 1e-5; i++) {
    sig = pmdl->sig[i];
    std += sig * sig * (pmdl->sigtms[i] - t);
    t = pmdl->sigtms[i];
  }
  std += sig * sig * (pde->times[pde->nt - 1] - t);
  std = sqrt(std);

  x_min = -3.0 * std;
  x_max = 5.0 * std;
  dx = (x_max - x_min) / (pde->nx - 1);
  pde->ix0 = (int)(-x_min / dx);
  pde->x[pde->ix0] = 0.0;
  for (i = pde->ix0 + 1; i < pde->nx; i++)
    pde->x[i] = pde->x[i - 1] + dx;
  for (i = pde->ix0 - 1; i >= 0; i--)
    pde->x[i] = pde->x[i + 1] - dx;
}

static void ExpandBounds(SCheyBetaPDE *pde, int start, int end) {
  double dt, v_min, v_max, v_mid, std, xcut_u, xcut_d;
  double sig, beta, lam = pde->pmdl->lambda;
  int i, isig = 0;

  for (i = start; i < end; i++) {
    dt = pde->times[i + 1] - pde->times[i];
    while (isig < pde->pmdl->nsig - 1 &&
           pde->pmdl->sigtms[isig] < pde->times[i] + 1e-5)
      isig++;
    sig = pde->pmdl->sig[isig];
    beta = pde->pmdl->beta[isig];

    std = sqrt(pde->cumxvar[i]);
    xcut_u = pde->cutcoef * pde->fr[i] * std;
    xcut_d = -pde->cutcoef * pde->fr[i] * std;
    v_mid = CheyBeta_Vol(0.0, pde->phi_mid[i], sig, beta, pde->gamT[i],
                         pde->fr[i], xcut_u, xcut_d);
    v_max = CheyBeta_Vol(xcut_u, 0.0, sig, beta, pde->gamT[i], pde->fr[i],
                         xcut_u, xcut_d);
    v_min = CheyBeta_Vol(xcut_d, 0.0, sig, beta, pde->gamT[i], pde->fr[i],
                         xcut_u, xcut_d);

    pde->cumxvar[i + 1] = pde->cumxvar[i] + v_mid * v_mid * dt;
    pde->phi_mid[i + 1] =
        pde->phi_mid[i] + (v_mid * v_mid - 2.0 * lam * pde->phi_mid[i]) * dt;
    pde->phi_max[i + 1] =
        pde->phi_max[i] + (v_max * v_max - 2.0 * lam * pde->phi_max[i]) * dt;
    pde->phi_min[i + 1] =
        pde->phi_min[i] + (v_min * v_min - 2.0 * lam * pde->phi_min[i]) * dt;
  }
}

static Err InitGridAt_i(SCheyBetaPDE *pde, int i) {
  int itmp, j;
  double dphi;

  if (pde->phi_max[i] - pde->phi_min[i] > 1e-12) {
    itmp = 3 + (int)((pde->nphi - 3) * pde->times[i] / pde->times[pde->nt - 1] +
                     1e-8);
    dphi = (pde->phi_max[i] - pde->phi_min[i]) / (itmp - 1);
    pde->grids[!pde->dir].phi = (double *)calloc(itmp, sizeof(double));
    if (!pde->grids[!pde->dir].phi)
      return serror("Memory failure");
    pde->grids[!pde->dir].nphi = itmp;
    for (j = 0; j < itmp; j++)
      pde->grids[!pde->dir].phi[j] = pde->phi_min[i] + j * dphi;
  } else // If phi is constant (i.e. LGM)
  {
    pde->grids[!pde->dir].phi = (double *)calloc(1, sizeof(double));
    if (!pde->grids[!pde->dir].phi)
      return serror("Memory failure");
    pde->grids[!pde->dir].nphi = 1;
    pde->grids[!pde->dir].phi[0] = pde->phi_mid[i];
  }

  pde->grids[!pde->dir].nx = pde->nx;
  if (pde->isfwd) {
    pde->grids[!pde->dir].ad =
        dmatrix(0, pde->grids[!pde->dir].nphi - 1, 0, pde->nx - 1);
    if (!pde->grids[!pde->dir].ad)
      return serror("Memory failure");

    memset(&(pde->grids[!pde->dir].ad[0][0]), 0,
           pde->grids[!pde->dir].nphi * pde->nx * sizeof(double));
  } else {
    pde->grids[!pde->dir].pv = f3tensor(0, pde->grids[!pde->dir].nphi - 1, 0,
                                        pde->nx - 1, 0, pde->g->ninst - 1);
    if (!pde->grids[!pde->dir].pv)
      return serror("Memory failure");
    pde->grids[!pde->dir].ninst = pde->g->ninst;

    memset(&(pde->grids[!pde->dir].pv[0][0][0]), 0,
           pde->grids[!pde->dir].nphi * pde->nx * pde->g->ninst *
               sizeof(double));
  }
  return NULL;
}

Err CheyBetaPDE_ExpandFwd(SCheyBetaPDE *pde, double T) {
  Err err = NULL;
  int new_last;
  double *var = NULL, *mu = NULL, *r = NULL, *ad = NULL;
  int i, isig = 0, xi, j, k;
  double sig, beta, lam, dt, std, xcut_u, xcut_d;
  double phi_k, phi_nxt, coef, sum;

  if (!pde->isfwd)
    return serror("The PDE was initialized as backward");

  // Find new last
  new_last = pde->last;
  while (new_last < pde->nt - 1 && pde->times[new_last] < T - 1e-5)
    new_last++;
  if (new_last == pde->last)
    return NULL;

  if (pde->last != pde->savetime) // Save grid
  {
    err = CopyGrid(&pde->grids[pde->dir], &pde->grids[2]);
    if (err)
      return err;
    pde->savetime = pde->last;
  }

  // Proceed from pde->last to new_last

  if (pde->last == 0)
    MakeGrid(pde);
  ExpandBounds(pde, pde->last, new_last);

  var = (double *)calloc(pde->nx, sizeof(double));
  mu = (double *)calloc(pde->nx, sizeof(double));
  r = (double *)calloc(pde->nx, sizeof(double));
  ad = (double *)calloc(pde->nx, sizeof(double));

  if (!var || !mu || !r || !ad) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }

  memset(r, 0, pde->nx * sizeof(double));
  lam = pde->pmdl->lambda;

  for (i = pde->last; i < new_last; i++) {
    dt = pde->times[i + 1] - pde->times[i];
    while (isig < pde->pmdl->nsig - 1 &&
           pde->pmdl->sigtms[isig] < pde->times[i] + 1e-5)
      isig++;
    sig = pde->pmdl->sig[isig];
    beta = pde->pmdl->beta[isig];

    std = sqrt(pde->cumxvar[i]);
    xcut_u = pde->cutcoef * pde->fr[i] * std;
    xcut_d = -pde->cutcoef * pde->fr[i] * std;

    // Init grid at the next period
    err = InitGridAt_i(pde, i + 1);
    if (err)
      goto FREE_RETURN;

    // Now calculate AD densities at the next period

    if (i == 0) // First step
    {
      // Calculate densities at time 1 (we have only one phi at time 1 !!!)
      std = sig * sqrt(dt);
      for (xi = 0; xi < pde->nx; xi++)
        pde->grids[!pde->dir].ad[0][xi] = gauss(pde->x[xi] / std);
    } else {
      // Calculate mu
      for (xi = 0; xi < pde->nx; xi++)
        mu[xi] = -lam * pde->x[xi] * dt;

      // Go phi slice by phi slice
      for (k = 0; k < pde->grids[pde->dir].nphi; k++) {
        phi_k = pde->grids[pde->dir].phi[k];

        // Calculate var
        for (xi = 0; xi < pde->nx; xi++) {
          std = CheyBeta_Vol(pde->x[xi], phi_k, sig, beta, pde->gamT[i],
                             pde->fr[i], xcut_u, xcut_d);
          var[xi] = std * std * dt;
        }

        // Convolve slice into the ad vector using Crank-Nicholson
        num_f_pde_one_step_forward(pde->tmppde, pde->nx, pde->x,
                                   pde->grids[pde->dir].ad[k], mu, var, r, 0.55,
                                   ad);

        // Share densities between next phis
        for (xi = 0; xi < pde->nx; xi++) {
          phi_nxt = phi_k + var[xi] - 2.0 * lam * phi_k * dt;

          if (phi_nxt < pde->grids[!pde->dir].phi[0])
            pde->grids[!pde->dir].ad[0][xi] += ad[xi];
          else if (phi_nxt >=
                   pde->grids[!pde->dir].phi[pde->grids[!pde->dir].nphi - 1])
            pde->grids[!pde->dir].ad[pde->grids[!pde->dir].nphi - 1][xi] +=
                ad[xi];
          else // find j: phi[j] <= phi_nxt < phi[j+1]
          {
            for (j = 0; pde->grids[!pde->dir].phi[j + 1] <= phi_nxt; j++)
              ;
            coef = (pde->grids[!pde->dir].phi[j + 1] - phi_nxt) /
                   (pde->grids[!pde->dir].phi[j + 1] -
                    pde->grids[!pde->dir].phi[j]);
            pde->grids[!pde->dir].ad[j][xi] += coef * ad[xi];
            pde->grids[!pde->dir].ad[j + 1][xi] += (1.0 - coef) * ad[xi];
          }
        }
      } // for (k=0; k < pde->grids[pde->dir].nphi; k++)
    }   // if (i == 0)

    // Renormalize densities
    sum = 0.0;
    for (k = 0; k < pde->grids[!pde->dir].nphi; k++)
      for (xi = 0; xi < pde->nx; xi++)
        sum += pde->grids[!pde->dir].ad[k][xi];

    // if (i==0) {	// Uncomment to switch renormalize off
    for (k = 0; k < pde->grids[!pde->dir].nphi; k++)
      for (xi = 0; xi < pde->nx; xi++)
        pde->grids[!pde->dir].ad[k][xi] /= sum;
    //}				// Uncomment to switch renormalize off
    // Switch grids
    FreeGrid(&pde->grids[pde->dir]);
    pde->dir = !pde->dir;

  } // for (i = pde->last; i < new_last; i++)
  pde->last = new_last;

FREE_RETURN:
  free(var);
  free(mu);
  free(r);
  free(ad);

  return err;
}

Err CheyBetaPDE_ResetFwd(SCheyBetaPDE *pde, double T) {
  Err err = NULL;
  int new_last;

  if (!pde->isfwd)
    return serror("The PDE was initialized as backward");

  // Find new last
  new_last = 0;
  while (new_last < pde->last && pde->times[new_last] < T - 1e-5)
    new_last++;
  if (new_last == pde->last)
    return NULL;

  if (pde->savetime == new_last) // Load grid
  {
    err = CopyGrid(&pde->grids[2], &pde->grids[pde->dir]);
    pde->last = new_last;
  } else {
    pde->last = pde->savetime = 0;
    FreeGrid(&pde->grids[pde->dir]);
    FreeGrid(&pde->grids[2]);
    err = CheyBetaPDE_ExpandFwd(pde, T);
  }
  return err;
}

#define MAXCPN 1000

static Err PayoffAt_i(SCheyBetaPDE *pde, int idx, int i, double ***pv) {
  Err err = NULL;
  int xi, j, k;
  long d1, d2;
  double x, phi, dffs[MAXCPN], gams[MAXCPN], dfs[MAXCPN], dff, gam, num,
      *dfs_ = dfs;
  SMktData mkt_data = {&dfs_, MKT_IRSTANDARD, NULL};

  // Precalculate dffs and gams
  d1 = (long)(pde->pmdl->today + pde->times[i] * DAYS_IN_YEAR + 1e-5);
  for (k = 0; k < pde->g->nmat[0][idx]; k++) {
    d2 =
        (long)(pde->pmdl->today + pde->g->mat[0][idx][k] * DAYS_IN_YEAR + 1e-5);
    dffs[k] = swp_f_df(d1, d2, pde->pmdl->ycname);
    gams[k] = (1.0 - exp(-pde->pmdl->lambda *
                         (pde->g->mat[0][idx][k] - pde->times[i]))) /
              pde->pmdl->lambda;
  }
  dff = pde->df[pde->nt - 1] / pde->df[i];
  gam = pde->gamT[i];

  // Calculate payoff at all nodes
  for (j = 0; j < pde->grids[pde->dir].nphi; j++) {
    phi = pde->grids[pde->dir].phi[j];
    for (xi = 0; xi < pde->nx; xi++) {
      x = pde->x[xi] - phi * gam;
      num = dff * exp(-gam * (x + 0.5 * gam * phi));

      for (k = 0; k < pde->g->nmat[0][idx]; k++)
        dfs[k] = dffs[k] * exp(-gams[k] * (x + 0.5 * gams[k] * phi));

      for (k = 0; k < pde->g->ninst; k++)
        pv[j][xi][k] *= num;

      err = pde->g->Payoff(pde->g, idx, pde->times[i], &mkt_data, pv[j][xi]);
      if (err)
        return err;

      for (k = 0; k < pde->g->ninst; k++)
        pv[j][xi][k] /= num;
    }
  }
  return NULL;
}

Err CheyBetaPDE_PriceSwaption(SCheyBetaPDE *pde, int iswp, double *result) {
  Err err = NULL;
  int xi, j;
  double ***payoff = NULL;

  if (!pde->isfwd)
    return serror("The PDE was initialized as backward");
  if (pde->g->type != PRODUCT_SWAPTIONS)
    return serror("Product is not a swaption");

  if (pde->times[pde->last] > pde->g->ex[iswp] + 1e-5)
    err = CheyBetaPDE_ResetFwd(pde, pde->g->ex[iswp]);
  else
    err = CheyBetaPDE_ExpandFwd(pde, pde->g->ex[iswp]);
  if (err)
    return err;

  payoff = f3tensor(0, pde->grids[pde->dir].nphi - 1, 0, pde->nx - 1, 0,
                    pde->g->ninst - 1);
  if (!payoff)
    return serror("Memory failure");

  memset(&(payoff[0][0][0]), 0,
         pde->grids[pde->dir].nphi * pde->nx * pde->g->ninst * sizeof(double));

  err = PayoffAt_i(pde, iswp, pde->last, payoff);
  if (err)
    goto FREE_RETURN;

  *result = 0.0;
  for (j = 0; j < pde->grids[pde->dir].nphi; j++)
    for (xi = 0; xi < pde->nx; xi++)
      *result += payoff[j][xi][iswp] * pde->grids[pde->dir].ad[j][xi];

  *result *= pde->df[pde->nt - 1];

FREE_RETURN:
  if (payoff)
    free_f3tensor(payoff, 0, pde->grids[pde->dir].nphi - 1, 0, pde->nx - 1, 0,
                  pde->g->ninst - 1);
  return err;
}

Err CheyBetaPDE_PriceBkwd(SCheyBetaPDE *pde, double *pv) {
  Err err = NULL;
  double *var = NULL, *mu = NULL, *r = NULL, **pvinst = NULL;
  int i, isig, iex, xi, j, k, m;
  double sig, beta, lam, dt, std, xcut_u, xcut_d;
  double phi_k, phi_nxt, coef;

  if (pde->isfwd)
    return serror("The PDE was initialized as forward");

  MakeGrid(pde);
  ExpandBounds(pde, 0, pde->nt - 1);

  var = (double *)calloc(pde->nx, sizeof(double));
  mu = (double *)calloc(pde->nx, sizeof(double));
  r = (double *)calloc(pde->nx, sizeof(double));
  pvinst = dmatrix(0, pde->nx - 1, 0, pde->g->ninst - 1);

  if (!var || !mu || !r || !pvinst) {
    err = serror("Memory failure");
    goto FREE_RETURN;
  }

  memset(r, 0, pde->nx * sizeof(double));
  lam = pde->pmdl->lambda;

  // Evaluate final payoff

  // Check that the final event is at the final time step
  if (fabs(pde->g->ex[iex = pde->g->nex - 1] - pde->times[pde->nt - 1]) >
      1e-5) {
    err = serror("Last event and last time step do not coincide");
    goto FREE_RETURN;
  }

  // Init grid at the last step in grids[!dir]
  err = InitGridAt_i(pde, pde->nt - 1);
  if (err)
    goto FREE_RETURN;
  pde->dir = !pde->dir;

  // Calculate payoffs at all nodes in grids[dir]
  err = PayoffAt_i(pde, iex, pde->nt - 1, pde->grids[pde->dir].pv);
  if (err)
    goto FREE_RETURN;

  // Go backwards to zero
  isig = pde->pmdl->nsig - 1;
  for (i = pde->nt - 2; i >= 0; i--) {
    dt = pde->times[i + 1] - pde->times[i];
    while (isig > 0 && pde->pmdl->sigtms[isig - 1] > pde->times[i] + 1e-5)
      isig--;
    sig = pde->pmdl->sig[isig];
    beta = pde->pmdl->beta[isig];

    while (iex >= 0 && pde->g->ex[iex] > pde->times[i] + 1e-5)
      iex--;

    std = sqrt(pde->cumxvar[i]);
    xcut_u = pde->cutcoef * pde->fr[i] * std;
    xcut_d = -pde->cutcoef * pde->fr[i] * std;

    // Init grid at period i in grids[!dir]
    err = InitGridAt_i(pde, i);
    if (err)
      goto FREE_RETURN;
    pde->dir = !pde->dir;

    // Calculate mu
    for (xi = 0; xi < pde->nx; xi++)
      mu[xi] = -lam * pde->x[xi] * dt;

    // Go phi slice by phi slice
    for (k = 0; k < pde->grids[pde->dir].nphi; k++) {
      phi_k = pde->grids[pde->dir].phi[k];

      // Calculate var
      for (xi = 0; xi < pde->nx; xi++) {
        std = CheyBeta_Vol(pde->x[xi], phi_k, sig, beta, pde->gamT[i],
                           pde->fr[i], xcut_u, xcut_d);
        var[xi] = std * std * dt;
      }

      // Convolve explicitly in phi into pvinst matrix
      for (xi = 0; xi < pde->nx; xi++) {
        phi_nxt = phi_k + var[xi] - 2.0 * lam * phi_k * dt;

        if (phi_nxt < pde->grids[!pde->dir].phi[0])
          memcpy(pvinst[xi], pde->grids[!pde->dir].pv[0][xi],
                 pde->g->ninst * sizeof(double));
        else if (phi_nxt >=
                 pde->grids[!pde->dir].phi[pde->grids[!pde->dir].nphi - 1])
          memcpy(pvinst[xi],
                 pde->grids[!pde->dir].pv[pde->grids[!pde->dir].nphi - 1][xi],
                 pde->g->ninst * sizeof(double));
        else // find j: phi[j] <= phi_nxt < phi[j+1]
        {
          for (j = 0; pde->grids[!pde->dir].phi[j + 1] <= phi_nxt; j++)
            ;
          coef =
              (pde->grids[!pde->dir].phi[j + 1] - phi_nxt) /
              (pde->grids[!pde->dir].phi[j + 1] - pde->grids[!pde->dir].phi[j]);

          for (m = 0; m < pde->g->ninst; m++)
            pvinst[xi][m] =
                coef * pde->grids[!pde->dir].pv[j][xi][m] +
                (1.0 - coef) * pde->grids[!pde->dir].pv[j + 1][xi][m];
        }
      }

      // Convolve slice using Crank-Nicholson
      num_f_pde_one_step_backward(pde->tmppde, pde->nx, pde->x, pde->g->ninst,
                                  pvinst, mu, var, r, 0.55,
                                  pde->grids[pde->dir].pv[k]);

    } // for (k=0; k < pde->grids[pde->dir].nphi; k++)

    // Evaluate events if any
    if (iex >= 0 && ((fabs(pde->g->ex[iex] - pde->times[i]) < 1e-5) ||
                     (pde->g->am && pde->g->am[iex]))) {
      err = PayoffAt_i(pde, iex, i, pde->grids[pde->dir].pv);
      if (err)
        goto FREE_RETURN;
    }
    FreeGrid(&pde->grids[!pde->dir]);
  }

  // Return results
  for (m = 0; m < pde->g->ninst; m++)
    pv[m] = pde->grids[pde->dir].pv[0][pde->ix0][m] * pde->df[pde->nt - 1];

FREE_RETURN:
  free(var);
  free(mu);
  free(r);
  if (pvinst)
    free_dmatrix(pvinst, 0, pde->nx - 1, 0, pde->g->ninst - 1);

  return err;
}
