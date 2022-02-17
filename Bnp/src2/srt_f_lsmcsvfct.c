/* ----------------------------------------------------------------------------------
   AUTHOR: E. FOURNIE

   DATE : SEPTEMBRE 99

   FILENAME: srt_f_lsmcsvfct.c

   PURPOSE:  utility functions for the implementation of LSM algo
                      for the Cheyette beta stoch vol

   MODIFICATION:

   ----------------------------------------------------------------------------------
 */
#include "math.h"
#include "num_h_gaussj.h"
#include "srt_h_all.h"

#include "srt_h_lsmcsvfct.h"

static Err cheybetastochvol_evolve_in_Grid(double **rndm_mat, SrtStpPtr step,
                                           MCTreePointCheyBetaStochVol **grid,
                                           long n) {
  Err err = NULL;
  SrtStpPtr top;
  long t;
  SrtIRMTmInf *tminf1;
  double dt, sqdt, c1, c2, dw1, dw2, sig, eta, rho, lam, r, x, phi, beta, vol,
      vovol, sqrvovol;

  /*--- Start from the first time step ---*/
  step = top = gototop(step);
  t = 0;
  tminf1 = step->tminf[0];

  /* Initialization for the first time step (today) */
  x = grid[0][n].x = 0.0;
  eta = grid[0][n].sigma = 1.0;
  phi = grid[0][n].phi = 0.0;
  r = grid[0][n].short_rate = sam_get(tminf1->fwd_sam, 0, F_0_t);
  grid[0][n].df = exp(-grid[0][n].short_rate * step->delta_t);
  grid[0][n].numeraire = 1.0;

  /* Loops on all the time steps */
  while (step->next != NULL) {
    dt = step->delta_t;
    sqdt = sqrt(dt);
    lam = tminf1->ev.onef.lambda;
    sig = tminf1->ev.onef.sig;
    beta = tminf1->ev.onef.beta;
    rho = tminf1->rho;
    vovol = tminf1->vovol;
    sqrvovol = vovol * vovol;

    /* moves to the next time step */
    step = step->next;
    t++;
    tminf1 = step->tminf[0];

    /* Compute the correlated random numbers */
    c1 = sqrt((1.0 + rho) / 2.0);
    c2 = sqrt((1.0 - rho) / 2.0);
    dw1 = (c1 * rndm_mat[0][t] + c2 * rndm_mat[1][t]) * sqdt;
    dw2 = (c1 * rndm_mat[0][t] - c2 * rndm_mat[1][t]) * sqdt;

    /* Populate the next time step with the new values */
    vol = sig * eta * pow(fabs(r), beta);
    x = grid[t][n].x = x + (phi - lam * x) * dt + DMIN(vol, 0.3) * dw1;
    phi = grid[t][n].phi = phi + (vol * vol - 2.0 * lam * phi) * dt;
    eta = grid[t][n].sigma = eta * exp(-sqrvovol * dt + vovol * dw2);

    r = grid[t][n].short_rate =
        sam_get(tminf1->fwd_sam, 0, F_0_t) + grid[t][n].x;
    grid[t][n].df = exp(-grid[t][n].short_rate * dt);
    grid[t][n].numeraire = grid[t - 1][n].numeraire / grid[t - 1][n].df;

  } /* END while step loop */

  return NULL;
}

/* ----------------------------------------------------------------------------------
 */

/* ---------------------------------------------------------------------------------
   FUNCTION:     Err srt_f_cheybeta1f1equtree_make_grid()
   PURPOSE:      Generate the McGrid
   ---------------------------------------------------------------------------------
 */
Err cheybetastochvol_simul(SrtGrfnParam *grfnparam, SrtStpPtr step,
                           SrtUndInfo *und_info,
                           MCTreePointCheyBetaStochVol **grid, long nt) {
  Err err = NULL;
  SrtUndPtr und = NULL;
  SrtStpPtr top = gototop(step);
  long n;
  SrtMCRanSrc mcrs;

  /*	Initialises mcrs according to inputs  , makes memory allocation or
          even generates fully the non correlated random numbers */
  err = SrtMCRanSrc_start(&mcrs, 0, grfnparam->num_MCarlo_paths - 1, 1, nt,
                          und_info->no_of_brownians, grfnparam->sample_type,
                          grfnparam->rand_seed, top);
  if (err)
    return err;

  /* ----------------- PATHS SIMULATION ------------------------ */
  /* Loop on all the paths used to generate the initial draw */
  for (n = 0; n < grfnparam->num_MCarlo_paths; n++) {
    /* Generates or gets the next MC sample path in mcrs.r[mcrs.cur_path] */
    err = SrtMCRanSrc_next(&mcrs);
    if (err) {
      SrtMCRanSrc_end(&mcrs);
      return err;
    }

    /* Generate the current path for the Cheybeta underlying in the Grid */
    err = cheybetastochvol_evolve_in_Grid(mcrs.r[mcrs.cur_path], top, grid, n);
    if (err) {
      SrtMCRanSrc_end(&mcrs);
      return err;
    }
  }

  /* Frees the random number for the MC simulation */
  SrtMCRanSrc_end(&mcrs);

  return NULL;
}

/* ---------------------------------------------------------------------------------------------
 */
/* ---------------------------------------------------------------------------------------------
 */

typedef double (*funcbasiscsv)(double, double, double);
static double **cov, **vec, *y;
funcbasiscsv *phi;

/* ---------------------------------------------------------------------------------------------
 */
static double phi1(double x, double phi, double sigma) { return (1.0); }

static double phi2(double x, double phi, double sigma) { return (x); }

static double phi3(double x, double phi, double sigma) { return (x * x); }

static double phi4(double x, double phi, double sigma) { return (x * x * x); }

static double phi5(double x, double phi, double sigma) {
  return (x * x * x * x * x);
}

static double phi6(double x, double phi, double sigma) {
  return (x * x * x * x * x * x);
}

static double phi7(double x, double phi, double sigma) {
  return (x * x * x * x * x * x * x);
}

static double phi8(double x, double phi, double sigma) { return (sigma); }

static double phi9(double x, double phi, double sigma) {
  return (sigma * sigma);
}

static double phi10(double x, double phi, double sigma) { return (phi); }

static double phi11(double x, double phi, double sigma) {

  return (sigma * sigma * sigma);
}

static double phi12(double x, double phi, double sigma) { return (x * sigma); }

static double phi13(double x, double phi, double sigma) { return (x * phi); }

/*
static double phi11(double x  , double phi  , double sigma)
{

        return( phi*phi );
}

static double phi12(double x  , double phi  , double sigma)
{
        return( x*sigma);

}

static double phi13(double x  , double phi  , double sigma)
{
        return( x*phi );
}
*/
/* ---------------------------------------------------------------------------------------------
 */

Err srt_f_alloc_lsmcsv_resources(int nreg, int node_dim) {
  /*--- REGR ---*/

  cov = dmatrix(1, nreg, 1, nreg);
  vec = dmatrix(1, nreg, 1, node_dim);
  y = dvector(1, nreg);

  phi = (funcbasiscsv *)malloc(nreg * sizeof(funcbasiscsv));
  phi[0] = phi1;
  phi[1] = phi2;
  phi[2] = phi3;
  phi[3] = phi4;
  phi[4] = phi5;
  phi[5] = phi6;
  phi[6] = phi7;
  phi[7] = phi8;
  phi[8] = phi9;
  phi[9] = phi10;
  phi[10] = phi11;
  phi[11] = phi12;
  phi[12] = phi13;

  return NULL;
}

Err srt_f_free_lsmcsv_resources(int nreg, int node_dim) {
  free_dmatrix(cov, 1, nreg, 1, nreg);
  free_dmatrix(vec, 1, nreg, 1, node_dim);
  free_dvector(y, 1, nreg);
  free(phi);

  return NULL;
}

/* --------------------------------------------------------- */
/* solve the quadratic minimization problem by direct method */
/* --------------------------------------------------------- */

Err srt_f_lsmcsv_regr(int nreg, int node_dim, int ntraj, int nt,
                      MCTreePointCheyBetaStochVol **grid, double ***bwd_assets,
                      double m1, double m2, double s11, double s12, double s21,
                      double s22) {
  Err err = NULL;
  int i, k, l;
  double x, sigma, psi;

  /*--- initialisation ---*/
  for (k = 1; k <= nreg; k++) {
    for (l = 1; l <= nreg; l++)
      cov[k][l] = 0.0;
    for (l = 1; l <= node_dim; l++)
      vec[k][l] = 0.0;
  }

  for (i = 0; i < ntraj; i++) {
    x = s11 * (grid[nt][i].x - m1) + s12 * (grid[nt][i].sigma - m2);
    sigma = s21 * (grid[nt][i].x - m1) + s22 * (grid[nt][i].sigma - m2);
    psi = grid[nt][i].phi;

    for (k = 1; k <= nreg; k++)
      y[k] = (phi[k - 1])(x, psi, sigma);

    for (k = 1; k <= nreg; k++) {
      for (l = 1; l <= node_dim; l++)
        vec[k][l] += y[k] * bwd_assets[nt + 1][i][l - 1] * grid[nt][i].df;
      for (l = 1; l <= nreg; l++)
        cov[k][l] += y[k] * y[l];
    }
  }

  gaussj(cov, nreg, vec, node_dim);

  /* approximated value at t */
  for (i = 0; i < ntraj; i++) {
    x = s11 * (grid[nt][i].x - m1) + s12 * (grid[nt][i].sigma - m2);
    sigma = s21 * (grid[nt][i].x - m1) + s22 * (grid[nt][i].sigma - m2);
    psi = grid[nt][i].phi;

    for (l = 1; l <= node_dim; l++) {
      bwd_assets[nt][i][l - 1] = 0.0;
      for (k = 1; k <= nreg; k++) {
        y[k] = (phi[k - 1])(x, psi, sigma);
        bwd_assets[nt][i][l - 1] += vec[k][l] * y[k];
      }
    }
  }

  return err;
}

/* ----------------------------------------------------------------------------------*/
/*--------------------------------- End of File
 * -------------------------------------*/
