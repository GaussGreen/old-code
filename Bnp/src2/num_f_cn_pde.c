/* -------------------------------------------------------------------------
   FILENAME		: num_f_cn_pde.c
   AUTHOR		: Toto 15Dec1999

   PURPOSE		: convolve backward/forward using Cranck-Nicholson
   equations
   ------------------------------------------------------------------------- */

#include "math.h"
#include "num_h_allhdr.h"

/*	Calculate transition probabilities	*/
static void calc_trans_prob(
    CNPDE_TEMP *tmp,
    /*	Size of the array	*/
    int n,
    /*	X values	*/
    double *x,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double *mu,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double *var,
    /*	Discount rates from t to t + 1	*/
    double *r,
    /*	Cranck-Nicholson theta	*/
    double theta,
    /*	Backward (-1) or forward (+1)	*/
    int for_back)
/*	Result: probabilities to go up  , central and down	*/
{
  int i;
  double du, dd;

  for (i = 1; i < n - 1; i++) {
    du = x[i + 1] - x[i];
    dd = x[i] - x[i - 1];

    if (for_back == -1)
    /*	Backward	*/
    {
      tmp->exp_pu[i] = (var[i] + mu[i] * dd) / (du * (du + dd));
      tmp->exp_pd[i] = (tmp->exp_pu[i] * du - mu[i]) / dd;
      tmp->exp_pc[i] = -tmp->exp_pu[i] - tmp->exp_pd[i];
    } else
    /*	Forward	*/
    {
      tmp->exp_pd[i] = (mu[i - 1] * du + var[i - 1]) / (dd * (du + dd));
      tmp->exp_pc[i] = (-mu[i] * (du - dd) - var[i]) / du / dd;
      tmp->exp_pu[i] = (-mu[i + 1] * dd + var[i + 1]) / (du * (du + dd));
    }
  }

  du = x[1] - x[0];
  tmp->exp_pd[0] = 0.0;
  if (for_back == -1)
  /*	Backward	*/
  {
    tmp->exp_pc[0] = -mu[0] / du;
    tmp->exp_pu[0] = mu[0] / du;
  } else
  /*	Forward	*/
  {
    tmp->exp_pc[0] = mu[0] / du;
    tmp->exp_pu[0] = -mu[1] / du;
  }

  dd = x[n - 1] - x[n - 2];
  tmp->exp_pu[n - 1] = 0.0;
  if (for_back == -1)
  /*	Backward	*/
  {
    tmp->exp_pc[n - 1] = mu[n - 1] / dd;
    tmp->exp_pd[n - 1] = -mu[n - 1] / dd;
  } else
  /*	Forward	*/
  {
    tmp->exp_pc[n - 1] = -mu[n - 1] / dd;
    tmp->exp_pd[n - 1] = mu[n - 2] / dd;
  }

  for (i = 0; i < n; i++) {
    tmp->imp_pu[i] = -theta * tmp->exp_pu[i];
    tmp->imp_pd[i] = -theta * tmp->exp_pd[i];
    tmp->imp_pc[i] = -theta * tmp->exp_pc[i];

    tmp->exp_pu[i] *= (1.0 - theta);
    tmp->exp_pd[i] *= (1.0 - theta);
    tmp->exp_pc[i] *= (1.0 - theta);

    tmp->exp_pc[i] += 1.0 - (1.0 - theta) * r[i];
    tmp->imp_pc[i] += 1.0 + theta * r[i];
  }
}

/*	Calculate right member of FDM equation	*/
/*	1 function	*/
static void calc_right_mem_1(CNPDE_TEMP *tmp,
                             /*	Size of the array	*/
                             int n,
                             /*	f values	*/
                             /*	i: state  , j: function number */
                             double *f) {
  int i;

  for (i = 1; i < n - 1; i++) {
    tmp->right_mem[i][0] = tmp->exp_pu[i] * f[i + 1] + tmp->exp_pc[i] * f[i] +
                           tmp->exp_pd[i] * f[i - 1];
  }

  tmp->right_mem[0][0] = tmp->exp_pu[0] * f[1] + tmp->exp_pc[0] * f[0];

  tmp->right_mem[n - 1][0] =
      tmp->exp_pc[n - 1] * f[n - 1] + tmp->exp_pd[n - 1] * f[n - 2];
}

/*	Calculate right member of FDM equation	*/
/*	m functions	*/
static void calc_right_mem_m(CNPDE_TEMP *tmp,
                             /*	Size of the array	*/
                             int n,
                             /*	Number of functions	*/
                             int m,
                             /*	f values	*/
                             /*	i: state  , j: function number */
                             double **f) {
  int i, j;

  for (j = 0; j < m; j++) {
    for (i = 1; i < n - 1; i++) {
      tmp->right_mem[i][j] = tmp->exp_pu[i] * f[i + 1][j] +
                             tmp->exp_pc[i] * f[i][j] +
                             tmp->exp_pd[i] * f[i - 1][j];
    }

    tmp->right_mem[0][j] = tmp->exp_pu[0] * f[1][j] + tmp->exp_pc[0] * f[0][j];

    tmp->right_mem[n - 1][j] =
        tmp->exp_pc[n - 1] * f[n - 1][j] + tmp->exp_pd[n - 1] * f[n - 2][j];
  }
}

/*	Solve FDM equation	*/
/*	1 function	*/
static void calc_solve_fdm_1(CNPDE_TEMP *tmp,
                             /*	Size of the array	*/
                             int n,
                             /*	Result: f values	*/
                             /*	i: state  , j: function number */
                             double *f) {
  int j;
  double bet;

  f[0] = tmp->right_mem[0][0] / (bet = tmp->imp_pc[0]);

  for (j = 1; j < n; j++) {
    tmp->gam[j] = tmp->imp_pu[j - 1] / bet;
    bet = tmp->imp_pc[j] - tmp->imp_pd[j] * tmp->gam[j];
    f[j] = (tmp->right_mem[j][0] - tmp->imp_pd[j] * f[j - 1]) / bet;
  }

  for (j = n - 2; j >= 0; j--) {
    f[j] -= tmp->gam[j + 1] * f[j + 1];
  }
}

/*	Solve FDM equation	*/
/*	m functions	*/
static void calc_solve_fdm_m(CNPDE_TEMP *tmp,
                             /*	Size of the array	*/
                             int n,
                             /*	Number of functions	*/
                             int m,
                             /*	Result: f values	*/
                             /*	i: state  , j: function number */
                             double **f) {
  int j, k;
  double bet;

  for (k = 0; k < m; k++) {
    f[0][k] = tmp->right_mem[0][k] / (bet = tmp->imp_pc[0]);
  }

  for (j = 1; j < n; j++) {
    tmp->gam[j] = tmp->imp_pu[j - 1] / bet;
    bet = tmp->imp_pc[j] - tmp->imp_pd[j] * tmp->gam[j];
    for (k = 0; k < m; k++) {
      f[j][k] = (tmp->right_mem[j][k] - tmp->imp_pd[j] * f[j - 1][k]) / bet;
    }
  }

  for (j = n - 2; j >= 0; j--) {
    for (k = 0; k < m; k++) {
      f[j][k] -= tmp->gam[j + 1] * f[j + 1][k];
    }
  }
}

/*	Call before PDE	*/
void num_f_pde_init(CNPDE_TEMP *tmp,
                    /*	Size of the array	*/
                    int n,
                    /*	Number of functions	*/
                    int m) {
  tmp->exp_pu = (double *)calloc(n, sizeof(double));
  tmp->exp_pc = (double *)calloc(n, sizeof(double));
  tmp->exp_pd = (double *)calloc(n, sizeof(double));

  tmp->imp_pu = (double *)calloc(n, sizeof(double));
  tmp->imp_pc = (double *)calloc(n, sizeof(double));
  tmp->imp_pd = (double *)calloc(n, sizeof(double));

  tmp->gam = (double *)calloc(n, sizeof(double));

  tmp->right_mem = dmatrix(0, n - 1, 0, m - 1);
}

/*	Call after PDE	*/
void num_f_pde_free(CNPDE_TEMP *tmp,
                    /*	Size of the array	*/
                    int n,
                    /*	Number of functions	*/
                    int m) {
  free(tmp->exp_pu);
  free(tmp->exp_pc);
  free(tmp->exp_pd);

  free(tmp->imp_pu);
  free(tmp->imp_pc);
  free(tmp->imp_pd);

  free(tmp->gam);

  free_dmatrix(tmp->right_mem, 0, n - 1, 0, m - 1);
}

/*	Convolve one step backward using Kolmogorov backward equation	*/
void num_f_pde_one_step_backward(
    CNPDE_TEMP *tmp,
    /*	Size of the array	*/
    int n,
    /*	X values	*/
    double *x,
    /*	Number of functions	*/
    int m,
    /*	f values @ t+1	*/
    /*	i: state (1..n)  , j:function number (1..m)	*/
    double **f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double *mu,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double *var,
    /*	Discount rates from t to t + 1	*/
    double *r,
    /*	Cranck-Nicholson theta	*/
    double theta,
    /*	Result: f values @ t	*/
    double **f_t) {
  calc_trans_prob(tmp, n, x, mu, var, r, theta, -1);
  calc_right_mem_m(tmp, n, m, f_t_plus_1);
  calc_solve_fdm_m(tmp, n, m, f_t);
}

/*	De-convolve one step forward using Kolmogorov forward equation	*/
void num_f_pde_one_step_forward(
    CNPDE_TEMP *tmp,
    /*	Size of the array	*/
    int n,
    /*	X values	*/
    double *x,
    /*	f values @ t	*/
    double *f_t,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double *mu,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double *var,
    /*	Discount rates from t to t + 1	*/
    double *r,
    /*	Cranck-Nicholson theta	*/
    double theta,
    /*	Result: f values @ t+1	*/
    double *f_t_plus_1) {
  calc_trans_prob(tmp, n, x, mu, var, r, theta, 1);
  calc_right_mem_1(tmp, n, f_t);
  calc_solve_fdm_1(tmp, n, f_t_plus_1);
}

/*	2D Scheme */

/*	Call before PDE	2D */
void num_f_pde_init_2d_lod(CNPDE_TEMP_2D_LOD *tmp,
                           /*	Size of the arrays	*/
                           int nx, int ny,
                           /*	Number of functions	*/
                           int m) {
  tmp->tmpx = malloc(sizeof(CNPDE_TEMP));
  tmp->tmpy = malloc(sizeof(CNPDE_TEMP));
  num_f_pde_init(tmp->tmpx, nx, m);
  num_f_pde_init(tmp->tmpy, ny, m);
}

/*	Call after PDE 2D */
void num_f_pde_free_2d_lod(CNPDE_TEMP_2D_LOD *tmp,
                           /*	Size of the arrays	*/
                           int nx, int ny,
                           /*	Number of functions	*/
                           int m) {
  num_f_pde_free(tmp->tmpx, nx, m);
  num_f_pde_free(tmp->tmpy, ny, m);

  free(tmp->tmpx);
  free(tmp->tmpy);
}

/*	Convolve one step backward using LOD scheme */
void num_f_pde_one_step_backward_2f_lod(
    CNPDE_TEMP_2D_LOD *tmp,
    /*	Size of the X array	*/
    int nx,
    /*	X values	*/
    double *x,
    /*	Size of the Y array	*/
    int ny,
    /*	X values	*/
    double *y,
    /*	Number of functions	*/
    int m,
    /*	f values @ t+1	*/
    /*	i: x state  , j: y state  , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Cranck-Nicholson theta	*/
    double theta,
    /*	Result: f values @ t	*/
    double ***f_t, int lx, int ux, int ly, int uy) {
  static int i, j, k;
  static double du, dd, bet;
  static double *muyi, *varyi, *ri, **fi;

  CNPDE_TEMP *tmpx = tmp->tmpx, *tmpy = tmp->tmpy;

  for (j = ly; j <= uy; j++) {
    /*	For each slice in y  , do pde in x */

    /*	calc trans prob */

    for (i = lx + 1; i <= ux - 1; i++) {
      du = x[i + 1] - x[i];
      dd = x[i] - x[i - 1];

      tmpx->exp_pu[i] = (varx[i][j] + mux[i][j] * dd) / (du * (du + dd));
      tmpx->exp_pd[i] = (tmpx->exp_pu[i] * du - mux[i][j]) / dd;
      tmpx->exp_pc[i] = -tmpx->exp_pu[i] - tmpx->exp_pd[i];
    }

    du = x[lx + 1] - x[lx];
    tmpx->exp_pd[lx] = 0.0;
    tmpx->exp_pc[lx] = -mux[lx][j] / du;
    tmpx->exp_pu[lx] = -tmpx->exp_pc[lx];

    dd = x[ux] - x[ux - 1];
    tmpx->exp_pu[ux] = 0.0;
    tmpx->exp_pc[ux] = mux[ux][j] / dd;
    tmpx->exp_pd[ux] = -tmpx->exp_pc[ux];

    for (i = lx; i <= ux; i++) {
      tmpx->imp_pu[i] = -theta * tmpx->exp_pu[i];
      tmpx->imp_pd[i] = -theta * tmpx->exp_pd[i];
      tmpx->imp_pc[i] = -theta * tmpx->exp_pc[i];

      tmpx->exp_pu[i] *= (1.0 - theta);
      tmpx->exp_pd[i] *= (1.0 - theta);
      tmpx->exp_pc[i] *= (1.0 - theta);

      tmpx->exp_pc[i] += 1.0 - (1.0 - theta) * r[i][j] * 0.50;
      tmpx->imp_pc[i] += 1.0 + theta * r[i][j] * 0.50;
    }

    /*	calc right mem */

    for (k = 0; k < m; k++) {
      for (i = lx + 1; i <= ux - 1; i++) {
        tmpx->right_mem[i][k] = tmpx->exp_pu[i] * f_t_plus_1[i + 1][j][k] +
                                tmpx->exp_pc[i] * f_t_plus_1[i][j][k] +
                                tmpx->exp_pd[i] * f_t_plus_1[i - 1][j][k];
      }

      tmpx->right_mem[lx][k] = tmpx->exp_pu[lx] * f_t_plus_1[lx + 1][j][k] +
                               tmpx->exp_pc[lx] * f_t_plus_1[lx][j][k];

      tmpx->right_mem[ux][k] = tmpx->exp_pc[ux] * f_t_plus_1[ux][j][k] +
                               tmpx->exp_pd[ux] * f_t_plus_1[ux - 1][j][k];
    }

    /*	solve fdm */

    for (k = 0; k < m; k++) {
      f_t[lx][j][k] = tmpx->right_mem[lx][k] / (bet = tmpx->imp_pc[lx]);
    }

    for (i = lx + 1; i <= ux; i++) {
      tmpx->gam[i] = tmpx->imp_pu[i - 1] / bet;
      bet = tmpx->imp_pc[i] - tmpx->imp_pd[i] * tmpx->gam[i];
      for (k = 0; k < m; k++) {
        f_t[i][j][k] =
            (tmpx->right_mem[i][k] - tmpx->imp_pd[i] * f_t[i - 1][j][k]) / bet;
      }
    }

    for (i = ux - 1; i >= lx; i--) {
      for (k = 0; k < m; k++) {
        f_t[i][j][k] -= tmpx->gam[i + 1] * f_t[i + 1][j][k];
      }
    }
  }

  for (i = lx; i <= ux; i++) {
    /*	For each slice in x  , do pde in y */
    muyi = muy[i];
    varyi = vary[i];
    ri = r[i];
    fi = f_t[i];

    /*	calc trans prob */

    for (j = ly + 1; j <= uy - 1; j++) {
      du = y[j + 1] - y[j];
      dd = y[j] - y[j - 1];

      tmpy->exp_pu[j] = (varyi[j] + muyi[j] * dd) / (du * (du + dd));
      tmpy->exp_pd[j] = (tmpy->exp_pu[j] * du - muyi[j]) / dd;
      tmpy->exp_pc[j] = -tmpy->exp_pu[j] - tmpy->exp_pd[j];
    }

    du = y[ly + 1] - y[ly];
    tmpy->exp_pd[ly] = 0.0;
    tmpy->exp_pc[ly] = -muyi[ly] / du;
    tmpy->exp_pu[ly] = -tmpy->exp_pc[ly];

    dd = y[uy] - y[uy - 1];
    tmpy->exp_pu[uy] = 0.0;
    tmpy->exp_pc[uy] = muyi[uy] / dd;
    tmpy->exp_pd[uy] = -tmpy->exp_pc[uy];

    for (j = ly; j <= uy; j++) {
      tmpy->imp_pu[j] = -theta * tmpy->exp_pu[j];
      tmpy->imp_pd[j] = -theta * tmpy->exp_pd[j];
      tmpy->imp_pc[j] = -theta * tmpy->exp_pc[j];

      tmpy->exp_pu[j] *= (1.0 - theta);
      tmpy->exp_pd[j] *= (1.0 - theta);
      tmpy->exp_pc[j] *= (1.0 - theta);

      tmpy->exp_pc[j] += 1.0 - (1.0 - theta) * ri[j] * 0.50;
      tmpy->imp_pc[j] += 1.0 + theta * ri[j] * 0.50;
    }

    /*	calc right mem */

    for (k = 0; k < m; k++) {
      for (j = ly + 1; j <= uy - 1; j++) {
        tmpy->right_mem[j][k] = tmpy->exp_pu[j] * fi[j + 1][k] +
                                tmpy->exp_pc[j] * fi[j][k] +
                                tmpy->exp_pd[j] * fi[j - 1][k];
      }

      tmpy->right_mem[ly][k] =
          tmpy->exp_pu[ly] * fi[ly + 1][k] + tmpy->exp_pc[ly] * fi[ly][k];

      tmpy->right_mem[uy][k] =
          tmpy->exp_pc[uy] * fi[uy][k] + tmpy->exp_pd[uy] * fi[uy - 1][k];
    }

    /*	solve fdm */

    for (k = 0; k < m; k++) {
      fi[ly][k] = tmpy->right_mem[ly][k] / (bet = tmpy->imp_pc[ly]);
    }

    for (j = ly + 1; j <= uy; j++) {
      tmpy->gam[j] = tmpy->imp_pu[j - 1] / bet;
      bet = tmpy->imp_pc[j] - tmpy->imp_pd[j] * tmpy->gam[j];
      for (k = 0; k < m; k++) {
        fi[j][k] =
            (tmpy->right_mem[j][k] - tmpy->imp_pd[j] * fi[j - 1][k]) / bet;
      }
    }

    for (j = uy - 1; j >= ly; j--) {
      for (k = 0; k < m; k++) {
        fi[j][k] -= tmpy->gam[j + 1] * fi[j + 1][k];
      }
    }
  }
}

/*	Call before PDE	2D */
void num_f_pde_init_2d_adi(CNPDE_TEMP_2D_ADI *tmp,
                           /*	Size of the arrays	*/
                           int nx, int ny,
                           /*	Number of functions	*/
                           int m) {
  tmp->xexp_pu = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->xexp_pc = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->xexp_pd = dmatrix(0, nx - 1, 0, ny - 1);

  tmp->ximp_pu = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->ximp_pc = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->ximp_pd = dmatrix(0, nx - 1, 0, ny - 1);

  tmp->yexp_pu = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->yexp_pc = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->yexp_pd = dmatrix(0, nx - 1, 0, ny - 1);

  tmp->yimp_pu = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->yimp_pc = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->yimp_pd = dmatrix(0, nx - 1, 0, ny - 1);

  tmp->f = f3tensor(0, nx - 1, 0, ny - 1, 0, m - 1);

  tmp->gamx = (double *)calloc(nx, sizeof(double));
  tmp->gamy = (double *)calloc(ny, sizeof(double));
}

/*	Call after PDE 2D */
void num_f_pde_free_2d_adi(CNPDE_TEMP_2D_ADI *tmp,
                           /*	Size of the arrays	*/
                           int nx, int ny,
                           /*	Number of functions	*/
                           int m) {
  free_dmatrix(tmp->xexp_pu, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->xexp_pc, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->xexp_pd, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->ximp_pu, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->ximp_pc, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->ximp_pd, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yexp_pu, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yexp_pc, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yexp_pd, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yimp_pu, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yimp_pc, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yimp_pd, 0, nx - 1, 0, ny - 1);

  free_f3tensor(tmp->f, 0, nx - 1, 0, ny - 1, 0, m - 1);

  free(tmp->gamx);
  free(tmp->gamy);
}

/*	Convolve one step backward using ADI scheme */
void num_f_pde_one_step_backward_2f_adi(
    CNPDE_TEMP_2D_ADI *tmp,
    /*	Size of the X array	*/
    int nx,
    /*	X values	*/
    double *x,
    /*	Size of the Y array	*/
    int ny,
    /*	X values	*/
    double *y,
    /*	Number of functions	*/
    int mstart, int mend,
    /*	f values @ t+1	*/
    /*	i: x state  , j: y state  , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy) {
  static int i, j, k;
  static double dux, ddx, duy, ddy, bet;
  static double *muxi, *varxi, *ri, **fi, **tmpfi, *xexp_pui, *xexp_pci,
      *xexp_pdi, *ximp_pui, *ximp_pci, *ximp_pdi, *yexp_pui, *yexp_pci,
      *yexp_pdi, *yimp_pui, *yimp_pci, *yimp_pdi;
  static double temp1, temp2;

  double **xexp_pu = tmp->xexp_pu, **xexp_pc = tmp->xexp_pc,
         **xexp_pd = tmp->xexp_pd,

         **ximp_pu = tmp->ximp_pu, **ximp_pc = tmp->ximp_pc,
         **ximp_pd = tmp->ximp_pd,

         **yexp_pu = tmp->yexp_pu, **yexp_pc = tmp->yexp_pc,
         **yexp_pd = tmp->yexp_pd,

         **yimp_pu = tmp->yimp_pu, **yimp_pc = tmp->yimp_pc,
         **yimp_pd = tmp->yimp_pd,

         ***tmpf = tmp->f,

         *gamx = tmp->gamx, *gamy = tmp->gamy;

  /*	Calculate all probabilities */

  /*	Interior probas on y */
  for (j = ly + 1; j <= uy - 1; j++) {
    duy = y[j + 1] - y[j];
    ddy = y[j] - y[j - 1];
    temp1 = 0.5 / (duy * (duy + ddy));

    for (i = lx; i <= ux; i++) {
      yexp_pu[i][j] = temp1 * (vary[i][j] + muy[i][j] * ddy);
      yexp_pd[i][j] = (yexp_pu[i][j] * duy - 0.5 * muy[i][j]) / ddy;
      yexp_pc[i][j] = -yexp_pu[i][j] - yexp_pd[i][j];

      yimp_pu[i][j] = -yexp_pu[i][j];
      yimp_pd[i][j] = -yexp_pd[i][j];
      yimp_pc[i][j] = -yexp_pc[i][j];

      temp2 = 0.25 * r[i][j];
      yexp_pc[i][j] += 1.0 - temp2;
      yimp_pc[i][j] += 1.0 + temp2;
    }
  }

  /*	Down probas on y */
  duy = y[ly + 1] - y[ly];

  for (i = lx; i <= ux; i++) {
    yexp_pd[i][ly] = 0.0;
    yexp_pc[i][ly] = -0.5 * muy[i][ly] / duy;
    yexp_pu[i][ly] = -yexp_pc[i][ly];

    yimp_pu[i][ly] = -yexp_pu[i][ly];
    yimp_pd[i][ly] = -yexp_pd[i][ly];
    yimp_pc[i][ly] = -yexp_pc[i][ly];

    temp2 = 0.25 * r[i][ly];
    yexp_pc[i][ly] += 1.0 - temp2;
    yimp_pc[i][ly] += 1.0 + temp2;
  }

  /*	Up probas on y */
  ddy = y[uy] - y[uy - 1];

  for (i = lx; i <= ux; i++) {
    yexp_pu[i][uy] = 0.0;
    yexp_pc[i][uy] = 0.5 * muy[i][uy] / ddy;
    yexp_pd[i][uy] = -yexp_pc[i][uy];

    yimp_pu[i][uy] = -yexp_pu[i][uy];
    yimp_pd[i][uy] = -yexp_pd[i][uy];
    yimp_pc[i][uy] = -yexp_pc[i][uy];

    temp2 = 0.25 * r[i][uy];
    yexp_pc[i][uy] += 1.0 - temp2;
    yimp_pc[i][uy] += 1.0 + temp2;
  }

  /*	Interior probas on x */
  for (i = lx + 1; i <= ux - 1; i++) {
    dux = x[i + 1] - x[i];
    ddx = x[i] - x[i - 1];
    temp1 = 0.5 / (dux * (dux + ddx));

    muxi = mux[i];
    varxi = varx[i];
    ri = r[i];
    xexp_pui = xexp_pu[i];
    xexp_pci = xexp_pc[i];
    xexp_pdi = xexp_pd[i];
    ximp_pui = ximp_pu[i];
    ximp_pci = ximp_pc[i];
    ximp_pdi = ximp_pd[i];

    for (j = ly; j <= uy; j++) {
      xexp_pui[j] = temp1 * (varxi[j] + muxi[j] * ddx);
      xexp_pdi[j] = (xexp_pui[j] * dux - 0.5 * muxi[j]) / ddx;
      xexp_pci[j] = -xexp_pui[j] - xexp_pdi[j];

      ximp_pui[j] = -xexp_pui[j];
      ximp_pdi[j] = -xexp_pdi[j];
      ximp_pci[j] = -xexp_pci[j];

      temp2 = 0.25 * ri[j];
      xexp_pci[j] += 1.0 - temp2;
      ximp_pci[j] += 1.0 + temp2;
    }
  }

  /*	Down probas on x */
  dux = x[lx + 1] - x[lx];

  muxi = mux[lx];
  varxi = varx[lx];
  ri = r[lx];
  xexp_pui = xexp_pu[lx];
  xexp_pci = xexp_pc[lx];
  xexp_pdi = xexp_pd[lx];
  ximp_pui = ximp_pu[lx];
  ximp_pci = ximp_pc[lx];
  ximp_pdi = ximp_pd[lx];

  for (j = ly; j <= uy; j++) {
    xexp_pdi[j] = 0.0;
    xexp_pci[j] = -0.5 * muxi[j] / dux;
    xexp_pui[j] = -xexp_pci[j];

    ximp_pui[j] = -xexp_pui[j];
    ximp_pdi[j] = -xexp_pdi[j];
    ximp_pci[j] = -xexp_pci[j];

    temp2 = 0.25 * ri[j];
    xexp_pci[j] += 1.0 - temp2;
    ximp_pci[j] += 1.0 + temp2;
  }

  /*	Up probas on x */
  ddx = x[ux] - x[ux - 1];

  muxi = mux[ux];
  varxi = varx[ux];
  ri = r[ux];
  xexp_pui = xexp_pu[ux];
  xexp_pci = xexp_pc[ux];
  xexp_pdi = xexp_pd[ux];
  ximp_pui = ximp_pu[ux];
  ximp_pci = ximp_pc[ux];
  ximp_pdi = ximp_pd[ux];

  for (j = ly; j <= uy; j++) {
    xexp_pui[j] = 0.0;
    xexp_pci[j] = 0.5 * muxi[j] / ddx;
    xexp_pdi[j] = -xexp_pci[j];

    ximp_pui[j] = -xexp_pui[j];
    ximp_pdi[j] = -xexp_pdi[j];
    ximp_pci[j] = -xexp_pci[j];

    temp2 = 0.25 * ri[j];
    xexp_pci[j] += 1.0 - temp2;
    ximp_pci[j] += 1.0 + temp2;
  }

  /*	For each slice in y  , implicit in x */

  for (j = ly; j <= uy; j++) {
    for (k = mstart; k <= mend; k++) {
      f_t[lx][j][k] = f_t_plus_1[lx][j][k] / (bet = ximp_pc[lx][j]);
    }

    for (i = lx + 1; i <= ux; i++) {
      gamx[i] = ximp_pu[i - 1][j] / bet;
      bet = ximp_pc[i][j] - ximp_pd[i][j] * gamx[i];
      for (k = mstart; k <= mend; k++) {
        f_t[i][j][k] =
            (f_t_plus_1[i][j][k] - ximp_pd[i][j] * f_t[i - 1][j][k]) / bet;
      }
    }

    for (i = ux - 1; i >= lx; i--) {
      for (k = mstart; k <= mend; k++) {
        f_t[i][j][k] -= gamx[i + 1] * f_t[i + 1][j][k];
      }
    }
  }

  /*	For each slice in x  , explicit in y */

  for (i = lx; i <= ux; i++) {
    fi = f_t[i];
    tmpfi = tmpf[i];
    yexp_pui = yexp_pu[i];
    yexp_pci = yexp_pc[i];
    yexp_pdi = yexp_pd[i];

    for (k = mstart; k <= mend; k++) {
      for (j = ly + 1; j <= uy - 1; j++) {
        tmpfi[j][k - mstart] = yexp_pui[j] * fi[j + 1][k] +
                               yexp_pci[j] * fi[j][k] +
                               yexp_pdi[j] * fi[j - 1][k];
      }

      tmpfi[ly][k - mstart] =
          yexp_pui[ly] * fi[ly + 1][k] + yexp_pci[ly] * fi[ly][k];

      tmpfi[uy][k - mstart] =
          yexp_pci[uy] * fi[uy][k] + yexp_pdi[uy] * fi[uy - 1][k];
    }
  }

  /*	For each slice in x  , implicit in y */

  for (i = lx; i <= ux; i++) {
    tmpfi = tmpf[i];
    yimp_pui = yimp_pu[i];
    yimp_pci = yimp_pc[i];
    yimp_pdi = yimp_pd[i];

    for (k = mstart; k <= mend; k++) {
      tmpfi[ly][k - mstart] /= (bet = yimp_pci[ly]);
    }

    for (j = ly + 1; j <= uy; j++) {
      gamy[j] = yimp_pui[j - 1] / bet;
      bet = yimp_pci[j] - yimp_pdi[j] * gamy[j];
      for (k = mstart; k <= mend; k++) {
        tmpfi[j][k - mstart] =
            (tmpfi[j][k - mstart] - yimp_pdi[j] * tmpfi[j - 1][k - mstart]) /
            bet;
      }
    }

    for (j = uy - 1; j >= ly; j--) {
      for (k = mstart; k <= mend; k++) {
        tmpfi[j][k - mstart] -= gamy[j + 1] * tmpfi[j + 1][k - mstart];
      }
    }
  }

  /*	For each slice in y  , explicit in x */

  for (j = ly; j <= uy; j++) {
    for (k = mstart; k <= mend; k++) {
      for (i = lx + 1; i <= ux - 1; i++) {
        f_t[i][j][k] = xexp_pu[i][j] * tmpf[i + 1][j][k - mstart] +
                       xexp_pc[i][j] * tmpf[i][j][k - mstart] +
                       xexp_pd[i][j] * tmpf[i - 1][j][k - mstart];
      }

      f_t[lx][j][k] = xexp_pu[lx][j] * tmpf[lx + 1][j][k - mstart] +
                      xexp_pc[lx][j] * tmpf[lx][j][k - mstart];

      f_t[ux][j][k] = xexp_pc[ux][j] * tmpf[ux][j][k - mstart] +
                      xexp_pd[ux][j] * tmpf[ux - 1][j][k - mstart];
    }
  }
}

/*	Convolve one step backward using ADI scheme */
void num_f_pde_one_step_backward_2f_save_adi(
    CNPDE_TEMP_2D_ADI *tmp,
    /*	Size of the X array	*/
    int nx,
    /*	X values	*/
    double *x,
    /*	Size of the Y array	*/
    int ny,
    /*	X values	*/
    double *y,
    /*	Number of functions	*/
    int mstart, int mend,
    /*	f values @ t+1	*/
    /*	i: x state  , j: y state  , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy, int recalc) {
  static int i, j, k;
  static double dux, ddx, duy, ddy, bet;
  static double *muxi, *varxi, *ri, **fi, **tmpfi, *xexp_pui, *xexp_pci,
      *xexp_pdi, *ximp_pui, *ximp_pci, *ximp_pdi, *yexp_pui, *yexp_pci,
      *yexp_pdi, *yimp_pui, *yimp_pci, *yimp_pdi;
  static double temp1, temp2;

  double **xexp_pu = tmp->xexp_pu, **xexp_pc = tmp->xexp_pc,
         **xexp_pd = tmp->xexp_pd,

         **ximp_pu = tmp->ximp_pu, **ximp_pc = tmp->ximp_pc,
         **ximp_pd = tmp->ximp_pd,

         **yexp_pu = tmp->yexp_pu, **yexp_pc = tmp->yexp_pc,
         **yexp_pd = tmp->yexp_pd,

         **yimp_pu = tmp->yimp_pu, **yimp_pc = tmp->yimp_pc,
         **yimp_pd = tmp->yimp_pd,

         ***tmpf = tmp->f,

         *gamx = tmp->gamx, *gamy = tmp->gamy;

  if (recalc) {

    /*	Calculate all probabilities */

    /*	Interior probas on y */
    for (j = ly + 1; j <= uy - 1; j++) {
      duy = y[j + 1] - y[j];
      ddy = y[j] - y[j - 1];
      temp1 = 0.5 / (duy * (duy + ddy));

      for (i = lx; i <= ux; i++) {
        yexp_pu[i][j] = temp1 * (vary[i][j] + muy[i][j] * ddy);
        yexp_pd[i][j] = (yexp_pu[i][j] * duy - 0.5 * muy[i][j]) / ddy;
        yexp_pc[i][j] = -yexp_pu[i][j] - yexp_pd[i][j];

        yimp_pu[i][j] = -yexp_pu[i][j];
        yimp_pd[i][j] = -yexp_pd[i][j];
        yimp_pc[i][j] = -yexp_pc[i][j];

        temp2 = 0.25 * r[i][j];
        yexp_pc[i][j] += 1.0 - temp2;
        yimp_pc[i][j] += 1.0 + temp2;
      }
    }

    /*	Down probas on y */
    duy = y[ly + 1] - y[ly];

    for (i = lx; i <= ux; i++) {
      yexp_pd[i][ly] = 0.0;
      yexp_pc[i][ly] = -0.5 * muy[i][ly] / duy;
      yexp_pu[i][ly] = -yexp_pc[i][ly];

      yimp_pu[i][ly] = -yexp_pu[i][ly];
      yimp_pd[i][ly] = -yexp_pd[i][ly];
      yimp_pc[i][ly] = -yexp_pc[i][ly];

      temp2 = 0.25 * r[i][ly];
      yexp_pc[i][ly] += 1.0 - temp2;
      yimp_pc[i][ly] += 1.0 + temp2;
    }

    /*	Up probas on y */
    ddy = y[uy] - y[uy - 1];

    for (i = lx; i <= ux; i++) {
      yexp_pu[i][uy] = 0.0;
      yexp_pc[i][uy] = 0.5 * muy[i][uy] / ddy;
      yexp_pd[i][uy] = -yexp_pc[i][uy];

      yimp_pu[i][uy] = -yexp_pu[i][uy];
      yimp_pd[i][uy] = -yexp_pd[i][uy];
      yimp_pc[i][uy] = -yexp_pc[i][uy];

      temp2 = 0.25 * r[i][uy];
      yexp_pc[i][uy] += 1.0 - temp2;
      yimp_pc[i][uy] += 1.0 + temp2;
    }

    /*	Interior probas on x */
    for (i = lx + 1; i <= ux - 1; i++) {
      dux = x[i + 1] - x[i];
      ddx = x[i] - x[i - 1];
      temp1 = 0.5 / (dux * (dux + ddx));

      muxi = mux[i];
      varxi = varx[i];
      ri = r[i];
      xexp_pui = xexp_pu[i];
      xexp_pci = xexp_pc[i];
      xexp_pdi = xexp_pd[i];
      ximp_pui = ximp_pu[i];
      ximp_pci = ximp_pc[i];
      ximp_pdi = ximp_pd[i];

      for (j = ly; j <= uy; j++) {
        xexp_pui[j] = temp1 * (varxi[j] + muxi[j] * ddx);
        xexp_pdi[j] = (xexp_pui[j] * dux - 0.5 * muxi[j]) / ddx;
        xexp_pci[j] = -xexp_pui[j] - xexp_pdi[j];

        ximp_pui[j] = -xexp_pui[j];
        ximp_pdi[j] = -xexp_pdi[j];
        ximp_pci[j] = -xexp_pci[j];

        temp2 = 0.25 * ri[j];
        xexp_pci[j] += 1.0 - temp2;
        ximp_pci[j] += 1.0 + temp2;
      }
    }

    /*	Down probas on x */
    dux = x[lx + 1] - x[lx];

    muxi = mux[lx];
    varxi = varx[lx];
    ri = r[lx];
    xexp_pui = xexp_pu[lx];
    xexp_pci = xexp_pc[lx];
    xexp_pdi = xexp_pd[lx];
    ximp_pui = ximp_pu[lx];
    ximp_pci = ximp_pc[lx];
    ximp_pdi = ximp_pd[lx];

    for (j = ly; j <= uy; j++) {
      xexp_pdi[j] = 0.0;
      xexp_pci[j] = -0.5 * muxi[j] / dux;
      xexp_pui[j] = -xexp_pci[j];

      ximp_pui[j] = -xexp_pui[j];
      ximp_pdi[j] = -xexp_pdi[j];
      ximp_pci[j] = -xexp_pci[j];

      temp2 = 0.25 * ri[j];
      xexp_pci[j] += 1.0 - temp2;
      ximp_pci[j] += 1.0 + temp2;
    }

    /*	Up probas on x */
    ddx = x[ux] - x[ux - 1];

    muxi = mux[ux];
    varxi = varx[ux];
    ri = r[ux];
    xexp_pui = xexp_pu[ux];
    xexp_pci = xexp_pc[ux];
    xexp_pdi = xexp_pd[ux];
    ximp_pui = ximp_pu[ux];
    ximp_pci = ximp_pc[ux];
    ximp_pdi = ximp_pd[ux];

    for (j = ly; j <= uy; j++) {
      xexp_pui[j] = 0.0;
      xexp_pci[j] = 0.5 * muxi[j] / ddx;
      xexp_pdi[j] = -xexp_pci[j];

      ximp_pui[j] = -xexp_pui[j];
      ximp_pdi[j] = -xexp_pdi[j];
      ximp_pci[j] = -xexp_pci[j];

      temp2 = 0.25 * ri[j];
      xexp_pci[j] += 1.0 - temp2;
      ximp_pci[j] += 1.0 + temp2;
    }
  }

  /*	For each slice in y  , implicit in x */

  for (j = ly; j <= uy; j++) {
    for (k = mstart; k <= mend; k++) {
      f_t[lx][j][k] = f_t_plus_1[lx][j][k] / (bet = ximp_pc[lx][j]);
    }

    for (i = lx + 1; i <= ux; i++) {
      gamx[i] = ximp_pu[i - 1][j] / bet;
      bet = ximp_pc[i][j] - ximp_pd[i][j] * gamx[i];
      for (k = mstart; k <= mend; k++) {
        f_t[i][j][k] =
            (f_t_plus_1[i][j][k] - ximp_pd[i][j] * f_t[i - 1][j][k]) / bet;
      }
    }

    for (i = ux - 1; i >= lx; i--) {
      for (k = mstart; k <= mend; k++) {
        f_t[i][j][k] -= gamx[i + 1] * f_t[i + 1][j][k];
      }
    }
  }

  /*	For each slice in x  , explicit in y */

  for (i = lx; i <= ux; i++) {
    fi = f_t[i];
    tmpfi = tmpf[i];
    yexp_pui = yexp_pu[i];
    yexp_pci = yexp_pc[i];
    yexp_pdi = yexp_pd[i];

    for (k = mstart; k <= mend; k++) {
      for (j = ly + 1; j <= uy - 1; j++) {
        tmpfi[j][k - mstart] = yexp_pui[j] * fi[j + 1][k] +
                               yexp_pci[j] * fi[j][k] +
                               yexp_pdi[j] * fi[j - 1][k];
      }

      tmpfi[ly][k - mstart] =
          yexp_pui[ly] * fi[ly + 1][k] + yexp_pci[ly] * fi[ly][k];

      tmpfi[uy][k - mstart] =
          yexp_pci[uy] * fi[uy][k] + yexp_pdi[uy] * fi[uy - 1][k];
    }
  }

  /*	For each slice in x  , implicit in y */

  for (i = lx; i <= ux; i++) {
    tmpfi = tmpf[i];
    yimp_pui = yimp_pu[i];
    yimp_pci = yimp_pc[i];
    yimp_pdi = yimp_pd[i];

    for (k = mstart; k <= mend; k++) {
      tmpfi[ly][k - mstart] /= (bet = yimp_pci[ly]);
    }

    for (j = ly + 1; j <= uy; j++) {
      gamy[j] = yimp_pui[j - 1] / bet;
      bet = yimp_pci[j] - yimp_pdi[j] * gamy[j];
      for (k = mstart; k <= mend; k++) {
        tmpfi[j][k - mstart] =
            (tmpfi[j][k - mstart] - yimp_pdi[j] * tmpfi[j - 1][k - mstart]) /
            bet;
      }
    }

    for (j = uy - 1; j >= ly; j--) {
      for (k = mstart; k <= mend; k++) {
        tmpfi[j][k - mstart] -= gamy[j + 1] * tmpfi[j + 1][k - mstart];
      }
    }
  }

  /*	For each slice in y  , explicit in x */

  for (j = ly; j <= uy; j++) {
    for (k = mstart; k <= mend; k++) {
      for (i = lx + 1; i <= ux - 1; i++) {
        f_t[i][j][k] = xexp_pu[i][j] * tmpf[i + 1][j][k - mstart] +
                       xexp_pc[i][j] * tmpf[i][j][k - mstart] +
                       xexp_pd[i][j] * tmpf[i - 1][j][k - mstart];
      }

      f_t[lx][j][k] = xexp_pu[lx][j] * tmpf[lx + 1][j][k - mstart] +
                      xexp_pc[lx][j] * tmpf[lx][j][k - mstart];

      f_t[ux][j][k] = xexp_pc[ux][j] * tmpf[ux][j][k - mstart] +
                      xexp_pd[ux][j] * tmpf[ux - 1][j][k - mstart];
    }
  }
}

/*	Convolve one step backward using ADI scheme */
void num_f_pde_one_step_backward_2f_adi_LGMSV(
    CNPDE_TEMP_2D_ADI *tmp,
    /*	Size of the X array	*/
    int nx,
    /*	X values	*/
    double *x,
    /*	Size of the Y array	*/
    int ny,
    /*	X values	*/
    double *y,
    /*	Number of functions	*/
    int mstart, int mend,
    /*	f values @ t+1	*/
    /*	i: x state  , j: y state  , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int *ly, int *uy,

    int *invlx, int *invux, int invly, int invuy) {
  static int i, j, k;
  static double dux, ddx, duy, ddy, bet;
  static double *muxi, *varxi, *ri, **fi, **tmpfi, *xexp_pui, *xexp_pci,
      *xexp_pdi, *ximp_pui, *ximp_pci, *ximp_pdi, *yexp_pui, *yexp_pci,
      *yexp_pdi, *yimp_pui, *yimp_pci, *yimp_pdi;
  static double temp1, temp2;

  double **xexp_pu = tmp->xexp_pu, **xexp_pc = tmp->xexp_pc,
         **xexp_pd = tmp->xexp_pd,

         **ximp_pu = tmp->ximp_pu, **ximp_pc = tmp->ximp_pc,
         **ximp_pd = tmp->ximp_pd,

         **yexp_pu = tmp->yexp_pu, **yexp_pc = tmp->yexp_pc,
         **yexp_pd = tmp->yexp_pd,

         **yimp_pu = tmp->yimp_pu, **yimp_pc = tmp->yimp_pc,
         **yimp_pd = tmp->yimp_pd,

         ***tmpf = tmp->f,

         *gamx = tmp->gamx, *gamy = tmp->gamy;

  /*	Calculate all probabilities */

  /*	Interior probas on y */
  for (j = invly + 1; j <= invuy - 1; j++) {
    duy = y[j + 1] - y[j];
    ddy = y[j] - y[j - 1];
    temp1 = 0.5 / (duy * (duy + ddy));

    for (i = invlx[j]; i <= invux[j]; i++) {
      yexp_pu[i][j] = temp1 * (vary[i][j] + muy[i][j] * ddy);
      yexp_pd[i][j] = (yexp_pu[i][j] * duy - 0.5 * muy[i][j]) / ddy;
      yexp_pc[i][j] = -yexp_pu[i][j] - yexp_pd[i][j];

      yimp_pu[i][j] = -yexp_pu[i][j];
      yimp_pd[i][j] = -yexp_pd[i][j];
      yimp_pc[i][j] = -yexp_pc[i][j];

      temp2 = 0.25 * r[i][j];
      yexp_pc[i][j] += 1.0 - temp2;
      yimp_pc[i][j] += 1.0 + temp2;
    }
  }

  /*	Down probas on y */
  duy = y[invly + 1] - y[invly];

  for (i = invlx[invly]; i <= invux[invly]; i++) {
    k = ly[i];
    yexp_pd[i][k] = 0.0;
    yexp_pc[i][k] = -0.5 * muy[i][k] / duy;
    yexp_pu[i][k] = -yexp_pc[i][k];

    yimp_pu[i][k] = -yexp_pu[i][k];
    yimp_pd[i][k] = -yexp_pd[i][k];
    yimp_pc[i][k] = -yexp_pc[i][k];

    temp2 = 0.25 * r[i][k];
    yexp_pc[i][k] += 1.0 - temp2;
    yimp_pc[i][k] += 1.0 + temp2;
  }

  /*	Up probas on y */
  ddy = y[invuy] - y[invuy - 1];

  for (i = invlx[invuy]; i <= invux[invuy]; i++) {
    k = uy[i];
    yexp_pu[i][k] = 0.0;
    yexp_pc[i][k] = 0.5 * muy[i][k] / ddy;
    yexp_pd[i][k] = -yexp_pc[i][k];

    yimp_pu[i][k] = -yexp_pu[i][k];
    yimp_pd[i][k] = -yexp_pd[i][k];
    yimp_pc[i][k] = -yexp_pc[i][k];

    temp2 = 0.25 * r[i][k];
    yexp_pc[i][k] += 1.0 - temp2;
    yimp_pc[i][k] += 1.0 + temp2;
  }

  /*	Interior probas on x */
  for (i = lx + 1; i <= ux - 1; i++) {
    dux = x[i + 1] - x[i];
    ddx = x[i] - x[i - 1];
    temp1 = 0.5 / (dux * (dux + ddx));

    muxi = mux[i];
    varxi = varx[i];
    ri = r[i];
    xexp_pui = xexp_pu[i];
    xexp_pci = xexp_pc[i];
    xexp_pdi = xexp_pd[i];
    ximp_pui = ximp_pu[i];
    ximp_pci = ximp_pc[i];
    ximp_pdi = ximp_pd[i];

    for (j = ly[i]; j <= uy[i]; j++) {
      xexp_pui[j] = temp1 * (varxi[j] + muxi[j] * ddx);
      xexp_pdi[j] = (xexp_pui[j] * dux - 0.5 * muxi[j]) / ddx;
      xexp_pci[j] = -xexp_pui[j] - xexp_pdi[j];

      ximp_pui[j] = -xexp_pui[j];
      ximp_pdi[j] = -xexp_pdi[j];
      ximp_pci[j] = -xexp_pci[j];

      temp2 = 0.25 * ri[j];
      xexp_pci[j] += 1.0 - temp2;
      ximp_pci[j] += 1.0 + temp2;
    }
  }

  /*	Down probas on x */
  dux = x[lx + 1] - x[lx];

  muxi = mux[lx];
  varxi = varx[lx];
  ri = r[lx];
  xexp_pui = xexp_pu[lx];
  xexp_pci = xexp_pc[lx];
  xexp_pdi = xexp_pd[lx];
  ximp_pui = ximp_pu[lx];
  ximp_pci = ximp_pc[lx];
  ximp_pdi = ximp_pd[lx];

  for (j = ly[lx]; j <= uy[lx]; j++) {
    xexp_pdi[j] = 0.0;
    xexp_pci[j] = -0.5 * muxi[j] / dux;
    xexp_pui[j] = -xexp_pci[j];

    ximp_pui[j] = -xexp_pui[j];
    ximp_pdi[j] = -xexp_pdi[j];
    ximp_pci[j] = -xexp_pci[j];

    temp2 = 0.25 * ri[j];
    xexp_pci[j] += 1.0 - temp2;
    ximp_pci[j] += 1.0 + temp2;
  }

  /*	Up probas on x */
  ddx = x[ux] - x[ux - 1];

  muxi = mux[ux];
  varxi = varx[ux];
  ri = r[ux];
  xexp_pui = xexp_pu[ux];
  xexp_pci = xexp_pc[ux];
  xexp_pdi = xexp_pd[ux];
  ximp_pui = ximp_pu[ux];
  ximp_pci = ximp_pc[ux];
  ximp_pdi = ximp_pd[ux];

  for (j = ly[ux]; j <= uy[ux]; j++) {
    xexp_pui[j] = 0.0;
    xexp_pci[j] = 0.5 * muxi[j] / ddx;
    xexp_pdi[j] = -xexp_pci[j];

    ximp_pui[j] = -xexp_pui[j];
    ximp_pdi[j] = -xexp_pdi[j];
    ximp_pci[j] = -xexp_pci[j];

    temp2 = 0.25 * ri[j];
    xexp_pci[j] += 1.0 - temp2;
    ximp_pci[j] += 1.0 + temp2;
  }

  /*	For each slice in y  , implicit in x */

  for (j = invly; j <= invuy; j++) {
    for (k = mstart; k <= mend; k++) {
      f_t[invlx[j]][j][k] =
          f_t_plus_1[invlx[j]][j][k] / (bet = ximp_pc[invlx[j]][j]);
    }

    for (i = invlx[j] + 1; i <= invux[j]; i++) {
      gamx[i] = ximp_pu[i - 1][j] / bet;
      bet = ximp_pc[i][j] - ximp_pd[i][j] * gamx[i];
      for (k = mstart; k <= mend; k++) {
        f_t[i][j][k] =
            (f_t_plus_1[i][j][k] - ximp_pd[i][j] * f_t[i - 1][j][k]) / bet;
      }
    }

    for (i = invux[j] - 1; i >= invlx[j]; i--) {
      for (k = mstart; k <= mend; k++) {
        f_t[i][j][k] -= gamx[i + 1] * f_t[i + 1][j][k];
      }
    }
  }

  /*	For each slice in x  , explicit in y */

  for (i = lx; i <= ux; i++) {
    fi = f_t[i];
    tmpfi = tmpf[i];
    yexp_pui = yexp_pu[i];
    yexp_pci = yexp_pc[i];
    yexp_pdi = yexp_pd[i];

    for (k = mstart; k <= mend; k++) {
      for (j = ly[i] + 1; j <= uy[i] - 1; j++) {
        tmpfi[j][k - mstart] = yexp_pui[j] * fi[j + 1][k] +
                               yexp_pci[j] * fi[j][k] +
                               yexp_pdi[j] * fi[j - 1][k];
      }

      tmpfi[ly[i]][k - mstart] =
          yexp_pui[ly[i]] * fi[ly[i] + 1][k] + yexp_pci[ly[i]] * fi[ly[i]][k];

      tmpfi[uy[i]][k - mstart] =
          yexp_pci[uy[i]] * fi[uy[i]][k] + yexp_pdi[uy[i]] * fi[uy[i] - 1][k];
    }
  }

  /*	For each slice in x  , implicit in y */

  for (i = lx; i <= ux; i++) {
    tmpfi = tmpf[i];
    yimp_pui = yimp_pu[i];
    yimp_pci = yimp_pc[i];
    yimp_pdi = yimp_pd[i];

    for (k = mstart; k <= mend; k++) {
      tmpfi[ly[i]][k - mstart] /= (bet = yimp_pci[ly[i]]);
    }

    for (j = ly[i] + 1; j <= uy[i]; j++) {
      gamy[j] = yimp_pui[j - 1] / bet;
      bet = yimp_pci[j] - yimp_pdi[j] * gamy[j];
      for (k = mstart; k <= mend; k++) {
        tmpfi[j][k - mstart] =
            (tmpfi[j][k - mstart] - yimp_pdi[j] * tmpfi[j - 1][k - mstart]) /
            bet;
      }
    }

    for (j = uy[i] - 1; j >= ly[i]; j--) {
      for (k = mstart; k <= mend; k++) {
        tmpfi[j][k - mstart] -= gamy[j + 1] * tmpfi[j + 1][k - mstart];
      }
    }
  }

  /*	For each slice in y  , explicit in x */

  for (j = invly; j <= invuy; j++) {
    for (k = mstart; k <= mend; k++) {
      for (i = invlx[j] + 1; i <= invux[j] - 1; i++) {
        f_t[i][j][k] = xexp_pu[i][j] * tmpf[i + 1][j][k - mstart] +
                       xexp_pc[i][j] * tmpf[i][j][k - mstart] +
                       xexp_pd[i][j] * tmpf[i - 1][j][k - mstart];
      }

      f_t[invlx[j]][j][k] =
          xexp_pu[invlx[j]][j] * tmpf[invlx[j] + 1][j][k - mstart] +
          xexp_pc[invlx[j]][j] * tmpf[invlx[j]][j][k - mstart];

      f_t[invux[j]][j][k] =
          xexp_pc[invux[j]][j] * tmpf[invux[j]][j][k - mstart] +
          xexp_pd[invux[j]][j] * tmpf[invux[j] - 1][j][k - mstart];
    }
  }
}

/*	Calculate all proba for backward ADI scheme */
void calculate_proba_backward_2f_adi(
    CNPDE_TEMP_2D_ADI *tmp,
    /*	Size of the X array	*/
    int nx,
    /*	X values	*/
    double *x,
    /*	Size of the Y array	*/
    int ny,
    /*	X values	*/
    double *y,
    /*	Number of functions	*/
    int m,
    /*	f values @ t+1	*/
    /*	i: x state  , j: y state  , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy) {
  static int i, j, k;
  static double dux, ddx, duy, ddy, bet;
  static double *muxi, *varxi, *ri, **fi, **tmpfi, *xexp_pui, *xexp_pci,
      *xexp_pdi, *ximp_pui, *ximp_pci, *ximp_pdi, *yexp_pui, *yexp_pci,
      *yexp_pdi, *yimp_pui, *yimp_pci, *yimp_pdi;
  static double temp1, temp2;

  double **xexp_pu = tmp->xexp_pu, **xexp_pc = tmp->xexp_pc,
         **xexp_pd = tmp->xexp_pd,

         **ximp_pu = tmp->ximp_pu, **ximp_pc = tmp->ximp_pc,
         **ximp_pd = tmp->ximp_pd,

         **yexp_pu = tmp->yexp_pu, **yexp_pc = tmp->yexp_pc,
         **yexp_pd = tmp->yexp_pd,

         **yimp_pu = tmp->yimp_pu, **yimp_pc = tmp->yimp_pc,
         **yimp_pd = tmp->yimp_pd,

         ***tmpf = tmp->f,

         *gamx = tmp->gamx, *gamy = tmp->gamy;

  /*	Calculate all probabilities */

  /*	Interior probas on y */
  for (j = ly + 1; j <= uy - 1; j++) {
    duy = y[j + 1] - y[j];
    ddy = y[j] - y[j - 1];
    temp1 = 0.5 / (duy * (duy + ddy));

    for (i = lx; i <= ux; i++) {
      yexp_pu[i][j] = temp1 * (vary[i][j] + muy[i][j] * ddy);
      yexp_pd[i][j] = (yexp_pu[i][j] * duy - 0.5 * muy[i][j]) / ddy;
      yexp_pc[i][j] = -yexp_pu[i][j] - yexp_pd[i][j];

      yimp_pu[i][j] = -yexp_pu[i][j];
      yimp_pd[i][j] = -yexp_pd[i][j];
      yimp_pc[i][j] = -yexp_pc[i][j];

      temp2 = 0.25 * r[i][j];
      yexp_pc[i][j] += 1.0 - temp2;
      yimp_pc[i][j] += 1.0 + temp2;
    }
  }

  /*	Down probas on y */
  duy = y[ly + 1] - y[ly];

  for (i = lx; i <= ux; i++) {
    yexp_pd[i][ly] = 0.0;
    yexp_pc[i][ly] = -0.5 * muy[i][ly] / duy;
    yexp_pu[i][ly] = -yexp_pc[i][ly];

    yimp_pu[i][ly] = -yexp_pu[i][ly];
    yimp_pd[i][ly] = -yexp_pd[i][ly];
    yimp_pc[i][ly] = -yexp_pc[i][ly];

    temp2 = 0.25 * r[i][ly];
    yexp_pc[i][ly] += 1.0 - temp2;
    yimp_pc[i][ly] += 1.0 + temp2;
  }

  /*	Up probas on y */
  ddy = y[uy] - y[uy - 1];

  for (i = lx; i <= ux; i++) {
    yexp_pu[i][uy] = 0.0;
    yexp_pc[i][uy] = 0.5 * muy[i][uy] / ddy;
    yexp_pd[i][uy] = -yexp_pc[i][uy];

    yimp_pu[i][uy] = -yexp_pu[i][uy];
    yimp_pd[i][uy] = -yexp_pd[i][uy];
    yimp_pc[i][uy] = -yexp_pc[i][uy];

    temp2 = 0.25 * r[i][uy];
    yexp_pc[i][uy] += 1.0 - temp2;
    yimp_pc[i][uy] += 1.0 + temp2;
  }

  /*	Interior probas on x */
  for (i = lx + 1; i <= ux - 1; i++) {
    dux = x[i + 1] - x[i];
    ddx = x[i] - x[i - 1];
    temp1 = 0.5 / (dux * (dux + ddx));

    muxi = mux[i];
    varxi = varx[i];
    ri = r[i];
    xexp_pui = xexp_pu[i];
    xexp_pci = xexp_pc[i];
    xexp_pdi = xexp_pd[i];
    ximp_pui = ximp_pu[i];
    ximp_pci = ximp_pc[i];
    ximp_pdi = ximp_pd[i];

    for (j = ly; j <= uy; j++) {
      xexp_pui[j] = temp1 * (varxi[j] + muxi[j] * ddx);
      xexp_pdi[j] = (xexp_pui[j] * dux - 0.5 * muxi[j]) / ddx;
      xexp_pci[j] = -xexp_pui[j] - xexp_pdi[j];

      ximp_pui[j] = -xexp_pui[j];
      ximp_pdi[j] = -xexp_pdi[j];
      ximp_pci[j] = -xexp_pci[j];

      temp2 = 0.25 * ri[j];
      xexp_pci[j] += 1.0 - temp2;
      ximp_pci[j] += 1.0 + temp2;
    }
  }

  /*	Down probas on x */
  dux = x[lx + 1] - x[lx];

  muxi = mux[lx];
  varxi = varx[lx];
  ri = r[lx];
  xexp_pui = xexp_pu[lx];
  xexp_pci = xexp_pc[lx];
  xexp_pdi = xexp_pd[lx];
  ximp_pui = ximp_pu[lx];
  ximp_pci = ximp_pc[lx];
  ximp_pdi = ximp_pd[lx];

  for (j = ly; j <= uy; j++) {
    xexp_pdi[j] = 0.0;
    xexp_pci[j] = -0.5 * muxi[j] / dux;
    xexp_pui[j] = -xexp_pci[j];

    ximp_pui[j] = -xexp_pui[j];
    ximp_pdi[j] = -xexp_pdi[j];
    ximp_pci[j] = -xexp_pci[j];

    temp2 = 0.25 * ri[j];
    xexp_pci[j] += 1.0 - temp2;
    ximp_pci[j] += 1.0 + temp2;
  }

  /*	Up probas on x */
  ddx = x[ux] - x[ux - 1];

  muxi = mux[ux];
  varxi = varx[ux];
  ri = r[ux];
  xexp_pui = xexp_pu[ux];
  xexp_pci = xexp_pc[ux];
  xexp_pdi = xexp_pd[ux];
  ximp_pui = ximp_pu[ux];
  ximp_pci = ximp_pc[ux];
  ximp_pdi = ximp_pd[ux];

  for (j = ly; j <= uy; j++) {
    xexp_pui[j] = 0.0;
    xexp_pci[j] = 0.5 * muxi[j] / ddx;
    xexp_pdi[j] = -xexp_pci[j];

    ximp_pui[j] = -xexp_pui[j];
    ximp_pdi[j] = -xexp_pdi[j];
    ximp_pci[j] = -xexp_pci[j];

    temp2 = 0.25 * ri[j];
    xexp_pci[j] += 1.0 - temp2;
    ximp_pci[j] += 1.0 + temp2;
  }
}

/*	Convolve one step backward using ADI scheme with precalculated probas */
void num_f_pde_one_step_backward_2f_adi_precalc_proba(
    CNPDE_TEMP_2D_ADI *tmp,
    /*	Size of the X array	*/
    int nx,
    /*	X values	*/
    double *x,
    /*	Size of the Y array	*/
    int ny,
    /*	X values	*/
    double *y,
    /*	Number of functions	*/
    int m,
    /*	f values @ t+1	*/
    /*	i: x state  , j: y state  , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy) {
  static int i, j, k;
  static double dux, ddx, duy, ddy, bet;
  static double *muxi, *varxi, ri, **fi, **tmpfi, *xexp_pui, *xexp_pci,
      *xexp_pdi, *ximp_pui, *ximp_pci, *ximp_pdi, *yexp_pui, *yexp_pci,
      *yexp_pdi, *yimp_pui, *yimp_pci, *yimp_pdi;
  static double temp1, temp2;

  double **xexp_pu = tmp->xexp_pu, **xexp_pc = tmp->xexp_pc,
         **xexp_pd = tmp->xexp_pd,

         **ximp_pu = tmp->ximp_pu, **ximp_pc = tmp->ximp_pc,
         **ximp_pd = tmp->ximp_pd,

         **yexp_pu = tmp->yexp_pu, **yexp_pc = tmp->yexp_pc,
         **yexp_pd = tmp->yexp_pd,

         **yimp_pu = tmp->yimp_pu, **yimp_pc = tmp->yimp_pc,
         **yimp_pd = tmp->yimp_pd,

         ***tmpf = tmp->f,

         *gamx = tmp->gamx, *gamy = tmp->gamy;

  /*	For each slice in y  , implicit in x */

  for (j = ly; j <= uy; j++) {
    for (k = 0; k < m; k++) {
      f_t[lx][j][k] = f_t_plus_1[lx][j][k] / (bet = ximp_pc[lx][j]);
    }

    for (i = lx + 1; i <= ux; i++) {
      gamx[i] = ximp_pu[i - 1][j] / bet;
      bet = ximp_pc[i][j] - ximp_pd[i][j] * gamx[i];
      for (k = 0; k < m; k++) {
        f_t[i][j][k] =
            (f_t_plus_1[i][j][k] - ximp_pd[i][j] * f_t[i - 1][j][k]) / bet;
      }
    }

    for (i = ux - 1; i >= lx; i--) {
      for (k = 0; k < m; k++) {
        f_t[i][j][k] -= gamx[i + 1] * f_t[i + 1][j][k];
      }
    }
  }

  /*	For each slice in x  , explicit in y */

  for (i = lx; i <= ux; i++) {
    fi = f_t[i];
    tmpfi = tmpf[i];
    yexp_pui = yexp_pu[i];
    yexp_pci = yexp_pc[i];
    yexp_pdi = yexp_pd[i];

    for (k = 0; k < m; k++) {
      for (j = ly + 1; j <= uy - 1; j++) {
        tmpfi[j][k] = yexp_pui[j] * fi[j + 1][k] + yexp_pci[j] * fi[j][k] +
                      yexp_pdi[j] * fi[j - 1][k];
      }

      tmpfi[ly][k] = yexp_pui[ly] * fi[ly + 1][k] + yexp_pci[ly] * fi[ly][k];

      tmpfi[uy][k] = yexp_pci[uy] * fi[uy][k] + yexp_pdi[uy] * fi[uy - 1][k];
    }
  }

  /*	For each slice in x  , implicit in y */

  for (i = lx; i <= ux; i++) {
    tmpfi = tmpf[i];
    yimp_pui = yimp_pu[i];
    yimp_pci = yimp_pc[i];
    yimp_pdi = yimp_pd[i];

    for (k = 0; k < m; k++) {
      tmpfi[ly][k] /= (bet = yimp_pci[ly]);
    }

    for (j = ly + 1; j <= uy; j++) {
      gamy[j] = yimp_pui[j - 1] / bet;
      bet = yimp_pci[j] - yimp_pdi[j] * gamy[j];
      for (k = 0; k < m; k++) {
        tmpfi[j][k] = (tmpfi[j][k] - yimp_pdi[j] * tmpfi[j - 1][k]) / bet;
      }
    }

    for (j = uy - 1; j >= ly; j--) {
      for (k = 0; k < m; k++) {
        tmpfi[j][k] -= gamy[j + 1] * tmpfi[j + 1][k];
      }
    }
  }

  /*	For each slice in y  , explicit in x */

  for (j = ly; j <= uy; j++) {
    for (k = 0; k < m; k++) {
      for (i = lx + 1; i <= ux - 1; i++) {
        f_t[i][j][k] = xexp_pu[i][j] * tmpf[i + 1][j][k] +
                       xexp_pc[i][j] * tmpf[i][j][k] +
                       xexp_pd[i][j] * tmpf[i - 1][j][k];
      }

      f_t[lx][j][k] =
          xexp_pu[lx][j] * tmpf[lx + 1][j][k] + xexp_pc[lx][j] * tmpf[lx][j][k];

      f_t[ux][j][k] =
          xexp_pc[ux][j] * tmpf[ux][j][k] + xexp_pd[ux][j] * tmpf[ux - 1][j][k];
    }
  }
}

/*	Convolve one step backward using ADI scheme */
void num_f_pde_one_step_backward_2f_adi_bar(
    CNPDE_TEMP_2D_ADI *tmp,
    /*	Size of the X array	*/
    int nx,
    /*	X values	*/
    double *x,
    /*	Size of the Y array	*/
    int ny,
    /*	X values	*/
    double *y,
    /*	Number of functions	*/
    int mstart, int mend,
    /*	f values @ t+1	*/
    /*	i: x state  , j: y state  , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy, int is_up, int is_down) {
  static int i, j, k;
  static double dux, ddx, duy, ddy, bet;
  static double *muxi, *varxi, *ri, **fi, **tmpfi, *xexp_pui, *xexp_pci,
      *xexp_pdi, *ximp_pui, *ximp_pci, *ximp_pdi, *yexp_pui, *yexp_pci,
      *yexp_pdi, *yimp_pui, *yimp_pci, *yimp_pdi;
  static double temp1, temp2;

  double **xexp_pu = tmp->xexp_pu, **xexp_pc = tmp->xexp_pc,
         **xexp_pd = tmp->xexp_pd,

         **ximp_pu = tmp->ximp_pu, **ximp_pc = tmp->ximp_pc,
         **ximp_pd = tmp->ximp_pd,

         **yexp_pu = tmp->yexp_pu, **yexp_pc = tmp->yexp_pc,
         **yexp_pd = tmp->yexp_pd,

         **yimp_pu = tmp->yimp_pu, **yimp_pc = tmp->yimp_pc,
         **yimp_pd = tmp->yimp_pd,

         ***tmpf = tmp->f,

         *gamx = tmp->gamx, *gamy = tmp->gamy;

  /*	Calculate all probabilities */

  /*	Interior probas on y */
  for (j = ly + 1; j <= uy - 1; j++) {
    duy = y[j + 1] - y[j];
    ddy = y[j] - y[j - 1];
    temp1 = 0.5 / (duy * (duy + ddy));

    for (i = lx; i <= ux; i++) {
      yexp_pu[i][j] = temp1 * (vary[i][j] + muy[i][j] * ddy);
      yexp_pd[i][j] = (yexp_pu[i][j] * duy - 0.5 * muy[i][j]) / ddy;
      yexp_pc[i][j] = -yexp_pu[i][j] - yexp_pd[i][j];

      yimp_pu[i][j] = -yexp_pu[i][j];
      yimp_pd[i][j] = -yexp_pd[i][j];
      yimp_pc[i][j] = -yexp_pc[i][j];

      temp2 = 0.25 * r[i][j];
      yexp_pc[i][j] += 1.0 - temp2;
      yimp_pc[i][j] += 1.0 + temp2;
    }
  }

  /*	Down probas on y */
  duy = y[ly + 1] - y[ly];

  for (i = lx; i <= ux; i++) {
    yexp_pd[i][ly] = 0.0;
    yexp_pc[i][ly] = -0.5 * muy[i][ly] / duy;
    yexp_pu[i][ly] = -yexp_pc[i][ly];

    yimp_pu[i][ly] = -yexp_pu[i][ly];
    yimp_pd[i][ly] = -yexp_pd[i][ly];
    yimp_pc[i][ly] = -yexp_pc[i][ly];

    temp2 = 0.25 * r[i][ly];
    yexp_pc[i][ly] += 1.0 - temp2;
    yimp_pc[i][ly] += 1.0 + temp2;
  }

  /*	Up probas on y */
  ddy = y[uy] - y[uy - 1];

  for (i = lx; i <= ux; i++) {
    yexp_pu[i][uy] = 0.0;
    yexp_pc[i][uy] = 0.5 * muy[i][uy] / ddy;
    yexp_pd[i][uy] = -yexp_pc[i][uy];

    yimp_pu[i][uy] = -yexp_pu[i][uy];
    yimp_pd[i][uy] = -yexp_pd[i][uy];
    yimp_pc[i][uy] = -yexp_pc[i][uy];

    temp2 = 0.25 * r[i][uy];
    yexp_pc[i][uy] += 1.0 - temp2;
    yimp_pc[i][uy] += 1.0 + temp2;
  }

  /*	Interior probas on x */
  for (i = lx + 1; i <= ux - 1; i++) {
    dux = x[i + 1] - x[i];
    ddx = x[i] - x[i - 1];
    temp1 = 0.5 / (dux * (dux + ddx));

    muxi = mux[i];
    varxi = varx[i];
    ri = r[i];

    xexp_pui = xexp_pu[i];
    xexp_pci = xexp_pc[i];
    xexp_pdi = xexp_pd[i];
    ximp_pui = ximp_pu[i];
    ximp_pci = ximp_pc[i];
    ximp_pdi = ximp_pd[i];

    for (j = ly; j <= uy; j++) {
      xexp_pui[j] = temp1 * (varxi[j] + muxi[j] * ddx);
      xexp_pdi[j] = (xexp_pui[j] * dux - 0.5 * muxi[j]) / ddx;
      xexp_pci[j] = -xexp_pui[j] - xexp_pdi[j];

      ximp_pui[j] = -xexp_pui[j];
      ximp_pdi[j] = -xexp_pdi[j];
      ximp_pci[j] = -xexp_pci[j];

      temp2 = 0.25 * ri[j];
      xexp_pci[j] += 1.0 - temp2;
      ximp_pci[j] += 1.0 + temp2;
    }
  }

  /*	Down probas on x */
  dux = x[lx + 1] - x[lx];

  muxi = mux[lx];
  varxi = varx[lx];
  ri = r[lx];

  xexp_pui = xexp_pu[lx];
  xexp_pci = xexp_pc[lx];
  xexp_pdi = xexp_pd[lx];
  ximp_pui = ximp_pu[lx];
  ximp_pci = ximp_pc[lx];
  ximp_pdi = ximp_pd[lx];

  if (is_down) {
    yexp_pui = yexp_pu[lx];
    yexp_pci = yexp_pc[lx];
    yexp_pdi = yexp_pd[lx];
    yimp_pui = yimp_pu[lx];
    yimp_pci = yimp_pc[lx];
    yimp_pdi = yimp_pd[lx];

    for (j = ly; j <= uy; j++) {
      xexp_pdi[j] = 0.0;
      xexp_pui[j] = 0.0;
      xexp_pci[j] = 1.0;

      ximp_pui[j] = 0.0;
      ximp_pdi[j] = 0.0;
      ximp_pci[j] = 1.0;

      yexp_pdi[j] = 0.0;
      yexp_pui[j] = 0.0;
      yexp_pci[j] = 1.0;

      yimp_pui[j] = 0.0;
      yimp_pdi[j] = 0.0;
      yimp_pci[j] = 1.0;
    }
  } else {
    for (j = ly; j <= uy; j++) {
      xexp_pdi[j] = 0.0;
      xexp_pci[j] = -0.5 * muxi[j] / dux;
      xexp_pui[j] = -xexp_pci[j];

      ximp_pui[j] = -xexp_pui[j];
      ximp_pdi[j] = -xexp_pdi[j];
      ximp_pci[j] = -xexp_pci[j];

      temp2 = 0.25 * ri[j];
      xexp_pci[j] += 1.0 - temp2;
      ximp_pci[j] += 1.0 + temp2;
    }
  }

  /*	Up probas on x */
  ddx = x[ux] - x[ux - 1];

  muxi = mux[ux];
  varxi = varx[ux];
  ri = r[ux];

  xexp_pui = xexp_pu[ux];
  xexp_pci = xexp_pc[ux];
  xexp_pdi = xexp_pd[ux];
  ximp_pui = ximp_pu[ux];
  ximp_pci = ximp_pc[ux];
  ximp_pdi = ximp_pd[ux];

  if (is_up) {
    yexp_pui = yexp_pu[ux];
    yexp_pci = yexp_pc[ux];
    yexp_pdi = yexp_pd[ux];
    yimp_pui = yimp_pu[ux];
    yimp_pci = yimp_pc[ux];
    yimp_pdi = yimp_pd[ux];

    for (j = ly; j <= uy; j++) {
      xexp_pui[j] = 0.0;
      xexp_pdi[j] = 0.0;
      xexp_pci[j] = 1.0;

      ximp_pui[j] = 0.0;
      ximp_pdi[j] = 0.0;
      ximp_pci[j] = 1.0;

      yexp_pui[j] = 0.0;
      yexp_pdi[j] = 0.0;
      yexp_pci[j] = 1.0;

      yimp_pui[j] = 0.0;
      yimp_pdi[j] = 0.0;
      yimp_pci[j] = 1.0;
    }
  } else {
    for (j = ly; j <= uy; j++) {
      xexp_pui[j] = 0.0;
      xexp_pci[j] = 0.5 * muxi[j] / ddx;
      xexp_pdi[j] = -xexp_pci[j];

      ximp_pui[j] = -xexp_pui[j];
      ximp_pdi[j] = -xexp_pdi[j];
      ximp_pci[j] = -xexp_pci[j];

      temp2 = 0.25 * ri[j];
      xexp_pci[j] += 1.0 - temp2;
      ximp_pci[j] += 1.0 + temp2;
    }
  }

  /*	For each slice in y  , implicit in x */

  for (j = ly; j <= uy; j++) {
    for (k = mstart; k <= mend; k++) {
      f_t[lx][j][k] = f_t_plus_1[lx][j][k] / (bet = ximp_pc[lx][j]);
    }

    for (i = lx + 1; i <= ux; i++) {
      gamx[i] = ximp_pu[i - 1][j] / bet;
      bet = ximp_pc[i][j] - ximp_pd[i][j] * gamx[i];
      for (k = mstart; k <= mend; k++) {
        f_t[i][j][k] =
            (f_t_plus_1[i][j][k] - ximp_pd[i][j] * f_t[i - 1][j][k]) / bet;
      }
    }

    for (i = ux - 1; i >= lx; i--) {
      for (k = mstart; k <= mend; k++) {
        f_t[i][j][k] -= gamx[i + 1] * f_t[i + 1][j][k];
      }
    }
  }

  /*	For each slice in x  , explicit in y */

  for (i = lx; i <= ux; i++) {
    fi = f_t[i];
    tmpfi = tmpf[i];
    yexp_pui = yexp_pu[i];
    yexp_pci = yexp_pc[i];
    yexp_pdi = yexp_pd[i];

    for (k = mstart; k <= mend; k++) {
      for (j = ly + 1; j <= uy - 1; j++) {
        tmpfi[j][k - mstart] = yexp_pui[j] * fi[j + 1][k] +
                               yexp_pci[j] * fi[j][k] +
                               yexp_pdi[j] * fi[j - 1][k];
      }

      tmpfi[ly][k - mstart] =
          yexp_pui[ly] * fi[ly + 1][k] + yexp_pci[ly] * fi[ly][k];

      tmpfi[uy][k - mstart] =
          yexp_pci[uy] * fi[uy][k] + yexp_pdi[uy] * fi[uy - 1][k];
    }
  }

  /*	For each slice in x  , implicit in y */

  for (i = lx; i <= ux; i++) {
    tmpfi = tmpf[i];
    yimp_pui = yimp_pu[i];
    yimp_pci = yimp_pc[i];
    yimp_pdi = yimp_pd[i];

    for (k = mstart; k <= mend; k++) {
      tmpfi[ly][k - mstart] /= (bet = yimp_pci[ly]);
    }

    for (j = ly + 1; j <= uy; j++) {
      gamy[j] = yimp_pui[j - 1] / bet;
      bet = yimp_pci[j] - yimp_pdi[j] * gamy[j];
      for (k = mstart; k <= mend; k++) {
        tmpfi[j][k - mstart] =
            (tmpfi[j][k - mstart] - yimp_pdi[j] * tmpfi[j - 1][k - mstart]) /
            bet;
      }
    }

    for (j = uy - 1; j >= ly; j--) {
      for (k = mstart; k <= mend; k++) {
        tmpfi[j][k - mstart] -= gamy[j + 1] * tmpfi[j + 1][k - mstart];
      }
    }
  }

  /*	For each slice in y  , explicit in x */

  for (j = ly; j <= uy; j++) {
    for (k = mstart; k <= mend; k++) {
      for (i = lx + 1; i <= ux - 1; i++) {
        f_t[i][j][k] = xexp_pu[i][j] * tmpf[i + 1][j][k - mstart] +
                       xexp_pc[i][j] * tmpf[i][j][k - mstart] +
                       xexp_pd[i][j] * tmpf[i - 1][j][k - mstart];
      }

      f_t[lx][j][k] = xexp_pu[lx][j] * tmpf[lx + 1][j][k - mstart] +
                      xexp_pc[lx][j] * tmpf[lx][j][k - mstart];

      f_t[ux][j][k] = xexp_pc[ux][j] * tmpf[ux][j][k - mstart] +
                      xexp_pd[ux][j] * tmpf[ux - 1][j][k - mstart];
    }
  }
}

/*	Call before PDE	2D */
void num_f_pde_init_2d_LGMSV_adi(CNPDE_TEMP_2D_LGMSV_ADI *tmp,
                                 /*	Size of the arrays	*/
                                 int npsi, int nx, int ny,
                                 /*	Number of functions	*/
                                 int m) {
  tmp->xexp_pu = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->xexp_pc = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->xexp_pd = dmatrix(0, nx - 1, 0, ny - 1);

  tmp->ximp_pu = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->ximp_pc = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->ximp_pd = dmatrix(0, nx - 1, 0, ny - 1);

  tmp->yexp_pu = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->yexp_pc = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->yexp_pd = dmatrix(0, nx - 1, 0, ny - 1);

  tmp->yimp_pu = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->yimp_pc = dmatrix(0, nx - 1, 0, ny - 1);
  tmp->yimp_pd = dmatrix(0, nx - 1, 0, ny - 1);

  tmp->f = f4tensor(0, nx - 1, 0, ny - 1, 0, npsi - 1, 0, m - 1);

  tmp->gamx = (double *)calloc(nx, sizeof(double));
  tmp->gamy = (double *)calloc(ny, sizeof(double));
}

/*	Call after PDE 2D */
void num_f_pde_free_2d_LGMSV_adi(CNPDE_TEMP_2D_LGMSV_ADI *tmp,
                                 /*	Size of the arrays	*/
                                 int npsi, int nx, int ny,
                                 /*	Number of functions	*/
                                 int m) {
  free_dmatrix(tmp->xexp_pu, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->xexp_pc, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->xexp_pd, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->ximp_pu, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->ximp_pc, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->ximp_pd, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yexp_pu, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yexp_pc, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yexp_pd, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yimp_pu, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yimp_pc, 0, nx - 1, 0, ny - 1);
  free_dmatrix(tmp->yimp_pd, 0, nx - 1, 0, ny - 1);

  free_f4tensor(tmp->f, 0, nx - 1, 0, ny - 1, 0, npsi - 1, 0, m - 1);

  free(tmp->gamx);
  free(tmp->gamy);
}

/*	Convolve one step backward using ADI scheme */
void num_f_pde_one_step_backward_2f_save_LGMSV_adi(
    CNPDE_TEMP_2D_LGMSV_ADI *tmp,
    /*	Size of the X array	*/
    int nx,
    /*	X values	*/
    double *x,
    /*	Size of the Y array	*/
    int ny,
    /*	X values	*/
    double *y,
    /*	Number of functions	*/
    int mstart, int mend,
    /*	f values @ t+1	*/
    /*	i: x state  , j: y state  , k: function number */
    double ****f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually  , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ****f_t,
    /*	Lower and upper Psi		*/
    int lpsi, int upsi,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy) {
  static int i, j, l, k;
  static double dux, ddx, duy, ddy, bet;
  static double *muxi, *varxi, *ri, *xexp_pui, *xexp_pci, *xexp_pdi, *ximp_pui,
      *ximp_pci, *ximp_pdi, *yexp_pui, *yexp_pci, *yexp_pdi, *yimp_pui,
      *yimp_pci, *yimp_pdi;

  static double ***fi, ***fi1, ***fi2, ***fi3, **fij, **fij1, **fij2, **fij3,
      ***fi_plus1, **fij_plus1, ***tmpfi, ***tmpfi1, **tmpfij, **tmpfij1,
      **tmpfij2, **tmpfij3;

  static double temp1, temp2;

  double **xexp_pu = tmp->xexp_pu, **xexp_pc = tmp->xexp_pc,
         **xexp_pd = tmp->xexp_pd,

         **ximp_pu = tmp->ximp_pu, **ximp_pc = tmp->ximp_pc,
         **ximp_pd = tmp->ximp_pd,

         **yexp_pu = tmp->yexp_pu, **yexp_pc = tmp->yexp_pc,
         **yexp_pd = tmp->yexp_pd,

         **yimp_pu = tmp->yimp_pu, **yimp_pc = tmp->yimp_pc,
         **yimp_pd = tmp->yimp_pd,

         ****tmpf = tmp->f,

         *gamx = tmp->gamx, *gamy = tmp->gamy;

  /*	Calculate all probabilities */

  /*	Interior probas on y */
  for (j = ly + 1; j <= uy - 1; j++) {
    duy = y[j + 1] - y[j];
    ddy = y[j] - y[j - 1];
    temp1 = 0.5 / (duy * (duy + ddy));

    for (i = lx; i <= ux; i++) {
      yexp_pu[i][j] = temp1 * (vary[i][j] + muy[i][j] * ddy);
      yexp_pd[i][j] = (yexp_pu[i][j] * duy - 0.5 * muy[i][j]) / ddy;
      yexp_pc[i][j] = -yexp_pu[i][j] - yexp_pd[i][j];

      yimp_pu[i][j] = -yexp_pu[i][j];
      yimp_pd[i][j] = -yexp_pd[i][j];
      yimp_pc[i][j] = -yexp_pc[i][j];

      temp2 = 0.25 * r[i][j];
      yexp_pc[i][j] += 1.0 - temp2;
      yimp_pc[i][j] += 1.0 + temp2;
    }
  }

  /*	Down probas on y */
  duy = y[ly + 1] - y[ly];

  for (i = lx; i <= ux; i++) {
    yexp_pd[i][ly] = 0.0;
    yexp_pc[i][ly] = -0.5 * muy[i][ly] / duy;
    yexp_pu[i][ly] = -yexp_pc[i][ly];

    yimp_pu[i][ly] = -yexp_pu[i][ly];
    yimp_pd[i][ly] = -yexp_pd[i][ly];
    yimp_pc[i][ly] = -yexp_pc[i][ly];

    temp2 = 0.25 * r[i][ly];
    yexp_pc[i][ly] += 1.0 - temp2;
    yimp_pc[i][ly] += 1.0 + temp2;
  }

  /*	Up probas on y */
  ddy = y[uy] - y[uy - 1];

  for (i = lx; i <= ux; i++) {
    yexp_pu[i][uy] = 0.0;
    yexp_pc[i][uy] = 0.5 * muy[i][uy] / ddy;
    yexp_pd[i][uy] = -yexp_pc[i][uy];

    yimp_pu[i][uy] = -yexp_pu[i][uy];
    yimp_pd[i][uy] = -yexp_pd[i][uy];
    yimp_pc[i][uy] = -yexp_pc[i][uy];

    temp2 = 0.25 * r[i][uy];
    yexp_pc[i][uy] += 1.0 - temp2;
    yimp_pc[i][uy] += 1.0 + temp2;
  }

  /*	Interior probas on x */
  for (i = lx + 1; i <= ux - 1; i++) {
    dux = x[i + 1] - x[i];
    ddx = x[i] - x[i - 1];
    temp1 = 0.5 / (dux * (dux + ddx));

    muxi = mux[i];
    varxi = varx[i];
    ri = r[i];
    xexp_pui = xexp_pu[i];
    xexp_pci = xexp_pc[i];
    xexp_pdi = xexp_pd[i];
    ximp_pui = ximp_pu[i];
    ximp_pci = ximp_pc[i];
    ximp_pdi = ximp_pd[i];

    for (j = ly; j <= uy; j++) {
      xexp_pui[j] = temp1 * (varxi[j] + muxi[j] * ddx);
      xexp_pdi[j] = (xexp_pui[j] * dux - 0.5 * muxi[j]) / ddx;
      xexp_pci[j] = -xexp_pui[j] - xexp_pdi[j];

      ximp_pui[j] = -xexp_pui[j];
      ximp_pdi[j] = -xexp_pdi[j];
      ximp_pci[j] = -xexp_pci[j];

      temp2 = 0.25 * ri[j];
      xexp_pci[j] += 1.0 - temp2;
      ximp_pci[j] += 1.0 + temp2;
    }
  }

  /*	Down probas on x */
  dux = x[lx + 1] - x[lx];

  muxi = mux[lx];
  varxi = varx[lx];
  ri = r[lx];
  xexp_pui = xexp_pu[lx];
  xexp_pci = xexp_pc[lx];
  xexp_pdi = xexp_pd[lx];
  ximp_pui = ximp_pu[lx];
  ximp_pci = ximp_pc[lx];
  ximp_pdi = ximp_pd[lx];

  for (j = ly; j <= uy; j++) {
    xexp_pdi[j] = 0.0;
    xexp_pci[j] = -0.5 * muxi[j] / dux;
    xexp_pui[j] = -xexp_pci[j];

    ximp_pui[j] = -xexp_pui[j];
    ximp_pdi[j] = -xexp_pdi[j];
    ximp_pci[j] = -xexp_pci[j];

    temp2 = 0.25 * ri[j];
    xexp_pci[j] += 1.0 - temp2;
    ximp_pci[j] += 1.0 + temp2;
  }

  /*	Up probas on x */
  ddx = x[ux] - x[ux - 1];

  muxi = mux[ux];
  varxi = varx[ux];
  ri = r[ux];
  xexp_pui = xexp_pu[ux];
  xexp_pci = xexp_pc[ux];
  xexp_pdi = xexp_pd[ux];
  ximp_pui = ximp_pu[ux];
  ximp_pci = ximp_pc[ux];
  ximp_pdi = ximp_pd[ux];

  for (j = ly; j <= uy; j++) {
    xexp_pui[j] = 0.0;
    xexp_pci[j] = 0.5 * muxi[j] / ddx;
    xexp_pdi[j] = -xexp_pci[j];

    ximp_pui[j] = -xexp_pui[j];
    ximp_pdi[j] = -xexp_pdi[j];
    ximp_pci[j] = -xexp_pci[j];

    temp2 = 0.25 * ri[j];
    xexp_pci[j] += 1.0 - temp2;
    ximp_pci[j] += 1.0 + temp2;
  }

#if 0

	/*	For each slice in y  , implicit in x */

	for (j=ly; j<=uy; j++)
	{
		bet = ximp_pc[lx][j];
		
		for (k=lpsi; k<=upsi; k++)
		{
			for (l=mstart; l<=mend; l++) 
			{
				f_t[lx][j][k][l] = f_t_plus_1[lx][j][k][l] / bet;
			}
		}

		for (i=lx+1; i<=ux; i++) 
		{
			gamx[i] = ximp_pu[i-1][j] / bet;
			bet = ximp_pc[i][j] - ximp_pd[i][j] * gamx[i];

			for (k=lpsi; k<=upsi; k++)
			{
				for (l=mstart; l<=mend; l++)
				{
					f_t[i][j][k][l] = (f_t_plus_1[i][j][k][l] - ximp_pd[i][j] * f_t[i-1][j][k][l]) / bet;
				}
			}
		}

		for (i=ux-1; i>=lx; i--)
		{			
			for (k=lpsi; k<=upsi; k++)
			{
				for (l=mstart; l<=mend; l++)
				{
					f_t[i][j][k][l] -= gamx[i+1] * f_t[i+1][j][k][l];
				}
			}
		}
	}

	/*	For each slice in x  , explicit in y */

	for (i=lx; i<=ux; i++)
	{
		fi = f_t[i];
		tmpfi = tmpf[i];

		yexp_pui = yexp_pu[i];
		yexp_pci = yexp_pc[i];
		yexp_pdi = yexp_pd[i];
		
		for (j=ly+1; j<=uy-1; j++)
		{
			fij = fi[j];
			fij1 = fi[j+1];
			fij2 = fi[j-1];
			tmpfij = tmpfi[j];

			for (k=lpsi; k<=upsi; k++)
			{
				for (l=mstart; l<=mend; l++)
				{
					tmpfij[k][l] = yexp_pui[j] * fij1[k][l]
								+ yexp_pci[j] * fij[k][l]
								+ yexp_pdi[j] * fij2[k][l];
				}
			}
		}

		for (k=lpsi; k<=upsi; k++)
		{
			for (l=mstart; l<=mend; l++)
			{
				tmpfi[ly][k][l] = yexp_pui[ly] * fi[ly+1][k][l]
							+ yexp_pci[ly] * fi[ly][k][l];
				
				tmpfi[uy][k][l] = yexp_pci[uy] * fi[uy][k][l]
								+ yexp_pdi[uy] * fi[uy-1][k][l];
			}
		}
	}
	
	/*	For each slice in x  , implicit in y */

	for (i=lx; i<=ux; i++)
	{
		tmpfi = tmpf[i];
		tmpfij = tmpfi[ly];

		yimp_pui = yimp_pu[i];
		yimp_pci = yimp_pc[i];
		yimp_pdi = yimp_pd[i];

		bet = yimp_pci[ly];

		for (k=lpsi; k<=upsi; k++)
		{
			for (l=mstart; l<=mend; l++) 
			{
				tmpfij[k][l] /= bet;
			}
		}

		for (j=ly+1; j<=uy; j++) 
		{
			tmpfij = tmpfi[j];
			tmpfij1 = tmpfi[j-1];

			gamy[j] = yimp_pui[j-1] / bet;
			bet = yimp_pci[j] - yimp_pdi[j] * gamy[j];
			
			for (k=lpsi; k<=upsi; k++)
			{
				for (l=mstart; l<=mend; l++)
				{
					tmpfij[k][l] = (tmpfij[k][l] - yimp_pdi[j] * tmpfij1[k][l]) / bet;
				}
			}
		}

		for (j=uy-1; j>=ly; j--)
		{
			tmpfij = tmpfi[j];
			tmpfij1 = tmpfi[j+1];

			for (k=lpsi; k<=upsi; k++)
			{
				for (l=mstart; l<=mend; l++)
				{				
					tmpfij[k][l] -= gamy[j+1] * tmpfij1[k][l];
				}
			}
		}
	}

	/*	For each slice in y  , explicit in x */		

	for (j=ly; j<=uy; j++)
	{
		fi = f_t[i];
		tmpfi = tmpf[i];

		for (k=lpsi; k<=upsi; k++)
		{
			for (l=mstart; l<=mend; l++)
			{
				for (i=lx+1; i<=ux-1; i++)
				{
					f_t[i][j][k][l] =	xexp_pu[i][j] * tmpf[i+1][j][k][l]
									+ xexp_pc[i][j] * tmpf[i][j][k][l]
									+ xexp_pd[i][j] * tmpf[i-1][j][k][l];
				}

				f_t[lx][j][k][l] = xexp_pu[lx][j] * tmpf[lx+1][j][k][l]
							+ xexp_pc[lx][j] * tmpf[lx][j][k][l];
				
				f_t[ux][j][k][l] = xexp_pc[ux][j] * tmpf[ux][j][k][l]
							+ xexp_pd[ux][j] * tmpf[ux-1][j][k][l];
			}
		}
	}

#else

  /*	For each slice in x  , implicit in y */

  for (i = lx; i <= ux; i++) {
    fi = f_t[i];
    fi_plus1 = f_t_plus_1[i];

    fij = fi[ly];
    fij_plus1 = fi_plus1[ly];

    yimp_pui = yimp_pu[i];
    yimp_pci = yimp_pc[i];
    yimp_pdi = yimp_pd[i];

    bet = yimp_pci[ly];

    for (k = lpsi; k <= upsi; k++) {
      for (l = mstart; l <= mend; l++) {
        fij[k][l] = fij_plus1[k][l] / bet;
      }
    }

    for (j = ly + 1; j <= uy; j++) {
      fij = fi[j];
      fij1 = fi[j - 1];
      fij_plus1 = fi_plus1[j];

      gamy[j] = yimp_pui[j - 1] / bet;
      bet = yimp_pci[j] - yimp_pdi[j] * gamy[j];

      for (k = lpsi; k <= upsi; k++) {
        for (l = mstart; l <= mend; l++) {
          fij[k][l] = (fij_plus1[k][l] - yimp_pdi[j] * fij1[k][l]) / bet;
        }
      }
    }

    for (j = uy - 1; j >= ly; j--) {
      fij = fi[j];
      fij1 = fi[j + 1];

      for (k = lpsi; k <= upsi; k++) {
        for (l = mstart; l <= mend; l++) {
          fij[k][l] -= gamy[j + 1] * fij1[k][l];
        }
      }
    }
  }

  /*	For each slice in y  , explicit in x */

  for (i = lx + 1; i <= ux - 1; i++) {
    fi = f_t[i];
    fi1 = f_t[i + 1];
    fi2 = f_t[i - 1];
    tmpfi = tmpf[i];

    for (j = ly; j <= uy; j++) {
      fij = fi[j];
      fij1 = fi1[j];
      fij2 = fi2[j];
      tmpfij = tmpfi[j];

      for (k = lpsi; k <= upsi; k++) {
        for (l = mstart; l <= mend; l++) {
          tmpfij[k][l] = xexp_pu[i][j] * fi1[j][k][l] +
                         xexp_pc[i][j] * fij[k][l] + xexp_pd[i][j] * fij2[k][l];
        }
      }
    }
  }

  fi = f_t[lx];
  fi1 = f_t[lx + 1];
  fi2 = f_t[ux];
  fi3 = f_t[ux - 1];

  tmpfi = tmpf[lx];
  tmpfi1 = tmpf[ux];

  for (j = ly; j <= uy; j++) {
    fij = fi[j];
    fij1 = fi1[j];
    fij2 = fi2[j];
    fij3 = fi3[j];

    tmpfij = tmpfi[j];
    tmpfij1 = tmpfi1[j];

    for (k = lpsi; k <= upsi; k++) {
      for (l = mstart; l <= mend; l++) {
        tmpfij[k][l] = xexp_pu[lx][j] * fij1[k][l] + xexp_pc[lx][j] * fij[k][l];

        tmpfij1[k][l] =
            xexp_pc[ux][j] * fij2[k][l] + xexp_pd[ux][j] * fij3[k][l];
      }
    }
  }

  /*	For each slice in y  , implicit in x */

  tmpfi = tmpf[lx];

  for (j = ly; j <= uy; j++) {
    tmpfij = tmpfi[j];

    bet = ximp_pc[lx][j];

    for (k = lpsi; k <= upsi; k++) {
      for (l = mstart; l <= mend; l++) {
        tmpfij[k][l] /= bet;
      }
    }

    for (i = lx + 1; i <= ux; i++) {
      tmpfij = tmpf[i][j];
      tmpfij1 = tmpf[i - 1][j];

      gamx[i] = ximp_pu[i - 1][j] / bet;
      bet = ximp_pc[i][j] - ximp_pd[i][j] * gamx[i];

      for (k = lpsi; k <= upsi; k++) {
        for (l = mstart; l <= mend; l++) {
          tmpfij[k][l] =
              (tmpf[i][j][k][l] - ximp_pd[i][j] * tmpfij1[k][l]) / bet;
        }
      }
    }

    for (i = ux - 1; i >= lx; i--) {
      tmpfij = tmpf[i][j];
      tmpfij1 = tmpf[i + 1][j];

      for (k = lpsi; k <= upsi; k++) {
        for (l = mstart; l <= mend; l++) {
          tmpfij[k][l] -= gamx[i + 1] * tmpfij1[k][l];
        }
      }
    }
  }

  /*	For each slice in x  , explicit in y */

  for (i = lx; i <= ux; i++) {
    fi = f_t[i];
    tmpfi = tmpf[i];
    yexp_pui = yexp_pu[i];
    yexp_pci = yexp_pc[i];
    yexp_pdi = yexp_pd[i];

    for (j = ly + 1; j <= uy - 1; j++) {
      fij = fi[j];
      tmpfij = tmpfi[j];
      tmpfij1 = tmpfi[j + 1];
      tmpfij2 = tmpfi[j - 1];

      for (k = lpsi; k <= upsi; k++) {
        for (l = mstart; l <= mend; l++) {
          fij[k][l] = yexp_pui[j] * tmpfij1[k][l] + yexp_pci[j] * tmpfij[k][l] +
                      yexp_pdi[j] * tmpfij2[k][l];
        }
      }
    }

    fij = fi[ly];
    fij1 = fi[uy];

    tmpfij = tmpfi[ly];
    tmpfij1 = tmpfi[ly + 1];
    tmpfij2 = tmpfi[uy];
    tmpfij3 = tmpfi[uy - 1];

    for (k = lpsi; k <= upsi; k++) {
      for (l = mstart; l <= mend; l++) {
        fij[k][l] = yexp_pui[ly] * tmpfij1[k][l] + yexp_pci[ly] * tmpfij[k][l];

        fij1[k][l] =
            yexp_pci[uy] * tmpfij2[k][l] + yexp_pdi[uy] * tmpfij3[k][l];
      }
    }
  }

#endif
}