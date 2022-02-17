
/* ==========================================================================
   FILE_NAME:	Fx3FBetaDLMTree.cxx

   PURPOSE:		Modelling of a 3 Facors Model:

                                - First one is LGM one factor on the Domestic
   Market        ,
                                - Second one is LGM one factor on the Foreign
   Market        ,
                                - Third one is Lognormal BS for the Spot Fx.

                                A volatility term structure is allowed for the
   three dynamics        , but for the LGM dynamics        , we assume a fixed
   tau. The correlations have also to be constant.

                                The pricing is done trough a trinomial tree.

   DATE:		08/01/00

   AUTHOR:		L.C.
   ========================================================================== */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"

/*	Global vectors and matrices */
double *temp_vec1, *temp_vec2, **temp_mat1, **temp_mat2, **temp_mat3;

static double Boundaries_Matrix[3][8] = {
    {1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0},
    {1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0},
    {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0}};

static double Zero_Vector[3] = {0.0, 0.0, 0.0};

/*	To calculate the boundaries at a given time t        , we first
 * calculate the global dual basis		*/
/*	according to the volatilities and expectations seen from today.
 */
/*	We apply NSTD standard deviation in this base (vol in this base in equal
 * to 1)				*/
/*	Eventually        , we put the resulting boundaries in the local dual
 * base (according to the		*/
/*	volatilities seen at date t).
 */

static Err calculate_boundaries(
    double sqt_inv, double *expect, double *var, double *correl,
    double **transf_matrix_inv, double **maxmin_dual, double **maxmin_orig,
    double *volt, double **maxmin_dual_global, double **maxmin_dual_local,
    double *expect_dual, double **transf_matrix_global,
    double **transf_matrix_inv_global, double **covar_matrix,
    double **eigen_vector, double **eigen_vector_transp, double *eigen_val) {
  Err err = NULL;
  static unsigned short i, j, gi, gj, gk;
  static double gsum, gmin, gmax;

  for (i = 0; i < 3; i++) {
    if (var[i] < 1.0e-16) {
      err =
          "Variance is NOT positive - check your volatilities and correlations";
      return err;
    }

    volt[i] = sqrt(var[i]);
  }

  /*	Get the global transformation matricies  */
  err = get_transf_matrix(volt, correl, transf_matrix_global,
                          transf_matrix_inv_global, covar_matrix, eigen_vector,
                          eigen_vector_transp, eigen_val);

  /*	Calculation of the expectation in the local dual base */
  prod_matrix_vector_constant(transf_matrix_inv, expect, sqt_inv, expect_dual,
                              gi, gk, gsum);

  /*	Calculate the deviations in the original base corresponding
          to -1 and +1 standard deviation in the global base
   */
  prod_matrix_special(transf_matrix_global, Boundaries_Matrix,
                      maxmin_dual_global, gi, gj, gk, gsum);

  /*	Find the corresponding deviation in the local dual base */
  prod_matrix_constant_special(transf_matrix_inv, maxmin_dual_global, sqt_inv,
                               maxmin_dual_local, gi, gj, gk, gsum);

  /*	Find in each direction the min and the max deviation        ,
          multiply it by NSTD        , and add it to the expected value
                                                                                                                                  */
  for (i = 0; i < 3; i++) {
    gmin = maxmin_dual_local[i][0];
    gmax = gmin;
    for (j = 1; j < 8; j++) {
      if (maxmin_dual_local[i][j] > gmax) {
        gmax = maxmin_dual_local[i][j];
      }
      if (maxmin_dual_local[i][j] < gmin) {
        gmin = maxmin_dual_local[i][j];
      }
    }
    maxmin_dual[i][0] = expect_dual[i] + NSTD * gmax;
    maxmin_dual[i][1] = expect_dual[i] + NSTD * gmin;
  }

  return err;
}

/*	Calculate the boundaries at each time step and save		*/
/*	the maximum of nodes that are used in each direction	*/
static Err calculate_all_boundaries(
    long nstp, double *time, int *vol_change, double *sig_dom, double *sig_for,
    double *sig_fx, double *dom_fwd, double *dom_var, double *for_fwd,
    double *for_var, double *fx_fwd, double *fx_var, double **correl,
    double **glob_correl, int *is_bar, double *bar_lvl, long *bar_k,
    long ***bar_idx, long nprod, long **nb_nodes_tab, long *nb_nodes_max,
    double ***maxmin_dual_tab, double *sigma, double **maxmin_orig,
    double **transf_matrix, double **transf_matrix_inv, double **covar_matrix,
    double **eigen_vector, double **eigen_vector_transp, double *eigen_val) {
  double **transf_matrix_global = NULL, **transf_matrix_inv_global = NULL,
         **maxmin_dual_global = NULL, **maxmin_dual_local = NULL,
         *expectation = NULL, *vart = NULL, *volt = NULL, *expect_dual = NULL,
         *transf_fx = NULL;

  double t, prev_t, sqt_prev_inv, k1_p;
  double a, b, c, k2_p;
  long step, i, j, k, k1;
  long idx_rel, idx1, idx2;
  Err err;

  /*	Memory allocation */
  transf_matrix_global = dmatrix(0, 2, 0, 2);
  transf_matrix_inv_global = dmatrix(0, 2, 0, 2);
  expectation = dvector(0, 2);
  vart = dvector(0, 2);
  maxmin_dual_global = dmatrix(0, 2, 0, 7);
  maxmin_dual_local = dmatrix(0, 2, 0, 7);
  expect_dual = dvector(0, 2);
  volt = dvector(0, 2);

  if (!(transf_matrix_global) || !(transf_matrix_inv_global) ||
      !(expectation) || !(vart) || !(volt) || !(maxmin_dual_global) ||
      !(maxmin_dual_local) || !(expect_dual)) {
    err = "Memory allocation failed in calculate_all_boundaries";
    goto FREE_RETURN_BOUNDARIES;
  }

  /*	Initialisation */
  t = time[nstp - 1];
  nb_nodes_max[0] = 0;
  nb_nodes_max[1] = 0;
  nb_nodes_max[2] = 0;

  /*	Loop on time steps */
  for (step = nstp - 1; step >= 1; step--) {
    prev_t = time[step - 1];
    sqt_prev_inv = 1.0 / (sqrt(t - prev_t));

    /*	Find transf matrix from prev_t to t */
    if (vol_change[step]) {
      /*	We need to recompute everything  */
      sigma[0] = sig_dom[step];
      sigma[1] = sig_for[step];
      sigma[2] = sig_fx[step];

      /*	Calculates transf_matrix and transf_matrix_inv */
      err = get_transf_matrix(sigma, correl[step], transf_matrix,
                              transf_matrix_inv, covar_matrix, eigen_vector,
                              eigen_vector_transp, eigen_val);

      transf_fx = transf_matrix[2];

      if (err) {
        goto FREE_RETURN_BOUNDARIES;
      }
    }

    /*	Find boundaries        , nbNodes        , maxmin_dual at date t */
    expectation[0] = dom_fwd[step];
    expectation[1] = for_fwd[step];
    expectation[2] = fx_fwd[step];
    vart[0] = dom_var[step];
    vart[1] = for_var[step];
    vart[2] = fx_var[step];

    err = calculate_boundaries(
        sqt_prev_inv, expectation, vart, glob_correl[step], transf_matrix_inv,
        maxmin_dual_tab[step], maxmin_orig, volt, maxmin_dual_global,
        maxmin_dual_local, expect_dual, transf_matrix_global,
        transf_matrix_inv_global, covar_matrix, eigen_vector,
        eigen_vector_transp, eigen_val);

    if (err) {
      goto FREE_RETURN_BOUNDARIES;
    }

    for (i = 0; i < 3; i++) {
      nb_nodes_tab[step][i] =
          (long)((maxmin_dual_tab[step][i][0] - maxmin_dual_tab[step][i][1]) /
                 MESH) +
          1;

      maxmin_dual_tab[step][i][0] =
          maxmin_dual_tab[step][i][1] + (nb_nodes_tab[step][i] - 1) * MESH;

      /*	Check if we update the maximum of nodes */
      if (nb_nodes_tab[step][i] > nb_nodes_max[i]) {
        nb_nodes_max[i] = nb_nodes_tab[step][i];
      }
      if (nb_nodes_tab[step][i] > MAXNODE3D) {
        err = serror("MAXNODE3D exceeded in tree_main_3dfx: %d > %d !!!",
                     nb_nodes_tab[step][i], MAXNODE3D);
        goto FREE_RETURN_BOUNDARIES;
      }
    }

    /*	Barrier stuff
            -------------	*/

    /*	1)	Determine which factor is the most relevant to Fx */

    if (vol_change[step]) {
      if (transf_fx[0] >= transf_fx[1]) {
        if (transf_fx[0] >= transf_fx[2]) {
          bar_k[step] = 0;
        } else {
          bar_k[step] = 2;
        }
      } else if (transf_fx[1] >= transf_fx[2]) {
        bar_k[step] = 1;
      } else {
        bar_k[step] = 2;
      }
    } else {
      bar_k[step] = bar_k[step + 1];
    }

    switch (bar_k[step]) {
    case 0:
      idx_rel = 0;
      idx1 = 1;
      idx2 = 2;
      break;
    case 1:
      idx_rel = 1;
      idx1 = 0;
      idx2 = 2;
      break;
    case 2:
      idx_rel = 2;
      idx1 = 0;
      idx2 = 1;
      break;
    default:
      break;
    }

    /*	2)	Modify the relevant min so as to force a node at the barrier
                            in the medium plane in (j        , k)
                            i.e. fx (i0) (j=n2/2        , k=n3/2) = barrier */

    a = (bar_lvl[step] * sqt_prev_inv -
         transf_fx[0] * maxmin_dual_tab[step][0][1] -
         transf_fx[1] * maxmin_dual_tab[step][1][1] -
         transf_fx[2] * maxmin_dual_tab[step][2][1]) /
        (transf_fx[idx_rel] * MESH);
    b = -transf_fx[idx1] / transf_fx[idx_rel];
    c = -transf_fx[idx2] / transf_fx[idx_rel];

    k1_p = a + b * ((int)(nb_nodes_tab[step][idx1] / 2 + 0.5) * 1.0) +
           c * ((int)(nb_nodes_tab[step][idx2] / 2 + 0.5) * 1.0);

    k1 = (int)(k1_p + 0.5);

    if (k1 >= 1 && k1 <= nb_nodes_tab[step][idx_rel] - 2) {
      maxmin_dual_tab[step][idx_rel][1] += (k1_p - k1) * MESH;
      maxmin_dual_tab[step][idx_rel][0] =
          maxmin_dual_tab[step][idx_rel][1] +
          (nb_nodes_tab[step][idx_rel] - 1) * MESH;
    } else {
      if (bar_lvl[step] == -100000000.00000) {
        is_bar[step] = 0;
      } else {
        if (k1 < 1) {
          is_bar[step] = -1; /* barrier is too low */
        } else {
          is_bar[step] = -2; /* barrier is too high */
        }
      }

      a = (fx_fwd[step] * sqt_prev_inv -
           transf_fx[0] * maxmin_dual_tab[step][0][1] -
           transf_fx[1] * maxmin_dual_tab[step][1][1] -
           transf_fx[2] * maxmin_dual_tab[step][2][1]) /
          (transf_fx[idx_rel] * MESH);
      b = -transf_fx[idx1] / transf_fx[idx_rel];
      c = -transf_fx[idx2] / transf_fx[idx_rel];

      k1_p = a + b * ((int)(nb_nodes_tab[step][idx1] / 2 + 0.5) * 1.0) +
             c * ((int)(nb_nodes_tab[step][idx2] / 2 + 0.5) * 1.0);

      k1 = (int)(k1_p + 0.5);

      maxmin_dual_tab[step][idx_rel][1] += (k1_p - k1) * MESH;
      maxmin_dual_tab[step][idx_rel][0] =
          maxmin_dual_tab[step][idx_rel][1] +
          (nb_nodes_tab[step][idx_rel] - 1) * MESH;
    }

    /*	3)	For each couple (j        , k) in the irrelevant directions        ,
       find the index i0 in the relevant direction such that the node i0 is the
       closest to the barrier	*/

    /*	Also check that 0 <= i0 <= n-2        , i.e. i0 is within the bounds */

    if (is_bar[step] == 1) {
      a -= (k1_p - k1);
      bar_idx[step] = lngmatrix(0, nb_nodes_tab[step][idx1] - 1, 0,
                                nb_nodes_tab[step][idx2] - 1);
      if (!bar_idx[step]) {
        err = "Memory allocation error in calculate_all_boundaries";
        goto FREE_RETURN_BOUNDARIES;
      }

      k1_p = a;
      k2_p = a;
      for (j = 0; j < nb_nodes_tab[step][idx1]; j++) {
        for (k = 0; k < nb_nodes_tab[step][idx2]; k++) {
          k1 = (long)(k1_p + 0.5);
          if (k1 <= 0 || k1 >= nb_nodes_tab[step][idx_rel] - 1) {
            k1 = -2;
          }
          bar_idx[step][j][k] = k1;

          k1_p += c;
        }
        k2_p += b;
        k1_p = k2_p;
      }
    }

    t = prev_t;
  }

FREE_RETURN_BOUNDARIES:

  /*	Free memory and return */

  if (maxmin_dual_global) {
    free_dmatrix(maxmin_dual_global, 0, 2, 0, 7);
  }

  if (maxmin_dual_local) {
    free_dmatrix(maxmin_dual_local, 0, 2, 0, 7);
  }

  if (transf_matrix_global) {
    free_dmatrix(transf_matrix_global, 0, 2, 0, 2);
  }

  if (transf_matrix_inv_global) {
    free_dmatrix(transf_matrix_inv_global, 0, 2, 0, 2);
  }

  if (expectation) {
    free_dvector(expectation, 0, 2);
  }

  if (vart) {
    free_dvector(vart, 0, 2);
  }

  if (volt) {
    free_dvector(volt, 0, 2);
  }

  if (expect_dual) {
    free_dvector(expect_dual, 0, 2);
  }

  return err;
}

/*	Main function under QBeta */
/*	------------------------- */
Err tree_main_3dBetaDLM_QBeta(
    /*	Time data */
    long nstp, double *time, double *date,

    /*	Model data */
    int *vol_change,

    double dom_lam, double for_lam,

    double *sig_dom, double *sig_for, double *sig_fx, double *corr_dom_for,
    double *corr_dom_fx, double *corr_for_fx,

    double *dr_const_dom, double *dr_coef_dom, double *dr_const_for,
    double *dr_coef1_for, double *dr_coef2_for, double *dr_coef3_for,
    double *dr_const_fx, double *dr_coef1_fx, double *dr_coef2_fx,
    double *dr_coef3_fx,

    /*	Distributions and Constants */
    double *dom_fwd, double *dom_var, double *for_fwd, double *for_var,
    double *fx_fwd, double *fx_var,

    double *glob_corr_dom_for, double *glob_corr_dom_fx,
    double *glob_corr_for_fx,

    /*	Market data */
    double spot_fx, char *dom_yc, char *for_yc, double *dom_ifr,

    /*	Product data */
    void **func_parm_tab, int *eval_evt, double *bar_lvl, int *bar_col,
    int *is_bar,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,

                       /* Market Data */
                       double spot_fx, char *dom_yc, char *for_yc,

                       /* Nodes data */
                       long n1, long n2, long n3, double ****sv,
                       /* Results */
                       long nprod, double ****prod_val,
                       /* Barrier details */
                       int is_bar, int bar_k, int **bar_idx, int bar_col,
                       double bar_lvl),
    /*	Result */
    int nprod, int discount, double *res) {

  long i, j, k, step, l, m, n;
  int r;
  double t, prev_t, next_t, sqt_prev, delta_t, sqt, sqt_inv;
  long gi, gj, gk;
  long index_fwd;
  double gsum, vol_fx;
  double F, min_value, gnu;
  double proba1, proba12, proba123;
  int flag;
  double rd;

  double ***link1 = NULL, **link2 = NULL, *link3 = NULL;

  double ****values = NULL, ****values_next = NULL, ****values_temp = NULL,
         ***valuesi = NULL, **valuesij = NULL, *valuesijk = NULL, ****sv = NULL,
         ***svi = NULL, **svij = NULL, *svijk = NULL, ****expect_dual = NULL,
         ***expect_duali = NULL, **expect_dualij = NULL, *expect_dualijk = NULL;

  double **transf_matrix = NULL, **transf_matrix_next = NULL,
         **transf_matrix_temp = NULL, **transf_matrix_inv = NULL,
         **transf_matrix_inv_next = NULL, **transf_matrix_inv_temp = NULL,
         **maxmin_dual = NULL, ***maxmin_dual_tab = NULL,
         **maxmin_dual_next = NULL, *sigma = NULL, *sigma_next = NULL,
         *volt = NULL, *sigma_temp = NULL;

  long **nb_nodes_tab = NULL, *nb_nodes = NULL, *nb_nodes_max = NULL,
       *nb_nodes_next = NULL;

  long *bar_k = NULL, ***bar_idx = NULL;

  double **correl = NULL, **glob_correl = NULL, **covar_matrix = NULL,
         **eigen_vector = NULL, **eigen_vector_transp = NULL, *eigen_val = NULL,
         **maxmin_orig = NULL, *expect_node = NULL, *expect_dual_node = NULL,
         **proba = NULL, *coos_dual = NULL;

  double temp1, temp2, temp3;

  double ca[3], cai[3], caj[3], cak[3], cajn[3], cakn[3];
  double ea[3], eai[3], eaj[3], eak[3], eajn[3], eakn[3];

  time_t t1, t2;

  size_t nprod_doubles = (nprod + 1) * sizeof(double);

  long fwd_coos[3], link_type[3];

  long n0, n1, n2;

  int do_dig;
  double dig_strike;
  double dig_dt;
  int dig_col;
  double dig_ivol;

  Err err = NULL;

  /*	Allocate memory once for all */

  t1 = clock();

  temp_vec1 = temp_vec2 = NULL;
  temp_mat1 = temp_mat2 = temp_mat3 = NULL;

  transf_matrix = dmatrix(0, 2, 0, 2);
  transf_matrix_next = dmatrix(0, 2, 0, 2);
  transf_matrix_inv = dmatrix(0, 2, 0, 2);
  transf_matrix_inv_next = dmatrix(0, 2, 0, 2);
  sigma = dvector(0, 2);
  sigma_next = dvector(0, 2);
  maxmin_dual_tab = f3tensor(0, nstp - 1, 0, 2, 0, 1);
  nb_nodes_tab = lngmatrix(0, nstp - 1, 0, 2);
  nb_nodes_max = lvector(0, 2);
  correl = dmatrix(0, nstp - 1, 0, 2);
  glob_correl = dmatrix(0, nstp - 1, 0, 2);
  covar_matrix = dmatrix(0, 2, 0, 2);
  eigen_vector = dmatrix(0, 2, 0, 2);
  eigen_vector_transp = dmatrix(0, 2, 0, 2);
  eigen_val = dvector(0, 2);
  maxmin_orig = dmatrix(0, 2, 0, 1);
  expect_node = dvector(0, 2);
  expect_dual_node = dvector(0, 2);
  proba = dmatrix(0, 2, 0, 2);
  coos_dual = dvector(0, 2);

  temp_vec1 = dvector(0, 2);
  temp_vec2 = dvector(0, 2);
  temp_mat1 = dmatrix(0, 2, 0, 2);
  temp_mat2 = dmatrix(0, 2, 0, 2);
  temp_mat3 = dmatrix(0, 2, 0, 2);

  if (!transf_matrix || !transf_matrix_next || !transf_matrix_inv ||
      !transf_matrix_inv_next || !maxmin_dual_tab || !sigma || !sigma_next ||
      !nb_nodes_tab || !nb_nodes_max || !correl || !glob_correl ||
      !covar_matrix || !eigen_vector || !eigen_vector_transp || !eigen_val ||
      !maxmin_orig || !expect_node || !expect_dual_node || !proba ||
      !coos_dual || !temp_vec1 || !temp_vec2 || !temp_mat1 || !temp_mat2 ||
      !temp_mat3) {
    err = "Memory allocation error (1) in tree_main_3dfx";
    goto FREE_RETURN;
  }

  /*	Calculates all boundaries and number of nodes */

  for (i = 0; i < nstp; i++) {
    correl[i][0] = corr_dom_for[i];
    correl[i][1] = corr_dom_fx[i];
    correl[i][2] = corr_for_fx[i];

    glob_correl[i][0] = glob_corr_dom_for[i];
    glob_correl[i][1] = glob_corr_dom_fx[i];
    glob_correl[i][2] = glob_corr_for_fx[i];
  }

  bar_k = lngvector(0, nstp - 1);
  bar_idx = (long ***)calloc(nstp, sizeof(long **));

  if (!bar_k || !bar_idx) {
    err = "Memory allocation (2) error in tree_main_3dfx";
    goto FREE_RETURN;
  }

  err = calculate_all_boundaries(
      nstp, time, vol_change, sig_dom, sig_for, sig_fx, dom_fwd, dom_var,
      for_fwd, for_var, fx_fwd, fx_var, correl, glob_correl, is_bar, bar_lvl,
      bar_k, bar_idx, nprod, nb_nodes_tab, nb_nodes_max, maxmin_dual_tab, sigma,
      maxmin_orig, transf_matrix, transf_matrix_inv, covar_matrix, eigen_vector,
      eigen_vector_transp, eigen_val);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Now we can allocate memory for the tensors */

  values = f4tensor(0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                    nb_nodes_max[2] - 1, 0, nprod);
  values_next = f4tensor(0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                         nb_nodes_max[2] - 1, 0, nprod);
  sv = f4tensor(0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                nb_nodes_max[2] - 1, 0, 2);
  expect_dual = f4tensor(0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                         nb_nodes_max[2] - 1, 0, 2);

  if (!values || !values_next || !sv || !expect_dual) {
    err = "Memory allocation (3) error in tree_main_3dfx";
    goto FREE_RETURN;
  }

  /*	Initialisation */

  transf_matrix_temp = transf_matrix;
  transf_matrix_inv_temp = transf_matrix_inv;
  sigma_temp = sigma;

  t = time[nstp - 1];
  prev_t = time[nstp - 2];
  sqt = sqrt(t - prev_t);

  sigma_next[0] = sig_dom[nstp - 1];
  sigma_next[1] = sig_for[nstp - 1];
  sigma_next[2] = sig_fx[nstp - 1];
  vol_fx = sigma_next[2];

  err = get_transf_matrix(sigma_next, correl[nstp - 1], transf_matrix_next,
                          transf_matrix_inv_next, covar_matrix, eigen_vector,
                          eigen_vector_transp, eigen_val);

  if (err) {
    goto FREE_RETURN;
  }

  nb_nodes_next = nb_nodes_tab[nstp - 1];
  maxmin_dual_next = maxmin_dual_tab[nstp - 1];

  t2 = clock();
  smessage("Phase 2 -initialisation        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  t1 = clock();

  /*	Final payoff valuation */

  if (!eval_evt[nstp - 1]) {
    err = "No event at last step in tree_main_3dfx";
    goto FREE_RETURN;
  }

  /*	Calculate state variables */
  coos_dual[0] = maxmin_dual_next[0][1];
  for (i = 0; i < nb_nodes_next[0]; i++) {
    coos_dual[1] = maxmin_dual_next[1][1];
    for (j = 0; j < nb_nodes_next[1]; j++) {
      coos_dual[2] = maxmin_dual_next[2][1];
      for (k = 0; k < nb_nodes_next[2]; k++) {
        prod_matrix_vector_constant(transf_matrix_next, coos_dual, sqt,
                                    sv[i][j][k], gj, gk, gsum);
        coos_dual[2] += MESH;
      } /* end for k	*/
      coos_dual[1] += MESH;
    } /* end for j	*/
    coos_dual[0] += MESH;
  } /* end for i	*/

  /*	Eval payoff */
  err = payoff_func(date[nstp - 1], time[nstp - 1], func_parm_tab[nstp - 1],
                    spot_fx, dom_yc, for_yc, nb_nodes_next[0], nb_nodes_next[1],
                    nb_nodes_next[2], sv, nprod, values_next, is_bar[nstp - 1],
                    bar_k[nstp - 1], bar_idx[nstp - 1], bar_col[nstp - 1],
                    bar_lvl[nstp - 1]);

  if (err) {
    goto FREE_RETURN;
  }

  /*	If there is a barrier        , precalculate digital parameters */
  if (is_bar[nstp - 1]) {
    do_dig = 1;
    dig_strike = bar_lvl[nstp - 1];
    dig_col = bar_col[nstp - 1];
    dig_dt = time[nstp - 1] - time[nstp - 2];
    dig_ivol = sig_fx[nstp - 1];
    dig_ivol *= sqrt(dig_dt);
  } else {
    do_dig = 0;
  }

  t2 = clock();
  smessage("Phase 3 -payoff evaluation        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  t1 = clock();

  /*	Discount */

  next_t = t;
  t = prev_t;

  for (step = nstp - 2; step >= 1; step--) {
    prev_t = time[step - 1];
    sqt_prev = sqrt(t - prev_t);
    delta_t = next_t - t;
    sqt_inv = 1.0 / sqt;

    /*	Update volatilities and related transformation matrices        , etc. */
    flag = vol_change[step];
    if (flag) {
      /*	We need to recompute everything  */
      sigma[0] = sig_dom[step];
      sigma[1] = sig_for[step];
      sigma[2] = sig_fx[step];

      /*	Calculate transformation and inverse matrices */
      err = get_transf_matrix(sigma, correl[step], transf_matrix,
                              transf_matrix_inv, covar_matrix, eigen_vector,
                              eigen_vector_transp, eigen_val);

      if (err) {
        goto FREE_RETURN;
      }

    } else {
      /*	We keep same vol and matrices */
      sigma_temp = sigma;
      transf_matrix_temp = transf_matrix;
      transf_matrix_inv_temp = transf_matrix_inv;

      sigma = sigma_next;
      transf_matrix = transf_matrix_next;
      transf_matrix_inv = transf_matrix_inv_next;
    }

    nb_nodes = nb_nodes_tab[step];
    maxmin_dual = maxmin_dual_tab[step];

    /*	Precalculate quantities used to compute
                    original state variables and dual expectations */

    n0 = nb_nodes[0];
    n1 = nb_nodes[1];
    n2 = nb_nodes[2];

    temp_vec1[0] = dr_const_dom[step] * delta_t;
    temp_vec1[1] = dr_const_for[step] * delta_t;
    temp_vec1[2] = dr_const_fx[step] * delta_t;

    prod_matrix_vector(transf_matrix_inv_next, temp_vec1, temp_vec2, i, j,
                       gsum);

    temp_mat3[0][0] = 1.0 + dr_coef_dom[step] * delta_t;
    temp_mat3[0][1] = 0.0;
    temp_mat3[0][2] = 0.0;
    temp_mat3[1][0] = dr_coef1_for[step] * delta_t;
    temp_mat3[1][1] = 1.0 + dr_coef2_for[step] * delta_t;
    temp_mat3[1][2] = dr_coef3_for[step] * delta_t;
    temp_mat3[2][0] = dr_coef1_fx[step] * delta_t;
    temp_mat3[2][1] = dr_coef2_fx[step] * delta_t;
    temp_mat3[2][2] = 1.0 + dr_coef3_fx[step] * delta_t;

    prod_matrix(temp_mat3, transf_matrix, temp_mat1, i, j, k, gsum);
    prod_matrix(transf_matrix_inv_next, temp_mat1, temp_mat2, i, j, k, gsum);

    temp1 = sqt_inv * sqt_prev;
    temp2 = sqt_prev * MESH;
    temp3 = temp1 * MESH;

    for (l = 0; l < 3; l++) {
      ca[l] = ea[l] = 0.0;
      for (j = 0; j < 3; j++) {
        ca[l] += transf_matrix[l][j] * maxmin_dual[j][1];
        ea[l] += temp_mat2[l][j] * maxmin_dual[j][1];
      }
      ca[l] *= sqt_prev;
      ea[l] *= temp1;

      ea[l] += sqt_inv * temp_vec2[l];

      cai[l] = temp2 * transf_matrix[l][0];
      caj[l] = temp2 * transf_matrix[l][1];
      cak[l] = temp2 * transf_matrix[l][2];
      cajn[l] = caj[l] * n1;
      cakn[l] = cak[l] * n2;

      eai[l] = temp3 * temp_mat2[l][0];
      eaj[l] = temp3 * temp_mat2[l][1];
      eak[l] = temp3 * temp_mat2[l][2];
      eajn[l] = eaj[l] * n1;
      eakn[l] = eak[l] * n2;
    }

    /*	Main loop on nodes */
    for (i = 0; i < n0; i++) {
      svi = sv[i];
      expect_duali = expect_dual[i];
      valuesi = values[i];
      for (j = 0; j < n1; j++) {
        svij = svi[j];
        expect_dualij = expect_duali[j];
        valuesij = valuesi[j];
        for (k = 0; k < n2; k++) {
          svijk = svij[k];
          expect_dualijk = expect_dualij[k];
          valuesijk = valuesij[k];

          /*	Calculate original state variables and dual expectations */
          for (l = 0; l < 3; l++) {
            svijk[l] = ca[l];
            expect_dualijk[l] = ea[l];
            ca[l] += cak[l];
            ea[l] += eak[l];
          }

          /*	Calculate discount rate */
          rd = 0.5 * (svijk[0] + dom_ifr[step]) * delta_t;
          rd = (1.0 - rd) / (1.0 + rd);

          /*	Recombine */
          for (gi = 0; gi < 3; gi++) {
            min_value = maxmin_dual_next[gi][1];
            F = expect_dualijk[gi];

            index_fwd = (long)((F - min_value) / MESH);
            gnu = F - min_value - index_fwd * MESH;
            if (gnu > HALF_MESH) {
              index_fwd += 1;
              gnu -= MESH;
            }
            if (index_fwd >= 1 && index_fwd <= nb_nodes_next[gi] - 2) {
              /*	Trinomial case */
              fwd_coos[gi] = index_fwd;
              proba[gi][2] = (1.0 + gnu * (gnu + MESH)) / DOUBLE_MESH_SQUARE;
              proba[gi][0] = -(gnu - MESH * proba[gi][2]) / MESH;
              proba[gi][1] = 1.0 - proba[gi][2] - proba[gi][0];
              link_type[gi] = 1;
            } else if (index_fwd <= 0) {
              /*	Binomial case down */
              fwd_coos[gi] = 1;
              link_type[gi] = 0;
              proba[gi][1] = (F - min_value) / MESH;
              proba[gi][0] = 1.0 - proba[gi][1];
            } else {
              /*	Binomial case up */
              fwd_coos[gi] = nb_nodes_next[gi] - 1;
              link_type[gi] = 0;
              proba[gi][0] = (maxmin_dual_next[gi][0] - F) / MESH;
              proba[gi][1] = 1.0 - proba[gi][0];
            }
          }

          /*	Discount */
          memset(res, 0, nprod_doubles);

          for (l = -1; l <= link_type[0]; l++) {
            proba1 = proba[0][l + 1];
            link1 = values_next[fwd_coos[0] + l];
            for (m = -1; m <= link_type[1]; m++) {
              proba12 = proba1 * proba[1][m + 1];
              link2 = link1[fwd_coos[1] + m];
              for (n = -1; n <= link_type[2]; n++) {
                link3 = link2[fwd_coos[2] + n];
                proba123 = proba12 * proba[2][n + 1];
                for (r = 0; r < nprod + do_dig; r++) {
                  res[r] += link3[r] * proba123;
                }
              }
            }
          }

          /*	Apply DF */
          for (r = 0; r < nprod + do_dig; r++) {
            valuesijk[r] = (discount ? res[r] * rd : res[r]);
          }

        } /* end for k	*/
        for (l = 0; l < 3; l++) {
          ca[l] += caj[l] - cakn[l];
          ea[l] += eaj[l] - eakn[l];
        }
      } /* end for j	*/
      for (l = 0; l < 3; l++) {
        ca[l] += cai[l] - cajn[l];
        ea[l] += eai[l] - eajn[l];
      }
    } /* end for i	*/

    /*	Add digital if relevant */
    if (do_dig) {
    }

    /*	Eval payoff */
    if (eval_evt[step]) {
      err = payoff_func(date[step], time[step], func_parm_tab[step], spot_fx,
                        dom_yc, for_yc, nb_nodes[0], nb_nodes[1], nb_nodes[2],
                        sv, nprod, values, is_bar[step], bar_k[step],
                        bar_idx[step], bar_col[step], bar_lvl[step]);

      if (err)
        goto FREE_RETURN;

      /*	If there is a barrier        , precalculate digital parameters
       */
      if (is_bar[step] == 1) {
        do_dig = 1;
        dig_strike = bar_lvl[step];
        dig_col = bar_col[step];
        dig_dt = time[step] - time[step - 1];
        dig_ivol = sig_fx[step];
        dig_ivol *= sqrt(dig_dt);
      } else {
        do_dig = 0;
      }
    }

    next_t = t;
    t = prev_t;
    sqt = sqt_prev;

    maxmin_dual_next = maxmin_dual;
    nb_nodes_next = nb_nodes;

    values_temp = values_next;
    values_next = values;
    values = values_temp;

    if (flag) {
      transf_matrix_temp = transf_matrix_next;
      transf_matrix_next = transf_matrix;
      transf_matrix = transf_matrix_temp;

      transf_matrix_inv_temp = transf_matrix_inv_next;
      transf_matrix_inv_next = transf_matrix_inv;
      transf_matrix_inv = transf_matrix_inv_temp;

      sigma_temp = sigma_next;
      sigma_next = sigma;
      sigma = sigma_temp;
      vol_fx = sigma_next[2];
    } else {
      sigma = sigma_temp;
      transf_matrix = transf_matrix_temp;
      transf_matrix_inv = transf_matrix_inv_temp;
    }
  } /* end for step */

  t2 = clock();
  smessage("Phase 4 -discounting        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  t1 = clock();

  /*	Last step: seperated because a lot of calculations are not done */

  delta_t = next_t;
  sqt = sqrt(delta_t);
  sqt_inv = 1.0 / sqt;

  sv[0][0][0][0] = sv[0][0][0][1] = sv[0][0][0][2] = 0.0;
  i = j = k = 0;

  rd = 0.5 * dom_ifr[0] * delta_t;
  rd = (1.0 - rd) / (1.0 + rd);

  expect_node[0] = dr_const_dom[0] * delta_t;
  expect_node[1] = dr_const_for[0] * delta_t;
  expect_node[2] = dr_const_fx[0] * delta_t;

  prod_matrix_vector_constant(transf_matrix_inv_next, expect_node, sqt_inv,
                              expect_dual_node, gj, gk, gsum);

  for (gi = 0; gi < 3; gi++) {
    min_value = maxmin_dual_next[gi][1];
    F = expect_dual_node[gi];

    index_fwd = (long)((F - min_value) / MESH);
    gnu = F - min_value - index_fwd * MESH;
    if (gnu > HALF_MESH) {
      index_fwd += 1;
      gnu -= MESH;
    }
    if (index_fwd >= 1 && index_fwd <= nb_nodes_next[gi] - 2) {
      fwd_coos[gi] = index_fwd;
      proba[gi][2] = (1 + gnu * (gnu + MESH)) / DOUBLE_MESH_SQUARE;
      proba[gi][0] = -(gnu - MESH * proba[gi][2]) / MESH;
      proba[gi][1] = 1.0 - proba[gi][2] - proba[gi][0];
      link_type[gi] = 1;
    } else if (index_fwd <= 0) {
      fwd_coos[gi] = 1;
      link_type[gi] = 0;
      proba[gi][1] = (F - min_value) / MESH;
      proba[gi][0] = 1.0 - proba[gi][1];
    } else {
      fwd_coos[gi] = nb_nodes_next[gi] - 1;
      link_type[gi] = 0;
      proba[gi][0] = (maxmin_dual_next[gi][0] - F) / MESH;
      proba[gi][1] = 1.0 - proba[gi][0];
    }
  }

  memset(res, 0, nprod_doubles);

  for (l = -1; l <= link_type[0]; l++) {
    proba1 = proba[0][l + 1];
    link1 = values_next[fwd_coos[0] + l];
    for (m = -1; m <= link_type[1]; m++) {
      proba12 = proba1 * proba[1][m + 1];
      link2 = link1[fwd_coos[1] + m];
      for (n = -1; n <= link_type[2]; n++) {
        proba123 = proba12 * proba[2][n + 1];
        link3 = link2[fwd_coos[2] + n];
        for (r = 0; r < nprod + do_dig; r++) {
          res[r] += link3[r] * proba123;
        }
      }
    }
  }

  if (discount) {
    for (r = 0; r < nprod + do_dig; r++) {
      res[r] *= rd;
    }
  }

  /*	Add digital if relevant */

  if (eval_evt[0]) {
    memcpy(values[0][0][0], res, nprod_doubles);

    err = payoff_func(date[0], time[0], func_parm_tab[0], spot_fx, dom_yc,
                      for_yc, 1, 1, 1, sv, nprod, values, is_bar[0], bar_k[0],
                      bar_idx[0], bar_col[0], bar_lvl[0]);

    memcpy(res, values[0][0][0], nprod_doubles);

    if (err) {
      goto FREE_RETURN;
    }
  }

  t2 = clock();
  smessage("Phase 5 -last step        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  /*	Free memory and return */

  if (bar_idx) {
    for (i = 0; i < nstp; i++) {
      if (bar_idx[i]) {
        switch (bar_k[i]) {
        case 0:
          free_lngmatrix(bar_idx[i], 0, nb_nodes_tab[i][1] - 1, 0,
                         nb_nodes_tab[i][2] - 1);
          break;
        case 1:
          free_lngmatrix(bar_idx[i], 0, nb_nodes_tab[i][0] - 1, 0,
                         nb_nodes_tab[i][2] - 1);
          break;
        case 2:
          free_lngmatrix(bar_idx[i], 0, nb_nodes_tab[i][0] - 1, 0,
                         nb_nodes_tab[i][1] - 1);
          break;
        default:
          break;
        }
      }
    }
    free(bar_idx);
  }

  if (bar_k) {
    free_lngvector(bar_k, 0, nstp - 1);
  }

  if (values) {
    free_f4tensor(values, 0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                  nb_nodes_max[2] - 1, 0, nprod - 1);
  }

  if (values_next) {
    free_f4tensor(values_next, 0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1,
                  0, nb_nodes_max[2] - 1, 0, nprod - 1);
  }

  if (sv) {
    free_f4tensor(sv, 0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                  nb_nodes_max[2] - 1, 0, 2);
  }

  if (expect_dual) {
    free_f4tensor(expect_dual, 0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1,
                  0, nb_nodes_max[2] - 1, 0, 2);
  }

  if (transf_matrix) {
    free_dmatrix(transf_matrix, 0, 2, 0, 2);
  }

  if (transf_matrix_next) {
    free_dmatrix(transf_matrix_next, 0, 2, 0, 2);
  }

  if (transf_matrix_inv) {
    free_dmatrix(transf_matrix_inv, 0, 2, 0, 2);
  }

  if (transf_matrix_inv_next) {
    free_dmatrix(transf_matrix_inv_next, 0, 2, 0, 2);
  }

  if (maxmin_dual_tab) {
    free_f3tensor(maxmin_dual_tab, 0, nstp - 1, 0, 2, 0, 1);
  }

  if (sigma) {
    free_dvector(sigma, 0, 2);
  }

  if (sigma_next) {
    free_dvector(sigma_next, 0, 2);
  }

  if (nb_nodes_tab) {
    free_lngmatrix(nb_nodes_tab, 0, nstp - 1, 0, 2);
  }

  if (nb_nodes_max) {
    free_lvector(nb_nodes_max, 0, 2);
  }

  if (correl) {
    free_dmatrix(correl, 0, nstp - 1, 0, 2);
  }

  if (glob_correl) {
    free_dmatrix(glob_correl, 0, nstp - 1, 0, 2);
  }

  if (covar_matrix) {
    free_dmatrix(covar_matrix, 0, 2, 0, 2);
  }

  if (eigen_vector) {
    free_dmatrix(eigen_vector, 0, 2, 0, 2);
  }

  if (eigen_vector_transp) {
    free_dmatrix(eigen_vector_transp, 0, 2, 0, 2);
  }

  if (eigen_val) {
    free_dvector(eigen_val, 0, 2);
  }

  if (maxmin_orig) {
    free_dmatrix(maxmin_orig, 0, 2, 0, 1);
  }

  if (expect_node) {
    free_dvector(expect_node, 0, 2);
  }

  if (expect_dual_node) {
    free_dvector(expect_dual_node, 0, 2);
  }

  if (proba) {
    free_dmatrix(proba, 0, 2, 0, 2);
  }

  if (coos_dual) {
    free_dvector(coos_dual, 0, 2);
  }

  if (temp_vec1) {
    free_dvector(temp_vec1, 0, 2);
  }

  if (temp_vec2) {
    free_dvector(temp_vec2, 0, 2);
  }

  if (temp_mat1) {
    free_dmatrix(temp_mat1, 0, 2, 0, 2);
  }

  if (temp_mat2) {
    free_dmatrix(temp_mat2, 0, 2, 0, 2);
  }

  if (temp_mat3) {
    free_dmatrix(temp_mat3, 0, 2, 0, 2);
  }

  return err;
}

/*	Main function under QTStar */
/*	-------------------------- */
Err tree_main_3dBetaDLM_QTStar(
    /*	Time data */
    long nstp, double *time, double *date,
    /*	Model data */
    int *vol_change, /*	1 if one of the volatlities has changed        ,
                                       0 otherwise */
    double *sig_dom, /*	Term structures */
    double *mu_quanto_const, double *mu_quanto_lin, double *sig_for,
    double *sig_fx, double *dom_fwd, double *dom_var, double *for_fwd,
    double *for_var, double *fx_fwd, double *fx_var, double dom_lam,
    double for_lam, double *corr_dom_for, double *corr_dom_fx,
    double *corr_for_fx, double *glob_corr_dom_for, double *glob_corr_dom_fx,
    double *glob_corr_for_fx,

    /*	Product data */
    void **func_parm_tab, int *eval_evt, double *bar_lvl, int *bar_col,
    int *is_bar,

    /*	Payoff function */
    Err (*payoff_func)(/* Event */
                       double evt_date, double evt_time, void *func_parm,
                       /* Nodes data */
                       long n1, long n2, long n3,
                       /* i: d1        , j: d2        , k: d3        ,
                                       l = {0: xDom        , 1: xFor        , 2:
                          log (Fx/Fx0)} */
                       double ****sv,
                       /* Vector of results to be updated */
                       long nprod, double ****prod_val,
                       /* Barrier details */
                       int is_bar, int bar_k, int **bar_idx, int bar_col,
                       double bar_lvl),
    /*	Result */
    int nprod, double *res) {

  long i, j, k, step, l, m, n;
  int r;
  double t, prev_t, next_t, sqt_prev, delta_t, sqt, sqt_inv;
  long gi, gj, gk;
  long index_fwd;
  double gsum, vol_fx;
  double F, min_value, gnu;
  double proba1, proba12, proba123;
  double a, b;
  int flag;

  double ***link1 = NULL, **link2 = NULL, *link3 = NULL;

  double ****values = NULL, ****values_next = NULL, ****values_temp = NULL,
         ***valuesi = NULL, **valuesij = NULL, *valuesijk = NULL, ****sv = NULL,
         ***svi = NULL, **svij = NULL, *svijk = NULL, ****expect_dual = NULL,
         ***expect_duali = NULL, **expect_dualij = NULL, *expect_dualijk = NULL;

  double **transf_matrix = NULL, **transf_matrix_next = NULL,
         **transf_matrix_temp = NULL, **transf_matrix_inv = NULL,
         **transf_matrix_inv_next = NULL, **transf_matrix_inv_temp = NULL,
         **maxmin_dual = NULL, ***maxmin_dual_tab = NULL,
         **maxmin_dual_next = NULL, *sigma = NULL, *sigma_next = NULL,
         *volt = NULL, *sigma_temp = NULL;

  long **nb_nodes_tab = NULL, *nb_nodes = NULL, *nb_nodes_max = NULL,
       *nb_nodes_next = NULL;

  long *bar_k = NULL, ***bar_idx = NULL;

  double **correl = NULL, **glob_correl = NULL, **covar_matrix = NULL,
         **eigen_vector = NULL, **eigen_vector_transp = NULL, *eigen_val = NULL,
         **maxmin_orig = NULL, *expect_node = NULL, *expect_dual_node = NULL,
         **proba = NULL, *coos_dual = NULL;

  double *dom_ifr = NULL, *for_ifr = NULL;

  double temp1, temp2, temp3;

  double ca[3], cai[3], caj[3], cak[3], cajn[3], cakn[3];
  double ea[3], eai[3], eaj[3], eak[3], eajn[3], eakn[3];

  time_t t1, t2;

  size_t nprod_doubles = (nprod + 1) * sizeof(double);

  long fwd_coos[3], link_type[3];

  long n0, n1, n2;

  Err err = NULL;

  /*	Allocate memory once for all */

  t1 = clock();

  temp_vec1 = temp_vec2 = NULL;
  temp_mat1 = temp_mat2 = temp_mat3 = NULL;

  transf_matrix = dmatrix(0, 2, 0, 2);
  transf_matrix_next = dmatrix(0, 2, 0, 2);
  transf_matrix_inv = dmatrix(0, 2, 0, 2);
  transf_matrix_inv_next = dmatrix(0, 2, 0, 2);
  sigma = dvector(0, 2);
  sigma_next = dvector(0, 2);
  maxmin_dual_tab = f3tensor(0, nstp - 1, 0, 2, 0, 1);
  nb_nodes_tab = lngmatrix(0, nstp - 1, 0, 2);
  nb_nodes_max = lvector(0, 2);
  correl = dmatrix(0, nstp - 1, 0, 2);
  glob_correl = dmatrix(0, nstp - 1, 0, 2);
  covar_matrix = dmatrix(0, 2, 0, 2);
  eigen_vector = dmatrix(0, 2, 0, 2);
  eigen_vector_transp = dmatrix(0, 2, 0, 2);
  eigen_val = dvector(0, 2);
  maxmin_orig = dmatrix(0, 2, 0, 1);
  expect_node = dvector(0, 2);
  expect_dual_node = dvector(0, 2);
  proba = dmatrix(0, 2, 0, 2);
  coos_dual = dvector(0, 2);

  temp_vec1 = dvector(0, 2);
  temp_vec2 = dvector(0, 2);
  temp_mat1 = dmatrix(0, 2, 0, 2);
  temp_mat2 = dmatrix(0, 2, 0, 2);
  temp_mat3 = dmatrix(0, 2, 0, 2);

  if (!transf_matrix || !transf_matrix_next || !transf_matrix_inv ||
      !transf_matrix_inv_next || !maxmin_dual_tab || !sigma || !sigma_next ||
      !nb_nodes_tab || !nb_nodes_max || !correl || !glob_correl ||
      !covar_matrix || !eigen_vector || !eigen_vector_transp || !eigen_val ||
      !maxmin_orig || !expect_node || !expect_dual_node || !proba ||
      !coos_dual || !temp_vec1 || !temp_vec2 || !temp_mat1 || !temp_mat2 ||
      !temp_mat3) {
    err = "Memory allocation error (1) in tree_main_3dfx";
    goto FREE_RETURN;
  }

  /*	Calculates all boundaries and number of nodes */

  for (i = 0; i < nstp; i++) {
    correl[i][0] = corr_dom_for[i];
    correl[i][1] = corr_dom_fx[i];
    correl[i][2] = corr_for_fx[i];

    glob_correl[i][0] = glob_corr_dom_for[i];
    glob_correl[i][1] = glob_corr_dom_fx[i];
    glob_correl[i][2] = glob_corr_for_fx[i];
  }

  bar_k = lngvector(0, nstp - 1);
  bar_idx = (long ***)calloc(nstp, sizeof(long **));

  if (!bar_k || !bar_idx) {
    err = "Memory allocation (2) error in tree_main_3dfx";
    goto FREE_RETURN;
  }

  err = calculate_all_boundaries(
      nstp, time, vol_change, sig_dom, sig_for, sig_fx, dom_fwd, dom_var,
      for_fwd, for_var, fx_fwd, fx_var, correl, glob_correl, is_bar, bar_lvl,
      bar_k, bar_idx, nprod, nb_nodes_tab, nb_nodes_max, maxmin_dual_tab, sigma,
      maxmin_orig, transf_matrix, transf_matrix_inv, covar_matrix, eigen_vector,
      eigen_vector_transp, eigen_val);

  if (err) {
    goto FREE_RETURN;
  }

  /*	Now we can allocate memory for the tensors */

  values = f4tensor(0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                    nb_nodes_max[2] - 1, 0, nprod);
  values_next = f4tensor(0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                         nb_nodes_max[2] - 1, 0, nprod);
  sv = f4tensor(0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                nb_nodes_max[2] - 1, 0, 2);
  expect_dual = f4tensor(0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                         nb_nodes_max[2] - 1, 0, 2);

  if (!values || !values_next || !sv || !expect_dual) {
    err = "Memory allocation (3) error in tree_main_3dfx";
    goto FREE_RETURN;
  }

  /*	Initialisation */

  transf_matrix_temp = transf_matrix;
  transf_matrix_inv_temp = transf_matrix_inv;
  sigma_temp = sigma;

  t = time[nstp - 1];
  prev_t = time[nstp - 2];
  sqt = sqrt(t - prev_t);

  sigma_next[0] = sig_dom[nstp - 1];
  sigma_next[1] = sig_for[nstp - 1];
  sigma_next[2] = sig_fx[nstp - 1];
  vol_fx = sigma_next[2];

  err = get_transf_matrix(sigma_next, correl[nstp - 1], transf_matrix_next,
                          transf_matrix_inv_next, covar_matrix, eigen_vector,
                          eigen_vector_transp, eigen_val);

  if (err) {
    goto FREE_RETURN;
  }

  nb_nodes_next = nb_nodes_tab[nstp - 1];
  maxmin_dual_next = maxmin_dual_tab[nstp - 1];

  t2 = clock();
  smessage("Phase 2 -initialisation        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  t1 = clock();

  /*	Final payoff valuation */

  if (!eval_evt[nstp - 1]) {
    err = "No event at last step in tree_main_3dfx";
    goto FREE_RETURN;
  }

  /*	Calculate state variables */
  coos_dual[0] = maxmin_dual_next[0][1];
  for (i = 0; i < nb_nodes_next[0]; i++) {
    coos_dual[1] = maxmin_dual_next[1][1];
    for (j = 0; j < nb_nodes_next[1]; j++) {
      coos_dual[2] = maxmin_dual_next[2][1];
      for (k = 0; k < nb_nodes_next[2]; k++) {
        prod_matrix_vector_constant(transf_matrix_next, coos_dual, sqt,
                                    sv[i][j][k], gj, gk, gsum);
        coos_dual[2] += MESH;
      } /* end for k	*/
      coos_dual[1] += MESH;
    } /* end for j	*/
    coos_dual[0] += MESH;
  } /* end for i	*/

  /*	Eval payoff */
  err = payoff_func(date[nstp - 1], time[nstp - 1], func_parm_tab[nstp - 1],
                    nb_nodes_next[0], nb_nodes_next[1], nb_nodes_next[2], sv,
                    nprod, values_next, is_bar[nstp - 1], bar_k[nstp - 1],
                    bar_idx[nstp - 1], bar_col[nstp - 1], bar_lvl[nstp - 1]);

  if (err)
    goto FREE_RETURN;

  t2 = clock();
  smessage("Phase 3 -payoff evaluation        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  t1 = clock();

  /*	Discount */

  next_t = t;
  t = prev_t;

  for (step = nstp - 2; step >= 1; step--) {
    prev_t = time[step - 1];
    sqt_prev = sqrt(t - prev_t);
    delta_t = next_t - t;
    sqt_inv = 1.0 / sqt;

    /*	Update volatilities and related transformation matrices        , etc. */
    flag = vol_change[step];
    if (flag) {
      /*	We need to recompute everything  */
      sigma[0] = sig_dom[step];
      sigma[1] = sig_for[step];
      sigma[2] = sig_fx[step];

      /*	Calculate transformation and inverse matrices */
      err = get_transf_matrix(sigma, correl[step], transf_matrix,
                              transf_matrix_inv, covar_matrix, eigen_vector,
                              eigen_vector_transp, eigen_val);

      if (err)
        goto FREE_RETURN;
    } else {
      /*	We keep same vol and matrices */
      sigma_temp = sigma;
      transf_matrix_temp = transf_matrix;
      transf_matrix_inv_temp = transf_matrix_inv;

      sigma = sigma_next;
      transf_matrix = transf_matrix_next;
      transf_matrix_inv = transf_matrix_inv_next;
    }

    nb_nodes = nb_nodes_tab[step];
    maxmin_dual = maxmin_dual_tab[step];

    /*	Precalculate quantities used to compute
                    original state variables and dual expectations */

    a = mu_quanto_const[step] * delta_t;
    b = mu_quanto_lin[step] * delta_t;

    n0 = nb_nodes[0];
    n1 = nb_nodes[1];
    n2 = nb_nodes[2];

    prod_matrix(transf_matrix_inv_next, transf_matrix, temp_mat2, i, j, k,
                gsum);

    for (l = 0; l < 3; l++) {
      temp_vec2[l] = transf_matrix_inv_next[l][1] * a;

      for (j = 0; j < 3; j++) {
        temp_mat2[l][j] +=
            transf_matrix_inv_next[l][1] * b * transf_matrix[2][j];
      }
    }

    temp1 = sqt_inv * sqt_prev;
    temp2 = sqt_prev * MESH;
    temp3 = temp1 * MESH;

    if (eval_evt[step]) {
      for (l = 0; l < 3; l++) {
        ca[l] = ea[l] = 0.0;
        for (j = 0; j < 3; j++) {
          ca[l] += transf_matrix[l][j] * maxmin_dual[j][1];
          ea[l] += temp_mat2[l][j] * maxmin_dual[j][1];
        }
        ca[l] *= sqt_prev;
        ea[l] *= temp1;

        ea[l] += sqt_inv * temp_vec2[l];

        cai[l] = temp2 * transf_matrix[l][0];
        caj[l] = temp2 * transf_matrix[l][1];
        cak[l] = temp2 * transf_matrix[l][2];
        cajn[l] = caj[l] * n1;
        cakn[l] = cak[l] * n2;

        eai[l] = temp3 * temp_mat2[l][0];
        eaj[l] = temp3 * temp_mat2[l][1];
        eak[l] = temp3 * temp_mat2[l][2];
        eajn[l] = eaj[l] * n1;
        eakn[l] = eak[l] * n2;
      }

      /*	Main loop on nodes */
      for (i = 0; i < n0; i++) {
        svi = sv[i];
        expect_duali = expect_dual[i];
        valuesi = values[i];
        for (j = 0; j < n1; j++) {
          svij = svi[j];
          expect_dualij = expect_duali[j];
          valuesij = valuesi[j];
          for (k = 0; k < n2; k++) {
            svijk = svij[k];
            expect_dualijk = expect_dualij[k];
            valuesijk = valuesij[k];

            /*	Calculate original state variables and dual expectations */
            for (l = 0; l < 3; l++) {
              svijk[l] = ca[l];
              expect_dualijk[l] = ea[l];
              ca[l] += cak[l];
              ea[l] += eak[l];
            }

            /*	Recombine */
            for (gi = 0; gi < 3; gi++) {
              min_value = maxmin_dual_next[gi][1];
              F = expect_dualijk[gi];

              index_fwd = (long)((F - min_value) / MESH);
              gnu = F - min_value - index_fwd * MESH;
              if (gnu > HALF_MESH) {
                index_fwd += 1;
                gnu -= MESH;
              }
              if (index_fwd >= 1 && index_fwd <= nb_nodes_next[gi] - 2) {
                /*	Trinomial case */
                fwd_coos[gi] = index_fwd;
                proba[gi][2] = (1.0 + gnu * (gnu + MESH)) / DOUBLE_MESH_SQUARE;
                proba[gi][0] = -(gnu - MESH * proba[gi][2]) / MESH;
                proba[gi][1] = 1.0 - proba[gi][2] - proba[gi][0];
                link_type[gi] = 1;
              } else if (index_fwd <= 0) {
                /*	Binomial case down */
                fwd_coos[gi] = 1;
                link_type[gi] = 0;
                proba[gi][1] = (F - min_value) / MESH;
                proba[gi][0] = 1.0 - proba[gi][1];
              } else {
                /*	Binomial case up */
                fwd_coos[gi] = nb_nodes_next[gi] - 1;
                link_type[gi] = 0;
                proba[gi][0] = (maxmin_dual_next[gi][0] - F) / MESH;
                proba[gi][1] = 1.0 - proba[gi][0];
              }
            }

            /*	Discount */
            memset(valuesijk, 0, nprod_doubles);

            for (l = -1; l <= link_type[0]; l++) {
              proba1 = proba[0][l + 1];
              link1 = values_next[fwd_coos[0] + l];
              for (m = -1; m <= link_type[1]; m++) {
                proba12 = proba1 * proba[1][m + 1];
                link2 = link1[fwd_coos[1] + m];
                for (n = -1; n <= link_type[2]; n++) {
                  link3 = link2[fwd_coos[2] + n];
                  proba123 = proba12 * proba[2][n + 1];

                  for (r = 0; r < nprod; r++) {
                    valuesijk[r] += link3[r] * proba123;
                  }
                }
              }
            }
          } /* end for k	*/

          for (l = 0; l < 3; l++) {
            ca[l] += caj[l] - cakn[l];
            ea[l] += eaj[l] - eakn[l];
          }
        } /* end for j	*/
        for (l = 0; l < 3; l++) {
          ca[l] += cai[l] - cajn[l];
          ea[l] += eai[l] - eajn[l];
        }
      } /* end for i */

      /*	Eval payoff */
      err =
          payoff_func(date[step], time[step], func_parm_tab[step], nb_nodes[0],
                      nb_nodes[1], nb_nodes[2], sv, nprod, values, is_bar[step],
                      bar_k[step], bar_idx[step], bar_col[step], bar_lvl[step]);

      if (err)
        goto FREE_RETURN;
    } else {
      for (l = 0; l < 3; l++) {
        ea[l] = 0.0;
        for (j = 0; j < 3; j++) {
          ea[l] += temp_mat2[l][j] * maxmin_dual[j][1];
        }

        ea[l] *= temp1;
        ea[l] += sqt_inv * temp_vec2[l];

        eai[l] = temp3 * temp_mat2[l][0];
        eaj[l] = temp3 * temp_mat2[l][1];
        eak[l] = temp3 * temp_mat2[l][2];
        eajn[l] = eaj[l] * n1;
        eakn[l] = eak[l] * n2;
      }

      /*	Main loop on nodes */
      for (i = 0; i < n0; i++) {
        expect_duali = expect_dual[i];
        valuesi = values[i];
        for (j = 0; j < n1; j++) {
          expect_dualij = expect_duali[j];
          valuesij = valuesi[j];
          for (k = 0; k < n2; k++) {
            expect_dualijk = expect_dualij[k];
            valuesijk = valuesij[k];

            /*	Calculate original state variables and dual expectations */
            for (l = 0; l < 3; l++) {
              expect_dualijk[l] = ea[l];
              ea[l] += eak[l];
            }

            /*	Recombine */
            for (gi = 0; gi < 3; gi++) {
              min_value = maxmin_dual_next[gi][1];
              F = expect_dualijk[gi];

              index_fwd = (long)((F - min_value) / MESH);
              gnu = F - min_value - index_fwd * MESH;
              if (gnu > HALF_MESH) {
                index_fwd += 1;
                gnu -= MESH;
              }
              if (index_fwd >= 1 && index_fwd <= nb_nodes_next[gi] - 2) {
                /*	Trinomial case */
                fwd_coos[gi] = index_fwd;
                proba[gi][2] = (1.0 + gnu * (gnu + MESH)) / DOUBLE_MESH_SQUARE;
                proba[gi][0] = -(gnu - MESH * proba[gi][2]) / MESH;
                proba[gi][1] = 1.0 - proba[gi][2] - proba[gi][0];
                link_type[gi] = 1;
              } else if (index_fwd <= 0) {
                /*	Binomial case down */
                fwd_coos[gi] = 1;
                link_type[gi] = 0;
                proba[gi][1] = (F - min_value) / MESH;
                proba[gi][0] = 1.0 - proba[gi][1];
              } else {
                /*	Binomial case up */
                fwd_coos[gi] = nb_nodes_next[gi] - 1;
                link_type[gi] = 0;
                proba[gi][0] = (maxmin_dual_next[gi][0] - F) / MESH;
                proba[gi][1] = 1.0 - proba[gi][0];
              }
            }

            /*	Discount */
            memset(valuesijk, 0, nprod_doubles);

            for (l = -1; l <= link_type[0]; l++) {
              proba1 = proba[0][l + 1];
              link1 = values_next[fwd_coos[0] + l];
              for (m = -1; m <= link_type[1]; m++) {
                proba12 = proba1 * proba[1][m + 1];
                link2 = link1[fwd_coos[1] + m];
                for (n = -1; n <= link_type[2]; n++) {
                  link3 = link2[fwd_coos[2] + n];
                  proba123 = proba12 * proba[2][n + 1];
                  for (r = 0; r < nprod; r++) {
                    valuesijk[r] += link3[r] * proba123;
                  }
                }
              }
            }

          } /* end for k	*/
          for (l = 0; l < 3; l++) {
            ea[l] += eaj[l] - eakn[l];
          }
        } /* end for j	*/
        for (l = 0; l < 3; l++) {
          ea[l] += eai[l] - eajn[l];
        }
      } /* end for i	*/
    }

    next_t = t;
    t = prev_t;
    sqt = sqt_prev;

    maxmin_dual_next = maxmin_dual;
    nb_nodes_next = nb_nodes;

    values_temp = values_next;
    values_next = values;
    values = values_temp;

    if (flag) {
      transf_matrix_temp = transf_matrix_next;
      transf_matrix_next = transf_matrix;
      transf_matrix = transf_matrix_temp;

      transf_matrix_inv_temp = transf_matrix_inv_next;
      transf_matrix_inv_next = transf_matrix_inv;
      transf_matrix_inv = transf_matrix_inv_temp;

      sigma_temp = sigma_next;
      sigma_next = sigma;
      sigma = sigma_temp;
      vol_fx = sigma_next[2];
    } else {
      sigma = sigma_temp;
      transf_matrix = transf_matrix_temp;
      transf_matrix_inv = transf_matrix_inv_temp;
    }
  } /* end for step */

  t2 = clock();
  smessage("Phase 4 -discounting        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);
  t1 = clock();

  /*	Last step: seperated because a lot of calculations are not done */

  delta_t = next_t;
  sqt = sqrt(delta_t);
  sqt_inv = 1.0 / sqt;

  sv[0][0][0][0] = sv[0][0][0][1] = sv[0][0][0][2] = 0.0;
  i = j = k = 0;

  a = mu_quanto_const[0] * delta_t;
  b = mu_quanto_lin[0] * delta_t;

  for (l = 0; l < 3; l++) {
    expect_dual_node[l] = a * transf_matrix_inv_next[l][1] * sqt_inv;
  }

  for (gi = 0; gi < 3; gi++) {
    min_value = maxmin_dual_next[gi][1];
    F = expect_dual_node[gi];

    index_fwd = (long)((F - min_value) / MESH);
    gnu = F - min_value - index_fwd * MESH;
    if (gnu > HALF_MESH) {
      index_fwd += 1;
      gnu -= MESH;
    }
    if (index_fwd >= 1 && index_fwd <= nb_nodes_next[gi] - 2) {
      fwd_coos[gi] = index_fwd;
      proba[gi][2] = (1 + gnu * (gnu + MESH)) / DOUBLE_MESH_SQUARE;
      proba[gi][0] = -(gnu - MESH * proba[gi][2]) / MESH;
      proba[gi][1] = 1.0 - proba[gi][2] - proba[gi][0];
      link_type[gi] = 1;
    } else if (index_fwd <= 0) {
      fwd_coos[gi] = 1;
      link_type[gi] = 0;
      proba[gi][1] = (F - min_value) / MESH;
      proba[gi][0] = 1.0 - proba[gi][1];
    } else {
      fwd_coos[gi] = nb_nodes_next[gi] - 1;
      link_type[gi] = 0;
      proba[gi][0] = (maxmin_dual_next[gi][0] - F) / MESH;
      proba[gi][1] = 1.0 - proba[gi][0];
    }
  }

  memset(res, 0, nprod_doubles);

  for (l = -1; l <= link_type[0]; l++) {
    proba1 = proba[0][l + 1];
    link1 = values_next[fwd_coos[0] + l];
    for (m = -1; m <= link_type[1]; m++) {
      proba12 = proba1 * proba[1][m + 1];
      link2 = link1[fwd_coos[1] + m];
      for (n = -1; n <= link_type[2]; n++) {
        proba123 = proba12 * proba[2][n + 1];
        link3 = link2[fwd_coos[2] + n];
        for (r = 0; r < nprod; r++) {
          res[r] += link3[r] * proba123;
        }
      }
    }
  }

  if (eval_evt[0]) {
    memcpy(values[0][0][0], res, nprod_doubles);

    err = payoff_func(date[0], time[0], func_parm_tab[0], 1, 1, 1, sv, nprod,
                      values, is_bar[0], bar_k[0], bar_idx[0], bar_col[0],
                      bar_lvl[0]);

    memcpy(res, values[0][0][0], nprod_doubles);

    if (err) {
      goto FREE_RETURN;
    }
  }

  t2 = clock();
  smessage("Phase 5 -last step        , time in sec: %.2f",
           (double)(t2 - t1) / CLOCKS_PER_SEC);

FREE_RETURN:

  /*	Free memory and return */

  if (bar_idx) {
    for (i = 0; i < nstp; i++) {
      if (bar_idx[i]) {
        switch (bar_k[i]) {
        case 0:
          free_lngmatrix(bar_idx[i], 0, nb_nodes_tab[i][1] - 1, 0,
                         nb_nodes_tab[i][2] - 1);
          break;
        case 1:
          free_lngmatrix(bar_idx[i], 0, nb_nodes_tab[i][0] - 1, 0,
                         nb_nodes_tab[i][2] - 1);
          break;
        case 2:
          free_lngmatrix(bar_idx[i], 0, nb_nodes_tab[i][0] - 1, 0,
                         nb_nodes_tab[i][1] - 1);
          break;
        default:
          break;
        }
      }
    }
    free(bar_idx);
  }

  if (bar_k) {
    free_lngvector(bar_k, 0, nstp - 1);
  }

  if (values) {
    free_f4tensor(values, 0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                  nb_nodes_max[2] - 1, 0, nprod - 1);
  }

  if (values_next) {
    free_f4tensor(values_next, 0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1,
                  0, nb_nodes_max[2] - 1, 0, nprod - 1);
  }

  if (sv) {
    free_f4tensor(sv, 0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1, 0,
                  nb_nodes_max[2] - 1, 0, 2);
  }

  if (expect_dual) {
    free_f4tensor(expect_dual, 0, nb_nodes_max[0] - 1, 0, nb_nodes_max[1] - 1,
                  0, nb_nodes_max[2] - 1, 0, 2);
  }

  if (transf_matrix) {
    free_dmatrix(transf_matrix, 0, 2, 0, 2);
  }

  if (transf_matrix_next) {
    free_dmatrix(transf_matrix_next, 0, 2, 0, 2);
  }

  if (transf_matrix_inv) {
    free_dmatrix(transf_matrix_inv, 0, 2, 0, 2);
  }

  if (transf_matrix_inv_next) {
    free_dmatrix(transf_matrix_inv_next, 0, 2, 0, 2);
  }

  if (maxmin_dual_tab) {
    free_f3tensor(maxmin_dual_tab, 0, nstp - 1, 0, 2, 0, 1);
  }

  if (sigma) {
    free_dvector(sigma, 0, 2);
  }

  if (sigma_next) {
    free_dvector(sigma_next, 0, 2);
  }

  if (nb_nodes_tab) {
    free_lngmatrix(nb_nodes_tab, 0, nstp - 1, 0, 2);
  }

  if (nb_nodes_max) {
    free_lvector(nb_nodes_max, 0, 2);
  }

  if (correl) {
    free_dmatrix(correl, 0, nstp - 1, 0, 2);
  }

  if (glob_correl) {
    free_dmatrix(glob_correl, 0, nstp - 1, 0, 2);
  }

  if (covar_matrix) {
    free_dmatrix(covar_matrix, 0, 2, 0, 2);
  }

  if (eigen_vector) {
    free_dmatrix(eigen_vector, 0, 2, 0, 2);
  }

  if (eigen_vector_transp) {
    free_dmatrix(eigen_vector_transp, 0, 2, 0, 2);
  }

  if (eigen_val) {
    free_dvector(eigen_val, 0, 2);
  }

  if (maxmin_orig) {
    free_dmatrix(maxmin_orig, 0, 2, 0, 1);
  }

  if (expect_node) {
    free_dvector(expect_node, 0, 2);
  }

  if (expect_dual_node) {
    free_dvector(expect_dual_node, 0, 2);
  }

  if (proba) {
    free_dmatrix(proba, 0, 2, 0, 2);
  }

  if (coos_dual) {
    free_dvector(coos_dual, 0, 2);
  }

  if (temp_vec1) {
    free_dvector(temp_vec1, 0, 2);
  }

  if (temp_vec2) {
    free_dvector(temp_vec2, 0, 2);
  }

  if (temp_mat1) {
    free_dmatrix(temp_mat1, 0, 2, 0, 2);
  }

  if (temp_mat2) {
    free_dmatrix(temp_mat2, 0, 2, 0, 2);
  }

  if (temp_mat3) {
    free_dmatrix(temp_mat3, 0, 2, 0, 2);
  }

  return err;
}