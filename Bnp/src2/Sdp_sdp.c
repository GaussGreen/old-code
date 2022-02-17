/*  Author: AdA */
/*  Date:   22/06/1999 */

/*  Associated header: sdp.h */
/*  Uses: some nr routines and some BLAS routines.  */

/*  BLAS from the Intel Math Kernel library can be linked directly  ,  */
/*  or the MiniBLAS (ToDo) provided for use in griffin can be used instead. */
/*  See sdp.h for BLAS library choice. */

/*  This defines the main routine for Semi Deinite Programming */
/*  It takes a set of matrix and boubles as constraints and  */
/*  a matrix for the optimisation direction. */

#include "sdp_sdp.h"
#include "math.h"
#include "sdp_sdplib.h"
#include "time.h"

Err srt_sdp(double ***const_mat, double *const_val, double **cmat, int dim,
            int num_const, double **xmat, double *obj_val, double *error,
            int *return_error, double toler, int printlevel,
            int max_iterations) {
  /*  General parameters for the algorithm. */
  int max_iter = max_iterations;
  double tolerance = toler;
  double start_gamma = 0.8;
  double frac_gamma_reduction = 0.65;
  double expon = 3;
  double first_kappa = 1.0;
  double first_tau = 1.0;
  double line_reduc = 0.6;
  int default_num_improve = 0;
  double condition_bound = 1e+35;
  double tolerance_infeasibility = 10e+03 * tolerance;

  /*  ---------- */
  /*  Algortihm variables. */
  Err err = NULL;
  int go_on = 0;
  int iter = 0;
  double **zmat = dmatrix(1, dim, 1, dim);
  double **dpzmat = dmatrix(1, dim, 1, dim);
  double **dpxmat = dmatrix(1, dim, 1, dim);
  double **dczmat = dmatrix(1, dim, 1, dim);
  double **dcxmat = dmatrix(1, dim, 1, dim);
  double **gram = dmatrix(1, num_const, 1, num_const);
  double **inverse_gram = dmatrix(1, num_const, 1, num_const);
  double **big_inv_E = dmatrix(1, dim * dim, 1, dim * dim);
  double *yvec = dvector(1, num_const);
  double *dpyvec = dvector(1, num_const);
  double *dcyvec = dvector(1, num_const);
  double gamma = start_gamma;
  double tau = first_tau, dptau, dctau;
  double primal_tau, dual_tau;
  double kappa = first_kappa, dpkappa, dckappa;
  double alphap, betap, alphac = 1.0, betac = 1.0;
  double mu, first_mu;
  double phi;
  double sigma;
  double eta;
  /* ----------- */
  /*  Working variables. */
  double **status_tracker = dmatrix(0, max_iter, 1, 11);
  double **matbuf1 = dmatrix(1, num_const + 1, 1, num_const + 1);
  double **systmat = dmatrix(1, num_const + 1, 1, num_const + 1);
  double *vecbuf1 = dvector(0, num_const + 1);
  double *vecbuf2 = dvector(0, num_const + 1);
  double *bufvec1 = dvector(0, dim);
  double numer_error, bufpinf, bufdinf;
  time_t start_time, finish_time;
  int feasibility_step = 0;
  static int num_calls = 0;

  /*  Initializing. */
  time(&start_time);
  if ((num_calls == 0) && (printlevel >= 1)) {
    smessage("SDP  , Test version. \n");
    smessage("Go ahead  , make my day... \n");
    smessage("\n");
  }
  /*  Starting point. Initialization. */
  num_calls++;
  if (printlevel >= 2)
    smessage("SDP: Initializing... Call number: %d \n", num_calls);
  init_sol(const_mat, num_const, const_val, cmat, dim, xmat, zmat);
  /*  Compute the inverse Gram matrix for the A operator  */
  /*  for use in the projection on the fesible set. */
  compute_gram(const_mat, num_const, dim, gram);
  err = inv_mat(gram, num_const, inverse_gram);
  mat_mult_nobuf(gram, inverse_gram, num_const, matbuf1);
  /*  Test for constraints independence. */
  if (err != NULL)
    go_on = 8;
  /* 	   ************************ Start of main loop. ******************* */
  while (go_on == 0) {
    iter++;
    mu = compute_mu(xmat, zmat, dim, tau, kappa);
    phi = compute_infeas(const_mat, const_val, num_const, xmat, zmat, cmat, dim,
                         yvec, tau, &bufpinf, &bufdinf);
    /*  Status display and tracking. In debug send the status_tracker to the
     * monitor. */
    status(iter, const_mat, const_val, num_const, systmat, xmat, zmat, cmat,
           dim, yvec, tau, kappa, sigma, eta, alphac, betac, status_tracker,
           printlevel, max_iter, numer_error);
    /*  Send matrix to monitor if in debug. **** DEBUG ****** */
    visu_mat(status_tracker, max_iter, 11, 1);
    /*  ---------------------------- */
    if (iter == 1)
      first_mu = mu;
    /*  Test for success. */
    if (global_error(const_mat, const_val, num_const, xmat, zmat, cmat, dim,
                     yvec, tau, kappa) < tolerance)
      go_on = 1;
    /*  Or for infeasibility... */
    if ((fabs(mu / first_mu) < tolerance_infeasibility) &&
        ((tau / kappa) / (first_tau / first_kappa) < tolerance_infeasibility))
      go_on = 2;
    /*  Then test for increasing primal infeasibility. */
    if ((status_tracker[iter - 1][2] <= status_tracker[iter][2]) &&
        (feasibility_step == 0) && (iter >= 10))
      feasibility_step = 0;
    else
      feasibility_step = 0;
    /*  *************************** FEASIBILITY STEP ******************* */
    if (feasibility_step == 1) {
      take_feasibility_step(const_mat, const_val, num_const, xmat, dcxmat, zmat,
                            dczmat, dim, yvec, dcyvec, &dctau, tau, &dckappa,
                            inverse_gram);
      alphac =
          step_length(xmat, dcxmat, dim, tau, dctau, kappa, dckappa, gamma);
      betac = step_length(zmat, dczmat, dim, tau, dctau, kappa, dckappa, gamma);
      set_steps(const_mat, const_val, num_const, yvec, dcyvec, xmat, dcxmat,
                zmat, dczmat, cmat, dim, kappa, dckappa, tau, dctau, eta,
                &alphac, &betac, line_reduc, "Corrector");
    }
    /*  *************************** PREDICTOR STEP ********************* */
    if ((go_on == 0) && (feasibility_step == 0)) {
      /*  Predictor step: sigma=0  , eta=1. */
      sigma = 0;
      eta = 1;
      /*  Solve for (dy  ,dtau). */
      err =
          system_matrix(const_mat, num_const, xmat, zmat, cmat, dim, "AHO",
                        systmat, big_inv_E, &numer_error, default_num_improve);
      /*  Test for singular system operator; */
      if (err != NULL)
        go_on = 7;
      left_matrix(systmat, num_const, const_val, kappa, tau, systmat);
      right_vector(const_mat, num_const, big_inv_E, xmat, zmat, cmat, dpxmat,
                   dpzmat, "Predictor", dim, yvec, const_val, eta, tau, sigma,
                   mu, kappa, dptau, dpkappa, vecbuf1, default_num_improve);
      /*  Solve the main system for (dy.dtau): */
      /*  Test the condition number  , if to bad  , stops. */
      if (condition_mat(systmat, num_const + 1) > condition_bound)
        go_on = 6;
      err = lsolve(systmat, vecbuf1, num_const + 1, vecbuf2, matbuf1, 0,
                   (10e-2) * tolerance, default_num_improve);
      /*  -------- */
      /*  Test for singular system matrix */
      if (err != NULL)
        go_on = 6;
      add_vec(0.0, yvec, 1.0, vecbuf2, num_const, dpyvec);
      dptau = vecbuf2[num_const + 1];
      /*  Get dx  ,dz  ,dkappa. */
      delta_z(const_mat, num_const, zmat, cmat, dim, yvec, dpyvec, tau, eta,
              dptau, dpzmat);
      delta_x(big_inv_E, xmat, dpxmat, zmat, dpzmat, dczmat, dim, sigma, mu,
              "Predictor", dpxmat, default_num_improve);
      dpkappa = (rc_val(sigma, mu, tau, kappa) - kappa * dptau) / tau;
      /*  Define the primal and dual step length. */
      alphap =
          step_length(xmat, dpxmat, dim, tau, dptau, kappa, dpkappa, gamma);
      betap = step_length(zmat, dpzmat, dim, tau, dptau, kappa, dpkappa, gamma);
      set_steps(const_mat, const_val, num_const, yvec, dpyvec, xmat, dpxmat,
                zmat, dpzmat, cmat, dim, kappa, dpkappa, tau, dptau, eta,
                &alphap, &betap, line_reduc, "Predictor");
      if ((alphap == 0) || (betap == 0))
        go_on = 3;
    }
    /*  *************************** CORRECTOR STEP ********************** */
    if ((go_on == 0) && (feasibility_step == 0)) {
      sigma = compute_sigma(xmat, dpxmat, zmat, dpzmat, dim, alphap, betap, tau,
                            dptau, kappa, dpkappa, expon);
      eta = 1 - sigma;
      /*  Update gamma */
      gamma = frac_gamma_reduction +
              0.1 * frac_gamma_reduction * (-d_max(-alphap, -betap));
      /*  Then apply corrector step as for predictor with a few add-ins. */
      /*  Solve for (dy  ,dtau). */
      right_vector(const_mat, num_const, big_inv_E, xmat, zmat, cmat, dpxmat,
                   dpzmat, "Corrector", dim, yvec, const_val, eta, tau, sigma,
                   mu, kappa, dptau, dpkappa, vecbuf1, default_num_improve);
      /*  System matrix is the same as for the predictor (stored in systmat)  ,
       * matbuf1 is the inverse for further usage. */
      /*  Solves the main system for (dy  ,dtau) */
      /*  If error small enough  , system is dimension num_cont instead of
       * num_const+1 */
      /*  for path-following. */
      err = lsolve(systmat, vecbuf1, num_const + 1, vecbuf2, matbuf1, 0,
                   (10e-2) * tolerance, default_num_improve);
      /*  ---------- */
      add_vec(0.0, yvec, 1.0, vecbuf2, num_const, dcyvec);
      dctau = vecbuf2[num_const + 1];
      /*  Get dx  ,dz  ,dkappa. */
      delta_z(const_mat, num_const, zmat, cmat, dim, yvec, dcyvec, tau, eta,
              dctau, dczmat);
      delta_x(big_inv_E, xmat, dpxmat, zmat, dpzmat, dczmat, dim, sigma, mu,
              "Corrector", dcxmat, default_num_improve);
      /*  For match with mathematica ******* DEBUG ********** */
      /* dckappa=(rc_val(sigma  ,mu  ,tau
       * ,kappa)-dptau*dpkappa-kappa*dctau)/tau; */
      dckappa = (rc_val(sigma, mu, tau, kappa) - kappa * dctau) / tau;
      /*  Define the primal and dual step length. */
      alphac =
          step_length(xmat, dcxmat, dim, tau, dctau, kappa, dckappa, gamma);
      betac = step_length(zmat, dczmat, dim, tau, dctau, kappa, dckappa, gamma);
      set_steps(const_mat, const_val, num_const, yvec, dcyvec, xmat, dcxmat,
                zmat, dczmat, cmat, dim, kappa, dckappa, tau, dctau, eta,
                &alphac, &betac, line_reduc, "Corrector");
      if ((alphac == 0) || (betac == 0))
        go_on = 4;
    }
    /*  ********* STEP ******** */
    if (go_on == 0) {
      /*  Update tau. */
      primal_tau = tau;
      dual_tau = tau;
      primal_tau += alphac * dctau;
      dual_tau += betac * dctau;
      if (primal_tau > dual_tau) {
        kappa += alphac * dckappa;
      } else {
        kappa += betac * dckappa;
      }
      tau = d_max(primal_tau, dual_tau);
      add_mat(tau / primal_tau, xmat, (tau / primal_tau) * alphac, dcxmat, dim,
              xmat, "No");
      add_mat(tau / dual_tau, zmat, (tau / dual_tau) * betac, dczmat, dim, zmat,
              "No");
      add_vec(tau / dual_tau, yvec, (tau / dual_tau) * betac, dcyvec, num_const,
              yvec);
    }
    /*  For DEBUG
     * ****************************************************************** */
    visu_mat(xmat, dim, dim, 2);
    visu_mat(dcxmat, dim, dim, 4);
    visu_mat(zmat, dim, dim, 5);
    visu_mat(dczmat, dim, dim, 7);
    /*  ****************************************************************************
     */
    /*  Test for max iteration. */
    if (iter >= max_iter)
      go_on = 5;
    /*  ***************************** End of Main loop ******************** */
  }
  /*  Error handling and return; */
  switch (go_on) {
  case 0:
    err = serror("completely fucked up! Wake up the the author /n");
    break;
  case 1:
    err = serror("succeeded  , it's enormous... \n");
    break;
  case 2:
    err = serror("infeasibility detected. ");
    break;
  case 3:
    err = serror("Predictor hit boundary. ");
    break;
  case 4:
    err = serror("Corrector hit boundary. ");
    break;
  case 5:
    err = serror("Max iterations reached. ");
    break;
  case 6:
    err = serror("System matrix singular. ");
    break;
  case 7:
    err = serror("System operator singular. ");
    break;
  case 8:
    err = serror("Constraints not independent. ");
    break;
  }
  time(&finish_time);
  if (printlevel >= 1) {
    if (printlevel == 2) {
      /*  Status display and tracking. In debug send the status_tracker to the
       * monitor. */
      status(iter, const_mat, const_val, num_const, systmat, xmat, zmat, cmat,
             dim, yvec, tau, kappa, sigma, eta, alphac, betac, status_tracker,
             4, max_iter, numer_error);
    }
    smessage("SDP:  %s", err);
    if (printlevel >= 2) {
      smessage("Computing time: %6.1f secs.\n",
               difftime(finish_time, start_time));
    }
    smessage("\n");
  }
  /*  Return result. */
  *obj_val = scal_mat(cmat, xmat, dim) / tau;
  *return_error = go_on;
  *error = global_error(const_mat, const_val, num_const, xmat, zmat, cmat, dim,
                        yvec, tau, kappa);
  /*  Free */
  free_dmatrix(zmat, 1, dim, 1, dim);
  free_dmatrix(dpzmat, 1, dim, 1, dim);
  free_dmatrix(dpxmat, 1, dim, 1, dim);
  free_dmatrix(dczmat, 1, dim, 1, dim);
  free_dmatrix(dcxmat, 1, dim, 1, dim);
  free_dmatrix(gram, 1, num_const, 1, num_const);
  free_dmatrix(inverse_gram, 1, num_const, 1, num_const);
  free_dmatrix(big_inv_E, 1, dim * dim, 1, dim * dim);
  free_dvector(yvec, 1, num_const);
  free_dvector(dpyvec, 1, num_const);
  free_dvector(dcyvec, 1, num_const);
  free_dmatrix(status_tracker, 0, max_iter, 1, 11);
  free_dmatrix(matbuf1, 1, num_const + 1, 1, num_const + 1);
  free_dmatrix(systmat, 1, num_const + 1, 1, num_const + 1);
  free_dvector(vecbuf1, 0, num_const + 1);
  free_dvector(vecbuf2, 0, num_const + 1);
  free_dvector(bufvec1, 0, dim);
  /*  --------- */
  return err;
}

Err normalizevarcovar(double **matrix, int dim, double tau) {
  Err err = NULL;
  int i, j;

  for (i = 1; i <= dim; i++) {
    for (j = 1; j <= dim; j++) {
      matrix[i][j] = matrix[i][j] / tau;
    }
  }

  return err;
}

Err srt_bgm_sdp(double ***const_mat, double *const_val, double **cmat, int dim,
                int num_const, double **xmat, double *obj_val, double *error,
                int *return_error, double toler, int printlevel,
                int max_iterations) {
  /*  General parameters for the algorithm. */
  int max_iter = max_iterations;
  double tolerance = toler;
  double start_gamma = 0.8;
  double frac_gamma_reduction = 0.65;
  double expon = 3;
  double first_kappa = 1.0;
  double first_tau = 1.0;
  double line_reduc = 0.6;
  int default_num_improve = 0;
  double condition_bound = 1e+35;
  double tolerance_infeasibility = 10e+03 * tolerance;

  /*  ---------- */
  /*  Algortihm variables. */
  Err err = NULL;
  int go_on = 0;
  int iter = 0;
  double **zmat = dmatrix(1, dim, 1, dim);
  double **dpzmat = dmatrix(1, dim, 1, dim);
  double **dpxmat = dmatrix(1, dim, 1, dim);
  double **dczmat = dmatrix(1, dim, 1, dim);
  double **dcxmat = dmatrix(1, dim, 1, dim);
  double **gram = dmatrix(1, num_const, 1, num_const);
  double **inverse_gram = dmatrix(1, num_const, 1, num_const);
  double **big_inv_E = dmatrix(1, dim * dim, 1, dim * dim);
  double *yvec = dvector(1, num_const);
  double *dpyvec = dvector(1, num_const);
  double *dcyvec = dvector(1, num_const);
  double gamma = start_gamma;
  double tau = first_tau, dptau, dctau;
  double primal_tau, dual_tau;
  double kappa = first_kappa, dpkappa, dckappa;
  double alphap, betap, alphac = 1.0, betac = 1.0;
  double mu, first_mu;
  double phi;
  double sigma;
  double eta;
  /* ----------- */
  /*  Working variables. */
  double **status_tracker = dmatrix(0, max_iter, 1, 11);
  double **matbuf1 = dmatrix(1, num_const + 1, 1, num_const + 1);
  double **systmat = dmatrix(1, num_const + 1, 1, num_const + 1);
  double *vecbuf1 = dvector(0, num_const + 1);
  double *vecbuf2 = dvector(0, num_const + 1);
  double *bufvec1 = dvector(0, dim);
  double numer_error, bufpinf, bufdinf;
  time_t start_time, finish_time;
  int feasibility_step = 0;
  static int num_calls = 0;

  /*  Initializing. */
  time(&start_time);
  if ((num_calls == 0) && (printlevel >= 1)) {
    smessage("SDP  , Test version. \n");
    smessage("Go ahead  , make my day... \n");
    smessage("\n");
  }
  /*  Starting point. Initialization. */
  num_calls++;
  if (printlevel >= 2)
    smessage("SDP: Initializing... Call number: %d \n", num_calls);
  init_sol(const_mat, num_const, const_val, cmat, dim, xmat, zmat);
  /*  Compute the inverse Gram matrix for the A operator  */
  /*  for use in the projection on the fesible set. */
  compute_gram(const_mat, num_const, dim, gram);
  err = inv_mat(gram, num_const, inverse_gram);
  mat_mult_nobuf(gram, inverse_gram, num_const, matbuf1);
  /*  Test for constraints independence. */
  if (err != NULL)
    go_on = 8;
  /* 	   ************************ Start of main loop. ******************* */
  while (go_on == 0) {
    iter++;
    mu = compute_mu(xmat, zmat, dim, tau, kappa);
    phi = compute_infeas(const_mat, const_val, num_const, xmat, zmat, cmat, dim,
                         yvec, tau, &bufpinf, &bufdinf);
    /*  Status display and tracking. In debug send the status_tracker to the
     * monitor. */
    status(iter, const_mat, const_val, num_const, systmat, xmat, zmat, cmat,
           dim, yvec, tau, kappa, sigma, eta, alphac, betac, status_tracker,
           printlevel, max_iter, numer_error);
    /*  Send matrix to monitor if in debug. **** DEBUG ****** */
    visu_mat(status_tracker, max_iter, 11, 1);
    /*  ---------------------------- */
    if (iter == 1)
      first_mu = mu;
    /*  Test for success. */
    if (global_error(const_mat, const_val, num_const, xmat, zmat, cmat, dim,
                     yvec, tau, kappa) < tolerance)
      go_on = 1;
    /*  Or for infeasibility... */
    if ((fabs(mu / first_mu) < tolerance_infeasibility) &&
        ((tau / kappa) / (first_tau / first_kappa) < tolerance_infeasibility))
      go_on = 2;
    /*  Then test for increasing primal infeasibility. */
    if ((status_tracker[iter - 1][2] <= status_tracker[iter][2]) &&
        (feasibility_step == 0) && (iter >= 10))
      feasibility_step = 0;
    else
      feasibility_step = 0;
    /*  *************************** FEASIBILITY STEP ******************* */
    if (feasibility_step == 1) {
      take_feasibility_step(const_mat, const_val, num_const, xmat, dcxmat, zmat,
                            dczmat, dim, yvec, dcyvec, &dctau, tau, &dckappa,
                            inverse_gram);
      alphac =
          step_length(xmat, dcxmat, dim, tau, dctau, kappa, dckappa, gamma);
      betac = step_length(zmat, dczmat, dim, tau, dctau, kappa, dckappa, gamma);
      set_steps(const_mat, const_val, num_const, yvec, dcyvec, xmat, dcxmat,
                zmat, dczmat, cmat, dim, kappa, dckappa, tau, dctau, eta,
                &alphac, &betac, line_reduc, "Corrector");
    }
    /*  *************************** PREDICTOR STEP ********************* */
    if ((go_on == 0) && (feasibility_step == 0)) {
      /*  Predictor step: sigma=0  , eta=1. */
      sigma = 0;
      eta = 1;
      /*  Solve for (dy  ,dtau). */
      err =
          system_matrix(const_mat, num_const, xmat, zmat, cmat, dim, "AHO",
                        systmat, big_inv_E, &numer_error, default_num_improve);
      /*  Test for singular system operator; */
      if (err != NULL)
        go_on = 7;
      left_matrix(systmat, num_const, const_val, kappa, tau, systmat);
      right_vector(const_mat, num_const, big_inv_E, xmat, zmat, cmat, dpxmat,
                   dpzmat, "Predictor", dim, yvec, const_val, eta, tau, sigma,
                   mu, kappa, dptau, dpkappa, vecbuf1, default_num_improve);
      /*  Solve the main system for (dy.dtau): */
      /*  Test the condition number  , if to bad  , stops. */
      if (condition_mat(systmat, num_const + 1) > condition_bound)
        go_on = 6;
      err = lsolve(systmat, vecbuf1, num_const + 1, vecbuf2, matbuf1, 0,
                   (10e-2) * tolerance, default_num_improve);
      /*  -------- */
      /*  Test for singular system matrix */
      if (err != NULL)
        go_on = 6;
      add_vec(0.0, yvec, 1.0, vecbuf2, num_const, dpyvec);
      dptau = vecbuf2[num_const + 1];
      /*  Get dx  ,dz  ,dkappa. */
      delta_z(const_mat, num_const, zmat, cmat, dim, yvec, dpyvec, tau, eta,
              dptau, dpzmat);
      delta_x(big_inv_E, xmat, dpxmat, zmat, dpzmat, dczmat, dim, sigma, mu,
              "Predictor", dpxmat, default_num_improve);
      dpkappa = (rc_val(sigma, mu, tau, kappa) - kappa * dptau) / tau;
      /*  Define the primal and dual step length. */
      alphap =
          step_length(xmat, dpxmat, dim, tau, dptau, kappa, dpkappa, gamma);
      betap = step_length(zmat, dpzmat, dim, tau, dptau, kappa, dpkappa, gamma);
      set_steps(const_mat, const_val, num_const, yvec, dpyvec, xmat, dpxmat,
                zmat, dpzmat, cmat, dim, kappa, dpkappa, tau, dptau, eta,
                &alphap, &betap, line_reduc, "Predictor");
      if ((alphap == 0) || (betap == 0))
        go_on = 3;
    }
    /*  *************************** CORRECTOR STEP ********************** */
    if ((go_on == 0) && (feasibility_step == 0)) {
      sigma = compute_sigma(xmat, dpxmat, zmat, dpzmat, dim, alphap, betap, tau,
                            dptau, kappa, dpkappa, expon);
      eta = 1 - sigma;
      /*  Update gamma */
      gamma = frac_gamma_reduction +
              0.1 * frac_gamma_reduction * (-d_max(-alphap, -betap));
      /*  Then apply corrector step as for predictor with a few add-ins. */
      /*  Solve for (dy  ,dtau). */
      right_vector(const_mat, num_const, big_inv_E, xmat, zmat, cmat, dpxmat,
                   dpzmat, "Corrector", dim, yvec, const_val, eta, tau, sigma,
                   mu, kappa, dptau, dpkappa, vecbuf1, default_num_improve);
      /*  System matrix is the same as for the predictor (stored in systmat)  ,
       * matbuf1 is the inverse for further usage. */
      /*  Solves the main system for (dy  ,dtau) */
      /*  If error small enough  , system is dimension num_cont instead of
       * num_const+1 */
      /*  for path-following. */
      err = lsolve(systmat, vecbuf1, num_const + 1, vecbuf2, matbuf1, 0,
                   (10e-2) * tolerance, default_num_improve);
      /*  ---------- */
      add_vec(0.0, yvec, 1.0, vecbuf2, num_const, dcyvec);
      dctau = vecbuf2[num_const + 1];
      /*  Get dx  ,dz  ,dkappa. */
      delta_z(const_mat, num_const, zmat, cmat, dim, yvec, dcyvec, tau, eta,
              dctau, dczmat);
      delta_x(big_inv_E, xmat, dpxmat, zmat, dpzmat, dczmat, dim, sigma, mu,
              "Corrector", dcxmat, default_num_improve);
      /*  For match with mathematica ******* DEBUG ********** */
      /* dckappa=(rc_val(sigma  ,mu  ,tau
       * ,kappa)-dptau*dpkappa-kappa*dctau)/tau; */
      dckappa = (rc_val(sigma, mu, tau, kappa) - kappa * dctau) / tau;
      /*  Define the primal and dual step length. */
      alphac =
          step_length(xmat, dcxmat, dim, tau, dctau, kappa, dckappa, gamma);
      betac = step_length(zmat, dczmat, dim, tau, dctau, kappa, dckappa, gamma);
      set_steps(const_mat, const_val, num_const, yvec, dcyvec, xmat, dcxmat,
                zmat, dczmat, cmat, dim, kappa, dckappa, tau, dctau, eta,
                &alphac, &betac, line_reduc, "Corrector");
      if ((alphac == 0) || (betac == 0))
        go_on = 4;
    }
    /*  ********* STEP ******** */
    if (go_on == 0) {
      /*  Update tau. */
      primal_tau = tau;
      dual_tau = tau;
      primal_tau += alphac * dctau;
      dual_tau += betac * dctau;
      if (primal_tau > dual_tau) {
        kappa += alphac * dckappa;
      } else {
        kappa += betac * dckappa;
      }
      tau = d_max(primal_tau, dual_tau);
      add_mat(tau / primal_tau, xmat, (tau / primal_tau) * alphac, dcxmat, dim,
              xmat, "No");
      add_mat(tau / dual_tau, zmat, (tau / dual_tau) * betac, dczmat, dim, zmat,
              "No");
      add_vec(tau / dual_tau, yvec, (tau / dual_tau) * betac, dcyvec, num_const,
              yvec);
    }
    /*  For DEBUG
     * ****************************************************************** */
    visu_mat(xmat, dim, dim, 2);
    visu_mat(dcxmat, dim, dim, 4);
    visu_mat(zmat, dim, dim, 5);
    visu_mat(dczmat, dim, dim, 7);
    /*  ****************************************************************************
     */
    /*  Test for max iteration. */
    if (iter >= max_iter)
      go_on = 5;
    /*  ***************************** End of Main loop ******************** */
  }
  /*  Error handling and return; */
  switch (go_on) {
  case 0:
    err = serror("completely fucked up! Wake up the the author /n");
    break;
  case 1:
    err = serror("succeeded  , it's enormous... \n");
    break;
  case 2:
    err = serror("infeasibility detected. ");
    break;
  case 3:
    err = serror("Predictor hit boundary. ");
    break;
  case 4:
    err = serror("Corrector hit boundary. ");
    break;
  case 5:
    err = serror("Max iterations reached. ");
    break;
  case 6:
    err = serror("System matrix singular. ");
    break;
  case 7:
    err = serror("System operator singular. ");
    break;
  case 8:
    err = serror("Constraints not independent. ");
    break;
  }
  time(&finish_time);
  if (printlevel >= 1) {
    if (printlevel == 2) {
      /*  Status display and tracking. In debug send the status_tracker to the
       * monitor. */
      status(iter, const_mat, const_val, num_const, systmat, xmat, zmat, cmat,
             dim, yvec, tau, kappa, sigma, eta, alphac, betac, status_tracker,
             4, max_iter, numer_error);
    }
    smessage("SDP:  %s", err);
    if (printlevel >= 2) {
      smessage("Computing time: %6.1f secs.\n",
               difftime(finish_time, start_time));
    }
    smessage("\n");
  }
  /*  Return result. */
  *obj_val = scal_mat(cmat, xmat, dim) / tau;
  *return_error = go_on;
  *error = global_error(const_mat, const_val, num_const, xmat, zmat, cmat, dim,
                        yvec, tau, kappa);
  /*****************************Added line from original
   * version***********************/
  /*Renormalize the variance covariance matrix obtained*/
  normalizevarcovar(xmat, dim, tau);
  /*End of changes*/
  /*  Free */
  free_dmatrix(zmat, 1, dim, 1, dim);
  free_dmatrix(dpzmat, 1, dim, 1, dim);
  free_dmatrix(dpxmat, 1, dim, 1, dim);
  free_dmatrix(dczmat, 1, dim, 1, dim);
  free_dmatrix(dcxmat, 1, dim, 1, dim);
  free_dmatrix(gram, 1, num_const, 1, num_const);
  free_dmatrix(inverse_gram, 1, num_const, 1, num_const);
  free_dmatrix(big_inv_E, 1, dim * dim, 1, dim * dim);
  free_dvector(yvec, 1, num_const);
  free_dvector(dpyvec, 1, num_const);
  free_dvector(dcyvec, 1, num_const);
  free_dmatrix(status_tracker, 0, max_iter, 1, 11);
  free_dmatrix(matbuf1, 1, num_const + 1, 1, num_const + 1);
  free_dmatrix(systmat, 1, num_const + 1, 1, num_const + 1);
  free_dvector(vecbuf1, 0, num_const + 1);
  free_dvector(vecbuf2, 0, num_const + 1);
  free_dvector(bufvec1, 0, dim);
  /*  --------- */
  return err;
}