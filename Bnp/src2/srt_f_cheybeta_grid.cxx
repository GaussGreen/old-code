/* ==============================================================================

   FILE NAME:      srt_f_cheybeta_grid.cxx

   OBJECT:         functions implementing the local vol Cheyette PDE grid

  ===============================================================================
*/

#include "math.h"
#include "srt_h_all.h"

/*	Grid definition	*/

/*	Initialise        , i.e. set pointers to NULL	*/
void chey_beta_grid_init(CHEYBETA_GRID *grid) {
  /*	Just set pointers to NULL	*/
  grid->nt = 0;
  grid->t = NULL;

  grid->maxvar = NULL;
  grid->minvar = NULL;

  grid->fwdrt = NULL;

  grid->nx = 0;
  grid->x = NULL;

  grid->maxnphi = 0;
  grid->nphi = NULL;
  grid->phi = NULL;

  grid->adt1 = NULL;
  grid->adt2 = NULL;

  grid->num_prod = 0;

  grid->vt1 = NULL;
  grid->vt2 = NULL;
}

/*	Free	*/
void chey_beta_grid_free(CHEYBETA_GRID *grid) {
  int i, j;

  if (grid->nt) {
    if (grid->nphi) {
      for (i = 0; i < grid->nt; i++) {
        if (grid->nphi + i) {
          if (grid->nphi[i]) {
            free(grid->phi[i]);
            grid->nphi[i] = 0;
            grid->phi[i] = NULL;
          }
        }
      }
      free(grid->phi);
      grid->phi = NULL;
    }

    free(grid->t);
    grid->nt = 0;
    grid->t = NULL;
  }

  if (grid->maxvar) {
    free(grid->maxvar);
    grid->maxvar = NULL;
  }

  if (grid->minvar) {
    free(grid->minvar);
    grid->minvar = NULL;
  }

  if (grid->fwdrt) {
    free(grid->fwdrt);
    grid->fwdrt = NULL;
  }

  if (grid->adt1) {
    for (i = 0; i < grid->maxnphi; i++) {
      if (grid->adt1 + i) {
        if (grid->adt1[i]) {
          free(grid->adt1[i]);
          grid->adt1[i] = NULL;
        }
      }
    }
    free(grid->adt1);
    grid->adt1 = NULL;
  }

  if (grid->adt2) {
    for (i = 0; i < grid->maxnphi; i++) {
      if (grid->adt2 + i) {
        if (grid->adt2[i]) {
          free(grid->adt2[i]);
          grid->adt2[i] = NULL;
        }
      }
    }
    free(grid->adt2);
    grid->adt2 = NULL;
  }

  if (grid->num_prod && grid->vt1) {
    for (i = 0; i < grid->maxnphi; i++) {
      if (grid->vt1 + i) {
        if (grid->vt1[i]) {
          for (j = 0; j < grid->nx; j++) {
            if (grid->vt1[i] + j) {
              if (grid->vt1[i][j]) {
                free(grid->vt1[i][j]);
                grid->vt1[i][j] = NULL;
              }
            }
          }
          free(grid->vt1[i]);
          grid->vt1[i] = NULL;
        }
      }
    }
    free(grid->vt1);
    grid->vt1 = NULL;
  }

  if (grid->num_prod && grid->vt2) {
    for (i = 0; i < grid->maxnphi; i++) {
      if (grid->vt2 + i) {
        if (grid->vt2[i]) {
          for (j = 0; j < grid->nx; j++) {
            if (grid->vt2[i] + j) {
              if (grid->vt2[i][j]) {
                free(grid->vt2[i][j]);
                grid->vt2[i][j] = NULL;
              }
            }
          }
          free(grid->vt2[i]);
          grid->vt2[i] = NULL;
        }
      }
    }
    free(grid->vt2);
    grid->vt2 = NULL;
  }

  grid->num_prod = 0;

  if (grid->nx) {
    free(grid->x);
    grid->nx = 0;
    grid->x = NULL;
  }

  grid->maxnphi = 0;
}

/*	Construct time grid	*/
void chey_beta_grid_time(
    /*	Model	*/
    CHEYBETA_MDL *mdl,
    /*	Vector of required stop times	*/
    int num_req_times, double *req_times,
    /*	Number of steps required	*/
    int nt,
    /*	Grid	*/
    CHEYBETA_GRID *grid) {
  /*	Copy required times	*/
  grid->nt = num_req_times;
  grid->t = (double *)calloc(grid->nt, sizeof(double));
  memcpy(grid->t, req_times, grid->nt * sizeof(double));

  /*	Add today	*/
  num_f_add_number(&(grid->nt), &(grid->t), 0.0);
  /*	Sort	*/
  num_f_sort_vector(grid->nt, grid->t);
  /*	Remove doubles	*/
  num_f_unique_vector(&(grid->nt), grid->t);
  /*	Fill evenly to the required number of steps	*/
  num_f_fill_vector(grid->nt, &(grid->t), nt, nt / 3);
  grid->nt = nt;

  /*	Allocate time dependant parameters in grid	*/
  grid->maxvar = (double *)calloc(grid->nt, sizeof(double));
  grid->minvar = (double *)calloc(grid->nt, sizeof(double));
  grid->fwdrt = (double *)calloc(grid->nt, sizeof(double));
  grid->nphi = (int *)calloc(grid->nt, sizeof(int));
  grid->phi = (double **)calloc(grid->nt, sizeof(double *));
}

/*	Construct phi and x grids	*/
void chey_beta_grid_phi_x(
    /*	Model	*/
    CHEYBETA_MDL *mdl,
    /*	Grid	*/
    CHEYBETA_GRID *grid,
    /*	Num x	*/
    int nx,
    /*	Num phi	*/
    int nphi) {
  int t1, t2;
  double phi_min, phi_mid, phi_max;
  double x_min, x_mid, x_max;
  double v_min, v_mid, v_max;
  CHEYBETA_PARAM p;
  double cum_var, delta_x;

  /*	Initialise phi and x	*/

  phi_min = 0.0;
  phi_mid = 0.0;
  phi_max = 0.0;

  x_min = 0.0;
  x_mid = 0.0;
  x_max = 0.0;

  /*	Phi grid at time 0	*/
  grid->nphi[0] = 1;
  grid->phi[0] = (double *)calloc(1, sizeof(double));
  grid->phi[0][0] = 0.0;
  grid->maxnphi = 1;

  /*	Initialise cum_var        , exp and std	*/
  cum_var = 0.0;

  /*	Recursively construct phi grid	*/
  for (t1 = 0, t2 = 1; t2 < grid->nt; t1++, t2++) {
    /*	Get parameters between t1 and t2	*/
    chey_beta_mdl_param(mdl, grid->t[t1], grid->t[t2], &p);

    /*	Calculate forward rate from t1 to t2	*/
    grid->fwdrt[t1] =
        swp_f_zr(mdl->today_date + DAYS_IN_YEAR * grid->t[t1],
                 mdl->today_date + DAYS_IN_YEAR * grid->t[t2], mdl->yc_name);

    /*	Update expectation and variances	*/

    v_mid = chey_beta_mdl_norm_var_at_t(&p, x_mid, phi_mid, grid->fwdrt[t1],
                                        DBL_MAX, 0.0);

    grid->maxvar[t1] = chey_beta_mdl_norm_var_at_t(
        &p, x_max, phi_max, grid->fwdrt[t1], DBL_MAX, 0.0);

    grid->minvar[t1] = 2 * v_mid - grid->maxvar[t1];

    if (grid->minvar[t1] < 1.0e-08) {
      grid->minvar[t1] = 1.0e-08;
    }

    cum_var += v_mid * p.dt;

    v_mid = chey_beta_mdl_var(mdl, &p, x_mid, phi_mid, v_mid);

    v_max = chey_beta_mdl_var(mdl, &p, x_max, phi_max, grid->maxvar[t1]);

    v_min = chey_beta_mdl_var(mdl, &p, x_min, phi_min, grid->minvar[t1]);

    /*	Update x_mid        , m_min = x_mid -2std        , x_max = m_mid+2std
     */
    /*	Also calculate phi along mid        , min and max lines	*/

    delta_x = chey_beta_mdl_drift(mdl, &p, x_mid, phi_mid);

    phi_mid = chey_beta_mdl_phi(mdl, &p, x_mid, phi_mid, v_mid);
    phi_max = chey_beta_mdl_phi(mdl, &p, x_max, phi_max, v_max);
    phi_min = chey_beta_mdl_phi(mdl, &p, x_min, phi_min, v_min);

    x_mid += delta_x;
    x_max = x_mid + 2 * sqrt(cum_var);
    x_min = x_mid - 2 * sqrt(cum_var);

    /*	Construct phi grid at t	*/

    if (phi_max - phi_mid > 1.0e-8) {
      grid->nphi[t2] =
          (int)(3 + (nphi - 3) * grid->t[t2] / grid->t[grid->nt - 1]);
      /*	Force nphi to be an odd number	*/
      grid->nphi[t2] = (int)(grid->nphi[t2] / 2) * 2 + 1;
      grid->phi[t2] = (double *)calloc(3, sizeof(double));
      grid->phi[t2][0] = phi_min;
      grid->phi[t2][1] = phi_mid;
      grid->phi[t2][2] = phi_max;
      if (grid->nphi[t2] < 3) {
        grid->nphi[t2] = 3;
      }
      num_f_fill_vector(3, &(grid->phi[t2]), grid->nphi[t2], 0);
    } else
    /*	If phi is constant (i.e. LGM)	*/
    {
      grid->nphi[t2] = 1;
      grid->phi[t2] = (double *)calloc(1, sizeof(double));
      grid->phi[t2][0] = phi_mid;
    }

    if (grid->maxnphi < grid->nphi[t2]) {
      grid->maxnphi = grid->nphi[t2];
    }
  }

  /*	Force nx to be an odd number	*/
  grid->nx = (int)(nx / 2) * 2 + 1;

  /*	Make x grid	*/
  grid->x = (double *)calloc(2, sizeof(double));
  grid->x[0] = -5.0 * sqrt(cum_var);
  grid->x[1] = 8.0 * sqrt(cum_var);
  num_f_fill_vector(2, &(grid->x), grid->nx, 0);
}

/*	For the forward PDE	*/

/*	Allocate AD grids	*/
void chey_beta_grid_ad(
    /*	Grid	*/
    CHEYBETA_GRID *grid) {
  int i;

  grid->adt1 = (double **)calloc(grid->maxnphi, sizeof(double *));
  for (i = 0; i < grid->maxnphi; i++) {
    grid->adt1[i] = (double *)calloc(grid->nx, sizeof(double));
  }

  grid->adt2 = (double **)calloc(grid->maxnphi, sizeof(double *));
  for (i = 0; i < grid->maxnphi; i++) {
    grid->adt2[i] = (double *)calloc(grid->nx, sizeof(double));
  }
}

/*	Copy density	*/
void chey_beta_density_at_t(
    /*	Grid	*/
    CHEYBETA_GRID *grid,
    /*	Time step index	*/
    int ti,
    /*	X grid	*/
    int *nx, double **x,
    /*	Phi grid	*/
    int *nphi, double **phi,
    /*	AD grid	idx1 = phi        , idx2 = r	*/
    double ***ad) {
  int i;

  *nx = grid->nx;
  *x = (double *)calloc(*nx, sizeof(double));
  memcpy(*x, grid->x, *nx * sizeof(double));

  *nphi = grid->nphi[ti];
  *phi = (double *)calloc(*nphi, sizeof(double));
  memcpy(*phi, grid->phi[ti], *nphi * sizeof(double));

  *ad = (double **)calloc(*nphi, sizeof(double *));
  for (i = 0; i < *nphi; i++) {
    (*ad)[i] = (double *)calloc(*nx, sizeof(double));
    memcpy((*ad)[i], grid->adt1[i], *nx * sizeof(double));
  }
}

/*	Free density	*/
void chey_beta_free_density_at_t(
    /*	X grid	*/
    int *nx, double **x,
    /*	Phi grid	*/
    int *nphi, double **phi,
    /*	AD grid	idx1 = phi        , idx2 = r	*/
    double ***ad) {
  int i;

  free(*x);
  free(*phi);

  for (i = 0; i < *nphi; i++) {
    free((*ad)[i]);
  }
  free(*ad);
}

/*	For the backward PDE	*/

/*	Allocate value grids	*/
void chey_beta_grid_prod_val(
    /*	Grid	*/
    CHEYBETA_GRID *grid,
    /*	Number of products to be valued	*/
    int num_prod) {
  int i, j;

  grid->num_prod = num_prod;

  grid->vt1 = (double ***)calloc(grid->maxnphi, sizeof(double **));
  for (i = 0; i < grid->maxnphi; i++) {
    grid->vt1[i] = (double **)calloc(grid->nx, sizeof(double *));
    for (j = 0; j < grid->nx; j++) {
      grid->vt1[i][j] = (double *)calloc(num_prod, sizeof(double));
    }
  }

  grid->vt2 = (double ***)calloc(grid->maxnphi, sizeof(double **));
  for (i = 0; i < grid->maxnphi; i++) {
    grid->vt2[i] = (double **)calloc(grid->nx, sizeof(double *));
    for (j = 0; j < grid->nx; j++) {
      grid->vt2[i][j] = (double *)calloc(num_prod, sizeof(double));
    }
  }
}

/*	Get	product values	*/
void chey_beta_prod_val_at_t(
    /*	Grid	*/
    CHEYBETA_GRID *grid,
    /*	Time step index	*/
    int ti,
    /*	X value required	*/
    double x,
    /*	Phi value required	*/
    double phi,
    /*	Product values	*/
    int *num_prod,
    /*	Allocated inside        , to be freed by caller	*/
    double **prod_val) {
  int i, j;

  if (x <= grid->x[0]) {
    i = 0;
  } else if (x >= grid->x[grid->nx - 1]) {
    i = grid->nx - 1;
  } else {
    i = 0;
    while (grid->x[i] < x) {
      i++;
    }
  }

  if (phi <= grid->phi[ti][0]) {
    j = 0;
  } else if (phi >= grid->phi[ti][grid->nphi[ti] - 1]) {
    j = grid->nphi[ti] - 1;
  } else {
    j = 0;
    while (grid->phi[ti][j] < phi) {
      j++;
    }
  }

  *num_prod = grid->num_prod;
  *prod_val = (double *)calloc(*num_prod, sizeof(double));
  memcpy(*prod_val, grid->vt2[i][j], *num_prod * sizeof(double));
}