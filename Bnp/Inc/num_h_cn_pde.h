/* -------------------------------------------------------------------------
   FILENAME		: num_h_cn_pde.h
   AUTHOR		: Toto 15Dec1999

   PURPOSE		: convolve backward/forward using Cranck-Nicholson
   equations
   ------------------------------------------------------------------------- */

#ifndef NUM_H_CN_PDE_H
#define NUM_H_CN_PDE_H

typedef struct {
  double *exp_pu;
  double *exp_pc;
  double *exp_pd;

  double *imp_pu;
  double *imp_pc;
  double *imp_pd;

  double **right_mem;

  double *gam;
} CNPDE_TEMP;

/*	Call before PDE	*/
void num_f_pde_init(CNPDE_TEMP *tmp,
                    /*	Size of the array	*/
                    int n,
                    /*	Number of functions	*/
                    int m);

/*	Call after PDE	*/
void num_f_pde_free(CNPDE_TEMP *tmp,
                    /*	Size of the array	*/
                    int n,
                    /*	Number of functions	*/
                    int m);

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
    /*	i: state      , j: function number */
    double **f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
    double *mu,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double *var,
    /*	Discount rates from t to t + 1	*/
    double *r,
    /*	Cranck-Nicholson theta	*/
    double theta,
    /*	Result: f values @ t	*/
    double **f_t);

/*	De-convolve one step forward using Kolmogorov forward equation	*/
void num_f_pde_one_step_forward(
    CNPDE_TEMP *tmp,
    /*	Size of the array	*/
    int n,
    /*	X values	*/
    double *x,
    /*	f values @ t	*/
    double *f_t,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
    double *mu,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double *var,
    /*	Discount rates from t to t + 1	*/
    double *r,
    /*	Cranck-Nicholson theta	*/
    double theta,
    /*	Result: f values @ t+1	*/
    double *f_t_plus_1);

/*	2D Scheme */

typedef struct {
  CNPDE_TEMP *tmpx;
  CNPDE_TEMP *tmpy;
} CNPDE_TEMP_2D_LOD;

/*	Call before PDE	2D */
void num_f_pde_init_2d_lod(CNPDE_TEMP_2D_LOD *tmp,
                           /*	Size of the arrays	*/
                           int nx, int ny,
                           /*	Number of functions	*/
                           int m);

/*	Call after PDE 2D */
void num_f_pde_free_2d_lod(CNPDE_TEMP_2D_LOD *tmp,
                           /*	Size of the arrays	*/
                           int nx, int ny,
                           /*	Number of functions	*/
                           int m);

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
    /*	i: x state      , j: y state      , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Cranck-Nicholson theta	*/
    double theta,
    /*	Result: f values @ t	*/
    double ***f_t, int lx, int ux, int ly, int uy);

typedef struct {
  double **xexp_pu;
  double **xexp_pc;
  double **xexp_pd;

  double **ximp_pu;
  double **ximp_pc;
  double **ximp_pd;

  double **yexp_pu;
  double **yexp_pc;
  double **yexp_pd;

  double **yimp_pu;
  double **yimp_pc;
  double **yimp_pd;

  double ***f;

  double *gamx;
  double *gamy;
} CNPDE_TEMP_2D_ADI;

/*	Call before PDE	2D */
void num_f_pde_init_2d_adi(CNPDE_TEMP_2D_ADI *tmp,
                           /*	Size of the arrays	*/
                           int nx, int ny,
                           /*	Number of functions	*/
                           int m);

/*	Call after PDE 2D */
void num_f_pde_free_2d_adi(CNPDE_TEMP_2D_ADI *tmp,
                           /*	Size of the arrays	*/
                           int nx, int ny,
                           /*	Number of functions	*/
                           int m);

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
    /*	i: x state      , j: y state      , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy);

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
    /*	i: x state      , j: y state      , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy, int recalc);

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
    /*	i: x state      , j: y state      , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int *ly, int *uy,

    int *invlx, int *invux, int invly, int invuy);

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
    /*	i: x state      , j: y state      , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy);

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
    /*	i: x state      , j: y state      , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy);

/*	Convolve one step backward using ADI scheme with barrier */
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
    /*	i: x state      , j: y state      , k: function number */
    double ***f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
    double **mux, double **muy,
    /*	Variances under Qt+1 * dt (actually Vt+1 [Xt+1 - Xt / Ft] )	*/
    double **varx, double **vary,
    /*	Discount rates from t to t + 1	*/
    double **r,
    /*	Result: f values @ t	*/
    double ***f_t,
    /*	Lower and upper x and y */
    int lx, int ux, int ly, int uy, int is_up, int is_down);

typedef struct {
  double **xexp_pu;
  double **xexp_pc;
  double **xexp_pd;

  double **ximp_pu;
  double **ximp_pc;
  double **ximp_pd;

  double **yexp_pu;
  double **yexp_pc;
  double **yexp_pd;

  double **yimp_pu;
  double **yimp_pc;
  double **yimp_pd;

  double ****f;

  double *gamx;
  double *gamy;
} CNPDE_TEMP_2D_LGMSV_ADI;

/*	Call before PDE	2D */
void num_f_pde_init_2d_LGMSV_adi(CNPDE_TEMP_2D_LGMSV_ADI *tmp,
                                 /*	Size of the arrays	*/
                                 int npsi, int nx, int ny,
                                 /*	Number of functions	*/
                                 int m);

/*	Call after PDE 2D */
void num_f_pde_free_2d_LGMSV_adi(CNPDE_TEMP_2D_LGMSV_ADI *tmp,
                                 /*	Size of the arrays	*/
                                 int npsi, int nx, int ny,
                                 /*	Number of functions	*/
                                 int m);

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
    /*	i: x state      , j: y state      , k: function number */
    double ****f_t_plus_1,
    /*	Drifts under Qt+1 * dt (actually      , Et+1 [Xt+1 - Xt / Ft] )	*/
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
    int lx, int ux, int ly, int uy);

#endif
