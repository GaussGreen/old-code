#ifndef IMSL_AT_H

#define IMSL_AT_H



#include "imsl_inc.h"



#ifdef DOUBLE

#define ANAME(FNAME, DNAME)     static char *aname[2]={DNAME,"?????"};

#else

#define ANAME(FNAME, DNAME)     static char *aname[2]={FNAME,"?????"};

#endif



#ifdef DOUBLE

#define     imsl_tchkr      imsl_dtchkr

#define     imsl_tchkc      imsl_dtchkc

#define     imsl_tchki      imsl_dtchki

#define     imsl_schkr      imsl_dschkr

#define     imsl_schkc      imsl_dschkc

#define     imsl_schki      imsl_dschki

#define     imsl_testr      imsl_dtestr

#define     imsl_testc      imsl_dtestc

#define     imsl_testn      imsl_dtestn

#define     imsl_ststr      imsl_dststr

#define     imsl_ststn      imsl_dststn



#define     imsl_diffs      imsl_ddiffs

#define     imsl_abnrm      imsl_dabnrm

#define     imsl_cbnrm      imsl_dcbnrm



/*      SPARSE MATRIX DATA      */



#define imsl_f_NOS1                     imsl_d_NOS1

#define imsl_f_662_BUS                  imsl_d_662_BUS1

#define imsl_f_IMPCOL_A                 imsl_d_IMPCOL_A

#define imsl_f_IMPCOL_B                 imsl_d_IMPCOL_B

#define imsl_f_IMPCOL_C                 imsl_d_IMPCOL_C

#define imsl_f_IMPCOL_D                 imsl_d_IMPCOL_D

#define imsl_f_IMPCOL_E                 imsl_d_IMPCOL_E



#endif



void 	PROTO(imsl_begin_test,(char*, ...));

void 	PROTO(imsl_end_test,());

Mdouble PROTO(imsl_get_chksum,());

void 	PROTO(imsl_put_chksum,(Mdouble));

void	PROTO(imsl_tchkr, (Mfloat[], Mint, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mfloat*, Mint*, Mint));

void	PROTO(imsl_tchki, (Mint[], Mint, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mfloat*, Mint*, Mint));

void	PROTO(imsl_tchkc, (Mf_complex[], Mint, Mint[], Mint, Mchar**, Mf_complex, Mint*, Mfloat, Mf_complex*, Mint*, Mint));

void	PROTO(imsl_schkr, (Mfloat, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mint*, Mint));

void	PROTO(imsl_schki, (Mint, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mint*, Mint));

void	PROTO(imsl_schkc, (Mf_complex, Mint[], Mint, Mchar**, Mf_complex, Mint*, Mfloat, Mint*, Mint));



int	PROTO(imsl_f_check_vector, (char *, int, float[], ...));

int	PROTO(imsl_d_check_vector, (char *, int, double[], ...));

int	PROTO(imsl_i_check_vector, (char *, int, int[], ...));

int	PROTO(imsl_c_check_vector, (char *, int, f_complex[], ...));

int	PROTO(imsl_z_check_vector, (char *, int, d_complex[], ...));

int	PROTO(imsl_f_check_matrix, (char *, int, int, float[], ...));

int	PROTO(imsl_d_check_matrix, (char *, int, int, double[], ...));

int	PROTO(imsl_i_check_matrix, (char *, int, int, int[], ...));

int	PROTO(imsl_c_check_matrix, (char *, int, int, f_complex[], ...));

int	PROTO(imsl_z_check_matrix, (char *, int, int, d_complex[], ...));

int	PROTO(imsl_f_check_scalar, (char *, float, ...));

int	PROTO(imsl_d_check_scalar, (char *, double, ...));

int	PROTO(imsl_i_check_scalar, (char *, int, ...));

int	PROTO(imsl_c_check_scalar, (char *, f_complex, ...));

int	PROTO(imsl_z_check_scalar, (char *, d_complex, ...));

int     PROTO(imsl_check_pointer, (char *, void *, void *, ...));





void    PROTO(imsl_ercheck, (char *, ...));

void    PROTO(imsl_chksum,(int, double));



void	PROTO(imsl_ststn,(Mint, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mint*, Mint));

void	PROTO(imsl_ststr,(Mfloat, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mint*, Mint));

void	PROTO(imsl_ststc,(Mf_complex, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mint*, Mint));

void	PROTO(imsl_testn,(Mint[], Mint, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mfloat*, Mint*, Mint));

void	PROTO(imsl_testr,(Mfloat[], Mint, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mfloat*, Mint*, Mint));

void    PROTO(imsl_testc, (Mf_complex[], Mint, Mint[], Mint, Mchar**, Mfloat, Mint*, Mfloat, Mfloat*, Mint*, Mint));



void    PROTO(imsl_diffs,(Mfloat x, Mfloat y, Mfloat ax, Mfloat *er));

void    PROTO(imsl_abnrm,(Mfloat stat[], Mint *k, Mfloat *resul));

void    PROTO(imsl_cbnrm,(Mf_complex imsl_cstat[], Mint *k, Mfloat *ctot));

void    PROTO(imsl_rititm,(Mint, Mint*));

void    PROTO(imsl_pasckd,(Mchar**, Mint));

void    PROTO(imsl_pascks,(Mchar**, Mint));

void    PROTO(imsl_pasckt,(Mchar**, Mint, Mint));

void    PROTO(imsl_e1csm,(Mint,Mdouble*));

void    PROTO(imsl_e1ucs,(Mint type, Mint code, Mchar *error_message));

void    PROTO(imsl_e1chk,(Mdouble));



	/* Chapter 1 */



void    at_lin_sol_gen();

void    at_f_lin_sol_gen();

void    at_d_lin_sol_gen();

void    at_c_lin_sol_gen();

void    at_z_lin_sol_gen();



void	at_lin_least_squares_gen();

void	at_f_lin_least_squares_gen();

void	at_d_lin_least_squares_gen();



void	at_lin_sol_nonnegdef();

void	at_f_lin_sol_nonnegdef();

void	at_d_lin_sol_nonnegdef();



void	at_f_lin_sol_posdef();

void	at_d_lin_sol_posdef();

void	at_c_lin_sol_posdef();

void	at_z_lin_sol_posdef();



void    at_lin_svd_gen();

void    at_f_lin_svd_gen();

void    at_d_lin_svd_gen();

void    at_c_lin_svd_gen();

void    at_z_lin_svd_gen();

void	at_lin_cgsol_posdef();



void	at_f_lin_cgsol_posdef();

void	at_d_lin_cgsol_posdef();



void    at_mat_mul_rect();

void    at_f_mat_mul_rect();

void    at_d_mat_mul_rect();

void    at_c_mat_mul_rect();

void    at_z_mat_mul_rect();



void	at_lin_lsq_lin_constraints();

void	at_f_lin_lsq_lin_constraints();

void	at_d_lin_lsq_lin_constraints();



	/* Chapter 2 */



void	at_eig_gen();

void	at_f_eig_gen();

void	at_d_eig_gen();

void	at_c_eig_gen();

void	at_z_eig_gen();



void	at_eig_sym();

void	at_f_eig_sym();

void	at_d_eig_sym();



void	at_eig_herm();

void	at_c_eig_herm();

void	at_z_eig_herm();



void	at_eig_symgen();

void	at_f_eig_symgen();

void	at_d_eig_symgen();



	/* Chapter 3 */

void    at_cub_spline_interp();         

void    at_f_cub_spline_interp();         

void    at_d_cub_spline_interp();



void    at_cub_spline_interp_e_cnd();

void    at_f_cub_spline_interp_e_cnd();

void    at_d_cub_spline_interp_e_cnd();



void    at_cub_spline_interp_shape();

void    at_f_cub_spline_interp_shape();

void    at_d_cub_spline_interp_shape();



void    at_cub_spline_value();

void    at_f_cub_spline_value();

void    at_d_cub_spline_value();



void    at_cub_spline_integral();       

void    at_f_cub_spline_integral();       

void    at_d_cub_spline_integral();



void    at_spline_interp();             

void    at_f_spline_interp();             

void    at_d_spline_interp(); 



void    at_spline_integral();           

void    at_f_spline_integral();           

void    at_d_spline_integral();



void    at_spline_value();              

void    at_f_spline_value();              

void    at_d_spline_value();



void    at_spline_knots();              

void    at_f_spline_knots();              

void    at_d_spline_knots();



void    at_spline_2d_interp();             

void    at_f_spline_2d_interp();             

void    at_d_spline_2d_interp(); 



void    at_spline_2d_integral();           

void    at_f_spline_2d_integral();           

void    at_d_spline_2d_integral();



void    at_spline_2d_value();              

void    at_f_spline_2d_value();              

void    at_d_spline_2d_value();



void    at_spline_2d_least_squares();              

void    at_f_spline_2d_least_squares();              

void    at_d_spline_2d_least_squares();              



void    at_spline_least_squares();              

void    at_f_spline_least_squares();              

void    at_d_spline_least_squares();              



void    at_poly_least_squares();

void    at_f_poly_least_squares();

void    at_d_poly_least_squares();



void    at_user_fcn_least_squares();

void    at_f_user_fcn_least_squares();

void    at_d_user_fcn_least_squares();



void    at_cub_spline_smooth();

void    at_f_cub_spline_smooth();

void    at_d_cub_spline_smooth();



void    at_scattered_2d_interp();

void    at_f_scattered_2d_interp();

void    at_d_scattered_2d_interp();



void    at_radial_scattered_fit();

void    at_f_radial_scattered_fit();

void    at_d_radial_scattered_fit();



void	at_radial_evaluate();

void	at_f_radial_evaluate();

void	at_d_radial_evaluate();



void	at_spline_lsq_constrained();

void	at_f_spline_lsq_constrained();

void	at_d_spline_lsq_constrained();



	/* Chapter 4 */



void	at_int_fcn_sing();

void	at_f_int_fcn_sing();

void	at_d_int_fcn_sing();



void    at_int_fcn();

void    at_f_int_fcn();

void    at_d_int_fcn();



void	at_int_fcn_sing_pts();

void	at_f_int_fcn_sing_pts();

void	at_d_int_fcn_sing_pts();



void	at_int_fcn_alg_log();

void	at_f_int_fcn_alg_log();

void	at_d_int_fcn_alg_log();



void	at_int_fcn_inf();

void	at_f_int_fcn_inf();

void	at_d_int_fcn_inf();



void	at_int_fcn_trig();

void	at_f_int_fcn_trig();

void	at_d_int_fcn_trig();



void	at_int_fcn_fourier();

void	at_f_int_fcn_fourier();

void	at_d_int_fcn_fourier();



void	at_int_fcn_cauchy();

void	at_f_int_fcn_cauchy();

void	at_d_int_fcn_cauchy();



void	at_int_fcn_smooth();

void	at_f_int_fcn_smooth();

void	at_d_int_fcn_smooth();



void    at_int_fcn_2d();

void    at_f_int_fcn_2d();

void    at_d_int_fcn_2d();



void	at_int_fcn_hyper_rect();

void	at_f_int_fcn_hyper_rect();

void	at_d_int_fcn_hyper_rect();



void    at_gauss_quad_rule();

void    at_f_gauss_quad_rule();

void    at_d_gauss_quad_rule();





	/* Chapter 5 */



void	at_ode_runge_kutta();

void	at_f_ode_runge_kutta();

void	at_d_ode_runge_kutta();



void	at_ode_adams_gear();

void	at_f_ode_adams_gear();

void	at_d_ode_adams_gear();



	/* Chapter 6 */



void    at_fft_real();

void    at_f_fft_real();

void    at_d_fft_real();



void    at_fft_real_init();

void    at_f_fft_real_init();

void    at_d_fft_real_init();



void	at_c_fft_complex();

void	at_z_fft_complex();



void	at_c_fft_complex_init();

void	at_z_fft_complex_init();



void	at_fft_2d_complex();

void	at_c_fft_2d_complex();

void	at_z_fft_2d_complex();



void	at_convolution();

void	at_f_convolution();

void	at_d_convolution();

void	at_c_convolution();

void	at_z_convolution();



	/* Chapter 7 */



void	at_zeros_poly();

void	at_f_zeros_poly();

void	at_d_zeros_poly();

void	at_c_zeros_poly();

void	at_z_zeros_poly();



void	at_zeros_fcn();

void	at_f_zeros_fcn();

void	at_d_zeros_fcn();



void	at_zeros_sys_eqn();

void	at_f_zeros_sys_eqn();

void	at_d_zeros_sys_eqn();



	/* Chapter 8 */



void	at_min_uncon();

void	at_f_min_uncon();

void	at_d_min_uncon();



void	at_min_uncon_deriv();

void	at_f_min_uncon_deriv();

void	at_d_min_uncon_deriv();



void	at_min_uncon_multivar();

void	at_f_min_uncon_multivar();

void	at_d_min_uncon_multivar();



void	at_nonlin_least_squares();

void	at_f_nonlin_least_squares();

void	at_d_nonlin_least_squares();



void	at_lin_prog();

void	at_f_lin_prog();

void	at_d_lin_prog();



void	at_quadratic_prog();

void	at_f_quadratic_prog();

void	at_d_quadratic_prog();



void	at_min_con_nonlin();

void	at_f_min_con_nonlin();

void	at_d_min_con_nonlin();



	/* Chapter 9 */



void	at_erf();

void	at_f_erf();

void	at_d_erf();



void	at_erf_inverse();

void	at_f_erf_inverse();

void	at_d_erf_inverse();



void	at_erfc();

void	at_f_erfc();

void	at_d_erfc();



void	at_erfc_inverse();

void	at_f_erfc_inverse();

void	at_d_erfc_inverse();



void	at_erfci();

void	at_f_erfci();

void	at_d_erfci();



void    at_chi_squared_inverse_cdf();

void    at_f_chi_squared_inverse_cdf();

void    at_d_chi_squared_inverse_cdf();



void    at_t_inverse_cdf();

void    at_f_t_inverse_cdf();

void    at_d_t_inverse_cdf();



void	at_F_cdf();

void	at_f_F_cdf();

void	at_d_F_cdf();



void    at_gamma();

void    at_f_gamma();

void    at_d_gamma();



void    at_beta();

void    at_f_beta();

void    at_d_beta();



void    at_log_gamma();

void    at_f_log_gamma();

void    at_d_log_gamma();



void    at_log_beta();

void    at_f_log_beta();

void    at_d_log_beta();



void    at_beta_incomplete();

void    at_f_beta_incomplete();

void    at_d_beta_incomplete();

 

void	at_binomial_cdf();

void	at_f_binomial_cdf();

void	at_d_binomial_cdf();

 

void	at_normal_cdf();

void	at_f_normal_cdf();

void	at_d_normal_cdf();

 

void	at_normal_inverse_cdf();

void	at_f_normal_inverse_cdf();

void	at_d_normal_inverse_cdf();

 

void	at_gamma_cdf();

void	at_f_gamma_cdf();

void	at_d_gamma_cdf();



void    at_gamma_incomplete();

void    at_f_gamma_incomplete();

void    at_d_gamma_incomplete();

 

void    at_t_cdf();

void    at_f_t_cdf();

void    at_d_t_cdf();



void    at_F_inverse_cdf();

void    at_f_F_inverse_cdf();

void    at_d_F_inverse_cdf();



void	at_hypergeometric_cdf();

void	at_f_hypergeometric_cdf();

void	at_d_hypergeometric_cdf();



void    at_poisson_cdf();

void    at_f_poisson_cdf();

void    at_d_poisson_cdf();



void    at_chi_squared_cdf();

void    at_f_chi_squared_cdf();

void    at_d_chi_squared_cdf();



void	at_bessel_J0();

void	at_f_bessel_J0();

void	at_d_bessel_J0();



void	at_bessel_J1();

void	at_f_bessel_J1();

void	at_d_bessel_J1();



void	at_bessel_Jx();

void	at_c_bessel_Jx();

void	at_z_bessel_Jx();



void	at_bessel_Y0();

void	at_f_bessel_Y0();

void	at_d_bessel_Y0();



void	at_bessel_Y1();

void	at_f_bessel_Y1();

void	at_d_bessel_Y1();



void	at_bessel_Yx();

void	at_c_bessel_Yx();

void	at_z_bessel_Yx();



void	at_bessel_I0();

void	at_f_bessel_I0();

void	at_d_bessel_I0();



void	at_bessel_I1();

void	at_f_bessel_I1();

void	at_d_bessel_I1();



void	at_bessel_Ix();

void	at_c_bessel_Ix();

void	at_z_bessel_Ix();



void	at_bessel_k0();

void	at_f_bessel_K0();

void	at_d_bessel_K0();



void	at_bessel_K1();

void	at_f_bessel_K1();

void	at_d_bessel_K1();



void	at_bessel_Kx();

void	at_c_bessel_Kx();

void	at_z_bessel_Kx();



	/* Chapter 10 */



void    at_i_write_matrix();

void    at_f_write_matrix();

void    at_f_write_matrix_1();

void    at_f_write_matrix_2();

void    at_d_write_matrix();

void    at_d_write_matrix_1();

void    at_d_write_matrix_2();

void    at_c_write_matrix();

void    at_c_write_matrix_1();

void    at_c_write_matrix_2();

void    at_z_write_matrix();

void    at_z_write_matrix_1();

void    at_z_write_matrix_2();



	/* Chapter 11 */



	/* Chapter 12 */



void    at_random_uniform();

void    at_f_random_uniform();

void    at_d_random_uniform();



void    at_random_normal();

void    at_f_random_normal();

void    at_d_random_normal();



void    at_random_exponential();

void    at_f_random_exponential();

void    at_d_random_exponential();





void    at_random_poisson();



void	at_random_gamma();

void	at_f_random_gamma();

void	at_d_random_gamma();



void	at_random_beta();

void	at_f_random_beta();

void	at_d_random_beta();



void	at_ranks();

void	at_f_ranks();

void	at_d_ranks();



void	at_simple_statistics();

void	at_f_simple_statistics();

void	at_d_simple_statistics();



void	at_chi_squared_test();

void	at_f_chi_squared_test();

void	at_d_chi_squared_test();

 

void	at_covariances();

void	at_f_covariances();

void	at_d_covariances();



void    at_regression();

void    at_f_regression();

void    at_d_regression();



void    at_poly_regression();

void    at_f_poly_regression();

void    at_d_poly_regression();



void    at_beta_cdf();

void    at_f_beta_cdf();

void    at_d_beta_cdf();



void    at_bivariate_normal_cdf();

void    at_f_bivariate_normal_cdf();

void    at_d_bivariate_normal_cdf();



void    at_beta_inverse_cdf();

void    at_f_beta_inverse_cdf();

void    at_d_beta_inverse_cdf();



void	at_table_oneway();

void	at_f_table_oneway();

void	at_d_table_oneway();



	/* Chapter 13 */



void    at_version();



void    at_output_file();


void    at_constant();

void    at_f_constant();

void    at_d_constant();



void    at_error_options();



void    at_error_code();



void    at_i_machine();



void    at_machine();

void    at_f_machine();

void    at_d_machine();



void	at_norm();

void	at_f_norm();

void	at_d_norm();

void	at_i_norm();



/*	SPARSE MATRIX MODULE 	*/



void	at_lin_sol_posdef_coordinate();

void	at_f_lin_sol_posdef_coordinate();

void	at_d_lin_sol_posdef_coordinate();

void	at_c_lin_sol_posdef_coordinate();

void	at_z_lin_sol_posdef_coordinate();



void    at_lin_sol_gen_coordinate ();

void    at_f_lin_sol_gen_coordinate ();

void    at_d_lin_sol_gen_coordinate ();

void    at_c_lin_sol_gen_coordinate ();

void    at_z_lin_sol_gen_coordinate ();



void    at_lin_sol_posdef_band ();

void    at_f_lin_sol_posdef_band ();

void    at_d_lin_sol_posdef_band ();

void    at_c_lin_sol_posdef_band ();

void    at_z_lin_sol_posdef_band ();



void	at_lin_sol_gen_band ();

void	at_f_lin_sol_gen_band ();

void	at_d_lin_sol_gen_band ();

void	at_c_lin_sol_gen_band ();

void	at_z_lin_sol_gen_band ();



void	at_lin_sol_def_cg ();

void	at_f_lin_sol_def_cg ();

void	at_d_lin_sol_def_cg ();



void  at_lin_sol_gen_min_residual ();

void  at_f_lin_sol_gen_min_residual ();

void  at_d_lin_sol_gen_min_residual ();



void	at_mat_mul_rect_band ();

void	at_f_mat_mul_rect_band ();

void	at_d_mat_mul_rect_band ();

void	at_c_mat_mul_rect_band ();

void	at_z_mat_mul_rect_band ();



void	at_mat_mul_rect_coordinate ();

void	at_f_mat_mul_rect_coordinate ();

void	at_d_mat_mul_rect_coordinate ();

void	at_c_mat_mul_rect_coordinate ();

void	at_z_mat_mul_rect_coordinate ();



/*      SPARSE MATRIX DATA      */

 

void PROTO (imsl_f_NOS1, (Mint*, Mint*, char*, Mint**, Mint**, Mfloat**));

void PROTO (imsl_f_662_BUS, (Mint*, Mint*, char*, Mint**, Mint**, Mfloat**));

void PROTO (imsl_f_IMPCOL_A, (Mint*, Mint*, char*, Mint**, Mint**, Mfloat**));

void PROTO (imsl_f_IMPCOL_B, (Mint*, Mint*, char*, Mint**, Mint**, Mfloat**));

void PROTO (imsl_f_IMPCOL_C, (Mint*, Mint*, char*, Mint**, Mint**, Mfloat**));

void PROTO (imsl_f_IMPCOL_D, (Mint*, Mint*, char*, Mint**, Mint**, Mfloat**));

void PROTO (imsl_f_IMPCOL_E, (Mint*, Mint*, char*, Mint**, Mint**, Mfloat**));



#ifdef DOUBLE

#define	    imsl_f_check_vector imsl_d_check_vector 

#define	    imsl_c_check_vector imsl_z_check_vector 

#define	    imsl_f_check_matrix imsl_d_check_matrix 

#define	    imsl_c_check_matrix imsl_z_check_matrix 

#define	    imsl_f_check_scalar imsl_d_check_scalar 

#define	    imsl_c_check_scalar imsl_z_check_scalar 



	/* Chapter 1 */



#define     at_f_lin_sol_gen                at_d_lin_sol_gen

#define     at_c_lin_sol_gen                at_z_lin_sol_gen

#define	    at_f_lin_least_squares_gen	    at_d_lin_least_squares_gen

#define     at_f_lin_sol_nonnegdef         at_d_lin_sol_nonnegdef

#define     at_f_lin_sol_posdef             at_d_lin_sol_posdef

#define     at_c_lin_sol_posdef		    at_z_lin_sol_posdef

#define     at_f_lin_svd_gen		    at_d_lin_svd_gen

#define     at_c_lin_svd_gen		    at_z_lin_svd_gen

#define	    at_f_lin_cgsol_posdef	    at_d_lin_cgsol_posdef

#define     at_f_mat_mul_rect               at_d_mat_mul_rect

#define     at_c_mat_mul_rect               at_z_mat_mul_rect

#define     at_f_dm_band_gen		    at_d_dm_band_gen

#define     at_f_lin_lsq_lin_constraints    at_d_lin_lsq_lin_constraints



	/* Chapter 2 */

#define at_f_eig_sym                    at_d_eig_sym

#define at_f_eig_gen                    at_d_eig_gen

#define at_c_eig_gen                    at_z_eig_gen

#define at_c_eig_herm                   at_z_eig_herm

#define	at_f_eig_symgen			at_d_eig_symgen

       /* Exhaustive Test */

#define	et_f_eig_symgen			et_d_eig_symgen



	/* Chapter 3 */

#define at_f_cub_spline_interp          at_d_cub_spline_interp

#define at_f_cub_spline_interp_e_cnd    at_d_cub_spline_interp_e_cnd

#define at_f_cub_spline_interp_shape    at_d_cub_spline_interp_shape

#define at_f_cub_spline_value           at_d_cub_spline_value

#define at_f_cub_spline_integral        at_d_cub_spline_integral

#define at_f_spline_interp              at_d_spline_interp 

#define at_f_spline_integral            at_d_spline_integral

#define at_f_spline_value               at_d_spline_value

#define at_f_spline_knots               at_d_spline_knots

#define at_f_spline_2d_interp           at_d_spline_2d_interp 

#define at_f_spline_2d_integral         at_d_spline_2d_integral

#define at_f_spline_2d_value            at_d_spline_2d_value

#define at_f_spline_2d_least_squares    at_d_spline_2d_least_squares

#define at_f_spline_least_squares       at_d_spline_least_squares

#define at_f_poly_least_squares         at_d_poly_least_squares

#define at_f_user_fcn_least_squares     at_d_user_fcn_least_squares

#define at_f_cub_spline_smooth          at_d_cub_spline_smooth

#define at_f_scattered_2d_interp        at_d_scattered_2d_interp

#define at_f_radial_scattered_fit	at_d_radial_scattered_fit

#define at_f_radial_evaluate		at_d_radial_evaluate

#define at_f_spline_lsq_constrained     at_d_spline_lsq_constrained





	/* Chapter 4 */



#define at_f_int_fcn_sing               at_d_int_fcn_sing

#define at_f_int_fcn_sing_pts           at_d_int_fcn_sing_pts

#define at_f_int_fcn_inf                at_d_int_fcn_inf

#define at_f_int_fcn_trig               at_d_int_fcn_trig

#define at_f_int_fcn_alg_log            at_d_int_fcn_alg_log

#define at_f_int_fcn_fourier            at_d_int_fcn_fourier

#define at_f_int_fcn_smooth             at_d_int_fcn_smooth

#define at_f_int_fcn_cauchy             at_d_int_fcn_cauchy

#define at_f_int_fcn_hyper_rect         at_d_int_fcn_hyper_rect

#define     at_f_int_fcn                    at_d_int_fcn

#define     at_f_int_fcn_2d                 at_d_int_fcn_2d

#define     at_f_gauss_quad_rule            at_d_gauss_quad_rule





	/* Chapter 5 */



#define at_f_ode_runge_kutta            at_d_ode_runge_kutta

#define at_f_ode_adams_gear             at_d_ode_adams_gear



	/* Chapter 6 */



#define     at_f_fft_real                   at_d_fft_real

#define at_f_fft_real_init              at_d_fft_real_init

#define at_c_fft_complex_init           at_z_fft_complex_init

#define at_c_fft_complex                at_z_fft_complex

#define at_c_fft_2d_complex             at_z_fft_2d_complex

#define	    at_f_convolution		    at_d_convolution

#define	    at_c_convolution		    at_z_convolution



	/* Chapter 7 */



#define at_f_zeros_poly                 at_d_zeros_poly

#define at_c_zeros_poly                 at_z_zeros_poly

#define at_f_zeros_fcn                  at_d_zeros_fcn

#define at_f_zeros_sys_eqn              at_d_zeros_sys_eqn





	/* Chapter 8 */



#define at_f_min_uncon                  at_d_min_uncon

#define at_f_min_uncon_deriv            at_d_min_uncon_deriv

#define at_f_min_uncon_multivar         at_d_min_uncon_multivar

#define at_f_lin_prog                   at_d_lin_prog

#define at_f_quadratic_prog             at_d_quadratic_prog

#define at_f_min_con_nonlin             at_d_min_con_nonlin

#define at_f_nonlin_least_squares       at_d_nonlin_least_squares



	/* Chapter 9 */



#define at_f_erf                        at_d_erf

#define     at_f_erfc                       at_d_erfc

#define     at_f_erf_inverse                at_d_erf_inverse

#define     at_f_erfc_inverse               at_d_erfc_inverse

#define     at_f_erfci                      at_d_erfci

#define     at_f_gamma_incomplete           at_d_gamma_incomplete

#define     at_f_gamma                      at_d_gamma

#define     at_f_beta                       at_d_beta

#define     at_f_log_gamma                  at_d_log_gamma

#define     at_f_log_beta                   at_d_log_beta

#define     at_f_beta_incomplete            at_d_beta_incomplete

#define     at_f_chi_squared_cdf	    at_d_chi_squared_cdf

#define     at_f_chi_squared_inverse_cdf    at_d_chi_squared_inverse_cdf

#define     at_f_t_cdf                      at_d_t_cdf

#define     at_f_t_inverse_cdf              at_d_t_inverse_cdf

#define     at_f_F_cdf                      at_d_F_cdf

#define     at_f_F_inverse_cdf              at_d_F_inverse_cdf

#define     at_f_normal_cdf		    at_d_normal_cdf

#define     at_f_normal_inverse_cdf	    at_d_normal_inverse_cdf

#define     at_f_gamma_cdf		    at_d_gamma_cdf

#define     at_f_hypergeometric_cdf         at_d_hypergeometric_cdf

#define     at_f_poisson_cdf                at_d_poisson_cdf

#define     at_f_binomial_cdf		    at_d_binomial_cdf

#define at_f_bessel_J0                  at_d_bessel_J0

#define at_f_bessel_J1                  at_d_bessel_J1

#define at_f_bessel_I0                  at_d_bessel_I0

#define at_f_bessel_I1                  at_d_bessel_I1

#define at_f_bessel_Y0                  at_d_bessel_Y0

#define at_f_bessel_Y1                  at_d_bessel_Y1

#define at_f_bessel_K0                  at_d_bessel_K0

#define at_f_bessel_K1                  at_d_bessel_K1

#define at_c_bessel_Ix                  at_z_bessel_Ix

#define at_c_bessel_Kx                  at_z_bessel_Kx

#define at_c_bessel_Jx                  at_z_bessel_Jx

#define at_c_bessel_Yx                  at_z_bessel_Yx







	/* Chapter 10 */



#define    at_f_write_matrix                at_d_write_matrix

#define    at_f_write_matrix_1              at_d_write_matrix_1

#define    at_f_write_matrix_2              at_d_write_matrix_2

#define    at_c_write_matrix                at_z_write_matrix

#define    at_c_write_matrix_1              at_z_write_matrix_1

#define    at_c_write_matrix_2              at_z_write_matrix_2



	/* Chapter 11 */



	/* Chapter 12 */



#define     at_f_random_uniform             at_d_random_uniform

#define     at_f_random_normal              at_d_random_normal

#define     at_f_random_exponential         at_d_random_exponential

#define	    at_f_random_gamma		    at_d_random_gamma

#define	    at_f_random_beta		    at_d_random_beta

#define     at_f_ranks                      at_d_ranks

#define     at_f_simple_statistics          at_d_simple_statistics

#define     at_f_chi_squared_test           at_d_chi_squared_test

#define     at_f_covariances                at_d_covariances

#define     at_f_regression                 at_d_regression

#define     at_f_poly_regression            at_d_poly_regression

#define     at_f_beta_cdf                   at_d_beta_cdf

#define     at_f_beta_inverse_cdf           at_d_beta_inverse_cdf

#define     at_f_bivariate_normal_cdf       at_d_bivariate_normal_cdf

#define     at_f_table_oneway		    at_d_table_oneway

	/* Chapter 13 */

#define     at_f_norm			    at_d_norm

#define     at_f_sort			    at_d_sort

#define     at_f_constant                   at_d_constant

#define     at_f_machine                    at_d_machine



/*	SPARSE MATRIX MODULE 	*/



#define at_f_lin_sol_posdef_coordinate	at_d_lin_sol_posdef_coordinate

#define at_c_lin_sol_posdef_coordinate	at_z_lin_sol_posdef_coordinate

#define at_f_lin_sol_gen_coordinate	at_d_lin_sol_gen_coordinate

#define at_c_lin_sol_gen_coordinate	at_z_lin_sol_gen_coordinate

#define at_f_lin_sol_posdef_band	at_d_lin_sol_posdef_band

#define at_c_lin_sol_posdef_band	at_z_lin_sol_posdef_band

#define at_f_lin_sol_gen_band		at_d_lin_sol_gen_band

#define at_c_lin_sol_gen_band		at_z_lin_sol_gen_band

#define at_f_lin_sol_def_cg		at_d_lin_sol_def_cg

#define at_f_lin_sol_gen_min_residual	at_d_lin_sol_gen_min_residual

#define at_f_mat_mul_rect_band		at_d_mat_mul_rect_band

#define at_c_mat_mul_rect_band		at_z_mat_mul_rect_band

#define at_f_mat_mul_rect_coordinate	at_d_mat_mul_rect_coordinate

#define at_c_mat_mul_rect_coordinate	at_z_mat_mul_rect_coordinate



#endif   /* DOUBLE */



enum Imsl_AT_keyword {

    IMSL_AT_NO_PRINT	    = 100001,

    IMSL_AT_ERMAX	    = 100002,

    IMSL_AT_OPTION	    = 100003,

    IMSL_AT_CHECKSUM_FILE   = 100004,

#if defined (COMPUTER_APLC)

    IMSL_AT_CHECK_SUM_FILE  = 100004,

#else

    IMSL_AT_CHECK_SUM_FILE  = IMSL_AT_CHECKSUM_FILE,

#endif

    IMSL_AT_CLEAR_CHECK_SUM = 100005,

    IMSL_AT_CHECK_SUM	    = 100006,

    IMSL_AT_OLD_TEST	    = 100007,

#if defined (COMPUTER_APLC)

    IMSL_AT_SET_ERMAX       = 100002,

#else

    IMSL_AT_SET_ERMAX	    = IMSL_AT_ERMAX,

#endif

    IMSL_AT_COL_DIM	    = IMSL_A_COL_DIM,

    IMSL_AT_COL_LABELS	    = IMSL_COL_LABELS,

    IMSL_AT_ROW_LABELS	    = IMSL_ROW_LABELS,

    IMSL_AT_FORMAT	    = IMSL_WRITE_FORMAT,

    IMSL_AT_TRANSPOSE	    = IMSL_TRANSPOSE

};





#endif

