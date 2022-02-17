#ifndef IMSL_FD_H

#define IMSL_FD_H



#if defined(COMPUTER_APLC) || defined(IMSL_MACHINE_80X86) || defined(IMSL_MACHINE_NT)

#define imsl_ifnan		imsl_difnan

#endif



#define imsl_f_lin_least_squares_gen imsl_d_lin_least_squares_gen

#define imsl_f_lin_sol_gen	    imsl_d_lin_sol_gen

#define imsl_f_lin_sol_posdef	    imsl_d_lin_sol_posdef

#define imsl_c_lin_sol_posdef       imsl_z_lin_sol_posdef

#define imsl_c_lin_sol_gen          imsl_z_lin_sol_gen

#define imsl_f_lin_sol_nonnegdef    imsl_d_lin_sol_nonnegdef

#define imsl_c_lin_sol_full	    imsl_z_lin_sol_full

#define imsl_f_lin_svd_gen	    imsl_d_lin_svd_gen

#define imsl_c_lin_svd_gen	    imsl_z_lin_svd_gen

#define imsl_f_mat_mul_rect         imsl_d_mat_mul_rect

#define imsl_c_mat_mul_rect         imsl_z_mat_mul_rect

#define imsl_f_lin_lsq_lin_constraints	imsl_d_lin_lsq_lin_constraints

#define imsl_crgrg		    imsl_dcrgrg

#define imsl_svrgn		    imsl_dsvrgn

#define imsl_lfsrg		    imsl_dlfsrg

#define imsl_l2trg		    imsl_dl2trg

#define imsl_l2rrr		    imsl_dl2rrr

#define imsl_l2qrr                  imsl_dl2qrr

#define imsl_lqrsl                  imsl_dlqrsl

#define imsl_cgeru		    imsl_zgeru



#define imsl_f_lin_sol_gen_coordinate 	imsl_d_lin_sol_gen_coordinate

#define imsl_c_lin_sol_gen_coordinate 	imsl_z_lin_sol_gen_coordinate

#define imsl_f_lin_sol_posdef_coordinate imsl_d_lin_sol_posdef_coordinate

#define imsl_c_lin_sol_posdef_coordinate imsl_z_lin_sol_posdef_coordinate

#define imsl_f_lin_sol_gen_band 	imsl_d_lin_sol_gen_band

#define imsl_c_lin_sol_gen_band 	imsl_z_lin_sol_gen_band

#define imsl_f_lin_sol_posdef_band 	imsl_d_lin_sol_posdef_band

#define imsl_c_lin_sol_posdef_band 	imsl_z_lin_sol_posdef_band

#define imsl_f_lin_sol_def_cg		imsl_d_lin_sol_def_cg

#define imsl_f_lin_sol_gen_min_residual imsl_d_lin_sol_gen_min_residual

#define imsl_f_free_numeric_factor	imsl_d_free_numeric_factor

#define imsl_c_free_numeric_factor	imsl_z_free_numeric_factor

#define imsl_f_generate_test_coordinate	imsl_d_generate_test_coordinate

#define imsl_f_generate_test_band	imsl_d_generate_test_band

#define imsl_c_generate_test_coordinate	imsl_z_generate_test_coordinate

#define imsl_c_generate_test_band	imsl_z_generate_test_band

#define imsl_f_mat_mul_rect_coordinate	imsl_d_mat_mul_rect_coordinate

#define imsl_c_mat_mul_rect_coordinate	imsl_z_mat_mul_rect_coordinate

#define imsl_f_mat_mul_rect_band	imsl_d_mat_mul_rect_band

#define imsl_c_mat_mul_rect_band	imsl_z_mat_mul_rect_band





#define imsl_f_eig_sym		imsl_d_eig_sym

#define imsl_f_eig_gen		imsl_d_eig_gen

#define imsl_c_eig_gen          imsl_z_eig_gen

#define imsl_c_eig_herm		imsl_z_eig_herm

#define imsl_f_eig_symgen       imsl_d_eig_symgen

#define imsl_f_geneig_sym_posdef	imsl_d_geneig_sym_posdef



#define imsl_c1sor                      imsl_dc1sor

#define imsl_c2dec                      imsl_dc2dec

#define imsl_c2int                      imsl_dc2int

#define imsl_csder                      imsl_dcsder

#define imsl_csitg                      imsl_dcsitg

#define imsl_ppder                      imsl_dppder

#define imsl_p3der                      imsl_dp3der

#define imsl_ppitg                      imsl_dppitg

#define imsl_f_cub_spline_interp        imsl_d_cub_spline_interp

#define imsl_f_cub_spline_interp_e_cnd  imsl_d_cub_spline_interp_e_cnd

#define imsl_f_cub_spline_value         imsl_d_cub_spline_value

#define imsl_f_cub_spline_integral      imsl_d_cub_spline_integral

#define imsl_f_cub_spline_interp_shape  imsl_d_cub_spline_interp_shape

#define imsl_f_ppoly_create             imsl_d_ppoly_create

#define imsl_f_ppoly_print              imsl_d_ppoly_print

#define imsl_f_spline_print             imsl_d_spline_print

#define imsl_b2der                      imsl_db2der

#define imsl_b3der                      imsl_db3der

#define imsl_b4der                      imsl_db4der

#define imsl_b2int                      imsl_db2int

#define imsl_b3int                      imsl_db3int

#define imsl_b4int                      imsl_db4int

#define imsl_b5int                      imsl_db5int

#define imsl_b2itg                      imsl_db2itg

#define imsl_b3itg                      imsl_db3itg

#define imsl_b4itg                      imsl_db4itg

#define imsl_b5itg                      imsl_db5itg

#define imsl_b2nak                      imsl_db2nak

#define imsl_c1not                      imsl_dc1not

#define imsl_crbrb                      imsl_dcrbrb

#define imsl_l2lrb                      imsl_dl2lrb

#define imsl_l2crb                      imsl_dl2crb

#define imsl_l2trb                      imsl_dl2trb

#define imsl_lfirb                      imsl_dlfirb

#define imsl_lfsrb                      imsl_dlfsrb

#define imsl_nr1rb                      imsl_dnr1rb

#define imsl_stbsv                      imsl_dstbsv

#define imsl_b2opk                      imsl_db2opk

#define imsl_b3opk                      imsl_db3opk

#define imsl_b4opk                      imsl_db4opk

#define imsl_f_spline_create            imsl_d_spline_create

#define imsl_f_spline_interp            imsl_d_spline_interp

#define imsl_f_spline_value             imsl_d_spline_value

#define imsl_f_spline_integral          imsl_d_spline_integral

#define imsl_f_spline_knots             imsl_d_spline_knots

#define imsl_f_spline_2d_interp         imsl_d_spline_2d_interp

#define imsl_f_spline_2d_value          imsl_d_spline_2d_value

#define imsl_f_spline_2d_integral       imsl_d_spline_2d_integral

#define imsl_b22dr                      imsl_db22dr

#define imsl_b32dr                      imsl_db32dr

#define imsl_b22ig                      imsl_db22ig

#define imsl_b32ig                      imsl_db32ig

#define imsl_b22in                      imsl_db22in

#define imsl_b32in                      imsl_db32in

#define imsl_b42in                      imsl_db42in

#define imsl_f_spline_least_squares     imsl_d_spline_least_squares

#define imsl_f_spline_2d_least_squares  imsl_d_spline_2d_least_squares

#define imsl_b2lsq                      imsl_db2lsq  

#define imsl_b3lsq                      imsl_db3lsq  

#define imsl_b4lsq                      imsl_db4lsq  

#define imsl_b5lsq                      imsl_db5lsq  

#define imsl_b6lsq                      imsl_db6lsq  

#define imsl_b2ls2                      imsl_db2ls2  

#define imsl_b3ls2                      imsl_db3ls2  

#define imsl_b4ls2                      imsl_db4ls2  

#define imsl_b5ls2                      imsl_db5ls2  

#define imsl_b6ls2                      imsl_db6ls2  

#define imsl_b7ls2                      imsl_db7ls2  

#define imsl_f_poly_least_squares       imsl_d_poly_least_squares

#define imsl_f_user_fcn_least_squares   imsl_d_user_fcn_least_squares

#define imsl_f2lsq                      imsl_df2lsq

#define imsl_b2vls                      imsl_db2vls  

#define imsl_b3vls                      imsl_db3vls  

#define imsl_b4vls                      imsl_db4vls  

#define imsl_b5vls                      imsl_db5vls  

#define imsl_b6vls                      imsl_db6vls  

#define imsl_b7vls                      imsl_db7vls  

#define imsl_b8vls                      imsl_db8vls  

#define imsl_b2cpp                      imsl_db2cpp  

#define imsl_b3cpp                      imsl_db3cpp  

#define imsl_c2scv                      imsl_dc2scv

#define imsl_c3scv                      imsl_dc3scv 

#define imsl_c4scv                      imsl_dc4scv

#define imsl_c5scv                      imsl_dc5scv

#define imsl_c2smh                      imsl_dc2smh

#define imsl_c3smh                      imsl_dc3smh

#define imsl_c4smh                      imsl_dc4smh

#define imsl_f_cub_spline_smooth        imsl_d_cub_spline_smooth

#define imsl_f_scattered_2d_interp      imsl_d_scattered_2d_interp

#define Imsl_f_radial_basis_fit     Imsl_d_radial_basis_fit 

#define imsl_f_radial_scattered_fit imsl_d_radial_scattered_fit

#define imsl_f_radial_evaluate      imsl_d_radial_evaluate     

#define imsl_f_radial_function      imsl_d_radial_function

#define imsl_f_spline_lsq_constrained imsl_d_spline_lsq_constrained



#define imsl_c21gd			imsl_dc21gd

#define imsl_b21gd			imsl_db21gd

#define imsl_b22gd			imsl_db22gd







#define imsl_f_ode_runge_kutta      imsl_d_ode_runge_kutta

#define imsl_f_ode_runge_kutta_mgr  imsl_d_ode_runge_kutta_mgr

#define imsl_f_ode_adams_gear     imsl_d_ode_adams_gear

#define imsl_f_ode_adams_gear_mgr   imsl_d_ode_adams_gear_mgr



#define imsl_isamax		imsl_idamax

#define imsl_icamax             imsl_izamax

#define imsl_sdot		imsl_ddot

#define imsl_strsv		imsl_dtrsv

#define imsl_sger		imsl_dger

#define imsl_csfrg              imsl_dcsfrg

#define imsl_sasum              imsl_dasum

#define imsl_sgemv              imsl_dgemv

#define imsl_ssyr2              imsl_dsyr2

#define imsl_svcal              imsl_dvcal

#define imsl_svrbn              imsl_dsvrbn

#define imsl_svrgp              imsl_dsvrgp

#define imsl_trnrr              imsl_dtrnrr

#define imsl_snrm2              imsl_dnrm2

#define imsl_scnrm2             imsl_sznrm2

#define imsl_trncr              imsl_dtrncr



#define imsl_f_m1ran		imsl_d_m1ran

#define imsl_c_m1ran		imsl_z_m1ran

#define imsl_f_write_matrix	imsl_d_write_matrix

#define imsl_f_wrrrl	        imsl_d_wrrrl

#define imsl_c_wrcrl	        imsl_z_wrcrl

#define imsl_w13rl_f		imsl_dw13rl_f

#define imsl_fmtx		imsl_dfmtx

#define imsl_w1iss		imsl_dw1iss

#define imsl_c_write_matrix	imsl_z_write_matrix



#define imsl_f_sort		imsl_d_sort

#define imsl_f_vector_norm	imsl_d_vector_norm

#define imsl_c_vector_norm	imsl_z_vector_norm

#define imsl_f_matrix_norm	imsl_d_matrix_norm

#define imsl_f_constant		imsl_d_constant



                /* defined(imsl_fc_convert) is taken to imply that

                   all of the macros in this group are defined. */

#ifndef imsl_fc_convert

#define imsl_fc_convert		imsl_dz_convert

#define imsl_cf_convert		imsl_zd_convert

#define imsl_c_neg		imsl_z_neg

#define imsl_c_add		imsl_z_add

#define imsl_c_sub		imsl_z_sub

#define imsl_c_mul		imsl_z_mul

#define imsl_c_eq		imsl_z_eq

#define imsl_c_aimag		imsl_z_aimag

#define imsl_c_conjg		imsl_z_conjg

#define imsl_c_abs		imsl_z_abs

#endif



#define imsl_c_div		imsl_z_div

#define imsl_c_sqrt		imsl_z_sqrt

#define imsl_c_log		imsl_z_log

#define imsl_c_exp		imsl_z_exp

#define imsl_c_sin		imsl_z_sin

#define imsl_c_cos		imsl_z_cos

#define imsl_c_arg		imsl_z_arg



#define imsl_amach		imsl_dmach

#define imsl_f_machine		imsl_d_machine

#define imsl_inits		imsl_initds

#define imsl_csevl		imsl_dcsevl

#define imsl_f_min		imsl_d_min

#define imsl_f_vmin		imsl_d_vmin

#define imsl_f_max		imsl_d_max

#define imsl_f_vmax		imsl_d_vmax

#define imsl_fi_power		imsl_di_power

#define imsl_ff_power		imsl_dd_power

#define imsl_cf_power		imsl_zd_power

#define imsl_cc_power		imsl_zz_power

#define imsl_f_erf		imsl_d_erf

#define imsl_f_erfi		imsl_d_erfi

#define imsl_f_erfc		imsl_d_erfc

#define imsl_f_erfci		imsl_d_erfci

#define imsl_f_erf_inverse    imsl_d_erf_inverse

#define imsl_f_erfc_inverse     imsl_d_erfc_inverse

#define imsl_f_random_normal    imsl_d_random_normal

#define imsl_f_normal_cdf	imsl_d_normal_cdf

#define imsl_f_binomial_cdf     imsl_d_binomial_cdf

#define imsl_f_chi_squared_cdf  imsl_d_chi_squared_cdf

#define imsl_f_simple_statistics imsl_d_simple_statistics

#define imsl_f_regression       imsl_d_regression

#define imsl_f_poly_regression  imsl_d_poly_regression

#define imsl_f_table_oneway	imsl_d_table_oneway

#define imsl_f_bivariate_normal_cdf	imsl_d_bivariate_normal_cdf

#define imsl_f_beta_cdf		imsl_d_beta_cdf

#define imsl_f_beta_inverse_cdf	imsl_d_beta_inverse_cdf

#define imsl_f_random_exponential	imsl_d_random_exponential

#define imsl_f_random_normal_multivariate	imsl_d_random_normal_multivariate



#define imsl_f_fft_real_init	imsl_d_fft_real_init

#define imsl_f_fft_real		imsl_d_fft_real

#define imsl_c_fft_complex_init	imsl_z_fft_complex_init

#define imsl_c_fft_complex	imsl_z_fft_complex

#define imsl_c_fft_2d_complex	imsl_z_fft_2d_complex

#define imsl_f_convolution	imsl_d_convolution

#define imsl_c_convolution	imsl_z_convolution

#define imsl_f2trf		imsl_df2trf

#define imsl_f2trb		imsl_df2trb

#define imsl_f2tcf		imsl_df2tcf

#define imsl_f3tcf		imsl_df3tcf

#define imsl_f5tcf		imsl_df5tcf

#define imsl_f6tcf		imsl_df6tcf

#define imsl_f7tcf		imsl_df7tcf

#define imsl_f8tcf		imsl_df8tcf

#define imsl_f2tcb		imsl_df2tcb

#define imsl_f3tcb		imsl_df3tcb

#define imsl_f5tcb		imsl_df5tcb

#define imsl_f6tcb		imsl_df6tcb

#define imsl_f7tcb		imsl_df7tcb

#define imsl_f8tcb		imsl_df8tcb



#define imsl_f_int_fcn		imsl_d_int_fcn

#define imsl_f_int_fcn_sing	imsl_d_int_fcn_sing

#define imsl_f_int_fcn_sing_pts	imsl_d_int_fcn_sing_pts

#define imsl_f_int_fcn_inf	imsl_d_int_fcn_inf

#define imsl_f_int_fcn_trig	imsl_d_int_fcn_trig

#define imsl_f_int_fcn_alg_log	imsl_d_int_fcn_alg_log

#define imsl_f_int_fcn_fourier	imsl_d_int_fcn_fourier

#define imsl_f_int_fcn_smooth	imsl_d_int_fcn_smooth

#define imsl_f_int_fcn_cauchy	imsl_d_int_fcn_cauchy

#define imsl_f_int_fcn_hyper_rect	imsl_d_int_fcn_hyper_rect

#define imsl_f_int_fcn_2d       imsl_d_int_fcn_2d

#define imsl_f_gauss_quad_rule  imsl_d_gauss_quad_rule

#define imsl_g2rcf              imsl_dg2rcf

#define imsl_g3rcf              imsl_dg3rcf

#define imsl_g4rcf              imsl_dg4rcf

#define imsl_g2rul              imsl_dg2rul

#define imsl_reccf              imsl_dreccf





#define imsl_qdag		imsl_dqdag

#define imsl_q2ag		imsl_dq2ag

#define imsl_q9ag		imsl_dq9ag

#define imsl_q10g		imsl_dq10g

#define imsl_q4ng		imsl_dq4ng

#define imsl_q3agi		imsl_dq3agi

#define imsl_q3awo		imsl_dq3awo

#define imsl_q4awo		imsl_dq4awo

#define imsl_q7awo		imsl_dq7awo

#define imsl_q8awo		imsl_dq8awo



#define imsl_f_zeros_poly           imsl_d_zeros_poly

#define imsl_c_zeros_poly           imsl_z_zeros_poly

#define imsl_f_zeros_fcn            imsl_d_zeros_fcn

#define imsl_f_zeros_sys_eqn        imsl_d_zeros_sys_eqn



#define imsl_f_min_uncon            imsl_d_min_uncon

#define imsl_f_min_uncon_deriv      imsl_d_min_uncon_deriv

#define imsl_f_min_uncon_multivar   imsl_d_min_uncon_multivar

#define imsl_f_min_con_nonlin       imsl_d_min_con_nonlin

#define imsl_f_quadratic_prog       imsl_d_quadratic_prog

#define imsl_f_lin_prog             imsl_d_lin_prog

#define imsl_f_nonlin_least_squares imsl_d_nonlin_least_squares

#define imsl_d2prs		    imsl_dd2prs

#define imsl_q2rog		    imsl_dq2rog



#define imsl_f_int_prog             imsl_d_int_prog



#define tf1tri_c		uf1tri_c

#define tf1trf_c		uf1trf_c

#define tf1trb_c		uf1trb_c

#define tf1tci_c		uf1tci_c

#define imsl_fftci		imsl_dfftci

#define imsl_f3tci		imsl_df3tci

#define tf1tcf_c		uf1tcf_c

#define imsl_fftcf		imsl_dfftcf

#define imsl_f2tcf		imsl_df2tcf

#define imsl_f3tcf		imsl_df3tcf

#define imsl_f4tcf		imsl_df4tcf

#define imsl_f5tcf		imsl_df5tcf

#define imsl_f6tcf		imsl_df6tcf

#define imsl_f7tcf		imsl_df7tcf

#define imsl_f8tcf		imsl_df8tcf

#define tf1tcb_c		uf1tcb_c

#define imsl_fftcb		imsl_dfftcb

#define imsl_f2tcb		imsl_df2tcb

#define imsl_f3tcb		imsl_df3tcb

#define imsl_f4tcb		imsl_df4tcb

#define imsl_f5tcb		imsl_df5tcb

#define imsl_f6tcb		imsl_df6tcb

#define imsl_f7tcb		imsl_df7tcb

#define imsl_f8tcb		imsl_df8tcb

#define tf1t2d_c		uf1t2d_c 

#define imsl_fft2d		imsl_dfft2d 

#define imsl_f2t2d		imsl_df2t2d

#define tf1t2b_c		uf1t2b_c 

#define imsl_fft2b		imsl_dfft2b 

#define imsl_f2t2b		imsl_df2t2b

#define aimag			daimag

#define cexp			dcexp

#define cfmul			cmul

#define cfdiv			cdiv

#define cfpowf			cpowf

#define ccopy			zcopy

#define imsl_sdot		imsl_ddot

#define imsl_wrcrn		imsl_dwrcrn



#define imsl_f_least_squares_fit    imsl_d_least_squares_fit

#define imsl_r2ivn		    imsl_dr2ivn

#define imsl_f_ranks		    imsl_d_ranks

#define imsl_f_normal_inverse_cdf   imsl_d_normal_inverse_cdf

#define imsl_f_random_uniform	    imsl_d_random_uniform

#define imsl_rnses		    imsl_drnses

#define imsl_f_chi_squared_test	    imsl_d_chi_squared_test

#define imsl_c1wfr		    imsl_dc1wfr

#define imsl_f_covariances	    imsl_d_covariances

#define imsl_f_random_beta	    imsl_d_random_beta

#define imsl_rcovb		    imsl_drcovb



#define imsl_f_bessel_J0	    	imsl_d_bessel_J0

#define imsl_f_bessel_J1		imsl_d_bessel_J1

#define imsl_f_bessel_I0	    	imsl_d_bessel_I0

#define imsl_f_bessel_I1		imsl_d_bessel_I1

#define imsl_f_bessel_Y0	    	imsl_d_bessel_Y0

#define imsl_f_bessel_Y1		imsl_d_bessel_Y1

#define imsl_f_bessel_K0	    	imsl_d_bessel_K0

#define imsl_f_bessel_K1		imsl_d_bessel_K1

#define imsl_f_bessel_exp_I0		imsl_d_bessel_exp_I0

#define imsl_f_bessel_exp_I1		imsl_d_bessel_exp_I1

#define imsl_f_bessel_exp_K0		imsl_d_bessel_exp_K0

#define imsl_f_bessel_exp_K1		imsl_d_bessel_exp_K1



#define imsl_c3is			imsl_dc3is

#define imsl_c_bessel_Ix		imsl_z_bessel_Ix

#define imsl_c_bessel_Kx		imsl_z_bessel_Kx

#define imsl_c_bessel_Jx		imsl_z_bessel_Jx

#define imsl_c_bessel_Yx		imsl_z_bessel_Yx

/* The following fout #define's are here for IMSL/IDL. */

#define imsl_c_bessel_Ix_adr		imsl_z_bessel_Ix_adr

#define imsl_c_bessel_Jx_adr		imsl_z_bessel_Jx_adr

#define imsl_c_bessel_Kx_adr		imsl_z_bessel_Kx_adr

#define imsl_c_bessel_Yx_adr		imsl_z_bessel_Yx_adr



#define imsl_f_log_gamma	    imsl_d_log_gamma

#define imsl_f_gamma_incomplete	    imsl_d_gamma_incomplete

#define imsl_f_gamma		    imsl_d_gamma

#define imsl_r9lgmc		    imsl_d9lgmc

#define imsl_f_chi_squared_inverse_cdf imsl_d_chi_squared_inverse_cdf

#define imsl_f_beta_incomplete	    imsl_d_beta_incomplete

#define imsl_f_log_beta             imsl_d_log_beta

#define imsl_f_t_inverse_cdf	    imsl_d_t_inverse_cdf

#define imsl_f_F_cdf		    imsl_d_F_cdf

#define imsl_f_gamma_cdf	    imsl_d_gamma_cdf

#define imsl_srotmg		    imsl_drotmg

#define imsl_srotm		    imsl_drotm

#define imsl_srot		    imsl_drot

#define imsl_srotg                  imsl_drotg

#define imsl_f_t_cdf		    imsl_d_t_cdf

#define imsl_f_F_inverse_cdf	    imsl_d_F_inverse_cdf

#define imsl_f_hypergeometric_cdf   imsl_d_hypergeometric_cdf

#define imsl_betin		    imsl_dbetin

#define imsl_f_poisson_cdf	    imsl_d_poisson_cdf

#define imsl_f_random_gamma	    imsl_d_random_gamma

#define imsl_sxyz		    imsl_dxyz

#define imsl_c1div		    imsl_dc1div

#define imsl_g1aov		    imsl_dg1aov

#define imsl_c1r		    imsl_dc1r

#define imsl_a1ot		    imsl_da1ot

#define imsl_c1ge0		    imsl_dc1ge0



#define imsl_scopy                      imsl_dcopy

#define imsl_ccopy                      imsl_zcopy 

#define imsl_cscal                      imsl_zscal  

#define imsl_csscal                     imsl_zsscal

#define imsl_cswap                      imsl_zswap

#define imsl_cadd                       imsl_zadd

#define imsl_caxpy                      imsl_zaxpy

#define imsl_cdotc                      imsl_zdotc

#define imsl_cdotu                      imsl_zdotu

#define imsl_cgemv                      imsl_zgemv

#define imsl_cgerc                      imsl_zgerc

#define imsl_cset                       imsl_zset

#define imsl_scasum                     imsl_szasum

#define imsl_ccgcg                      imsl_dccgcg

#define imsl_ccbcb			imsl_dccbcb

#define imsl_ctbsv			imsl_dctbsv

#define imsl_ctrsv			imsl_dctrsv



#ifdef saxpy

#undef saxpy

#define saxpy(N,A,X,INCX,Y,INCY)    AXPY(double,N,A,X,INCX,Y,INCY)

#endif



#ifdef scopy

#undef scopy

#define scopy(N,X,INCX,Y,INCY)	COPY(double,N,X,INCX,Y,INCY)

#endif



#ifdef sswap

#undef sswap

#define sswap(N,X,INCX,Y,INCY)	SWAP(double,N,X,INCX,Y,INCY)

#endif



#ifdef sscal

#undef sscal

#define sscal(N,A,X,INCX)	SCAL(double,N,A,X,INCX)

#endif



#ifdef sset

#undef sset

#define sset(N,A,X,INCX)	SET(double,N,A,X,INCX)

#endif



#ifdef sadd

#undef sadd

#define sadd(N,A,X,INCX)	ADD(double,N,A,X,INCX)

#endif



/*  ET's */

#define et_f_eig_gen                    et_d_eig_gen



/*  Special structures  */

 

#define Imsl_f_ppoly                    Imsl_d_ppoly

#define Imsl_f_spline                   Imsl_d_spline



#endif

