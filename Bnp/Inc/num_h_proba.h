/* ======================================================
   FILENAME:  num_h_proba.h
   
   PURPOSE:   Broydn 
   ====================================================== */

#ifndef NUM_H_PROBA_H
#define NUM_H_PROBA_H

#include "srt_h_types.h"

double norm(double d);
double norm_accurate(double x);
double bivar(double x, double y, double rho);      

double gauss(double x);

double inv_cumnorm_newton(double start,double target);
double inv_cumnorm_fast(double u);
double normal_dist_inv(double p);

double srt_f_intnrm_i1(double a0, double b0, double a1, double b1);
double srt_f_intnrm_i2(double a0, double b0, double a1, double b1);
double srt_f_intnrm_j1(double alpha, double a0, double b0, double a1, double b1);
double srt_f_intnrm_j2(double alpha, double a0, double b0, double a1, double b1);
double 	 srt_f_intnrm_l(double a, double b, double ap, double bp);
double srt_f_intnrm_m1	(double a, double b, double ap, double bp);
double srt_f_intnrm_m2	(double a, double b, double ap, double bp); 

Err num_f_Gaussian_copulas
					(double **Datas, /* Datas[0..nDim - 1][0..nNumDates - 1] */
					 long   nDim,
					 long   nNumDates,  
					 double ***CopulaCorrelation); /* no allocation inside ! */

////////////////////////// Student distribution /////////////////////////////////////


double Student_Dis(
					const int n,
					const double x

				  );

double Student_Integrand( const double x, va_list argptr );
double Student_Distribution( const double x, const double m, const double sigma );
double Std_Student_Distribution( const double x, const double m );
Err Inverse_Student_Diff( double x, double *f, double *df, double *data );
Err Inverse_Std_Student_Distribution2( double *x, const double t, const double m, const double acc );

double Inverse_Student_Distribution( const double t, const unsigned int m, const double acc );

double Student_Dis_Inv(
						const double xacc,
						const double tVal,
						const unsigned int degree);

void  TryStudent_Dis_Inv		(
					    const double x,
						const double tVal,
						const unsigned int degree,
						double *fVal,
						double *f1Val
						);


//////////////////////// Get Linear Correlation from Rank Correlation  //////////////

double RankToLinCorr  (	 
						const double SpRho,
						const int degree,
						const int acc_degree);

////////////////////////////// Student copulas  /////////////////////////////////////////

double CopulaRhoSpearman(
						const double r,
						const double SpRho,
						const int degree,
						const int acc_degree
					);

double IntegrateCopula ( 
						const double r,
						const int n,
						const int acc_degree
					  );

double EvalStudCpl(
					const double		xmin,
					const double		xmax,
					const double		ymin,
					const double		ymax,
					const double		r,
					const int	degree,
					const int	acc_degree);

double EvalStudCplIntegrand(	
						const double x,
						const double y,
						const double r,
						const int n
						);



//--------------------------------------------------------------------------------------------------------------------------------------
//
// Functions for fitting Student Copula
//
//--------------------------------------------------------------------------------------------------------------------------------------
void Stud_Deg_Fit(			 
					unsigned long *m,
					double *linear_corr,				
					double *data1,
					const double m1,
					const double sigma1,
					double *data2,						
					const double m2,
					const double sigma2,
					const unsigned short N,
					double Inv_Stud_acc,
					unsigned long Max_Stud_Degree,		
					unsigned long RankToLin_acc );

void Stud_Bivar_Copula(			 
					unsigned long *m,
					double *llhd,
					double *linear_corr,				
					double *data1,
					const double m1,
					const double sigma1,
					double *data2,						
					const double m2,
					const double sigma2,
					const unsigned short N,
					double Inv_Stud_acc,
					unsigned long Max_Stud_Degree,		
					unsigned long RankToLin_acc );

double log_Stud_ML( unsigned long degree, double linear_corr, double *x1, double *x2, unsigned long NumData );
double Stud_cop_llhd( unsigned long m, double rho, double *CumDensity1, double *CumDensity2, unsigned short N, double acc );
double log_Student_copula_dist( double x1, double x2, unsigned long m, double rho, double acc );

double Student_Density( double x, double m, double sigma );
double Std_Student_Density( double x, double m );
double log_Student_dist( double x, unsigned long m );


double Hill_estimator( double *data, unsigned long k, unsigned long N );
Err Stud_Log_Likelihd( double sigma, double *f, double *df, double *data );
Err Stud_Deg_Diff( double m, double *f, double *df, double *data );


void Student_Simluation( double *data1, 
					    double *data2, 
						const unsigned int N, 
						const unsigned int m, 
						const double rho, 
						const unsigned int m1, 
						const double sigma1, 
						const unsigned int m2, 
						const double sigma2, 
						const double acc,
						long seed );

void BiStudCpl_Variates( double *u1, double *u2, const unsigned int m, const double rho, long *seed );

	void GetStudentCplDev (
						const int degree,
	  					const double *mean_v,
						double **corr_mtx,   
						const long p,
						const int d, 
						double **xa,
						double **ya,
						const long n_pts,
						int n_conv,
						SrtMCSamType  MCType,
						double **res
					   );

double *GetStudentCplDev_Rand (
								const int degree,
	  							const double *mean_v,
								double **sqrt_corr_mtx,   
								const long p,
								const int d, 
								double **xa,
								double **ya,
								const long n_pts,
								int n_conv,
								long *idum,
								double *GaussSample
					          );

#endif



